import logging
import os

import pandas as pd


SEQUENCE_SUFFIXES = (".fasta", ".fna", ".fa", ".fas")

METADATA_ALIAS_SOURCES = {
    "Geographic Location": ["Country", "geographic_location", "country"],
    "Collection Date": ["Collection_Year", "Collection Year", "collection_year"],
    "Host": ["Host_SD", "Host_Cleaned", "Host_Original", "host"],
    "Isolation Source": ["Isolation_Source_SD", "Sample_Type_SD", "Environment_Medium_SD", "isolation_source"],
}


def convert_tab_to_csv(tab_file, csv_file):
    """Convert a .tab file to .csv format.

    Args:
        tab_file (str): Path to the .tab file.
        csv_file (str): Path to save the .csv file.
    """
    try:
        df = pd.read_csv(tab_file, sep="\t")
        df.to_csv(csv_file, index=False)
        logging.info(f"Converted {tab_file} to {csv_file}")
    except Exception as e:
        logging.error(f"Error converting {tab_file} to CSV: {e}")
        raise

def read_table_auto(path):
    """Read CSV or TAB input based on file extension."""
    sep = "\t" if path.endswith((".tab", ".tsv")) else ","
    try:
        return pd.read_csv(path, sep=sep)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def normalize_metadata_aliases(df):
    """Add legacy PanR2 metadata columns from richer FetchM2-style fields when needed."""
    if df is None or df.empty:
        return df
    metadata = df.copy()
    for target, sources in METADATA_ALIAS_SOURCES.items():
        if target not in metadata.columns:
            metadata[target] = ""
        target_values = metadata[target].fillna("").astype(str).str.strip()
        for source in sources:
            if source not in metadata.columns:
                continue
            source_values = metadata[source].fillna("").astype(str).str.strip()
            missing = target_values.eq("") | target_values.str.lower().isin({"nan", "none", "0"})
            if missing.any():
                metadata.loc[missing, target] = source_values.loc[missing]
                target_values = metadata[target].fillna("").astype(str).str.strip()
    if "Collection Date" in metadata.columns:
        extracted = metadata["Collection Date"].astype(str).str.extract(r"(\d{4})", expand=False)
        metadata["Collection Date"] = extracted.fillna(metadata["Collection Date"].astype(str))
    if "Collection_Year" in metadata.columns and "Collection Year" not in metadata.columns:
        metadata["Collection Year"] = metadata["Collection_Year"]
    if "Country" in metadata.columns and "Geographic Location" not in metadata.columns:
        metadata["Geographic Location"] = metadata["Country"]
    if "Organism Taxonomic ID" in metadata.columns:
        if "TaxID" not in metadata.columns:
            metadata["TaxID"] = metadata["Organism Taxonomic ID"]
        else:
            missing_taxid = metadata["TaxID"].fillna("").astype(str).str.strip().eq("")
            metadata.loc[missing_taxid, "TaxID"] = metadata.loc[missing_taxid, "Organism Taxonomic ID"]
    if "Organism Name" in metadata.columns:
        organism = metadata["Organism Name"].fillna("").astype(str).str.strip()
        genus = organism.str.split().str[0].replace({"nan": "", "None": "", "0": ""})
        species = organism.str.split().str[:2].str.join(" ").replace({"nan": "", "None": "", "0": ""})
        if "Genus" not in metadata.columns:
            metadata["Genus"] = genus
        else:
            missing_genus = metadata["Genus"].fillna("").astype(str).str.strip().eq("")
            metadata.loc[missing_genus, "Genus"] = genus.loc[missing_genus]
        if "Species" not in metadata.columns:
            metadata["Species"] = species
        else:
            missing_species = metadata["Species"].fillna("").astype(str).str.strip().eq("")
            metadata.loc[missing_species, "Species"] = species.loc[missing_species]
    return metadata


def _sample_key_variants(value):
    text = str(value or "").strip()
    if not text:
        return []
    base = os.path.basename(text)
    if base.endswith(".gz"):
        base = base[:-3]
    variants = [text, base]
    lower_base = base.lower()
    for suffix in SEQUENCE_SUFFIXES:
        if lower_base.endswith(suffix):
            variants.append(base[:-len(suffix)])
            break
    stem, _ = os.path.splitext(base)
    variants.append(stem)
    seen = []
    for item in variants:
        if item and item not in seen:
            seen.append(item)
    return seen


def read_sample_map(sample_map_path):
    """Read optional sample-name to Assembly Accession mapping."""
    if not sample_map_path:
        return {}
    if not os.path.exists(sample_map_path):
        raise FileNotFoundError(f"Sample map not found: {sample_map_path}")
    df = read_table_auto(sample_map_path)
    lower_cols = {str(col).strip().lower().replace("_", " "): col for col in df.columns}
    sample_col = lower_cols.get("sample id") or lower_cols.get("sample") or lower_cols.get("sample name") or lower_cols.get("file")
    accession_col = lower_cols.get("assembly accession") or lower_cols.get("accession")
    if not sample_col or not accession_col:
        raise ValueError("Sample map must contain columns `sample_id` and `Assembly Accession`.")
    mapping = {}
    for _, row in df.iterrows():
        accession = str(row.get(accession_col, "")).strip()
        if not accession:
            continue
        for key in _sample_key_variants(row.get(sample_col, "")):
            mapping.setdefault(key, accession)
    return mapping

def normalize_assembly_accession(value, preserve_version=True):
    """Return a normalized GCF/GCA accession while preserving versions by default."""
    match = pd.Series([value]).astype(str).str.extract(r"((?:GCF|GCA)_\d+(?:\.\d+)?)")[0].iloc[0]
    if pd.isna(match):
        return None
    if preserve_version:
        return match
    return str(match).split(".")[0]

def extract_assembly_accessions(file_series, preserve_version=True):
    """Extract GCF/GCA assembly accessions from file paths or sample names.

    Version suffixes such as `.1` are retained by default. Use
    ``preserve_version=False`` only when intentionally performing
    version-insensitive matching.
    """
    pattern = r"((?:GCF|GCA)_\d+(?:\.\d+)?)"
    accessions = file_series.astype(str).str.extract(pattern)[0]
    if not preserve_version:
        accessions = accessions.astype(str).str.split(".").str[0].where(accessions.notna())
    return accessions


def resolve_sample_accessions(sample_series, sample_map=None, preserve_version=True):
    """Resolve sample/file names to assembly accessions using GCF/GCA extraction and an optional sample map."""
    sample_map = sample_map or {}
    extracted = extract_assembly_accessions(sample_series, preserve_version=preserve_version)
    resolved = []
    for original, accession in zip(sample_series.astype(str), extracted):
        mapped = None
        for key in _sample_key_variants(original):
            if key in sample_map:
                mapped = sample_map[key]
                break
        resolved.append(mapped or accession)
    return pd.Series(resolved, index=sample_series.index)


def write_sample_map_qc(output_dir, source, sample_series, resolved_series, sample_map=None):
    """Write mapping diagnostics for a source that used accession resolution."""
    if not output_dir:
        return None
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    sample_map = sample_map or {}
    rows = []
    for original, resolved in zip(sample_series.astype(str), resolved_series):
        matched_keys = [key for key in _sample_key_variants(original) if key in sample_map]
        extracted = normalize_assembly_accession(original)
        rows.append({
            "source": source,
            "sample_id": original,
            "resolved_assembly_accession": resolved if pd.notna(resolved) else "",
            "mapping_status": "mapped_by_sample_map" if matched_keys else "extracted_accession" if extracted else "unmapped",
            "matched_sample_map_key": matched_keys[0] if matched_keys else "",
        })
    path = os.path.join(qc_dir, f"sample_map_qc_{source}.csv")
    pd.DataFrame(rows).to_csv(path, index=False)
    summary_path = os.path.join(qc_dir, f"sample_map_qc_{source}_summary.txt")
    counts = pd.Series([row["mapping_status"] for row in rows]).value_counts().to_dict()
    duplicates = pd.Series([row["sample_id"] for row in rows]).value_counts()
    duplicate_count = int((duplicates > 1).sum())
    with open(summary_path, "w") as handle:
        handle.write(f"Sample-map QC for {source}\n")
        for status, count in sorted(counts.items()):
            handle.write(f"{status}: {count}\n")
        handle.write(f"duplicate_sample_ids: {duplicate_count}\n")
    return path

def unique_input_files(files):
    """Return one file per basename, preferring CSV over TAB when both exist."""
    selected = {}
    for path in sorted(files):
        base, ext = os.path.splitext(path)
        current = selected.get(base)
        if current is None or ext.lower() == ".csv":
            selected[base] = path
    return sorted(selected.values())

def load_and_merge_data(ncbi_clean_path, abricate_summary_file, sample_map=None, sample_map_qc_dir=None, sample_map_source="abricate"):
    """Load and merge NCBI and Abricate data."""
    try:
        ncbi_clean_df = normalize_metadata_aliases(pd.read_csv(ncbi_clean_path))
        abricate_summary_df = pd.read_csv(abricate_summary_file)
        
        # Extract Assembly Accession from the ABRicate file column.
        if "#File" in abricate_summary_df.columns:
            file_col = "#File"
        elif "#FILE" in abricate_summary_df.columns:
            file_col = "#FILE"
        else:
            raise ValueError("ABRicate summary must contain a #FILE or #File column.")

        abricate_summary_df["Assembly Accession"] = resolve_sample_accessions(abricate_summary_df[file_col], sample_map=sample_map)
        if sample_map_qc_dir:
            write_sample_map_qc(sample_map_qc_dir, sample_map_source, abricate_summary_df[file_col], abricate_summary_df["Assembly Accession"], sample_map)
        if abricate_summary_df["Assembly Accession"].isna().all():
            raise ValueError("No GCF_ or GCA_ assembly accessions were found in the ABRicate summary file column and no sample-map entries matched.")
        
        # Merge dataframes
        merged_df = pd.merge(ncbi_clean_df, abricate_summary_df, on="Assembly Accession", how="left")
        merged_df.fillna("0", inplace=True)
        
        return merged_df
    except Exception as e:
        logging.error(f"Error loading or merging data: {e}")
        raise

def save_merged_data(merged_df, output_dir, output_filename):
    """Save the merged dataframe to a CSV file."""
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, output_filename)
    merged_df.to_csv(output_file, index=False)
    logging.info(f"Merged file saved to: {output_file}")

def get_gene_columns(df):
    """Return ABRicate gene columns from a merged dataframe."""
    if "NUM_FOUND" in df.columns:
        gene_start = df.columns.get_loc("NUM_FOUND") + 1
        return list(df.columns[gene_start:])
    if df.shape[1] > 15:
        return list(df.columns[15:])
    return []

def convert_to_tidy_format(df):
    """Convert the merged dataframe to tidy format."""
    try:
        if df.shape[1] < 14:
            raise ValueError("Not enough columns in the dataframe.")

        # ABRicate summary gene columns begin after NUM_FOUND when present.
        value_vars = get_gene_columns(df)
        if "NUM_FOUND" in df.columns:
            id_vars = list(df.columns[:df.columns.get_loc("NUM_FOUND") + 1])
        else:
            id_vars = list(df.columns[:15])

        if not value_vars:
            raise ValueError("No gene columns found in the merged dataframe.")
        
        # Melt the dataframe
        df_tidy = df.melt(id_vars=id_vars, value_vars=value_vars, var_name="Gene", value_name="Identity")
        df_tidy["Identity"] = pd.to_numeric(df_tidy["Identity"], errors="coerce").fillna(0)
        
        # Add Presence column
        df_tidy["Presence"] = df_tidy["Identity"].apply(lambda x: 0 if x == 0 else 1)
        
        return df_tidy
    except Exception as e:
        logging.error(f"Error converting to tidy format: {e}")
        raise
