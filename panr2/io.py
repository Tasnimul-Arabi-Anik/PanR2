import logging
import os

import pandas as pd


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
    sep = "\t" if path.endswith(".tab") else ","
    return pd.read_csv(path, sep=sep)

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

def unique_input_files(files):
    """Return one file per basename, preferring CSV over TAB when both exist."""
    selected = {}
    for path in sorted(files):
        base, ext = os.path.splitext(path)
        current = selected.get(base)
        if current is None or ext.lower() == ".csv":
            selected[base] = path
    return sorted(selected.values())

def load_and_merge_data(ncbi_clean_path, abricate_summary_file):
    """Load and merge NCBI and Abricate data."""
    try:
        ncbi_clean_df = pd.read_csv(ncbi_clean_path)
        abricate_summary_df = pd.read_csv(abricate_summary_file)
        
        # Extract Assembly Accession from the ABRicate file column.
        if "#File" in abricate_summary_df.columns:
            file_col = "#File"
        elif "#FILE" in abricate_summary_df.columns:
            file_col = "#FILE"
        else:
            raise ValueError("ABRicate summary must contain a #FILE or #File column.")

        abricate_summary_df["Assembly Accession"] = extract_assembly_accessions(abricate_summary_df[file_col])
        if abricate_summary_df["Assembly Accession"].isna().all():
            raise ValueError("No GCF_ or GCA_ assembly accessions were found in the ABRicate summary file column.")
        
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
