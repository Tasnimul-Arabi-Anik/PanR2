import logging
import os

import pandas as pd

from panr2.io import extract_assembly_accessions, read_table_auto, unique_input_files


def _find_first(files, token):
    matches = [path for path in files if token in os.path.basename(path).lower()]
    return sorted(matches)[0] if matches else None


def find_abricate_feature_files(feature_dir):
    """Find one ABRicate summary/results pair in a feature directory."""
    if not feature_dir:
        return None, None
    if not os.path.isdir(feature_dir):
        raise FileNotFoundError(f"Feature directory not found: {feature_dir}")
    files = unique_input_files([
        os.path.join(feature_dir, name)
        for name in os.listdir(feature_dir)
        if name.endswith((".csv", ".tab"))
    ])
    summary = _find_first(files, "summary")
    results = _find_first(files, "results")
    if not summary:
        raise FileNotFoundError(f"No ABRicate summary CSV/TAB file found in {feature_dir}")
    if not results:
        logging.warning(f"No ABRicate results CSV/TAB file found in {feature_dir}; annotation summaries will be limited.")
    return summary, results


def _read_feature_table(path):
    if not path:
        return pd.DataFrame()
    return read_table_auto(path)


def _feature_columns(summary_df):
    if "NUM_FOUND" in summary_df.columns:
        start = summary_df.columns.get_loc("NUM_FOUND") + 1
        return list(summary_df.columns[start:])
    file_cols = {"#FILE", "#File", "Assembly Accession"}
    return [col for col in summary_df.columns if col not in file_cols]


def _feature_category(results_df, feature_id, mode):
    if results_df.empty or "GENE" not in results_df.columns:
        return "Unknown"
    sub_df = results_df[results_df["GENE"].astype(str) == str(feature_id)]
    if sub_df.empty:
        return "Unknown"
    if mode == "virulence":
        for col in ["PRODUCT", "SEQUENCE", "DATABASE", "ACCESSION"]:
            if col in sub_df.columns and sub_df[col].notna().any():
                value = str(sub_df[col].dropna().iloc[0]).strip()
                if value:
                    return value
    if mode == "plasmid":
        # PlasmidFinder feature IDs usually carry the useful replicon/Inc label.
        return str(feature_id)
    return "Unknown"


def analyze_abricate_feature_database(ncbi_clean_path, feature_dir, output_dir, feature_type, mode, min_identity=0.0):
    """Analyze an optional ABRicate-style database such as VFDB or PlasmidFinder."""
    summary_path, results_path = find_abricate_feature_files(feature_dir)
    ncbi_df = pd.read_csv(ncbi_clean_path)
    summary_df = _read_feature_table(summary_path)
    results_df = _read_feature_table(results_path) if results_path else pd.DataFrame()

    if "#File" in summary_df.columns:
        file_col = "#File"
    elif "#FILE" in summary_df.columns:
        file_col = "#FILE"
    else:
        raise ValueError(f"{feature_type} summary must contain #FILE or #File column")

    summary_df["Assembly Accession"] = extract_assembly_accessions(summary_df[file_col])
    merged_df = ncbi_df.merge(summary_df, on="Assembly Accession", how="left").fillna("0")
    feature_cols = _feature_columns(merged_df)
    if not feature_cols:
        raise ValueError(f"No feature columns found in {summary_path}")

    output_feature_dir = os.path.join(output_dir, feature_type)
    os.makedirs(output_feature_dir, exist_ok=True)

    numeric = merged_df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    if min_identity > 0:
        numeric = numeric.where(numeric >= min_identity, 0)
        merged_df.loc[:, feature_cols] = numeric
    if "NUM_FOUND" in merged_df.columns:
        merged_df["NUM_FOUND"] = (numeric > 0).sum(axis=1).astype(int)

    merged_path = os.path.join(output_feature_dir, f"{feature_type}_merged.csv")
    merged_df.to_csv(merged_path, index=False)

    id_vars = list(merged_df.columns[:merged_df.columns.get_loc("NUM_FOUND") + 1]) if "NUM_FOUND" in merged_df.columns else list(ncbi_df.columns)
    tidy = merged_df.melt(id_vars=id_vars, value_vars=feature_cols, var_name="feature_id", value_name="identity")
    tidy["identity"] = pd.to_numeric(tidy["identity"], errors="coerce").fillna(0)
    tidy["presence"] = (tidy["identity"] > 0).astype(int)
    tidy["feature_type"] = feature_type
    tidy["feature_category"] = tidy["feature_id"].apply(lambda value: _feature_category(results_df, value, mode))
    tidy_path = os.path.join(output_feature_dir, f"{feature_type}_tidy.csv")
    tidy.to_csv(tidy_path, index=False)

    sample_col = "Assembly BioSample Accession" if "Assembly BioSample Accession" in tidy.columns else "Assembly Accession"
    total_samples = max(tidy[sample_col].nunique(), 1)
    present = tidy[tidy["presence"] == 1].copy()

    if present.empty:
        feature_summary = pd.DataFrame(columns=["feature_id", "feature_category", "present_samples", "prevalence_percentage", "mean_identity", "max_identity"])
        category_summary = pd.DataFrame(columns=["feature_category", "present_samples", "prevalence_percentage", "unique_features"])
    else:
        feature_summary = present.groupby(["feature_id", "feature_category"], as_index=False).agg(
            present_samples=(sample_col, "nunique"),
            mean_identity=("identity", "mean"),
            max_identity=("identity", "max"),
        )
        feature_summary["prevalence_percentage"] = (feature_summary["present_samples"] / total_samples) * 100
        feature_summary = feature_summary.sort_values(["prevalence_percentage", "present_samples", "feature_id"], ascending=[False, False, True])

        category_summary = present.groupby("feature_category", as_index=False).agg(
            present_samples=(sample_col, "nunique"),
            unique_features=("feature_id", "nunique"),
        )
        category_summary["prevalence_percentage"] = (category_summary["present_samples"] / total_samples) * 100
        category_summary = category_summary.sort_values(["prevalence_percentage", "unique_features", "feature_category"], ascending=[False, False, True])

    feature_summary_path = os.path.join(output_feature_dir, f"{feature_type}_feature_summary.csv")
    category_summary_path = os.path.join(output_feature_dir, f"{feature_type}_category_summary.csv")
    feature_summary.to_csv(feature_summary_path, index=False)
    category_summary.to_csv(category_summary_path, index=False)

    sample_meta_cols = [col for col in ["Assembly Accession", "Assembly BioSample Accession", "Geographic Location", "Continent", "Collection Date"] if col in merged_df.columns]
    sample_burden = merged_df[sample_meta_cols].copy() if sample_meta_cols else merged_df[["Assembly Accession"]].copy()
    sample_burden[f"{feature_type}_feature_count"] = (numeric > 0).sum(axis=1).astype(int)
    sample_burden_path = os.path.join(output_feature_dir, f"{feature_type}_sample_burden.csv")
    sample_burden.to_csv(sample_burden_path, index=False)

    logging.info(f"{feature_type} analysis outputs saved to {output_feature_dir}")
    return {
        "feature_type": feature_type,
        "summary": summary_path,
        "results": results_path or "not available",
        "merged": merged_path,
        "tidy": tidy_path,
        "feature_summary": feature_summary_path,
        "category_summary": category_summary_path,
        "sample_burden": sample_burden_path,
    }
