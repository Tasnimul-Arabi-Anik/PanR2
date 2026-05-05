import os

import pandas as pd


METADATA_REPORT_COLUMNS = [
    "Assembly Accession",
    "Assembly BioSample Accession",
    "Organism Name",
    "Genus",
    "Species",
    "TaxID",
    "Geographic Location",
    "Continent",
    "Subcontinent",
    "Collection Date",
    "Host",
    "Host_SD",
    "Host_Rank",
    "Host_Genus",
    "Host_Species",
    "Host_Class",
    "Host_Order",
    "Host_Anatomical_Site_SD",
    "Host_Age_Group_SD",
    "Host_Production_Context_SD",
    "Isolation Source",
    "Isolation_Source_SD",
    "Sample_Type_SD",
    "Environment_Medium_SD",
    "Environment_Broad_Scale_SD",
    "Environment_Local_Scale_SD",
]


def _clean_series(series):
    return series.fillna("").astype(str).str.strip().replace({"nan": "", "None": "", "0": ""})


def write_metadata_reports(ncbi_clean_path, output_dir, min_group_size=5):
    """Write metadata completeness and underpowered-group reports."""
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    completeness_path = os.path.join(qc_dir, "metadata_completeness_report.csv")
    group_sizes_path = os.path.join(qc_dir, "metadata_group_sample_sizes.csv")
    warning_path = os.path.join(qc_dir, "metadata_bias_warning.txt")

    if not ncbi_clean_path or not os.path.exists(ncbi_clean_path):
        pd.DataFrame().to_csv(completeness_path, index=False)
        pd.DataFrame().to_csv(group_sizes_path, index=False)
        with open(warning_path, "w") as handle:
            handle.write("Metadata completeness could not be assessed because ncbi_clean.csv was not available.\n")
        return {"metadata_completeness": completeness_path, "metadata_group_sample_sizes": group_sizes_path, "metadata_bias_warning": warning_path}

    df = pd.read_csv(ncbi_clean_path)
    total = len(df)
    rows = []
    group_rows = []
    warnings = []
    for column in [col for col in METADATA_REPORT_COLUMNS if col in df.columns]:
        clean = _clean_series(df[column])
        present = int((clean != "").sum())
        missing = int(total - present)
        completeness = (present / max(total, 1)) * 100
        if completeness < 50:
            status = "FAIL"
            warnings.append(f"{column} is available for only {completeness:.1f}% of samples. Interpret analyses using this metadata cautiously.")
        elif completeness < 80:
            status = "WARN"
            warnings.append(f"{column} is partially complete ({completeness:.1f}% of samples). Group comparisons may be biased.")
        else:
            status = "PASS"
        rows.append({
            "metadata_column": column,
            "total_samples": total,
            "present_samples": present,
            "missing_samples": missing,
            "completeness_percentage": completeness,
            "unique_nonmissing_values": int(clean[clean != ""].nunique()),
            "status": status,
        })

        if column not in {"Assembly Accession", "Assembly BioSample Accession"}:
            for value, count in clean[clean != ""].value_counts().items():
                group_rows.append({
                    "metadata_column": column,
                    "group": value,
                    "sample_count": int(count),
                    "underpowered": bool(count < min_group_size),
                    "min_group_size": int(min_group_size),
                })

    pd.DataFrame(rows).to_csv(completeness_path, index=False)
    pd.DataFrame(group_rows).to_csv(group_sizes_path, index=False)
    if not warnings:
        warnings.append("No major metadata completeness warnings were detected using the configured thresholds.")
    with open(warning_path, "w") as handle:
        handle.write("\n".join(warnings) + "\n")
    return {"metadata_completeness": completeness_path, "metadata_group_sample_sizes": group_sizes_path, "metadata_bias_warning": warning_path}
