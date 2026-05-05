import logging
import os

import pandas as pd

from panr2.io import extract_assembly_accessions, normalize_metadata_aliases, read_table_auto


def write_input_qc_report(ncbi_clean_path, abricate_summary_file, abricate_results_file, output_dir):
    """Validate PanR2 inputs and write machine-readable QC reports."""
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)

    checks = []

    def add_check(check, status, detail, value=None):
        checks.append({
            "check": check,
            "status": status,
            "detail": detail,
            "value": "" if value is None else value,
        })

    required_ncbi_cols = [
        "Assembly Accession",
        "Assembly BioSample Accession",
        "Geographic Location",
        "Continent",
        "Subcontinent",
        "Collection Date",
    ]
    required_results_cols = ["GENE", "RESISTANCE"]

    ncbi_df = normalize_metadata_aliases(read_table_auto(ncbi_clean_path))
    summary_df = read_table_auto(abricate_summary_file)
    results_df = read_table_auto(abricate_results_file)

    missing_ncbi_cols = [col for col in required_ncbi_cols if col not in ncbi_df.columns]
    if missing_ncbi_cols:
        add_check("ncbi_required_columns", "FAIL", "Missing required NCBI metadata columns.", ";".join(missing_ncbi_cols))
    else:
        add_check("ncbi_required_columns", "PASS", "All required NCBI metadata columns are present.")

    if "#File" in summary_df.columns:
        summary_file_col = "#File"
    elif "#FILE" in summary_df.columns:
        summary_file_col = "#FILE"
    else:
        summary_file_col = None
        add_check("abricate_summary_file_column", "FAIL", "ABRicate summary is missing #FILE/#File column.")

    if summary_file_col:
        add_check("abricate_summary_file_column", "PASS", "ABRicate summary file column is present.", summary_file_col)
        summary_accessions = extract_assembly_accessions(summary_df[summary_file_col])
        missing_summary_accessions = int(summary_accessions.isna().sum())
        if missing_summary_accessions:
            add_check(
                "abricate_summary_accession_parse",
                "WARN",
                "Some ABRicate summary rows do not contain parseable GCF/GCA accessions.",
                missing_summary_accessions,
            )
        else:
            add_check("abricate_summary_accession_parse", "PASS", "All ABRicate summary rows have parseable GCF/GCA accessions.")
    else:
        summary_accessions = pd.Series(dtype=object)

    missing_results_cols = [col for col in required_results_cols if col not in results_df.columns]
    if missing_results_cols:
        add_check("abricate_results_required_columns", "FAIL", "ABRicate results is missing required columns.", ";".join(missing_results_cols))
    else:
        add_check("abricate_results_required_columns", "PASS", "ABRicate results contains GENE and RESISTANCE columns.")

    ncbi_accessions = set(ncbi_df.get("Assembly Accession", pd.Series(dtype=object)).dropna().astype(str))
    summary_accession_set = set(summary_accessions.dropna().astype(str))

    unmatched_ncbi = sorted(ncbi_accessions - summary_accession_set)
    unmatched_summary = sorted(summary_accession_set - ncbi_accessions)
    if unmatched_ncbi:
        add_check("ncbi_without_abricate_summary", "WARN", "NCBI assemblies without matching ABRicate summary rows.", len(unmatched_ncbi))
    else:
        add_check("ncbi_without_abricate_summary", "PASS", "Every NCBI assembly has an ABRicate summary row.")

    if unmatched_summary:
        add_check("abricate_summary_without_ncbi", "WARN", "ABRicate summary rows without matching NCBI metadata.", len(unmatched_summary))
    else:
        add_check("abricate_summary_without_ncbi", "PASS", "Every ABRicate summary accession has NCBI metadata.")

    if "NUM_FOUND" in summary_df.columns:
        num_found = pd.to_numeric(summary_df["NUM_FOUND"], errors="coerce").fillna(0)
        zero_hit_count = int((num_found == 0).sum())
        if zero_hit_count:
            add_check("zero_arg_hit_samples", "WARN", "Assemblies with zero ARG hits in ABRicate summary.", zero_hit_count)
        else:
            add_check("zero_arg_hit_samples", "PASS", "All ABRicate summary rows report at least one ARG hit.")
    else:
        add_check("abricate_summary_num_found", "WARN", "ABRicate summary has no NUM_FOUND column; zero-hit sample check skipped.")

    if not missing_results_cols:
        missing_resistance = int(results_df["RESISTANCE"].isna().sum() + (results_df["RESISTANCE"].astype(str).str.strip() == "").sum())
        unique_genes = int(results_df["GENE"].dropna().nunique())
        unique_resistance = int(results_df["RESISTANCE"].dropna().nunique())
        add_check("unique_arg_genes", "INFO", "Unique ARG genes in ABRicate results.", unique_genes)
        add_check("unique_resistance_classes", "INFO", "Unique resistance classes in ABRicate results.", unique_resistance)
        if missing_resistance:
            add_check("missing_resistance_annotations", "WARN", "ABRicate result rows missing resistance class annotation.", missing_resistance)
        else:
            add_check("missing_resistance_annotations", "PASS", "All ABRicate result rows have resistance class annotations.")

    add_check("ncbi_assembly_count", "INFO", "Assemblies in ncbi_clean.csv.", len(ncbi_df))
    add_check("abricate_summary_row_count", "INFO", "Rows in ABRicate summary file.", len(summary_df))
    add_check("abricate_results_row_count", "INFO", "Rows in ABRicate results file.", len(results_df))

    checks_df = pd.DataFrame(checks)
    checks_path = os.path.join(qc_dir, "panr2_input_qc.csv")
    checks_df.to_csv(checks_path, index=False)

    unmatched_path = os.path.join(qc_dir, "panr2_unmatched_accessions.csv")
    pd.DataFrame(
        [{"source": "ncbi_without_abricate_summary", "Assembly Accession": acc} for acc in unmatched_ncbi]
        + [{"source": "abricate_summary_without_ncbi", "Assembly Accession": acc} for acc in unmatched_summary]
    ).to_csv(unmatched_path, index=False)

    summary_path = os.path.join(qc_dir, "panr2_input_qc_summary.txt")
    with open(summary_path, "w") as handle:
        handle.write("PanR2 input QC summary\n")
        handle.write(f"NCBI file: {ncbi_clean_path}\n")
        handle.write(f"ABRicate summary: {abricate_summary_file}\n")
        handle.write(f"ABRicate results: {abricate_results_file}\n\n")
        for status in ["FAIL", "WARN", "PASS", "INFO"]:
            subset = checks_df[checks_df["status"] == status]
            handle.write(f"{status}: {len(subset)}\n")
            for row in subset.itertuples(index=False):
                value = f" ({row.value})" if str(row.value) else ""
                handle.write(f"- {row.check}: {row.detail}{value}\n")
            handle.write("\n")

    if (checks_df["status"] == "FAIL").any():
        failed = ", ".join(checks_df.loc[checks_df["status"] == "FAIL", "check"].tolist())
        raise ValueError(f"Input QC failed: {failed}. See {checks_path}")

    logging.info(f"Input QC report saved to {checks_path}")
    return checks_path
