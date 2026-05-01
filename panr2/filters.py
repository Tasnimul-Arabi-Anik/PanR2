import logging
import os

import pandas as pd

from panr2.io import get_gene_columns


def apply_analysis_filters(merged_df, min_identity=0.0, drop_unmatched_accessions=False, output_dir=None, base_name="panr2"):
    """Apply optional analysis filters and write a filter report."""
    filtered_df = merged_df.copy()
    report = []

    def add_filter(name, enabled, before, after, detail):
        report.append({
            "filter": name,
            "enabled": bool(enabled),
            "before": before,
            "after": after,
            "removed": before - after,
            "detail": detail,
        })

    before_rows = len(filtered_df)
    if drop_unmatched_accessions:
        file_col = "#File" if "#File" in filtered_df.columns else "#FILE" if "#FILE" in filtered_df.columns else None
        if file_col:
            filtered_df = filtered_df[filtered_df[file_col].astype(str) != "0"].copy()
            detail = f"Kept assemblies with a matching ABRicate summary row using {file_col}."
        else:
            detail = "No ABRicate #FILE/#File column found after merge; no rows dropped."
    else:
        detail = "Filter disabled."
    add_filter("drop_unmatched_accessions", drop_unmatched_accessions, before_rows, len(filtered_df), detail)

    gene_cols = get_gene_columns(filtered_df)
    before_present_calls = 0
    after_present_calls = 0
    if gene_cols:
        numeric_genes = filtered_df[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        before_present_calls = int((numeric_genes > 0).sum().sum())
        if min_identity > 0:
            numeric_genes = numeric_genes.where(numeric_genes >= min_identity, 0)
            filtered_df.loc[:, gene_cols] = numeric_genes
            if "NUM_FOUND" in filtered_df.columns:
                filtered_df["NUM_FOUND"] = (numeric_genes > 0).sum(axis=1).astype(int)
            detail = f"Set ARG calls below {min_identity:g}% identity to absent."
        else:
            detail = "Filter disabled."
        after_present_calls = int((numeric_genes > 0).sum().sum())
    else:
        detail = "No gene columns found; identity filter skipped."
    add_filter("min_identity", min_identity > 0, before_present_calls, after_present_calls, detail)

    if output_dir:
        qc_dir = os.path.join(output_dir, "qc")
        os.makedirs(qc_dir, exist_ok=True)
        report_path = os.path.join(qc_dir, f"{base_name}_filter_report.csv")
        pd.DataFrame(report).to_csv(report_path, index=False)
        logging.info(f"Filter report saved to {report_path}")

    return filtered_df

