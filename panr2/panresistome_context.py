import os
from pathlib import Path

import pandas as pd


SAMPLE_COLUMNS = [
    "Assembly Accession",
    "assembly_accession",
    "sample_id",
    "sequence_accession",
    "sequence_file",
    "genome",
]


def _norm(value):
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    name = Path(text).name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in [".fna", ".fa", ".fasta", ".fas"]:
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
            break
    if name.endswith("_genomic"):
        name = name[:-8]
    return name.lower().replace(".", "_").replace("-", "_")


def _sample_key(row):
    for column in SAMPLE_COLUMNS:
        if column in row and _norm(row.get(column)):
            return _norm(row.get(column))
    return ""


def _read_csv(path):
    if path and os.path.exists(path):
        try:
            return pd.read_csv(path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    return pd.DataFrame()


def _first_existing(output_dir, rel_paths):
    for rel_path in rel_paths:
        path = os.path.join(output_dir, rel_path)
        if os.path.exists(path):
            return path
    return ""


def _write_csv(df, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)
    return path


def _numeric(df, columns):
    for column in columns:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def _load_burden(output_dir, base_name, cross_database_outputs=None):
    cross_database_outputs = cross_database_outputs or {}
    burden = _read_csv(cross_database_outputs.get("sample_integrated_feature_burden", ""))
    if burden.empty:
        burden = _read_csv(os.path.join(output_dir, "ncbi", "analysis", f"{base_name}_sample_resistome_burden.csv"))
    if burden.empty:
        return burden
    burden = burden.copy()
    burden["_sample_key"] = burden.apply(_sample_key, axis=1)
    rename_map = {
        "unique_arg_genes": "amr_feature_count",
        "resistance_class_count": "amr_class_count",
    }
    for old, new in rename_map.items():
        if old in burden.columns and new not in burden.columns:
            burden[new] = burden[old]
    return burden


def _load_context_tables(output_dir):
    paths = {
        "qc_master": _first_existing(output_dir, [
            "qc/qc_master_report.csv",
            "panr2_inputs/qc/qc_master_report.csv",
        ]),
        "excluded": _first_existing(output_dir, [
            "qc/excluded_for_panr2.csv",
            "panr2_inputs/qc/excluded_for_panr2.csv",
        ]),
        "quast": _first_existing(output_dir, [
            "quast/analysis/assembly_qc.csv",
        ]),
        "quast_contract": _first_existing(output_dir, [
            "quast/analysis/panr2_quast_summary.csv",
            "assembly_qc/analysis/panr2_quast_summary.csv",
            "panr2_inputs/assembly_qc/analysis/panr2_quast_summary.csv",
        ]),
        "ani_closest": _first_existing(output_dir, [
            "ani/analysis/closest_genome.csv",
        ]),
        "ani_clusters": _first_existing(output_dir, [
            "ani/analysis/duplicate_clusters.csv",
        ]),
        "ani_contract": _first_existing(output_dir, [
            "ani/analysis/panr2_ani_summary.csv",
            "panr2_inputs/ani/analysis/panr2_ani_summary.csv",
        ]),
        "ani_outliers": _first_existing(output_dir, [
            "ani/analysis/ani_outliers.csv",
        ]),
    }
    tables = {key: _read_csv(path) for key, path in paths.items()}
    return paths, tables


def _merge_context(qc_master, burden):
    if qc_master.empty:
        return pd.DataFrame()
    merged = qc_master.copy()
    merged["_sample_key"] = merged.apply(_sample_key, axis=1)
    if not burden.empty and "_sample_key" in burden.columns:
        burden_columns = [
            column for column in [
                "_sample_key",
                "amr_feature_count",
                "amr_class_count",
                "vfdb_feature_count",
                "plasmid_feature_count",
                "mge_feature_count",
                "is_feature_count",
                "integron_feature_count",
                "ice_feature_count",
                "defense_feature_count",
                "prophage_feature_count",
                "mob_feature_count",
                "kleborate_feature_count",
                "kaptive_feature_count",
                "ectyper_feature_count",
                "serotype_feature_count",
                "sccmec_feature_count",
                "total_mobileome_count",
                "total_feature_count",
                "mobility_associated_amr_score",
            ] if column in burden.columns
        ]
        if burden_columns:
            burden_subset = burden[burden_columns]
            burden_subset = burden_subset[burden_subset["_sample_key"].astype(str) != ""]
            merged = merged.merge(burden_subset.drop_duplicates("_sample_key"), on="_sample_key", how="left")
    numeric_columns = [
        "amr_feature_count",
        "amr_class_count",
        "vfdb_feature_count",
        "plasmid_feature_count",
        "mge_feature_count",
        "is_feature_count",
        "integron_feature_count",
        "ice_feature_count",
        "defense_feature_count",
        "prophage_feature_count",
        "mob_feature_count",
        "kleborate_feature_count",
        "kaptive_feature_count",
        "ectyper_feature_count",
        "serotype_feature_count",
        "sccmec_feature_count",
        "total_mobileome_count",
        "total_feature_count",
        "mobility_associated_amr_score",
        "checkm2_completeness",
        "checkm2_contamination",
        "quast_n50",
        "quast_num_contigs",
        "quast_total_length",
        "ani_closest_ani",
    ]
    merged = _numeric(merged, numeric_columns)
    return merged


def _summarize_duplicate_clusters(clusters, merged):
    if clusters.empty:
        return pd.DataFrame()
    clusters = clusters.copy()
    clusters["_sample_key"] = clusters.apply(_sample_key, axis=1)
    if not merged.empty:
        extra = merged[[
            column for column in [
                "_sample_key",
                "qc_master_status",
                "amr_feature_count",
                "vfdb_feature_count",
                "plasmid_feature_count",
                "total_mobileome_count",
                "mob_feature_count",
                "sccmec_feature_count",
            ] if column in merged.columns
        ]].drop_duplicates("_sample_key")
        extra = extra[extra["_sample_key"].astype(str) != ""]
        clusters = clusters.merge(extra, on="_sample_key", how="left")
    rows = []
    for cluster_id, group in clusters.groupby("ani_cluster", dropna=False):
        row = {
            "ani_cluster": cluster_id,
            "representative": group["representative"].dropna().astype(str).iloc[0] if "representative" in group.columns and group["representative"].notna().any() else "",
            "cluster_size": len(group),
            "members": ";".join(sorted(group["genome"].dropna().astype(str).unique())) if "genome" in group.columns else "",
        }
        if "qc_master_status" in group.columns:
            row["qc_pass_samples"] = int((group["qc_master_status"] == "PASS").sum())
            row["qc_warn_samples"] = int((group["qc_master_status"] == "WARN").sum())
            row["qc_fail_samples"] = int((group["qc_master_status"] == "FAIL").sum())
        for metric in ["amr_feature_count", "vfdb_feature_count", "plasmid_feature_count", "mge_feature_count", "integron_feature_count", "mob_feature_count", "sccmec_feature_count", "total_mobileome_count"]:
            if metric in group.columns:
                row[f"mean_{metric}"] = group[metric].mean()
        rows.append(row)
    return pd.DataFrame(rows)


def _burden_by_ani_cluster(merged):
    if merged.empty or "ani_cluster" not in merged.columns:
        return pd.DataFrame()
    metric_columns = [
        column for column in [
            "amr_feature_count",
            "amr_class_count",
            "vfdb_feature_count",
            "plasmid_feature_count",
            "mge_feature_count",
            "integron_feature_count",
            "mob_feature_count",
            "sccmec_feature_count",
            "kleborate_feature_count",
            "kaptive_feature_count",
            "ectyper_feature_count",
            "serotype_feature_count",
            "total_mobileome_count",
            "total_feature_count",
        ] if column in merged.columns
    ]
    rows = []
    for cluster_id, group in merged.groupby("ani_cluster", dropna=False):
        row = {"ani_cluster": cluster_id, "sample_count": len(group)}
        if "ani_cluster_representative" in group.columns:
            representatives = [value for value in group["ani_cluster_representative"].dropna().astype(str).unique() if value]
            row["representative"] = representatives[0] if representatives else ""
        if "qc_master_status" in group.columns:
            row["qc_pass_samples"] = int((group["qc_master_status"] == "PASS").sum())
            row["qc_warn_samples"] = int((group["qc_master_status"] == "WARN").sum())
            row["qc_fail_samples"] = int((group["qc_master_status"] == "FAIL").sum())
        for metric in metric_columns:
            row[f"mean_{metric}"] = group[metric].mean()
            row[f"median_{metric}"] = group[metric].median()
            row[f"max_{metric}"] = group[metric].max()
        rows.append(row)
    return pd.DataFrame(rows).sort_values("sample_count", ascending=False)


def _qc_feature_correlations(merged):
    if merged.empty:
        return pd.DataFrame()
    qc_metrics = [
        "checkm2_completeness",
        "checkm2_contamination",
        "quast_n50",
        "quast_num_contigs",
        "quast_total_length",
        "ani_closest_ani",
    ]
    burden_metrics = [
        "amr_feature_count",
        "amr_class_count",
        "vfdb_feature_count",
        "plasmid_feature_count",
        "total_mobileome_count",
        "total_feature_count",
    ]
    rows = []
    for qc_metric in [column for column in qc_metrics if column in merged.columns]:
        for burden_metric in [column for column in burden_metrics if column in merged.columns]:
            subset = merged[[qc_metric, burden_metric]].dropna()
            if len(subset) < 3 or subset[qc_metric].nunique() < 2 or subset[burden_metric].nunique() < 2:
                continue
            rows.append({
                "qc_metric": qc_metric,
                "burden_metric": burden_metric,
                "n_samples": len(subset),
                "pearson_r": subset[qc_metric].corr(subset[burden_metric], method="pearson"),
                "spearman_r": subset[qc_metric].corr(subset[burden_metric], method="spearman"),
            })
    return pd.DataFrame(rows)


def _species_summary(merged, closest, outliers):
    source = merged if "ani_species_consistency_status" in merged.columns else closest
    if source.empty or "ani_species_consistency_status" not in source.columns:
        return pd.DataFrame()
    summary = source["ani_species_consistency_status"].fillna("not_recorded").value_counts().reset_index()
    summary.columns = ["ani_species_consistency_status", "sample_count"]
    if not outliers.empty:
        summary["ani_outlier_count"] = len(outliers)
    return summary


def _representatives(clusters, merged):
    if clusters.empty or "representative" not in clusters.columns or "genome" not in clusters.columns:
        return pd.DataFrame()
    reps = clusters[clusters["representative"].astype(str) == clusters["genome"].astype(str)].copy()
    reps["_sample_key"] = reps.apply(_sample_key, axis=1)
    if not merged.empty and "_sample_key" in merged.columns:
        extra = merged[[
            column for column in [
                "_sample_key",
                "qc_master_status",
                "qc_master_fail_reasons",
                "qc_master_warning_reasons",
                "amr_feature_count",
                "total_mobileome_count",
            ] if column in merged.columns
        ]].drop_duplicates("_sample_key")
        extra = extra[extra["_sample_key"].astype(str) != ""]
        reps = reps.merge(extra, on="_sample_key", how="left")
    return reps


def generate_panresistome_context_outputs(output_dir, base_name, cross_database_outputs=None):
    """Generate lightweight PanR2 panels from PanResistome heavy-tool outputs."""
    paths, tables = _load_context_tables(output_dir)
    if all(table.empty for table in tables.values()):
        return {}

    analysis_dir = os.path.join(output_dir, "panresistome_context", "analysis")
    os.makedirs(analysis_dir, exist_ok=True)
    outputs = {}

    burden = _load_burden(output_dir, base_name, cross_database_outputs)
    merged = _merge_context(tables["qc_master"], burden)
    if not merged.empty:
        outputs["qc_context_sample_burden"] = _write_csv(
            merged.drop(columns=[column for column in ["_sample_key"] if column in merged.columns]),
            os.path.join(analysis_dir, "qc_context_sample_burden.csv"),
        )
        outputs["qc_feature_correlation_summary"] = _write_csv(
            _qc_feature_correlations(merged),
            os.path.join(analysis_dir, "qc_feature_correlation_summary.csv"),
        )
        if "qc_master_status" in merged.columns:
            counts = merged["qc_master_status"].fillna("not_recorded").value_counts().reset_index()
            counts.columns = ["qc_master_status", "sample_count"]
            outputs["qc_master_status_summary"] = _write_csv(counts, os.path.join(analysis_dir, "qc_master_status_summary.csv"))

    species_summary = _species_summary(merged, tables["ani_closest"], tables["ani_outliers"])
    if not species_summary.empty:
        outputs["species_consistency_summary"] = _write_csv(species_summary, os.path.join(analysis_dir, "species_consistency_summary.csv"))

    cluster_summary = _summarize_duplicate_clusters(tables["ani_clusters"], merged)
    if not cluster_summary.empty:
        outputs["duplicate_cluster_summary"] = _write_csv(cluster_summary, os.path.join(analysis_dir, "duplicate_cluster_summary.csv"))

    reps = _representatives(tables["ani_clusters"], merged)
    if not reps.empty:
        outputs["representative_samples"] = _write_csv(reps.drop(columns=[column for column in ["_sample_key"] if column in reps.columns]), os.path.join(analysis_dir, "representative_samples.csv"))

    cluster_burden = _burden_by_ani_cluster(merged)
    if not cluster_burden.empty:
        outputs["burden_by_ani_cluster"] = _write_csv(cluster_burden, os.path.join(analysis_dir, "burden_by_ani_cluster.csv"))

    detected = []
    for key, path in paths.items():
        if path:
            detected.append({"context": key, "path": os.path.relpath(path, output_dir)})
    outputs["panresistome_context_manifest"] = _write_csv(pd.DataFrame(detected), os.path.join(analysis_dir, "panresistome_context_manifest.csv"))
    return outputs
