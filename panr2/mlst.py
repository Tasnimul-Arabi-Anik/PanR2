import os
import shutil
import subprocess
from datetime import datetime

import pandas as pd
import plotly.express as px

from panr2.io import normalize_metadata_aliases, resolve_sample_accessions, write_sample_map_qc
from panr2.runners import find_sequence_files, write_tool_manifest


def _read_mlst_table(path):
    delimiter = "\t" if path.lower().endswith((".tab", ".tsv")) else ","
    try:
        df = pd.read_csv(path, sep=delimiter, dtype=str)
    except Exception:
        df = pd.read_csv(path, sep="\t", dtype=str, header=None)
    lower_cols = {str(col).lower(): col for col in df.columns}
    if {"scheme", "st"}.issubset(lower_cols):
        sample_col = next((lower_cols[key] for key in ["assembly", "assembly_accession", "sample", "file", "#file"] if key in lower_cols), df.columns[0])
        out = pd.DataFrame({
            "sample": df[sample_col],
            "scheme": df[lower_cols["scheme"]],
            "st": df[lower_cols["st"]],
        })
        allele_cols = [col for col in df.columns if col not in {sample_col, lower_cols["scheme"], lower_cols["st"]}]
        out["allele_profile"] = df[allele_cols].fillna("").astype(str).agg(";".join, axis=1) if allele_cols else ""
        return out

    raw = pd.read_csv(path, sep="\t", dtype=str, header=None)
    if raw.shape[1] < 3:
        raw = pd.read_csv(path, sep=",", dtype=str, header=None)
    if raw.shape[1] < 3:
        return pd.DataFrame(columns=["sample", "scheme", "st", "allele_profile"])
    out = pd.DataFrame({
        "sample": raw.iloc[:, 0],
        "scheme": raw.iloc[:, 1],
        "st": raw.iloc[:, 2],
    })
    out["allele_profile"] = raw.iloc[:, 3:].fillna("").astype(str).agg(";".join, axis=1) if raw.shape[1] > 3 else ""
    return out


def read_mlst_tables(mlst_dir, sample_map=None):
    if not os.path.isdir(mlst_dir):
        raise FileNotFoundError(f"MLST directory not found: {mlst_dir}")
    frames = []
    for name in sorted(os.listdir(mlst_dir)):
        if name.lower().endswith((".csv", ".tsv", ".tab")):
            frame = _read_mlst_table(os.path.join(mlst_dir, name))
            if not frame.empty:
                frames.append(frame)
    if not frames:
        raise FileNotFoundError(f"No MLST CSV/TSV/TAB outputs found in {mlst_dir}")
    mlst = pd.concat(frames, ignore_index=True)
    mlst["Assembly Accession"] = resolve_sample_accessions(mlst["sample"], sample_map=sample_map).fillna(mlst["sample"].astype(str))
    mlst["scheme"] = mlst["scheme"].fillna("unknown_scheme").astype(str)
    mlst["st"] = mlst["st"].fillna("unknown").astype(str)
    mlst["sequence_type"] = mlst["scheme"] + ":ST" + mlst["st"].str.replace("^ST", "", regex=True)
    return mlst.drop_duplicates(subset=["Assembly Accession", "scheme", "st"])


def run_mlst(sequence_dir, output_dir, mlst_bin="mlst", force=False):
    executable = shutil.which(mlst_bin)
    if not executable:
        raise FileNotFoundError(
            f"MLST executable not found: {mlst_bin}. Run `panr doctor` or provide --mlst-dir with existing MLST output."
        )
    sequence_files = find_sequence_files(sequence_dir)
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")
    raw_dir = os.path.join(output_dir, "tool_results", "mlst", "raw")
    os.makedirs(raw_dir, exist_ok=True)
    raw_path = os.path.join(raw_dir, "mlst.tsv")
    if force or not os.path.exists(raw_path):
        completed = subprocess.run([executable, *sequence_files], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())
        with open(raw_path, "w") as handle:
            handle.write(completed.stdout)
    version = "available"
    try:
        version_result = subprocess.run([executable, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        version = (version_result.stdout or version_result.stderr or "available").strip().splitlines()[0]
    except Exception:
        pass
    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sequence_dir": sequence_dir,
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "mlst",
            "executable": executable,
            "version": version,
            "runs": [{
                "database": "pubmlst",
                "database_sequences": "",
                "database_date": "",
                "results": raw_path,
                "summary": raw_path,
                "status": "completed",
            }],
        }],
    }
    manifest_paths = write_tool_manifest(output_dir, manifest)
    return {"mlst_dir": raw_dir, "raw_table": raw_path, "manifest": manifest_paths}


def analyze_mlst(ncbi_clean_path, mlst_dir, output_dir, fig_format="png", sample_map=None):
    mlst = read_mlst_tables(mlst_dir, sample_map=sample_map)
    write_sample_map_qc(output_dir, "mlst", mlst["sample"], mlst["Assembly Accession"], sample_map)
    ncbi = normalize_metadata_aliases(pd.read_csv(ncbi_clean_path, dtype=str))
    merged = ncbi.merge(mlst, on="Assembly Accession", how="left")

    base_dir = os.path.join(output_dir, "mlst")
    analysis_dir = os.path.join(base_dir, "analysis")
    merged_dir = os.path.join(base_dir, "merged_output")
    figures_dir = os.path.join(base_dir, "figures")
    html_dir = os.path.join(figures_dir, "html_files")
    for path in [analysis_dir, merged_dir, figures_dir, html_dir]:
        os.makedirs(path, exist_ok=True)

    merged_path = os.path.join(merged_dir, "mlst_merged.csv")
    summary_path = os.path.join(analysis_dir, "sample_mlst_summary.csv")
    metadata_path = os.path.join(analysis_dir, "mlst_by_metadata.csv")
    feature_summary_path = os.path.join(analysis_dir, "mlst_feature_summary.csv")
    category_summary_path = os.path.join(analysis_dir, "mlst_category_summary.csv")
    qc_path = os.path.join(analysis_dir, "mlst_qc_summary.csv")
    tidy_path = os.path.join(merged_dir, "mlst_tidy.csv")
    burden_path = os.path.join(analysis_dir, "mlst_sample_burden.csv")

    merged.to_csv(merged_path, index=False)
    sample_summary = merged[[
        col for col in [
            "Assembly Accession", "Assembly BioSample Accession", "scheme", "st", "sequence_type",
            "Geographic Location", "Country", "Continent", "Subcontinent", "Collection Date", "Collection_Year",
            "Host_SD", "Host_Rank", "Host_Genus", "Host_Species", "Sample_Type_SD", "Isolation_Source_SD",
            "Environment_Medium_SD", "Organism Name", "Genus", "Species"
        ]
        if col in merged.columns
    ]].copy()
    sample_summary.to_csv(summary_path, index=False)

    rows = []
    for meta_col in [
        "Geographic Location", "Country", "Continent", "Subcontinent", "Collection Date", "Collection_Year",
        "Host", "Host_SD", "Host_Rank", "Host_Genus", "Host_Species", "Isolation Source",
        "Isolation_Source_SD", "Isolation_Source_SD_Broad", "Sample_Type_SD", "Sample_Type_SD_Broad",
        "Environment_Medium_SD", "Environment_Medium_SD_Broad", "Environment_Broad_Scale_SD",
        "Environment_Local_Scale_SD", "Host_Disease_SD", "Host_Health_State_SD", "Organism Name",
        "Genus", "Species"
    ]:
        if meta_col not in merged.columns:
            continue
        for (group, sequence_type), sub_df in merged.dropna(subset=["sequence_type"]).groupby([meta_col, "sequence_type"], dropna=False):
            rows.append({
                "metadata_variable": meta_col,
                "group": group,
                "sequence_type": sequence_type,
                "sample_count": int(len(sub_df)),
            })
    pd.DataFrame(rows).to_csv(metadata_path, index=False)

    tidy_cols = [
        col for col in [
            "Assembly Accession", "Assembly BioSample Accession", "sequence_type", "Geographic Location",
            "Country", "Continent", "Subcontinent", "Collection Date", "Collection_Year", "Host_SD",
            "Host_Rank", "Host_Genus", "Host_Species", "Sample_Type_SD", "Isolation_Source_SD",
            "Environment_Medium_SD", "Organism Name", "Genus", "Species"
        ] if col in merged.columns
    ]
    tidy = merged.dropna(subset=["sequence_type"])[tidy_cols].copy()
    tidy["feature_id"] = tidy["sequence_type"]
    tidy["presence"] = 1
    tidy["identity"] = 100.0
    tidy.to_csv(tidy_path, index=False)

    burden = sample_summary.copy()
    burden["mlst_feature_count"] = burden["sequence_type"].notna().astype(int)
    burden.to_csv(burden_path, index=False)

    typed = merged.dropna(subset=["sequence_type"]).copy()
    if typed.empty:
        pd.DataFrame(columns=["feature_id", "feature_category", "present_samples", "prevalence_percentage"]).to_csv(feature_summary_path, index=False)
        pd.DataFrame(columns=["feature_category", "present_samples", "prevalence_percentage", "unique_features"]).to_csv(category_summary_path, index=False)
    else:
        feature_summary = typed.groupby("sequence_type")["Assembly Accession"].nunique().reset_index(name="present_samples")
        feature_summary["feature_id"] = feature_summary["sequence_type"]
        feature_summary["feature_category"] = typed.groupby("sequence_type")["scheme"].first().reindex(feature_summary["sequence_type"]).values
        feature_summary["prevalence_percentage"] = (feature_summary["present_samples"] / max(len(merged), 1)) * 100
        feature_summary[["feature_id", "feature_category", "present_samples", "prevalence_percentage"]].to_csv(feature_summary_path, index=False)
        category_summary = typed.groupby("scheme").agg(
            present_samples=("Assembly Accession", "nunique"),
            unique_features=("sequence_type", "nunique"),
        ).reset_index().rename(columns={"scheme": "feature_category"})
        category_summary["prevalence_percentage"] = (category_summary["present_samples"] / max(len(merged), 1)) * 100
        category_summary.to_csv(category_summary_path, index=False)
    pd.DataFrame([
        {"metric": "total_metadata_samples", "value": int(len(merged)), "detail": "Samples in ncbi_clean.csv"},
        {"metric": "typed_samples", "value": int(merged["sequence_type"].notna().sum()), "detail": "Samples with an MLST assignment"},
        {"metric": "untyped_samples", "value": int(merged["sequence_type"].isna().sum()), "detail": "Samples without an MLST assignment"},
    ]).to_csv(qc_path, index=False)

    if not sample_summary.empty and "sequence_type" in sample_summary.columns:
        plot_df = sample_summary.dropna(subset=["sequence_type"])
        if not plot_df.empty:
            fig = px.histogram(plot_df, y="sequence_type", color="Continent" if "Continent" in plot_df.columns else None,
                               title="MLST sequence type distribution")
            fig.write_html(os.path.join(html_dir, "ST_distribution_by_country.html"))

    return {
        "merged": merged_path,
        "tidy": tidy_path,
        "sample_burden": burden_path,
        "feature_summary": feature_summary_path,
        "category_summary": category_summary_path,
        "qc_summary": qc_path,
        "sample_mlst_summary": summary_path,
        "mlst_by_metadata": metadata_path,
        "st_feature_burden_summary": "",
        "html_index": os.path.join(figures_dir, "index.html"),
    }


def finalize_mlst_analysis(output_dir, mlst_outputs, cross_database_outputs=None):
    if not mlst_outputs:
        return mlst_outputs
    sample_summary = pd.read_csv(mlst_outputs["sample_mlst_summary"])
    burden_path = (cross_database_outputs or {}).get("sample_integrated_feature_burden")
    analysis_dir = os.path.join(output_dir, "mlst", "analysis")
    html_dir = os.path.join(output_dir, "mlst", "figures", "html_files")
    os.makedirs(analysis_dir, exist_ok=True)
    os.makedirs(html_dir, exist_ok=True)
    out_path = os.path.join(analysis_dir, "st_feature_burden_summary.csv")
    if burden_path and os.path.exists(burden_path):
        burden = pd.read_csv(burden_path)
        merged = sample_summary.merge(burden, on="Assembly Accession", how="left", suffixes=("", "_feature"))
        count_cols = [col for col in merged.columns if col.endswith("_feature_count") or col in {"total_mobileome_count", "total_feature_count", "mobility_associated_amr_score"}]
        summary = merged.dropna(subset=["sequence_type"]).groupby(["scheme", "st", "sequence_type"], dropna=False)[count_cols].agg(["count", "mean", "median", "max"])
        summary.columns = ["_".join(col).strip("_") for col in summary.columns]
        summary = summary.reset_index()
        summary.to_csv(out_path, index=False)
        if "amr_feature_count" in merged.columns:
            fig = px.box(merged.dropna(subset=["sequence_type"]), x="sequence_type", y="amr_feature_count",
                         points="all", title="AMR burden by MLST sequence type")
            fig.write_html(os.path.join(html_dir, "ST_vs_AMR_burden.html"))
    else:
        pd.DataFrame().to_csv(out_path, index=False)
    mlst_outputs["st_feature_burden_summary"] = out_path
    return mlst_outputs
