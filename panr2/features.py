import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

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


def _safe_figsize(row_count, base_height=2.0, per_row=0.35, width=9, max_height=14):
    return (width, max(4, min(max_height, base_height + per_row * max(row_count, 1))))


def _write_feature_plots(feature_type, output_feature_dir, fig_format, tidy, feature_summary, category_summary, sample_burden, sample_col):
    plot_dir = os.path.join(output_feature_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    plot_paths = {}
    label = "Virulence" if feature_type == "virulence" else "Plasmid"

    if not feature_summary.empty:
        plot_df = feature_summary.head(25).sort_values("prevalence_percentage", ascending=True)
        plt.figure(figsize=_safe_figsize(len(plot_df), width=10))
        sns.barplot(data=plot_df, x="prevalence_percentage", y="feature_id", color="#4C78A8")
        plt.xlabel("Samples carrying feature (%)")
        plt.ylabel(f"{label} feature")
        plt.title(f"Top {label.lower()} features by prevalence")
        plt.xlim(0, 100)
        plt.tight_layout()
        path = os.path.join(plot_dir, f"{feature_type}_feature_prevalence.{fig_format}")
        plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()
        plot_paths["feature_prevalence_plot"] = path

    if not category_summary.empty:
        plot_df = category_summary.head(25).sort_values("prevalence_percentage", ascending=True)
        plt.figure(figsize=_safe_figsize(len(plot_df), width=10))
        sns.barplot(data=plot_df, x="prevalence_percentage", y="feature_category", color="#59A14F")
        plt.xlabel("Samples carrying category (%)")
        ylabel = "VFDB product/category" if feature_type == "virulence" else "Plasmid replicon/category"
        plt.ylabel(ylabel)
        plt.title(f"Top {label.lower()} categories by prevalence")
        plt.xlim(0, 100)
        plt.tight_layout()
        path = os.path.join(plot_dir, f"{feature_type}_category_prevalence.{fig_format}")
        plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()
        plot_paths["category_prevalence_plot"] = path

    count_col = f"{feature_type}_feature_count"
    if count_col in sample_burden.columns and not sample_burden.empty:
        burden_df = sample_burden.copy()
        burden_df[count_col] = pd.to_numeric(burden_df[count_col], errors="coerce").fillna(0)
        if "Continent" in burden_df.columns and burden_df["Continent"].nunique() > 1:
            plt.figure(figsize=(9, 5))
            order = burden_df.groupby("Continent")[count_col].median().sort_values(ascending=False).index
            sns.boxplot(data=burden_df, x="Continent", y=count_col, order=order, color="#D8DEE9")
            sns.stripplot(data=burden_df, x="Continent", y=count_col, order=order, color="#2F4B7C", size=4, jitter=True)
            plt.xlabel("Continent")
            plt.ylabel(f"{label} feature count per sample")
            plt.title(f"{label} feature burden by continent")
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout()
            path = os.path.join(plot_dir, f"{feature_type}_burden_by_continent.{fig_format}")
            plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
            plt.close()
            plot_paths["burden_by_continent_plot"] = path

    present = tidy[tidy["presence"] == 1].copy()
    if not present.empty and tidy["feature_id"].nunique() >= 2 and tidy[sample_col].nunique() >= 2:
        top_features = feature_summary.head(40)["feature_id"].tolist() if not feature_summary.empty else sorted(tidy["feature_id"].unique())[:40]
        matrix = (
            tidy[tidy["feature_id"].isin(top_features)]
            .pivot_table(index=sample_col, columns="feature_id", values="presence", fill_value=0, aggfunc="max")
        )
        if matrix.shape[0] >= 2 and matrix.shape[1] >= 2:
            width = max(6, min(16, 0.35 * matrix.shape[1] + 4))
            height = max(4, min(16, 0.28 * matrix.shape[0] + 3))
            plt.figure(figsize=(width, height))
            sns.heatmap(matrix, cmap="Blues", linewidths=0.2, linecolor="white", cbar_kws={"label": "Presence"})
            plt.xlabel(f"{label} feature")
            plt.ylabel("Sample")
            plt.title(f"{label} feature presence matrix")
            plt.xticks(rotation=45, ha="right")
            plt.yticks(rotation=0)
            plt.tight_layout()
            path = os.path.join(plot_dir, f"{feature_type}_presence_heatmap.{fig_format}")
            plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
            plt.close()
            plot_paths["presence_heatmap"] = path

    return plot_paths


def _write_feature_geographic_summary(feature_type, output_feature_dir, tidy, sample_burden):
    count_col = f"{feature_type}_feature_count"
    rows = []
    for geo_col in ["Geographic Location", "Continent", "Subcontinent"]:
        if geo_col not in sample_burden.columns or count_col not in sample_burden.columns:
            continue
        for region, sub_df in sample_burden.groupby(geo_col, dropna=False):
            region_label = str(region) if str(region) not in ["", "nan", "None"] else "Unknown"
            counts = pd.to_numeric(sub_df[count_col], errors="coerce").fillna(0)
            sample_ids = set(sub_df["Assembly BioSample Accession"] if "Assembly BioSample Accession" in sub_df.columns else sub_df["Assembly Accession"])
            region_tidy = tidy[tidy["presence"].eq(1)]
            sample_key = "Assembly BioSample Accession" if "Assembly BioSample Accession" in region_tidy.columns and "Assembly BioSample Accession" in sub_df.columns else "Assembly Accession"
            region_tidy = region_tidy[region_tidy[sample_key].isin(sample_ids)] if sample_key in region_tidy.columns else pd.DataFrame()
            top_features = ""
            if not region_tidy.empty:
                top_features = ";".join(region_tidy.groupby("feature_id")[sample_key].nunique().sort_values(ascending=False).head(5).index.astype(str))
            rows.append({
                "geographic_level": geo_col,
                "region": region_label,
                "sample_count": int(len(sub_df)),
                "samples_with_feature": int((counts > 0).sum()),
                "mean_feature_count": float(counts.mean()) if len(counts) else 0.0,
                "median_feature_count": float(counts.median()) if len(counts) else 0.0,
                "max_feature_count": int(counts.max()) if len(counts) else 0,
                "top_features": top_features,
            })
    summary = pd.DataFrame(rows)
    path = os.path.join(output_feature_dir, f"{feature_type}_geographic_summary.csv")
    summary.to_csv(path, index=False)
    return path


def analyze_abricate_feature_database(ncbi_clean_path, feature_dir, output_dir, feature_type, mode, min_identity=0.0, fig_format="png"):
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

    sample_meta_cols = [col for col in ["Assembly Accession", "Assembly BioSample Accession", "Geographic Location", "Continent", "Subcontinent", "Collection Date"] if col in merged_df.columns]
    sample_burden = merged_df[sample_meta_cols].copy() if sample_meta_cols else merged_df[["Assembly Accession"]].copy()
    sample_burden[f"{feature_type}_feature_count"] = (numeric > 0).sum(axis=1).astype(int)
    sample_burden_path = os.path.join(output_feature_dir, f"{feature_type}_sample_burden.csv")
    sample_burden.to_csv(sample_burden_path, index=False)

    geographic_summary_path = _write_feature_geographic_summary(feature_type, output_feature_dir, tidy, sample_burden)
    plot_paths = _write_feature_plots(feature_type, output_feature_dir, fig_format, tidy, feature_summary, category_summary, sample_burden, sample_col)

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
        "geographic_summary": geographic_summary_path,
        **plot_paths,
    }
