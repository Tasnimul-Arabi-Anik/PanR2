import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu

from panr2.io import extract_assembly_accessions, read_table_auto, unique_input_files


FEATURE_LABELS = {
    "vfdb": "VFDB",
    "plasmidfinder": "PlasmidFinder",
    "mobileelementfinder": "MobileElementFinder",
    "isfinder": "ISfinder",
    "integronfinder": "IntegronFinder",
    "iceberg": "ICEberg",
}


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
    if mode == "plasmid":
        # PlasmidFinder feature IDs usually carry the useful replicon/Inc label.
        return str(feature_id)
    for col in ["PRODUCT", "SEQUENCE", "DATABASE", "ACCESSION"]:
        if col in sub_df.columns and sub_df[col].notna().any():
            value = str(sub_df[col].dropna().iloc[0]).strip()
            if value:
                return value
    return "Unknown"


def _safe_figsize(row_count, base_height=2.0, per_row=0.35, width=9, max_height=14):
    return (width, max(4, min(max_height, base_height + per_row * max(row_count, 1))))


def _feature_label(feature_type, mode):
    if feature_type.lower() in FEATURE_LABELS:
        return FEATURE_LABELS[feature_type.lower()]
    if mode == "virulence" or feature_type.lower() == "vfdb":
        return "VFDB"
    if mode == "plasmid" or feature_type.lower() == "plasmidfinder":
        return "PlasmidFinder"
    return feature_type.replace("_", " ").title()


def _write_feature_plots(feature_type, mode, figures_dir, fig_format, tidy, feature_summary, category_summary, sample_burden, sample_col):
    plot_dir = figures_dir
    html_dir = os.path.join(figures_dir, "html_files")
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(html_dir, exist_ok=True)
    plot_paths = {}
    label = _feature_label(feature_type, mode)

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
        fig = px.bar(plot_df.sort_values("prevalence_percentage", ascending=False), x="prevalence_percentage", y="feature_id", orientation="h",
                     hover_data=[col for col in ["present_samples", "mean_identity", "min_identity", "max_identity", "feature_category"] if col in plot_df.columns],
                     labels={"prevalence_percentage": "Samples carrying feature (%)", "feature_id": f"{label} feature"},
                     title=f"Top {label.lower()} features by prevalence")
        fig.update_layout(yaxis={"categoryorder": "total ascending"})
        html_path = os.path.join(html_dir, f"{feature_type}_feature_prevalence.html")
        fig.write_html(html_path)
        plot_paths["feature_prevalence_html"] = html_path

    if not category_summary.empty:
        plot_df = category_summary.head(25).sort_values("prevalence_percentage", ascending=True)
        plt.figure(figsize=_safe_figsize(len(plot_df), width=10))
        sns.barplot(data=plot_df, x="prevalence_percentage", y="feature_category", color="#59A14F")
        plt.xlabel("Samples carrying category (%)")
        ylabel = "VFDB product/category" if mode == "virulence" else "PlasmidFinder replicon/category"
        plt.ylabel(ylabel)
        plt.title(f"Top {label.lower()} categories by prevalence")
        plt.xlim(0, 100)
        plt.tight_layout()
        path = os.path.join(plot_dir, f"{feature_type}_category_prevalence.{fig_format}")
        plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()
        plot_paths["category_prevalence_plot"] = path
        fig = px.bar(plot_df.sort_values("prevalence_percentage", ascending=False), x="prevalence_percentage", y="feature_category", orientation="h",
                     hover_data=[col for col in ["present_samples", "unique_features"] if col in plot_df.columns],
                     labels={"prevalence_percentage": "Samples carrying category (%)", "feature_category": ylabel},
                     title=f"Top {label.lower()} categories by prevalence")
        fig.update_layout(yaxis={"categoryorder": "total ascending"})
        html_path = os.path.join(html_dir, f"{feature_type}_category_prevalence.html")
        fig.write_html(html_path)
        plot_paths["category_prevalence_html"] = html_path

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
            fig = px.box(burden_df, x="Continent", y=count_col, points="all",
                         labels={count_col: f"{label} feature count per sample"},
                         title=f"{label} feature burden by continent")
            html_path = os.path.join(html_dir, f"{feature_type}_burden_by_continent.html")
            fig.write_html(html_path)
            plot_paths["burden_by_continent_html"] = html_path

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
            fig = px.imshow(matrix, color_continuous_scale="Blues", aspect="auto",
                            labels={"x": f"{label} feature", "y": "Sample", "color": "Presence"},
                            title=f"{label} feature presence matrix")
            html_path = os.path.join(html_dir, f"{feature_type}_presence_heatmap.html")
            fig.write_html(html_path)
            plot_paths["presence_heatmap_html"] = html_path

    identity_df = tidy[tidy["presence"] == 1].copy()
    if not identity_df.empty:
        top_features = feature_summary.head(25)["feature_id"].tolist() if not feature_summary.empty else sorted(identity_df["feature_id"].unique())[:25]
        identity_df = identity_df[identity_df["feature_id"].isin(top_features)]
        if not identity_df.empty:
            plt.figure(figsize=_safe_figsize(identity_df["feature_id"].nunique(), width=11))
            order = identity_df.groupby("feature_id")["identity"].median().sort_values(ascending=True).index
            sns.boxplot(data=identity_df, x="identity", y="feature_id", order=order, color="#E5E7EB")
            sns.stripplot(data=identity_df, x="identity", y="feature_id", order=order, color="#2F4B7C", size=3, jitter=True)
            plt.xlabel("ABRicate identity (%)")
            plt.ylabel(f"{label} feature")
            plt.title(f"{label} feature identity distribution")
            plt.xlim(0, 100)
            plt.tight_layout()
            path = os.path.join(plot_dir, f"{feature_type}_identity_distribution.{fig_format}")
            plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
            plt.close()
            plot_paths["identity_distribution_plot"] = path
            fig = px.box(identity_df, x="identity", y="feature_id", points="all",
                         labels={"identity": "ABRicate identity (%)", "feature_id": f"{label} feature"},
                         title=f"{label} feature identity distribution")
            html_path = os.path.join(html_dir, f"{feature_type}_identity_distribution.html")
            fig.write_html(html_path)
            plot_paths["identity_distribution_html"] = html_path

    return plot_paths


def _write_feature_geographic_summary(feature_type, mode, analysis_dir, figures_dir, tidy, sample_burden):
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
    path = os.path.join(analysis_dir, f"{feature_type}_geographic_summary.csv")
    summary.to_csv(path, index=False)
    html_path = None
    if not summary.empty:
        html_dir = os.path.join(figures_dir, "html_files")
        os.makedirs(html_dir, exist_ok=True)
        label = _feature_label(feature_type, mode)
        fig = px.bar(
            summary.sort_values(["geographic_level", "mean_feature_count"], ascending=[True, False]),
            x="mean_feature_count",
            y="region",
            color="geographic_level",
            orientation="h",
            hover_data=["sample_count", "samples_with_feature", "median_feature_count", "max_feature_count", "top_features"],
            labels={"mean_feature_count": "Mean feature count", "region": "Region", "geographic_level": "Geographic level"},
            title=f"{label} geographic feature burden",
        )
        html_path = os.path.join(html_dir, f"{feature_type}_geographic_burden.html")
        fig.write_html(html_path)
    return path, html_path


def _write_feature_cooccurrence(feature_type, mode, analysis_dir, figures_dir, fig_format, tidy, total_samples, sample_col, top_n=25):
    present = tidy[tidy["presence"] == 1].copy()
    matrix_path = os.path.join(analysis_dir, f"{feature_type}_feature_cooccurrence_matrix.csv")
    pairs_path = os.path.join(analysis_dir, f"{feature_type}_top_feature_pairs.csv")
    plot_paths = {}
    if present.empty:
        pd.DataFrame().to_csv(matrix_path)
        pd.DataFrame(columns=[
            "feature_1", "feature_2", "cooccurring_samples", "cooccurrence_percentage",
            "feature_1_samples", "feature_2_samples", "jaccard_index"
        ]).to_csv(pairs_path, index=False)
        return matrix_path, pairs_path, plot_paths

    presence = (
        present[[sample_col, "feature_id"]]
        .dropna()
        .drop_duplicates()
        .assign(value=1)
        .pivot_table(index=sample_col, columns="feature_id", values="value", fill_value=0, aggfunc="max")
    )
    if presence.empty:
        pd.DataFrame().to_csv(matrix_path)
        pd.DataFrame().to_csv(pairs_path, index=False)
        return matrix_path, pairs_path, plot_paths

    cooc_matrix = presence.T.dot(presence).astype(int)
    cooc_matrix.to_csv(matrix_path)

    sample_counts = presence.sum(axis=0).astype(int).to_dict()
    pair_rows = []
    features = list(cooc_matrix.columns)
    for i, first in enumerate(features):
        for second in features[i + 1:]:
            co_count = int(cooc_matrix.loc[first, second])
            union_count = int(((presence[first] + presence[second]) > 0).sum())
            pair_rows.append({
                "feature_1": first,
                "feature_2": second,
                "cooccurring_samples": co_count,
                "cooccurrence_percentage": (co_count / max(total_samples, 1)) * 100,
                "feature_1_samples": sample_counts[first],
                "feature_2_samples": sample_counts[second],
                "jaccard_index": co_count / union_count if union_count else 0,
            })
    pair_df = pd.DataFrame(pair_rows)
    if not pair_df.empty:
        pair_df = pair_df[pair_df["cooccurring_samples"] > 0]
        pair_df = pair_df.sort_values(["cooccurring_samples", "jaccard_index", "feature_1", "feature_2"], ascending=[False, False, True, True]).head(top_n)
    pair_df.to_csv(pairs_path, index=False)

    if cooc_matrix.shape[0] >= 2:
        top_features = presence.sum(axis=0).sort_values(ascending=False).head(top_n).index.tolist()
        plot_matrix = cooc_matrix.loc[top_features, top_features]
        size = max(5, min(14, 0.45 * len(top_features) + 3))
        label = _feature_label(feature_type, mode)
        plt.figure(figsize=(size, size))
        sns.heatmap(plot_matrix, cmap="Blues", square=True, linewidths=0.3, linecolor="white", cbar_kws={"label": "Co-occurring samples"})
        plt.title(f"{label} feature co-occurrence")
        plt.xlabel(f"{label} feature")
        plt.ylabel(f"{label} feature")
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
        plt.tight_layout()
        path = os.path.join(figures_dir, f"{feature_type}_feature_cooccurrence_heatmap.{fig_format}")
        plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()
        plot_paths["feature_cooccurrence_heatmap"] = path
        html_dir = os.path.join(figures_dir, "html_files")
        os.makedirs(html_dir, exist_ok=True)
        fig = px.imshow(plot_matrix, color_continuous_scale="Blues", aspect="auto",
                        labels={"x": f"{label} feature", "y": f"{label} feature", "color": "Co-occurring samples"},
                        title=f"{label} feature co-occurrence")
        html_path = os.path.join(html_dir, f"{feature_type}_feature_cooccurrence_heatmap.html")
        fig.write_html(html_path)
        plot_paths["feature_cooccurrence_heatmap_html"] = html_path

    return matrix_path, pairs_path, plot_paths


def _write_feature_temporal_summary(feature_type, mode, analysis_dir, figures_dir, sample_burden):
    count_col = f"{feature_type}_feature_count"
    path = os.path.join(analysis_dir, f"{feature_type}_temporal_summary.csv")
    columns = ["collection_year", "sample_count", "samples_with_feature", "mean_feature_count", "median_feature_count", "max_feature_count"]
    if "Collection Date" not in sample_burden.columns or count_col not in sample_burden.columns:
        pd.DataFrame(columns=columns).to_csv(path, index=False)
        return path, None
    df = sample_burden.copy()
    df["collection_year"] = pd.to_numeric(df["Collection Date"], errors="coerce")
    df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(0)
    df = df.dropna(subset=["collection_year"])
    if df.empty:
        pd.DataFrame(columns=columns).to_csv(path, index=False)
        return path, None
    df["collection_year"] = df["collection_year"].astype(int)
    summary = df.groupby("collection_year", as_index=False).agg(
        sample_count=(count_col, "size"),
        samples_with_feature=(count_col, lambda values: int((values > 0).sum())),
        mean_feature_count=(count_col, "mean"),
        median_feature_count=(count_col, "median"),
        max_feature_count=(count_col, "max"),
    )
    summary.to_csv(path, index=False)
    html_dir = os.path.join(figures_dir, "html_files")
    os.makedirs(html_dir, exist_ok=True)
    label = _feature_label(feature_type, mode)
    fig = px.line(summary, x="collection_year", y="mean_feature_count", markers=True,
                  hover_data=["sample_count", "samples_with_feature", "median_feature_count", "max_feature_count"],
                  labels={"collection_year": "Collection year", "mean_feature_count": "Mean feature count"},
                  title=f"{label} temporal feature burden")
    html_path = os.path.join(html_dir, f"{feature_type}_temporal_burden.html")
    fig.write_html(html_path)
    return path, html_path




def _write_feature_group_analysis(feature_type, mode, analysis_dir, figures_dir, fig_format, sample_burden, min_group_size=2):
    count_col = f"{feature_type}_feature_count"
    summary_path = os.path.join(analysis_dir, f"{feature_type}_group_burden_summary.csv")
    overall_path = os.path.join(analysis_dir, f"{feature_type}_group_overall_tests.csv")
    pairwise_path = os.path.join(analysis_dir, f"{feature_type}_group_pairwise_tests.csv")
    plot_paths = {}
    if count_col not in sample_burden.columns:
        pd.DataFrame().to_csv(summary_path, index=False)
        pd.DataFrame().to_csv(overall_path, index=False)
        pd.DataFrame().to_csv(pairwise_path, index=False)
        return summary_path, overall_path, pairwise_path, plot_paths

    df = sample_burden.copy()
    df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(0)
    if "Collection Date" in df.columns:
        df["Collection Year"] = pd.to_numeric(df["Collection Date"], errors="coerce").astype("Int64").astype(str).replace("<NA>", "Unknown")

    group_cols = [col for col in ["Geographic Location", "Continent", "Subcontinent", "Collection Year"] if col in df.columns]
    summary_rows = []
    overall_rows = []
    pairwise_rows = []
    html_dir = os.path.join(figures_dir, "html_files")
    os.makedirs(html_dir, exist_ok=True)
    label = _feature_label(feature_type, mode)

    for group_col in group_cols:
        work = df[[group_col, count_col]].copy()
        work[group_col] = work[group_col].fillna("Unknown").astype(str).replace({"": "Unknown", "nan": "Unknown", "None": "Unknown"})
        for group_name, sub_df in work.groupby(group_col):
            values = sub_df[count_col]
            summary_rows.append({
                "grouping_variable": group_col,
                "group": group_name,
                "sample_count": int(len(values)),
                "samples_with_feature": int((values > 0).sum()),
                "mean_feature_count": float(values.mean()) if len(values) else 0.0,
                "median_feature_count": float(values.median()) if len(values) else 0.0,
                "min_feature_count": float(values.min()) if len(values) else 0.0,
                "max_feature_count": float(values.max()) if len(values) else 0.0,
            })

        valid_groups = [(name, sub[count_col].astype(float).tolist()) for name, sub in work.groupby(group_col) if len(sub) >= min_group_size]
        if len(valid_groups) >= 2:
            try:
                statistic, p_value = kruskal(*[values for _, values in valid_groups])
                overall_rows.append({
                    "grouping_variable": group_col,
                    "test": "Kruskal-Wallis",
                    "groups_tested": len(valid_groups),
                    "min_group_size": min_group_size,
                    "statistic": statistic,
                    "p_value": p_value,
                    "significant_0_05": bool(p_value < 0.05),
                })
            except Exception as exc:
                overall_rows.append({"grouping_variable": group_col, "test": "Kruskal-Wallis", "error": str(exc)})

            for i, (first_name, first_values) in enumerate(valid_groups):
                for second_name, second_values in valid_groups[i + 1:]:
                    try:
                        statistic, p_value = mannwhitneyu(first_values, second_values, alternative="two-sided")
                        pairwise_rows.append({
                            "grouping_variable": group_col,
                            "group_1": first_name,
                            "group_2": second_name,
                            "sample_size_1": len(first_values),
                            "sample_size_2": len(second_values),
                            "test": "Mann-Whitney U",
                            "statistic": statistic,
                            "p_value": p_value,
                            "significant_0_05": bool(p_value < 0.05),
                        })
                    except Exception as exc:
                        pairwise_rows.append({"grouping_variable": group_col, "group_1": first_name, "group_2": second_name, "error": str(exc)})

        plot_df = pd.DataFrame([row for row in summary_rows if row["grouping_variable"] == group_col])
        if not plot_df.empty:
            plot_df = plot_df.sort_values("mean_feature_count", ascending=True)
            plt.figure(figsize=_safe_figsize(len(plot_df), width=10))
            sns.barplot(data=plot_df, x="mean_feature_count", y="group", color="#4C78A8")
            plt.xlabel("Mean feature count")
            plt.ylabel(group_col)
            plt.title(f"{label} mean feature burden by {group_col}")
            plt.tight_layout()
            safe_group = group_col.lower().replace(" ", "_")
            path = os.path.join(figures_dir, f"{feature_type}_mean_burden_by_{safe_group}.{fig_format}")
            plt.savefig(path, dpi=300, bbox_inches="tight", format=fig_format)
            plt.close()
            plot_paths[f"mean_burden_by_{safe_group}_plot"] = path

            fig = px.bar(
                plot_df.sort_values("mean_feature_count", ascending=False),
                x="mean_feature_count",
                y="group",
                orientation="h",
                hover_data=["sample_count", "samples_with_feature", "median_feature_count", "min_feature_count", "max_feature_count"],
                labels={"mean_feature_count": "Mean feature count", "group": group_col},
                title=f"{label} mean feature burden by {group_col}",
            )
            fig.update_layout(yaxis={"categoryorder": "total ascending"})
            html_path = os.path.join(html_dir, f"{feature_type}_mean_burden_by_{safe_group}.html")
            fig.write_html(html_path)
            plot_paths[f"mean_burden_by_{safe_group}_html"] = html_path

    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    pd.DataFrame(overall_rows).to_csv(overall_path, index=False)
    pd.DataFrame(pairwise_rows).to_csv(pairwise_path, index=False)
    return summary_path, overall_path, pairwise_path, plot_paths

def _write_feature_html_index(feature_type, mode, figures_dir):
    html_dir = os.path.join(figures_dir, "html_files")
    if not os.path.isdir(html_dir):
        return None
    html_files = sorted(name for name in os.listdir(html_dir) if name.endswith(".html"))
    if not html_files:
        return None
    label = _feature_label(feature_type, mode)
    options = "\n".join([f'<option value="html_files/{name}">{name}</option>' for name in html_files])
    index = f'''<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{label} Feature Figures</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 0; color: #1f2933; }}
    header {{ padding: 1rem 1.25rem; border-bottom: 1px solid #d9e2ec; display: flex; gap: 1rem; align-items: center; }}
    h1 {{ font-size: 1.1rem; margin: 0; }}
    select {{ min-width: 20rem; padding: 0.4rem; }}
    iframe {{ width: 100%; height: calc(100vh - 70px); border: 0; }}
  </style>
</head>
<body>
  <header>
    <h1>{label} Feature Figures</h1>
    <select id="plotSelect" onchange="document.getElementById('plotFrame').src=this.value">
      {options}
    </select>
  </header>
  <iframe id="plotFrame" src="html_files/{html_files[0]}"></iframe>
</body>
</html>
'''
    path = os.path.join(figures_dir, "index.html")
    with open(path, "w") as handle:
        handle.write(index)
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
    analysis_dir = os.path.join(output_feature_dir, "analysis")
    figures_dir = os.path.join(output_feature_dir, "figures")
    merged_output_dir = os.path.join(output_feature_dir, "merged_output")
    os.makedirs(analysis_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(merged_output_dir, exist_ok=True)

    numeric = merged_df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    if min_identity > 0:
        numeric = numeric.where(numeric >= min_identity, 0)
        merged_df.loc[:, feature_cols] = numeric
    if "NUM_FOUND" in merged_df.columns:
        merged_df["NUM_FOUND"] = (numeric > 0).sum(axis=1).astype(int)

    merged_path = os.path.join(merged_output_dir, f"{feature_type}_merged.csv")
    merged_df.to_csv(merged_path, index=False)

    id_vars = list(merged_df.columns[:merged_df.columns.get_loc("NUM_FOUND") + 1]) if "NUM_FOUND" in merged_df.columns else list(ncbi_df.columns)
    tidy = merged_df.melt(id_vars=id_vars, value_vars=feature_cols, var_name="feature_id", value_name="identity")
    tidy["identity"] = pd.to_numeric(tidy["identity"], errors="coerce").fillna(0)
    tidy["presence"] = (tidy["identity"] > 0).astype(int)
    tidy["feature_type"] = feature_type
    tidy["feature_category"] = tidy["feature_id"].apply(lambda value: _feature_category(results_df, value, mode))
    tidy_path = os.path.join(merged_output_dir, f"{feature_type}_tidy.csv")
    tidy.to_csv(tidy_path, index=False)

    sample_col = "Assembly BioSample Accession" if "Assembly BioSample Accession" in tidy.columns else "Assembly Accession"
    total_samples = max(tidy[sample_col].nunique(), 1)
    present = tidy[tidy["presence"] == 1].copy()

    if present.empty:
        feature_summary = pd.DataFrame(columns=["feature_id", "feature_category", "present_samples", "prevalence_percentage", "mean_identity", "min_identity", "max_identity"])
        category_summary = pd.DataFrame(columns=["feature_category", "present_samples", "prevalence_percentage", "unique_features"])
    else:
        feature_summary = present.groupby(["feature_id", "feature_category"], as_index=False).agg(
            present_samples=(sample_col, "nunique"),
            mean_identity=("identity", "mean"),
            min_identity=("identity", "min"),
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

    feature_summary_path = os.path.join(analysis_dir, f"{feature_type}_feature_summary.csv")
    category_summary_path = os.path.join(analysis_dir, f"{feature_type}_category_summary.csv")
    feature_summary.to_csv(feature_summary_path, index=False)
    category_summary.to_csv(category_summary_path, index=False)

    sample_meta_cols = [col for col in ["Assembly Accession", "Assembly BioSample Accession", "Geographic Location", "Continent", "Subcontinent", "Collection Date"] if col in merged_df.columns]
    sample_burden = merged_df[sample_meta_cols].copy() if sample_meta_cols else merged_df[["Assembly Accession"]].copy()
    sample_burden[f"{feature_type}_feature_count"] = (numeric > 0).sum(axis=1).astype(int)
    sample_burden_path = os.path.join(analysis_dir, f"{feature_type}_sample_burden.csv")
    sample_burden.to_csv(sample_burden_path, index=False)

    geographic_summary_path, geographic_html_path = _write_feature_geographic_summary(feature_type, mode, analysis_dir, figures_dir, tidy, sample_burden)
    temporal_summary_path, temporal_html_path = _write_feature_temporal_summary(feature_type, mode, analysis_dir, figures_dir, sample_burden)
    plot_paths = _write_feature_plots(feature_type, mode, figures_dir, fig_format, tidy, feature_summary, category_summary, sample_burden, sample_col)
    cooc_matrix_path, top_pairs_path, cooc_plot_paths = _write_feature_cooccurrence(feature_type, mode, analysis_dir, figures_dir, fig_format, tidy, total_samples, sample_col)
    group_summary_path, group_overall_tests_path, group_pairwise_tests_path, group_plot_paths = _write_feature_group_analysis(feature_type, mode, analysis_dir, figures_dir, fig_format, sample_burden)
    html_index_path = _write_feature_html_index(feature_type, mode, figures_dir)

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
        "temporal_summary": temporal_summary_path,
        "geographic_burden_html": geographic_html_path,
        "temporal_burden_html": temporal_html_path,
        "feature_cooccurrence_matrix": cooc_matrix_path,
        "top_feature_pairs": top_pairs_path,
        "group_burden_summary": group_summary_path,
        "group_overall_tests": group_overall_tests_path,
        "group_pairwise_tests": group_pairwise_tests_path,
        "html_index": html_index_path,
        **plot_paths,
        **group_plot_paths,
        **cooc_plot_paths,
    }
