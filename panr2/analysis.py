import logging
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def get_resistance_class(series):
    """Return a stable resistance class label for grouped ARG rows."""
    clean = series.dropna().astype(str)
    clean = clean[~clean.isin(["", "0", "nan", "None"])]
    if clean.empty:
        return "Unknown"
    return clean.mode().iloc[0]

def generate_comprehensive_analysis_outputs(
    tidy_df,
    output_dir,
    base_name,
    fig_format,
    core_threshold=95.0,
    rare_threshold=5.0,
    top_n=25,
    cooccurrence_min_prevalence=0.0,
    cooccurrence_top_n=25,
):
    """Generate clean panresistome summary tables and compact plots."""
    analysis_dir = os.path.join(output_dir, "analysis")
    plot_dir = os.path.join(analysis_dir, "plots")
    os.makedirs(analysis_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    df = tidy_df.copy()
    df["Identity"] = pd.to_numeric(df["Identity"], errors="coerce").fillna(0)
    df["Presence"] = pd.to_numeric(df["Presence"], errors="coerce").fillna(0).astype(int)
    if "RESISTANCE" not in df.columns:
        df["RESISTANCE"] = "Unknown"
    df["RESISTANCE"] = df["RESISTANCE"].fillna("Unknown").replace({"": "Unknown", "0": "Unknown"})

    sample_col = "Assembly BioSample Accession" if "Assembly BioSample Accession" in df.columns else "Assembly Accession"
    sample_meta_cols = [
        col for col in [
            "Assembly Accession", "Assembly BioSample Accession", "Organism Name", "Geographic Location",
            "Continent", "Subcontinent", "Collection Date", "Host", "Isolation Source", "NUM_FOUND"
        ]
        if col in df.columns
    ]
    total_samples = max(df[sample_col].nunique(), 1)
    present_df = df[df["Presence"] == 1].copy()

    sample_meta = df[sample_meta_cols].drop_duplicates(subset=[sample_col]) if sample_meta_cols else df[[sample_col]].drop_duplicates()
    if present_df.empty:
        burden = sample_meta.copy()
        burden["unique_arg_genes"] = 0
        burden["resistance_class_count"] = 0
        burden["mean_arg_identity"] = 0.0
        burden["max_arg_identity"] = 0.0
    else:
        burden_metrics = present_df.groupby(sample_col).agg(
            unique_arg_genes=("Gene", "nunique"),
            resistance_class_count=("RESISTANCE", "nunique"),
            mean_arg_identity=("Identity", "mean"),
            max_arg_identity=("Identity", "max"),
        ).reset_index()
        burden = sample_meta.merge(burden_metrics, on=sample_col, how="left")
        for col in ["unique_arg_genes", "resistance_class_count", "mean_arg_identity", "max_arg_identity"]:
            burden[col] = burden[col].fillna(0)
        for col in ["unique_arg_genes", "resistance_class_count"]:
            burden[col] = burden[col].astype(int)
    burden_path = os.path.join(analysis_dir, f"{base_name}_sample_resistome_burden.csv")
    burden.to_csv(burden_path, index=False)

    gene_group = df.groupby("Gene")
    gene_prevalence = gene_group.agg(
        present_samples=("Presence", "sum"),
        mean_identity_all=("Identity", "mean"),
    ).reset_index()
    gene_prevalence["prevalence_percentage"] = (gene_prevalence["present_samples"] / total_samples) * 100

    if present_df.empty:
        gene_present_metrics = pd.DataFrame(columns=["Gene", "mean_identity_present", "min_identity_present", "max_identity_present", "countries_detected", "continents_detected", "resistance_class"])
    else:
        present_aggs = {
            "mean_identity_present": ("Identity", "mean"),
            "min_identity_present": ("Identity", "min"),
            "max_identity_present": ("Identity", "max"),
            "resistance_class": ("RESISTANCE", get_resistance_class),
        }
        if "Geographic Location" in present_df.columns:
            present_aggs["countries_detected"] = ("Geographic Location", "nunique")
        if "Continent" in present_df.columns:
            present_aggs["continents_detected"] = ("Continent", "nunique")
        gene_present_metrics = present_df.groupby("Gene").agg(**present_aggs).reset_index()

    gene_prevalence = gene_prevalence.merge(gene_present_metrics, on="Gene", how="left")
    gene_prevalence["resistance_class"] = gene_prevalence["resistance_class"].fillna("Not detected")
    for col in ["mean_identity_present", "min_identity_present", "max_identity_present", "countries_detected", "continents_detected"]:
        if col in gene_prevalence.columns:
            gene_prevalence[col] = gene_prevalence[col].fillna(0)
    gene_prevalence = gene_prevalence.sort_values(["prevalence_percentage", "present_samples", "Gene"], ascending=[False, False, True])
    gene_path = os.path.join(analysis_dir, f"{base_name}_gene_prevalence_summary.csv")
    gene_prevalence.to_csv(gene_path, index=False)

    def categorize(prevalence):
        if prevalence >= core_threshold:
            return "core"
        if prevalence <= rare_threshold:
            return "rare"
        return "accessory"

    category_summary = gene_prevalence.copy()
    category_summary["resistome_category"] = category_summary["prevalence_percentage"].apply(categorize)
    category_path = os.path.join(analysis_dir, f"{base_name}_resistome_category_summary.csv")
    category_summary.to_csv(category_path, index=False)

    if present_df.empty:
        class_summary = pd.DataFrame(columns=["resistance_class", "samples_with_class", "prevalence_percentage", "unique_genes", "top_genes"])
    else:
        class_rows = []
        for resistance_class, sub_df in present_df.groupby("RESISTANCE"):
            gene_counts = sub_df.groupby("Gene")[sample_col].nunique().sort_values(ascending=False)
            top_genes = ";".join(gene_counts.head(10).index.astype(str))
            samples_with_class = sub_df[sample_col].nunique()
            class_rows.append({
                "resistance_class": resistance_class,
                "samples_with_class": samples_with_class,
                "prevalence_percentage": (samples_with_class / total_samples) * 100,
                "unique_genes": sub_df["Gene"].nunique(),
                "top_genes": top_genes,
            })
        class_summary = pd.DataFrame(class_rows).sort_values(["prevalence_percentage", "unique_genes", "resistance_class"], ascending=[False, False, True])
    class_path = os.path.join(analysis_dir, f"{base_name}_resistance_class_summary.csv")
    class_summary.to_csv(class_path, index=False)

    def build_cooccurrence_outputs(entity_col, label, prevalence_col, matrix_filename, pairs_filename, plot_filename):
        if present_df.empty:
            empty_matrix = pd.DataFrame()
            empty_pairs = pd.DataFrame(columns=[
                f"{label}_1", f"{label}_2", "cooccurring_samples", "cooccurrence_percentage",
                f"{label}_1_samples", f"{label}_2_samples", "jaccard_index"
            ])
            matrix_path = os.path.join(analysis_dir, matrix_filename)
            pairs_path = os.path.join(analysis_dir, pairs_filename)
            empty_matrix.to_csv(matrix_path)
            empty_pairs.to_csv(pairs_path, index=False)
            return matrix_path, pairs_path

        presence = (
            present_df[[sample_col, entity_col]]
            .dropna()
            .drop_duplicates()
            .assign(value=1)
            .pivot_table(index=sample_col, columns=entity_col, values="value", fill_value=0, aggfunc="max")
        )
        if presence.empty:
            matrix_path = os.path.join(analysis_dir, matrix_filename)
            pairs_path = os.path.join(analysis_dir, pairs_filename)
            pd.DataFrame().to_csv(matrix_path)
            pd.DataFrame().to_csv(pairs_path, index=False)
            return matrix_path, pairs_path

        prevalence = (presence.sum(axis=0) / total_samples) * 100
        kept_entities = prevalence[prevalence >= cooccurrence_min_prevalence].index.tolist()
        presence = presence[kept_entities]
        if presence.empty:
            matrix_path = os.path.join(analysis_dir, matrix_filename)
            pairs_path = os.path.join(analysis_dir, pairs_filename)
            pd.DataFrame().to_csv(matrix_path)
            pd.DataFrame(columns=[
                f"{label}_1", f"{label}_2", "cooccurring_samples", "cooccurrence_percentage",
                f"{label}_1_samples", f"{label}_2_samples", "jaccard_index"
            ]).to_csv(pairs_path, index=False)
            return matrix_path, pairs_path

        cooc_matrix = presence.T.dot(presence).astype(int)
        matrix_path = os.path.join(analysis_dir, matrix_filename)
        cooc_matrix.to_csv(matrix_path)

        pair_rows = []
        entities = list(cooc_matrix.columns)
        sample_counts = presence.sum(axis=0).astype(int).to_dict()
        for i, first in enumerate(entities):
            for second in entities[i + 1:]:
                co_count = int(cooc_matrix.loc[first, second])
                union_count = int(((presence[first] + presence[second]) > 0).sum())
                pair_rows.append({
                    f"{label}_1": first,
                    f"{label}_2": second,
                    "cooccurring_samples": co_count,
                    "cooccurrence_percentage": (co_count / total_samples) * 100,
                    f"{label}_1_samples": sample_counts[first],
                    f"{label}_2_samples": sample_counts[second],
                    "jaccard_index": co_count / union_count if union_count else 0,
                })
        pair_df = pd.DataFrame(pair_rows)
        if not pair_df.empty:
            pair_df = pair_df[pair_df["cooccurring_samples"] > 0]
            pair_df = pair_df.sort_values(["cooccurring_samples", "jaccard_index", f"{label}_1", f"{label}_2"], ascending=[False, False, True, True]).head(cooccurrence_top_n)
        pairs_path = os.path.join(analysis_dir, pairs_filename)
        pair_df.to_csv(pairs_path, index=False)

        plot_entities = prevalence.loc[kept_entities].sort_values(ascending=False).head(cooccurrence_top_n).index.tolist()
        if len(plot_entities) >= 2:
            plot_matrix = cooc_matrix.loc[plot_entities, plot_entities]
            size = max(5, min(14, 0.45 * len(plot_entities) + 3))
            plt.figure(figsize=(size, size))
            sns.heatmap(plot_matrix, cmap="Blues", square=True, linewidths=0.3, linecolor="white", cbar_kws={"label": "Co-occurring samples"})
            plt.title(f"{label.replace('_', ' ').title()} co-occurrence")
            plt.xlabel(label.replace("_", " ").title())
            plt.ylabel(label.replace("_", " ").title())
            plt.xticks(rotation=45, ha="right")
            plt.yticks(rotation=0)
            plt.tight_layout()
            plot_path = os.path.join(plot_dir, plot_filename)
            plt.savefig(plot_path, dpi=300, bbox_inches="tight", format=fig_format)
            plt.close()

        return matrix_path, pairs_path

    gene_cooc_matrix_path, gene_cooc_pairs_path = build_cooccurrence_outputs(
        "Gene",
        "gene",
        "prevalence_percentage",
        f"{base_name}_gene_cooccurrence_matrix.csv",
        f"{base_name}_top_gene_pairs.csv",
        f"{base_name}_gene_cooccurrence_heatmap.{fig_format}",
    )
    class_cooc_matrix_path, class_cooc_pairs_path = build_cooccurrence_outputs(
        "RESISTANCE",
        "class",
        "prevalence_percentage",
        f"{base_name}_class_cooccurrence_matrix.csv",
        f"{base_name}_top_class_pairs.csv",
        f"{base_name}_class_cooccurrence_heatmap.{fig_format}",
    )

    # Clean, compact plots for the summary layer.
    if not gene_prevalence.empty:
        top_genes_df = gene_prevalence.head(top_n).sort_values("prevalence_percentage", ascending=True)
        height = max(4, min(14, 0.35 * len(top_genes_df) + 1.5))
        plt.figure(figsize=(10, height))
        sns.barplot(data=top_genes_df, x="prevalence_percentage", y="Gene", hue="resistance_class", dodge=False, palette="tab20")
        plt.xlabel("Prevalence (%)")
        plt.ylabel("ARG")
        plt.title(f"Top {len(top_genes_df)} ARGs by prevalence")
        plt.xlim(0, 100)
        plt.legend(title="Resistance class", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
        plt.tight_layout()
        plot_path = os.path.join(plot_dir, f"{base_name}_top_gene_prevalence.{fig_format}")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()

    if not class_summary.empty:
        class_plot_df = class_summary.head(top_n).sort_values("prevalence_percentage", ascending=True)
        height = max(4, min(12, 0.4 * len(class_plot_df) + 1.5))
        plt.figure(figsize=(9, height))
        sns.barplot(data=class_plot_df, x="prevalence_percentage", y="resistance_class", color="#4C78A8")
        plt.xlabel("Samples carrying class (%)")
        plt.ylabel("Resistance class")
        plt.title("Resistance class prevalence")
        plt.xlim(0, 100)
        plt.tight_layout()
        plot_path = os.path.join(plot_dir, f"{base_name}_resistance_class_prevalence.{fig_format}")
        plt.savefig(plot_path, dpi=300, bbox_inches="tight", format=fig_format)
        plt.close()

    logging.info(f"Comprehensive analysis outputs saved to {analysis_dir}")
    return {
        "sample_resistome_burden": burden_path,
        "gene_prevalence_summary": gene_path,
        "resistome_category_summary": category_path,
        "resistance_class_summary": class_path,
        "gene_cooccurrence_matrix": gene_cooc_matrix_path,
        "top_gene_pairs": gene_cooc_pairs_path,
        "class_cooccurrence_matrix": class_cooc_matrix_path,
        "top_class_pairs": class_cooc_pairs_path,
    }

