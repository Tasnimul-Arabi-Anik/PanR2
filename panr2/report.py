import html
import logging
import os
from datetime import datetime

import pandas as pd

from panr2.dashboard import write_dashboard


def _read_csv(path):
    if os.path.exists(path):
        try:
            return pd.read_csv(path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    return pd.DataFrame()


def _fmt_pct(value):
    try:
        return f"{float(value):.1f}%"
    except Exception:
        return "not available"


def _fmt_num(value):
    try:
        if float(value).is_integer():
            return str(int(value))
        return f"{float(value):.2f}"
    except Exception:
        return "not available"


def _sentence_list(items):
    items = [str(item) for item in items if str(item)]
    if not items:
        return "none"
    if len(items) == 1:
        return items[0]
    if len(items) == 2:
        return f"{items[0]} and {items[1]}"
    return ", ".join(items[:-1]) + f", and {items[-1]}"


def _markdown_table(df, columns=None, max_rows=10):
    if df.empty:
        return "No records available.\n"
    table = df.copy()
    if columns:
        table = table[[col for col in columns if col in table.columns]]
    table = table.head(max_rows).fillna("")
    if table.empty:
        return "No records available.\n"
    headers = list(table.columns)
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for _, row in table.iterrows():
        values = []
        for col in headers:
            value = row[col]
            if str(value).lower() == "nan":
                value = ""
            values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines) + "\n"


def _status_counts(qc_df):
    if qc_df.empty or "status" not in qc_df.columns:
        return {}
    return qc_df["status"].value_counts().to_dict()


def _top_names(df, name_col, value_col, n=5):
    if df.empty or name_col not in df.columns or value_col not in df.columns:
        return []
    ranked = df.sort_values(value_col, ascending=False).head(n)
    return [f"{row[name_col]} ({_fmt_pct(row[value_col])})" for _, row in ranked.iterrows()]


def _feature_label(feature_type):
    labels = {
        "vfdb": "VFDB",
        "virulence": "VFDB",
        "plasmidfinder": "PlasmidFinder",
        "plasmid": "PlasmidFinder",
        "mobileelementfinder": "MobileElementFinder",
        "isfinder": "ISfinder",
        "integronfinder": "IntegronFinder",
        "iceberg": "ICEberg",
        "mlst": "MLST",
        "defensefinder": "DefenseFinder",
        "prophage": "Prophage",
    }
    return labels.get(feature_type.lower(), feature_type.replace("_", " ").title())


def write_report(output_dir, base_name, ncbi_output_dir=None, options=None, panr2_version="unknown", input_files=None, feature_outputs=None, cross_database_outputs=None, citation_outputs=None, temporal_outputs=None):
    """Write a deterministic journal-style PanR2 report."""
    options = options or {}
    input_files = input_files or {}
    feature_outputs = feature_outputs or {}
    cross_database_outputs = cross_database_outputs or {}
    citation_outputs = citation_outputs or {}
    temporal_outputs = temporal_outputs or {}
    report_dir = os.path.join(output_dir, "report")
    os.makedirs(report_dir, exist_ok=True)
    ncbi_output_dir = ncbi_output_dir or output_dir

    qc_dir = os.path.join(output_dir, "qc")
    analysis_dir = os.path.join(ncbi_output_dir, "analysis")
    figures_dir = os.path.join(ncbi_output_dir, "figures")
    stat_dir = os.path.join(figures_dir, "Stat_analysis")

    qc_df = _read_csv(os.path.join(qc_dir, "panr2_input_qc.csv"))
    tool_manifest = _read_csv(os.path.join(qc_dir, "panr2_tool_manifest.csv"))
    filter_df = _read_csv(os.path.join(qc_dir, f"{base_name}_filter_report.csv"))
    burden = _read_csv(os.path.join(analysis_dir, f"{base_name}_sample_resistome_burden.csv"))
    genes = _read_csv(os.path.join(analysis_dir, f"{base_name}_gene_prevalence_summary.csv"))
    categories = _read_csv(os.path.join(analysis_dir, f"{base_name}_resistome_category_summary.csv"))
    classes = _read_csv(os.path.join(analysis_dir, f"{base_name}_resistance_class_summary.csv"))
    gene_pairs = _read_csv(os.path.join(analysis_dir, f"{base_name}_top_gene_pairs.csv"))
    class_pairs = _read_csv(os.path.join(analysis_dir, f"{base_name}_top_class_pairs.csv"))
    correlation = _read_csv(os.path.join(stat_dir, "combined_geographic_correlation_summary.csv"))

    total_samples = int(len(burden)) if not burden.empty else 0
    total_genes = int(len(genes)) if not genes.empty else 0
    total_classes = int(len(classes)) if not classes.empty else 0
    mean_burden = burden["unique_arg_genes"].mean() if "unique_arg_genes" in burden.columns and not burden.empty else 0
    median_burden = burden["unique_arg_genes"].median() if "unique_arg_genes" in burden.columns and not burden.empty else 0
    max_burden = burden["unique_arg_genes"].max() if "unique_arg_genes" in burden.columns and not burden.empty else 0

    qc_counts = _status_counts(qc_df)
    filter_removed = 0
    if not filter_df.empty and "removed" in filter_df.columns:
        filter_removed = int(pd.to_numeric(filter_df["removed"], errors="coerce").fillna(0).sum())

    top_gene_text = _sentence_list(_top_names(genes, "Gene", "prevalence_percentage", 5))
    top_class_text = _sentence_list(_top_names(classes, "resistance_class", "prevalence_percentage", 5))

    category_counts = {}
    if not categories.empty and "resistome_category" in categories.columns:
        category_counts = categories["resistome_category"].value_counts().to_dict()
    core_count = int(category_counts.get("core", 0))
    accessory_count = int(category_counts.get("accessory", 0))
    rare_count = int(category_counts.get("rare", 0))

    lines = []
    lines.append("# PanR2 Panresistome Analysis Report")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"**Dataset prefix:** `{base_name}`")
    lines.append(f"**PanR2 version:** `{panr2_version}`")
    lines.append("")

    lines.append("## Executive Summary")
    lines.append("")
    lines.append(
        f"A total of {total_samples} assemblies or samples were included in the panresistome analysis. "
        f"The analysis identified {total_genes} antimicrobial resistance gene(s) distributed across {total_classes} resistance class(es). "
        f"The mean ARG burden was {_fmt_num(mean_burden)} gene(s) per sample, with a median of {_fmt_num(median_burden)} and a maximum of {_fmt_num(max_burden)}."
    )
    lines.append("")
    lines.append(
        f"The most prevalent ARGs were {top_gene_text}. The dominant resistance classes were {top_class_text}. "
        f"Using the configured thresholds, the resistome contained {core_count} core, {accessory_count} accessory, and {rare_count} rare ARG(s)."
    )
    lines.append("")

    lines.append("## Input Data and Quality Control")
    lines.append("")
    lines.append(f"NCBI metadata input: `{input_files.get('ncbi_clean', 'not recorded')}`")
    lines.append(f"ABRicate summary input: `{input_files.get('abricate_summary', 'not recorded')}`")
    lines.append(f"ABRicate results input: `{input_files.get('abricate_results', 'not recorded')}`")
    lines.append("")
    lines.append(
        f"Input validation produced {qc_counts.get('FAIL', 0)} failure(s), {qc_counts.get('WARN', 0)} warning(s), "
        f"{qc_counts.get('PASS', 0)} passed check(s), and {qc_counts.get('INFO', 0)} informational record(s)."
    )
    lines.append("")
    lines.append(_markdown_table(qc_df, ["check", "status", "detail", "value"], max_rows=25))

    if not tool_manifest.empty:
        lines.append("## Integrated Tool Runs")
        lines.append("")
        lines.append(
            "PanR2 recorded integrated upstream annotation runs for reproducibility. "
            "Database sequence counts and dates are reported when the upstream tool exposes them."
        )
        lines.append("")
        lines.append(_markdown_table(tool_manifest, ["tool", "version", "database", "database_sequences", "database_date", "results", "summary", "status"], max_rows=30))

    lines.append("## Filtering and Analysis Parameters")
    lines.append("")
    lines.append("The analysis used the following command options and thresholds:")
    lines.append("")
    for key in sorted(options):
        lines.append(f"- `{key}`: `{options[key]}`")
    lines.append("")
    lines.append(f"Filtering removed {filter_removed} row(s) or ARG call(s) across recorded filters.")
    lines.append("")
    lines.append(_markdown_table(filter_df, ["filter", "enabled", "before", "after", "removed", "detail"], max_rows=20))

    lines.append("## Resistome Burden")
    lines.append("")
    lines.append(
        f"Per-sample ARG burden ranged from {_fmt_num(burden['unique_arg_genes'].min() if 'unique_arg_genes' in burden.columns and not burden.empty else 0)} "
        f"to {_fmt_num(max_burden)} unique ARG(s). Resistance class richness ranged from "
        f"{_fmt_num(burden['resistance_class_count'].min() if 'resistance_class_count' in burden.columns and not burden.empty else 0)} to "
        f"{_fmt_num(burden['resistance_class_count'].max() if 'resistance_class_count' in burden.columns and not burden.empty else 0)} class(es) per sample."
    )
    lines.append("")
    lines.append(_markdown_table(burden, ["Assembly Accession", "Assembly BioSample Accession", "Geographic Location", "Continent", "Collection Date", "unique_arg_genes", "resistance_class_count", "mean_arg_identity"], max_rows=15))

    lines.append("## Gene Prevalence")
    lines.append("")
    lines.append(
        "Gene-level prevalence was calculated as the percentage of samples in which each ARG was detected after filtering. "
        "Identity summaries describe detected calls only where available."
    )
    lines.append("")
    lines.append(_markdown_table(genes, ["Gene", "present_samples", "prevalence_percentage", "resistance_class", "mean_identity_present", "min_identity_present", "max_identity_present", "countries_detected", "continents_detected"], max_rows=20))

    lines.append("## Core, Accessory, and Rare Resistome")
    lines.append("")
    lines.append(
        f"ARGs were classified using the configured prevalence thresholds: core >= {options.get('core_threshold', 'not recorded')}%, "
        f"rare <= {options.get('rare_threshold', 'not recorded')}%, and accessory between those thresholds."
    )
    lines.append("")
    lines.append(_markdown_table(categories, ["Gene", "prevalence_percentage", "resistance_class", "resistome_category"], max_rows=25))

    lines.append("## Resistance Class Composition")
    lines.append("")
    lines.append(
        "Resistance class prevalence was calculated from samples carrying at least one ARG assigned to each class. "
        "The top gene list records the highest-support genes within each class."
    )
    lines.append("")
    lines.append(_markdown_table(classes, ["resistance_class", "samples_with_class", "prevalence_percentage", "unique_genes", "top_genes"], max_rows=20))

    lines.append("## Co-occurrence Analysis")
    lines.append("")
    if gene_pairs.empty and class_pairs.empty:
        lines.append("No positive ARG or resistance-class co-occurrence pairs were detected after filtering and prevalence thresholding.")
    else:
        lines.append(
            "Co-occurrence was calculated from shared sample presence. Pair tables report shared sample count, "
            "co-occurrence percentage, marginal sample counts, and Jaccard index."
        )
    lines.append("")
    lines.append("### Top ARG Pairs")
    lines.append(_markdown_table(gene_pairs, ["gene_1", "gene_2", "cooccurring_samples", "cooccurrence_percentage", "gene_1_samples", "gene_2_samples", "jaccard_index"], max_rows=20))
    lines.append("### Top Resistance Class Pairs")
    lines.append(_markdown_table(class_pairs, ["class_1", "class_2", "cooccurring_samples", "cooccurrence_percentage", "class_1_samples", "class_2_samples", "jaccard_index"], max_rows=20))

    if feature_outputs:
        lines.append("## Optional Database Feature Analysis")
        lines.append("")
        lines.append(
            "Optional non-AMR feature analyses were kept separate from resistance-class analysis because virulence factors, plasmid replicons, and mobile genetic element features do not share the same biological class structure as antimicrobial resistance genes. "
            "Prevalence values below describe feature presence after the same identity filtering threshold used for this run."
        )
        lines.append("")
        feature_order = ["vfdb", "plasmidfinder", "mobileelementfinder", "isfinder", "integronfinder", "iceberg", "mlst", "defensefinder", "prophage", "virulence", "plasmid"]
        for feature_type in feature_order:
            if feature_type not in feature_outputs:
                continue
            outputs = feature_outputs.get(feature_type, {})
            feature_summary = _read_csv(outputs.get("feature_summary", ""))
            category_summary = _read_csv(outputs.get("category_summary", ""))
            sample_burden = _read_csv(outputs.get("sample_burden", ""))
            feature_qc = _read_csv(outputs.get("qc_summary", ""))
            geographic_summary = _read_csv(outputs.get("geographic_summary", ""))
            temporal_summary = _read_csv(outputs.get("temporal_summary", ""))
            top_pairs = _read_csv(outputs.get("top_feature_pairs", ""))
            group_burden = _read_csv(outputs.get("group_burden_summary", ""))
            group_overall = _read_csv(outputs.get("group_overall_tests", ""))
            feature_count = len(feature_summary) if not feature_summary.empty else 0
            carrying_samples = 0
            count_col = f"{feature_type}_feature_count"
            if count_col in sample_burden.columns:
                carrying_samples = int((pd.to_numeric(sample_burden[count_col], errors="coerce").fillna(0) > 0).sum())
            label = _feature_label(feature_type)
            lines.append(f"### {label} Features")
            lines.append("")
            lines.append(
                f"The {label} module detected {feature_count} feature(s), with {carrying_samples} sample(s) carrying at least one {label} feature."
            )
            if feature_type == "mlst":
                lines.append("")
                lines.append(
                    "MLST sequence types are treated as typing features for metadata-linked comparisons. They are not interpreted as AMR, virulence, plasmid, or mobile-element annotations."
                )
            if feature_type == "defensefinder":
                lines.append("")
                lines.append(
                    "DefenseFinder features describe predicted genome-defense systems and are interpreted separately from AMR and mobile-element features."
                )
            if feature_type == "prophage":
                lines.append("")
                lines.append(
                    "Prophage features describe user-provided viral-region or prophage annotations. They indicate sample-level viral-region burden and do not by themselves prove active phage, transfer, or AMR linkage."
                )
            if feature_type == "iceberg":
                lines.append("")
                lines.append(
                    "The ICEberg module is treated as ICEberg-style feature analysis. PanR2 converts user-provided ICE/IME/CIME annotation tables or ABRicate-style ICEberg inputs into database-specific summaries; it does not run an ICEberg annotation program directly."
                )
            lines.append("")
            lines.append(f"### {label} QC Summary")
            lines.append(_markdown_table(feature_qc, ["metric", "value", "detail"], max_rows=25))
            lines.append(f"### {label} Top Features")
            lines.append(_markdown_table(feature_summary, ["feature_id", "feature_category", "present_samples", "prevalence_percentage", "mean_identity", "min_identity", "max_identity"], max_rows=20))
            lines.append(f"### {label} Categories")
            lines.append(_markdown_table(category_summary, ["feature_category", "present_samples", "prevalence_percentage", "unique_features"], max_rows=20))
            lines.append(f"### {label} Geographic Summary")
            lines.append(
                "Feature-burden summaries are descriptive and are not treated as antibiotic resistance classes. "
                "They report sample counts, samples carrying at least one feature, mean feature count, and leading detected features by metadata group."
            )
            lines.append("")
            lines.append(_markdown_table(geographic_summary, ["geographic_level", "region", "sample_count", "samples_with_feature", "mean_feature_count", "median_feature_count", "max_feature_count", "top_features"], max_rows=20))
            lines.append(f"### {label} Temporal Summary")
            lines.append(_markdown_table(temporal_summary, ["collection_year", "sample_count", "samples_with_feature", "mean_feature_count", "median_feature_count", "max_feature_count"], max_rows=20))
            lines.append(f"### {label} Feature Co-occurrence")
            if top_pairs.empty:
                lines.append("No positive feature co-occurrence pairs were detected after identity filtering.")
            else:
                lines.append("Feature co-occurrence was calculated from shared sample presence and is reported descriptively; it does not imply physical linkage or causal association.")
            lines.append("")
            lines.append(_markdown_table(top_pairs, ["feature_1", "feature_2", "cooccurring_samples", "cooccurrence_percentage", "feature_1_samples", "feature_2_samples", "jaccard_index"], max_rows=20))
            lines.append(f"### {label} Group Burden Comparisons")
            lines.append(
                "Mean feature burden was summarized across available metadata groups. Nonparametric tests are reported only when at least two groups meet the minimum sample-size requirement."
            )
            lines.append("")
            lines.append(_markdown_table(group_burden, ["grouping_variable", "group", "sample_count", "samples_with_feature", "mean_feature_count", "median_feature_count", "min_feature_count", "max_feature_count"], max_rows=25))
            lines.append(f"### {label} Group-Level Statistical Tests")
            lines.append(_markdown_table(group_overall, ["grouping_variable", "test", "groups_tested", "min_group_size", "statistic", "p_value", "significant_0_05", "error"], max_rows=20))

    if cross_database_outputs:
        associations = _read_csv(cross_database_outputs.get("cross_database_top_associations", ""))
        prevalence = _read_csv(cross_database_outputs.get("cross_database_feature_prevalence", ""))
        burden = _read_csv(cross_database_outputs.get("sample_integrated_feature_burden", ""))
        enrichment = _read_csv(cross_database_outputs.get("feature_enrichment_by_metadata", ""))
        amr_mge = _read_csv(cross_database_outputs.get("amr_mge_associations", ""))
        amr_plasmid = _read_csv(cross_database_outputs.get("amr_plasmid_associations", ""))
        amr_virulence = _read_csv(cross_database_outputs.get("amr_virulence_associations", ""))

        lines.append("## Cross-Database Comparative Genomics")
        lines.append("")
        lines.append(
            "PanR2 built a unified sample-by-feature matrix across AMR, virulence, plasmid, and mobile genetic element feature families. "
            "Pairwise associations report co-occurring sample count, Jaccard index, phi coefficient, Fisher exact-test odds ratio, nominal p-value, and Benjamini-Hochberg FDR-adjusted q-value. "
            "These associations describe shared sample/genome presence only and should not be interpreted as physical linkage, plasmid localization, horizontal transfer, shared regulation, clinical phenotype, or causality without genomic-context validation."
        )
        lines.append("")
        if not prevalence.empty:
            database_counts = prevalence.groupby("database")["feature"].nunique().to_dict()
            lines.append(
                "The integrated matrix included "
                + _sentence_list([f"{count} {database} feature(s)" for database, count in sorted(database_counts.items())])
                + "."
            )
            lines.append("")
        lines.append("### Top Cross-Database Associations")
        cross_only = associations[associations["association_scope"].eq("cross_database")] if "association_scope" in associations.columns else associations
        lines.append(_markdown_table(cross_only, ["feature_1", "feature_2", "database_1", "database_2", "cooccurring_samples", "jaccard_index", "phi_coefficient", "odds_ratio", "p_value", "q_value"], max_rows=25))
        lines.append("### AMR-Mobileome Associations")
        lines.append(_markdown_table(amr_mge, ["feature_1", "feature_2", "database_1", "database_2", "cooccurring_samples", "jaccard_index", "phi_coefficient", "odds_ratio", "q_value"], max_rows=15))
        lines.append("### AMR-Plasmid Associations")
        lines.append(_markdown_table(amr_plasmid, ["feature_1", "feature_2", "database_1", "database_2", "cooccurring_samples", "jaccard_index", "phi_coefficient", "odds_ratio", "q_value"], max_rows=15))
        lines.append("### AMR-Virulence Associations")
        lines.append(_markdown_table(amr_virulence, ["feature_1", "feature_2", "database_1", "database_2", "cooccurring_samples", "jaccard_index", "phi_coefficient", "odds_ratio", "q_value"], max_rows=15))
        lines.append("### Integrated Sample Burden")
        lines.append(
            "Genome-level composite burden metrics summarize AMR, virulence, plasmid, and mobileome feature counts per sample. "
            "The mobility-associated AMR score is a screening metric that highlights AMR-bearing samples with any detected mobileome feature."
        )
        lines.append("")
        lines.append(_markdown_table(burden, ["Assembly Accession", "Geographic Location", "Continent", "amr_feature_count", "vfdb_feature_count", "plasmid_feature_count", "total_mobileome_count", "total_feature_count", "mobility_associated_amr_score"], max_rows=15))
        lines.append("### Metadata Enrichment")
        lines.append(
            "Feature-level enrichment by metadata group used Fisher exact tests with FDR correction where metadata groups met the configured minimum sample-size threshold."
        )
        lines.append("")
        lines.append(_markdown_table(enrichment, ["feature", "database", "metadata_variable", "group", "feature_present_in_group", "feature_present_outside_group", "odds_ratio", "p_value", "q_value", "direction"], max_rows=25))

    if temporal_outputs:
        feature_trends = _read_csv(temporal_outputs.get("temporal_feature_trends", ""))
        burden_trends = _read_csv(temporal_outputs.get("temporal_burden_trends", ""))
        lines.append("## Advanced Temporal Trends")
        lines.append("")
        lines.append(
            "Temporal trend summaries use available collection years only. Feature trends include Mann-Kendall trend statistics, Spearman correlation across yearly prevalence, and a logistic feature-presence slope when enough variation is available, with FDR-adjusted q-values for multiple testing. Burden trends use per-sample linear regression and Mann-Kendall summaries. These temporal tests are screening-level summaries; missing years, uneven sampling, repeated tied prevalence values, and changing surveillance intensity can strongly affect interpretation."
        )
        lines.append("")
        lines.append("### Feature-Level Temporal Trends")
        lines.append(_markdown_table(feature_trends, ["feature", "database", "years_observed", "first_year", "last_year", "mann_kendall_trend", "mann_kendall_p_value", "mann_kendall_q_value", "spearman_r", "spearman_q_value", "logistic_odds_ratio_per_year", "logistic_q_value"], max_rows=20))
        lines.append("### Burden-Level Temporal Trends")
        lines.append(_markdown_table(burden_trends, ["burden_metric", "years_observed", "linear_slope_per_year", "linear_r_value", "linear_p_value", "linear_burden_q_value", "mann_kendall_trend", "mann_kendall_q_value"], max_rows=20))

    lines.append("## Geographic and Temporal Patterns")
    lines.append("")
    if correlation.empty:
        lines.append("Correlation summaries were not available, most commonly because one or more geographic groupings did not meet the minimum sample threshold.")
    else:
        lines.append(
            "Correlation summaries describe associations between collection year and ARG burden (`NUM_FOUND`) within geographic groupings. "
            "Interpretation should consider sample size, missing dates, and constant-input warnings."
        )
        lines.append("")
        lines.append(_markdown_table(correlation, ["Geographic_Level", "Geographic_Region", "n_samples", "pearson_r", "pearson_p", "spearman_r", "spearman_p"], max_rows=25))
    lines.append("")

    lines.append("## Output Files")
    lines.append("")
    lines.append("Primary machine-readable outputs generated by this run include:")
    lines.append("")
    primary_outputs = [
        f"qc/panr2_input_qc.csv",
        f"qc/{base_name}_filter_report.csv",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'merged_output', f'ncbi_{base_name}_summary.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'merged_output', f'ncbi_{base_name}_tidy_summary.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_sample_resistome_burden.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_gene_prevalence_summary.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_resistome_category_summary.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_resistance_class_summary.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_gene_cooccurrence_matrix.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_top_gene_pairs.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_class_cooccurrence_matrix.csv'), output_dir)}",
        f"{os.path.relpath(os.path.join(ncbi_output_dir, 'analysis', f'{base_name}_top_class_pairs.csv'), output_dir)}",
    ]
    for manifest_rel_path in ["qc/panr2_tool_manifest.csv", "qc/panr2_tool_manifest.json"]:
        if os.path.exists(os.path.join(output_dir, manifest_rel_path)):
            primary_outputs.insert(2, manifest_rel_path)
    for rel_path in primary_outputs:
        lines.append(f"- `{rel_path}`")
    if cross_database_outputs:
        lines.append("")
        lines.append("Cross-database comparative genomics outputs include:")
        lines.append("")
        for key in [
            "cross_database_feature_matrix",
            "cross_database_feature_prevalence",
            "cross_database_cooccurrence_matrix",
            "cross_database_jaccard_matrix",
            "cross_database_phi_correlation_matrix",
            "cross_database_top_associations",
            "amr_mge_associations",
            "amr_plasmid_associations",
            "amr_integron_associations",
            "amr_virulence_associations",
            "amr_defense_associations",
            "amr_prophage_associations",
            "plasmid_mge_associations",
            "defense_mge_associations",
            "prophage_mge_associations",
            "sample_integrated_feature_burden",
            "feature_enrichment_by_metadata",
            "network_edges",
            "network_nodes",
            "plot_readability_warnings",
            "figure_manifest",
        ]:
            output_path = cross_database_outputs.get(key)
            if output_path and isinstance(output_path, str) and os.path.exists(output_path):
                lines.append(f"- `{os.path.relpath(output_path, output_dir)}`")
    lines.append("")

    lines.append("## Interpretive Notes and Limitations")
    lines.append("")
    lines.append(
        "PanR2 reports descriptive panresistome patterns from the provided metadata and ABRicate calls. "
        "The report should not be interpreted as evidence of clinical resistance phenotype without appropriate phenotypic validation, database curation, and epidemiological context. "
        "Small group sizes, incomplete collection dates, uneven geographic sampling, and constant ARG burdens can limit correlation and temporal analyses. "
        "Co-occurrence indicates shared presence in samples/genomes and does not imply physical linkage, plasmid localization, horizontal transfer, transfer direction, or causal association."
    )
    lines.append("")

    if feature_outputs:
        lines.append("Optional feature output files generated by this run include:")
        lines.append("")
        for feature_type in sorted(feature_outputs):
            for key in [
                "merged",
                "tidy",
                "feature_summary",
                "category_summary",
                "sample_burden",
                "qc_summary",
                "unmatched_samples",
                "geographic_summary",
                "temporal_summary",
                "feature_cooccurrence_matrix",
                "top_feature_pairs",
                "group_burden_summary",
                "group_overall_tests",
                "group_pairwise_tests",
                "feature_prevalence_plot",
                "category_prevalence_plot",
                "burden_by_continent_plot",
                "presence_heatmap",
                "identity_distribution_plot",
                "feature_cooccurrence_heatmap",
                "feature_prevalence_html",
                "category_prevalence_html",
                "burden_by_continent_html",
                "presence_heatmap_html",
                "identity_distribution_html",
                "feature_cooccurrence_heatmap_html",
                "geographic_burden_html",
                "temporal_burden_html",
                "mean_burden_by_geographic_location_plot",
                "mean_burden_by_continent_plot",
                "mean_burden_by_subcontinent_plot",
                "mean_burden_by_collection_year_plot",
                "mean_burden_by_geographic_location_html",
                "mean_burden_by_continent_html",
                "mean_burden_by_subcontinent_html",
                "mean_burden_by_collection_year_html",
                "html_index",
            ]:
                output_path = feature_outputs.get(feature_type, {}).get(key)
                if output_path and os.path.exists(output_path):
                    lines.append(f"- `{os.path.relpath(output_path, output_dir)}`")
        lines.append("")

    lines.append("## Methods Summary")
    lines.append("")
    lines.append(
        "PanR2 merged NCBI-derived sample metadata with ABRicate antimicrobial resistance gene summary output by assembly accession. "
        "ABRicate gene identity values were converted into a tidy presence/absence table after optional identity filtering. "
        "Per-sample burden, per-gene prevalence, resistance class summaries, core/accessory/rare categories, and co-occurrence statistics were then calculated from the filtered tidy table. "
        "Static and interactive visualizations were generated from the same filtered data tables. Optional VFDB, PlasmidFinder, MobileElementFinder, ISfinder, IntegronFinder, ICEberg-style, MLST, DefenseFinder, and prophage-style analyses, when provided, were written to database-named output folders and parsed as separate feature families and summarized independently from AMR resistance classes using feature prevalence, category or replicon prevalence, sample burden, database-specific QC summaries, geographic and temporal feature-burden summaries, feature presence heatmaps, identity-distribution plots where applicable, descriptive feature co-occurrence tables, and database-specific interactive HTML figures, group-burden summaries, and nonparametric group-comparison tests where sample sizes permit. ICEberg-style analysis uses user-provided ICE/IME/CIME annotations or ABRicate-style inputs converted to PanR2 feature tables; PanR2 does not run an ICEberg annotation program directly."
    )
    lines.append("")
    if cross_database_outputs:
        lines.append(
            "For cross-database comparative genomics, PanR2 joined AMR, virulence, plasmid, insertion-sequence, integron, ICE/IME/CIME, and mobile genetic element feature calls into a unified binary sample-feature matrix. "
            "Pairwise binary associations were summarized using co-occurrence counts, Jaccard index, phi coefficient, Fisher exact tests, odds ratios, and Benjamini-Hochberg FDR-adjusted q-values. "
            "Metadata enrichment used feature presence/absence against metadata-defined groups with the configured minimum group-size threshold."
        )
        lines.append("")

    if citation_outputs:
        lines.append("## Citations and Software Versions")
        lines.append("")
        lines.append("PanR2 wrote run-specific citation and software-version files for reproducibility and manuscript preparation:")
        lines.append("")
        for key in ["citations_md", "citations_bib", "software_versions"]:
            output_path = citation_outputs.get(key)
            if output_path and os.path.exists(output_path):
                lines.append(f"- `{os.path.relpath(output_path, output_dir)}`")
        lines.append("")

    lines.append("## Reproducibility")
    lines.append("")
    lines.append(f"PanR2 version: `{panr2_version}`")
    lines.append("")
    lines.append("Recorded options:")
    lines.append("")
    for key in sorted(options):
        lines.append(f"- `{key}`: `{options[key]}`")
    lines.append("")
    lines.append(
        "This report is generated deterministically from PanR2 output tables. Values should be interpreted in the context of input metadata completeness, database version, ABRicate settings, and sample size."
    )

    markdown = "\n".join(lines).rstrip() + "\n"
    md_path = os.path.join(report_dir, f"{base_name}_panr2_report.md")
    html_path = os.path.join(report_dir, f"{base_name}_panr2_report.html")
    methods_path = os.path.join(report_dir, f"{base_name}_methods.txt")

    with open(md_path, "w") as handle:
        handle.write(markdown)

    html_body = html.escape(markdown)
    html_doc = f"""<!doctype html>
<html>
<head>
  <meta charset=\"utf-8\">
  <title>PanR2 Report - {html.escape(base_name)}</title>
  <style>
    body {{ font-family: Arial, sans-serif; line-height: 1.5; max-width: 1100px; margin: 2rem auto; padding: 0 1rem; color: #1f2933; }}
    pre {{ white-space: pre-wrap; background: #f7f9fb; border: 1px solid #d9e2ec; padding: 1rem; overflow-x: auto; }}
  </style>
</head>
<body>
<pre>{html_body}</pre>
</body>
</html>
"""
    with open(html_path, "w") as handle:
        handle.write(html_doc)

    with open(methods_path, "w") as handle:
        handle.write(
            "PanR2 methods summary\n\n"
            "PanR2 merged NCBI metadata with ABRicate antimicrobial resistance gene outputs by assembly accession. "
            "ARG calls were filtered according to user-specified identity and sample inclusion thresholds, converted to tidy presence/absence format, and summarized by sample, gene, resistance class, resistome prevalence category, and co-occurrence. "
            "Core, accessory, and rare ARG classifications used the prevalence thresholds recorded in the report. "
            "Geographic and temporal summaries were generated when sample-size requirements were met.\n"
        )

    logging.info(f"PanR2 report saved to {md_path}")
    dashboard_outputs = write_dashboard(
        output_dir,
        base_name,
        feature_outputs=feature_outputs,
        cross_database_outputs=cross_database_outputs,
        citation_outputs=citation_outputs,
        temporal_outputs=temporal_outputs,
        panr2_version=panr2_version,
    )
    return {"markdown": md_path, "html": html_path, "methods": methods_path, **dashboard_outputs}
