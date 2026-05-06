import html
import os
from datetime import datetime


MODULE_LABELS = {
    "amr": "AMR / NCBI",
    "vfdb": "VFDB virulence",
    "plasmidfinder": "PlasmidFinder",
    "mobileelementfinder": "MobileElementFinder",
    "isfinder": "ISfinder",
    "integronfinder": "IntegronFinder",
    "iceberg": "ICEberg-style inputs",
    "mlst": "MLST",
    "defensefinder": "DefenseFinder",
    "prophage": "Prophage / viral-region inputs",
    "mobsuite": "MOB-suite",
    "kleborate": "Kleborate",
    "kaptive": "Kaptive",
    "ectyper": "ECTyper",
    "serotypefinder": "SerotypeFinder",
    "sccmecfinder": "SCCmecFinder",
}


def _exists(output_dir, rel_path):
    return os.path.exists(os.path.join(output_dir, rel_path))


def _link(rel_path, label=None):
    label = label or rel_path
    return f'<li><a href="../{html.escape(rel_path)}">{html.escape(label)}</a></li>'


def _optional_link(output_dir, rel_path, label=None):
    return _link(rel_path, label) if _exists(output_dir, rel_path) else ""


def _section_links(output_dir, links):
    items = [_optional_link(output_dir, rel, label) for rel, label in links]
    items = [item for item in items if item]
    if not items:
        return "<p>No files were generated for this section.</p>"
    return "<ul>" + "\n".join(items) + "</ul>"


def _module_status(feature_outputs, cross_database_outputs, temporal_outputs, panresistome_context_outputs=None):
    panresistome_context_outputs = panresistome_context_outputs or {}
    modules = [{"module": "amr", "label": MODULE_LABELS["amr"], "status": "used"}]
    for key in ["vfdb", "plasmidfinder", "mobileelementfinder", "isfinder", "integronfinder", "iceberg", "mlst", "defensefinder", "prophage", "mobsuite", "kleborate", "kaptive", "ectyper", "serotypefinder", "sccmecfinder"]:
        modules.append({
            "module": key,
            "label": MODULE_LABELS[key],
            "status": "used" if key in feature_outputs else "not provided",
        })
    modules.append({
        "module": "cross_database",
        "label": "Cross-database associations",
        "status": "used" if cross_database_outputs else "not generated",
    })
    modules.append({
        "module": "metadata_enrichment",
        "label": "Metadata enrichment",
        "status": "used" if cross_database_outputs.get("feature_enrichment_by_metadata") else "not generated",
    })
    modules.append({
        "module": "temporal_trends",
        "label": "Advanced temporal trends",
        "status": "used" if temporal_outputs else "not generated",
    })
    modules.append({
        "module": "panresistome_context",
        "label": "PanResistome QC / ANI / QUAST context",
        "status": "used" if panresistome_context_outputs else "not detected",
    })
    return modules


def write_dashboard(output_dir, base_name, feature_outputs=None, cross_database_outputs=None, citation_outputs=None, temporal_outputs=None, panresistome_context_outputs=None, panr2_version="unknown"):
    """Write the top-level HTML dashboard users should open first."""
    feature_outputs = feature_outputs or {}
    cross_database_outputs = cross_database_outputs or {}
    citation_outputs = citation_outputs or {}
    temporal_outputs = temporal_outputs or {}
    panresistome_context_outputs = panresistome_context_outputs or {}
    report_dir = os.path.join(output_dir, "report")
    os.makedirs(report_dir, exist_ok=True)
    dashboard_path = os.path.join(report_dir, "index.html")

    module_rows = []
    for row in _module_status(feature_outputs, cross_database_outputs, temporal_outputs, panresistome_context_outputs):
        status_class = "used" if row["status"] == "used" else "missing"
        module_rows.append(
            f"<tr><td>{html.escape(row['label'])}</td><td class='{status_class}'>{html.escape(row['status'])}</td></tr>"
        )

    recommended = [
        ("qc/panr2_input_qc_summary.txt", "QC summary"),
        ("qc/metadata_bias_warning.txt", "Metadata bias warning"),
        ("panr2_inputs/metadata/fetchm2_clean.csv", "FetchM2 clean metadata"),
        ("panr2_inputs/metadata_analysis/metadata_analysis_report.md", "FetchM2 metadata analysis report"),
        ("panr2_inputs/metadata_audit/production_readiness_gate.md", "FetchM2 production readiness gate"),
        ("qc/qc_master_report.csv", "PanResistome QC master report"),
        ("panresistome_context/analysis/qc_context_sample_burden.csv", "QC context with feature burden"),
        ("panresistome_context/analysis/burden_by_ani_cluster.csv", "Feature burden by ANI cluster"),
        ("cross_database/analysis/amr_mge_associations.csv", "AMR-MGE associations"),
        ("cross_database/analysis/amr_plasmid_associations.csv", "AMR-plasmid associations"),
        ("cross_database/analysis/cross_database_top_associations.csv", "Top cross-database associations"),
        ("cross_database/figures/html_files/cross_database_feature_network.html", "Cross-database network"),
        ("temporal/analysis/temporal_feature_trends.csv", "Temporal feature trends"),
        ("mlst/analysis/st_feature_burden_summary.csv", "MLST feature-burden summary"),
        ("panresistome_context/analysis/duplicate_cluster_summary.csv", "Representative and duplicate clusters"),
        ("report/citations.md", "Citations"),
        ("report/software_versions.csv", "Software versions"),
    ]

    db_sections = []
    for key in ["ncbi", "vfdb", "plasmidfinder", "mobileelementfinder", "isfinder", "integronfinder", "iceberg", "mlst", "defensefinder", "prophage", "mobsuite", "kleborate", "kaptive", "ectyper", "serotypefinder", "sccmecfinder"]:
        label = MODULE_LABELS.get(key, key)
        links = [
            (f"{key}/analysis/{key}_feature_summary.csv", "Feature summary"),
            (f"{key}/analysis/{key}_sample_burden.csv", "Sample burden"),
            (f"{key}/analysis/{key}_qc_summary.csv", "Feature QC"),
            (f"{key}/analysis/{key}_group_burden_summary.csv", "Group burden"),
            (f"{key}/figures/index.html", "Figure index"),
            (f"{key}/merged_output/{key}_tidy.csv", "Tidy feature table"),
        ]
        if key == "ncbi":
            links = [
                (f"ncbi/analysis/{base_name}_gene_prevalence_summary.csv", "AMR gene prevalence"),
                (f"ncbi/analysis/{base_name}_sample_resistome_burden.csv", "AMR sample burden"),
                (f"ncbi/analysis/{base_name}_resistance_class_summary.csv", "Resistance classes"),
                ("ncbi/figures/index.html", "AMR figure index"),
            ]
        section_html = _section_links(output_dir, links)
        if "No files" not in section_html or key == "ncbi":
            db_sections.append(f"<h3>{html.escape(label)}</h3>{section_html}")

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>PanR2 Dashboard - {html.escape(base_name)}</title>
  <style>
    body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; color: #17202a; background: #f8fafc; }}
    header {{ background: #243447; color: white; padding: 24px 32px; }}
    main {{ max-width: 1180px; margin: 0 auto; padding: 24px; }}
    section {{ background: white; border: 1px solid #d8dee9; border-radius: 6px; padding: 18px 20px; margin-bottom: 18px; }}
    h1, h2, h3 {{ margin-top: 0; }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ padding: 8px 10px; border-bottom: 1px solid #e5e7eb; text-align: left; }}
    .used {{ color: #116329; font-weight: 700; }}
    .missing {{ color: #6b7280; }}
    .warning {{ border-left: 5px solid #b45309; background: #fff7ed; }}
    a {{ color: #155e75; }}
    ul {{ line-height: 1.65; }}
  </style>
</head>
<body>
  <header>
    <h1>PanR2 Analysis Dashboard</h1>
    <p>Dataset prefix: <strong>{html.escape(base_name)}</strong> | PanR2 {html.escape(panr2_version)} | Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
  </header>
  <main>
    <section>
      <h2>Run Summary</h2>
      <table><thead><tr><th>Module</th><th>Status</th></tr></thead><tbody>{''.join(module_rows)}</tbody></table>
    </section>
    <section>
      <h2>Recommended Outputs To Inspect First</h2>
      {_section_links(output_dir, recommended)}
    </section>
    <section class="warning">
      <h2>Biological Interpretation Warning</h2>
      <p>Cross-database co-occurrence is calculated at the sample/genome level. It does not prove physical linkage, plasmid localization, horizontal transfer, shared regulation, clinical phenotype, or causality. Treat AMR-MGE, AMR-plasmid, and AMR-virulence links as screening hypotheses unless genomic-context validation is performed.</p>
    </section>
    <section>
      <h2>QC, Metadata, And Reports</h2>
      {_section_links(output_dir, [
        ('qc/panr2_input_qc.csv', 'Input QC table'),
        ('qc/panr2_input_qc_summary.txt', 'Input QC summary'),
        ('qc/metadata_completeness_report.csv', 'Metadata completeness'),
        ('qc/metadata_group_sample_sizes.csv', 'Metadata group sample sizes'),
        ('qc/metadata_bias_warning.txt', 'Metadata bias warning'),
        ('panr2_inputs/metadata/fetchm2_clean.csv', 'FetchM2 clean metadata'),
        ('panr2_inputs/metadata_analysis/metadata_analysis_report.md', 'FetchM2 metadata analysis report'),
        ('panr2_inputs/metadata_audit/standardization_audit.md', 'FetchM2 standardization audit'),
        ('panr2_inputs/metadata_audit/production_readiness_gate.md', 'FetchM2 production readiness gate'),
        ('qc/qc_master_report.csv', 'PanResistome QC master report'),
        ('qc/excluded_for_panr2.csv', 'Samples excluded before PanR2'),
        (f'report/{base_name}_panr2_report.md', 'Journal-style Markdown report'),
        (f'report/{base_name}_panr2_report.html', 'Journal-style HTML report'),
      ])}
    </section>
    <section>
      <h2>PanResistome QC, ANI, And Assembly Context</h2>
      <p>These panels summarize heavy-tool outputs generated upstream by PanResistome. PanR2 uses them for comparative reporting without requiring FastANI, skani, QUAST, CheckM2, GTDB-Tk, Mash, or other external tools to be installed in the PanR2 environment.</p>
      {_section_links(output_dir, [
        ('panresistome_context/analysis/qc_context_sample_burden.csv', 'QC context with PanR2 feature burden'),
        ('panresistome_context/analysis/qc_master_status_summary.csv', 'QC master status summary'),
        ('panresistome_context/analysis/qc_feature_correlation_summary.csv', 'QC metric vs feature-burden correlations'),
        ('panresistome_context/analysis/species_consistency_summary.csv', 'ANI species-consistency summary'),
        ('panresistome_context/analysis/duplicate_cluster_summary.csv', 'Duplicate cluster summary'),
        ('panresistome_context/analysis/representative_samples.csv', 'Representative samples'),
        ('panresistome_context/analysis/burden_by_ani_cluster.csv', 'Feature burden by ANI cluster'),
        ('panresistome_context/analysis/panresistome_context_manifest.csv', 'Detected PanResistome context files'),
        ('ani/analysis/closest_genome.csv', 'Closest genome by ANI'),
        ('ani/analysis/duplicate_clusters.csv', 'Raw ANI duplicate clusters'),
        ('ani/analysis/ani_outliers.csv', 'ANI outliers'),
        ('quast/analysis/assembly_qc.csv', 'QUAST assembly QC'),
        ('assembly_qc/analysis/panr2_quast_summary.csv', 'PanR2-compatible QUAST summary'),
        ('ani/analysis/panr2_ani_summary.csv', 'PanR2-compatible ANI summary'),
      ])}
    </section>
    <section>
      <h2>Database-Specific Outputs</h2>
      {''.join(db_sections)}
    </section>
    <section>
      <h2>Cross-Database Outputs</h2>
      {_section_links(output_dir, [
        ('cross_database/analysis/cross_database_feature_matrix.csv', 'Unified feature matrix'),
        ('cross_database/analysis/cross_database_top_associations.csv', 'Top associations'),
        ('cross_database/analysis/amr_mge_associations.csv', 'AMR-MGE associations'),
        ('cross_database/analysis/amr_plasmid_associations.csv', 'AMR-plasmid associations'),
        ('cross_database/analysis/amr_virulence_associations.csv', 'AMR-virulence associations'),
        ('cross_database/analysis/amr_defense_associations.csv', 'AMR-defense associations'),
        ('cross_database/analysis/amr_prophage_associations.csv', 'AMR-prophage associations'),
        ('cross_database/analysis/sample_integrated_feature_burden.csv', 'Integrated sample burden'),
        ('cross_database/figures/html_files/global_feature_association_heatmap.html', 'Association heatmap'),
        ('cross_database/figures/html_files/cross_database_feature_network.html', 'Feature network'),
        ('cross_database/figures/figure_manifest.csv', 'Figure manifest'),
        ('cross_database/figures/plot_readability_warnings.csv', 'Plot readability warnings'),
      ])}
    </section>
    <section>
      <h2>Temporal And Typing Outputs</h2>
      {_section_links(output_dir, [
        ('temporal/analysis/temporal_feature_trends.csv', 'Temporal feature trends'),
        ('temporal/analysis/temporal_burden_trends.csv', 'Temporal burden trends'),
        ('temporal/figures/html_files/temporal_top_feature_trends.html', 'Temporal trend plot'),
        ('mlst/analysis/sample_mlst_summary.csv', 'MLST sample summary'),
        ('mlst/analysis/mlst_by_metadata.csv', 'MLST by metadata'),
        ('mlst/analysis/st_feature_burden_summary.csv', 'ST feature-burden summary'),
      ])}
    </section>
    <section>
      <h2>Citations And Versions</h2>
      {_section_links(output_dir, [
        ('report/citations.md', 'Citation report'),
        ('report/citations.bib', 'BibTeX citations'),
        ('report/software_versions.csv', 'Software versions'),
      ])}
    </section>
  </main>
</body>
</html>
"""
    with open(dashboard_path, "w") as handle:
        handle.write(html_text)
    return {"dashboard": dashboard_path}
