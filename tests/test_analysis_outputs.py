import os
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


class PanRAnalysisOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from panr2 import cli as panr_cli
        cls.panr = panr_cli
        cls.ncbi_clean = REPO_ROOT / "tests" / "fixtures" / "ncbi" / "ncbi_clean.csv"
        cls.summary_tab = REPO_ROOT / "tests" / "fixtures" / "abricate" / "ncbi_summary.tab"
        cls.results_tab = REPO_ROOT / "tests" / "fixtures" / "abricate" / "ncbi_results.tab"

    def build_tidy_df(self, tmpdir, min_identity=80):
        tmpdir = Path(tmpdir)
        summary_csv = tmpdir / "ncbi_summary.csv"
        results_csv = tmpdir / "ncbi_results.csv"
        self.panr.convert_tab_to_csv(str(self.summary_tab), str(summary_csv))
        self.panr.convert_tab_to_csv(str(self.results_tab), str(results_csv))

        merged = self.panr.load_and_merge_data(str(self.ncbi_clean), str(summary_csv))
        merged = self.panr.apply_analysis_filters(
            merged,
            min_identity=min_identity,
            output_dir=str(tmpdir),
            base_name="ncbi",
        )
        tidy = self.panr.convert_to_tidy_format(merged)

        results = pd.read_csv(results_csv, dtype=str)[["GENE", "RESISTANCE"]].drop_duplicates()
        tidy = tidy.merge(results, left_on="Gene", right_on="GENE", how="left").drop(columns=["GENE"])
        return tidy

    def test_filter_report_tracks_removed_identity_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.build_tidy_df(tmp, min_identity=90)
            report = pd.read_csv(Path(tmp) / "qc" / "ncbi_filter_report.csv")

        min_identity = report[report["filter"] == "min_identity"].iloc[0]
        self.assertTrue(bool(min_identity["enabled"]))
        self.assertEqual(int(min_identity["before"]), 6)
        self.assertEqual(int(min_identity["after"]), 3)
        self.assertEqual(int(min_identity["removed"]), 3)

    def test_comprehensive_analysis_outputs_expected_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            tidy = self.build_tidy_df(tmp, min_identity=80)
            outputs = self.panr.generate_comprehensive_analysis_outputs(
                tidy,
                tmp,
                "ncbi",
                "png",
                core_threshold=75,
                rare_threshold=25,
                top_n=10,
                cooccurrence_min_prevalence=0,
                cooccurrence_top_n=10,
            )
            for output_path in outputs.values():
                self.assertTrue(Path(output_path).exists(), output_path)

            analysis_dir = Path(tmp) / "analysis"
            burden = pd.read_csv(analysis_dir / "ncbi_sample_resistome_burden.csv")
            genes = pd.read_csv(analysis_dir / "ncbi_gene_prevalence_summary.csv")
            categories = pd.read_csv(analysis_dir / "ncbi_resistome_category_summary.csv")
            classes = pd.read_csv(analysis_dir / "ncbi_resistance_class_summary.csv")

        self.assertEqual(len(burden), 4)
        self.assertIn("unique_arg_genes", burden.columns)
        self.assertEqual(set(genes["Gene"]), {"blaA", "tetB"})
        self.assertIn("resistome_category", categories.columns)
        self.assertEqual(set(classes["resistance_class"]), {"beta-lactam", "tetracycline"})

    def test_cooccurrence_outputs_have_expected_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            tidy = self.build_tidy_df(tmp, min_identity=80)
            self.panr.generate_comprehensive_analysis_outputs(
                tidy, tmp, "ncbi", "png", core_threshold=75, rare_threshold=25,
                top_n=10, cooccurrence_min_prevalence=0, cooccurrence_top_n=10
            )
            analysis_dir = Path(tmp) / "analysis"
            gene_matrix = pd.read_csv(analysis_dir / "ncbi_gene_cooccurrence_matrix.csv", index_col=0)
            gene_pairs = pd.read_csv(analysis_dir / "ncbi_top_gene_pairs.csv")
            class_pairs = pd.read_csv(analysis_dir / "ncbi_top_class_pairs.csv")

        self.assertEqual(int(gene_matrix.loc["blaA", "blaA"]), 3)
        self.assertEqual(int(gene_matrix.loc["tetB", "tetB"]), 3)
        self.assertEqual(int(gene_matrix.loc["blaA", "tetB"]), 2)
        self.assertEqual(int(gene_pairs.iloc[0]["cooccurring_samples"]), 2)
        self.assertAlmostEqual(float(gene_pairs.iloc[0]["jaccard_index"]), 0.5)
        self.assertEqual(int(class_pairs.iloc[0]["cooccurring_samples"]), 2)

    def test_report_contains_journal_style_summary_sections(self):
        with tempfile.TemporaryDirectory() as tmp:
            tidy = self.build_tidy_df(tmp, min_identity=80)
            self.panr.generate_comprehensive_analysis_outputs(
                tidy, tmp, "ncbi", "png", core_threshold=75, rare_threshold=25,
                top_n=10, cooccurrence_min_prevalence=0, cooccurrence_top_n=10
            )
            from panr2.report import write_report
            outputs = write_report(
                tmp,
                "ncbi",
                options={"core_threshold": 75, "rare_threshold": 25, "min_identity": 80},
                panr2_version="test-version",
                input_files={"ncbi_clean": str(self.ncbi_clean), "abricate_summary": str(self.summary_tab), "abricate_results": str(self.results_tab)},
            )
            report_text = Path(outputs["markdown"]).read_text()

        self.assertIn("# PanR2 Panresistome Analysis Report", report_text)
        self.assertIn("## Executive Summary", report_text)
        self.assertIn("A total of 4 assemblies or samples", report_text)
        self.assertIn("blaA", report_text)
        self.assertIn("## Methods Summary", report_text)
        self.assertIn("## Reproducibility", report_text)

    def test_optional_vfdb_and_plasmidfinder_feature_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            from panr2.features import analyze_abricate_feature_database
            vfdb = analyze_abricate_feature_database(
                str(self.ncbi_clean), str(REPO_ROOT / "tests" / "fixtures" / "vfdb"), tmp, "vfdb", "virulence", min_identity=90
            )
            plasmid = analyze_abricate_feature_database(
                str(self.ncbi_clean), str(REPO_ROOT / "tests" / "fixtures" / "plasmidfinder"), tmp, "plasmidfinder", "plasmid", min_identity=90
            )
            vf_summary = pd.read_csv(vfdb["feature_summary"])
            plasmid_summary = pd.read_csv(plasmid["feature_summary"])
            self.assertIn("/vfdb/analysis/", vfdb["feature_summary"])
            self.assertIn("/vfdb/merged_output/", vfdb["merged"])
            self.assertIn("/vfdb/figures/", vfdb["feature_prevalence_plot"])
            self.assertIn("/plasmidfinder/analysis/", plasmid["feature_summary"])
            self.assertIn("/plasmidfinder/merged_output/", plasmid["merged"])
            self.assertIn("/plasmidfinder/figures/", plasmid["feature_prevalence_plot"])
            vf_geo = pd.read_csv(vfdb["geographic_summary"])
            plasmid_geo = pd.read_csv(plasmid["geographic_summary"])
            vf_pairs = pd.read_csv(vfdb["top_feature_pairs"])
            plasmid_pairs = pd.read_csv(plasmid["top_feature_pairs"])
            vf_group = pd.read_csv(vfdb["group_burden_summary"])
            plasmid_group = pd.read_csv(plasmid["group_burden_summary"])

            for outputs in [vfdb, plasmid]:
                self.assertTrue(Path(outputs["feature_prevalence_plot"]).exists())
                self.assertTrue(Path(outputs["category_prevalence_plot"]).exists())
                self.assertTrue(Path(outputs["presence_heatmap"]).exists())
                self.assertTrue(Path(outputs["identity_distribution_plot"]).exists())
                self.assertTrue(Path(outputs["feature_cooccurrence_heatmap"]).exists())
                self.assertTrue(Path(outputs["feature_prevalence_html"]).exists())
                self.assertTrue(Path(outputs["presence_heatmap_html"]).exists())
                self.assertTrue(Path(outputs["html_index"]).exists())
                self.assertTrue(Path(outputs["group_burden_summary"]).exists())
                self.assertTrue(Path(outputs["group_overall_tests"]).exists())
                self.assertTrue(Path(outputs["mean_burden_by_continent_plot"]).exists())
                self.assertTrue(Path(outputs["mean_burden_by_continent_html"]).exists())

        self.assertEqual(set(vf_summary["feature_id"]), {"espA", "stx2"})
        self.assertEqual(set(plasmid_summary["feature_id"]), {"IncFIB", "IncI1"})
        self.assertIn("min_identity", vf_summary.columns)
        self.assertIn("min_identity", plasmid_summary.columns)
        self.assertIn("geographic_level", vf_geo.columns)
        self.assertIn("geographic_level", plasmid_geo.columns)
        self.assertIn("jaccard_index", vf_pairs.columns)
        self.assertIn("jaccard_index", plasmid_pairs.columns)
        self.assertIn("grouping_variable", vf_group.columns)
        self.assertIn("grouping_variable", plasmid_group.columns)

    def test_optional_mobile_element_feature_analysis(self):
        databases = [
            ("mobileelementfinder", {"Tn3", "IS26"}),
            ("isfinder", {"IS26", "ISEcp1"}),
            ("integronfinder", {"intI1", "attC"}),
            ("iceberg", {"ICEKp1", "Tn916"}),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            from panr2.features import analyze_abricate_feature_database
            for feature_type, expected_features in databases:
                with self.subTest(feature_type=feature_type):
                    outputs = analyze_abricate_feature_database(
                        str(self.ncbi_clean),
                        str(REPO_ROOT / "tests" / "fixtures" / feature_type),
                        tmp,
                        feature_type,
                        "mge",
                        min_identity=90,
                    )
                    summary = pd.read_csv(outputs["feature_summary"])
                    categories = pd.read_csv(outputs["category_summary"])
                    group = pd.read_csv(outputs["group_burden_summary"])

                    self.assertEqual(set(summary["feature_id"]), expected_features)
                    self.assertIn("feature_category", summary.columns)
                    self.assertTrue((summary["feature_category"] != "Unknown").any())
                    self.assertFalse(categories.empty)
                    self.assertIn("grouping_variable", group.columns)
                    self.assertIn(f"/{feature_type}/analysis/", outputs["feature_summary"])
                    self.assertIn(f"/{feature_type}/merged_output/", outputs["merged"])
                    self.assertIn(f"/{feature_type}/figures/", outputs["feature_prevalence_plot"])
                    self.assertTrue(Path(outputs["feature_prevalence_plot"]).exists())
                    self.assertTrue(Path(outputs["category_prevalence_plot"]).exists())
                    self.assertTrue(Path(outputs["presence_heatmap"]).exists())
                    self.assertTrue(Path(outputs["identity_distribution_plot"]).exists())
                    self.assertTrue(Path(outputs["feature_cooccurrence_heatmap"]).exists())
                    self.assertTrue(Path(outputs["feature_prevalence_html"]).exists())
                    self.assertTrue(Path(outputs["presence_heatmap_html"]).exists())
                    self.assertTrue(Path(outputs["html_index"]).exists())


if __name__ == "__main__":
    unittest.main()
