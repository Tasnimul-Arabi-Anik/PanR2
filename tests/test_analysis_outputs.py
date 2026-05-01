import importlib.util
from importlib.machinery import SourceFileLoader
import os
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
PANR_PATH = REPO_ROOT / "bin" / "panr"


def load_panr_module():
    loader = SourceFileLoader("panr_cli", str(PANR_PATH))
    spec = importlib.util.spec_from_loader("panr_cli", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


class PanRAnalysisOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panr = load_panr_module()
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


if __name__ == "__main__":
    unittest.main()
