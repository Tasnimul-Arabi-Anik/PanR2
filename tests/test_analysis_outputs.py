import os
import io
import tempfile
import unittest
from pathlib import Path
from contextlib import redirect_stdout

os.environ.setdefault("MPLBACKEND", "Agg")

import pandas as pd

from panr2.io import extract_assembly_accessions, normalize_assembly_accession


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

    def write_larger_fixture(self, root):
        root = Path(root)
        accessions = [f"GCF_10000{i}.1" for i in range(1, 11)]
        sample_rows = [
            (accessions[0], "SAMN100001", "United States", "North America", "Northern America", "2020", "clinical"),
            (accessions[1], "SAMN100002", "United States", "North America", "Northern America", "2021", "clinical"),
            (accessions[2], "SAMN100003", "Canada", "North America", "Northern America", "2022", "environmental"),
            (accessions[3], "SAMN100004", "Mexico", "North America", "Central America", "2023", "food"),
            (accessions[4], "SAMN100005", "Mexico", "North America", "Central America", "2024", "clinical"),
            (accessions[5], "SAMN100006", "Bangladesh", "Asia", "Southern Asia", "2020", "clinical"),
            (accessions[6], "SAMN100007", "Bangladesh", "Asia", "Southern Asia", "2021", "food"),
            (accessions[7], "SAMN100008", "India", "Asia", "Southern Asia", "2022", "clinical"),
            (accessions[8], "SAMN100009", "India", "Asia", "Southern Asia", "2023", "environmental"),
            (accessions[9], "SAMN100010", "Thailand", "Asia", "South-eastern Asia", "2024", "clinical"),
        ]
        ncbi_dir = root / "ncbi"
        ncbi_dir.mkdir(parents=True)
        ncbi_lines = [
            "Assembly Accession,Assembly BioSample Accession,Organism Name,Geographic Location,Continent,Subcontinent,Collection Date,Host,Isolation Source,Submitter,Assembly Level,Genome Representation,RefSeq Category"
        ]
        for accession, biosample, country, continent, subcontinent, year, source in sample_rows:
            ncbi_lines.append(
                f"{accession},{biosample},Example bacterium,{country},{continent},{subcontinent},{year},Homo sapiens,{source},Example Lab,Complete Genome,Full,na"
            )
        (ncbi_dir / "ncbi_clean.csv").write_text("\n".join(ncbi_lines) + "\n")

        def write_abricate_pair(db, features, hit_map, metadata, include_resistance=False):
            db_dir = root / db
            db_dir.mkdir(parents=True, exist_ok=True)
            summary_lines = ["#FILE\tNUM_FOUND\t" + "\t".join(features)]
            result_header = "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT"
            if include_resistance:
                result_header += "\tRESISTANCE"
            result_lines = [result_header]
            for sample_index, accession in enumerate(accessions, start=1):
                calls = hit_map.get(accession, {})
                values = [str(calls.get(feature, 0)) for feature in features]
                found = sum(1 for value in values if float(value) > 0)
                file_path = f"/synthetic/{accession}.fna"
                summary_lines.append(f"{file_path}\t{found}\t" + "\t".join(values))
                for feature_index, feature in enumerate(features, start=1):
                    identity = calls.get(feature, 0)
                    if float(identity) <= 0:
                        continue
                    product, category = metadata[feature]
                    start = 100 + feature_index * 100
                    end = start + 799
                    fields = [
                        file_path,
                        f"contig{sample_index}",
                        str(start),
                        str(end),
                        feature,
                        "1-800/800",
                        "=",
                        "0/0",
                        "100.00",
                        str(identity),
                        db,
                        f"{db.upper()}{feature_index:03d}",
                        product,
                    ]
                    if include_resistance:
                        fields.append(category)
                    result_lines.append("\t".join(fields))
            (db_dir / f"{db}_summary.tab").write_text("\n".join(summary_lines) + "\n")
            (db_dir / f"{db}_results.tab").write_text("\n".join(result_lines) + "\n")
            return db_dir

        def calls(patterns):
            return {accessions[index]: values for index, values in patterns.items()}

        ncbi_hits = calls({
            0: {"blaA": 97.5, "tetB": 95.0},
            1: {"blaA": 96.1, "sul1": 94.0},
            2: {"tetB": 93.5},
            3: {"blaA": 95.2, "qnrS": 92.0},
            4: {"blaA": 98.0, "tetB": 91.0, "sul1": 93.0},
            5: {"tetB": 96.0, "sul1": 95.5},
            6: {"blaA": 94.0, "qnrS": 93.2},
            7: {"tetB": 97.1, "qnrS": 94.4},
            8: {"sul1": 96.3},
            9: {"blaA": 95.8, "tetB": 92.6, "qnrS": 91.5},
        })
        write_abricate_pair(
            "abricate",
            ["blaA", "tetB", "sul1", "qnrS"],
            ncbi_hits,
            {
                "blaA": ("beta-lactamase", "beta-lactam"),
                "tetB": ("tetracycline efflux pump", "tetracycline"),
                "sul1": ("sulfonamide resistance protein", "sulfonamide"),
                "qnrS": ("quinolone resistance protein", "quinolone"),
            },
            include_resistance=True,
        )

        feature_specs = {
            "vfdb": (
                ["espA", "stx2", "fimH", "iutA"],
                {
                    "espA": ("type III secretion protein", "virulence"),
                    "stx2": ("Shiga toxin", "virulence"),
                    "fimH": ("type 1 fimbrial adhesin", "virulence"),
                    "iutA": ("aerobactin receptor", "virulence"),
                },
            ),
            "plasmidfinder": (
                ["IncFIB", "IncI1", "IncX3", "ColRNAI"],
                {
                    "IncFIB": ("IncFIB plasmid replicon", "plasmid"),
                    "IncI1": ("IncI1 plasmid replicon", "plasmid"),
                    "IncX3": ("IncX3 plasmid replicon", "plasmid"),
                    "ColRNAI": ("ColRNAI plasmid replicon", "plasmid"),
                },
            ),
            "mobileelementfinder": (
                ["Tn3", "IS26", "IS6100", "Tn21"],
                {
                    "Tn3": ("Tn3-family transposon", "mge"),
                    "IS26": ("IS6-family insertion sequence", "mge"),
                    "IS6100": ("IS6-family insertion sequence", "mge"),
                    "Tn21": ("Tn21-family transposon", "mge"),
                },
            ),
            "isfinder": (
                ["IS26", "ISEcp1", "IS6100", "IS903"],
                {
                    "IS26": ("IS6-family insertion sequence", "mge"),
                    "ISEcp1": ("ISEcp1 insertion sequence", "mge"),
                    "IS6100": ("IS6-family insertion sequence", "mge"),
                    "IS903": ("IS5-family insertion sequence", "mge"),
                },
            ),
            "integronfinder": (
                ["complete_integron_intI1", "attC1", "CALIN", "In0"],
                {
                    "complete_integron_intI1": ("class 1 complete integron", "mge"),
                    "attC1": ("attC recombination site", "mge"),
                    "CALIN": ("cluster of attC sites lacking integrase", "mge"),
                    "In0": ("integron integrase without cassette", "mge"),
                },
            ),
            "iceberg": (
                ["ICEKp1", "IME1", "CIME1", "Tn916"],
                {
                    "ICEKp1": ("integrative conjugative element", "mge"),
                    "IME1": ("integrative mobilizable element", "mge"),
                    "CIME1": ("cis-mobilizable element", "mge"),
                    "Tn916": ("conjugative transposon", "mge"),
                },
            ),
        }
        shared_patterns = [
            {0: 96.0, 1: 94.0},
            {0: 95.0, 2: 92.5},
            {1: 93.0},
            {0: 94.5, 3: 91.0},
            {0: 96.5, 1: 92.0, 2: 94.0},
            {1: 95.5, 2: 93.0},
            {0: 93.0, 3: 92.5},
            {1: 96.0, 3: 94.0},
            {2: 95.0},
            {0: 94.0, 1: 92.0, 3: 91.5},
        ]
        dirs = {"ncbi": ncbi_dir, "abricate": root / "abricate"}
        for db, (features, metadata) in feature_specs.items():
            hit_map = {}
            for sample_index, pattern in enumerate(shared_patterns):
                hit_map[accessions[sample_index]] = {
                    features[feature_index]: identity for feature_index, identity in pattern.items()
                }
            dirs[db] = write_abricate_pair(db, features, hit_map, metadata)
        return dirs

    def test_filter_report_tracks_removed_identity_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.build_tidy_df(tmp, min_identity=90)
            report = pd.read_csv(Path(tmp) / "qc" / "ncbi_filter_report.csv")

        min_identity = report[report["filter"] == "min_identity"].iloc[0]
        self.assertTrue(bool(min_identity["enabled"]))
        self.assertEqual(int(min_identity["before"]), 6)
        self.assertEqual(int(min_identity["after"]), 3)
        self.assertEqual(int(min_identity["removed"]), 3)

    def test_accession_normalization_preserves_versions_by_default(self):
        values = pd.Series(["/data/GCF_000123456.1_genomic.fna", "GCA_000987654"])
        accessions = extract_assembly_accessions(values)
        self.assertEqual(accessions.tolist(), ["GCF_000123456.1", "GCA_000987654"])
        self.assertEqual(normalize_assembly_accession("GCF_000123456.1"), "GCF_000123456.1")
        self.assertEqual(normalize_assembly_accession("GCF_000123456.1", preserve_version=False), "GCF_000123456")

    def test_doctor_reports_missing_external_tools_without_failing_analysis_only_install(self):
        stream = io.StringIO()
        with redirect_stdout(stream):
            exit_code = self.panr.run_doctor(
                abricate_bin="missing_abricate_for_panr2_test",
                mobileelementfinder_bin="missing_mefinder_for_panr2_test",
                integronfinder_bin="missing_integronfinder_for_panr2_test",
            )
        report = stream.getvalue()

        self.assertEqual(exit_code, 0)
        self.assertIn("PanR2 doctor report", report)
        self.assertIn("analysis-only", report)
        self.assertIn("not found on PATH", report)
        self.assertIn("PanR2 does not run ICEberg directly", report)

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
            vf_qc = pd.read_csv(vfdb["qc_summary"])
            plasmid_qc = pd.read_csv(plasmid["qc_summary"])

            for outputs in [vfdb, plasmid]:
                self.assertTrue(Path(outputs["feature_prevalence_plot"]).exists())
                self.assertTrue(Path(outputs["qc_summary"]).exists())
                self.assertTrue(Path(outputs["unmatched_samples"]).exists())
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
        self.assertIn("samples_with_at_least_one_feature", set(vf_qc["metric"]))
        self.assertIn("unmatched_metadata_samples", set(plasmid_qc["metric"]))

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

    def test_integrated_abricate_runner_feeds_panr2_analysis(self):
        script = """#!/usr/bin/env python3
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("abricate 1.0-test")
    sys.exit(0)
if args == ["--list"]:
    print("DATABASE\\tSEQUENCES\\tDBTYPE\\tDATE")
    print("ncbi\\t2\\tnucl\\t2026-May-02")
    sys.exit(0)
if args[:2] == ["--db", "ncbi"]:
    print("#FILE\\tSEQUENCE\\tSTART\\tEND\\tSTRAND\\tGENE\\tCOVERAGE\\tCOVERAGE_MAP\\tGAPS\\t%COVERAGE\\t%IDENTITY\\tDATABASE\\tACCESSION\\tPRODUCT\\tRESISTANCE")
    for path in args[2:]:
        if "GCF_000001.1" in path:
            print(f"{path}\\tcontig1\\t1\\t100\\t+\\tblaA\\t1-100/100\\t=\\t0/0\\t100.00\\t96.00\\tncbi\\tACC001\\tbeta lactamase\\tbeta-lactam")
        elif "GCA_000002.1" in path:
            print(f"{path}\\tcontig1\\t1\\t100\\t+\\ttetB\\t1-100/100\\t=\\t0/0\\t100.00\\t94.00\\tncbi\\tACC002\\ttetracycline efflux\\ttetracycline")
        elif "GCF_000003.1" in path:
            pass
        elif "GCA_000004.1" in path:
            print(f"{path}\\tcontig1\\t1\\t100\\t+\\tblaA\\t1-100/100\\t=\\t0/0\\t100.00\\t95.00\\tncbi\\tACC001\\tbeta lactamase\\tbeta-lactam")
            print(f"{path}\\tcontig2\\t1\\t100\\t+\\ttetB\\t1-100/100\\t=\\t0/0\\t100.00\\t93.00\\tncbi\\tACC002\\ttetracycline efflux\\ttetracycline")
    sys.exit(0)
if args and args[0] == "--summary":
    print("#FILE\\tNUM_FOUND\\tblaA\\ttetB")
    print("GCF_000001.1_genomic.fna\\t1\\t96.0\\t0")
    print("GCA_000002.1_genomic.fna\\t1\\t0\\t94.0")
    print("GCF_000003.1_genomic.fna\\t0\\t0\\t0")
    print("GCA_000004.1_genomic.fna\\t2\\t95.0\\t93.0")
    sys.exit(0)
print("unexpected arguments: " + " ".join(args), file=sys.stderr)
sys.exit(2)
"""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fake_abricate = tmp_path / "abricate"
            fake_abricate.write_text(script)
            fake_abricate.chmod(0o755)
            sequence_dir = tmp_path / "sequence"
            sequence_dir.mkdir()
            for accession in ["GCF_000001.1", "GCA_000002.1", "GCF_000003.1", "GCA_000004.1"]:
                (sequence_dir / f"{accession}_genomic.fna").write_text(">contig1\nATGC\n")

            output_dir = tmp_path / "panr2_output"
            self.panr.main(
                str(REPO_ROOT / "tests" / "fixtures" / "ncbi"),
                None,
                str(output_dir),
                "png",
                1,
                0,
                min_identity=90,
                min_samples_per_group=2,
                core_threshold=75,
                rare_threshold=25,
                top_n=10,
                cooccurrence_min_prevalence=0,
                cooccurrence_top_n=10,
                sequence_dir=str(sequence_dir),
                run_abricate=True,
                abricate_dbs=["ncbi"],
                abricate_bin=str(fake_abricate),
            )

            manifest = pd.read_csv(output_dir / "qc" / "panr2_tool_manifest.csv")
            report_text = (output_dir / "report" / "ncbi_panr2_report.md").read_text()
            self.assertEqual(manifest.iloc[0]["tool"], "abricate")
            self.assertEqual(manifest.iloc[0]["database"], "ncbi")
            self.assertTrue((output_dir / "tool_results" / "abricate" / "ncbi" / "ncbi_results.tab").exists())
            self.assertTrue((output_dir / "ncbi" / "analysis" / "ncbi_gene_prevalence_summary.csv").exists())
            self.assertIn("## Integrated Tool Runs", report_text)

    def test_integrated_mobileelementfinder_runner_feeds_feature_analysis(self):
        script = """#!/usr/bin/env python3
import os
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("MobileElementFinder 1.1-test")
    sys.exit(0)
if args and args[0] == "find":
    contig = args[args.index("--contig") + 1]
    output_prefix = args[-1]
    accession = os.path.basename(contig)
    with open(output_prefix + ".csv", "w") as handle:
        handle.write("# MobileElementFinder fake output\\n")
        handle.write("contig,start,end,mge_id,identity,coverage,type,accession\\n")
        if "GCF_000001.1" in accession:
            handle.write("contig1,10,900,Tn3,96.5,99.0,Tn3-family transposon,ME001\\n")
        elif "GCA_000002.1" in accession:
            handle.write("contig2,20,820,IS26,94.0,97.0,IS6-family insertion sequence,ME002\\n")
        elif "GCA_000004.1" in accession:
            handle.write("contig3,30,930,Tn3,95.0,98.0,Tn3-family transposon,ME001\\n")
            handle.write("contig3,940,1700,IS26,93.0,96.0,IS6-family insertion sequence,ME002\\n")
    sys.exit(0)
print("unexpected arguments: " + " ".join(args), file=sys.stderr)
sys.exit(2)
"""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fake_mefinder = tmp_path / "mefinder"
            fake_mefinder.write_text(script)
            fake_mefinder.chmod(0o755)
            sequence_dir = tmp_path / "sequence"
            sequence_dir.mkdir()
            for accession in ["GCF_000001.1", "GCA_000002.1", "GCF_000003.1", "GCA_000004.1"]:
                (sequence_dir / f"{accession}_genomic.fna").write_text(">contig1\nATGC\n")

            output_dir = tmp_path / "panr2_output"
            self.panr.main(
                str(REPO_ROOT / "tests" / "fixtures" / "ncbi"),
                str(REPO_ROOT / "tests" / "fixtures" / "abricate"),
                str(output_dir),
                "png",
                1,
                0,
                min_identity=90,
                min_samples_per_group=2,
                core_threshold=75,
                rare_threshold=25,
                top_n=10,
                cooccurrence_min_prevalence=0,
                cooccurrence_top_n=10,
                sequence_dir=str(sequence_dir),
                run_mobileelementfinder_tool=True,
                mobileelementfinder_bin=str(fake_mefinder),
            )

            manifest = pd.read_csv(output_dir / "qc" / "panr2_tool_manifest.csv")
            feature_summary = pd.read_csv(output_dir / "mobileelementfinder" / "analysis" / "mobileelementfinder_feature_summary.csv")
            report_text = (output_dir / "report" / "ncbi_panr2_report.md").read_text()
            self.assertEqual(manifest.iloc[0]["tool"], "mobileelementfinder")
            self.assertEqual(manifest.iloc[0]["database"], "mobileelementfinder")
            self.assertEqual(set(feature_summary["feature_id"]), {"Tn3", "IS26"})
            self.assertTrue((output_dir / "tool_results" / "mobileelementfinder" / "panr2_inputs" / "mobileelementfinder_results.tab").exists())
            self.assertTrue((output_dir / "mobileelementfinder" / "figures" / "index.html").exists())
            self.assertIn("MobileElementFinder", report_text)

    def test_integrated_integronfinder_runner_feeds_feature_analysis(self):
        script = """#!/usr/bin/env python3
import os
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("IntegronFinder 2.0-test")
    sys.exit(0)
if args:
    sequence = args[0]
    outdir = args[args.index("--outdir") + 1]
    os.makedirs(outdir, exist_ok=True)
    prefix = os.path.basename(sequence).replace(".fna", "")
    path = os.path.join(outdir, prefix + ".integrons")
    with open(path, "w") as handle:
        handle.write("ID_replicon\\tpos_beg\\tpos_end\\ttype\\tid\\tmodel\\n")
        if "GCF_000001.1" in sequence:
            handle.write("contig1\\t100\\t1100\\tcomplete_integron\\tintI1\\tclass 1 integron\\n")
        elif "GCA_000002.1" in sequence:
            handle.write("contig2\\t200\\t300\\tattC\\tattC1\\tattC recombination site\\n")
        elif "GCA_000004.1" in sequence:
            handle.write("contig3\\t300\\t1300\\tcomplete_integron\\tintI1\\tclass 1 integron\\n")
            handle.write("contig3\\t1400\\t1500\\tattC\\tattC1\\tattC recombination site\\n")
    sys.exit(0)
print("unexpected arguments", file=sys.stderr)
sys.exit(2)
"""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fake_integronfinder = tmp_path / "integron_finder"
            fake_integronfinder.write_text(script)
            fake_integronfinder.chmod(0o755)
            sequence_dir = tmp_path / "sequence"
            sequence_dir.mkdir()
            for accession in ["GCF_000001.1", "GCA_000002.1", "GCF_000003.1", "GCA_000004.1"]:
                (sequence_dir / f"{accession}_genomic.fna").write_text(">contig1\nATGC\n")

            output_dir = tmp_path / "panr2_output"
            self.panr.main(
                str(REPO_ROOT / "tests" / "fixtures" / "ncbi"),
                str(REPO_ROOT / "tests" / "fixtures" / "abricate"),
                str(output_dir),
                "png",
                1,
                0,
                min_identity=90,
                min_samples_per_group=2,
                core_threshold=75,
                rare_threshold=25,
                top_n=10,
                cooccurrence_min_prevalence=0,
                cooccurrence_top_n=10,
                sequence_dir=str(sequence_dir),
                run_integronfinder_tool=True,
                integronfinder_bin=str(fake_integronfinder),
            )

            manifest = pd.read_csv(output_dir / "qc" / "panr2_tool_manifest.csv")
            feature_summary = pd.read_csv(output_dir / "integronfinder" / "analysis" / "integronfinder_feature_summary.csv")
            report_text = (output_dir / "report" / "ncbi_panr2_report.md").read_text()
            self.assertEqual(manifest.iloc[0]["tool"], "integronfinder")
            self.assertEqual(manifest.iloc[0]["database"], "integronfinder")
            self.assertEqual(set(feature_summary["feature_id"]), {"attC1", "complete_integron_intI1"})
            self.assertTrue((output_dir / "tool_results" / "integronfinder" / "panr2_inputs" / "integronfinder_results.tab").exists())
            self.assertTrue((output_dir / "integronfinder" / "figures" / "index.html").exists())
            self.assertIn("IntegronFinder", report_text)

    def test_iceberg_table_converter_feeds_feature_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            iceberg_tables = tmp_path / "iceberg_tables"
            iceberg_tables.mkdir()
            (iceberg_tables / "iceberg_hits.tsv").write_text(
                "assembly_accession\tcontig\tstart\tend\tice_id\ttype\tidentity\tcoverage\tdescription\n"
                "GCF_000001.1\tcontig1\t100\t5000\tICEKp1\tICE\t95.0\t90.0\tKlebsiella integrative conjugative element\n"
                "GCA_000002.1\tcontig2\t200\t4200\tIME1\tIME\t93.0\t88.0\tintegrative mobilizable element\n"
                "GCA_000004.1\tcontig3\t300\t5200\tICEKp1\tICE\t94.0\t91.0\tKlebsiella integrative conjugative element\n"
                "GCA_000004.1\tcontig3\t5400\t8500\tCIME1\tCIME\t92.0\t87.0\tcis-mobilizable element\n"
            )

            output_dir = tmp_path / "panr2_output"
            self.panr.main(
                str(REPO_ROOT / "tests" / "fixtures" / "ncbi"),
                str(REPO_ROOT / "tests" / "fixtures" / "abricate"),
                str(output_dir),
                "png",
                1,
                0,
                min_identity=90,
                min_samples_per_group=2,
                core_threshold=75,
                rare_threshold=25,
                top_n=10,
                cooccurrence_min_prevalence=0,
                cooccurrence_top_n=10,
                iceberg_table_dir=str(iceberg_tables),
            )

            manifest = pd.read_csv(output_dir / "qc" / "panr2_tool_manifest.csv")
            feature_summary = pd.read_csv(output_dir / "iceberg" / "analysis" / "iceberg_feature_summary.csv")
            report_text = (output_dir / "report" / "ncbi_panr2_report.md").read_text()
            self.assertEqual(manifest.iloc[0]["tool"], "iceberg_table_converter")
            self.assertEqual(manifest.iloc[0]["database"], "iceberg")
            self.assertEqual(set(feature_summary["feature_id"]), {"ICEKp1", "IME1", "CIME1"})
            self.assertTrue((output_dir / "tool_results" / "iceberg" / "panr2_inputs" / "iceberg_results.tab").exists())
            self.assertTrue((output_dir / "iceberg" / "figures" / "index.html").exists())
            self.assertIn("ICEberg", report_text)
            self.assertIn("does not run an ICEberg annotation program directly", report_text)

    def test_full_fixture_workflow_outputs_all_database_reports(self):
        database_dirs = {
            "vfdb": "vfdb",
            "plasmidfinder": "plasmidfinder",
            "mobileelementfinder": "mobileelementfinder",
            "isfinder": "isfinder",
            "integronfinder": "integronfinder",
            "iceberg": "iceberg",
        }
        required_feature_outputs = [
            "analysis/{db}_feature_summary.csv",
            "analysis/{db}_category_summary.csv",
            "analysis/{db}_sample_burden.csv",
            "analysis/{db}_qc_summary.csv",
            "analysis/{db}_unmatched_samples.csv",
            "analysis/{db}_feature_cooccurrence_matrix.csv",
            "analysis/{db}_top_feature_pairs.csv",
            "analysis/{db}_group_burden_summary.csv",
            "analysis/{db}_group_overall_tests.csv",
            "merged_output/{db}_merged.csv",
            "merged_output/{db}_tidy.csv",
            "figures/{db}_feature_prevalence.png",
            "figures/{db}_category_prevalence.png",
            "figures/{db}_presence_heatmap.png",
            "figures/{db}_identity_distribution.png",
            "figures/{db}_feature_cooccurrence_heatmap.png",
            "figures/html_files/{db}_feature_prevalence.html",
            "figures/html_files/{db}_category_prevalence.html",
            "figures/html_files/{db}_presence_heatmap.html",
            "figures/html_files/{db}_identity_distribution.html",
            "figures/html_files/{db}_feature_cooccurrence_heatmap.html",
            "figures/index.html",
        ]

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "panr2_output"
            self.panr.main(
                str(REPO_ROOT / "tests" / "fixtures" / "ncbi"),
                str(REPO_ROOT / "tests" / "fixtures" / "abricate"),
                str(output_dir),
                "png",
                1,
                0,
                min_identity=90,
                min_samples_per_group=2,
                core_threshold=75,
                rare_threshold=25,
                top_n=10,
                cooccurrence_min_prevalence=0,
                cooccurrence_top_n=10,
                vfdb_dir=str(REPO_ROOT / "tests" / "fixtures" / "vfdb"),
                plasmidfinder_dir=str(REPO_ROOT / "tests" / "fixtures" / "plasmidfinder"),
                mobileelementfinder_dir=str(REPO_ROOT / "tests" / "fixtures" / "mobileelementfinder"),
                isfinder_dir=str(REPO_ROOT / "tests" / "fixtures" / "isfinder"),
                integronfinder_dir=str(REPO_ROOT / "tests" / "fixtures" / "integronfinder"),
                iceberg_dir=str(REPO_ROOT / "tests" / "fixtures" / "iceberg"),
            )

            report_text = (output_dir / "report" / "ncbi_panr2_report.md").read_text()
            self.assertTrue((output_dir / "report" / "ncbi_panr2_report.html").exists())
            self.assertTrue((output_dir / "qc" / "panr2_input_qc.csv").exists())
            self.assertTrue((output_dir / "ncbi" / "analysis" / "ncbi_gene_prevalence_summary.csv").exists())
            self.assertTrue((output_dir / "ncbi" / "figures" / "index.html").exists())
            self.assertIn("## Optional Database Feature Analysis", report_text)

            for db in database_dirs:
                with self.subTest(database=db):
                    self.assertIn(db, report_text)
                    for relative_path in required_feature_outputs:
                        expected = output_dir / db / relative_path.format(db=db)
                        self.assertTrue(expected.exists(), str(expected))

    def test_larger_fixture_workflow_exercises_group_statistics(self):
        optional_databases = [
            "vfdb",
            "plasmidfinder",
            "mobileelementfinder",
            "isfinder",
            "integronfinder",
            "iceberg",
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fixture_dirs = self.write_larger_fixture(tmp_path / "fixture")
            output_dir = tmp_path / "panr2_output"
            self.panr.main(
                str(fixture_dirs["ncbi"]),
                str(fixture_dirs["abricate"]),
                str(output_dir),
                "png",
                1,
                0,
                min_identity=90,
                min_samples_per_group=5,
                core_threshold=70,
                rare_threshold=30,
                top_n=10,
                cooccurrence_min_prevalence=0,
                cooccurrence_top_n=10,
                vfdb_dir=str(fixture_dirs["vfdb"]),
                plasmidfinder_dir=str(fixture_dirs["plasmidfinder"]),
                mobileelementfinder_dir=str(fixture_dirs["mobileelementfinder"]),
                isfinder_dir=str(fixture_dirs["isfinder"]),
                integronfinder_dir=str(fixture_dirs["integronfinder"]),
                iceberg_dir=str(fixture_dirs["iceberg"]),
            )

            report_path = sorted((output_dir / "report").glob("*_panr2_report.md"))[0]
            report_text = report_path.read_text()
            ncbi_burden = pd.read_csv(output_dir / "ncbi" / "analysis" / "abricate_sample_resistome_burden.csv")
            self.assertEqual(len(ncbi_burden), 10)
            self.assertIn("North America", set(ncbi_burden["Continent"]))
            self.assertIn("Asia", set(ncbi_burden["Continent"]))
            self.assertTrue((output_dir / "ncbi" / "figures" / "html_files" / "Continent_correlation_plot.html").exists())
            self.assertIn("## Optional Database Feature Analysis", report_text)

            for db in optional_databases:
                with self.subTest(database=db):
                    feature_summary = pd.read_csv(output_dir / db / "analysis" / f"{db}_feature_summary.csv")
                    group_summary = pd.read_csv(output_dir / db / "analysis" / f"{db}_group_burden_summary.csv")
                    overall_tests = pd.read_csv(output_dir / db / "analysis" / f"{db}_group_overall_tests.csv")
                    qc_summary = pd.read_csv(output_dir / db / "analysis" / f"{db}_qc_summary.csv")
                    self.assertGreaterEqual(len(feature_summary), 4)
                    self.assertGreaterEqual(int(feature_summary["present_samples"].max()), 5)
                    self.assertIn("top_20_features", set(qc_summary["metric"]))
                    self.assertIn("Continent", set(group_summary["grouping_variable"]))
                    self.assertIn("Continent", set(overall_tests["grouping_variable"]))
                    self.assertTrue((output_dir / db / "figures" / f"{db}_mean_burden_by_continent.png").exists())
                    self.assertTrue((output_dir / db / "figures" / "html_files" / f"{db}_mean_burden_by_continent.html").exists())


if __name__ == "__main__":
    unittest.main()
