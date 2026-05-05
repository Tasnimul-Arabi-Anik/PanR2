import argparse
import glob
import importlib.util
import json
import logging
import os
import platform
import shutil
import subprocess
import sys
from importlib import metadata as importlib_metadata

import pandas as pd

from panr2.analysis import generate_comprehensive_analysis_outputs
from panr2.associations import generate_cross_database_associations
from panr2.citations import write_citation_outputs
from panr2.features import analyze_abricate_feature_database
from panr2.filters import apply_analysis_filters
from panr2.io import (
    convert_tab_to_csv,
    convert_to_tidy_format,
    load_and_merge_data,
    read_sample_map,
    save_merged_data,
    unique_input_files,
)
from panr2.metadata import write_metadata_reports
from panr2.mlst import analyze_mlst, finalize_mlst_analysis, run_mlst
from panr2.panresistome_context import generate_panresistome_context_outputs
from panr2.plots import (
    analyze_gene_presence,
    generate_comparison_heatmap,
    generate_comparison_heatmap_plotly,
    generate_gene_identity_boxplot,
    generate_gene_identity_boxplot_plotly,
    generate_geographic_resistance_map_plotly,
    generate_index_html,
    generate_mean_arg_lollipop,
    generate_mean_arg_lollipop_plotly,
    generate_resistance_barplot,
    mean_Arg_resistance_analysis_plotly,
)
from panr2.qc import write_input_qc_report
from panr2.report import write_report
from panr2.runners import convert_iceberg_tables, run_abricate_databases, run_integronfinder, run_mobileelementfinder
from panr2.stats import combined_correlation_analysis, correlation_scatterplot_analysis
from panr2.table_features import convert_defensefinder_tables, convert_prophage_tables, run_defensefinder
from panr2.temporal import write_temporal_trends


PANR2_VERSION = "0.1.3-dev"


# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def _python_package_status(package):
    if not importlib.util.find_spec(package):
        return {"status": "missing", "version": ""}
    try:
        version = importlib_metadata.version(package)
    except Exception:
        version = "available"
    return {"status": "ok", "version": version}


def _command_version(command, version_args):
    executable = shutil.which(command)
    if not executable:
        return "missing", "", ""
    for args in version_args:
        try:
            completed = subprocess.run(
                [executable, *args],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                timeout=20,
            )
            output = (completed.stdout or completed.stderr or "").strip().splitlines()
            if completed.returncode == 0 and output:
                return "ok", executable, output[0]
            if completed.returncode == 0:
                return "ok", executable, "available"
        except Exception as exc:
            return "error", executable, str(exc)
    return "found", executable, "version unavailable"


def _list_abricate_databases(abricate_path):
    databases = []
    if not abricate_path:
        return databases
    try:
        completed = subprocess.run(
            [abricate_path, "--list"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
            timeout=30,
        )
        if completed.returncode == 0:
            for line in completed.stdout.splitlines()[1:]:
                if not line.strip():
                    continue
                parts = line.split("\t")
                databases.append({
                    "name": parts[0],
                    "sequences": parts[1] if len(parts) > 1 else "",
                    "type": parts[2] if len(parts) > 2 else "",
                    "date": parts[3] if len(parts) > 3 else "",
                })
    except Exception:
        return []
    return databases


def collect_doctor_report(abricate_bin="abricate", mobileelementfinder_bin="mefinder", integronfinder_bin="integron_finder", mlst_bin="mlst", defensefinder_bin="defense-finder", include_system=False):
    """Collect dependency status for analysis-only and integrated-runner modes."""
    python_packages = ["pandas", "numpy", "matplotlib", "seaborn", "scipy", "plotly"]
    package_rows = {package: _python_package_status(package) for package in python_packages}

    tool_specs = [
        ("abricate", abricate_bin, [["--version"]]),
        ("mobileelementfinder", mobileelementfinder_bin, [["--version"], ["version"]]),
        ("integronfinder", integronfinder_bin, [["--version"], ["--help"]]),
        ("mlst", mlst_bin, [["--version"]]),
        ("defensefinder", defensefinder_bin, [["--version"], ["--help"]]),
    ]
    tool_rows = {}
    abricate_path = None
    for label, command, version_args in tool_specs:
        status, path, version = _command_version(command, version_args)
        if label == "abricate":
            abricate_path = path
        tool_rows[label] = {
            "command": command,
            "status": status,
            "path": path,
            "version": version,
        }

    report = {
        "panr2_version": PANR2_VERSION,
        "python_packages": package_rows,
        "tools": tool_rows,
        "abricate_databases": _list_abricate_databases(abricate_path),
        "install_modes": {
            "analysis_only": "Requires PanR2 Python dependencies plus existing ABRicate-style result folders.",
            "integrated_runners": "Additionally requires ABRicate, MobileElementFinder, IntegronFinder, MLST, and/or DefenseFinder installed with their databases/schemes.",
            "iceberg_style": "Requires user-provided ICE/IME/CIME annotation tables or ABRicate-style ICEberg inputs; PanR2 does not run ICEberg directly.",
        },
    }
    if include_system:
        report["system"] = {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "executable": sys.executable,
            "cwd": os.getcwd(),
            "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
        }
    return report


def _print_doctor_report(report):
    print(f"PanR2 doctor report (PanR2 {report['panr2_version']})")
    if "system" in report:
        print("")
        print("System:")
        for key, value in report["system"].items():
            print(f"- {key}: {value}")

    print("")
    print("Python package dependencies:")
    for package, details in report["python_packages"].items():
        version = f" ({details['version']})" if details.get("version") else ""
        print(f"- {package}: {details['status']}{version}")
    print("")
    print("External annotation tools:")
    for label, details in report["tools"].items():
        command = details["command"]
        status = details["status"]
        path = details["path"]
        version = details["version"]
        detail = f"{path} ({version})" if path else "not found on PATH"
        print(f"- {label} [{command}]: {status}; {detail}")
    if report["tools"]["abricate"]["path"]:
        databases = report["abricate_databases"]
        if databases:
            names = [row["name"] for row in databases[:10]]
            print(f"- abricate databases: {len(databases)} available ({', '.join(names)})")
        else:
            print("- abricate databases: none detected; run `panr setup-db` or `abricate --setupdb` before integrated ABRicate analysis")
    print("")
    print("Install modes:")
    print(f"- analysis-only: {report['install_modes']['analysis_only']}")
    print(f"- integrated runners: {report['install_modes']['integrated_runners']}")
    print(f"- ICEberg-style analysis: {report['install_modes']['iceberg_style']}")


def setup_abricate_databases(abricate_bin="abricate", dbs=None, check_only=False, json_output=False):
    executable = shutil.which(abricate_bin)
    result = {"tool": "abricate", "command": abricate_bin, "path": executable, "check_only": check_only, "requested_databases": dbs or []}
    if not executable:
        result["status"] = "missing"
        result["message"] = "ABRicate was not found. Install with environment.yml/Docker or run `panr doctor` for details."
        if json_output:
            print(json.dumps(result, indent=2, sort_keys=True))
        else:
            print(result["message"])
        return 1

    before = _list_abricate_databases(executable)
    result["databases_before"] = before
    if not check_only:
        completed = subprocess.run([executable, "--setupdb"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        result["setupdb_returncode"] = completed.returncode
        result["setupdb_stdout"] = completed.stdout.strip()
        result["setupdb_stderr"] = completed.stderr.strip()
        if completed.returncode != 0:
            result["status"] = "failed"
            result["message"] = "ABRicate database setup failed. Check network access and ABRicate installation."
            if json_output:
                print(json.dumps(result, indent=2, sort_keys=True))
            else:
                print(result["message"])
                if result["setupdb_stderr"]:
                    print(result["setupdb_stderr"])
            return completed.returncode

    after = _list_abricate_databases(executable)
    result["databases_after"] = after
    available_names = {row["name"] for row in after}
    requested = set(dbs or [])
    missing = sorted(requested - available_names)
    result["missing_requested_databases"] = missing
    result["status"] = "ok" if not missing else "missing_requested_databases"
    if json_output:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"ABRicate path: {executable}")
        print(f"ABRicate databases available: {len(after)}")
        if after:
            print("Databases: " + ", ".join(row["name"] for row in after[:20]))
        if missing:
            print("Missing requested databases: " + ", ".join(missing))
    return 0 if not missing else 1


def _has_abricate_feature_files(feature_dir):
    if not feature_dir or not os.path.isdir(feature_dir):
        return False
    names = [name.lower() for name in os.listdir(feature_dir) if name.lower().endswith((".csv", ".tab"))]
    return any("summary" in name for name in names) and any("results" in name for name in names)


def run_doctor(abricate_bin="abricate", mobileelementfinder_bin="mefinder", integronfinder_bin="integron_finder", mlst_bin="mlst", defensefinder_bin="defense-finder", json_output=False, fix=False, include_system=False):
    """Print dependency status for analysis-only and integrated-runner modes."""
    report = collect_doctor_report(abricate_bin, mobileelementfinder_bin, integronfinder_bin, mlst_bin=mlst_bin, defensefinder_bin=defensefinder_bin, include_system=include_system)
    fixes = []
    if fix and report["tools"]["abricate"]["path"] and not report["abricate_databases"]:
        executable = report["tools"]["abricate"]["path"]
        completed = subprocess.run([executable, "--setupdb"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        fixes.append({
            "action": "abricate --setupdb",
            "returncode": completed.returncode,
            "status": "ok" if completed.returncode == 0 else "failed",
            "stderr": completed.stderr.strip(),
        })
        report = collect_doctor_report(abricate_bin, mobileelementfinder_bin, integronfinder_bin, mlst_bin=mlst_bin, defensefinder_bin=defensefinder_bin, include_system=include_system)
    if fixes:
        report["fixes"] = fixes

    if json_output:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        for fix_result in fixes:
            print(f"Fix action `{fix_result['action']}`: {fix_result['status']}")
            if fix_result.get("stderr") and fix_result["status"] != "ok":
                print(fix_result["stderr"])
        _print_doctor_report(report)

    missing_packages = [name for name, details in report["python_packages"].items() if details["status"] != "ok"]
    return 1 if missing_packages else 0


def _run_subcommand(argv):
    command = argv[0]
    if command in {"doctor", "install-info", "setup"}:
        parser = argparse.ArgumentParser(prog=f"panr {command}", description="Check PanR2 installation readiness.")
        parser.add_argument("--json", action="store_true", help="Write machine-readable JSON output.")
        parser.add_argument("--fix", action="store_true", help="Run safe fixes when possible, currently ABRicate database setup if ABRicate is installed.")
        parser.add_argument("--mode", choices=["analysis-only", "integrated", "databases-only"], default="integrated", help="Setup/readiness mode label for user guidance.")
        parser.add_argument("--abricate-bin", default="abricate", help="ABRicate executable name or path.")
        parser.add_argument("--mobileelementfinder-bin", default="mefinder", help="MobileElementFinder executable name or path.")
        parser.add_argument("--integronfinder-bin", default="integron_finder", help="IntegronFinder executable name or path.")
        parser.add_argument("--mlst-bin", default="mlst", help="MLST executable name or path.")
        parser.add_argument("--defensefinder-bin", default="defense-finder", help="DefenseFinder executable name or path.")
        args = parser.parse_args(argv[1:])
        include_system = command == "install-info"
        if command == "setup" and args.mode == "databases-only":
            return setup_abricate_databases(abricate_bin=args.abricate_bin, check_only=not args.fix, json_output=args.json)
        return run_doctor(
            abricate_bin=args.abricate_bin,
            mobileelementfinder_bin=args.mobileelementfinder_bin,
            integronfinder_bin=args.integronfinder_bin,
            mlst_bin=args.mlst_bin,
            defensefinder_bin=args.defensefinder_bin,
            json_output=args.json,
            fix=args.fix,
            include_system=include_system,
        )
    if command == "setup-db":
        parser = argparse.ArgumentParser(prog="panr setup-db", description="Check or initialize ABRicate databases for integrated PanR2 runs.")
        parser.add_argument("--abricate-bin", default="abricate", help="ABRicate executable name or path.")
        parser.add_argument("--dbs", default="", help="Comma-separated ABRicate database names expected after setup, for example ncbi,vfdb,plasmidfinder.")
        parser.add_argument("--check-only", action="store_true", help="Only report database visibility without running abricate --setupdb.")
        parser.add_argument("--json", action="store_true", help="Write machine-readable JSON output.")
        args = parser.parse_args(argv[1:])
        dbs = [db.strip() for db in args.dbs.split(",") if db.strip()]
        return setup_abricate_databases(abricate_bin=args.abricate_bin, dbs=dbs, check_only=args.check_only, json_output=args.json)
    if command == "citations":
        parser = argparse.ArgumentParser(prog="panr citations", description="Write PanR2 citation and software-version files for an output directory.")
        parser.add_argument("--output-dir", required=True, help="PanR2 output directory.")
        args = parser.parse_args(argv[1:])
        outputs = write_citation_outputs(args.output_dir, panr2_version=PANR2_VERSION)
        print(f"Citation report written to {outputs['citations_md']}")
        print(f"BibTeX file written to {outputs['citations_bib']}")
        print(f"Software versions written to {outputs['software_versions']}")
        return 0
    if command == "validate-demo":
        parser = argparse.ArgumentParser(prog="panr validate-demo", description="Run PanR2 on the bundled small validation dataset.")
        parser.add_argument("--output-dir", required=True, help="Directory where validation outputs will be written.")
        parser.add_argument("--format", default="png", choices=["tiff", "svg", "png", "pdf"], help="Output format for figures.")
        args = parser.parse_args(argv[1:])
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
        fixtures = os.path.join(repo_root, "tests", "fixtures")
        if not os.path.isdir(fixtures):
            raise FileNotFoundError("Bundled test fixtures were not found. Run validate-demo from a source checkout.")
        main(
            os.path.join(fixtures, "ncbi"),
            os.path.join(fixtures, "abricate"),
            args.output_dir,
            args.format,
            1,
            0,
            min_identity=90,
            min_samples_per_group=2,
            core_threshold=75,
            rare_threshold=25,
            top_n=10,
            cooccurrence_min_prevalence=0,
            cooccurrence_top_n=10,
            vfdb_dir=os.path.join(fixtures, "vfdb"),
            plasmidfinder_dir=os.path.join(fixtures, "plasmidfinder"),
            mobileelementfinder_dir=os.path.join(fixtures, "mobileelementfinder"),
            isfinder_dir=os.path.join(fixtures, "isfinder"),
            integronfinder_dir=os.path.join(fixtures, "integronfinder"),
            iceberg_dir=os.path.join(fixtures, "iceberg"),
            mlst_dir=os.path.join(fixtures, "mlst"),
            defensefinder_dir=os.path.join(fixtures, "defensefinder_tables"),
            prophage_dir=os.path.join(fixtures, "prophage_tables"),
            sample_map_path=os.path.join(fixtures, "sample_map.csv"),
        )
        print(f"Validation demo written to {args.output_dir}")
        print(f"Open {os.path.join(args.output_dir, 'report', 'index.html')}")
        return 0
    raise ValueError(f"Unknown PanR2 subcommand: {command}")


def main(ncbi_dir, abricate_dir, output_dir, fig_format, nseq, genep, min_identity=0.0, drop_unmatched_accessions=False, min_samples_per_group=5, core_threshold=95.0, rare_threshold=5.0, top_n=25, cooccurrence_min_prevalence=0.0, cooccurrence_top_n=25, vfdb_dir=None, plasmidfinder_dir=None, mobileelementfinder_dir=None, isfinder_dir=None, integronfinder_dir=None, iceberg_dir=None, mlst_dir=None, defensefinder_dir=None, prophage_dir=None, sample_map_path=None, sequence_dir=None, run_abricate=False, abricate_dbs=None, abricate_bin="abricate", abricate_summary_metric="identity", run_mobileelementfinder_tool=False, mobileelementfinder_bin="mefinder", mobileelementfinder_threads=1, run_integronfinder_tool=False, integronfinder_bin="integron_finder", integronfinder_threads=1, run_mlst_tool=False, mlst_bin="mlst", run_defensefinder_tool=False, defensefinder_bin="defense-finder", iceberg_table_dir=None, force_tool_run=False, run_cross_database=True, cross_database_max_features=300, plot_style="publication", label_max_length=None, run_temporal_trends=True):
    """Main function to process data and generate outputs."""
    logging.info("Starting the script.")

    if min_identity < 0 or min_identity > 100:
        raise ValueError("--min-identity must be between 0 and 100.")
    if min_samples_per_group < 1:
        raise ValueError("--min-samples-per-group must be at least 1.")
    if core_threshold < 0 or core_threshold > 100:
        raise ValueError("--core-threshold must be between 0 and 100.")
    if rare_threshold < 0 or rare_threshold > 100:
        raise ValueError("--rare-threshold must be between 0 and 100.")
    if rare_threshold > core_threshold:
        raise ValueError("--rare-threshold cannot be greater than --core-threshold.")
    if top_n < 1:
        raise ValueError("--top-n must be at least 1.")
    if cooccurrence_min_prevalence < 0 or cooccurrence_min_prevalence > 100:
        raise ValueError("--cooccurrence-min-prevalence must be between 0 and 100.")
    if cooccurrence_top_n < 1:
        raise ValueError("--cooccurrence-top-n must be at least 1.")

    sample_map = read_sample_map(sample_map_path) if sample_map_path else {}

    tool_manifest = {}
    if run_abricate:
        if not sequence_dir:
            raise ValueError("--sequence-dir is required when --run-abricate is used.")
        selected_dbs = abricate_dbs or ["ncbi"]
        abricate_run = run_abricate_databases(
            sequence_dir,
            output_dir,
            selected_dbs,
            abricate_bin=abricate_bin,
            summary_metric=abricate_summary_metric,
            force=force_tool_run,
        )
        tool_manifest = abricate_run.get("manifest", {})
        generated_dirs = abricate_run["database_dirs"]
        if not abricate_dir and "ncbi" in generated_dirs:
            abricate_dir = generated_dirs["ncbi"]
        if not vfdb_dir and "vfdb" in generated_dirs:
            vfdb_dir = generated_dirs["vfdb"]
        if not plasmidfinder_dir and "plasmidfinder" in generated_dirs:
            plasmidfinder_dir = generated_dirs["plasmidfinder"]
        if not isfinder_dir and "isfinder" in generated_dirs:
            isfinder_dir = generated_dirs["isfinder"]
    if run_mobileelementfinder_tool:
        if not sequence_dir:
            raise ValueError("--sequence-dir is required when --run-mobileelementfinder is used.")
        mobileelementfinder_run = run_mobileelementfinder(
            sequence_dir,
            output_dir,
            mefinder_bin=mobileelementfinder_bin,
            threads=mobileelementfinder_threads,
            force=force_tool_run,
        )
        tool_manifest = mobileelementfinder_run.get("manifest", tool_manifest)
        if not mobileelementfinder_dir:
            mobileelementfinder_dir = mobileelementfinder_run["feature_dir"]
    if run_integronfinder_tool:
        if not sequence_dir:
            raise ValueError("--sequence-dir is required when --run-integronfinder is used.")
        integronfinder_run = run_integronfinder(
            sequence_dir,
            output_dir,
            integronfinder_bin=integronfinder_bin,
            cpu=integronfinder_threads,
            force=force_tool_run,
        )
        tool_manifest = integronfinder_run.get("manifest", tool_manifest)
        if not integronfinder_dir:
            integronfinder_dir = integronfinder_run["feature_dir"]
    if run_mlst_tool:
        if not sequence_dir:
            raise ValueError("--sequence-dir is required when --run-mlst is used.")
        mlst_run = run_mlst(
            sequence_dir,
            output_dir,
            mlst_bin=mlst_bin,
            force=force_tool_run,
        )
        tool_manifest = mlst_run.get("manifest", tool_manifest)
        if not mlst_dir:
            mlst_dir = mlst_run["mlst_dir"]
    if run_defensefinder_tool:
        if not sequence_dir:
            raise ValueError("--sequence-dir is required when --run-defensefinder is used.")
        defensefinder_run = run_defensefinder(
            sequence_dir,
            output_dir,
            defensefinder_bin=defensefinder_bin,
            force=force_tool_run,
        )
        tool_manifest = defensefinder_run.get("manifest", tool_manifest)
        if not defensefinder_dir:
            defensefinder_dir = defensefinder_run["feature_dir"]
    if iceberg_table_dir:
        iceberg_run = convert_iceberg_tables(iceberg_table_dir, output_dir, sample_map=sample_map)
        tool_manifest = iceberg_run.get("manifest", tool_manifest)
        if not iceberg_dir:
            iceberg_dir = iceberg_run["feature_dir"]

    if not abricate_dir:
        raise ValueError("--abricate-dir is required unless --run-abricate includes the ncbi database.")
    
    # Define paths
    ncbi_clean_path = os.path.join(ncbi_dir, "ncbi_clean.csv")
    abricate_summary_files = unique_input_files(
        glob.glob(os.path.join(abricate_dir, "*summary.[ct][sa][bv]"))
    )  # Match .csv and .tab files
    abricate_results_files = unique_input_files(
        glob.glob(os.path.join(abricate_dir, "*results.[ct][sa][bv]"))
    )
    
    if not os.path.exists(ncbi_clean_path):
        raise FileNotFoundError(f"ncbi_clean.csv not found in {ncbi_dir}.")
    
    if not abricate_summary_files:
        raise FileNotFoundError(f"No CSV or TAB summary files (abricate) found in {abricate_dir}.")
        
    if not abricate_results_files:
        raise FileNotFoundError(f"No CSV or TAB results files (abricate) found in {abricate_dir}.")
    
    # Create database-specific output subdirectories. Shared QC and reports stay top-level.
    ncbi_output_dir = os.path.join(output_dir, "ncbi")
    merged_output_dir = os.path.join(ncbi_output_dir, "merged_output")
    figures_dir = os.path.join(ncbi_output_dir, "figures")
    os.makedirs(merged_output_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    first_summary_file = sorted(abricate_summary_files)[0]
    first_results_file = sorted(abricate_results_files)[0]
    write_input_qc_report(ncbi_clean_path, first_summary_file, first_results_file, output_dir)
    metadata_report_outputs = write_metadata_reports(ncbi_clean_path, output_dir, min_group_size=min_samples_per_group)

    optional_feature_outputs = {}
    mlst_outputs = {}
    if mlst_dir:
        mlst_outputs = analyze_mlst(ncbi_clean_path, mlst_dir, output_dir, fig_format=fig_format, sample_map=sample_map)
        optional_feature_outputs["mlst"] = mlst_outputs
    if defensefinder_dir and not _has_abricate_feature_files(defensefinder_dir):
        defensefinder_dir = convert_defensefinder_tables(defensefinder_dir, output_dir, sample_map=sample_map)["feature_dir"]
    if prophage_dir and not _has_abricate_feature_files(prophage_dir):
        prophage_dir = convert_prophage_tables(prophage_dir, output_dir, sample_map=sample_map)["feature_dir"]
    optional_database_specs = [
        ("vfdb", vfdb_dir, "virulence"),
        ("plasmidfinder", plasmidfinder_dir, "plasmid"),
        ("mobileelementfinder", mobileelementfinder_dir, "mge"),
        ("isfinder", isfinder_dir, "mge"),
        ("integronfinder", integronfinder_dir, "mge"),
        ("iceberg", iceberg_dir, "mge"),
        ("defensefinder", defensefinder_dir, "defense"),
        ("prophage", prophage_dir, "prophage"),
    ]
    for feature_type, feature_dir, mode in optional_database_specs:
        if feature_dir:
            optional_feature_outputs[feature_type] = analyze_abricate_feature_database(
                ncbi_clean_path,
                feature_dir,
                output_dir,
                feature_type,
                mode,
                min_identity=min_identity,
                fig_format=fig_format,
                sample_map=sample_map,
            )
    
    for abricate_summary_file in abricate_summary_files:
        try:
            # Get base name without extension
            base = os.path.splitext(abricate_summary_file)[0]
            csv_file = base + ".csv"
            tab_file = base + ".tab"

            if os.path.exists(csv_file):
                abricate_summary_file = csv_file
            elif os.path.exists(tab_file):
                convert_tab_to_csv(tab_file, csv_file)
                abricate_summary_file = csv_file
            else:
                print(f"⚠️ No .csv or .tab found for {os.path.basename(base)}, skipping.")
                continue

            # Load, filter, and merge data
            merged_df = load_and_merge_data(ncbi_clean_path, abricate_summary_file, sample_map=sample_map, sample_map_qc_dir=output_dir, sample_map_source="abricate")
            output_filename = f"ncbi_{os.path.basename(abricate_summary_file)}"
            base_name = os.path.basename(abricate_summary_file).replace("_summary.csv", "")
            merged_df = apply_analysis_filters(
                merged_df,
                min_identity=min_identity,
                drop_unmatched_accessions=drop_unmatched_accessions,
                output_dir=output_dir,
                base_name=base_name,
            )
            
            # Save merged data
            save_merged_data(merged_df, merged_output_dir, output_filename)
            
            # Convert to tidy format
            tidy_df = convert_to_tidy_format(merged_df)

            # Add RESISTANCE column from corresponding results file
            basename = os.path.basename(abricate_summary_file).replace("_summary.csv", "")
            expected_results_file = os.path.join(abricate_dir, f"{basename}_results.csv")

            # Convert .tab to .csv if needed
            if not os.path.exists(expected_results_file):
                tab_version = expected_results_file.replace(".csv", ".tab")
                if os.path.exists(tab_version):
                    convert_tab_to_csv(tab_version, expected_results_file)
                    logging.info(f"Converted {tab_version} to CSV.")
                else:
                    logging.warning(f"Results file not found for {basename}. Skipping RESISTANCE enrichment.")
                    expected_results_file = None
            # Merge RESISTANCE and GENE data
            if expected_results_file and os.path.exists(expected_results_file):
                try:
                    results_df = pd.read_csv(expected_results_file, dtype=str)
                    results_df = results_df[['GENE', 'RESISTANCE']].drop_duplicates()

                    # Perform the merge using 'Gene' from tidy_df and 'GENE' from results_df
                    tidy_df = tidy_df.merge(results_df, left_on='Gene', right_on='GENE', how='left')

                    # Drop redundant GENE column
                    tidy_df.drop(columns=['GENE'], inplace=True)

                    logging.info(f"Successfully added RESISTANCE info for {basename}.")
                except Exception as e:
                    logging.error(f"Error merging RESISTANCE data from {expected_results_file}: {e}")

            tidy_file = os.path.join(merged_output_dir, output_filename.replace("_summary.csv", "_tidy_summary.csv"))
            tidy_df.to_csv(tidy_file, index=False)
            logging.info(f"Tidied file saved to {tidy_file}")

            generate_comprehensive_analysis_outputs(
                tidy_df,
                ncbi_output_dir,
                base_name,
                fig_format,
                core_threshold=core_threshold,
                rare_threshold=rare_threshold,
                top_n=top_n,
                cooccurrence_min_prevalence=cooccurrence_min_prevalence,
                cooccurrence_top_n=cooccurrence_top_n,
            )
            cross_database_outputs = {}
            if run_cross_database and optional_feature_outputs:
                cross_database_outputs = generate_cross_database_associations(
                    output_dir,
                    base_name,
                    tidy_df,
                    optional_feature_outputs,
                    fig_format=fig_format,
                    top_n=top_n,
                    min_prevalence=cooccurrence_min_prevalence,
                    max_features=cross_database_max_features,
                    min_group_size=min_samples_per_group,
                    plot_style=plot_style,
                    label_max_length=label_max_length,
                )
            temporal_outputs = {}
            if run_temporal_trends:
                temporal_outputs = write_temporal_trends(
                    output_dir,
                    base_name,
                    tidy_df,
                    optional_feature_outputs,
                    fig_format=fig_format,
                )
            if mlst_outputs:
                mlst_outputs = finalize_mlst_analysis(output_dir, mlst_outputs, cross_database_outputs)
                optional_feature_outputs["mlst"] = mlst_outputs
            panresistome_context_outputs = generate_panresistome_context_outputs(
                output_dir,
                base_name,
                cross_database_outputs=cross_database_outputs,
            )
            
            # Analyze gene prevalence and generate figures
            analyze_gene_presence(tidy_df, figures_dir, base_name, fig_format)
            

            # Generate boxplot for gene identity
            generate_gene_identity_boxplot(tidy_file, figures_dir, fig_format)
            # Generate interactive boxplot for gene identity
            generate_gene_identity_boxplot_plotly(tidy_file, figures_dir)

            
            # Create subdirectory 'mean_ARG' inside figures_dir
            mean_ARG_dir = os.path.join(figures_dir, "mean_ARG")
            os.makedirs(mean_ARG_dir, exist_ok=True)
                # Generate mean ARG based on different groups
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Geographic Location")
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Collection Date")
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Continent")
            generate_mean_arg_lollipop(tidy_file, mean_ARG_dir, fig_format, group_by="Subcontinent")


            # Generate barplots for resistance compariosn
            generate_resistance_barplot(tidy_file, figures_dir, fig_format)

            # Create subdirectory 'heatmap' inside figures_dir
            heatmap_dir = os.path.join(figures_dir, "heatmap")
            os.makedirs(heatmap_dir, exist_ok=True)
                # Generate country comparison heatmap
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Geographic Location", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Collection Date", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Continent", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            generate_comparison_heatmap(tidy_file, heatmap_dir, fig_format, group_col="Subcontinent", resistance_col="RESISTANCE", genep_threshold=genep, nseq_threshold=nseq)
            

            # Plotly plots
            # Generate interactive heatmap with multiple grouping options
            generate_comparison_heatmap_plotly(tidy_file, figures_dir)
            # Generate geographic resistance map
            generate_geographic_resistance_map_plotly(tidy_file, figures_dir)
            # Generate mean ARG resistance analysis box plot
            mean_Arg_resistance_analysis_plotly(tidy_file, figures_dir)
            # Generate interactive lollipop plot with dropdown for group selection
            generate_mean_arg_lollipop_plotly(tidy_file, figures_dir)


            # Generate correlation scatterplot analysis
            print("Generating Geographic Location analysis...")
            correlation_scatterplot_analysis(tidy_file, figures_dir, group_col="Geographic Location", min_samples_per_group=min_samples_per_group)
            print("Generating Continent analysis...")
            correlation_scatterplot_analysis(tidy_file, figures_dir, group_col="Continent", min_samples_per_group=min_samples_per_group)
            print("Generating Subcontinent analysis...")
            correlation_scatterplot_analysis(tidy_file, figures_dir, group_col="Subcontinent", min_samples_per_group=min_samples_per_group)
             # Combine the three CSV files
            print("Combining correlation summary CSV files...")
            combined_correlation_analysis(figures_dir)


            # Keeping all the html files in html_dir inside figures_dir
            html_dir = os.path.join(figures_dir, "html_files")
            os.makedirs(html_dir, exist_ok=True)
            # Move all HTML files to html_dir
            for file in os.listdir(figures_dir):
                if file.endswith(".html"):
                    src_path = os.path.join(figures_dir, file)
                    dest_path = os.path.join(html_dir, file)
                    os.rename(src_path, dest_path)
                    logging.info(f"Moved {file} to {html_dir}")
            # Generate index.html for easy navigation
            generate_index_html(html_dir, figures_dir)

            # Moves every .csv file in figures_dir into that subdirectory.
            stat_analysis_dir = os.path.join(figures_dir, "Stat_analysis")
            os.makedirs(stat_analysis_dir, exist_ok=True)

            # Move all .csv files from figures_dir to figures_dir/Stat_analysis
            for file in os.listdir(figures_dir):
                if file.endswith(".csv"):
                    src = os.path.join(figures_dir, file)
                    dst = os.path.join(stat_analysis_dir, file)
                    shutil.move(src, dst)

            report_options = {
                "fig_format": fig_format,
                "nseq": nseq,
                "genep": genep,
                "min_identity": min_identity,
                "drop_unmatched_accessions": drop_unmatched_accessions,
                "min_samples_per_group": min_samples_per_group,
                "core_threshold": core_threshold,
                "rare_threshold": rare_threshold,
                "top_n": top_n,
                "cooccurrence_min_prevalence": cooccurrence_min_prevalence,
                "cooccurrence_top_n": cooccurrence_top_n,
                "run_cross_database": run_cross_database,
                "cross_database_max_features": cross_database_max_features,
                "plot_style": plot_style,
                "label_max_length": label_max_length or "default",
                "run_temporal_trends": run_temporal_trends,
                "vfdb_dir": vfdb_dir or "not provided",
                "plasmidfinder_dir": plasmidfinder_dir or "not provided",
                "mobileelementfinder_dir": mobileelementfinder_dir or "not provided",
                "isfinder_dir": isfinder_dir or "not provided",
                "integronfinder_dir": integronfinder_dir or "not provided",
                "iceberg_dir": iceberg_dir or "not provided",
                "mlst_dir": mlst_dir or "not provided",
                "defensefinder_dir": defensefinder_dir or "not provided",
                "prophage_dir": prophage_dir or "not provided",
                "sample_map": sample_map_path or "not provided",
                "sequence_dir": sequence_dir or "not provided",
                "run_abricate": run_abricate,
                "abricate_dbs": ",".join(abricate_dbs or []) if abricate_dbs else "not provided",
                "abricate_summary_metric": abricate_summary_metric,
                "run_mobileelementfinder": run_mobileelementfinder_tool,
                "mobileelementfinder_bin": mobileelementfinder_bin,
                "mobileelementfinder_threads": mobileelementfinder_threads,
                "run_integronfinder": run_integronfinder_tool,
                "integronfinder_bin": integronfinder_bin,
                "integronfinder_threads": integronfinder_threads,
                "run_mlst": run_mlst_tool,
                "mlst_bin": mlst_bin,
                "run_defensefinder": run_defensefinder_tool,
                "defensefinder_bin": defensefinder_bin,
                "iceberg_table_dir": iceberg_table_dir or "not provided",
                "tool_manifest_json": tool_manifest.get("json", "not available"),
                "tool_manifest_csv": tool_manifest.get("csv", "not available"),
                "metadata_completeness": metadata_report_outputs.get("metadata_completeness", "not available"),
                "metadata_group_sample_sizes": metadata_report_outputs.get("metadata_group_sample_sizes", "not available"),
                "metadata_bias_warning": metadata_report_outputs.get("metadata_bias_warning", "not available"),
            }
            input_files = {
                "ncbi_clean": ncbi_clean_path,
                "abricate_summary": abricate_summary_file,
                "abricate_results": expected_results_file or "not available",
            }
            citation_outputs = write_citation_outputs(
                output_dir,
                options=report_options,
                feature_outputs=optional_feature_outputs,
                input_files=input_files,
                panr2_version=PANR2_VERSION,
            )
            write_report(
                output_dir,
                base_name,
                ncbi_output_dir=ncbi_output_dir,
                options=report_options,
                panr2_version=PANR2_VERSION,
                feature_outputs=optional_feature_outputs,
                input_files=input_files,
                cross_database_outputs=cross_database_outputs,
                citation_outputs=citation_outputs,
                temporal_outputs=temporal_outputs,
                panresistome_context_outputs=panresistome_context_outputs,
            )

        except Exception as e:
            logging.error(f"Error processing {abricate_summary_file}: {e}")
    
    logging.info("panr run successfully.")


def run_cli(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "run-all":
        return run_cli(argv[1:] + ["--run-all"])
    if argv and argv[0] in {"doctor", "setup", "setup-db", "install-info", "citations", "validate-demo"}:
        return _run_subcommand(argv)

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process NCBI and Abricate data.")
    parser.add_argument("--doctor", action="store_true", help="Check PanR2 Python dependencies and optional external annotation tools, then exit.")
    parser.add_argument("--json", action="store_true", help="With --doctor, write machine-readable JSON output.")
    parser.add_argument("--fix", action="store_true", help="With --doctor, run safe fixes when possible, currently ABRicate database setup if ABRicate is installed.")
    parser.add_argument("--install-info", action="store_true", help="Print PanR2, Python, system, tool, and database readiness information, then exit.")
    parser.add_argument("--setup-db", action="store_true", help="Run ABRicate database setup, then exit.")
    parser.add_argument("--check-only", action="store_true", help="With --setup-db, only report database visibility.")
    parser.add_argument("--citations", action="store_true", help="Write citation/software-version files for --output-dir, then exit.")
    parser.add_argument("--ncbi-dir", help="Directory containing ncbi_clean.csv.")
    parser.add_argument("--abricate-dir", help="Directory containing Abricate summary CSV or TAB files. Required unless --run-abricate is used with the ncbi database.")
    parser.add_argument("--output-dir", help="Base output directory.")
    parser.add_argument("--sequence-dir", help="Directory containing assembly FASTA files used by integrated tool runners.")
    parser.add_argument("--run-all", action="store_true", help="Run all currently integrated annotation runners from --sequence-dir, then perform all PanR2 analyses.")
    parser.add_argument("--run-abricate", action="store_true", help="Run ABRicate internally before PanR2 analysis.")
    parser.add_argument("--abricate-dbs", default="ncbi", help="Comma-separated ABRicate databases to run when --run-abricate is used, for example ncbi,vfdb,plasmidfinder.")
    parser.add_argument("--abricate-bin", default="abricate", help="ABRicate executable name or path.")
    parser.add_argument("--abricate-summary-metric", default="identity", choices=["identity", "coverage"], help="Metric used in generated ABRicate summary matrices.")
    parser.add_argument("--run-mobileelementfinder", action="store_true", help="Run MobileElementFinder internally before PanR2 feature analysis.")
    parser.add_argument("--mobileelementfinder-bin", default="mefinder", help="MobileElementFinder executable name or path.")
    parser.add_argument("--mobileelementfinder-threads", type=int, default=1, help="Threads passed to MobileElementFinder.")
    parser.add_argument("--run-integronfinder", action="store_true", help="Run IntegronFinder internally before PanR2 feature analysis.")
    parser.add_argument("--integronfinder-bin", default="integron_finder", help="IntegronFinder executable name or path.")
    parser.add_argument("--integronfinder-threads", type=int, default=1, help="CPU threads passed to IntegronFinder.")
    parser.add_argument("--run-mlst", action="store_true", help="Run mlst internally before PanR2 typing analysis.")
    parser.add_argument("--mlst-bin", default="mlst", help="MLST executable name or path.")
    parser.add_argument("--run-defensefinder", action="store_true", help="Run DefenseFinder internally before PanR2 feature analysis.")
    parser.add_argument("--defensefinder-bin", default="defense-finder", help="DefenseFinder executable name or path.")
    parser.add_argument("--iceberg-table-dir", help="Directory containing ICE/IME/CIME CSV/TSV/TAB tables to convert into PanR2 ICEberg analysis inputs.")
    parser.add_argument("--force-tool-run", action="store_true", help="Re-run integrated tools even when result files already exist.")
    parser.add_argument("--genep", type=float, default=10.0, help="Minimum %% gene presence to include in heatmap.")
    parser.add_argument("--nseq", type=int, default=1, help="Minimum number of sequences required per group in heatmaps.")
    parser.add_argument("--format", default="tiff", choices=["tiff", "svg", "png", "pdf"], help="Output format for figures (tiff, svg, png, pdf).")
    parser.add_argument("--min-identity", type=float, default=0.0, help="Minimum ABRicate identity percentage to treat a gene call as present. Values below this are set to absent.")
    parser.add_argument("--drop-unmatched-accessions", action="store_true", help="Drop NCBI assemblies that do not have a matching ABRicate summary row after merge.")
    parser.add_argument("--min-samples-per-group", type=int, default=5, help="Minimum samples per group required for correlation analyses.")
    parser.add_argument("--core-threshold", type=float, default=95.0, help="Prevalence percentage used to classify core ARGs.")
    parser.add_argument("--rare-threshold", type=float, default=5.0, help="Prevalence percentage used to classify rare ARGs.")
    parser.add_argument("--top-n", type=int, default=25, help="Number of top genes/classes to include in compact summary plots.")
    parser.add_argument("--cooccurrence-min-prevalence", type=float, default=0.0, help="Minimum prevalence percentage for genes/classes included in co-occurrence matrices.")
    parser.add_argument("--cooccurrence-top-n", type=int, default=25, help="Number of top genes/classes or pairs to include in co-occurrence plots and pair tables.")
    parser.add_argument("--plot-style", default="publication", choices=["publication", "dashboard", "compact"], help="Plot readability preset for integrated figures.")
    parser.add_argument("--label-max-length", type=int, help="Maximum displayed feature-label length in crowded integrated figures.")
    parser.add_argument("--no-cross-database", action="store_true", help="Disable integrated cross-database association outputs.")
    parser.add_argument("--cross-database-max-features", type=int, default=300, help="Maximum most-prevalent features used for pairwise cross-database statistics; use 0 for no limit.")
    parser.add_argument("--no-temporal-trends", action="store_true", help="Disable advanced temporal trend outputs.")
    parser.add_argument("--vfdb-dir", help="Optional directory containing ABRicate VFDB summary/results files.")
    parser.add_argument("--plasmidfinder-dir", help="Optional directory containing ABRicate PlasmidFinder summary/results files.")
    parser.add_argument("--mobileelementfinder-dir", help="Optional directory containing ABRicate MobileElementFinder summary/results files.")
    parser.add_argument("--isfinder-dir", help="Optional directory containing ABRicate ISfinder summary/results files.")
    parser.add_argument("--integronfinder-dir", help="Optional directory containing IntegronFinder or ABRicate-style integron summary/results files.")
    parser.add_argument("--iceberg-dir", help="Optional directory containing ABRicate ICEberg summary/results files.")
    parser.add_argument("--mlst-dir", help="Optional directory containing mlst TSV/CSV output.")
    parser.add_argument("--defensefinder-dir", help="Optional directory containing DefenseFinder tables or PanR2-compatible DefenseFinder summary/results files.")
    parser.add_argument("--prophage-dir", help="Optional directory containing prophage/viral-region tables or PanR2-compatible prophage summary/results files.")
    parser.add_argument("--sample-map", help="Optional CSV/TSV mapping sample_id values to Assembly Accession for filenames or tool outputs without GCF/GCA accessions.")
    parser.add_argument('--version', action='version', version=f'PanR2 {PANR2_VERSION}')

    args = parser.parse_args(argv)
    if args.doctor or args.install_info:
        return run_doctor(
            abricate_bin=args.abricate_bin,
            mobileelementfinder_bin=args.mobileelementfinder_bin,
            integronfinder_bin=args.integronfinder_bin,
            mlst_bin=args.mlst_bin,
            defensefinder_bin=args.defensefinder_bin,
            json_output=args.json,
            fix=args.fix,
            include_system=args.install_info,
        )
    if args.setup_db:
        dbs = [db.strip() for db in args.abricate_dbs.split(",") if db.strip()]
        return setup_abricate_databases(abricate_bin=args.abricate_bin, dbs=dbs, check_only=args.check_only, json_output=args.json)
    if args.citations:
        if not args.output_dir:
            parser.error("--output-dir is required when --citations is used.")
        write_citation_outputs(args.output_dir, panr2_version=PANR2_VERSION)
        return 0

    if not args.ncbi_dir:
        parser.error("--ncbi-dir is required unless --doctor, --install-info, --setup-db, or --citations is used.")
    if not args.output_dir:
        parser.error("--output-dir is required unless --doctor, --install-info, --setup-db, or --citations is used.")

    if args.run_all:
        args.run_abricate = True
        args.run_mobileelementfinder = True
        args.run_integronfinder = True
        args.run_mlst = True
        args.run_defensefinder = True
        if args.abricate_dbs == "ncbi":
            args.abricate_dbs = "ncbi,vfdb,plasmidfinder,isfinder"

    abricate_dbs = [db.strip() for db in args.abricate_dbs.split(",") if db.strip()]
    
    # Run the main function
    main(
        args.ncbi_dir,
        args.abricate_dir,
        args.output_dir,
        args.format,
        args.nseq,
        args.genep,
        min_identity=args.min_identity,
        drop_unmatched_accessions=args.drop_unmatched_accessions,
        min_samples_per_group=args.min_samples_per_group,
        core_threshold=args.core_threshold,
        rare_threshold=args.rare_threshold,
        top_n=args.top_n,
        cooccurrence_min_prevalence=args.cooccurrence_min_prevalence,
        cooccurrence_top_n=args.cooccurrence_top_n,
        vfdb_dir=args.vfdb_dir,
        plasmidfinder_dir=args.plasmidfinder_dir,
        mobileelementfinder_dir=args.mobileelementfinder_dir,
        isfinder_dir=args.isfinder_dir,
        integronfinder_dir=args.integronfinder_dir,
        iceberg_dir=args.iceberg_dir,
        mlst_dir=args.mlst_dir,
        defensefinder_dir=args.defensefinder_dir,
        prophage_dir=args.prophage_dir,
        sample_map_path=args.sample_map,
        sequence_dir=args.sequence_dir,
        run_abricate=args.run_abricate,
        abricate_dbs=abricate_dbs,
        abricate_bin=args.abricate_bin,
        abricate_summary_metric=args.abricate_summary_metric,
        run_mobileelementfinder_tool=args.run_mobileelementfinder,
        mobileelementfinder_bin=args.mobileelementfinder_bin,
        mobileelementfinder_threads=args.mobileelementfinder_threads,
        run_integronfinder_tool=args.run_integronfinder,
        integronfinder_bin=args.integronfinder_bin,
        integronfinder_threads=args.integronfinder_threads,
        run_mlst_tool=args.run_mlst,
        mlst_bin=args.mlst_bin,
        run_defensefinder_tool=args.run_defensefinder,
        defensefinder_bin=args.defensefinder_bin,
        iceberg_table_dir=args.iceberg_table_dir,
        force_tool_run=args.force_tool_run,
        run_cross_database=not args.no_cross_database,
        cross_database_max_features=args.cross_database_max_features,
        plot_style=args.plot_style,
        label_max_length=args.label_max_length,
        run_temporal_trends=not args.no_temporal_trends,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli())
