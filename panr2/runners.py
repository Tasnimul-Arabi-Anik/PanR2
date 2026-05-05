import csv
import json
import logging
import os
import shutil
import subprocess
from datetime import datetime

import pandas as pd

from panr2.io import resolve_sample_accessions, write_sample_map_qc

SEQUENCE_EXTENSIONS = (".fa", ".fna", ".fasta", ".fas")


def find_sequence_files(sequence_dir):
    """Return FASTA-like assembly files from a sequence directory."""
    if not sequence_dir:
        return []
    if not os.path.isdir(sequence_dir):
        raise FileNotFoundError(f"Sequence directory not found: {sequence_dir}")
    files = []
    for name in sorted(os.listdir(sequence_dir)):
        lower = name.lower()
        if lower.endswith(SEQUENCE_EXTENSIONS) or lower.endswith(tuple(ext + ".gz" for ext in SEQUENCE_EXTENSIONS)):
            files.append(os.path.join(sequence_dir, name))
    return files


def _run_command(command, stdout_path=None):
    logging.info("Running: %s", " ".join(command))
    if stdout_path:
        with open(stdout_path, "w") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
    else:
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or getattr(completed, "stdout", "").strip()
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}\n{detail}")
    return completed.stdout if not stdout_path else ""


def _capture_command(command):
    try:
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        output = (completed.stdout or completed.stderr or "").strip()
        return output if completed.returncode == 0 else f"unavailable: {output}"
    except Exception as exc:
        return f"unavailable: {exc}"


def _parse_abricate_list(output):
    rows = {}
    reader = csv.DictReader(output.splitlines(), delimiter="\t")
    for row in reader:
        db = row.get("DATABASE")
        if db:
            rows[db] = row
    return rows


def write_tool_manifest(output_dir, manifest):
    """Write tool run metadata for reproducibility."""
    qc_dir = os.path.join(output_dir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    json_path = os.path.join(qc_dir, "panr2_tool_manifest.json")
    if os.path.exists(json_path):
        with open(json_path) as handle:
            existing = json.load(handle)
        existing.setdefault("tools", []).extend(manifest.get("tools", []))
        existing["generated_at"] = manifest.get("generated_at", existing.get("generated_at", ""))
        existing["sequence_dir"] = manifest.get("sequence_dir", existing.get("sequence_dir", ""))
        existing["sequence_count"] = manifest.get("sequence_count", existing.get("sequence_count", ""))
        manifest = existing
    with open(json_path, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    rows = []
    for tool in manifest.get("tools", []):
        for run in tool.get("runs", []):
            rows.append({
                "tool": tool.get("name", ""),
                "version": tool.get("version", ""),
                "database": run.get("database", ""),
                "database_sequences": run.get("database_sequences", ""),
                "database_date": run.get("database_date", ""),
                "results": run.get("results", ""),
                "summary": run.get("summary", ""),
                "status": run.get("status", ""),
            })
    csv_path = os.path.join(qc_dir, "panr2_tool_manifest.csv")
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "tool", "version", "database", "database_sequences", "database_date",
            "results", "summary", "status",
        ])
        writer.writeheader()
        writer.writerows(rows)
    return {"json": json_path, "csv": csv_path}


def _strip_gzip_suffix(path):
    base = os.path.basename(path)
    return base[:-3] if base.endswith(".gz") else base


def _sample_prefix(path):
    name = _strip_gzip_suffix(path)
    for ext in SEQUENCE_EXTENSIONS:
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return os.path.splitext(name)[0]


def _clean_feature_id(value):
    value = str(value or "").strip()
    return value if value else "unknown_feature"


def _float_or_default(value, default=0.0):
    try:
        if value is None or str(value).strip() == "":
            return default
        return float(str(value).strip().rstrip("%"))
    except Exception:
        return default


def _first_value(row, candidates, default=""):
    lower_map = {str(key).lower().replace(" ", "_").replace("-", "_"): key for key in row}
    for candidate in candidates:
        key = lower_map.get(candidate.lower().replace(" ", "_").replace("-", "_"))
        if key is not None and str(row.get(key, "")).strip():
            return row.get(key, "")
    return default


def _read_mobileelementfinder_csv(path):
    with open(path, newline="") as handle:
        lines = [line for line in handle if not line.startswith("#") and line.strip()]
    if not lines:
        return []
    return list(csv.DictReader(lines))


def convert_mobileelementfinder_outputs(sequence_files, raw_csv_paths, output_dir):
    """Convert MobileElementFinder CSV outputs into ABRicate-style files."""
    feature_dir = os.path.join(output_dir, "tool_results", "mobileelementfinder", "panr2_inputs")
    os.makedirs(feature_dir, exist_ok=True)
    results_path = os.path.join(feature_dir, "mobileelementfinder_results.tab")
    summary_path = os.path.join(feature_dir, "mobileelementfinder_summary.tab")

    results_rows = []
    feature_ids = []
    by_sample = {}
    file_lookup = {_sample_prefix(path): path for path in sequence_files}

    for csv_path in raw_csv_paths:
        prefix = os.path.basename(csv_path).replace(".csv", "")
        source_file = file_lookup.get(prefix)
        if not source_file:
            for sample_prefix, path in file_lookup.items():
                if sample_prefix in prefix or prefix in sample_prefix:
                    source_file = path
                    break
        source_file = source_file or prefix
        sample_key = _sample_prefix(source_file)
        by_sample.setdefault(sample_key, {"file": source_file, "features": {}})

        for row in _read_mobileelementfinder_csv(csv_path):
            feature_id = _clean_feature_id(_first_value(row, [
                "mge_id", "mge", "mobile_element", "element", "template", "gene", "name", "id"
            ]))
            identity = _float_or_default(_first_value(row, [
                "identity", "perc_identity", "percent_identity", "%identity", "%_identity"
            ], 100.0), 100.0)
            coverage = _float_or_default(_first_value(row, [
                "coverage", "perc_coverage", "percent_coverage", "%coverage", "%_coverage"
            ], 100.0), 100.0)
            contig = _first_value(row, ["contig", "sequence", "reference", "qseqid"], "contig")
            start = _first_value(row, ["start", "query_start", "qstart"], "0")
            end = _first_value(row, ["end", "query_end", "qend"], "0")
            accession = _first_value(row, ["accession", "reference_accession", "template_id"], feature_id)
            product = _first_value(row, ["type", "family", "description", "product", "element_type"], feature_id)
            feature_ids.append(feature_id)
            by_sample[sample_key]["features"][feature_id] = max(identity, by_sample[sample_key]["features"].get(feature_id, 0.0))
            results_rows.append({
                "#FILE": source_file,
                "SEQUENCE": contig,
                "START": start,
                "END": end,
                "GENE": feature_id,
                "COVERAGE": f"{coverage:.2f}",
                "%COVERAGE": f"{coverage:.2f}",
                "%IDENTITY": f"{identity:.2f}",
                "DATABASE": "mobileelementfinder",
                "ACCESSION": accession,
                "PRODUCT": product,
            })

    feature_ids = sorted(set(feature_ids), key=str.lower)
    with open(results_path, "w", newline="") as handle:
        fieldnames = ["#FILE", "SEQUENCE", "START", "END", "GENE", "COVERAGE", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results_rows)

    with open(summary_path, "w", newline="") as handle:
        fieldnames = ["#FILE", "NUM_FOUND"] + feature_ids
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sample_key in sorted(file_lookup):
            features = by_sample.get(sample_key, {"file": file_lookup[sample_key], "features": {}})
            row = {"#FILE": features["file"], "NUM_FOUND": sum(1 for value in features["features"].values() if value > 0)}
            for feature_id in feature_ids:
                row[feature_id] = features["features"].get(feature_id, 0)
            writer.writerow(row)

    return {"results": results_path, "summary": summary_path, "feature_dir": feature_dir}


def _find_mobileelementfinder_csv(raw_dir, prefix):
    exact = os.path.join(raw_dir, f"{prefix}.csv")
    if os.path.exists(exact):
        return exact
    candidates = []
    for name in os.listdir(raw_dir):
        if name.endswith(".csv") and prefix in name:
            candidates.append(os.path.join(raw_dir, name))
    return sorted(candidates)[0] if candidates else None


def run_mobileelementfinder(
    sequence_dir,
    output_dir,
    mefinder_bin="mefinder",
    threads=1,
    force=False,
):
    """Run MobileElementFinder per assembly and create PanR2-compatible inputs."""
    executable = shutil.which(mefinder_bin)
    if not executable:
        raise FileNotFoundError(
            f"MobileElementFinder executable not found: {mefinder_bin}. "
            "Run `panr doctor` to inspect the environment, or install integrated dependencies with environment.yml/Docker."
        )

    sequence_files = find_sequence_files(sequence_dir)
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")

    raw_dir = os.path.join(output_dir, "tool_results", "mobileelementfinder", "raw")
    os.makedirs(raw_dir, exist_ok=True)
    version = _capture_command([executable, "--version"]).splitlines()[0]
    raw_csv_paths = []

    for sequence_file in sequence_files:
        prefix = _sample_prefix(sequence_file)
        output_prefix = os.path.join(raw_dir, prefix)
        expected_csv = _find_mobileelementfinder_csv(raw_dir, prefix)
        if force or not expected_csv:
            _run_command([executable, "find", "--contig", sequence_file, "--threads", str(threads), output_prefix])
            expected_csv = _find_mobileelementfinder_csv(raw_dir, prefix)
        if expected_csv:
            raw_csv_paths.append(expected_csv)
        else:
            logging.warning("MobileElementFinder did not create a CSV output for %s", sequence_file)

    converted = convert_mobileelementfinder_outputs(sequence_files, raw_csv_paths, output_dir)
    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sequence_dir": sequence_dir,
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "mobileelementfinder",
            "executable": executable,
            "version": version,
            "runs": [{
                "database": "mobileelementfinder",
                "database_sequences": "",
                "database_date": "",
                "results": converted["results"],
                "summary": converted["summary"],
                "status": "completed",
            }],
        }],
    }
    manifest_paths = write_tool_manifest(output_dir, manifest)
    logging.info("MobileElementFinder tool manifest saved to %s", manifest_paths["json"])
    return {"feature_dir": converted["feature_dir"], "manifest": manifest_paths, "raw_csv": raw_csv_paths}


def _read_integronfinder_table(path):
    with open(path, newline="") as handle:
        lines = [line for line in handle if not line.startswith("#") and line.strip()]
    if not lines:
        return []
    sample = lines[0]
    delimiter = "\t" if "\t" in sample else ","
    return list(csv.DictReader(lines, delimiter=delimiter))


def _find_integronfinder_table(raw_dir, prefix):
    candidates = []
    for root, _, files in os.walk(raw_dir):
        for name in files:
            lower = name.lower()
            if prefix in name and (lower.endswith(".integrons") or lower.endswith(".summary") or lower.endswith(".tsv") or lower.endswith(".tab") or lower.endswith(".csv")):
                candidates.append(os.path.join(root, name))
    if not candidates:
        return None

    def priority(path):
        lower = path.lower()
        if lower.endswith(".integrons"):
            return 0
        if lower.endswith((".tsv", ".tab", ".csv")):
            return 1
        if lower.endswith(".summary"):
            return 2
        return 3

    selected = sorted(candidates, key=lambda path: (priority(path), path))[0]
    logging.info("Selected IntegronFinder table for %s: %s", prefix, selected)
    if priority(selected) == 2:
        logging.warning("Using IntegronFinder summary file for %s because no detailed .integrons/TSV/TAB/CSV table was found.", prefix)
    return selected


def _integron_feature_id(row):
    feature_type = _first_value(row, ["type", "element", "annotation", "model", "name"], "integron")
    raw_id = _first_value(row, ["id", "ID", "integron_id", "attc_id", "protein_id", "gene"], "")
    feature_type = _clean_feature_id(feature_type)
    if raw_id:
        raw_id = _clean_feature_id(raw_id)
        if raw_id.lower().startswith(feature_type.lower()):
            return raw_id
        return f"{feature_type}_{raw_id}"
    return feature_type


def convert_integronfinder_outputs(sequence_files, raw_table_paths, output_dir):
    """Convert IntegronFinder tabular outputs into ABRicate-style files."""
    feature_dir = os.path.join(output_dir, "tool_results", "integronfinder", "panr2_inputs")
    os.makedirs(feature_dir, exist_ok=True)
    results_path = os.path.join(feature_dir, "integronfinder_results.tab")
    summary_path = os.path.join(feature_dir, "integronfinder_summary.tab")

    results_rows = []
    feature_ids = []
    by_sample = {}
    file_lookup = {_sample_prefix(path): path for path in sequence_files}

    for table_path in raw_table_paths:
        prefix = os.path.basename(table_path)
        source_file = None
        for sample_prefix, path in file_lookup.items():
            if sample_prefix in prefix or sample_prefix in table_path:
                source_file = path
                break
        source_file = source_file or prefix
        sample_key = _sample_prefix(source_file)
        by_sample.setdefault(sample_key, {"file": source_file, "features": {}})

        for row in _read_integronfinder_table(table_path):
            feature_id = _integron_feature_id(row)
            contig = _first_value(row, ["contig", "replicon", "sequence", "ID_replicon"], "contig")
            start = _first_value(row, ["start", "pos_beg", "begin", "left"], "0")
            end = _first_value(row, ["end", "pos_end", "stop", "right"], "0")
            product = _first_value(row, ["type", "element", "annotation", "model"], feature_id)
            accession = _first_value(row, ["id", "integron_id", "protein_id", "gene"], feature_id)
            identity = _float_or_default(_first_value(row, ["identity", "score", "%identity"], 100.0), 100.0)
            coverage = _float_or_default(_first_value(row, ["coverage", "%coverage"], 100.0), 100.0)
            feature_ids.append(feature_id)
            by_sample[sample_key]["features"][feature_id] = max(identity, by_sample[sample_key]["features"].get(feature_id, 0.0))
            results_rows.append({
                "#FILE": source_file,
                "SEQUENCE": contig,
                "START": start,
                "END": end,
                "GENE": feature_id,
                "COVERAGE": f"{coverage:.2f}",
                "%COVERAGE": f"{coverage:.2f}",
                "%IDENTITY": f"{identity:.2f}",
                "DATABASE": "integronfinder",
                "ACCESSION": accession,
                "PRODUCT": product,
            })

    feature_ids = sorted(set(feature_ids), key=str.lower)
    with open(results_path, "w", newline="") as handle:
        fieldnames = ["#FILE", "SEQUENCE", "START", "END", "GENE", "COVERAGE", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results_rows)

    with open(summary_path, "w", newline="") as handle:
        fieldnames = ["#FILE", "NUM_FOUND"] + feature_ids
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sample_key in sorted(file_lookup):
            features = by_sample.get(sample_key, {"file": file_lookup[sample_key], "features": {}})
            row = {"#FILE": features["file"], "NUM_FOUND": sum(1 for value in features["features"].values() if value > 0)}
            for feature_id in feature_ids:
                row[feature_id] = features["features"].get(feature_id, 0)
            writer.writerow(row)

    return {"results": results_path, "summary": summary_path, "feature_dir": feature_dir}


def run_integronfinder(
    sequence_dir,
    output_dir,
    integronfinder_bin="integron_finder",
    cpu=1,
    force=False,
):
    """Run IntegronFinder per assembly and create PanR2-compatible inputs."""
    executable = shutil.which(integronfinder_bin)
    if not executable:
        raise FileNotFoundError(
            f"IntegronFinder executable not found: {integronfinder_bin}. "
            "Run `panr doctor` to inspect the environment, or install integrated dependencies with environment.yml/Docker."
        )

    sequence_files = find_sequence_files(sequence_dir)
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")

    raw_dir = os.path.join(output_dir, "tool_results", "integronfinder", "raw")
    os.makedirs(raw_dir, exist_ok=True)
    version = _capture_command([executable, "--version"]).splitlines()[0]
    raw_table_paths = []

    for sequence_file in sequence_files:
        prefix = _sample_prefix(sequence_file)
        sample_out_dir = os.path.join(raw_dir, prefix)
        os.makedirs(sample_out_dir, exist_ok=True)
        expected_table = _find_integronfinder_table(raw_dir, prefix)
        if force or not expected_table:
            _run_command([executable, sequence_file, "--outdir", sample_out_dir, "--cpu", str(cpu)])
            expected_table = _find_integronfinder_table(raw_dir, prefix)
        if expected_table:
            raw_table_paths.append(expected_table)
        else:
            logging.warning("IntegronFinder did not create a tabular output for %s", sequence_file)

    converted = convert_integronfinder_outputs(sequence_files, raw_table_paths, output_dir)
    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sequence_dir": sequence_dir,
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "integronfinder",
            "executable": executable,
            "version": version,
            "runs": [{
                "database": "integronfinder",
                "database_sequences": "",
                "database_date": "",
                "results": converted["results"],
                "summary": converted["summary"],
                "status": "completed",
            }],
        }],
    }
    manifest_paths = write_tool_manifest(output_dir, manifest)
    logging.info("IntegronFinder tool manifest saved to %s", manifest_paths["json"])
    return {"feature_dir": converted["feature_dir"], "manifest": manifest_paths, "raw_tables": raw_table_paths}


def _read_delimited_table(path):
    with open(path, newline="") as handle:
        lines = [line for line in handle if not line.startswith("#") and line.strip()]
    if not lines:
        return []
    sample = lines[0]
    delimiter = "\t" if "\t" in sample else ","
    return list(csv.DictReader(lines, delimiter=delimiter))


def _find_table_files(table_dir):
    if not os.path.isdir(table_dir):
        raise FileNotFoundError(f"Table directory not found: {table_dir}")
    paths = []
    for root, _, files in os.walk(table_dir):
        for name in sorted(files):
            lower = name.lower()
            if lower.endswith((".csv", ".tsv", ".tab")):
                paths.append(os.path.join(root, name))
    return paths


def convert_iceberg_tables(table_dir, output_dir, sample_map=None):
    """Convert user-provided ICE/IME/CIME tables into ABRicate-style files."""
    table_paths = _find_table_files(table_dir)
    if not table_paths:
        raise FileNotFoundError(f"No CSV/TSV/TAB ICEberg tables found in {table_dir}")

    feature_dir = os.path.join(output_dir, "tool_results", "iceberg", "panr2_inputs")
    os.makedirs(feature_dir, exist_ok=True)
    results_path = os.path.join(feature_dir, "iceberg_results.tab")
    summary_path = os.path.join(feature_dir, "iceberg_summary.tab")

    results_rows = []
    feature_ids = []
    by_sample = {}
    original_samples = []
    resolved_samples = []

    for table_path in table_paths:
        for row in _read_delimited_table(table_path):
            sample = _first_value(row, [
                "file", "#file", "assembly_file", "assembly", "assembly_accession",
                "genome", "genome_id", "sample", "sample_id", "isolate"
            ], os.path.basename(table_path))
            resolved_sample = resolve_sample_accessions(pd.Series([str(sample)]), sample_map=sample_map).iloc[0]
            sample_value = str(resolved_sample or sample)
            original_samples.append(str(sample))
            resolved_samples.append(sample_value)
            feature_id = _clean_feature_id(_first_value(row, [
                "ice_id", "element_id", "element", "mge_id", "name", "id", "accession", "gene"
            ], "iceberg_feature"))
            element_type = _first_value(row, ["type", "element_type", "category", "class"], "ICE/IME/CIME")
            identity = _float_or_default(_first_value(row, [
                "identity", "perc_identity", "percent_identity", "%identity", "%_identity", "score"
            ], 100.0), 100.0)
            coverage = _float_or_default(_first_value(row, [
                "coverage", "perc_coverage", "percent_coverage", "%coverage", "%_coverage"
            ], 100.0), 100.0)
            contig = _first_value(row, ["contig", "sequence", "replicon", "chromosome"], "contig")
            start = _first_value(row, ["start", "pos_beg", "begin", "left"], "0")
            end = _first_value(row, ["end", "pos_end", "stop", "right"], "0")
            accession = _first_value(row, ["accession", "ice_accession", "reference", "reference_id"], feature_id)
            product = _first_value(row, ["description", "product", "annotation", "element_type", "type"], element_type)

            sample_key = _sample_prefix(sample_value)
            by_sample.setdefault(sample_key, {"file": sample_value, "features": {}})
            feature_ids.append(feature_id)
            by_sample[sample_key]["features"][feature_id] = max(identity, by_sample[sample_key]["features"].get(feature_id, 0.0))
            results_rows.append({
                "#FILE": sample_value,
                "SEQUENCE": contig,
                "START": start,
                "END": end,
                "GENE": feature_id,
                "COVERAGE": f"{coverage:.2f}",
                "%COVERAGE": f"{coverage:.2f}",
                "%IDENTITY": f"{identity:.2f}",
                "DATABASE": "iceberg",
                "ACCESSION": accession,
                "PRODUCT": product,
            })

    feature_ids = sorted(set(feature_ids), key=str.lower)
    if not feature_ids:
        raise ValueError(f"No ICEberg features found in {table_dir}")

    with open(results_path, "w", newline="") as handle:
        fieldnames = ["#FILE", "SEQUENCE", "START", "END", "GENE", "COVERAGE", "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results_rows)

    with open(summary_path, "w", newline="") as handle:
        fieldnames = ["#FILE", "NUM_FOUND"] + feature_ids
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sample_key in sorted(by_sample):
            features = by_sample[sample_key]
            row = {"#FILE": features["file"], "NUM_FOUND": sum(1 for value in features["features"].values() if value > 0)}
            for feature_id in feature_ids:
                row[feature_id] = features["features"].get(feature_id, 0)
            writer.writerow(row)

    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sequence_dir": "",
        "sequence_count": "",
        "tools": [{
            "name": "iceberg_table_converter",
            "executable": "panr2",
            "version": "",
            "runs": [{
                "database": "iceberg",
                "database_sequences": "",
                "database_date": "",
                "results": results_path,
                "summary": summary_path,
                "status": "completed",
            }],
        }],
    }
    manifest_paths = write_tool_manifest(output_dir, manifest)
    write_sample_map_qc(output_dir, "iceberg", pd.Series(original_samples), pd.Series(resolved_samples), sample_map)
    logging.info("ICEberg table converter manifest saved to %s", manifest_paths["json"])
    return {"feature_dir": feature_dir, "manifest": manifest_paths, "source_tables": table_paths}


def run_abricate_databases(
    sequence_dir,
    output_dir,
    databases,
    abricate_bin="abricate",
    summary_metric="identity",
    force=False,
):
    """Run ABRicate for selected databases and return generated output directories."""
    if not databases:
        raise ValueError("At least one ABRicate database must be provided.")
    executable = shutil.which(abricate_bin)
    if not executable:
        raise FileNotFoundError(
            f"ABRicate executable not found: {abricate_bin}. "
            "Run `panr doctor` to inspect the environment, or install integrated dependencies with environment.yml/Docker."
        )

    sequence_files = find_sequence_files(sequence_dir)
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")

    list_output = _capture_command([executable, "--list"])
    available = _parse_abricate_list(list_output)
    if not available:
        raise ValueError(
            "ABRicate did not report any available databases. "
            "Run `panr setup-db` or `abricate --setupdb`, then verify with `panr doctor`."
        )
    missing = [db for db in databases if available and db not in available]
    if missing:
        raise ValueError(
            f"ABRicate database(s) not available: {', '.join(missing)}. "
            "Run `panr setup-db --dbs " + ",".join(databases) + "` or check the installed ABRicate databases with `panr doctor`."
        )

    base_dir = os.path.join(output_dir, "tool_results", "abricate")
    os.makedirs(base_dir, exist_ok=True)
    version = _capture_command([executable, "--version"]).splitlines()[0]
    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sequence_dir": sequence_dir,
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "abricate",
            "executable": executable,
            "version": version,
            "runs": [],
        }],
    }
    output_dirs = {}

    for db in databases:
        db_dir = os.path.join(base_dir, db)
        os.makedirs(db_dir, exist_ok=True)
        results_path = os.path.join(db_dir, f"{db}_results.tab")
        summary_path = os.path.join(db_dir, f"{db}_summary.tab")

        if force or not os.path.exists(results_path):
            _run_command([executable, "--db", db] + sequence_files, stdout_path=results_path)
        else:
            logging.info("Reusing existing ABRicate results: %s", results_path)

        if force or not os.path.exists(summary_path):
            summary_cmd = [executable, "--summary"]
            if summary_metric == "identity":
                summary_cmd.append("--identity")
            summary_cmd.append(results_path)
            _run_command(summary_cmd, stdout_path=summary_path)
        else:
            logging.info("Reusing existing ABRicate summary: %s", summary_path)

        db_info = available.get(db, {})
        manifest["tools"][0]["runs"].append({
            "database": db,
            "database_sequences": db_info.get("SEQUENCES", ""),
            "database_date": db_info.get("DATE", ""),
            "database_type": db_info.get("DBTYPE", ""),
            "results": results_path,
            "summary": summary_path,
            "status": "completed",
        })
        output_dirs[db] = db_dir

    manifest_paths = write_tool_manifest(output_dir, manifest)
    logging.info("ABRicate tool manifest saved to %s", manifest_paths["json"])
    return {"database_dirs": output_dirs, "manifest": manifest_paths}
