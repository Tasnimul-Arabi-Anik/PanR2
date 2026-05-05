import csv
import os
import shutil
import subprocess
from datetime import datetime

from panr2.runners import (
    _clean_feature_id,
    _find_table_files,
    _first_value,
    _float_or_default,
    _read_delimited_table,
    _sample_prefix,
    find_sequence_files,
    write_tool_manifest,
)


def _write_abricate_style(feature_type, database_name, rows, output_dir):
    feature_dir = os.path.join(output_dir, "tool_results", feature_type, "panr2_inputs")
    os.makedirs(feature_dir, exist_ok=True)
    results_path = os.path.join(feature_dir, f"{feature_type}_results.tab")
    summary_path = os.path.join(feature_dir, f"{feature_type}_summary.tab")

    feature_ids = sorted({row["GENE"] for row in rows}, key=str.lower)
    by_sample = {}
    for row in rows:
        sample_key = _sample_prefix(row["#FILE"])
        by_sample.setdefault(sample_key, {"file": row["#FILE"], "features": {}})
        identity = _float_or_default(row["%IDENTITY"], 100.0)
        by_sample[sample_key]["features"][row["GENE"]] = max(
            identity,
            by_sample[sample_key]["features"].get(row["GENE"], 0.0),
        )

    with open(results_path, "w", newline="") as handle:
        fieldnames = [
            "#FILE", "SEQUENCE", "START", "END", "GENE", "COVERAGE",
            "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with open(summary_path, "w", newline="") as handle:
        fieldnames = ["#FILE", "NUM_FOUND"] + feature_ids
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for sample_key in sorted(by_sample):
            sample = by_sample[sample_key]
            row = {
                "#FILE": sample["file"],
                "NUM_FOUND": sum(1 for value in sample["features"].values() if value > 0),
            }
            for feature_id in feature_ids:
                row[feature_id] = sample["features"].get(feature_id, 0)
            writer.writerow(row)

    return {
        "feature_dir": feature_dir,
        "results": results_path,
        "summary": summary_path,
        "database": database_name,
    }


def _convert_generic_tables(table_dir, output_dir, feature_type, database_name, feature_candidates, category_candidates):
    table_paths = _find_table_files(table_dir)
    if not table_paths:
        raise FileNotFoundError(f"No CSV/TSV/TAB tables found in {table_dir}")

    rows = []
    for table_path in table_paths:
        for source_row in _read_delimited_table(table_path):
            sample = _first_value(source_row, [
                "file", "#file", "assembly_file", "assembly", "assembly_accession",
                "genome", "genome_id", "sample", "sample_id", "isolate", "replicon",
            ], os.path.basename(table_path))
            feature_id = _clean_feature_id(_first_value(source_row, feature_candidates, f"{feature_type}_feature"))
            category = _first_value(source_row, category_candidates, feature_id)
            identity = _float_or_default(_first_value(source_row, [
                "identity", "perc_identity", "percent_identity", "%identity", "%_identity", "score", "confidence",
            ], 100.0), 100.0)
            coverage = _float_or_default(_first_value(source_row, [
                "coverage", "perc_coverage", "percent_coverage", "%coverage", "%_coverage",
            ], 100.0), 100.0)
            rows.append({
                "#FILE": str(sample),
                "SEQUENCE": _first_value(source_row, ["contig", "sequence", "replicon", "chromosome"], "contig"),
                "START": _first_value(source_row, ["start", "pos_beg", "begin", "left"], "0"),
                "END": _first_value(source_row, ["end", "pos_end", "stop", "right"], "0"),
                "GENE": feature_id,
                "COVERAGE": f"{coverage:.2f}",
                "%COVERAGE": f"{coverage:.2f}",
                "%IDENTITY": f"{identity:.2f}",
                "DATABASE": database_name,
                "ACCESSION": _first_value(source_row, ["accession", "id", "hit_id", "protein_id", "region_id"], feature_id),
                "PRODUCT": category,
            })

    if not rows:
        raise ValueError(f"No {feature_type} features found in {table_dir}")
    return _write_abricate_style(feature_type, database_name, rows, output_dir)


def convert_defensefinder_tables(table_dir, output_dir):
    """Convert DefenseFinder-style TSV/CSV outputs into PanR2 feature inputs."""
    return _convert_generic_tables(
        table_dir,
        output_dir,
        "defensefinder",
        "defensefinder",
        [
            "system", "system_type", "subtype", "type", "gene_name", "name",
            "hit_id", "sys_id", "protein_in_syst", "feature", "id",
        ],
        ["system_type", "subtype", "type", "system", "product", "description"],
    )


def convert_prophage_tables(table_dir, output_dir):
    """Convert prophage/viral-region TSV/CSV outputs into PanR2 feature inputs."""
    return _convert_generic_tables(
        table_dir,
        output_dir,
        "prophage",
        "prophage",
        [
            "prophage_id", "region_id", "viral_region", "phage", "phage_id",
            "name", "type", "completeness", "id", "feature",
        ],
        ["type", "completeness", "category", "tool", "product", "description"],
    )


def run_defensefinder(sequence_dir, output_dir, defensefinder_bin="defense-finder", force=False):
    """Run DefenseFinder per assembly when installed, then convert tabular outputs."""
    executable = shutil.which(defensefinder_bin)
    if not executable:
        raise FileNotFoundError(
            f"DefenseFinder executable not found: {defensefinder_bin}. Run `panr doctor` or provide --defensefinder-dir with existing tables."
        )
    sequence_files = find_sequence_files(sequence_dir)
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")
    raw_dir = os.path.join(output_dir, "tool_results", "defensefinder", "raw")
    os.makedirs(raw_dir, exist_ok=True)
    for sequence_file in sequence_files:
        prefix = _sample_prefix(sequence_file)
        sample_out = os.path.join(raw_dir, prefix)
        os.makedirs(sample_out, exist_ok=True)
        existing = [name for name in os.listdir(sample_out) if name.lower().endswith((".tsv", ".tab", ".csv"))]
        if existing and not force:
            continue
        completed = subprocess.run([executable, "run", sequence_file, "-o", sample_out], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if completed.returncode != 0:
            raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())
    converted = convert_defensefinder_tables(raw_dir, output_dir)
    version = "available"
    try:
        completed = subprocess.run([executable, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        version = (completed.stdout or completed.stderr or "available").strip().splitlines()[0]
    except Exception:
        pass
    manifest = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sequence_dir": sequence_dir,
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "defensefinder",
            "executable": executable,
            "version": version,
            "runs": [{
                "database": "defensefinder",
                "database_sequences": "",
                "database_date": "",
                "results": converted["results"],
                "summary": converted["summary"],
                "status": "completed",
            }],
        }],
    }
    manifest_paths = write_tool_manifest(output_dir, manifest)
    return {"feature_dir": converted["feature_dir"], "manifest": manifest_paths, "raw_dir": raw_dir}
