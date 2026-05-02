import csv
import json
import logging
import os
import shutil
import subprocess
from datetime import datetime


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
        raise FileNotFoundError(f"ABRicate executable not found: {abricate_bin}")

    sequence_files = find_sequence_files(sequence_dir)
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")

    list_output = _capture_command([executable, "--list"])
    available = _parse_abricate_list(list_output)
    missing = [db for db in databases if available and db not in available]
    if missing:
        raise ValueError(f"ABRicate database(s) not available: {', '.join(missing)}")

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
