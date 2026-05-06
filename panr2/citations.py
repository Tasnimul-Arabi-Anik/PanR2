import csv
import os
import platform
import sys
from importlib import metadata as importlib_metadata

import pandas as pd


CITATION_REGISTRY = {
    "panr2": {
        "name": "PanR2",
        "citation": "PanR2: Panresistome Analysis Tool. DOI: 10.1101/2025.04.08.647722",
        "url": "https://github.com/Tasnimul-Arabi-Anik/PanR2",
    },
    "fetchm": {
        "name": "FetchM",
        "citation": "FetchM metadata standardization workflow.",
        "url": "https://github.com/Tasnimul-Arabi-Anik/FetchM",
    },
    "fetchm2": {
        "name": "FetchM2",
        "citation": "FetchM2 metadata standardization, audit, and sequence-download workflow.",
        "url": "https://github.com/Tasnimul-Arabi-Anik/FetchM2",
    },
    "abricate": {
        "name": "ABRicate",
        "citation": "Seemann T. ABRicate: mass screening of contigs for antimicrobial resistance or virulence genes. GitHub repository.",
        "url": "https://github.com/tseemann/abricate",
    },
    "ncbi": {
        "name": "NCBI AMR resources",
        "citation": "NCBI antimicrobial resistance reference resources used through ABRicate-compatible databases.",
        "url": "https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/",
    },
    "vfdb": {
        "name": "VFDB",
        "citation": "Chen L, Yang J, Yu J, et al. VFDB: a reference database for bacterial virulence factors. Nucleic Acids Research. 2005;33(Database issue):D325-D328.",
        "url": "http://www.mgc.ac.cn/VFs/",
    },
    "plasmidfinder": {
        "name": "PlasmidFinder",
        "citation": "Carattoli A, Zankari E, Garcia-Fernandez A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrobial Agents and Chemotherapy. 2014;58(7):3895-3903.",
        "url": "https://cge.food.dtu.dk/services/PlasmidFinder/",
    },
    "mobileelementfinder": {
        "name": "MobileElementFinder",
        "citation": "Johansson MHK, Bortolaia V, Tansirichaiya S, et al. Detection of mobile genetic elements associated with antibiotic resistance in Salmonella enterica using a newly developed web tool: MobileElementFinder. Journal of Antimicrobial Chemotherapy. 2021;76(1):101-109.",
        "url": "https://bitbucket.org/genomicepidemiology/mobileelementfinder/",
    },
    "isfinder": {
        "name": "ISfinder",
        "citation": "ISfinder database for insertion sequence annotation.",
        "url": "https://isfinder.biotoul.fr/",
    },
    "integronfinder": {
        "name": "IntegronFinder",
        "citation": "Cury J, Jove T, Touchon M, Neron B, Rocha EPC. Identification and analysis of integrons and cassette arrays in bacterial genomes. Nucleic Acids Research. 2016;44(10):4539-4550.",
        "url": "https://github.com/gem-pasteur/Integron_Finder",
    },
    "iceberg": {
        "name": "ICEberg",
        "citation": "ICEberg database or ICE/IME/CIME annotation source used for ICEberg-style inputs. PanR2 does not run ICEberg directly.",
        "url": "https://db-mml.sjtu.edu.cn/ICEberg/",
    },
    "mlst": {
        "name": "MLST",
        "citation": "Seemann T. mlst: scan contig files against PubMLST typing schemes. GitHub repository.",
        "url": "https://github.com/tseemann/mlst",
    },
    "pubmlst": {
        "name": "PubMLST",
        "citation": "Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Research. 2018;3:124.",
        "url": "https://pubmlst.org/",
    },
    "defensefinder": {
        "name": "DefenseFinder",
        "citation": "Tesson F, Herve A, Mordret E, et al. Systematic and quantitative view of the antiviral arsenal of prokaryotes. Nature Communications. 2022;13:2561.",
        "url": "https://github.com/mdmparis/defense-finder",
    },
    "prophage": {
        "name": "Prophage or viral-region annotation source",
        "citation": "User-provided prophage or viral-region annotation tables parsed by PanR2. Cite the upstream prophage detection tool or database used to produce these inputs.",
        "url": "https://github.com/Tasnimul-Arabi-Anik/PanR2",
    },
    "mobsuite": {
        "name": "MOB-suite",
        "citation": "Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. Microbial Genomics. 2018;4(8).",
        "url": "https://github.com/phac-nml/mob-suite",
    },
    "kleborate": {
        "name": "Kleborate",
        "citation": "Kleborate was used or parsed for Klebsiella species-complex typing, virulence, resistance, and locus summaries. Cite the current Kleborate publication/documentation used for the run.",
        "url": "https://github.com/klebgenomics/Kleborate",
    },
    "kaptive": {
        "name": "Kaptive",
        "citation": "Kaptive was used or parsed for bacterial capsule/O-antigen locus typing. Cite the current Kaptive publication/documentation and database used for the run.",
        "url": "https://github.com/klebgenomics/Kaptive",
    },
    "ectyper": {
        "name": "ECTyper",
        "citation": "ECTyper was used or parsed for Escherichia coli serotype prediction. Cite the current ECTyper publication/documentation used for the run.",
        "url": "https://github.com/phac-nml/ecoli_serotyping",
    },
    "serotypefinder": {
        "name": "SerotypeFinder",
        "citation": "SerotypeFinder was used or parsed for serotype prediction. Cite the current CGE SerotypeFinder publication/documentation and database used for the run.",
        "url": "https://bitbucket.org/genomicepidemiology/serotypefinder/",
    },
    "sccmecfinder": {
        "name": "SCCmecFinder",
        "citation": "SCCmecFinder was used or parsed for Staphylococcus SCCmec typing. Cite the current CGE SCCmecFinder publication/documentation and database used for the run.",
        "url": "https://bitbucket.org/genomicepidemiology/sccmecfinder/",
    },
    "python": {
        "name": "Python scientific stack",
        "citation": "pandas, NumPy, SciPy, matplotlib, seaborn, and Plotly were used for analysis and visualization.",
        "url": "https://www.python.org/",
    },
}


def _package_version(package):
    try:
        return importlib_metadata.version(package)
    except Exception:
        return "not installed"


def _read_tool_manifest(output_dir):
    path = os.path.join(output_dir, "qc", "panr2_tool_manifest.csv")
    if os.path.exists(path):
        try:
            return pd.read_csv(path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    return pd.DataFrame()


def _metadata_engine_from_input(input_files):
    ncbi_clean = input_files.get("ncbi_clean")
    if not ncbi_clean:
        return ""
    metadata_dir = os.path.dirname(ncbi_clean)
    engine_path = os.path.join(metadata_dir, "metadata_engine.txt")
    if os.path.exists(engine_path):
        try:
            with open(engine_path, encoding="utf-8") as handle:
                return handle.read().strip().lower()
        except OSError:
            return ""
    if os.path.exists(os.path.join(metadata_dir, "fetchm2_clean.csv")):
        return "fetchm2"
    try:
        cols = pd.read_csv(ncbi_clean, nrows=1).columns
        if any(col in cols for col in ["Host_SD", "Collection_Year", "Metadata Fetch Status", "Host_Review_Status"]):
            return "fetchm2"
    except Exception:
        return ""
    return ""


def _selected_citation_keys(options=None, feature_outputs=None, input_files=None, output_dir=None):
    options = options or {}
    feature_outputs = feature_outputs or {}
    input_files = input_files or {}
    keys = {"panr2", "python"}
    if input_files.get("ncbi_clean") or options.get("ncbi_dir"):
        engine = _metadata_engine_from_input(input_files)
        keys.add("fetchm2" if engine == "fetchm2" else "fetchm")
    if options.get("abricate_dir") or options.get("run_abricate") or input_files.get("abricate_summary"):
        keys.update({"abricate", "ncbi"})
    for feature_type in feature_outputs:
        if feature_type in CITATION_REGISTRY:
            keys.add(feature_type)
        if feature_type == "mlst":
            keys.add("pubmlst")
        if feature_type in {"vfdb", "plasmidfinder", "isfinder", "iceberg"}:
            keys.add("abricate")
    if output_dir:
        manifest = _read_tool_manifest(output_dir)
        if not manifest.empty and "tool" in manifest.columns:
            for tool in manifest["tool"].dropna().astype(str).str.lower().unique():
                if tool in CITATION_REGISTRY:
                    keys.add(tool)
                if tool == "iceberg_table_converter":
                    keys.add("iceberg")
        if not manifest.empty and "database" in manifest.columns:
            for database in manifest["database"].dropna().astype(str).str.lower().unique():
                if database in CITATION_REGISTRY:
                    keys.add(database)
    return [key for key in CITATION_REGISTRY if key in keys]


def _bibtex_key(key):
    return "".join(ch for ch in key.title() if ch.isalnum())


def write_citation_outputs(output_dir, options=None, feature_outputs=None, input_files=None, panr2_version="unknown"):
    """Write citation and software-version files for a PanR2 run."""
    report_dir = os.path.join(output_dir, "report")
    os.makedirs(report_dir, exist_ok=True)
    keys = _selected_citation_keys(options, feature_outputs, input_files, output_dir)

    citations_md = os.path.join(report_dir, "citations.md")
    citations_bib = os.path.join(report_dir, "citations.bib")
    software_versions = os.path.join(report_dir, "software_versions.csv")

    lines = ["# PanR2 Citation Report", ""]
    lines.append("Cite PanR2 and the tools/databases actually used in this run. Verify database-specific citation wording against the current upstream documentation when preparing a manuscript.")
    lines.append("")
    for key in keys:
        entry = CITATION_REGISTRY[key]
        lines.append(f"## {entry['name']}")
        lines.append("")
        lines.append(entry["citation"])
        lines.append("")
        lines.append(f"URL: {entry['url']}")
        lines.append("")
    with open(citations_md, "w") as handle:
        handle.write("\n".join(lines).rstrip() + "\n")

    bib_lines = []
    for key in keys:
        entry = CITATION_REGISTRY[key]
        bib_lines.extend([
            f"@misc{{{_bibtex_key(key)},",
            f"  title = {{{entry['name']}}},",
            f"  howpublished = {{{entry['url']}}},",
            f"  note = {{{entry['citation']}}}",
            "}",
            "",
        ])
    with open(citations_bib, "w") as handle:
        handle.write("\n".join(bib_lines).rstrip() + "\n")

    rows = [
        {"component": "PanR2", "version": panr2_version, "source": "package"},
        {"component": "Python", "version": sys.version.split()[0], "source": platform.platform()},
    ]
    for package in ["pandas", "numpy", "scipy", "matplotlib", "seaborn", "plotly"]:
        rows.append({"component": package, "version": _package_version(package), "source": "python_package"})
    manifest = _read_tool_manifest(output_dir)
    if not manifest.empty:
        for _, row in manifest.iterrows():
            rows.append({
                "component": row.get("tool", ""),
                "version": row.get("version", ""),
                "source": row.get("database", ""),
                "database_sequences": row.get("database_sequences", ""),
                "database_date": row.get("database_date", ""),
            })
    with open(software_versions, "w", newline="") as handle:
        fieldnames = ["component", "version", "source", "database_sequences", "database_date"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    return {"citations_md": citations_md, "citations_bib": citations_bib, "software_versions": software_versions}
