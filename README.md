# PanR2: Panresistome Analysis Tool

## Overview
PanR2 is a comprehensive Python-based tool for analyzing panresistome data (Global Resistance Pattern). It processes NCBI and Abricate summary files, merges the data, and generates a wide range of visualizations including heatmaps, bar plots, boxplots, and interactive HTML plots. The tool is designed to help researchers analyze and visualize antibiotic resistance gene presence, prevalence, and distribution patterns across different geographic locations and temporal scales.

It enables robust statistical analysis, including group-wise comparisons, summary statistics, and correlation assessments, to support meaningful interpretation of resistome data. Given sufficient sequencing data, PanR2 can help identify local resistance spread, compare resistance patterns between regions or over time, and detect signals of emerging resistance epidemics. This makes it a valuable tool for surveillance, epidemiology, and public health research.

**Prerequisites:**
- `ncbi_clean.csv` from [FetchM2](https://github.com/Tasnimul-Arabi-Anik/FetchM2), PanResistome, or a FetchM-compatible metadata workflow
- Summary files in `.tab` (preferred) or `.csv` format from [Abricate](https://github.com/tseemann/abricate)

### Key Features:
- Merges and analyzes standardized metadata and ABRicate outputs
- Calculates gene presence/absence across multiple categories, including geography, collection year/date, host, sample type, isolation source, environment, and organism metadata when present
- Performs comprehensive statistical tests and correlation analyses
- Generates multiple visualization types: heatmaps, bar plots, boxplots, lollipop plots, and correlation plots
- Creates interactive HTML visualizations for enhanced data exploration
- Generates an interactive HTML index for easy navigation of all results
- Provides detailed statistical analysis outputs

---

## Installation

PanR2 has three practical install modes:

- Analysis-only install: use PanR2 with existing ABRicate-style result folders. This does not install external annotation tools or databases.
- Source checkout with `environment.yml`: installs PanR2 Python dependencies plus ABRicate, IntegronFinder, MLST, DefenseFinder, BLAST/KMA support, and MobileElementFinder. ABRicate databases still need to be initialized with `abricate --setupdb`.
- Container install: builds a reproducible command-line image for users who want fewer local dependency conflicts.

### Method 1: Analysis-only `pip` Install
```bash
conda create -n panr_env python=3.9
conda activate panr_env
pip install panR2
```

### Method 2: Source Checkout with Integrated Tool Dependencies
```bash
git clone https://github.com/Tasnimul-Arabi-Anik/PanR2.git
cd PanR2
conda env create -f environment.yml
conda activate panr2
pip install -e .
abricate --setupdb
```

### Method 3: Direct Installation from GitHub
```bash
conda create -n panr_env python=3.9
conda activate panr_env
pip install git+https://github.com/Tasnimul-Arabi-Anik/PanR2.git
```

### Method 4: Docker/Container Install
Build the container from a source checkout:

```bash
docker build -t panr2 .
```

To also initialize ABRicate databases during the image build:

```bash
docker build --build-arg SETUP_ABRICATE_DB=true -t panr2 .
```

Run PanR2 from the container by mounting your data directory:

```bash
docker run --rm -v "$PWD":/work -w /work panr2 --help
```

### Confirm Installation
```bash
panr --help
panr --doctor
panr doctor --json
panr install-info
```

`panr --doctor` reports whether the Python dependencies, optional external annotation tools, and ABRicate databases are visible in the current environment. Missing ABRicate, MobileElementFinder, IntegronFinder, MLST, or DefenseFinder is not fatal for analysis-only mode, but those tools are required when using the corresponding integrated runner flags.

For integrated ABRicate runs, initialize or check ABRicate databases with:

```bash
panr setup-db --dbs ncbi,vfdb,plasmidfinder --check-only
panr setup-db --dbs ncbi,vfdb,plasmidfinder
```

`panr doctor --fix` can run safe setup fixes when possible. At present, that means running ABRicate database setup when ABRicate is installed but no databases are visible.

### Run the Included Smoke Test
From a source checkout, run the bundled small fixture to confirm the CLI, merge step, tidy output, and plotting path work locally:

```bash
panr validate-demo --output-dir test_output --format png
```

The validation command writes a complete small multi-database run, including MLST, DefenseFinder, prophage table inputs, sample-map matching, temporal trends, and the main dashboard at `test_output/report/index.html`.

You can also run the equivalent explicit command:

```bash
python bin/panr \
  --ncbi-dir tests/fixtures/ncbi \
  --abricate-dir tests/fixtures/abricate \
  --output-dir test_output \
  --format png \
  --genep 0 \
  --nseq 1 \
  --min-identity 90 \
  --min-samples-per-group 2 \
  --core-threshold 75 \
  --rare-threshold 25 \
  --top-n 10 \
  --cooccurrence-min-prevalence 0 \
  --cooccurrence-top-n 10 \
  --vfdb-dir tests/fixtures/vfdb \
  --plasmidfinder-dir tests/fixtures/plasmidfinder \
  --mobileelementfinder-dir tests/fixtures/mobileelementfinder \
  --isfinder-dir tests/fixtures/isfinder \
  --integronfinder-dir tests/fixtures/integronfinder \
  --iceberg-dir tests/fixtures/iceberg
```

The fixture is intentionally small. Correlation plots that require at least five samples per group may be skipped, while the rest of the outputs should be generated.

---

## Source Layout

PanR2 keeps the installed command as `panr`, while the implementation lives in the `panr2` Python package. The `bin/panr` file is a compatibility wrapper around `panr2.cli`. Core code is split into focused modules for I/O, QC, filtering, AMR analysis, optional feature analysis, integrated runners, MLST, temporal trends, cross-database associations, reports, dashboards, plots, and statistics.

---

## Usage

### Command-Line Interface
```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

Utility subcommands are also available:

```bash
panr doctor [--json] [--fix]
panr setup-db [--dbs ncbi,vfdb,plasmidfinder] [--check-only] [--json]
panr install-info
panr citations --output-dir <OUTPUT_DIRECTORY>
panr validate-demo --output-dir <OUTPUT_DIRECTORY>
panr run-all --ncbi-dir <NCBI_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY>
```

PanR2 can also run ABRicate internally from assembly FASTA files:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-abricate --abricate-dbs ncbi,vfdb,plasmidfinder --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

For a broad integrated run, use `panr run-all` or `--run-all`. This enables ABRicate, MobileElementFinder, IntegronFinder, MLST, and DefenseFinder runners, then performs AMR, optional database, temporal trend, and cross-database comparative analyses:

```bash
panr run-all --ncbi-dir <NCBI_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> --min-identity 90
```

MobileElementFinder can also be run internally when installed:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-mobileelementfinder --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

IntegronFinder can also be run internally when installed:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-integronfinder --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

MLST can be run internally when installed, or existing `mlst` TSV/CSV output can be supplied:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-mlst --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --mlst-dir <MLST_OUTPUT_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

DefenseFinder can be run internally when installed, or existing DefenseFinder TSV/CSV outputs can be converted into PanR2-compatible feature inputs:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-defensefinder --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --defensefinder-dir <DEFENSEFINDER_TABLE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

User-provided prophage or viral-region tables can be converted into PanR2 feature inputs. PanR2 does not currently run a prophage caller directly; cite the upstream tool used to create the supplied tables:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --prophage-dir <PROPHAGE_TABLE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

MOB-suite plasmid reconstruction tables and organism-specific typing outputs can also be supplied as table inputs. PanR2 converts them into the same internal feature contract used for AMR, VFDB, plasmid, MGE, and integron analyses:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> \
  --mobsuite-dir <MOBSUITE_TABLE_DIRECTORY> \
  --kleborate-dir <KLEBORATE_TABLE_DIRECTORY> \
  --kaptive-dir <KAPTIVE_TABLE_DIRECTORY> \
  --ectyper-dir <ECTYPER_TABLE_DIRECTORY> \
  --serotypefinder-dir <SEROTYPEFINDER_TABLE_DIRECTORY> \
  --sccmecfinder-dir <SCCMECFINDER_TABLE_DIRECTORY> \
  --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

### Recommended Stable Workflow

For publication-oriented analyses, the most reliable workflow is:

1. Use FetchM2/PanResistome or another external workflow to download assemblies, run QC, run annotation tools, and capture tool/database versions.
2. Keep heavy tools and databases outside PanR2 where possible.
3. Convert each tool output into PanR2-compatible database folders or table inputs.
4. Provide `--sample-map` whenever tool outputs use local sample names instead of GCF/GCA assembly accessions.
5. Run PanR2 for metadata-linked summaries, cross-database associations, temporal trends, figures, dashboard output, citations, and journal-style reporting.

PanR2 and PanResistome share a formal feature-table contract documented in [`docs/panr2_input_contract.md`](docs/panr2_input_contract.md). In short, PanResistome should run heavyweight tools such as CheckM2, GTDB-Tk, QUAST, FastANI/skani, Mash, ABRicate, MLST, DefenseFinder, MobileElementFinder, and IntegronFinder, then export standardized tables. PanR2 should remain the lightweight statistical/reporting layer that reads those tables.

### Sample Naming And Sample Map

PanR2 preserves assembly accession versions such as `GCF_000123456.1`. If all inputs contain matching GCF/GCA accessions, no sample map is needed. If files or tool outputs use local names such as `sample_001.fna`, provide:

```bash
panr --sample-map sample_map.csv ...
```

The sample map may be CSV or TSV and must contain:

```csv
sample_id,Assembly Accession
sample_001,GCF_000123456.1
sample_002,GCA_000987654.1
```

PanR2 uses this map for ABRicate summaries, MLST outputs, DefenseFinder/prophage tables, ICEberg-style tables, and other PanR2-compatible feature inputs. Mapping diagnostics are written to `qc/sample_map_qc_<source>.csv` and `qc/sample_map_qc_<source>_summary.txt`.

Example templates are available under `examples/input_templates/`.

### Biological Interpretation Limits

Cross-database association outputs are sample/genome-level screening analyses. They do not prove physical linkage, plasmid localization, horizontal transfer, transfer direction, phenotype, shared regulation, or causality. Temporal trend outputs are also screening-level summaries and depend on metadata completeness, collection-year balance, sampling intensity, and repeated tied prevalence values.

User-provided ICE/IME/CIME tables can be converted into ICEberg-style PanR2 analysis inputs. PanR2 does not run an ICEberg annotation program directly; it converts existing ICE/IME/CIME annotation tables into PanR2-compatible feature tables:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --iceberg-table-dir <ICE_TABLE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

### Required Arguments
| Argument | Description |
|----------|-------------|
| `--ncbi-dir` | Directory containing standardized `ncbi_clean.csv` from FetchM2, PanResistome, or a FetchM-compatible workflow |
| `--output-dir` | Directory to store merged results and visualizations |

`--abricate-dir` is required for analysis-only mode. It is not required when `--run-abricate` is used with the `ncbi` database.

### Optional Arguments
| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--doctor` | flag | off | Report installed Python dependencies, optional external tools, and ABRicate database visibility, then exit |
| `--json` | flag | off | With `--doctor` or setup subcommands, write machine-readable JSON output |
| `--fix` | flag | off | With `--doctor`, run safe fixes when possible, currently ABRicate database setup |
| `--install-info` | flag | off | Print PanR2/Python/system/tool/database readiness information, then exit |
| `--setup-db` | flag | off | Run ABRicate database setup, then exit |
| `--check-only` | flag | off | With `--setup-db`, only report ABRicate database visibility |
| `--citations` | flag | off | Write citation/software-version files for `--output-dir`, then exit |
| `--genep` | float | `10.0` | Minimum % gene presence to include in heatmap |
| `--nseq` | int | `1` | Minimum number of sequences required per group in heatmaps |
| `--format` | str | `tiff` | Output format for figures (`tiff`, `svg`, `png`, `pdf`) |
| `--sequence-dir` | path | optional | Directory containing assembly FASTA files for integrated tool runners |
| `--run-all` | flag | off | Run all currently integrated annotation runners from `--sequence-dir`, then perform all PanR2 analyses |
| `--run-abricate` | flag | off | Run ABRicate internally before PanR2 analysis |
| `--abricate-dbs` | str | `ncbi` | Comma-separated ABRicate databases to run, for example `ncbi,vfdb,plasmidfinder` |
| `--abricate-bin` | path | `abricate` | ABRicate executable name or path |
| `--abricate-summary-metric` | str | `identity` | Metric used in generated ABRicate summary matrices (`identity` or `coverage`) |
| `--run-mobileelementfinder` | flag | off | Run MobileElementFinder internally before PanR2 feature analysis |
| `--mobileelementfinder-bin` | path | `mefinder` | MobileElementFinder executable name or path |
| `--mobileelementfinder-threads` | int | `1` | Threads passed to MobileElementFinder |
| `--run-integronfinder` | flag | off | Run IntegronFinder internally before PanR2 feature analysis |
| `--integronfinder-bin` | path | `integron_finder` | IntegronFinder executable name or path |
| `--integronfinder-threads` | int | `1` | CPU threads passed to IntegronFinder |
| `--run-mlst` | flag | off | Run `mlst` internally before PanR2 typing analysis |
| `--mlst-bin` | path | `mlst` | MLST executable name or path |
| `--run-defensefinder` | flag | off | Run DefenseFinder internally before PanR2 feature analysis |
| `--defensefinder-bin` | path | `defense-finder` | DefenseFinder executable name or path |
| `--iceberg-table-dir` | path | optional | Directory containing ICE/IME/CIME CSV/TSV/TAB tables to convert into PanR2 ICEberg analysis inputs |
| `--force-tool-run` | flag | off | Re-run integrated tools even when raw result files already exist |
| `--min-identity` | float | `0.0` | Minimum ABRicate identity percentage to treat a gene call as present |
| `--drop-unmatched-accessions` | flag | off | Drop NCBI rows with no matching ABRicate summary row |
| `--min-samples-per-group` | int | `5` | Minimum samples per group required for correlation analyses |
| `--core-threshold` | float | `95.0` | Prevalence percentage used to classify core ARGs |
| `--rare-threshold` | float | `5.0` | Prevalence percentage used to classify rare ARGs |
| `--top-n` | int | `25` | Number of top genes/classes to include in compact summary plots |
| `--cooccurrence-min-prevalence` | float | `0.0` | Minimum prevalence percentage for genes/classes included in co-occurrence matrices |
| `--cooccurrence-top-n` | int | `25` | Number of top genes/classes or pairs to include in co-occurrence plots and pair tables |
| `--plot-style` | str | `publication` | Integrated plot readability preset: `publication`, `dashboard`, or `compact` |
| `--label-max-length` | int | style default | Maximum displayed feature-label length in crowded integrated figures |
| `--no-cross-database` | flag | off | Disable integrated cross-database association outputs |
| `--cross-database-max-features` | int | `300` | Maximum most-prevalent features used for pairwise cross-database statistics; use `0` for no limit |
| `--no-temporal-trends` | flag | off | Disable advanced temporal trend outputs |
| `--vfdb-dir` | path | optional | Directory containing ABRicate VFDB summary/results files |
| `--plasmidfinder-dir` | path | optional | Directory containing ABRicate PlasmidFinder summary/results files |
| `--mobileelementfinder-dir` | path | optional | Directory containing ABRicate MobileElementFinder summary/results files |
| `--isfinder-dir` | path | optional | Directory containing ABRicate ISfinder summary/results files |
| `--integronfinder-dir` | path | optional | Directory containing IntegronFinder or ABRicate-style integron summary/results files |
| `--iceberg-dir` | path | optional | Directory containing ABRicate ICEberg summary/results files |
| `--mlst-dir` | path | optional | Directory containing MLST TSV/CSV output |
| `--defensefinder-dir` | path | optional | Directory containing DefenseFinder tables or PanR2-compatible summary/results files |
| `--prophage-dir` | path | optional | Directory containing prophage/viral-region tables or PanR2-compatible summary/results files |
| `--mobsuite-dir` | path | optional | Directory containing MOB-suite tables or PanR2-compatible summary/results files |
| `--kleborate-dir` | path | optional | Directory containing Kleborate tables or PanR2-compatible summary/results files |
| `--kaptive-dir` | path | optional | Directory containing Kaptive tables or PanR2-compatible summary/results files |
| `--ectyper-dir` | path | optional | Directory containing ECTyper tables or PanR2-compatible summary/results files |
| `--serotypefinder-dir` | path | optional | Directory containing SerotypeFinder tables or PanR2-compatible summary/results files |
| `--sccmecfinder-dir` | path | optional | Directory containing SCCmecFinder tables or PanR2-compatible summary/results files |
| `--sample-map` | path | optional | CSV/TSV mapping local sample IDs or filenames to `Assembly Accession` |
| `--version` | - | - | Show program's version number and exit |
| `-h, --help` | - | - | Show help message and exit |

### Example Usage
```bash
# Basic usage
panr --ncbi-dir ./data/ncbi_clean.csv --abricate-dir ./data/abricate --output-dir ./output

# With optional parameters
panr --ncbi-dir ./data/ncbi_clean.csv --abricate-dir ./data/abricate --output-dir ./output --format pdf --genep 10 --nseq 5

# Apply analysis filtering before plotting
panr --ncbi-dir ./data/ncbi --abricate-dir ./data/abricate --output-dir ./output_filtered --min-identity 90 --drop-unmatched-accessions

# Run ABRicate inside PanR2, then analyze AMR, VFDB, and PlasmidFinder outputs
panr --ncbi-dir ./data/metadata_output --sequence-dir ./data/sequence --run-abricate --abricate-dbs ncbi,vfdb,plasmidfinder --output-dir ./output_integrated --min-identity 90

# Run MobileElementFinder inside PanR2 and include MGE feature analysis
panr --ncbi-dir ./data/metadata_output --abricate-dir ./data/abricate --sequence-dir ./data/sequence --run-mobileelementfinder --output-dir ./output_mge --min-identity 90

# Run IntegronFinder inside PanR2 and include integron feature analysis
panr --ncbi-dir ./data/metadata_output --abricate-dir ./data/abricate --sequence-dir ./data/sequence --run-integronfinder --output-dir ./output_integrons --min-identity 90

# Add MLST typing and DefenseFinder/prophage table inputs to the same comparative workflow
panr --ncbi-dir ./data/metadata_output --abricate-dir ./data/abricate --mlst-dir ./data/mlst --defensefinder-dir ./data/defensefinder --prophage-dir ./data/prophage --output-dir ./output_comparative --min-identity 90

# Convert ICE/IME/CIME annotation tables and include ICEberg-style feature analysis
panr --ncbi-dir ./data/metadata_output --abricate-dir ./data/abricate --iceberg-table-dir ./data/iceberg_tables --output-dir ./output_iceberg --min-identity 90
```

---

## Output Structure

PanR2 first writes an input QC report under `qc/`. Review this report before interpreting downstream figures, especially for unmatched accessions, zero-hit samples, or missing resistance annotations.

PanR2 generates a comprehensive set of outputs organized in the following directory structure:

```
output/
├── tool_results/                      # Raw outputs from integrated upstream tools
│   └── abricate/
│       ├── ncbi/
│       ├── vfdb/
│       └── plasmidfinder/
│   └── mobileelementfinder/
│       ├── raw/
│       └── panr2_inputs/
│   └── integronfinder/
│       ├── raw/
│       └── panr2_inputs/
│   └── mlst/
│       └── raw/
│   └── defensefinder/
│       ├── raw/
│       └── panr2_inputs/
│   └── iceberg/
│       └── panr2_inputs/
├── ncbi/                              # AMR/resistome tables, merged data, and figures
│   ├── analysis/                      # Panresistome summary tables and compact plots
│   ├── figures/                       # AMR static and interactive figures
│   └── merged_output/                 # Merged and tidy AMR tables
├── vfdb/                              # Optional VFDB feature tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── plasmidfinder/                     # Optional PlasmidFinder feature tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── mobileelementfinder/               # Optional mobile genetic element tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── isfinder/                          # Optional insertion sequence tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── integronfinder/                    # Optional integron feature tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── iceberg/                           # Optional ICE/IME feature tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── mlst/                              # Optional sequence typing summaries and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── defensefinder/                     # Optional defense-system feature tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── prophage/                          # Optional prophage/viral-region feature tables and figures
│   ├── analysis/
│   ├── figures/
│   └── merged_output/
├── mobsuite/                          # Optional MOB-suite plasmid reconstruction/typing summaries
├── kleborate/                         # Optional Klebsiella typing/virulence/resistance summaries
├── kaptive/                           # Optional capsule/O-locus typing summaries
├── ectyper/                           # Optional E. coli serotype summaries
├── serotypefinder/                    # Optional serotype summaries
├── sccmecfinder/                      # Optional SCCmec typing summaries
├── cross_database/                    # Integrated AMR/VFDB/plasmid/MGE comparative genomics
│   ├── analysis/
│   └── figures/
├── temporal/                          # Advanced temporal trend tables and HTML plots
│   ├── analysis/
│   └── figures/
├── panresistome_context/              # Optional PanResistome QC/ANI/QUAST summaries if upstream files are present
│   └── analysis/
├── qc/                                # Shared input validation and filter reports
└── report/                            # Dashboard, journal-style report, citations, and methods text
```

### Output Files Description

#### 1. Input QC (`qc/` directory)
- **`panr2_input_qc.csv`** - Machine-readable input checks with `PASS`, `WARN`, `FAIL`, and `INFO` statuses
- **`panr2_input_qc_summary.txt`** - Human-readable summary of input checks
- **`panr2_tool_manifest.csv`** and **`panr2_tool_manifest.json`** - Tool versions, selected databases, database dates/counts where available, and raw output paths for integrated runs
- **`panr2_unmatched_accessions.csv`** - Accessions present in only one of the NCBI or ABRicate inputs
- **`*_filter_report.csv`** - Row and ARG-call counts before and after optional analysis filters
- **`metadata_completeness_report.csv`** - Completeness, missingness, and status for geography, host/source, environment, disease/health, organism, and FetchM2-style standardized metadata fields
- **`metadata_group_sample_sizes.csv`** - Per-group sample counts and underpowered-group flags for metadata-driven analyses
- **`metadata_bias_warning.txt`** - Human-readable warnings for incomplete or biased metadata fields
- **`sample_map_qc_<source>.csv`** and **`sample_map_qc_<source>_summary.txt`** - Sample-map matching diagnostics for each mapped input source

#### 2. NCBI/AMR Panresistome Analysis (`ncbi/` directory)
- **`*_sample_resistome_burden.csv`** - Per-sample ARG burden, resistance class count, and identity summary
- **`*_gene_prevalence_summary.csv`** - Per-gene prevalence, identity range, resistance class, and geographic spread
- **`*_resistome_category_summary.csv`** - Core/accessory/rare ARG categories using configurable thresholds
- **`*_resistance_class_summary.csv`** - Resistance class prevalence, gene diversity, and top genes
- **`*_gene_cooccurrence_matrix.csv`** - Gene-by-gene co-occurrence counts across samples
- **`*_class_cooccurrence_matrix.csv`** - Resistance-class co-occurrence counts across samples
- **`*_top_gene_pairs.csv`** - Ranked ARG pairs with support, prevalence, and Jaccard index
- **`*_top_class_pairs.csv`** - Ranked resistance-class pairs with support, prevalence, and Jaccard index
- **`ncbi/analysis/plots/`** - Compact prevalence and co-occurrence plots intended for quick review, not exhaustive visualization

#### 3. Optional Database-Named Feature Analysis
- **`vfdb/analysis/vfdb_feature_summary.csv`** - VFDB feature prevalence and identity summaries
- **`vfdb/analysis/vfdb_qc_summary.csv`** - VFDB sample/feature QC metrics, zero-feature sample counts, unmatched sample counts, and top feature list
- **`vfdb/analysis/vfdb_category_summary.csv`** - Virulence product/category summaries
- **`vfdb/analysis/vfdb_geographic_summary.csv`** - Descriptive VFDB feature burden by metadata group
- **`vfdb/analysis/vfdb_temporal_summary.csv`** - Descriptive VFDB feature burden by collection year
- **`vfdb/analysis/vfdb_feature_cooccurrence_matrix.csv`** - VFDB feature co-occurrence counts across samples
- **`vfdb/analysis/vfdb_top_feature_pairs.csv`** - Ranked VFDB feature pairs with support and Jaccard index
- **`vfdb/analysis/vfdb_group_burden_summary.csv`** - VFDB feature burden summaries by location, continent, subcontinent, and year
- **`vfdb/analysis/vfdb_group_overall_tests.csv`** - VFDB nonparametric group-level tests where sample sizes permit
- **`vfdb/analysis/vfdb_group_pairwise_tests.csv`** - VFDB pairwise group comparisons where sample sizes permit
- **`vfdb/figures/`** - VFDB static figures, burden-by-group plots, interactive HTML files under `html_files/`, and `index.html` navigation
- **`plasmidfinder/analysis/plasmidfinder_feature_summary.csv`** - PlasmidFinder replicon prevalence and identity summaries
- **`plasmidfinder/analysis/plasmidfinder_qc_summary.csv`** - PlasmidFinder sample/feature QC metrics, zero-feature sample counts, unmatched sample counts, and top feature list
- **`plasmidfinder/analysis/plasmidfinder_category_summary.csv`** - Plasmid replicon/category summaries
- **`plasmidfinder/analysis/plasmidfinder_geographic_summary.csv`** - Descriptive plasmid replicon burden by metadata group
- **`plasmidfinder/analysis/plasmidfinder_temporal_summary.csv`** - Descriptive plasmid replicon burden by collection year
- **`plasmidfinder/analysis/plasmidfinder_feature_cooccurrence_matrix.csv`** - Plasmid replicon co-occurrence counts across samples
- **`plasmidfinder/analysis/plasmidfinder_top_feature_pairs.csv`** - Ranked plasmid replicon pairs with support and Jaccard index
- **`plasmidfinder/analysis/plasmidfinder_group_burden_summary.csv`** - PlasmidFinder feature burden summaries by location, continent, subcontinent, and year
- **`plasmidfinder/analysis/plasmidfinder_group_overall_tests.csv`** - PlasmidFinder nonparametric group-level tests where sample sizes permit
- **`plasmidfinder/analysis/plasmidfinder_group_pairwise_tests.csv`** - PlasmidFinder pairwise group comparisons where sample sizes permit
- **`plasmidfinder/figures/`** - PlasmidFinder static figures, burden-by-group plots, interactive HTML files under `html_files/`, and `index.html` navigation

The same database-named output pattern is used for `mobileelementfinder/`, `isfinder/`, `integronfinder/`, `iceberg/`, `defensefinder/`, `prophage/`, `mobsuite/`, `kleborate/`, `kaptive/`, `ectyper/`, `serotypefinder/`, and `sccmecfinder/`. Each optional feature database writes `analysis/`, `figures/`, and `merged_output/` folders with feature summaries, product/category summaries, sample burden, per-database QC summaries, unmatched-sample reports, geographic and temporal summaries, co-occurrence tables, group-burden comparisons, static plots, interactive HTML plots, and an `index.html` figure browser when applicable.

VFDB, PlasmidFinder, MobileElementFinder, ISfinder, IntegronFinder, ICEberg, DefenseFinder, prophage, MOB-suite, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder inputs are handled as separate feature families, not as antibiotic resistance classes. PanR2 therefore uses feature prevalence, product/replicon/typing categories, sample burden, geography and temporal summaries, identity distributions where applicable, feature co-occurrence, feature presence heatmaps, and interactive HTML figures instead of resistance-class composition plots. MLST is treated as typing metadata/features, not as AMR, virulence, plasmid, or mobile-element annotation.

MLST-specific outputs include `mlst/analysis/sample_mlst_summary.csv`, `mlst/analysis/mlst_by_metadata.csv`, and `mlst/analysis/st_feature_burden_summary.csv`.

#### 4. Cross-Database Comparative Genomics (`cross_database/` directory)
- **`cross_database_feature_matrix.csv`** - Unified sample-by-feature matrix with prefixed features such as `AMR:blaA`, `VFDB:fimH`, `PLASMID:IncFIB`, `MGE:IS26`, `INTEGRON:intI1`, `ICE:ICEKp1`, `MLST:scheme:ST11`, `DEFENSE:RM_type_I`, `PROPHAGE:pp1`, `MOB:rep_type:IncFIB`, `KLEBORATE:st:11`, and `ECTYPER:serotype:O157:H7`
- **`cross_database_top_associations.csv`** - Pairwise feature associations with co-occurrence count, Jaccard index, phi coefficient, Fisher exact-test odds ratio, p-value, and FDR-adjusted q-value
- **`cross_database_cooccurrence_matrix.csv`**, **`cross_database_jaccard_matrix.csv`**, and **`cross_database_phi_correlation_matrix.csv`** - Global association matrices
- **`amr_mge_associations.csv`**, **`amr_plasmid_associations.csv`**, **`amr_integron_associations.csv`**, **`amr_virulence_associations.csv`**, **`amr_defense_associations.csv`**, **`amr_prophage_associations.csv`**, **`plasmid_mge_associations.csv`**, **`defense_mge_associations.csv`**, and **`prophage_mge_associations.csv`** - Biologically focused cross-database association tables
- **`sample_integrated_feature_burden.csv`** - Per-sample AMR, virulence, plasmid, mobileome, total feature, and mobility-associated AMR burden metrics
- **`cross_database_feature_enrichment_by_metadata.csv`** - Feature enrichment by metadata group using Fisher exact tests and FDR correction
- **`global_feature_association_heatmap.*`**, **`integrated_feature_presence_heatmap.*`**, and **`cross_database_feature_network.html`** - Static and interactive integrated comparative figures
- **`figure_manifest.csv`** - Figure inventory with descriptions and recommended use
- **`plot_readability_warnings.csv`** - Warnings when labels or dense networks/heatmaps were shortened or limited for readability

Cross-database co-occurrence is sample/genome-level only. It does not prove physical linkage, plasmid localization, horizontal transfer, shared regulation, phenotype, or causality.

#### 5. Advanced Temporal Trends (`temporal/` directory)
- **`temporal_feature_trends.csv`** - Feature-level Mann-Kendall, Spearman, and logistic presence trend summaries by collection year, including FDR q-values
- **`temporal_burden_trends.csv`** - Linear and Mann-Kendall trend summaries for per-sample feature burdens, including FDR q-values
- **`temporal_top_feature_trends.html`** - Interactive yearly prevalence trends for selected features

#### 6. PanResistome QC/ANI/QUAST Context (`panresistome_context/` directory)
When PanResistome outputs are present in the PanR2 output directory, PanR2 writes lightweight comparative context summaries without requiring those heavy tools to be installed:

- **`qc_context_sample_burden.csv`** - PanResistome QC master metrics merged with PanR2 feature-burden values
- **`qc_master_status_summary.csv`** - PASS/WARN/FAIL counts from the upstream combined QC decision engine
- **`qc_feature_correlation_summary.csv`** - Correlations between QC metrics such as completeness, contamination, N50, contig count, ANI, and feature-burden metrics
- **`species_consistency_summary.csv`** - FastANI/skani species-consistency screening summary
- **`duplicate_cluster_summary.csv`** - Near-duplicate ANI cluster and representative-genome summary
- **`representative_samples.csv`** - Representative genomes selected from near-duplicate ANI clusters
- **`burden_by_ani_cluster.csv`** - AMR, virulence, plasmid, and mobileome burden by ANI cluster
- **`panresistome_context_manifest.csv`** - PanResistome QC/ANI/QUAST files detected by PanR2

These outputs support checks such as AMR burden versus completeness/contamination, plasmid or MGE burden by ANI cluster, species-mismatch review, duplicate-cluster review, and representative-genome selection. They are descriptive screening summaries and should be interpreted with the upstream PanResistome QC reports.

#### 7. Written Report (`report/` directory)
- **`index.html`** - Top-level dashboard linking QC, metadata completeness, AMR, optional database outputs, cross-database outputs, temporal trends, figures, citations, and software versions
- **`*_panr2_report.md`** - Comprehensive journal-style narrative report generated from output tables
- **`*_panr2_report.html`** - Simple HTML rendering of the Markdown report
- **`*_methods.txt`** - Reusable methods description for manuscript drafting
- **`citations.md`** and **`citations.bib`** - Run-specific citation files for PanR2 and detected tools/databases
- **`software_versions.csv`** - PanR2, Python package, and integrated-tool versions when available

#### 8. NCBI/AMR Static Visualizations
- **`Resistance_gene_presence.{format}`** - Bar plot showing gene presence across samples
- **`Resistance_gene_percentage.{format}`** - Lollipop plot showing gene percentage distribution
- **`Resistance_gene_identity_boxplot.{format}`** - Boxplot of resistance gene variation across sequences
- **`Resistance_percentage_by_Antibiotics.{format}`** - Bar plot of resistance by antibiotic classes

#### 9. Heatmaps (`ncbi/figures/heatmap/` directory)
- **`Resistance gene distribution by continent, geographic location, subcontinent, and year.`**

#### 10. Mean ARG Analysis (`ncbi/figures/mean_ARG/` directory)
- **`Average antibiotic resistance genes by continent, geographic location, subcontinent, and year.`** - 

#### 11. Interactive HTML Visualizations (`ncbi/figures/html_files/` directory)
- **`Resistance_gene_distribution_heatmap.html`** - Interactive heatmap of gene distribution
- **`Resistance_gene_geographic_distribution.html`** - Geographic distribution map
- **`Resistance_gene_frequency_boxplot.html`** - Interactive frequency analysis
- **`Resistance_gene_identity_boxplot.html`** - Interactive identity score analysis
- **`Resistance_gene_presence.html`** - Interactive presence/absence visualization
- **`Resistance_gene_percentage.html`** - Interactive percentage analysis
- **`Resistance_percentage_by_Antibiotics.html`** - Interactive antibiotic class analysis
- **`Mean_Frequency_Antibiotic_Resistance_genes.html`** - Mean frequency analysis
- **`Continent_correlation_plot.html`** - Continental correlation analysis
- **`Geographic_Location_correlation_plot.html`** - Location-based correlations
- **`Subcontinent_correlation_plot.html`** - Subcontinental correlation patterns

#### 12. Statistical Analysis (`ncbi/figures/Stat_analysis/` directory)
- **`combined_geographic_correlation_summary.csv`** - Geographic correlation statistics
- **`combined_overall_tests.csv`** - Overall statistical test results
- **`combined_pairwise_comparisons.csv`** - Pairwise comparison results
- **`combined_summary_statistics.csv`** - Comprehensive summary statistics
- **`ncbi_gene_presence_count_percentage.csv`** - Gene presence counts and percentages

#### 12. Navigation
- **It generates an interactive combined index.html file** 
- **[View the interactive HTML report](https://tasnimul-arabi-anik.github.io/PanR2/)** – Interactive HTML index page for easy navigation of all generated visualizations

## Example Visualizations

### Static Plots

**Mean ARG by Geographic Location**
![Mean Antibiotic Resistance Genes by Country](ncbi/figures/mean_ARG/Mean_ARG_by_Geographic%20Location.png)

**Resistance Gene Presence Analysis:**
![Resistance Gene Presence](ncbi/figures/Resistance_gene_presence.png)
*Bar plot showing the presence of resistance genes across samples*

**Gene Percentage Distribution:**
![Resistance Gene Percentage](ncbi/figures/Resistance_gene_percentage.png) 
*Lollipop plot displaying gene percentage distribution*

**Geographic Distribution Heatmap:**
![Geographic Heatmap](ncbi/figures/heatmap/ncbi_ncbi_Continent_heatmap.png)
*Heatmap showing resistance gene distribution across continents*

**Antibiotic Resistance by Classes:**
![Resistance by Antibiotics](ncbi/figures/Resistance_percentage_by_Antibiotics.png)
*Bar plot showing resistance patterns by antibiotic classes*

**Gene Identity Analysis:**
![Gene Identity Boxplot](ncbi/figures/Resistance_gene_identity_boxplot.png)
*Boxplot analysis of resistance gene identity scores*

### Sample Output Directory Structure
```
ncbi/figures/
├── Resistance_gene_presence.png
├── Resistance_gene_percentage.png  
├── Resistance_gene_identity_boxplot.png
├── Resistance_percentage_by_Antibiotics.png
├── heatmap/
│   ├── ncbi_ncbi_Continent_heatmap.png
│   ├── ncbi_ncbi_Geographic_Location_heatmap.png
│   ├── ncbi_ncbi_Subcontinent_heatmap.png
│   └── ncbi_ncbi_Collection_Date_heatmap.png
├── mean_ARG/
│   ├── Mean_ARG_by_Continent.png
│   ├── Mean_ARG_by_Geographic Location.png
│   ├── Mean_ARG_by_Subcontinent.png
│   └── Mean_ARG_by_Collection Date.png
└── html_files/
    ├── index.html (Main navigation page)
    └── [Interactive HTML plots - open in browser]
```

### Interactive Features
The tool generates interactive HTML visualizations that provide enhanced data exploration capabilities:

- **Dynamic Interaction**: Zooming, panning, and selection tools
- **Detailed Tooltips**: Hover for comprehensive data information  
- **Filtering Options**: Dynamic data filtering and highlighting
- **Export Capabilities**: High-quality image export functionality
- **Responsive Design**: Optimized for different screen sizes

**Available Interactive Visualizations:**
- Geographic distribution maps with zoom capabilities
- Interactive heatmaps with data filtering
- Dynamic correlation plots with hover details
- Responsive boxplots and bar charts
- Time-series analysis with date selection

## Statistical Analysis Features

PanR2 provides comprehensive statistical analysis including:
- **Correlation Analysis**: Geographic and temporal correlations
- **Comparative Statistics**: Between-group comparisons
- **Summary Statistics**: Descriptive statistics for all variables
- **Pairwise Comparisons**: Detailed pairwise statistical tests
- **Geographic Patterns**: Spatial distribution analysis
---

## Requirements

- Python 3.9+
- Required Python packages (automatically installed):
  - pandas
  - numpy
  - matplotlib
  - seaborn
  - plotly
  - scipy
  - Other dependencies listed in `requirements.txt`

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request or open an Issue for bugs, feature requests, or improvements.

---

## License

This tool is provided under the MIT License. See `LICENSE` file for details.

---

## Citation

If you use PanR2 in your research, please cite: If you use PanR2 in your research, please cite: DOI: 10.1101/2025.04.08.647722 

## Support

For questions, issues, or feature requests, please:
1. Check the existing [Issues](https://github.com/Tasnimul-Arabi-Anik/PanR2/issues)
2. Create a new issue with detailed information
3. Contact the author: Tasnimul Arabi Anik (arabianik987@gmail.com)

---
