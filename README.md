# PanR2: Panresistome Analysis Tool

## Overview
PanR2 is a comprehensive Python-based tool for analyzing panresistome data (Global Resistance Pattern). It processes NCBI and Abricate summary files, merges the data, and generates a wide range of visualizations including heatmaps, bar plots, boxplots, and interactive HTML plots. The tool is designed to help researchers analyze and visualize antibiotic resistance gene presence, prevalence, and distribution patterns across different geographic locations and temporal scales.

It enables robust statistical analysis, including group-wise comparisons, summary statistics, and correlation assessments, to support meaningful interpretation of resistome data. Given sufficient sequencing data, PanR2 can help identify local resistance spread, compare resistance patterns between regions or over time, and detect signals of emerging resistance epidemics. This makes it a valuable tool for surveillance, epidemiology, and public health research.

**Prerequisites:**
- `ncbi_clean.csv` from [FetchM](https://github.com/Tasnimul-Arabi-Anik/FetchM)
- Summary files in `.tab` (preferred) or `.csv` format from [Abricate](https://github.com/tseemann/abricate)

### Key Features:
- Merges and analyzes NCBI and Abricate outputs
- Calculates gene presence/absence across multiple categories (continent, location, subcontinent, collection date)
- Performs comprehensive statistical tests and correlation analyses
- Generates multiple visualization types: heatmaps, bar plots, boxplots, lollipop plots, and correlation plots
- Creates interactive HTML visualizations for enhanced data exploration
- Generates an interactive HTML index for easy navigation of all results
- Provides detailed statistical analysis outputs

---

## Installation

### Method 1: Using `pip` with `conda` (Recommended)
```bash
conda create -n panr_env python=3.9
conda activate panr_env
pip install panR2
```

### Method 2: Direct installation from GitHub
```bash
conda create -n panr_env python=3.8
conda activate panr_env
pip install git+https://github.com/Tasnimul-Arabi-Anik/PanR2.git
```

### Method 3: Manual Installation from Source
```bash
git clone https://github.com/Tasnimul-Arabi-Anik/PanR2.git
cd PanR2
pip install -r requirements.txt
```

### Confirm Installation
```bash
panr --help
```

### Run the Included Smoke Test
From a source checkout, run the bundled small fixture to confirm the CLI, merge step, tidy output, and plotting path work locally:

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

PanR2 keeps the installed command as `panr`, while the implementation lives in the `panr2` Python package. The `bin/panr` file is a compatibility wrapper around `panr2.cli`. Core code is split into focused modules: `io`, `qc`, `filters`, `analysis`, `plots`, and `stats`.

---

## Usage

### Command-Line Interface
```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

PanR2 can also run ABRicate internally from assembly FASTA files:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-abricate --abricate-dbs ncbi,vfdb,plasmidfinder --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

MobileElementFinder can also be run internally when installed:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-mobileelementfinder --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

IntegronFinder can also be run internally when installed:

```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --sequence-dir <SEQUENCE_DIRECTORY> --run-integronfinder --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

### Required Arguments
| Argument | Description |
|----------|-------------|
| `--ncbi-dir` | Directory containing `ncbi_clean.csv` from FetchM |
| `--output-dir` | Directory to store merged results and visualizations |

`--abricate-dir` is required for analysis-only mode. It is not required when `--run-abricate` is used with the `ncbi` database.

### Optional Arguments
| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--genep` | float | `10.0` | Minimum % gene presence to include in heatmap |
| `--nseq` | int | `1` | Minimum number of sequences required per group in heatmaps |
| `--format` | str | `tiff` | Output format for figures (`tiff`, `svg`, `png`, `pdf`) |
| `--sequence-dir` | path | optional | Directory containing assembly FASTA files for integrated tool runners |
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
| `--force-tool-run` | flag | off | Re-run integrated tools even when raw result files already exist |
| `--min-identity` | float | `0.0` | Minimum ABRicate identity percentage to treat a gene call as present |
| `--drop-unmatched-accessions` | flag | off | Drop NCBI rows with no matching ABRicate summary row |
| `--min-samples-per-group` | int | `5` | Minimum samples per group required for correlation analyses |
| `--core-threshold` | float | `95.0` | Prevalence percentage used to classify core ARGs |
| `--rare-threshold` | float | `5.0` | Prevalence percentage used to classify rare ARGs |
| `--top-n` | int | `25` | Number of top genes/classes to include in compact summary plots |
| `--cooccurrence-min-prevalence` | float | `0.0` | Minimum prevalence percentage for genes/classes included in co-occurrence matrices |
| `--cooccurrence-top-n` | int | `25` | Number of top genes/classes or pairs to include in co-occurrence plots and pair tables |
| `--vfdb-dir` | path | optional | Directory containing ABRicate VFDB summary/results files |
| `--plasmidfinder-dir` | path | optional | Directory containing ABRicate PlasmidFinder summary/results files |
| `--mobileelementfinder-dir` | path | optional | Directory containing ABRicate MobileElementFinder summary/results files |
| `--isfinder-dir` | path | optional | Directory containing ABRicate ISfinder summary/results files |
| `--integronfinder-dir` | path | optional | Directory containing IntegronFinder or ABRicate-style integron summary/results files |
| `--iceberg-dir` | path | optional | Directory containing ABRicate ICEberg summary/results files |
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
├── qc/                                # Shared input validation and filter reports
└── report/                            # Journal-style narrative report and methods text
```

### Output Files Description

#### 1. Input QC (`qc/` directory)
- **`panr2_input_qc.csv`** - Machine-readable input checks with `PASS`, `WARN`, `FAIL`, and `INFO` statuses
- **`panr2_input_qc_summary.txt`** - Human-readable summary of input checks
- **`panr2_tool_manifest.csv`** and **`panr2_tool_manifest.json`** - Tool versions, selected databases, database dates/counts where available, and raw output paths for integrated runs
- **`panr2_unmatched_accessions.csv`** - Accessions present in only one of the NCBI or ABRicate inputs
- **`*_filter_report.csv`** - Row and ARG-call counts before and after optional analysis filters

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
- **`plasmidfinder/analysis/plasmidfinder_category_summary.csv`** - Plasmid replicon/category summaries
- **`plasmidfinder/analysis/plasmidfinder_geographic_summary.csv`** - Descriptive plasmid replicon burden by metadata group
- **`plasmidfinder/analysis/plasmidfinder_temporal_summary.csv`** - Descriptive plasmid replicon burden by collection year
- **`plasmidfinder/analysis/plasmidfinder_feature_cooccurrence_matrix.csv`** - Plasmid replicon co-occurrence counts across samples
- **`plasmidfinder/analysis/plasmidfinder_top_feature_pairs.csv`** - Ranked plasmid replicon pairs with support and Jaccard index
- **`plasmidfinder/analysis/plasmidfinder_group_burden_summary.csv`** - PlasmidFinder feature burden summaries by location, continent, subcontinent, and year
- **`plasmidfinder/analysis/plasmidfinder_group_overall_tests.csv`** - PlasmidFinder nonparametric group-level tests where sample sizes permit
- **`plasmidfinder/analysis/plasmidfinder_group_pairwise_tests.csv`** - PlasmidFinder pairwise group comparisons where sample sizes permit
- **`plasmidfinder/figures/`** - PlasmidFinder static figures, burden-by-group plots, interactive HTML files under `html_files/`, and `index.html` navigation

The same database-named output pattern is used for `mobileelementfinder/`, `isfinder/`, `integronfinder/`, and `iceberg/`. Each optional mobile genetic element database writes `analysis/`, `figures/`, and `merged_output/` folders with feature summaries, product/category summaries, sample burden, geographic and temporal summaries, co-occurrence tables, group-burden comparisons, static plots, interactive HTML plots, and an `index.html` figure browser.

VFDB, PlasmidFinder, MobileElementFinder, ISfinder, IntegronFinder, and ICEberg are handled as separate feature families, not as antibiotic resistance classes. PanR2 therefore uses feature prevalence, product/replicon categories, sample burden, geography and temporal summaries, identity distributions, feature co-occurrence, feature presence heatmaps, and interactive HTML figures instead of resistance-class composition plots.

#### 4. Written Report (`report/` directory)
- **`*_panr2_report.md`** - Comprehensive journal-style narrative report generated from output tables
- **`*_panr2_report.html`** - Simple HTML rendering of the Markdown report
- **`*_methods.txt`** - Reusable methods description for manuscript drafting

#### 5. NCBI/AMR Static Visualizations
- **`Resistance_gene_presence.{format}`** - Bar plot showing gene presence across samples
- **`Resistance_gene_percentage.{format}`** - Lollipop plot showing gene percentage distribution
- **`Resistance_gene_identity_boxplot.{format}`** - Boxplot of resistance gene variation across sequences
- **`Resistance_percentage_by_Antibiotics.{format}`** - Bar plot of resistance by antibiotic classes

#### 6. Heatmaps (`ncbi/figures/heatmap/` directory)
- **`Resistance gene distribution by continent, geographic location, subcontinent, and year.`**

#### 7. Mean ARG Analysis (`ncbi/figures/mean_ARG/` directory)
- **`Average antibiotic resistance genes by continent, geographic location, subcontinent, and year.`** - 

#### 8. Interactive HTML Visualizations (`ncbi/figures/html_files/` directory)
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

#### 9. Statistical Analysis (`ncbi/figures/Stat_analysis/` directory)
- **`combined_geographic_correlation_summary.csv`** - Geographic correlation statistics
- **`combined_overall_tests.csv`** - Overall statistical test results
- **`combined_pairwise_comparisons.csv`** - Pairwise comparison results
- **`combined_summary_statistics.csv`** - Comprehensive summary statistics
- **`ncbi_gene_presence_count_percentage.csv`** - Gene presence counts and percentages

#### 10. Navigation
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
