# PanR2: Panresistome Analysis Tool

## Overview
PanR2 is a comprehensive Python-based tool for analyzing panresistome data. It processes NCBI and Abricate summary files, merges the data, and generates a wide range of visualizations including heatmaps, bar plots, boxplots, and interactive HTML plots. The tool is designed to help researchers analyze and visualize antibiotic resistance gene presence, prevalence, and distribution patterns across different geographic locations and temporal scales.

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

---

## Usage

### Command-Line Interface
```bash
panr --ncbi-dir <NCBI_DIRECTORY> --abricate-dir <ABRICATE_DIRECTORY> --output-dir <OUTPUT_DIRECTORY> [OPTIONS]
```

### Required Arguments
| Argument | Description |
|----------|-------------|
| `--ncbi-dir` | Path to the `ncbi_clean.csv` file from FetchM |
| `--abricate-dir` | Directory containing Abricate summary `.tab` or `.csv` files |
| `--output-dir` | Directory to store merged results and visualizations |

### Optional Arguments
| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--genep` | float | - | Minimum % gene presence to include in heatmap |
| `--nseq` | int | - | Minimum number of sequences required per group in heatmaps |
| `--format` | str | `png` | Output format for figures (`tiff`, `svg`, `png`, `pdf`) |
| `--version` | - | - | Show program's version number and exit |
| `-h, --help` | - | - | Show help message and exit |

### Example Usage
```bash
# Basic usage
panr --ncbi-dir ./data/ncbi_clean.csv --abricate-dir ./data/abricate --output-dir ./output

# With optional parameters
panr --ncbi-dir ./data/ncbi_clean.csv --abricate-dir ./data/abricate --output-dir ./output --format pdf --genep 10 --nseq 5
```

---

## Output Structure

PanR2 generates a comprehensive set of outputs organized in the following directory structure:

```
output/
├── figures/
│   ├── heatmap/                          # Geographic heatmaps
│   ├── html_files/                       # Interactive HTML plots
│   ├── mean_ARG/                         # Mean antibiotic resistance gene plots
│   ├── Stat_analysis/                    # Statistical analysis results
│   ├── index.html                        # Main navigation page
│   └── [Various static plots]
└── [Processed CSV files]
```

### Output Files Description

#### 1. Static Visualizations
- **`Resistance_gene_presence.{format}`** - Bar plot showing gene presence across samples
- **`Resistance_gene_percentage.{format}`** - Lollipop plot showing gene percentage distribution
- **`Resistance_gene_identity_boxplot.{format}`** - Boxplot of resistance gene identity scores
- **`Resistance_percentage_by_Antibiotics.{format}`** - Bar plot of resistance by antibiotic classes

#### 2. Heatmaps (`heatmap/` directory)
- **`ncbi_ncbi_Continent_heatmap.{format}`** - Resistance gene distribution by continent
- **`ncbi_ncbi_Geographic_Location_heatmap.{format}`** - Distribution by geographic location
- **`ncbi_ncbi_Subcontinent_heatmap.{format}`** - Distribution by subcontinent
- **`ncbi_ncbi_Collection_Date_heatmap.{format}`** - Temporal distribution patterns

#### 3. Mean ARG Analysis (`mean_ARG/` directory)
- **`Mean_ARG_by_Continent.{format}`** - Average antibiotic resistance genes by continent
- **`Mean_ARG_by_Geographic Location.{format}`** - Average ARGs by geographic location
- **`Mean_ARG_by_Subcontinent.{format}`** - Average ARGs by subcontinent  
- **`Mean_ARG_by_Collection Date.{format}`** - Temporal trends in ARG abundance

#### 4. Interactive HTML Visualizations (`html_files/` directory)
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

#### 5. Statistical Analysis (`Stat_analysis/` directory)
- **`combined_geographic_correlation_summary.csv`** - Geographic correlation statistics
- **`combined_overall_tests.csv`** - Overall statistical test results
- **`combined_pairwise_comparisons.csv`** - Pairwise comparison results
- **`combined_summary_statistics.csv`** - Comprehensive summary statistics
- **`ncbi_gene_presence_count_percentage.csv`** - Gene presence counts and percentages

#### 6. Navigation
- **`index.html`** - Interactive HTML index page for easy navigation of all generated visualizations

---

## Example Visualizations

### Static Plots

**Resistance Gene Presence Analysis:**
![Resistance Gene Presence](figures/Resistance_gene_presence.png)
*Bar plot showing the presence of resistance genes across samples*

**Gene Percentage Distribution:**
![Resistance Gene Percentage](figures/Resistance_gene_percentage.png) 
*Lollipop plot displaying gene percentage distribution*

**Geographic Distribution Heatmap:**
![Geographic Heatmap](figures/heatmap/ncbi_ncbi_Continent_heatmap.png)
*Heatmap showing resistance gene distribution across continents*

**Antibiotic Resistance by Classes:**
![Resistance by Antibiotics](figures/Resistance_percentage_by_Antibiotics.png)
*Bar plot showing resistance patterns by antibiotic classes*

**Gene Identity Analysis:**
![Gene Identity Boxplot](figures/Resistance_gene_identity_boxplot.png)
*Boxplot analysis of resistance gene identity scores*

### Sample Output Directory Structure
```
figures/
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

**To Access Interactive Plots:**
1. Navigate to your output directory
2. Open `figures/index.html` in a web browser
3. Click on any visualization link to explore interactive plots

**Available Interactive Visualizations:**
- Geographic distribution maps with zoom capabilities
- Interactive heatmaps with data filtering
- Dynamic correlation plots with hover details
- Responsive boxplots and bar charts
- Time-series analysis with date selection

> **Note**: Interactive HTML files must be opened in a web browser to view. GitHub README cannot display interactive HTML content directly.

---

## Statistical Analysis Features

PanR2 provides comprehensive statistical analysis including:
- **Correlation Analysis**: Geographic and temporal correlations
- **Comparative Statistics**: Between-group comparisons
- **Summary Statistics**: Descriptive statistics for all variables
- **Pairwise Comparisons**: Detailed pairwise statistical tests
- **Geographic Patterns**: Spatial distribution analysis

---

## Requirements

- Python 3.8+
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
3. Contact the author: Tasnimul Arabi Anik

---

## Changelog

### Latest Updates
- Added interactive HTML visualizations
- Enhanced statistical analysis capabilities  
- Improved output organization and navigation
- Added support for multiple figure formats
- Enhanced correlation analysis features