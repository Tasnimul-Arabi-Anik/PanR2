# Changelog

## Unreleased

### Added
- Reproducible install support with `environment.yml`, `Dockerfile`, and `.dockerignore` for source and container-based PanR2 setups.
- `panr --doctor` environment check reporting Python dependencies, optional external annotation tools, ABRicate database visibility, and supported install modes.
- Per-database optional feature QC summaries and unmatched-sample reports for VFDB, PlasmidFinder, MobileElementFinder, ISfinder, IntegronFinder, and ICEberg-style analyses.
- Explicit assembly accession normalization helpers that preserve GCF/GCA version suffixes by default.
- Full multi-database workflow validation tests covering NCBI AMR, VFDB, PlasmidFinder, MobileElementFinder, ISfinder, IntegronFinder, ICEberg, QC outputs, reports, static figures, and interactive HTML outputs.
- Larger 10-sample synthetic multi-database integration fixture to exercise group-level summaries and continent-level analyses without requiring external tools or large datasets.
- ICEberg table-converter workflow for user-provided ICE/IME/CIME CSV/TSV/TAB annotations with PanR2-compatible feature table conversion and manifest tracking.
- Integrated IntegronFinder runner mode with raw output preservation, PanR2-compatible feature table conversion, and manifest tracking.
- Integrated MobileElementFinder runner mode with raw CSV preservation, PanR2-compatible feature table conversion, and manifest tracking.
- Integrated ABRicate runner mode with `--run-abricate`, `--sequence-dir`, database selection, reusable raw outputs, and tool/database manifest tracking.
- Optional MobileElementFinder, ISfinder, IntegronFinder, and ICEberg ABRicate-style feature analyses with database-named outputs, figures, HTML indexes, and report sections.
- VFDB/PlasmidFinder group-burden summaries, burden-by-metadata plots, and nonparametric group-comparison outputs.
- VFDB/PlasmidFinder feature co-occurrence, identity distribution plots, and temporal feature-burden summaries.
- Optional VFDB and PlasmidFinder ABRicate-style feature analysis with separate virulence/plasmid summaries and report sections.
- Deterministic Markdown/HTML report generation with methods and reproducibility summaries.
- Unit tests for filtering, comprehensive summaries, and co-occurrence outputs.
- Gene and resistance-class co-occurrence matrices, ranked pair tables, and compact co-occurrence heatmaps.
- Comprehensive panresistome analysis tables for sample burden, gene prevalence, core/accessory/rare categories, and resistance class summaries.
- Compact summary plots for top ARG prevalence and resistance class prevalence.
- Analysis filtering options: `--min-identity`, `--drop-unmatched-accessions`, and `--min-samples-per-group`.
- Per-run filter reports under `output/qc/`.
- Input QC reports under `output/qc/` for required columns, accession matching, zero-hit samples, and resistance annotations.
- Lightweight CI validation for the PanR2 CLI.
- Small ABRicate/NCBI fixtures for end-to-end smoke testing.

### Changed
- IntegronFinder table selection now prefers detailed `.integrons` outputs, then tabular files, and uses `.summary` files only as a logged fallback.
- Documentation and reports now describe ICEberg support as an ICE/IME/CIME table-conversion workflow rather than a direct ICEberg runner.
- VFDB and PlasmidFinder outputs now mirror the `ncbi/` folder shape with `analysis/`, `figures/`, `merged_output/`, and interactive HTML figure indexes.
- AMR/resistome outputs now live under `ncbi/`, leaving top-level output organized as `ncbi/`, `vfdb/`, `plasmidfinder/`, `qc/`, and `report/`.
- Optional VFDB and PlasmidFinder outputs now use database-named folders (`vfdb/` and `plasmidfinder/`).
- Split plotting and statistical helper code into `panr2.plots` and `panr2.stats`.
- Split table I/O, QC, filtering, and comprehensive analysis code into dedicated `panr2` modules.
- Moved the implementation into the `panr2` package while keeping `bin/panr` as a compatibility wrapper.
- Improve ABRicate accession parsing to support both GCF and GCA assembly accessions.
- Report missing required inputs as command failures instead of successful no-op runs.
- Avoid depending on global CLI arguments inside `main()`.
- Make tidy conversion identify gene columns after `NUM_FOUND` when available.
