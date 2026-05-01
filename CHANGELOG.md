# Changelog

## Unreleased

### Added
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
- Moved the implementation into the `panr2` package while keeping `bin/panr` as a compatibility wrapper.
- Improve ABRicate accession parsing to support both GCF and GCA assembly accessions.
- Report missing required inputs as command failures instead of successful no-op runs.
- Avoid depending on global CLI arguments inside `main()`.
- Make tidy conversion identify gene columns after `NUM_FOUND` when available.
