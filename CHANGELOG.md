# Changelog

## Unreleased

### Added
- Lightweight CI validation for the PanR2 CLI.
- Small ABRicate/NCBI fixtures for end-to-end smoke testing.

### Changed
- Improve ABRicate accession parsing to support both GCF and GCA assembly accessions.
- Report missing required inputs as command failures instead of successful no-op runs.
- Avoid depending on global CLI arguments inside `main()`.
- Make tidy conversion identify gene columns after `NUM_FOUND` when available.
