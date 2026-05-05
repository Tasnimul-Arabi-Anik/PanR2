# PanR2 Input Contract

PanR2 is the lightweight comparative analysis and reporting layer. Heavy external tools and databases should run in PanResistome, Nextflow, Snakemake, shell scripts, HPC jobs, or other user workflows, then be exported as standardized PanR2-compatible tables.

## Required Feature Columns

Every feature-like tool output should be convertible to these columns:

| column | description |
| --- | --- |
| `sample_id` | Local sample name or FASTA stem used by the upstream tool |
| `assembly_accession` | Stable assembly accession, preferably versioned `GCF_*.*` or `GCA_*.*` |
| `database` | Feature family or database name, such as `ncbi`, `vfdb`, `plasmidfinder`, `mobileelementfinder`, `isfinder`, `integronfinder`, `iceberg`, `mlst`, `defensefinder`, `prophage`, `ani`, or `assembly_qc` |
| `feature_id` | Gene, replicon, element, ST, cluster, or metric identifier |
| `feature_category` | Biological category, product, class, system type, or summary category |
| `presence` | `1` for detected/present feature; `0` can be used for explicit absence records |
| `identity` | Percent identity or tool-specific primary score when available |
| `coverage` | Percent coverage or tool-specific secondary score when available |
| `contig` | Contig or sequence identifier when available |
| `start` | Feature start coordinate when available |
| `end` | Feature end coordinate when available |
| `tool` | Tool that produced the call |
| `tool_version` | Tool version, if captured |
| `database_version` | Database version/date/path, if captured |

## Design Rule

PanResistome should run heavy tools and export PanR2-ready tables. PanR2 should read these tables, merge them with metadata, and generate statistics, figures, dashboards, citations, and manuscript-style reports without requiring the external tools to be installed.

This split keeps the analysis layer installable with Python while still allowing comprehensive genomics workflows when users run the full PanResistome execution layer.

## Metadata Input

PanR2 accepts the legacy `ncbi_clean.csv` name for compatibility, but the preferred upstream source is FetchM2 through PanResistome. When FetchM2 columns are present, PanR2 preserves and analyzes richer standardized fields such as `Country`, `Continent`, `Subcontinent`, `Collection_Year`, `Host_SD`, host taxonomy columns, `Sample_Type_SD`, `Isolation_Source_SD`, `Environment_Medium_SD`, `Environment_Broad_Scale_SD`, `Environment_Local_Scale_SD`, `Host_Disease_SD`, and metadata audit/status fields.

PanR2 also fills older aliases such as `Geographic Location`, `Collection Date`, `Host`, and `Isolation Source` from FetchM2 fields so existing workflows and figures remain compatible.
