# Early pig validation module

**GitHub module folder:** `02_pig_early/`

**Original local package folder:** `UPLOAD_minimal_pig_early_submission_20260524`

This folder contains the final early pig validation analysis package for the manuscript:

**Cross-species transcriptomic analysis identifies a persistent senescence-associated synovial remodeling program after anterior cruciate ligament injury**

## Purpose

This module is designed for manuscript submission, GitHub upload, and Zenodo archiving of the final early pig validation branch.

It retains the scripts, processed data, figure source files, result tables, QC summaries, and logs required to document and audit the early pig validation analysis.

## Public dataset

- **Early pig validation cohort:** E-MTAB-6664

## Main contents

- Final R scripts for the early pig branch.
- Sample manifests and QC summaries.
- Reference/index manifests.
- featureCounts gene-level count matrix and assignment summaries.
- edgeR QLF differential expression outputs.
- Current signature and ortholog validation tables.
- TMM-logCPM signature score source data.
- Full Hallmark GSEA outputs.
- Figure 4 source data and final figure panels.
- Early pig PCA/MDS supplementary figures, source data, method parameters, and logs.
- Singscore sensitivity outputs.
- Logs, software/package version records, and manifest files.

## Files intentionally excluded

The following large or regenerable files are intentionally excluded from the GitHub upload package:

- raw FASTQ files
- BAM files
- reference FASTA/GTF files
- Subread index files
- large `.RData` workspaces

These files are public, regenerable, or too large for routine GitHub upload. Their provenance and regeneration routes should be documented through manifests and scripts where applicable.

## Important final-vs-old decisions

- Keep the Step20 symbol-only v4 zero-row-fix script as the final GSEA implementation.
- Do not present `Step20_current78_pig_early_hallmark_gsea_add_MSigDB_db_version.R` as the final GSEA script.
- Treat any uploaded file name containing `old step16` carefully: by content, it may be the Step16 v3 featureCounts script and should be retained only if it is the full final featureCounts script.

## Recommended audit before release

Before GitHub release or Zenodo archiving, inspect:

- `07_logs_versions_and_manifest/copied_file_manifest.csv`
- `07_logs_versions_and_manifest/missing_required_or_optional_files.csv`
- `07_logs_versions_and_manifest/script_qc_report.csv`

Any row marked `required = TRUE` and `copied = FALSE` should be resolved before release.

## PCA/MDS note

This package checks multiple possible early-pig PCA/MDS output paths and prioritizes:

```text
rebuild_submission/02_pig_early/supplement/MDS_PCA
```

It also remains compatible with the path written in the uploaded PCA/MDS script:

```text
rebuild_submission/supplementary files/MDS_PCA
```

If PCA/MDS files are missing from all checked paths, run:

```text
StepPCA_MDS_02_pig_early_unified_limma_MDS.R
```

and then rerun the organizer.

## Relationship to the root README

The repository-level `README.md` provides the overall project structure and cross-module reproducibility notes. This module-level README documents the early pig validation branch specifically.
