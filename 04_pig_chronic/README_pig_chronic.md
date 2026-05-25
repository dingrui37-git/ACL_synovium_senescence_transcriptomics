# Chronic pig extension module

**GitHub module folder:** `04_pig_chronic/`

**Original local package folder:** `UPLOAD_minimal_chronic_figure5_submission_with_PNOM_20260525_104430`

This folder contains the chronic pig extension analysis package supporting Figure 5 and related chronic-stage supplementary audit outputs for the manuscript:

**Cross-species transcriptomic analysis identifies a persistent senescence-associated synovial remodeling program after anterior cruciate ligament injury**

## Purpose

This module provides the chronic pig synovium reproducibility package. It includes scripts, processed matrices, audit files, figure source data, final figure panels, logs, version records, and manifests used for the chronic pig extension analysis.

The original public sequencing dataset is not redistributed here. The chronic pig analysis was based on public Salmon `quant.sf` processed quantification files from **GSE228848**, summarized and audited through the scripts and processed outputs retained in this module.

## Original package-generation metadata

- Created at: 2026-05-25 10:44:51 CEST
- Original project root: `E:/R/ACLsenescence2`
- Original package root: `E:/R/ACLsenescence2/rebuild_submission/UPLOAD_minimal_chronic_figure5_submission_with_PNOM_20260525_104430`

These original local paths are retained for provenance only. In the GitHub repository, this module is stored under:

```text
04_pig_chronic/
```

## Locked chronic branch logic

The chronic pig extension follows the locked analysis/order logic below:

```text
B  Step23 quant.sf-level QC
C  Step24 tximport from confirmed quant files
P  StepCHRONQC_01 Salmon/tximport quantification audit supplement
D  Step24B clean tximport matrix audit
N  StepCHRONLOCK_01 synovium-only 48 / main-24 matrix lock supplement
O  StepCHRONMDS_01 chronic main-24 PCA/MDS supplement
E  Step25v2 synovium manifest/subset rebuild
F  Step25v3 final group-label rebuild and main-comparison matrices
A  Step22 methods/data-source evidence archive
G  Chronic Step2 current78 TMM-aligned signature score / Figure 5A
M  Chronic Step3 current78 TMM-aligned singscore sensitivity supplement
H  Chronic Step4 edgeR QLF DE + full Hallmark GSEA / Figure 5B
J  Locked Figure 5C core-gene chronic DE audit plot
K  Figure 5D early-versus-chronic summary panel
```

## Main Figure 5 panel policy

The locked Figure 5 panel policy is:

- **Figure 5A:** G TMM-aligned chronic signature score figure.
- **Figure 5B:** H target Hallmark GSEA running curves.
- **Figure 5C:** J Step4-locked core-gene chronic DE audit lollipop plot.
- **Figure 5D:** K early-versus-chronic summary tile plot.

## Important consistency notes

- Uploaded L reads the old non-TMM folder and is **not required** in this chronic package.
- Uploaded K reads `tables/chronic_step5B_current78_core_DE_audit/`.
- Uploaded J reads `tables/chronic_step5B_current78_core_DE_audit_Step4_locked/`.
- Both exact paths are retained in the package so this historical mismatch remains visible in the manifest.
- P/N/O/M were added as required supplementary outputs with exact source paths.
- N/O are retained as supplementary matrix-lock and PCA/MDS branches and do not replace Step25v3 unless explicitly stated in the manuscript.

## Recommended audit before upload or review

Before creating a GitHub release or Zenodo archive, check:

1. `copy_manifest_exact.csv`
2. Confirm all required files have `source_exists = TRUE` and `copied = TRUE`.
3. `missing_required_files.csv`
4. Confirm `missing_required_files.csv` has zero rows.

Original package summary:

- Number of manifest rows: 223
- Copied files: 162
- Missing required files: 0
- Missing optional files: 24

Missing optional files are retained as an audit record and do not indicate incompleteness of the required chronic submission package.

## Large-file and GitHub upload notes

Large public or regenerable files are not redistributed in this GitHub module unless needed as processed source data. Where large CSV matrices are required, compressed `.csv.gz` files should be used where possible. If any individual processed matrix remains too large for GitHub web upload, it should be deposited in the associated Zenodo archive and documented in the root README or a local `README_large_matrices.md`.

## Relationship to the root README

The repository-level `README.md` provides the overall project structure, public dataset accessions, general methodological anchors, and cross-module reproducibility notes. This module-level README provides chronic pig-specific analysis order, Figure 5 panel policy, and chronic branch audit notes.

