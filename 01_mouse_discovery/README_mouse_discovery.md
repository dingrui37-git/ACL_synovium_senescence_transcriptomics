# Mouse discovery module

**GitHub module folder:** `01_mouse_discovery/`

**Original local package folder:** `ACL_synovium_senescence_minimal_submission_package_20260523_195310`

This folder contains the mouse discovery analysis package for the manuscript:

**Cross-species transcriptomic analysis identifies a persistent senescence-associated synovial remodeling program after anterior cruciate ligament injury**

## Purpose

This module provides the mouse discovery reproducibility package. It includes final R scripts, figure source data, key processed result tables, final figures, logs, version records, and public raw-data accession records.

The package is intended to help reviewers and readers verify the source data behind the mouse discovery figures and rerun the mouse analysis from public data where applicable.

## Public dataset and reference resource

- **Mouse discovery cohort:** GSE271903
- **Senescence reference:** CellAge `cellage3.tsv`

Raw-data provenance is documented in `public_dataset_accessions.csv` and/or `raw_data_manifest/README_raw_data.md`, if present.

## Main contents

- `scripts/`: final manuscript R scripts for the mouse discovery branch.
- `source_data/`: source data used for mouse-related manuscript figure panels.
- `tables/`: processed result tables, including differential expression, strict DEG, persistent gene, CellAge, ortholog mapping, and GSEA outputs.
- `figures/`: final exported mouse discovery and CellAge-related figure files.
- `logs_versions/`: run logs, `sessionInfo`, software/package version files, fgsea settings, and gprofiler/gorth audit files where available.
- `external_data/CellAge/`: CellAge raw table if locally copied; otherwise, use the cleaned/audited CellAge files under `tables/` or `source_data/`.

## Recommended rerun order

```text
1.  step01_initialize_project_structure.R
2.  step02_download_GSE271903.R
3.  step03_load_expression.R
4.  step04_construct_paired_samples.R
5.  step05_mouse_DE_limma_voom_volcano.R
6.  step06_mouse_PCA_from_step05_voom.R
7.  step07_strict_DEG_upset_persistent.R
8.  step08_logFC_scatter.R
9.  step09_persistent_heatmap.R
10. step10_strict_gene_set_summary.R
11. step11A_CellAge_raw_to_clean.R
12. step11B_current_gorth_strict_1to1_remap.R
13. step11C_intersect_current_gorth_mapping_with_CellAge_clean.R
14. step12A_whole_CellAge_projected_score.R
15. step13_mouse_Hallmark_GSEA.R
16. step14_Figure3C_mapping_framework_barplot.R
17. step15_Figure3D_top_CellAge_overlap_genes.R
```

## Important methodological anchors

- Mouse differential expression uses a paired limma-voom workflow with `duplicateCorrelation`.
- Mouse PCA uses the second voom-transformed logCPM matrix from the locked differential-expression workflow.
- Strict DEGs are defined as genes with FDR < 0.05 and |logFC| > 1.
- Persistent genes are defined as genes that are strict at both early time points and show consistent logFC direction.
- CellAge and ortholog analyses use cleaned CellAge gene symbols and strict one-to-one gorth ortholog mapping.
- Mouse Hallmark GSEA uses the full gene-level ranked list with rank statistic `sign(logFC) * -log10(P value)`.
- fgsea settings and gprofiler/gorth audit files should be checked under `logs_versions/` when available.

## Package audit

Before GitHub release or Zenodo archiving, check:

- `MANIFEST_all_copied_files.csv`: every copied/generated file with target path, size, and MD5.
- `MANIFEST_missing_expected_files.csv`: expected files or folders that were not found.

Generated at: 2026-05-23 19:53:17 CEST  
Original source mouse base directory: `E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery`

## Relationship to the root README

The repository-level `README.md` provides the overall project structure and cross-module reproducibility notes. This module-level README documents the mouse discovery branch specifically.

