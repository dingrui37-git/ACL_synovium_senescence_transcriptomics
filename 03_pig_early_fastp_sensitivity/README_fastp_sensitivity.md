# Early pig fastp sensitivity module

**GitHub module folder:** `03_pig_early_fastp_sensitivity/`

**Original local package folder:** `H2w_fastp_minpkg_20260524_223735`

This folder contains the fastp-trimmed sensitivity analysis package for the early pig validation branch of the manuscript:

**Cross-species transcriptomic analysis identifies a persistent senescence-associated synovial remodeling program after anterior cruciate ligament injury**

## Purpose

This module documents the sensitivity analysis comparing the primary untrimmed early pig workflow with a fastp-trimmed branch. It is intended to support the robustness of the early pig differential expression, signature scoring, GSEA, and core-gene results.

## Public dataset

- **Early pig validation cohort:** E-MTAB-6664

## Original package-generation metadata

- Generated: 2026-05-24 22:37:39 CEST
- Original project root: `E:/R/ACLsenescence2`
- Original package root: `E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/10_submission/H2w_fastp_minpkg_20260524_223735`

These original local paths are retained for provenance only. In the GitHub repository, this module is stored under:

```text
03_pig_early_fastp_sensitivity/
```

## Main contents

Included files and outputs cover:

1. H2w1–H2w13 scripts, with exact A–M-prefixed filename alternatives supported.
2. Primary untrimmed inputs used by sensitivity comparisons.
3. Fastp QC summaries.
4. Alignment and featureCounts summaries.
5. Trimmed count matrix and differential expression outputs.
6. Differential-expression stability summaries.
7. Signature-score stability summaries.
8. GSEA stability summaries.
9. Core-gene audit source data.
10. Final supplementary figure PDFs/PNGs and corresponding source data/statistics.
11. Method, version, `sessionInfo`, and log records where present.

## Files intentionally excluded

The following large or regenerable files are intentionally excluded from this module:

- raw FASTQ files
- fastp-trimmed FASTQ files
- BAM files
- BAM summary files
- full per-sample fastp HTML reports
- full per-sample fastp JSON reports

## Audit and copy logic

The copy manifest records every exact source path and package destination. No fuzzy search was used.

This v4 package keeps the short package root name and fixes copy-log column consistency before `rbind`.

## Relationship to the root README

The repository-level `README.md` provides the overall project structure and cross-module reproducibility notes. This module-level README documents the early pig fastp sensitivity branch specifically.

