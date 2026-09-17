# Cross-species transcriptomic analysis of persistent senescence-associated synovial remodeling after ACL injury

## Overview

This repository contains the analysis scripts, source tables, figures and reproducibility records for the mouse discovery, early-pig validation and chronic-pig extension cohorts. The revision package also contains the additional analyses requested during peer review.

Repository: <https://github.com/dingrui37-git/ACL_synovium_senescence_transcriptomics>

## Recommended analysis order

1. `mouse`
2. `early pig main analysis`
3. `early pig fastp sensitivity check`
4. `chronic pig`
5. `revision_submission_materials_20260911`

The first four folders contain the original analysis package. The revision folder contains the updated manuscript-facing figures, supplementary tables and data, complete result tables, audit records and reviewer-response materials.

## Revision analyses included

The revision package includes:

- the main senescence-related and mechanistic enrichment panel (Figure 6);
- the matched-gene-set directional-score null control (Supplementary Figure 5 and Supplementary Data 5);
- broader Reactome/KEGG over-representation analysis (Supplementary Figure 6 and Supplementary Data 6);
- the sex-adjusted mouse sensitivity analysis (Supplementary Table 1);
- persistent, 1-week-only and 4-week-only temporal-component analyses (Supplementary Table 2);
- the Figure 6 p53 leading-edge audit (Supplementary Data 7);
- the Fridman-signature leading-edge audit (Supplementary Data 8);
- complete tested result tables and analysis-specific QA/provenance records.

The original Hallmark inflammatory-response and EMT analyses remain associated with Figures 3–5. They were not duplicated as new supplementary figures.

## Figure 5C and retired Figure 5D

Figure 5C has been replaced with the wide-layout export. The revised script uses a 23.6-inch plotting width, enlarged value labels and manuscript-scale typography while retaining the locked 24-gene chronic DE audit and its classifications.

Figure 5D was retired during revision. Its dedicated figures, script, source tables and logs have been removed from the repository. Shared upstream files that remain necessary for Figures 5A or 5C are retained.

## Reproducibility policy

The repository excludes large raw FASTQ/BAM files, reference indexes and other non-essential intermediate files. Included scripts and tables retain the exact source-table relationships used to generate the reported figures. The revision-package index is authoritative for the current submission materials. Historical chronic-pig copy manifests are retained for provenance and are explicitly described in `chronic pig/manifest/REVISION_STATUS.md`; they are not the current Figure 5D/file inventory. Analysis-specific QA files record filtering, orthology, statistical and provenance checks.

The complete tested results are retained rather than only the significant results. Display-only operations, such as heatmap scaling or clipping, do not replace the underlying statistical tables.

## Data availability

The study reuses three public source cohorts: mouse synovium data from GEO accession `GSE271903`, early-pig synovium data from ArrayExpress accession `E-MTAB-6664`, and chronic-pig synovium data from GEO accession `GSE228848`. Processed data, source data underlying the figures, complete tested result tables, analysis scripts, audit manifests, logs and software/version records are available in this GitHub repository and in the current revision archive at [Zenodo DOI 10.5281/zenodo.22819688](https://doi.org/10.5281/zenodo.22819688); the all-versions concept record is [10.5281/zenodo.20384121](https://doi.org/10.5281/zenodo.20384121). No new primary sequencing data were generated for this study. Large raw sequencing files, BAM files, reference genome/index files and other regenerable intermediate files are not redistributed here; their public sources and regeneration procedures are documented in the analysis scripts and manifests.

## Citation and license

Please cite the associated manuscript, the versioned Zenodo release (DOI: 10.5281/zenodo.22819688) and the three public source datasets when using these data or scripts. Citation metadata are provided in `CITATION.cff`. The repository is released under the MIT License.
