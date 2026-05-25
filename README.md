# Cross-species transcriptomic analysis identifies a persistent senescence-associated synovial remodeling program after ACL injury

This repository contains the manuscript reproducibility package for the cross-species synovial transcriptomic analysis of anterior cruciate ligament (ACL) injury-associated senescence and remodeling.

The package is intended to help reviewers and readers inspect the figure source data, processed result tables, scripts, logs, and version records used in the manuscript. The original public sequencing datasets are not redistributed here; instead, public accession numbers and provenance records are provided so that the analyses can be rerun from the original sources where applicable.

## Manuscript

**Title:** Cross-species transcriptomic analysis identifies a persistent senescence-associated synovial remodeling program after anterior cruciate ligament injury

## Public datasets and reference resources

- **Mouse discovery cohort:** GSE271903
- **Early pig validation cohort:** E-MTAB-6664
- **Chronic pig extension cohort:** GSE228848
- **Senescence reference:** CellAge `cellage3.tsv`

Raw-data provenance and accession information are provided in `public_dataset_accessions.csv` and/or `raw_data_manifest/README_raw_data.md`, if present.

## Repository organization

The repository is organized as a manuscript-level reproducibility archive. Folder names may be adapted to the final GitHub upload structure, but the intended contents are:

```text
ACL_synovium_senescence_transcriptomics/
├─ README.md
├─ LICENSE
├─ CITATION.cff
├─ public_dataset_accessions.csv
├─ raw_data_manifest/
├─ 01_mouse_discovery/
├─ 02_pig_early/
├─ 03_pig_early_fastp_sensitivity/
├─ 04_pig_chronic/
├─ scripts/
├─ source_data/
├─ tables/
├─ results/
├─ figures/
├─ logs_versions/
└─ external_data/
```

### Main folders

- `01_mouse_discovery/`: mouse discovery analysis, including paired limma-voom differential expression, PCA/MDS-related source data where applicable, persistent gene definition, CellAge analysis, ortholog mapping, Hallmark GSEA, and Figure 1–3 source files.
- `02_pig_early/`: early pig validation analysis, including sample manifests, QC summaries, featureCounts gene-level count matrix, edgeR QLF differential expression results, current signature/ortholog validation tables, TMM-logCPM signature score source data, full Hallmark GSEA outputs, Figure 4 source data/final panels, PCA/MDS supplements, singscore sensitivity outputs, and logs/version records.
- `03_pig_early_fastp_sensitivity/`: fastp-trimmed sensitivity branch for the early pig analysis, including H2w scripts, fastp QC summaries, trimmed count matrix and differential expression outputs, DE/signature/GSEA/core-gene stability summaries, supplementary figures, source data, statistics, logs, and version records.
- `04_pig_chronic/`: chronic pig extension analysis for Figure 5, including Salmon/tximport audit files, synovium-only matrix lock supplements, chronic PCA/MDS supplements, main 24-sample comparison matrices, chronic signature scoring, singscore sensitivity, edgeR QLF differential expression, full Hallmark GSEA, chronic core-gene audit outputs, early-versus-chronic summaries, source data, and figure panels.
- `scripts/`: final manuscript R scripts. If scripts are also stored inside module-specific folders, the module-level folder should be treated as the authoritative location for that branch.
- `source_data/`: source data used to generate manuscript figure panels and statistical summaries.
- `tables/`: processed result tables, including differential expression, persistent genes, CellAge overlaps, ortholog mappings, signature validation, core-gene audits, and GSEA outputs.
- `results/`: additional processed analysis outputs not classified as figure source data or final tables.
- `figures/`: final exported figures and supplementary figure files.
- `logs_versions/`: run logs, sessionInfo files, software/package version records, fgsea settings, gprofiler/gorth audit files, and copy/audit manifests.
- `external_data/`: small external reference files copied locally where permitted, such as cleaned or audited CellAge reference files. Large public raw data and regenerable reference files are not redistributed.

## Module-specific documentation

Detailed branch-level README files are retained inside the corresponding module folders. These files document exact package-generation logic, branch-specific warnings, audit manifests, and final-vs-old script decisions.

Recommended module README names:

- `01_mouse_discovery/README_mouse_discovery.md`
- `02_pig_early/README_pig_early.md`
- `03_pig_early_fastp_sensitivity/README_fastp_sensitivity.md`
- `04_pig_chronic/README_pig_chronic.md`

The root `README.md` provides the project-level overview. The module README files provide detailed branch-level notes and should be checked before rerunning or auditing each analysis branch.

## Important methodological anchors

### Mouse discovery

- Differential expression was performed using a paired limma-voom workflow with `duplicateCorrelation`.
- PCA used the second voom-transformed logCPM matrix from the locked differential expression workflow.
- Strict DEGs were defined as genes with FDR < 0.05 and |logFC| > 1.
- Persistent genes were defined as strict DEGs at both early time points with direction-consistent logFC values.
- CellAge and ortholog analyses used cleaned CellAge gene symbols and strict one-to-one ortholog mapping.
- Mouse Hallmark GSEA used full gene-level ranked lists, with the rank statistic defined as `sign(logFC) * -log10(P value)`.

### Early pig validation

- Early pig validation used raw FASTQ-derived gene-level counts generated through the final alignment/featureCounts workflow.
- Differential expression was performed using an edgeR quasi-likelihood framework.
- Cross-species validation used the current CellAge-overlapping persistent signature after strict one-to-one ortholog projection.
- Signature scoring used TMM-normalized logCPM-based source data and direction-aware scoring logic.
- Figure 4 source data, final panels, PCA/MDS supplements, singscore sensitivity outputs, and logs/version records are included where available.
- The final GSEA implementation should use the locked final script documented in the module-level README, rather than older superseded Step20 versions.

### Early pig fastp sensitivity branch

- The fastp sensitivity branch compares the primary untrimmed early pig analysis with a trimmed-analysis branch.
- Included outputs cover fastp QC summaries, trimmed alignment/featureCounts summaries, trimmed count matrix, trimmed DE outputs, DE stability, signature-score stability, GSEA stability, core-gene audit source data, supplementary figures, statistics, logs, and version records.
- Raw FASTQ, fastp-trimmed FASTQ, BAM files, BAM summary files, and full per-sample fastp HTML/JSON reports are intentionally excluded from the GitHub/Zenodo package.

### Chronic pig extension

- The chronic pig extension used public Salmon `quant.sf` processed quantification files.
- Gene-level matrices were generated using tximport, with quantification audit, matrix-lock, and main-comparison subset files retained.
- The main chronic comparison focuses on Control_52W versus ACLT_alone_52W synovium samples.
- Chronic analyses include TMM-aligned signature scoring, singscore sensitivity, edgeR QLF differential expression, full Hallmark GSEA, chronic core-gene audit outputs, and early-versus-chronic summary panels.
- Figure 5 panel policy follows the locked chronic branch:
  - Figure 5A: TMM-aligned chronic signature score.
  - Figure 5B: target Hallmark GSEA running curves.
  - Figure 5C: Step4-locked chronic core-gene differential expression audit.
  - Figure 5D: early-versus-chronic summary panel.
- The module-level chronic README should be consulted for exact B/C/P/D/N/O/E/F/A/G/M/H/J/K branch order and consistency warnings.

## Files intentionally excluded

To keep the GitHub/Zenodo archive lightweight and submission-ready, the following large or regenerable files are intentionally excluded unless required by a journal editor or reviewer:

- raw FASTQ files
- fastp-trimmed FASTQ files
- BAM files
- BAM summary files
- reference FASTA files
- reference GTF files
- Subread index files
- full per-sample fastp HTML/JSON reports
- large `.RData` or workspace files
- other public or regenerable raw-data files

The original datasets remain available through their public repositories. Processed source data, scripts, logs, version records, and audit manifests are provided for manuscript-level reproducibility.

## Recommended audit before GitHub/Zenodo upload

Before creating a GitHub release and archiving the repository in Zenodo, check the available manifest/audit files in each module. Typical files include:

- `MANIFEST_all_copied_files.csv`
- `MANIFEST_missing_expected_files.csv`
- `copied_file_manifest.csv`
- `missing_required_or_optional_files.csv`
- `script_qc_report.csv`
- `copy_manifest_exact.csv`
- `missing_required_files.csv`

Recommended checks:

1. Confirm that all required files are present.
2. Confirm that all required rows in copy manifests were successfully copied.
3. Confirm that missing-required-file tables have zero rows.
4. Confirm that final scripts, source data, figure panels, logs, and version records are included.
5. Confirm that old or superseded scripts are not presented as final analysis scripts.
6. Confirm that public raw-data accessions are listed and that large raw/regenerable files are not redistributed.

## Reproducibility notes

- Analyses were conducted using public mouse and pig synovial RNA-seq datasets.
- Mouse discovery used paired modeling where applicable.
- Pig early validation used raw FASTQ-derived featureCounts gene-level counts.
- Pig chronic extension used public Salmon quantification files summarized by tximport.
- GSEA analyses used a consistent rank statistic across datasets: `sign(logFC) * -log10(P value)`.
- Multiple-testing correction used FDR/adjusted P values where applicable.
- Software versions, session information, fgsea settings, and gprofiler/gorth audit files are retained in the logs/version folders when available.

## Data and code availability statement

The R scripts, processed source data, analysis outputs, figure source files, logs, and version records generated for this study are provided in this repository and will be archived in Zenodo. The original public RNA-seq datasets are available from GSE271903, E-MTAB-6664, and GSE228848.

After Zenodo archiving, replace the placeholder below with the final DOI:

```text
Zenodo DOI: https://doi.org/10.5281/zenodo.XXXXXXX
```

## Citation

If using this repository, please cite the associated manuscript and the Zenodo DOI for the archived release.

A `CITATION.cff` file should be added or updated after the final author list, manuscript title, repository URL, and Zenodo DOI are finalized.

## License

A license file should be added before public release. For manuscript reproducibility packages, common options include:

- MIT License for code
- CC BY 4.0 for documentation and processed non-sensitive source data

Choose the license according to journal policy, institutional requirements, and co-author agreement.
