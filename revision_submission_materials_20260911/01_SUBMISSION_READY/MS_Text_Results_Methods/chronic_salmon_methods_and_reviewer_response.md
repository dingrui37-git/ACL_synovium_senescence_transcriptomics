# GSE228848 chronic-pig Salmon workflow: audit and replacement text

## Audit conclusion

The original Salmon quantification and the present tximport summarization are distinct stages and should be described separately.

| Item | Original GSE228848 processing | This reanalysis |
|---|---|---|
| Quantifier | Salmon 1.8.0 | Salmon was not rerun |
| Input | Untrimmed FASTQ reads | Deposited per-sample `quant.sf` files |
| Mapping/reference information reported by the source | Quasi-mapping against the Sscrofa11.1 porcine reference hosted by Ensembl | Not reconstructed or inferred |
| Exact transcript FASTA, Ensembl release, index identifier and complete Salmon command line | Not reported in the GEO record or original paper and not included in the supplementary archive | Therefore cannot be assigned from the `quant.sf` files alone |
| Annotation used for downstream aggregation | Not applicable | `Sus_scrofa.Sscrofa11.1.115.gtf.gz` (Ensembl release 115), used only to make `tx2gene` |
| Summarization | Not applicable | tximport 1.38.2, `type = "salmon"`, `ignoreTxVersion = TRUE`, `countsFromAbundance = "no"` |

The GEO record reports Salmon 1.8.0, quasi-mapping of untrimmed FASTQ reads, the Sscrofa11.1 assembly, and the five standard fields in the deposited files (`Name`, `Length`, `EffectiveLength`, `TPM`, and `NumReads`). It also states that raw reads are available in SRA and that processed quantification files are provided as supplementary files. Thus, the manuscript should not say that raw FASTQ files were unavailable. The accurate statement is that raw reads were not re-quantified in this analysis because the original index and complete Salmon parameters were not available for a reproducible reconstruction.

## Compatibility-audit result

The compatibility audit was performed after removing trailing transcript-version suffixes, matching the behavior of `ignoreTxVersion = TRUE`. Each of the 96 files contained the same 46,295 unique normalized transcript identifiers. Of these, 45,845 (99.02797%) were present in the Ensembl release 115 `tx2gene` table, leaving 450 (0.97203%) without a release-115 transcript-to-gene entry. Thus, the sentence proposed for the manuscript can be completed as follows:

> “After removal of transcript-version suffixes, 45,845 of 46,295 unique deposited transcript identifiers (99.03%) were successfully mapped to Ensembl release 115 gene identifiers.”

As an optional supplementary QC, these mapped transcripts represented 95.98% of estimated `NumReads` on average across the 96 files (range, 90.64–98.63%; all-file count-weighted proportion, 96.12%). This abundance-weighted value is distinct from the requested identifier-level compatibility percentage and should not replace the 99.03% primary statement.

## Recommended replacement for Methods 2.1

For the chronic pig cohort (GSE228848), we used the processed per-sample Salmon transcript-quantification files deposited in the GEO supplementary archive rather than re-quantifying the linked SRA reads. The original study reports Salmon version 1.8.0, quasi-mapping of untrimmed FASTQ reads against the Sscrofa11.1 porcine reference hosted by Ensembl. However, the GEO record and the original publication do not specify the exact transcript FASTA, Ensembl release, Salmon-index identifier or complete command line; these parameters therefore were not inferred from the deposited files. The 96 deposited files (48 synovium and 48 cartilage) were checked for readability and for the standard Salmon columns `Name`, `Length`, `EffectiveLength`, `TPM` and `NumReads`; all files passed these checks. After removal of transcript-version suffixes, 45,845 of 46,295 unique deposited transcript identifiers (99.03%) matched a transcript-linked gene identifier in the release-115 `tx2gene` table. For downstream gene-level summarization, we constructed a separate transcript-to-gene map from the `Sus scrofa` Sscrofa11.1 Ensembl release 115 GTF (`Sus_scrofa.Sscrofa11.1.115.gtf.gz`), retaining unique `transcript_id`–`gene_id` pairs. The deposited transcript-level estimates were then aggregated with tximport (version 1.38.2; `type = "salmon"`, `ignoreTxVersion = TRUE`, `countsFromAbundance = "no"`) to obtain gene-level estimated counts, TPM abundance and effective-length matrices. `ignoreTxVersion = TRUE` was used because the deposited transcript identifiers contain Ensembl version suffixes (for example, `.4`). Salmon quantification itself was not rerun, and no unreported index parameter was reconstructed. The synovium subset comprised Control_52-week, ACLT_alone_52-week, Reconstruction_52-week and Repair_52-week groups (12 samples per group); the primary chronic comparison was Control_52-week versus ACLT_alone_52-week.

## Recommended citations

Use the original GSE228848 paper for the cohort-specific workflow and cite the Salmon method paper for the quantification algorithm:

1. Donnenfield et al. *Bioengineering* (2023), doi:10.3390/bioengineering10050527.
2. Patro et al. *Nature Methods* (2017), doi:10.1038/nmeth.4197.
3. Cite the tximport software paper for transcript-to-gene summarization: Soneson C, Love MI and Robinson MD. *Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences*. F1000Research (2015), 4:1521, doi:10.12688/f1000research.7563.1 (tximport version 1.38.2 in the present environment).

## Recommended response to the reviewer

We thank the reviewer for identifying this ambiguity. We have revised Section 2.1 to distinguish the original transcript quantification from our downstream gene-level summarization. The GEO record and the source publication specify that the original study used Salmon v1.8.0 to quasi-map untrimmed FASTQ reads against the Sscrofa11.1 porcine reference hosted by Ensembl. They do not report the exact transcript FASTA, Ensembl release, Salmon-index identifier or complete command-line parameters, and these details are not included in the GEO supplementary archive. We therefore do not assign Ensembl release 115 to the original Salmon index and did not attempt to infer unreported parameters from the `quant.sf` files. Instead, we used the 96 processed per-sample quantification files deposited by the study, verified their readability and required Salmon fields, and constructed an independent `tx2gene` table from the Sscrofa11.1 Ensembl release 115 GTF solely for tximport aggregation. After transcript-version normalization, 45,845 of 46,295 unique deposited transcript identifiers (99.03%) matched the release-115 `tx2gene` table. We now report tximport version 1.38.2 and its parameters (`type = "salmon"`, `ignoreTxVersion = TRUE`, and `countsFromAbundance = "no"`), explain the version-suffix handling, correct the wording concerning raw FASTQ availability, and cite both the original GSE228848 study and the established Salmon method (Patro et al., 2017).

## Sources checked

- [GEO sample record GSM7140500 / GSE228848](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7140500)
- [Original GSE228848 study (PubMed)](https://pubmed.ncbi.nlm.nih.gov/37237597/)
- [Salmon method paper (PubMed)](https://pubmed.ncbi.nlm.nih.gov/28263959/)
- [tximport Bioconductor documentation](https://bioconductor.posit.co/packages/3.23/bioc/html/tximport.html)

## Audit files

- `audit_transcript_id_compatibility.R` reproduces the identifier-level audit and the optional `NumReads`-weighted QC.
- `transcript_mapping_compatibility_summary.csv` reports 45,845/46,295 mapped identifiers (99.02797%), 450 unmapped identifiers, and the abundance-weighted QC values.
