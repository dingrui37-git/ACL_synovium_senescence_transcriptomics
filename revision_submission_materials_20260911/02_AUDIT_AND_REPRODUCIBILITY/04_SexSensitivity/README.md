# Mouse sex-adjusted sensitivity analysis

This folder contains a read-only sensitivity analysis of the frozen mouse discovery inputs. No file under `E:/R/ACLsenescence2/rebuild_submission` was modified.

## Analysis design

- Primary model reproduced for QA: `~ treatment`
- Requested sensitivity model: `~ sex + treatment`
- Repeated-measures handling: `duplicateCorrelation(block = mouse_id)` and `lmFit(block = mouse_id)`
- Time points: 1W and 4W fitted separately
- Animals per time point: 12 (six female and six male)
- Samples per time point: 24 paired ACLR/Contra samples
- Gene universe: the primary expression-filtered gene universe was retained at each time point
- Multiple testing: BH within each time point
- Strict DEG threshold: FDR < 0.05 and absolute log2 fold change > 1

## Headline results

| Metric | 1W | 4W |
|---|---:|---:|
| Pearson correlation of treatment logFC | 0.9996 | 0.9990 |
| Treatment-direction concordance | 99.87% | 99.76% |
| Primary strict DEGs retained | 3,393/3,429 (98.95%) | 1,680/1,716 (97.90%) |
| Primary FDR-supported genes retained | 10,595/10,704 (98.98%) | 9,179/9,269 (99.03%) |

Of the 1,416 primary direction-consistent persistent genes, 1,384 (97.74%) were retained after sex adjustment. The sex-adjusted analysis identified 1,413 persistent genes in total, comprising 957 upregulated and 456 downregulated genes.

## Interpretation boundary

The analysis demonstrates robustness to adjustment for the sex main effect. It does not test a `sex × treatment` interaction and therefore does not establish that the injury response is identical between females and males.

## Main files

- `reviewer_minor1_sex_covariate_response_package.md`: reviewer response and manuscript text.
- `run_sex_adjusted_sensitivity.R`: complete reproducible analysis.
- `sex_adjusted_robustness_summary.csv`: main numerical summary.
- `persistent_gene_robustness_summary.csv`: persistent-gene retention summary.
- `sex_adjusted_QA_report.txt`: QA status.
