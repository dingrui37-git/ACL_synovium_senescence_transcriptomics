# Reviewer Minor Comment 1 — sex-adjusted sensitivity analysis

Package readiness: `draft_with_placeholders`

The statistical analysis and response wording are complete. Manuscript line numbers and the final Supplementary Table/Source Data designation remain placeholders because the revised manuscript file and supplement numbering were not supplied.

## Preserved reviewer comment

> Please report whether sex was modeled or tested as a covariate in the mouse cohort, given paired same-animal design includes both sexes.

## Action classification

- Comment type: methodological/statistical clarification
- Severity: minor
- Actions: `CLARIFY_EXISTING` + `ACCEPT_ANALYSIS` + `ACCEPT_TEXT`
- Primary model status: unchanged
- Added analysis: sex-adjusted sensitivity analysis
- Interaction analysis: not performed; this analysis evaluates adjustment for the sex main effect, not sex-dependent injury responses

## Draft response to the reviewer

We thank the reviewer for raising this point. Sex was recorded for every mouse but was not included as an explicit covariate in the primary differential-expression model. The 1-week and 4-week cohorts each comprised 12 mice (six females and six males), and every animal contributed matched ACLR and contralateral synovial samples. The primary analysis was performed separately at each time point using limma–voom with treatment as the fixed effect and mouse identity as the blocking factor through `duplicateCorrelation`. Consequently, the ACLR-versus-contralateral contrast was estimated within the same animal, so sex, which is constant within each pair, could not confound the treatment contrast.

To confirm that between-animal expression differences associated with sex did not materially affect the results, we repeated the analysis at each time point using a sex-adjusted design (`~ sex + treatment`) while retaining mouse identity as the blocking factor and the same expression-filtered gene universe as in the primary analysis. Treatment-effect estimates were nearly identical between the primary and sex-adjusted models (Pearson correlation of log2 fold changes: r = 0.9996 at 1 week and r = 0.9990 at 4 weeks), with directional concordance of 99.87% and 99.76%, respectively. Of the genes meeting the primary strict differential-expression criteria (BH-FDR < 0.05 and |log2 fold change| > 1), 3,393 of 3,429 (98.95%) at 1 week and 1,680 of 1,716 (97.90%) at 4 weeks retained strict differential-expression status after adjustment. Likewise, 1,384 of the 1,416 direction-consistent persistent genes (97.74%) were retained, including 953 of 970 persistently upregulated genes (98.25%) and 431 of 446 persistently downregulated genes (96.64%). These results indicate that the principal differential-expression and persistent-gene findings are robust to adjustment for sex. We have clarified the paired design and added the sex-adjusted sensitivity analysis in the Methods and Results [Methods, lines XX–XX; Results, lines XX–XX; Supplementary Table/Source Data XX].

This sensitivity analysis addresses adjustment for the sex main effect; it was not designed to establish the absence of a sex-by-injury interaction.

## Proposed Methods insertion

### Sex-adjusted sensitivity analysis

The primary mouse differential-expression analysis was conducted separately at 1 and 4 weeks using limma–voom, with treatment (ACLR versus contralateral synovium) included as a fixed effect and mouse identity specified as the blocking factor using `duplicateCorrelation`. Each time-point cohort contained 12 mice (six female and six male), and each mouse contributed one ACLR and one contralateral sample. As a sensitivity analysis, the models were refitted using a design that included sex and treatment (`~ sex + treatment`) while retaining mouse identity as the blocking factor. The expression-filtered gene universe from the corresponding primary model was held fixed so that the comparison isolated the effect of adding sex to the design matrix. The treatment coefficient was evaluated using empirical Bayes moderation, and P values were adjusted within each time point using the Benjamini–Hochberg method. This analysis assessed robustness to adjustment for the sex main effect and did not test a sex-by-treatment interaction.

## Proposed Results insertion

Adjustment for sex had minimal influence on the estimated ACLR-versus-contralateral transcriptional response. Genome-wide log2 fold changes from the primary and sex-adjusted models were highly correlated at both 1 week (r = 0.9996) and 4 weeks (r = 0.9990), with directional concordance exceeding 99.7%. Among strict differentially expressed genes, 98.95% at 1 week and 97.90% at 4 weeks retained strict significance after adjustment. The persistent program was similarly stable: 1,384 of 1,416 direction-consistent persistent genes (97.74%) were retained in the sex-adjusted analysis. Thus, the main mouse differential-expression and persistent-gene findings were not materially altered by adjustment for sex [Supplementary Table/Source Data XX].

## 中文核对

- 原主分析没有显式加入sex，但同一只小鼠同时提供ACLR和Contra样本，且每个时间点均为6雌、6雄，因此sex不会与同鼠内的treatment主效应混杂。
- 新增的`~ sex + treatment`分析是稳健性分析，不替换原主模型，也不需要据此重新定义正式persistent signature。
- 结果支持“对sex主效应调整稳健”：1W和4W的logFC相关性分别为0.9996和0.9990，原persistent genes保留97.74%。
- 不能将结果写成“没有性别差异”，因为没有检验sex-by-treatment interaction；每性别每时间点仅有6只小鼠，也不适合做强结论。
- 将本段并入稿件时，需要补充真实Methods/Results行号，并决定把详细结果放入哪一个Supplementary Table或Source Data文件。

## Manuscript change checklist

- [ ] Methods：报告每时间点12只小鼠、6 female/6 male及同鼠ACLR–Contra配对结构。
- [ ] Methods：说明原模型为`~ treatment`并以`mouse_id`作为duplicateCorrelation/lmFit block。
- [ ] Methods：加入sex-adjusted sensitivity design、固定基因检验集合及BH校正范围。
- [ ] Results：加入logFC相关性、严格DEG保留率和persistent-gene保留率。
- [ ] Supplementary material：指定稳健性汇总表和逐基因比较表的正式编号。
- [ ] Response letter：用真实稿件行号替换`XX–XX`。

## Evidence files

- `sex_adjusted_robustness_summary.csv`: time-point-specific effect concordance and DEG retention.
- `persistent_gene_robustness_summary.csv`: persistent-program retention.
- `primary_vs_sex_adjusted_gene_comparison_1W.csv` and `primary_vs_sex_adjusted_gene_comparison_4W.csv`: complete gene-level audit.
- `sex_adjusted_DE_1W.csv` and `sex_adjusted_DE_4W.csv`: complete sex-adjusted differential-expression results.
- `primary_model_reproduction_QA.csv`: numerical reproduction of the frozen primary model.
- `sex_adjusted_QA_report.txt`: analysis integrity checks.
