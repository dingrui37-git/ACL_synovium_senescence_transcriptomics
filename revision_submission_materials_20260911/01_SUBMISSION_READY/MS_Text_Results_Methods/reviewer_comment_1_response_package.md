# Reviewer comment 1: matched-gene-set negative control

## Response strategy summary

- Decision type: unclear (not supplied)
- Task mode: draft plus completed analysis
- Package readiness: `draft_with_placeholders`
- Classification: statistical/methodological; major; `ACCEPT_ANALYSIS`
- Overall posture: accept the concern and add a signature-specific matched random-gene-set null analysis
- Remaining placeholders: manuscript section, page/line numbers, and final supplementary item number

## Comment-response tracker

| ID | Reviewer concern | Type | Severity | Proposed action | Readiness | Missing author input |
|---|---|---|---|---|---|---|
| R1.1 | The mouse-defined directional score may reward direction consistency by construction and requires a negative control demonstrating signature specificity in pig | Statistical / methodological | Major | `ACCEPT_ANALYSIS`; add 10,000 expression- and variability-matched random gene sets preserving the 75-gene size and 65/10 direction composition | `draft_with_placeholders` | Confirm manuscript and supplementary locations |

## Draft point-by-point response

**Reviewer comment R1.1**

> Directional score requires a negative control. The primary readout is a mouse-defined directional z-score projected into pig. Because the score is constructed to reward direction consistency, it is vulnerable to self-confirmation. The singscore sensitivity analysis helps, but a permutation/null control (e.g., randomly drawn gene sets of matched size, or non-orthologous control sets) is needed to show the activation is specific to this signature rather than a generic injury response. Please add this.

**Response**

We thank the reviewer for identifying this important specificity control. We agree that the existing singscore analysis tests robustness to the scoring method but does not establish that the projected 75-gene pig orthologue set performs beyond comparably expressed background gene sets.

To address this concern, we added a matched random-gene-set analysis in both pig datasets. For each dataset, we generated 10,000 control gene sets from measurable pig genes satisfying the same pig-side one-to-one pig-to-human orthology eligibility criterion used to construct the projected signature (exactly one g:Profiler pig-to-human orthologue), after excluding the 75 signature genes. Each control set contained 75 unique genes. Matching was performed jointly within the observed mean-expression × across-sample-SD strata: the number of genes assigned to the up and down components was preserved in every stratum, yielding the observed 65:10 directional composition without assigning arbitrary directions after an unmatched draw. Genes were sampled without replacement within each permutation, while resampling across independent permutations was allowed. Each random set was scored using exactly the same gene-wise z transformation and directional-score formula as the observed signature. We then compared the observed case-minus-control median score difference with the corresponding null distribution. The empirical one-sided P value was calculated as `(1 + number of null effects greater than or equal to the observed effect)/(10,000 + 1)`; two-sided empirical P values were also evaluated, and multiple testing was controlled across the three prespecified pig comparisons using the Benjamini-Hochberg procedure.

The observed median directional-score differences were 0.970 at 1 week and 0.691 at 4 weeks in the early pig dataset, compared with matched-null 95% intervals of -0.258 to 0.382 and -0.357 to 0.260, respectively. The corresponding one-sided empirical P values were both 0.00010 (BH-FDR = 0.00015 for each comparison). In the chronic pig dataset, the observed 52-week difference was 0.302, compared with a matched-null 95% interval of -0.154 to 0.219 (empirical P = 0.00290; BH-FDR = 0.00290). Two-sided empirical P values were 0.00010, 0.00010 and 0.00330, respectively (BH-FDR = 0.00015, 0.00015 and 0.00330). These negative-control results indicate that the projected directional activation is not reproduced by random gene sets with the same size, directional composition, expression distribution, and variability distribution within an orthology-eligible pig background. We have added the analysis to [Methods section], the results to [Results section], and the complete null distributions and matching audit to [Supplementary Table/Figure placeholder].

## Proposed Methods text

### Matched-gene-set negative-control analysis

To test whether the cross-species directional score was specific to the mouse-defined signature rather than a generic property of the pig injury transcriptome, we generated 10,000 matched random gene sets independently in the early and chronic pig datasets. The fixed mouse-derived signature comprised 75 measurable one-to-one pig orthologues, including 65 genes assigned to the up component and 10 assigned to the down component. Candidate control genes were drawn from measurable pig genes satisfying the same pig-side one-to-one pig-to-human orthologue eligibility criterion (exactly one g:Profiler pig-to-human orthologue) in the corresponding TMM-normalized logCPM matrix, excluding the 75 signature genes. Gene-level variability was calculated as the standard deviation across all samples within the corresponding cohort. Genes were jointly stratified into five quantile bins for mean expression and five quantile bins for variability. Within each mean-expression × variability stratum, each random set sampled the required number of unique controls without replacement and retained the observed numbers assigned to the up and down components; across independent permutations, genes could recur. Random sets were processed using the same gene-wise z standardization and directional-score definition as the observed signature. For each of the prespecified comparisons (early pig ACLT at 1 week versus control, early pig ACLT at 4 weeks versus control, and chronic pig ACLT at 52 weeks versus control), the test statistic was the difference in median directional score between cases and controls. Empirical one-sided P values for activation were calculated as `(1 + sum(null statistic >= observed statistic))/(10,000 + 1)` because the injury-concordant direction was prespecified as positive. Two-sided empirical P values based on absolute statistics were examined as a sensitivity analysis. Benjamini-Hochberg correction was applied across the three prespecified comparisons.

## Proposed Results text

The mouse-derived directional signature remained specific relative to orthology-eligible, expression- and variability-matched random gene sets. In the early pig dataset, the observed median case-control score differences were 0.970 at 1 week and 0.691 at 4 weeks; both exceeded all 10,000 matched-null effects (one-sided empirical P = 0.00010 for each; BH-FDR = 0.00015). The corresponding matched-null 95% intervals were -0.258 to 0.382 and -0.357 to 0.260. In the chronic pig dataset, the observed 52-week difference was 0.302, exceeding 9,972 of 10,000 matched-null effects (one-sided empirical P = 0.00290; BH-FDR = 0.00290), with a null 95% interval of -0.154 to 0.219. Two-sided empirical tests yielded the same conclusion (BH-FDR <= 0.00330). Thus, the cross-species directional signal was not reproduced by random gene sets matched for size, directional composition, mean expression, and variability within a pig-side one-to-one orthology-eligible background.

## Manuscript change checklist

- Add the matched-gene-set negative-control procedure to Methods.
- Add the three observed effects, null intervals, empirical P values, and BH-FDR values to Results.
- Provide the complete 30,000 null statistics and the stratum-capacity audit as Source Data or a Supplementary Table.
- Add a compact null-distribution panel only if the target figure layout permits; do not assign a figure number until the current manuscript is inspected.
- State clearly that singscore tests scoring-method robustness, whereas the matched random-gene-set analysis tests signature specificity.
- Add the analysis script and session information to the reproducibility package.

## Missing information / risk flags

- Current manuscript, supplementary-material version, and line numbering have not yet been supplied, so all locations remain placeholders.
- The analysis supports specificity relative to matched random pig gene sets; it does not by itself prove cell-type specificity, causal senescence biology, or absence of all generic injury contributions.
- The final manuscript should retain calibrated wording such as “specific relative to matched random gene sets” rather than “uniquely specific to senescence.”

## 中文核对

- 已完成真实补分析，不是仅提供方案。
- 早期猪 1 周、4 周及慢性猪 52 周三个比较均通过匹配随机基因集检验。
- 最终回复前需要确认当前稿件的 Methods、Results、补充图表和行号位置。
- 该分析能回答“是否优于同规模、同表达/变异分布的随机基因集”，不能扩大解释为“已证明衰老特异性或因果机制”。
