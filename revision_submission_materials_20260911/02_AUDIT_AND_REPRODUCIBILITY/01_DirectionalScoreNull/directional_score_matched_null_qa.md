# Directional-score matched-gene-set null: frozen QA

## Frozen design

- 10,000 independent permutations for each prespecified pig comparison (early 1W, early 4W, chronic 52W).
- Each null set has 75 unique genes and retains the observed 65-up/10-down composition.
- Controls are sampled from measured pig genes with exactly one g:Profiler pig-to-human orthologue, excluding the 75 signature genes.
- Mean TMM-normalized logCPM and gene-level SD are calculated within the corresponding cohort; matching uses their joint five-by-five quantile strata.
- Sampling is without replacement within a permutation and may recur across independent permutations.
- One-sided empirical P tests the prespecified injury-concordant direction (case-minus-control directional-score difference > 0); two-sided P is a sensitivity analysis.
- Fixed seeds: pig early = `20260901`; pig chronic = `20260902`.

## Orthology eligibility audit

The frozen g:Profiler pig-to-human query covered 35,682 unique measured pig genes: 18,448 had exactly one human orthologue, 866 had multiple human orthologues, and 16,368 were unmapped. After excluding the 75 signature genes and requiring finite expression/SD with SD > 0, the eligible control universes contained 18,373 genes in the early cohort and 17,469 genes in the chronic cohort. The control eligibility criterion is the same pig-side one-to-one pig-to-human mapping step used for projected-signature construction; a reciprocal human-to-mouse filter was not imposed because the current orthology snapshot would exclude two otherwise valid projected signature orthologues.

## Reproducibility checks

- `directional_score_matched_null_all_permutations.csv`: 30,000 rows (10,000 for each of the three comparisons).
- `directional_score_matched_null_strata_audit.csv`: 33 stratum rows; within each dataset, signature counts sum to 75, up counts to 65 and down counts to 10.
- The same row-wise z-score and signed directional-score formula is used for observed and null sets.
- Complete source code, raw g:Profiler mappings and session metadata are stored in this directory.

## Frozen results

| Comparison | Observed effect | Null 95% interval | One-sided empirical P | BH-FDR |
|---|---:|---:|---:|---:|
| Early pig, 1W | 0.970 | −0.258 to 0.382 | 0.00010 | 0.00015 |
| Early pig, 4W | 0.691 | −0.357 to 0.260 | 0.00010 | 0.00015 |
| Chronic pig, 52W | 0.302 | −0.154 to 0.219 | 0.00290 | 0.00290 |
