# Figure 6 FDR-scope and directional-score audit

## Prespecified testing families

1. **New senescence-focused family:** the eight new entries are adjusted by BH within each comparison across the eight testable entries.
2. **Existing Hallmark contextual family:** HALLMARK_INFLAMMATORY_RESPONSE and HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION retain the `padj` values from the original complete 50-set Hallmark analysis for each comparison. They are not re-adjusted with the eight new entries.
3. **Directional-score family:** mouse score P values are retained as descriptive derivation references and receive no inferential FDR; pig score P values retain the original dataset-specific primary score-analysis BH correction, based on all 12 score-comparison rows in early pig and all six score definitions in chronic pig.

The heatmap combines these rows for visual comparison only; it does not combine their FDR families.

## Contextual Hallmark values retained from the complete 50-set analysis

| Comparison | Inflammatory-response NES | Inflammatory-response FDR | EMT NES | EMT FDR |
|---|---:|---:|---:|---:|
| Mouse 1W | 2.481 | 2.91 × 10⁻¹³ | 3.128 | 7.84 × 10⁻³² |
| Mouse 4W | 2.408 | 1.43 × 10⁻¹³ | 3.126 | 7.38 × 10⁻³⁶ |
| Pig 1W | 1.362 | 0.0305 | 1.633 | 0.00112 |
| Pig 4W | 1.012 | 0.528 | 1.513 | 0.00348 |
| Pig 52W | −1.026 | 0.415 | 2.148 | 3.36 × 10⁻⁹ |

## Directional-score values retained from the original score analyses

| Comparison | Median score difference | Primary score-analysis BH-FDR | Source family |
|---|---:|---:|---|
| Pig 1W | 0.970 | 0.00507 | Early pig: 12 score-comparison rows |
| Pig 4W | 0.691 | 0.00507 | Early pig: 12 score-comparison rows |
| Pig 52W | 0.302 | 0.00402 | Chronic pig: six score definitions |

## Directional-score definition

The mouse row is a descriptive reference in the derivation species: the 75 mouse genes corresponding to the measurable one-to-one pig orthologues and their preassigned up/down labels are applied to mouse voom logCPM values. Because membership and direction were derived from the mouse data, the paired ACLR-versus-contralateral Wilcoxon P values are not interpreted as independent validation evidence and no inferential FDR is assigned to mouse cells. The pig rows are the cross-species validation scores based on the projected 75 pig orthologues and the prespecified rank-sum comparisons; their FDR values retain the original dataset-specific primary score-analysis families shown above. The score row is therefore not an additional GSEA pathway and is not part of either pathway FDR family.
