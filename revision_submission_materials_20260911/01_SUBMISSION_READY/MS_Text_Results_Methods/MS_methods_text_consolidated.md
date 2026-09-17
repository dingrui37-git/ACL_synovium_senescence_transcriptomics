# Methods text for the added analyses

## Senescence-focused enrichment panel

Gene-set enrichment was evaluated for the prespecified panel comprising REACTOME_CELLULAR_SENESCENCE, REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP, FRIDMAN_SENESCENCE_UP, FRIDMAN_SENESCENCE_DN, SenMayo, HALLMARK_P53_PATHWAY, HALLMARK_E2F_TARGETS and HALLMARK_G2M_CHECKPOINT. Human-derived gene sets were mapped to the corresponding species before testing. HALLMARK_INFLAMMATORY_RESPONSE and HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION were retained as contextual readouts from the analyses presented in Figures 3–5.

## Directional-score matched-gene-set null model

Specificity of the projected directional score was evaluated using 10,000 matched-gene-set permutations. Candidate genes were restricted to measurable genes eligible for the same one-to-one orthologue framework, excluding the observed signature genes. Each null set matched the observed size and up/down composition, with sampling performed within the observed joint mean-expression and variability strata.

## Reactome/KEGG over-representation analysis

Persistent-up and persistent-down genes were analyzed separately using one-sided hypergeometric over-representation tests against a common eligible, pathway-mappable background universe. Reactome and KEGG pathways were tested together within each direction and Benjamini–Hochberg adjusted.

## Sex-adjusted sensitivity analysis

Mouse differential expression was refit with the model `~ sex + treatment` while retaining the same-animal blocking structure used in the primary analysis. Concordance was assessed using treatment-associated logFC correlations, genome-wide direction agreement, strict-DEG retention and retention of the persistent direction-consistent program.
