# Reviewer comment 3 / Reviewer 2, comment 7

Package readiness: **draft_with_placeholders**. The analysis is complete and statistically verified. Final manuscript line numbers and Supplementary Table numbering remain to be assigned after the manuscript file is supplied.

## Response strategy summary

- Decision type: unclear from the supplied excerpt.
- Comment classification: major methodological and biological-interpretation request.
- Action: `ACCEPT_ANALYSIS` plus targeted Methods, Results and Discussion revisions.
- Overall posture: retain whole-ranked-transcriptome Hallmark GSEA for continuity, and add Reactome/KEGG over-representation analysis (ORA) of the unfiltered direction-consistent persistent mouse gene set.
- Important boundary: report pathway representation and remodeling context, not causal molecular drivers.
- Circularity control: do not use the 78 CellAge-overlap genes as the ORA input because they were already selected using a senescence database. The ORA input is the complete 1,416-gene persistent set identified before CellAge filtering.

## Comment-response tracker

| ID | Reviewer concern | Type | Severity | Action completed | Remaining input |
|---|---|---|---|---|---|
| R2.7 | Add complementary KEGG/Reactome pathway analyses to clarify the molecular processes represented by direction-consistent persistent genes, especially in Sections 3.4–3.5. | Methodological; evidence/interpretation | Major | Reactome and KEGG ORA completed for 970 upregulated and 446 downregulated persistent genes, with a 1,416-gene combined-list secondary analysis, database-aware backgrounds, complete BH correction and full source tables. | Final manuscript line numbers and Supplementary Table identifier. |

## Draft point-by-point response letter

> **Reviewer 2, Comment 7:** “Complementary pathway analyses using databases such as KEGG or Reactome.”

**Response:** We thank the reviewer for this suggestion. To provide broader biological context for the direction-consistent persistent genes, we retained the whole-ranked-transcriptome Hallmark GSEA and added an over-representation analysis using mouse Reactome and KEGG pathways. The ORA used the complete set of 1,416 genes that met the persistent-gene definition—FDR < 0.05 and |logFC| > 1 at both 1 and 4 weeks, with the same direction at both time points—before CellAge-based filtering. We analyzed the 970 persistently upregulated and 446 persistently downregulated genes separately to preserve biological direction; the combined 1,416-gene analysis was retained as secondary context.

The common eligibility universe consisted of 14,184 unique mouse Entrez genes that passed the relevant expression filtering, were formally tested at both time points, and had valid identifiers. For each database, we used the intersection of this common universe with the genes represented in that database (7,327 Reactome genes and 6,421 KEGG genes); the same database-specific annotated universe was used for the persistent-up and persistent-down analyses, with no direction-specific background. After restricting to pathways containing 10–500 genes in the corresponding annotated universe, 937 Reactome and 350 KEGG pathways were retained for testing. Enrichment was assessed with the one-sided hypergeometric upper-tail test implemented by `stats::phyper` in R 4.5.2, using Entrez identifiers and deduplicated gene–pathway mappings. Reactome definitions were obtained from mouse-native MSigDB M2:CP:REACTOME (release 2026.1.Mm), and KEGG definitions were retrieved for Mus musculus through KEGGREST (snapshot dated 31 August 2026). For each direction, Benjamini–Hochberg correction was applied jointly across all 1,287 tested Reactome and KEGG pathways. This yielded 69 FDR-supported gene-set records for the persistently upregulated genes and 32 for the persistently downregulated genes. Because many Reactome and KEGG pathways are nested or share genes, these 101 records are not interpreted as 101 independent biological mechanisms.

The upregulated persistent genes were dominated by extracellular-matrix synthesis and turnover. Reactome extracellular matrix organization included 89 of 225 pathway genes (enrichment ratio = 5.76, FDR = 1.13 × 10^-43), collagen formation included 36 of 71 genes (enrichment ratio = 7.39, FDR = 5.93 × 10^-21), and extracellular-matrix degradation included 38 of 81 genes (enrichment ratio = 6.83, FDR = 8.07 × 10^-21). Concordant KEGG results included ECM–receptor interaction (26/80 genes; enrichment ratio = 4.41; FDR = 3.25 × 10^-9), focal adhesion (32/191 genes; enrichment ratio = 2.27; FDR = 3.43 × 10^-4), and cytokine–cytokine receptor interaction (32/173 genes; enrichment ratio = 2.51; FDR = 4.71 × 10^-5). Reactome neutrophil degranulation and hemostasis were also enriched (FDR = 0.026 and 0.030, respectively), providing additional inflammatory and vascular-remodeling context.

In contrast, the downregulated persistent genes were enriched for metabolic and cytoprotective processes. KEGG PPAR signaling included 19 of 62 genes (enrichment ratio = 7.72; FDR = 1.47 × 10^-9), while AMPK signaling included 18 of 115 genes (enrichment ratio = 3.94; FDR = 1.28 × 10^-4). Reactome lipid metabolism included 41 of 496 genes (enrichment ratio = 2.60; FDR = 6.21 × 10^-6), and cytoprotection by HMOX1 included 9 of 43 genes (enrichment ratio = 6.58; FDR = 1.01 × 10^-3).

These results refine the broad Hallmark interpretation by showing that the persistent response is concentrated in coordinated matrix production and degradation, adhesion and inflammatory signaling, together with reduced lipid/metabolic and cytoprotective programs. Descriptively, most persistent effects were larger at 1 week than at 4 weeks: only 205/970 upregulated genes and 43/446 downregulated genes showed a larger absolute effect at 4 weeks. Thus, the data support a sustained but quantitatively attenuated remodeling program rather than a uniformly increasing response. We have added these analyses to **[Methods—new Reactome/KEGG ORA subsection]**, the principal findings to **[Results, Sections 3.4–3.5]**, and the complete all-tested pathway results, together with the FDR-supported subset and overlap genes, to **[Supplementary Table X/source-data file]**. The supplementary dot plot shows up to eight lowest-FDR supported pathways per database-by-direction stratum rather than padding strata with nonsignificant pathways. We also revised the Discussion to clarify that ORA identifies molecular processes represented among the persistent genes but does not establish that these pathways causally drive remodeling.

## Proposed Methods insertion

### Reactome and KEGG over-representation analysis of persistent genes

Direction-consistent persistent genes were defined as genes meeting FDR < 0.05 and |logFC| > 1 in the ACLR-versus-contralateral comparison at both 1 and 4 weeks, with concordant logFC signs at the two time points. The complete persistent set contained 1,416 unique mouse Entrez genes, including 970 persistently upregulated and 446 persistently downregulated genes. To preserve direction-specific interpretation, ORA was performed separately for the upregulated and downregulated sets; an analysis of all 1,416 genes was retained as secondary context. The 78-gene CellAge-overlap subset was not used as the ORA input because its prior senescence-based selection would introduce circular pathway interpretation.

The common eligibility universe comprised 14,184 unique Entrez genes that passed the relevant expression filtering, were formally tested at both time points, and had valid identifiers. For each database, the testing universe was the intersection of this common universe with database-annotated genes, and this same database-specific universe was used for both directions. Mouse-native Reactome gene sets were obtained from MSigDB 2026.1.Mm (M2:CP:REACTOME), and mouse KEGG pathway definitions were retrieved through KEGGREST from the KEGG pathway release dated 31 August 2026. After restricting to pathways containing 10–500 genes in the corresponding annotated universe, 937 Reactome and 350 KEGG pathways were retained for testing. Over-representation was assessed with a one-sided hypergeometric upper-tail test using `stats::phyper`; Entrez identifiers were used and duplicate gene–pathway mappings were removed. For each direction, P values were adjusted jointly across all 1,287 tested Reactome and KEGG pathways using the Benjamini–Hochberg method. FDR < 0.05 defined statistical support. Database-specific FDR values and a combined-direction ORA were retained as secondary analyses. Effect-size summaries across overlapping genes were descriptive and were not treated as formal time-by-injury interaction tests.

### Final Methods wording after reviewer audit

For the final manuscript, define the common ORA eligibility universe as the 14,184 unique mouse Entrez genes that passed the relevant expression filtering, were formally tested at both 1 and 4 weeks, and had valid identifiers. For each database, use the intersection of this common universe with database-annotated genes as the same background for persistent-up and persistent-down lists. After restricting pathway membership to 10–500 genes in the corresponding annotated universe, retain 937 Reactome and 350 KEGG pathways. Apply the one-sided hypergeometric upper-tail test with `stats::phyper`, deduplicate gene–pathway mappings, and adjust P values by BH jointly across all 1,287 Reactome and KEGG pathways within each direction. This paragraph supersedes the earlier concise draft wording in this response package.

## Proposed Results insertion for Sections 3.4–3.5

To resolve the molecular processes represented by the direction-consistent response, we performed Reactome and KEGG ORA on the 970 persistently upregulated and 446 persistently downregulated genes. After joint BH correction across 937 Reactome and 350 KEGG pathways within each direction, 69 upregulated and 32 downregulated gene-set records met FDR < 0.05. The number of enriched records should not be interpreted as the number of independent mechanisms because the databases contain nested and overlapping pathways.

Persistently upregulated genes showed a dominant extracellular-matrix remodeling program encompassing both matrix production and turnover. Reactome extracellular matrix organization (89/225 genes; enrichment ratio = 5.76; FDR = 1.13 × 10^-43), collagen formation (36/71; enrichment ratio = 7.39; FDR = 5.93 × 10^-21), and extracellular-matrix degradation (38/81; enrichment ratio = 6.83; FDR = 8.07 × 10^-21) were among the strongest results. Concordant KEGG annotations grouped overlapping genes into ECM–receptor interaction (26/80; enrichment ratio = 4.41; FDR = 3.25 × 10^-9), focal adhesion (32/191; enrichment ratio = 2.27; FDR = 3.43 × 10^-4), and cytokine–cytokine receptor interaction (32/173; enrichment ratio = 2.51; FDR = 4.71 × 10^-5). Reactome neutrophil degranulation and hemostasis provided additional inflammatory and vascular-remodeling context.

Persistently downregulated genes were instead concentrated in metabolic and protective programs, including KEGG PPAR signaling (19/62; enrichment ratio = 7.72; FDR = 1.47 × 10^-9), AMPK signaling (18/115; enrichment ratio = 3.94; FDR = 1.28 × 10^-4), Reactome lipid metabolism (41/496; enrichment ratio = 2.60; FDR = 6.21 × 10^-6), and cytoprotection by HMOX1 (9/43; enrichment ratio = 6.58; FDR = 1.01 × 10^-3). Although all 1,416 genes retained direction and strict differential-expression support at both time points, only 21.1% of the upregulated genes and 9.6% of the downregulated genes had a larger absolute effect at 4 weeks than at 1 week. Together, these findings describe a persistent but generally attenuated transcriptional program coupling matrix turnover and inflammatory/adhesive remodeling with reduced metabolic-homeostatic expression.

## Proposed Discussion insertion

The Reactome/KEGG analysis adds process-level resolution to the Hallmark results. The simultaneous enrichment of collagen formation, matrix degradation, proteoglycan processing, integrin interactions and focal adhesion indicates active matrix turnover rather than simply increased matrix deposition. The concurrent representation of cytokine-receptor signaling, neutrophil degranulation and hemostasis links this structural response to inflammatory and vascular components of synovial remodeling. Conversely, the persistent reduction of PPAR–AMPK, lipid-handling and HMOX1-associated programs is consistent with reduced metabolic and cytoprotective transcriptional programs. These annotations describe coordinated processes represented among the persistent genes; they do not demonstrate pathway activity at the protein level or establish that any individual pathway causally drives remodeling. Disease-labelled KEGG terms should likewise be interpreted through their overlapping genes rather than as evidence for the named disease.

## Manuscript change checklist

- [ ] Add the ORA Methods paragraph and confirm the exact section heading.
- [ ] Add the Results text to Sections 3.4–3.5, adjusting transitions to the existing Hallmark GSEA paragraph.
- [ ] Add the Discussion boundary on pathway representation versus causal pathway activity.
- [ ] Assign a final Supplementary Table number to the complete ORA table and overlap-gene table.
- [ ] Add database citations in the journal's required reference format.
- [ ] Insert final page/line numbers into the response letter after typesetting.
- [ ] Do not describe the 101 significant gene-set records as 101 independent pathways or mechanisms.
- [ ] Do not interpret KEGG disease-labelled terms literally when enrichment is driven by shared ECM, immune or metabolic genes.

## Missing information / risk flags

- `AUTHOR_INPUT_NEEDED`: final manuscript file, section transitions, line numbers and Supplementary Table numbering.
- The ORA is based on the mouse discovery persistent set. It does not by itself establish pathway-level replication in pig; cross-species claims should remain supported by the separate whole-transcriptome GSEA analyses.
- The timepoint summaries are descriptive because they compare pathway-overlap gene effect sizes rather than fitting a formal time-by-injury interaction.

## 中文核对

- 主输入是未经 CellAge 筛选的1,416个小鼠方向一致 persistent genes，不是78个 CellAge重叠基因，也不是猪的75基因评分集合。
- 主分析分别检验970个持续上调和446个持续下调基因；合并1,416基因的结果只作补充背景。
- 主FDR是在每个方向内，把937个Reactome和350个KEGG通路共1,287项联合做BH，不按结果挑选通路。
- 结果支持“持续但总体减弱的基质周转—炎症/黏附重塑，并伴随代谢和细胞保护程序降低”，不支持直接写成这些通路因果性地驱动重塑。
