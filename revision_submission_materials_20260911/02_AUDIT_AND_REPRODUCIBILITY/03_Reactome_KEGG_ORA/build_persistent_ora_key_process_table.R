options(stringsAsFactors = FALSE)

out_dir <- "D:/workspace/ACLsenescence2_reviewer3_persistent_ORA"
ora_file <- file.path(out_dir, "persistent_ORA_all_Reactome_KEGG_results.csv")
ora <- utils::read.csv(ora_file, stringsAsFactors = FALSE, check.names = FALSE)

key_terms <- data.frame(
  process_axis = c(
    rep("matrix synthesis, turnover and adhesion", 7),
    rep("inflammatory and vascular remodeling", 3),
    rep("metabolic and cytoprotective programs", 6)
  ),
  gene_list = c(
    rep("persistent_up", 10),
    rep("persistent_down", 6)
  ),
  database = c(
    "Reactome", "Reactome", "Reactome", "Reactome", "KEGG", "KEGG", "KEGG",
    "Reactome", "Reactome", "KEGG",
    "KEGG", "KEGG", "KEGG", "Reactome", "Reactome", "Reactome"
  ),
  pathway_name = c(
    "Extracellular matrix organization",
    "Collagen formation",
    "Degradation of the extracellular matrix",
    "Integrin cell surface interactions",
    "ECM-receptor interaction",
    "Focal adhesion",
    "PI3K-Akt signaling pathway",
    "Neutrophil degranulation",
    "Hemostasis",
    "Cytokine-cytokine receptor interaction",
    "PPAR signaling pathway",
    "AMPK signaling pathway",
    "Glycolysis / Gluconeogenesis",
    "Metabolism of lipids",
    "Metabolism of amino acids and derivatives",
    "Cytoprotection by HMOX1"
  ),
  stringsAsFactors = FALSE
)

key <- merge(key_terms, ora, by = c("gene_list", "database", "pathway_name"), all.x = TRUE, sort = FALSE)
if (nrow(key) != nrow(key_terms) || any(is.na(key$pval))) {
  stop("Failed to recover every prespecified key process row.", call. = FALSE)
}

key$selection_role <- "representative biological summary; inference remains based on the complete ORA table"
keep <- c(
  "process_axis", "gene_list", "database", "database_release", "pathway_id", "pathway", "pathway_name",
  "universe_size", "query_size_in_universe", "pathway_size_in_universe", "overlap_size",
  "expected_overlap", "enrichment_ratio", "odds_ratio_haldane", "pval",
  "FDR_within_database", "FDR_joint_Reactome_KEGG", "primary_FDR", "primary_FDR_status",
  "median_logFC_1W_overlap", "median_logFC_4W_overlap",
  "median_abs_logFC_1W_overlap", "median_abs_logFC_4W_overlap",
  "median_delta_logFC_4W_minus_1W", "proportion_overlap_stronger_abs_effect_at_4W",
  "overlap_gene_symbols", "selection_role"
)
key <- key[, keep, drop = FALSE]
key <- key[order(match(key$process_axis, unique(key_terms$process_axis)),
                 match(key$gene_list, c("persistent_up", "persistent_down")),
                 key$primary_FDR), ]

utils::write.csv(
  key,
  file.path(out_dir, "persistent_ORA_representative_processes_for_manuscript.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8",
  na = ""
)
cat("Saved representative process table with", nrow(key), "rows.\n")
