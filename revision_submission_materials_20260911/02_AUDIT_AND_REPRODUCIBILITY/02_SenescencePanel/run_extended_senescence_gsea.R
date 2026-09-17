options(stringsAsFactors = FALSE)

## Reviewer 2: prespecified senescence/SASP MSigDB panel enrichment.
## Primary multiplicity family: all testable panel terms within each comparison.
## Reuses frozen ranked vectors and writes only to D:/workspace.

out_dir <- "D:/workspace/ACLsenescence2_reviewer2_senescence_pathways"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("fgsea", "msigdbr")
for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p, call. = FALSE)
}

min_size <- 10L
max_size <- 500L
fgsea_eps <- 0
nproc <- 1L
panel_name_regex <- "SENESCEN|SASP"
panel_fdr_method <- "BH"

fdr_status <- function(fdr, pval) {
  ifelse(is.na(fdr), "not_available",
         ifelse(fdr < 0.05, "FDR_lt_0.05",
                ifelse(!is.na(pval) & pval < 0.05, "nominal_only", "not_significant")))
}

rank_specs <- data.frame(
  dataset = c("mouse_discovery", "mouse_discovery", "pig_early", "pig_early", "pig_chronic"),
  comparison = c("mouse_1W", "mouse_4W", "pig_early_1W", "pig_early_4W", "pig_chronic_52W"),
  species = c("Mus musculus", "Mus musculus", "Sus scrofa", "Sus scrofa", "Sus scrofa"),
  path = c(
    "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step13_Figure3B_mouse_Hallmark_GSEA/step13_mouse_GSEA_rank_1W.csv",
    "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step13_Figure3B_mouse_Hallmark_GSEA/step13_mouse_GSEA_rank_4W.csv",
    "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/step20_current78_ranked_genes_t7_vs_Control.csv",
    "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/step20_current78_ranked_genes_t28_vs_Control.csv",
    "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step4_current78_DE_GSEA/chronic_step4_ranked_genes_ACLT_alone_52W_vs_Control_52W.csv"
  ),
  symbol_col = "gene_symbol",
  rank_col = c("rank_stat", "rank_stat", "rank_stat", "rank_stat", "rank_statistic"),
  seed = c(22001L, 22002L, 22003L, 22004L, 22005L),
  stringsAsFactors = FALSE
)

collections <- data.frame(
  collection = c("C2", "C5", "C2", "C2"),
  subcollection = c("CP:REACTOME", "GO:BP", "CGP", "CP:WIKIPATHWAYS"),
  family = c("Reactome", "GO_Biological_Process", "Published_Gene_Signatures", "WikiPathways"),
  stringsAsFactors = FALSE
)

read_rank <- function(spec) {
  if (!file.exists(spec$path)) stop("Missing rank file: ", spec$path, call. = FALSE)
  x <- utils::read.csv(spec$path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c(spec$symbol_col, spec$rank_col) %in% names(x))) {
    stop("Rank columns missing in ", spec$path, call. = FALSE)
  }
  symbol <- trimws(as.character(x[[spec$symbol_col]]))
  score <- suppressWarnings(as.numeric(x[[spec$rank_col]]))
  keep <- !is.na(symbol) & nzchar(symbol) & is.finite(score)
  symbol <- symbol[keep]
  score <- score[keep]
  ord <- order(-abs(score), symbol)
  symbol <- symbol[ord]
  score <- score[ord]
  keep_unique <- !duplicated(symbol)
  score <- score[keep_unique]
  names(score) <- symbol[keep_unique]
  sort(score, decreasing = TRUE)
}

clean_fgsea <- function(x, spec, family, collection, subcollection, n_loaded, n_tested) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (nrow(x) == 0L) return(x)
  if ("leadingEdge" %in% names(x)) {
    x$leadingEdge <- vapply(x$leadingEdge, paste, collapse = ";", character(1))
  }
  x$dataset <- spec$dataset
  x$comparison <- spec$comparison
  x$species <- spec$species
  x$family <- family
  x$collection <- collection
  x$subcollection <- subcollection
  x$n_gene_sets_loaded <- n_loaded
  x$n_gene_sets_tested <- n_tested
  x[, c("dataset", "comparison", "species", "family", "collection", "subcollection",
        "n_gene_sets_loaded", "n_gene_sets_tested", setdiff(names(x),
        c("dataset", "comparison", "species", "family", "collection", "subcollection",
          "n_gene_sets_loaded", "n_gene_sets_tested"))), drop = FALSE]
}

species_cache <- list()
all_results <- list()
run_audit <- list()
membership_audit <- list()
result_i <- 0L
audit_i <- 0L
member_i <- 0L

for (sp in unique(rank_specs$species)) {
  message("Loading MSigDB for ", sp)
  species_cache[[sp]] <- msigdbr::msigdbr(species = sp)
}

for (i in seq_len(nrow(rank_specs))) {
  spec <- rank_specs[i, , drop = FALSE]
  stats_vec <- read_rank(spec)
  db <- species_cache[[spec$species]]
  for (j in seq_len(nrow(collections))) {
    cc <- collections[j, , drop = FALSE]
    keep_db <- db$gs_collection == cc$collection & db$gs_subcollection == cc$subcollection
    subdb <- db[keep_db, c("gs_name", "gene_symbol"), drop = FALSE]
    subdb <- subdb[!is.na(subdb$gene_symbol) & nzchar(subdb$gene_symbol), , drop = FALSE]
    pathways <- split(subdb$gene_symbol, subdb$gs_name)
    pathways <- lapply(pathways, unique)
    overlap_n <- vapply(pathways, function(g) sum(g %in% names(stats_vec)), integer(1))
    testable <- overlap_n >= min_size & overlap_n <= max_size
    pathways_test <- pathways[testable]

    set.seed(spec$seed + j * 1000L)
    fg <- fgsea::fgsea(
      pathways = pathways_test,
      stats = stats_vec,
      minSize = min_size,
      maxSize = max_size,
      eps = fgsea_eps,
      nproc = nproc
    )
    result_i <- result_i + 1L
    all_results[[result_i]] <- clean_fgsea(
      fg, spec, cc$family, cc$collection, cc$subcollection,
      length(pathways), length(pathways_test)
    )

    audit_i <- audit_i + 1L
    run_audit[[audit_i]] <- data.frame(
      dataset = spec$dataset,
      comparison = spec$comparison,
      species = spec$species,
      family = cc$family,
      collection = cc$collection,
      subcollection = cc$subcollection,
      ranked_genes = length(stats_vec),
      gene_sets_loaded = length(pathways),
      gene_sets_tested = length(pathways_test),
      minSize = min_size,
      maxSize = max_size,
      eps = fgsea_eps,
      nproc = nproc,
      seed = spec$seed + j * 1000L,
      rank_file = spec$path,
      stringsAsFactors = FALSE
    )

    sen_names <- names(pathways)[grepl(panel_name_regex, names(pathways), ignore.case = TRUE)]
    if (length(sen_names) > 0L) {
      member_i <- member_i + 1L
      membership_audit[[member_i]] <- do.call(rbind, lapply(sen_names, function(nm) {
        data.frame(species = spec$species, family = cc$family, pathway = nm,
                   database_genes = length(pathways[[nm]]),
                   ranked_genes_overlap = sum(pathways[[nm]] %in% names(stats_vec)),
                   tested = nm %in% names(pathways_test), stringsAsFactors = FALSE)
      }))
    }
    message(spec$comparison, " / ", cc$family, ": ", length(pathways_test), " sets tested")
  }
}

full <- do.call(rbind, all_results)
rownames(full) <- NULL
target <- full[grepl(panel_name_regex, full$pathway, ignore.case = TRUE), , drop = FALSE]
target$direct_senescence_term <- TRUE
target$context_flag <- ifelse(grepl("COVID19", target$pathway, ignore.case = TRUE),
                              "disease_context_specific", "general_or_mechanistic_senescence")

## Preserve fgsea's within-subcollection adjustment as a sensitivity result.
target$full_subcollection_FDR <- target$padj
target$full_subcollection_FDR_status <- fdr_status(target$full_subcollection_FDR, target$pval)

## Formal primary analysis: one prespecified senescence/SASP family per comparison,
## pooling all testable terms from Reactome, GO:BP, CGP, and WikiPathways.
target$panel_FDR <- NA_real_
for (cmp in rank_specs$comparison) {
  idx <- which(target$comparison == cmp & is.finite(target$pval))
  target$panel_FDR[idx] <- stats::p.adjust(target$pval[idx], method = panel_fdr_method)
}
target$primary_FDR <- target$panel_FDR
target$primary_FDR_status <- fdr_status(target$primary_FDR, target$pval)
target$primary_FDR_scope <- "BH within comparison across all testable predefined SENESCEN|SASP terms"

## Sensitivity analysis: a single BH family across every testable panel term and comparison.
target$global_135_FDR <- NA_real_
global_idx <- which(is.finite(target$pval))
target$global_135_FDR[global_idx] <- stats::p.adjust(target$pval[global_idx], method = panel_fdr_method)
target$global_135_FDR_status <- fdr_status(target$global_135_FDR, target$pval)

## Backward-compatible status field now follows the formal primary panel FDR.
target$FDR_status <- target$primary_FDR_status
target <- target[order(target$family, target$pathway, target$comparison), , drop = FALSE]

audit <- do.call(rbind, run_audit)
membership <- unique(do.call(rbind, membership_audit))
membership <- membership[order(membership$species, membership$family, membership$pathway), , drop = FALSE]

coverage_summary <- aggregate(
  tested ~ species + family + pathway,
  data = membership,
  FUN = function(x) all(x)
)
names(coverage_summary)[names(coverage_summary) == "tested"] <- "tested_in_all_comparisons_for_species"

n_predefined_terms <- length(unique(c(membership$pathway, target$pathway)))
primary_summary <- do.call(rbind, lapply(rank_specs$comparison, function(cmp) {
  d <- target[target$comparison == cmp, , drop = FALSE]
  data.frame(
    comparison = cmp,
    n_predefined_terms = n_predefined_terms,
    n_testable_terms = nrow(d),
    n_not_tested_terms = n_predefined_terms - nrow(d),
    n_primary_FDR_lt_0.05 = sum(d$primary_FDR < 0.05, na.rm = TRUE),
    n_primary_FDR_positive = sum(d$primary_FDR < 0.05 & d$NES > 0, na.rm = TRUE),
    n_primary_FDR_negative = sum(d$primary_FDR < 0.05 & d$NES < 0, na.rm = TRUE),
    n_nominal_only = sum(d$pval < 0.05 & d$primary_FDR >= 0.05, na.rm = TRUE),
    n_global_135_FDR_lt_0.05 = sum(d$global_135_FDR < 0.05, na.rm = TRUE),
    n_full_subcollection_FDR_lt_0.05 = sum(d$full_subcollection_FDR < 0.05, na.rm = TRUE),
    primary_correction = "BH within comparison across all testable predefined SENESCEN|SASP terms",
    sensitivity_corrections = "BH across all 135 testable panel rows; fgsea BH within each full subcollection",
    stringsAsFactors = FALSE
  )
}))

cross_comparison <- do.call(rbind, lapply(split(target, interaction(target$family, target$pathway, drop = TRUE)), function(d) {
  data.frame(
    family = d$family[1],
    pathway = d$pathway[1],
    n_comparisons_tested = nrow(d),
    n_positive_NES = sum(d$NES > 0, na.rm = TRUE),
    n_negative_NES = sum(d$NES < 0, na.rm = TRUE),
    n_primary_FDR_lt_0.05 = sum(d$primary_FDR < 0.05, na.rm = TRUE),
    n_FDR_lt_0.05 = sum(d$primary_FDR < 0.05, na.rm = TRUE),
    n_nominal_p_lt_0.05 = sum(d$pval < 0.05, na.rm = TRUE),
    min_primary_FDR = if (all(is.na(d$primary_FDR))) NA_real_ else min(d$primary_FDR, na.rm = TRUE),
    min_FDR = if (all(is.na(d$primary_FDR))) NA_real_ else min(d$primary_FDR, na.rm = TRUE),
    comparisons_primary_FDR_supported = paste(d$comparison[!is.na(d$primary_FDR) & d$primary_FDR < 0.05], collapse = ";"),
    comparisons_FDR_supported = paste(d$comparison[!is.na(d$primary_FDR) & d$primary_FDR < 0.05], collapse = ";"),
    n_global_135_FDR_lt_0.05 = sum(d$global_135_FDR < 0.05, na.rm = TRUE),
    min_global_135_FDR = if (all(is.na(d$global_135_FDR))) NA_real_ else min(d$global_135_FDR, na.rm = TRUE),
    n_full_subcollection_FDR_lt_0.05 = sum(d$full_subcollection_FDR < 0.05, na.rm = TRUE),
    min_full_subcollection_FDR = if (all(is.na(d$full_subcollection_FDR))) NA_real_ else min(d$full_subcollection_FDR, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
cross_comparison <- cross_comparison[order(-cross_comparison$n_primary_FDR_lt_0.05,
                                           cross_comparison$family, cross_comparison$pathway), , drop = FALSE]

utils::write.csv(full, file.path(out_dir, "extended_msigdb_full_subcollection_gsea_results.csv"), row.names = FALSE)
utils::write.csv(target, file.path(out_dir, "extended_msigdb_direct_senescence_results.csv"), row.names = FALSE)
utils::write.csv(cross_comparison, file.path(out_dir, "extended_msigdb_senescence_cross_comparison_summary.csv"), row.names = FALSE)
utils::write.csv(primary_summary, file.path(out_dir, "extended_msigdb_primary_panel_fdr_summary.csv"), row.names = FALSE)
utils::write.csv(audit, file.path(out_dir, "extended_msigdb_run_audit.csv"), row.names = FALSE)
utils::write.csv(membership, file.path(out_dir, "extended_msigdb_senescence_membership_coverage.csv"), row.names = FALSE)
utils::write.csv(coverage_summary, file.path(out_dir, "extended_msigdb_senescence_testability_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "extended_msigdb_sessionInfo.txt"))
writeLines(c(
  "Formal primary analysis: prespecified direct senescence/SASP panel.",
  paste0("Panel selection rule: official MSigDB gene-set name matches /", panel_name_regex, "/i; selection is independent of enrichment results."),
  paste0("Total unique panel terms: ", n_predefined_terms, "."),
  "Primary multiplicity family: all testable panel terms from Reactome, GO:BP, CGP, and WikiPathways within each prespecified comparison.",
  paste0("Primary correction: ", panel_fdr_method, "."),
  "Sensitivity 1: BH across all 135 testable panel term-by-comparison rows.",
  "Sensitivity 2: fgsea BH adjustment within each complete MSigDB subcollection.",
  "The five comparisons are treated as distinct prespecified biological contrasts; significance counts are not formal tests of between-contrast differences."
), file.path(out_dir, "extended_msigdb_primary_analysis_specification.txt"))

cat("Full GSEA rows:", nrow(full), "\n")
cat("Direct senescence/SASP rows:", nrow(target), "\n")
cat("Direct terms:", length(unique(target$pathway)), "\n")
cat("Primary panel-FDR-supported rows:", sum(target$primary_FDR < 0.05, na.rm = TRUE), "\n")
cat("Global-135-FDR-supported rows:", sum(target$global_135_FDR < 0.05, na.rm = TRUE), "\n")
cat("Full-subcollection-FDR-supported rows:", sum(target$full_subcollection_FDR < 0.05, na.rm = TRUE), "\n")
print(primary_summary, row.names = FALSE)
print(cross_comparison, row.names = FALSE)
