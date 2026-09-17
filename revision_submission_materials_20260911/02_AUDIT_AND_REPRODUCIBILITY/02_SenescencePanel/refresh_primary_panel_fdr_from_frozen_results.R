options(stringsAsFactors = FALSE)

## Deterministic postprocessing for Reviewer 2 senescence/SASP panel.
## Uses the frozen full-subcollection GSEA results and changes only multiplicity fields.

out_dir <- "D:/workspace/ACLsenescence2_reviewer2_senescence_pathways"
full_file <- file.path(out_dir, "extended_msigdb_full_subcollection_gsea_results.csv")
membership_file <- file.path(out_dir, "extended_msigdb_senescence_membership_coverage.csv")

if (!file.exists(full_file)) stop("Missing frozen full GSEA results: ", full_file, call. = FALSE)
if (!file.exists(membership_file)) stop("Missing panel membership audit: ", membership_file, call. = FALSE)

panel_name_regex <- "SENESCEN|SASP"
panel_fdr_method <- "BH"
comparison_order <- c("mouse_1W", "mouse_4W", "pig_early_1W", "pig_early_4W", "pig_chronic_52W")

fdr_status <- function(fdr, pval) {
  ifelse(is.na(fdr), "not_available",
         ifelse(fdr < 0.05, "FDR_lt_0.05",
                ifelse(!is.na(pval) & pval < 0.05, "nominal_only", "not_significant")))
}

full <- utils::read.csv(full_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("comparison", "family", "pathway", "pval", "padj", "NES")
if (!all(required_columns %in% names(full))) {
  stop("Frozen full GSEA table lacks required columns: ",
       paste(setdiff(required_columns, names(full)), collapse = ", "), call. = FALSE)
}

target <- full[grepl(panel_name_regex, full$pathway, ignore.case = TRUE), , drop = FALSE]
target$direct_senescence_term <- TRUE
target$context_flag <- ifelse(grepl("COVID19", target$pathway, ignore.case = TRUE),
                              "disease_context_specific", "general_or_mechanistic_senescence")

target$full_subcollection_FDR <- target$padj
target$full_subcollection_FDR_status <- fdr_status(target$full_subcollection_FDR, target$pval)

target$panel_FDR <- NA_real_
for (cmp in comparison_order) {
  idx <- which(target$comparison == cmp & is.finite(target$pval))
  target$panel_FDR[idx] <- stats::p.adjust(target$pval[idx], method = panel_fdr_method)
}
target$primary_FDR <- target$panel_FDR
target$primary_FDR_status <- fdr_status(target$primary_FDR, target$pval)
target$primary_FDR_scope <- "BH within comparison across all testable predefined SENESCEN|SASP terms"

target$global_135_FDR <- NA_real_
global_idx <- which(is.finite(target$pval))
target$global_135_FDR[global_idx] <- stats::p.adjust(target$pval[global_idx], method = panel_fdr_method)
target$global_135_FDR_status <- fdr_status(target$global_135_FDR, target$pval)

## Backward-compatible status field now follows the formal primary result.
target$FDR_status <- target$primary_FDR_status
target <- target[order(target$family, target$pathway,
                       match(target$comparison, comparison_order)), , drop = FALSE]

membership <- utils::read.csv(membership_file, check.names = FALSE, stringsAsFactors = FALSE)
n_predefined_terms <- length(unique(c(membership$pathway, target$pathway)))

primary_summary <- do.call(rbind, lapply(comparison_order, function(cmp) {
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

cross_comparison <- do.call(rbind, lapply(
  split(target, interaction(target$family, target$pathway, drop = TRUE)),
  function(d) {
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
  }
))
cross_comparison <- cross_comparison[order(-cross_comparison$n_primary_FDR_lt_0.05,
                                           cross_comparison$family,
                                           cross_comparison$pathway), , drop = FALSE]

utils::write.csv(target, file.path(out_dir, "extended_msigdb_direct_senescence_results.csv"), row.names = FALSE)
utils::write.csv(primary_summary, file.path(out_dir, "extended_msigdb_primary_panel_fdr_summary.csv"), row.names = FALSE)
utils::write.csv(cross_comparison, file.path(out_dir, "extended_msigdb_senescence_cross_comparison_summary.csv"), row.names = FALSE)

writeLines(c(
  "Formal primary analysis: prespecified direct senescence/SASP panel.",
  paste0("Panel selection rule: official MSigDB gene-set name matches /", panel_name_regex,
         "/i; selection is independent of enrichment results."),
  paste0("Total unique panel terms: ", n_predefined_terms, "."),
  "Primary multiplicity family: all testable panel terms from Reactome, GO:BP, CGP, and WikiPathways within each prespecified comparison.",
  paste0("Primary correction: ", panel_fdr_method, "."),
  "Sensitivity 1: BH across all 135 testable panel term-by-comparison rows.",
  "Sensitivity 2: fgsea BH adjustment within each complete MSigDB subcollection.",
  "The five comparisons are distinct prespecified biological contrasts; significance counts are not formal tests of between-contrast differences."
), file.path(out_dir, "extended_msigdb_primary_analysis_specification.txt"))
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "extended_msigdb_primary_panel_fdr_sessionInfo.txt"))

stopifnot(nrow(target) == 135L, n_predefined_terms == 33L)
stopifnot(identical(primary_summary$n_testable_terms, c(28L, 28L, 27L, 26L, 26L)))
stopifnot(identical(primary_summary$n_primary_FDR_lt_0.05, c(2L, 3L, 5L, 4L, 3L)))

cat("Frozen full GSEA rows:", nrow(full), "\n")
cat("Direct panel rows:", nrow(target), "\n")
cat("Unique predefined terms:", n_predefined_terms, "\n")
cat("Primary panel-FDR-supported rows:", sum(target$primary_FDR < 0.05, na.rm = TRUE), "\n")
cat("Global-135-FDR-supported rows:", sum(target$global_135_FDR < 0.05, na.rm = TRUE), "\n")
cat("Full-subcollection-FDR-supported rows:", sum(target$full_subcollection_FDR < 0.05, na.rm = TRUE), "\n")
print(primary_summary, row.names = FALSE)
