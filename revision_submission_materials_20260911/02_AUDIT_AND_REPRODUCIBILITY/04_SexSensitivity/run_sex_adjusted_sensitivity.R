options(stringsAsFactors = FALSE)

project_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"
output_dir <- "D:/workspace/ACLsenescence2_reviewer_minor1_sex_adjusted"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

strict_fdr <- 0.05
strict_abs_logfc <- 1

input_paths <- c(
  counts_1W = file.path(project_dir, "02_expression_matrix", "step04_expr_1W.csv"),
  counts_4W = file.path(project_dir, "02_expression_matrix", "step04_expr_4W.csv"),
  annotation_1W = file.path(project_dir, "01_metadata", "step04_anno_1W.csv"),
  annotation_4W = file.path(project_dir, "01_metadata", "step04_anno_4W.csv"),
  primary_DE_1W = file.path(project_dir, "03_DE_analysis", "step05_DE_1W_ACLR_vs_Contra_limma_voom.csv"),
  primary_DE_4W = file.path(project_dir, "03_DE_analysis", "step05_DE_4W_ACLR_vs_Contra_limma_voom.csv"),
  primary_persistent = file.path(
    project_dir, "07_tables", "step07_strict_DEG_upset_persistent",
    "step07_persistent_direction_consistent_genes.csv"
  )
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0) {
  stop("Missing input files: ", paste(missing_inputs, collapse = "; "))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

standardize_sex <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- ifelse(z %in% c("F", "FEMALE"), "Female",
                ifelse(z %in% c("M", "MALE"), "Male", NA_character_))
  factor(out, levels = c("Female", "Male"))
}

read_counts <- function(path, timepoint) {
  raw <- read.csv(path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  numeric_df <- as.data.frame(lapply(raw, function(z) suppressWarnings(as.numeric(as.character(z)))))
  rownames(numeric_df) <- rownames(raw)
  mat <- as.matrix(numeric_df)
  storage.mode(mat) <- "numeric"
  na_rows <- rowSums(is.na(mat)) > 0
  mat <- mat[!na_rows, , drop = FALSE]
  if (any(mat < 0)) stop(timepoint, ": negative counts detected")
  list(
    counts = mat,
    audit = data.frame(
      timepoint = timepoint,
      original_genes = nrow(raw),
      genes_removed_due_to_any_NA = sum(na_rows),
      genes_after_NA_removal = nrow(mat),
      samples = ncol(mat)
    )
  )
}

prepare_annotation <- function(path, counts, timepoint) {
  anno <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("sample_id", "mouse_id", "sex", "treatment")
  if (!all(required %in% colnames(anno))) {
    stop(timepoint, ": annotation missing columns: ", paste(setdiff(required, colnames(anno)), collapse = ", "))
  }
  anno$sex <- standardize_sex(anno$sex)
  anno$treatment <- factor(anno$treatment, levels = c("Contra", "ACLR"))
  anno$mouse_id <- factor(anno$mouse_id)
  anno <- anno[match(colnames(counts), anno$sample_id), , drop = FALSE]
  if (any(is.na(anno$sample_id))) stop(timepoint, ": annotation does not match count columns")
  if (any(is.na(anno$sex))) stop(timepoint, ": unrecognized or missing sex")
  pair_table <- table(anno$mouse_id, anno$treatment)
  if (!all(pair_table[, "Contra"] == 1L & pair_table[, "ACLR"] == 1L)) {
    stop(timepoint, ": each mouse must have exactly one Contra and one ACLR sample")
  }
  sex_per_mouse <- tapply(as.character(anno$sex), anno$mouse_id, function(z) length(unique(z)))
  if (!all(sex_per_mouse == 1L)) stop(timepoint, ": sex is not constant within mouse")
  anno
}

fit_limma_voom <- function(counts, anno, formula_text, timepoint, model_label, forced_gene_ids = NULL) {
  design <- model.matrix(as.formula(formula_text), data = anno)
  rownames(design) <- anno$sample_id
  if (qr(design)$rank != ncol(design)) stop(timepoint, " ", model_label, ": design is not full rank")

  dge <- edgeR::DGEList(counts = counts, samples = anno)
  if (is.null(forced_gene_ids)) {
    keep <- edgeR::filterByExpr(dge, design = design)
    filter_rule <- paste0("filterByExpr using ", formula_text)
  } else {
    missing_forced <- setdiff(forced_gene_ids, rownames(dge))
    if (length(missing_forced) > 0) {
      stop(timepoint, " ", model_label, ": forced gene universe contains missing genes")
    }
    keep <- rownames(dge) %in% forced_gene_ids
    filter_rule <- "frozen primary filterByExpr gene universe"
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  v1 <- limma::voom(dge, design = design, plot = FALSE)
  corfit <- limma::duplicateCorrelation(v1, design = design, block = anno$mouse_id)
  v2 <- limma::voom(
    dge, design = design, plot = FALSE,
    block = anno$mouse_id, correlation = corfit$consensus
  )
  fit <- limma::lmFit(
    v2, design = design,
    block = anno$mouse_id, correlation = corfit$consensus
  )
  fit <- limma::eBayes(fit)
  if (!("treatmentACLR" %in% colnames(fit$coefficients))) {
    stop(timepoint, " ", model_label, ": treatmentACLR coefficient missing")
  }
  tt <- limma::topTable(
    fit, coef = "treatmentACLR", number = Inf,
    sort.by = "P", adjust.method = "BH"
  )
  tt$gene_id <- rownames(tt)
  rownames(tt) <- NULL
  colnames(tt)[colnames(tt) == "adj.P.Val"] <- "FDR"
  tt$strict_DEG <- tt$FDR < strict_fdr & abs(tt$logFC) > strict_abs_logfc
  tt$direction <- ifelse(tt$strict_DEG & tt$logFC > 0, "Up",
                         ifelse(tt$strict_DEG & tt$logFC < 0, "Down", "Not_strict"))
  tt <- tt[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "FDR", "B", "strict_DEG", "direction")]

  design_out <- data.frame(sample_id = rownames(design), design, check.names = FALSE)
  list(
    table = tt,
    design = design_out,
    kept_gene_ids = rownames(v2$E),
    filter_rule = filter_rule,
    consensus_correlation = unname(corfit$consensus),
    n_genes_after_filter = nrow(v2$E)
  )
}

add_annotation <- function(tt, annotation_map) {
  idx <- match(tt$gene_id, annotation_map$gene_id)
  tt$ENTREZID <- annotation_map$ENTREZID[idx]
  tt$SYMBOL <- annotation_map$SYMBOL[idx]
  tt[, c("gene_id", "ENTREZID", "SYMBOL", setdiff(colnames(tt), c("gene_id", "ENTREZID", "SYMBOL")))]
}

compare_models <- function(primary, adjusted, timepoint) {
  p <- primary[, c("gene_id", "ENTREZID", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "FDR", "strict_DEG", "direction")]
  colnames(p)[4:10] <- paste0("primary_", colnames(p)[4:10])
  a <- adjusted[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "FDR", "strict_DEG", "direction")]
  colnames(a)[2:8] <- paste0("sex_adjusted_", colnames(a)[2:8])
  z <- merge(p, a, by = "gene_id", all = TRUE, sort = FALSE)
  z$timepoint <- timepoint
  z$delta_logFC_adjusted_minus_primary <- z$sex_adjusted_logFC - z$primary_logFC
  z$abs_delta_logFC <- abs(z$delta_logFC_adjusted_minus_primary)
  z$same_logFC_direction <- ifelse(
    is.na(z$primary_logFC) | is.na(z$sex_adjusted_logFC), NA,
    sign(z$primary_logFC) == sign(z$sex_adjusted_logFC)
  )
  z$primary_FDR_lt_0.05 <- !is.na(z$primary_FDR) & z$primary_FDR < strict_fdr
  z$sex_adjusted_FDR_lt_0.05 <- !is.na(z$sex_adjusted_FDR) & z$sex_adjusted_FDR < strict_fdr
  z$FDR_supported_in_both <- z$primary_FDR_lt_0.05 & z$sex_adjusted_FDR_lt_0.05
  z$strict_in_both <- (!is.na(z$primary_strict_DEG) & z$primary_strict_DEG) &
    (!is.na(z$sex_adjusted_strict_DEG) & z$sex_adjusted_strict_DEG)
  z
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok], method = method)
}

summarize_comparison <- function(comp, primary_fit, adjusted_fit, timepoint) {
  common <- !is.na(comp$primary_logFC) & !is.na(comp$sex_adjusted_logFC)
  p_fdr <- comp$primary_FDR_lt_0.05
  a_fdr <- comp$sex_adjusted_FDR_lt_0.05
  p_strict <- !is.na(comp$primary_strict_DEG) & comp$primary_strict_DEG
  a_strict <- !is.na(comp$sex_adjusted_strict_DEG) & comp$sex_adjusted_strict_DEG
  fdr_union <- sum(p_fdr | a_fdr)
  strict_union <- sum(p_strict | a_strict)
  data.frame(
    timepoint = timepoint,
    animals = 12L,
    female_animals = 6L,
    male_animals = 6L,
    samples = 24L,
    primary_design = "~ treatment",
    sensitivity_design = "~ sex + treatment",
    block_factor = "mouse_id",
    primary_genes_tested = sum(!is.na(comp$primary_logFC)),
    adjusted_genes_tested = sum(!is.na(comp$sex_adjusted_logFC)),
    common_genes = sum(common),
    pearson_logFC = safe_cor(comp$primary_logFC, comp$sex_adjusted_logFC, "pearson"),
    spearman_logFC = safe_cor(comp$primary_logFC, comp$sex_adjusted_logFC, "spearman"),
    pearson_t_statistic = safe_cor(comp$primary_t, comp$sex_adjusted_t, "pearson"),
    median_abs_delta_logFC = median(comp$abs_delta_logFC[common], na.rm = TRUE),
    p95_abs_delta_logFC = unname(quantile(comp$abs_delta_logFC[common], 0.95, na.rm = TRUE)),
    max_abs_delta_logFC = max(comp$abs_delta_logFC[common], na.rm = TRUE),
    all_gene_direction_concordance = mean(comp$same_logFC_direction[common], na.rm = TRUE),
    primary_FDR_lt_0.05 = sum(p_fdr),
    adjusted_FDR_lt_0.05 = sum(a_fdr),
    FDR_lt_0.05_intersection = sum(p_fdr & a_fdr),
    primary_FDR_retention = sum(p_fdr & a_fdr) / sum(p_fdr),
    FDR_jaccard = ifelse(fdr_union > 0, sum(p_fdr & a_fdr) / fdr_union, NA_real_),
    primary_strict_DEG = sum(p_strict),
    adjusted_strict_DEG = sum(a_strict),
    strict_DEG_intersection = sum(p_strict & a_strict),
    primary_strict_retention = sum(p_strict & a_strict) / sum(p_strict),
    strict_DEG_jaccard = ifelse(strict_union > 0, sum(p_strict & a_strict) / strict_union, NA_real_),
    primary_strict_up = sum(p_strict & comp$primary_logFC > 0, na.rm = TRUE),
    adjusted_strict_up = sum(a_strict & comp$sex_adjusted_logFC > 0, na.rm = TRUE),
    primary_strict_down = sum(p_strict & comp$primary_logFC < 0, na.rm = TRUE),
    adjusted_strict_down = sum(a_strict & comp$sex_adjusted_logFC < 0, na.rm = TRUE),
    primary_duplicateCorrelation = primary_fit$consensus_correlation,
    adjusted_duplicateCorrelation = adjusted_fit$consensus_correlation
  )
}

## Load inputs.
count_1w <- read_counts(input_paths[["counts_1W"]], "1W")
count_4w <- read_counts(input_paths[["counts_4W"]], "4W")
anno_1w <- prepare_annotation(input_paths[["annotation_1W"]], count_1w$counts, "1W")
anno_4w <- prepare_annotation(input_paths[["annotation_4W"]], count_4w$counts, "4W")
primary_1w <- read.csv(input_paths[["primary_DE_1W"]], stringsAsFactors = FALSE, check.names = FALSE)
primary_4w <- read.csv(input_paths[["primary_DE_4W"]], stringsAsFactors = FALSE, check.names = FALSE)
original_persistent <- read.csv(input_paths[["primary_persistent"]], stringsAsFactors = FALSE, check.names = FALSE)

annotation_map <- unique(rbind(
  primary_1w[, c("gene_id", "ENTREZID", "SYMBOL")],
  primary_4w[, c("gene_id", "ENTREZID", "SYMBOL")]
))
annotation_map <- annotation_map[!duplicated(annotation_map$gene_id), , drop = FALSE]

## Reproduce the primary model as a QA control, then fit the requested sensitivity model.
fit_primary_1w <- fit_limma_voom(count_1w$counts, anno_1w, "~ treatment", "1W", "primary reproduction")
fit_primary_4w <- fit_limma_voom(count_4w$counts, anno_4w, "~ treatment", "4W", "primary reproduction")
fit_adjusted_1w <- fit_limma_voom(
  count_1w$counts, anno_1w, "~ sex + treatment", "1W", "sex-adjusted",
  forced_gene_ids = fit_primary_1w$kept_gene_ids
)
fit_adjusted_4w <- fit_limma_voom(
  count_4w$counts, anno_4w, "~ sex + treatment", "4W", "sex-adjusted",
  forced_gene_ids = fit_primary_4w$kept_gene_ids
)

adjusted_1w <- add_annotation(fit_adjusted_1w$table, annotation_map)
adjusted_4w <- add_annotation(fit_adjusted_4w$table, annotation_map)

comp_1w <- compare_models(primary_1w, adjusted_1w, "1W")
comp_4w <- compare_models(primary_4w, adjusted_4w, "4W")
summary_table <- rbind(
  summarize_comparison(comp_1w, fit_primary_1w, fit_adjusted_1w, "1W"),
  summarize_comparison(comp_4w, fit_primary_4w, fit_adjusted_4w, "4W")
)

## Persistent-gene sensitivity using the same strict thresholds at both timepoints.
p1 <- adjusted_1w[, c("gene_id", "ENTREZID", "SYMBOL", "logFC", "FDR", "strict_DEG", "direction")]
p4 <- adjusted_4w[, c("gene_id", "logFC", "FDR", "strict_DEG", "direction")]
colnames(p1)[4:7] <- paste0(colnames(p1)[4:7], "_1W_adjusted")
colnames(p4)[2:5] <- paste0(colnames(p4)[2:5], "_4W_adjusted")
adjusted_cross_time <- merge(p1, p4, by = "gene_id", all = FALSE, sort = FALSE)
adjusted_cross_time$is_adjusted_persistent <- adjusted_cross_time$strict_DEG_1W_adjusted &
  adjusted_cross_time$strict_DEG_4W_adjusted &
  adjusted_cross_time$direction_1W_adjusted == adjusted_cross_time$direction_4W_adjusted
adjusted_cross_time$adjusted_persistent_direction <- ifelse(
  adjusted_cross_time$is_adjusted_persistent,
  adjusted_cross_time$direction_1W_adjusted,
  "Not_persistent"
)
adjusted_persistent <- adjusted_cross_time[adjusted_cross_time$is_adjusted_persistent, , drop = FALSE]

orig_ids <- unique(as.character(original_persistent$gene_id))
adj_ids <- unique(as.character(adjusted_persistent$gene_id))
orig_up <- unique(as.character(original_persistent$gene_id[original_persistent$direction_persistent == "Up"]))
orig_down <- unique(as.character(original_persistent$gene_id[original_persistent$direction_persistent == "Down"]))
adj_up <- unique(as.character(adjusted_persistent$gene_id[adjusted_persistent$adjusted_persistent_direction == "Up"]))
adj_down <- unique(as.character(adjusted_persistent$gene_id[adjusted_persistent$adjusted_persistent_direction == "Down"]))

original_persistent_audit <- merge(
  original_persistent,
  adjusted_cross_time,
  by = "gene_id", all.x = TRUE, sort = FALSE,
  suffixes = c("_primary", "")
)
original_persistent_audit$retained_as_adjusted_persistent <-
  original_persistent_audit$gene_id %in% adj_ids
original_persistent_audit$retained_same_direction <-
  original_persistent_audit$retained_as_adjusted_persistent &
  original_persistent_audit$direction_persistent == original_persistent_audit$adjusted_persistent_direction

persistent_summary <- data.frame(
  metric = c(
    "primary_persistent_total", "primary_persistent_up", "primary_persistent_down",
    "sex_adjusted_persistent_total", "sex_adjusted_persistent_up", "sex_adjusted_persistent_down",
    "persistent_intersection_total", "primary_persistent_retention_fraction",
    "persistent_jaccard", "primary_up_retained_as_adjusted_up",
    "primary_up_retention_fraction", "primary_down_retained_as_adjusted_down",
    "primary_down_retention_fraction"
  ),
  value = c(
    length(orig_ids), length(orig_up), length(orig_down),
    length(adj_ids), length(adj_up), length(adj_down),
    length(intersect(orig_ids, adj_ids)), length(intersect(orig_ids, adj_ids)) / length(orig_ids),
    length(intersect(orig_ids, adj_ids)) / length(union(orig_ids, adj_ids)),
    length(intersect(orig_up, adj_up)), length(intersect(orig_up, adj_up)) / length(orig_up),
    length(intersect(orig_down, adj_down)), length(intersect(orig_down, adj_down)) / length(orig_down)
  )
)

## Primary-model reproduction QA against the frozen Step05 tables.
primary_repro_compare <- function(frozen, reproduced, timepoint) {
  m <- merge(
    frozen[, c("gene_id", "logFC", "t", "P.Value", "FDR")],
    reproduced$table[, c("gene_id", "logFC", "t", "P.Value", "FDR")],
    by = "gene_id", all = TRUE, suffixes = c("_frozen", "_reproduced")
  )
  common <- complete.cases(m[, c("logFC_frozen", "logFC_reproduced")])
  data.frame(
    timepoint = timepoint,
    frozen_genes = nrow(frozen),
    reproduced_genes = nrow(reproduced$table),
    common_genes = sum(common),
    identical_gene_set = setequal(frozen$gene_id, reproduced$table$gene_id),
    max_abs_logFC_difference = max(abs(m$logFC_frozen[common] - m$logFC_reproduced[common])),
    max_abs_t_difference = max(abs(m$t_frozen[common] - m$t_reproduced[common])),
    max_abs_P_difference = max(abs(m$P.Value_frozen[common] - m$P.Value_reproduced[common])),
    max_abs_FDR_difference = max(abs(m$FDR_frozen[common] - m$FDR_reproduced[common]))
  )
}
reproduction_audit <- rbind(
  primary_repro_compare(primary_1w, fit_primary_1w, "1W"),
  primary_repro_compare(primary_4w, fit_primary_4w, "4W")
)

## Write outputs.
write_csv(adjusted_1w, file.path(output_dir, "sex_adjusted_DE_1W.csv"))
write_csv(adjusted_4w, file.path(output_dir, "sex_adjusted_DE_4W.csv"))
write_csv(fit_adjusted_1w$design, file.path(output_dir, "sex_adjusted_design_matrix_1W.csv"))
write_csv(fit_adjusted_4w$design, file.path(output_dir, "sex_adjusted_design_matrix_4W.csv"))
write_csv(comp_1w, file.path(output_dir, "primary_vs_sex_adjusted_gene_comparison_1W.csv"))
write_csv(comp_4w, file.path(output_dir, "primary_vs_sex_adjusted_gene_comparison_4W.csv"))
write_csv(summary_table, file.path(output_dir, "sex_adjusted_robustness_summary.csv"))
write_csv(adjusted_persistent, file.path(output_dir, "sex_adjusted_persistent_direction_consistent_genes.csv"))
write_csv(original_persistent_audit, file.path(output_dir, "primary_persistent_gene_retention_audit.csv"))
write_csv(persistent_summary, file.path(output_dir, "persistent_gene_robustness_summary.csv"))
write_csv(reproduction_audit, file.path(output_dir, "primary_model_reproduction_QA.csv"))
write_csv(rbind(count_1w$audit, count_4w$audit), file.path(output_dir, "input_count_matrix_audit.csv"))

analysis_specification <- data.frame(
  item = c(
    "analysis_role", "primary_design", "sensitivity_design", "coefficient_tested",
    "block_factor", "within_mouse_correlation", "gene_filter_scope", "normalization",
    "expression_model", "multiple_testing", "strict_DEG_definition",
    "timepoint_handling", "sex_interaction_tested"
  ),
  value = c(
    "sex-adjusted robustness sensitivity analysis", "~ treatment", "~ sex + treatment",
    "treatmentACLR", "mouse_id", "limma duplicateCorrelation consensus",
    "The frozen primary filterByExpr gene universe was retained at each timepoint so the comparison isolates addition of sex to the design matrix.",
    "edgeR TMM", "two-pass limma-voom followed by lmFit and eBayes", "BH within each timepoint",
    "FDR < 0.05 and absolute logFC > 1", "1W and 4W fitted separately", "No"
  )
)
write_csv(analysis_specification, file.path(output_dir, "analysis_specification.csv"))

sample_audit <- do.call(rbind, lapply(list(`1W` = anno_1w, `4W` = anno_4w), function(a) {
  unique_mouse <- a[!duplicated(a$mouse_id), c("mouse_id", "sex"), drop = FALSE]
  data.frame(
    samples = nrow(a), animals = nrow(unique_mouse),
    female_animals = sum(unique_mouse$sex == "Female"),
    male_animals = sum(unique_mouse$sex == "Male"),
    female_samples = sum(a$sex == "Female"), male_samples = sum(a$sex == "Male")
  )
}))
sample_audit$timepoint <- rownames(sample_audit)
rownames(sample_audit) <- NULL
sample_audit <- sample_audit[, c("timepoint", setdiff(colnames(sample_audit), "timepoint"))]
write_csv(sample_audit, file.path(output_dir, "sample_sex_and_pairing_audit.csv"))

input_info <- file.info(input_paths)
input_audit <- data.frame(
  input_name = names(input_paths),
  path = unname(input_paths),
  bytes = input_info$size,
  modified = format(input_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  md5 = unname(tools::md5sum(input_paths))
)
write_csv(input_audit, file.path(output_dir, "input_file_audit.csv"))

qa <- data.frame(
  check = c(
    "1W_sample_balance", "4W_sample_balance", "1W_pair_completeness", "4W_pair_completeness",
    "1W_adjusted_design_columns", "4W_adjusted_design_columns",
    "1W_primary_gene_set_reproduced", "4W_primary_gene_set_reproduced",
    "1W_primary_logFC_reproduced", "4W_primary_logFC_reproduced",
    "1W_adjusted_unique_gene_ids", "4W_adjusted_unique_gene_ids",
    "adjusted_probability_ranges", "original_persistent_unique",
    "persistent_audit_complete"
  ),
  pass = c(
    sample_audit$female_animals[sample_audit$timepoint == "1W"] == 6 & sample_audit$male_animals[sample_audit$timepoint == "1W"] == 6,
    sample_audit$female_animals[sample_audit$timepoint == "4W"] == 6 & sample_audit$male_animals[sample_audit$timepoint == "4W"] == 6,
    all(table(anno_1w$mouse_id, anno_1w$treatment) == 1),
    all(table(anno_4w$mouse_id, anno_4w$treatment) == 1),
    identical(colnames(fit_adjusted_1w$design), c("sample_id", "(Intercept)", "sexMale", "treatmentACLR")),
    identical(colnames(fit_adjusted_4w$design), c("sample_id", "(Intercept)", "sexMale", "treatmentACLR")),
    reproduction_audit$identical_gene_set[reproduction_audit$timepoint == "1W"],
    reproduction_audit$identical_gene_set[reproduction_audit$timepoint == "4W"],
    reproduction_audit$max_abs_logFC_difference[reproduction_audit$timepoint == "1W"] < 1e-10,
    reproduction_audit$max_abs_logFC_difference[reproduction_audit$timepoint == "4W"] < 1e-10,
    !anyDuplicated(adjusted_1w$gene_id),
    !anyDuplicated(adjusted_4w$gene_id),
    all(adjusted_1w$P.Value >= 0 & adjusted_1w$P.Value <= 1 & adjusted_1w$FDR >= 0 & adjusted_1w$FDR <= 1) &
      all(adjusted_4w$P.Value >= 0 & adjusted_4w$P.Value <= 1 & adjusted_4w$FDR >= 0 & adjusted_4w$FDR <= 1),
    !anyDuplicated(original_persistent$gene_id),
    nrow(original_persistent_audit) == nrow(original_persistent)
  )
)
qa$status <- ifelse(qa$pass, "PASS", "FAIL")
write_csv(qa, file.path(output_dir, "sex_adjusted_QA_checks.csv"))

writeLines(
  c(
    "Mouse sex-adjusted sensitivity analysis QA report",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Checks passed: ", sum(qa$pass), "/", nrow(qa)),
    paste0("Overall status: ", ifelse(all(qa$pass), "PASS", "FAIL")),
    "",
    paste(capture.output(print(qa[, c("check", "status")], row.names = FALSE)), collapse = "\n")
  ),
  file.path(output_dir, "sex_adjusted_QA_report.txt")
)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

if (!all(qa$pass)) stop("QA failed; inspect sex_adjusted_QA_checks.csv")

cat("Sex-adjusted sensitivity analysis completed successfully.\n")
print(summary_table)
cat("\nPersistent-gene robustness:\n")
print(persistent_summary)
cat("\nQA: ", sum(qa$pass), "/", nrow(qa), " PASS\n", sep = "")
