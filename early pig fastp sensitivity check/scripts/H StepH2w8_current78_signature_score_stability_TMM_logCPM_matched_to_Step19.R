# StepH2w8_current78: pig early fastp sensitivity - current78 signature score stability
# Corrected version: matches current Step19 TMM-normalized logCPM scoring exactly
# Purpose:
# Recompute current78-derived pig early signature scores using both untrimmed main counts
# and fastp-trimmed branch counts, then evaluate whether fastp trimming materially changes
# signature score conclusions.
#
# Current formal signature:
#   75 detected pig ortholog signature genes
#   65 Up_in_ACLR genes
#   10 Down_in_ACLR genes
#
# Important methodological choices:
#   - The signature gene table is fixed to the current Step18_current78 output.
#   - No automatic search for signature files is used.
#   - Scores follow the current Step19 definition:
#       directional_score = (sum z_up + sum -z_down) / (n_up + n_down)
#     not the older unweighted mean of up_score and down_score_reoriented.
#   - logCPM is computed exactly as the current Step19 implementation:
#       edgeR::DGEList(full count matrix) -> edgeR::calcNormFactors(method = "TMM") ->
#       edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE).
#     No filterByExpr() filtering is applied before signature scoring.
#   - No R script auto-archiving is performed.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths and packages
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

main_pig_dir <- file.path(rebuild_root, "02_pig_early")
fastp_dir <- file.path(rebuild_root, "02_pig_early_fastp_sensitivity")

out_dir <- file.path(fastp_dir, "tables", "stepH2w8_current78_signature_score_stability")
log_dir <- file.path(out_dir, "logs")
fig_dir <- file.path(out_dir, "figures")
obj_dir <- file.path(out_dir, "objects")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "stepH2w8_current78_signature_score_stability_full_log.txt")
summary_log_file <- file.path(log_dir, "stepH2w8_current78_signature_score_stability_summary_to_send_me.txt")

needed_pkgs <- c("ggplot2", "edgeR")
for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(edgeR)
})

sink(full_log_file, split = TRUE)

cat("===== StepH2w8_current78 signature score stability =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")

## =========================
## 1. Helper functions
## =========================

stop_if_missing <- function(file, label) {
  if (!file.exists(file)) stop("Missing ", label, ": ", file, call. = FALSE)
}

safe_read_csv <- function(file) {
  read.csv(file, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(x, file) {
  write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", file, "\n", sep = "")
}

clean_string <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}

first_existing_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop("Could not find column for ", label, ". Checked: ", paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  hit[1]
}

read_count_matrix <- function(file) {
  df <- safe_read_csv(file)
  if (ncol(df) < 2) stop("Count matrix has too few columns: ", file)

  gene_col_candidates <- c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id")
  gene_col <- intersect(gene_col_candidates, names(df))[1]
  if (is.na(gene_col) || length(gene_col) == 0) gene_col <- names(df)[1]

  gene_id <- clean_string(df[[gene_col]])
  count_df <- df[, setdiff(names(df), gene_col), drop = FALSE]

  numeric_ok <- vapply(count_df, function(x) {
    suppressWarnings(all(!is.na(as.numeric(as.character(x)))))
  }, logical(1))
  count_df <- count_df[, numeric_ok, drop = FALSE]

  mat <- as.matrix(data.frame(lapply(count_df, function(x) as.numeric(as.character(x))),
                              check.names = FALSE))
  colnames(mat) <- names(count_df)
  rownames(mat) <- gene_id

  keep <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep, , drop = FALSE]

  if (any(duplicated(rownames(mat)))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }

  mat
}

normalize_group <- function(x) {
  z <- tolower(clean_string(x))
  out <- rep(NA_character_, length(z))
  out[grepl("control|ctrl|con|t0", z)] <- "Control"
  out[grepl("t7|7d|7 day|7_day|1w|1 week|1_week", z) & grepl("acl|aclt|transection|untreated|inj", z)] <- "ACLT_t7"
  out[grepl("t28|28d|28 day|28_day|4w|4 week|4_week", z) & grepl("acl|aclt|transection|untreated|inj", z)] <- "ACLT_t28"
  out[is.na(out) & grepl("t7|7d|7 day|7_day|1w|1 week|1_week", z)] <- "ACLT_t7"
  out[is.na(out) & grepl("t28|28d|28 day|28_day|4w|4 week|4_week", z)] <- "ACLT_t28"
  out
}

row_z <- function(mat) {
  m <- rowMeans(mat, na.rm = TRUE)
  s <- apply(mat, 1, sd, na.rm = TRUE)
  z <- sweep(mat, 1, m, "-")
  z <- sweep(z, 1, s, "/")
  z[!is.finite(z)] <- NA_real_
  z
}

compute_current78_scores <- function(count_mat, signature_df, sample_info, branch_name) {
  sample_ids <- sample_info$sample_id

  missing_samples <- setdiff(sample_ids, colnames(count_mat))
  if (length(missing_samples) > 0) {
    stop(branch_name, " count matrix is missing samples: ", paste(missing_samples, collapse = ", "))
  }

  pig_ensg <- signature_df$pig_ensg
  up_genes <- signature_df$pig_ensg[signature_df$signature_direction == "Up_in_ACLR"]
  down_genes <- signature_df$pig_ensg[signature_df$signature_direction == "Down_in_ACLR"]

  missing_sig <- setdiff(pig_ensg, rownames(count_mat))
  if (length(missing_sig) > 0) {
    return(list(
      scores = NULL,
      logcpm = NULL,
      z = NULL,
      norm_factors = NULL,
      missing_signature_genes = missing_sig
    ))
  }

  # Match the corrected current Step19 implementation exactly:
  # estimate TMM normalization factors from the full count matrix for the 18 samples,
  # then extract the 75 signature genes from edgeR TMM-normalized logCPM.
  # No filterByExpr() filtering is applied before scoring.
  count_mat <- count_mat[, sample_ids, drop = FALSE]
  lib_size <- colSums(count_mat, na.rm = TRUE)
  if (any(lib_size <= 0 | !is.finite(lib_size))) {
    stop("Invalid library size detected in ", branch_name, " count matrix.")
  }

  dge <- edgeR::DGEList(counts = count_mat)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm_all <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 1,
    normalized.lib.sizes = TRUE
  )
  logcpm_mat <- logcpm_all[pig_ensg, , drop = FALSE]
  z_mat <- row_z(logcpm_mat)

  if (anyNA(z_mat)) {
    bad <- rownames(z_mat)[apply(is.na(z_mat), 1, any)]
    stop("NA z-scores detected in ", branch_name, " for: ", paste(bad, collapse = ", "))
  }

  up_score <- colMeans(z_mat[up_genes, , drop = FALSE])
  down_score_raw <- colMeans(z_mat[down_genes, , drop = FALSE])
  down_score_reoriented <- -down_score_raw
  directional_score <- colSums(rbind(z_mat[up_genes, , drop = FALSE], -z_mat[down_genes, , drop = FALSE])) / nrow(z_mat)
  directional_score_sum_means <- up_score + down_score_reoriented
  total_score_unoriented <- colMeans(z_mat)

  scores <- data.frame(
    branch = branch_name,
    sample_id = sample_ids,
    group = sample_info$group,
    total_score_unoriented = as.numeric(total_score_unoriented[sample_ids]),
    up_score = as.numeric(up_score[sample_ids]),
    down_score_raw = as.numeric(down_score_raw[sample_ids]),
    down_score_reoriented = as.numeric(down_score_reoriented[sample_ids]),
    directional_score = as.numeric(directional_score[sample_ids]),
    directional_score_sum_means = as.numeric(directional_score_sum_means[sample_ids]),
    stringsAsFactors = FALSE
  )
  scores$group <- factor(scores$group, levels = c("Control", "ACLT_t7", "ACLT_t28"))
  scores <- scores[order(scores$group, scores$sample_id), , drop = FALSE]

  norm_factors <- data.frame(
    branch = branch_name,
    sample_id = rownames(dge$samples),
    raw_library_size = as.numeric(dge$samples$lib.size),
    TMM_norm_factor = as.numeric(dge$samples$norm.factors),
    effective_library_size = as.numeric(dge$samples$lib.size * dge$samples$norm.factors),
    stringsAsFactors = FALSE
  )

  list(
    scores = scores,
    logcpm = logcpm_mat,
    z = z_mat,
    norm_factors = norm_factors,
    missing_signature_genes = character(0)
  )
}

compare_one <- function(df, score_col, case_group) {
  control_vals <- df[df$group == "Control", score_col]
  case_vals <- df[df$group == case_group, score_col]
  wt <- suppressWarnings(wilcox.test(case_vals, control_vals, exact = FALSE, alternative = "two.sided"))
  data.frame(
    comparison = paste0(case_group, "_vs_Control"),
    score = score_col,
    n_control = length(control_vals),
    n_case = length(case_vals),
    median_control = median(control_vals, na.rm = TRUE),
    median_case = median(case_vals, na.rm = TRUE),
    mean_control = mean(control_vals, na.rm = TRUE),
    mean_case = mean(case_vals, na.rm = TRUE),
    median_difference_case_minus_control = median(case_vals, na.rm = TRUE) - median(control_vals, na.rm = TRUE),
    mean_difference_case_minus_control = mean(case_vals, na.rm = TRUE) - mean(control_vals, na.rm = TRUE),
    wilcox_p_value = as.numeric(wt$p.value),
    stringsAsFactors = FALSE
  )
}

flag_stability <- function(spearman, max_abs_delta) {
  if (is.na(spearman)) return("needs_review")
  if (spearman >= 0.98 && max_abs_delta <= 0.06) return("highly_stable")
  if (spearman >= 0.95 && max_abs_delta <= 0.10) return("stable")
  if (spearman >= 0.90) return("moderately_stable")
  "needs_review"
}

## =========================
## 2. Input files
## =========================

signature_file <- file.path(main_pig_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_pig_signature_gene_table.csv")
untrimmed_counts_file <- file.path(main_pig_dir, "tables", "step16_pig_early_gene_count_matrix.csv")
trimmed_counts_file <- file.path(fastp_dir, "tables", "stepH2w5_v2_trimmed_gene_level_counts_matrix.csv")
trimmed_sample_info_file <- file.path(fastp_dir, "tables", "stepH2w5_v2_trimmed_sample_info_for_DE.csv")
official_untrimmed_score_file <- file.path(main_pig_dir, "tables", "step19_current78_pig_signature_score", "step19_current78_pig_signature_scores_by_sample.csv")

stop_if_missing(signature_file, "current78 Step18 pig signature gene table")
stop_if_missing(untrimmed_counts_file, "untrimmed Step16 pig early count matrix")
stop_if_missing(trimmed_counts_file, "trimmed StepH2w5 count matrix")
stop_if_missing(trimmed_sample_info_file, "trimmed StepH2w5 sample info")

cat("Signature file: ", signature_file, "\n", sep = "")
cat("Untrimmed count matrix: ", untrimmed_counts_file, "\n", sep = "")
cat("Trimmed count matrix: ", trimmed_counts_file, "\n", sep = "")
cat("Trimmed sample info: ", trimmed_sample_info_file, "\n", sep = "")
cat("Official untrimmed Step19 score file: ", ifelse(file.exists(official_untrimmed_score_file), official_untrimmed_score_file, "not found; audit skipped"), "\n\n", sep = "")

## =========================
## 3. Load inputs
## =========================

sig_raw <- safe_read_csv(signature_file)
untrimmed_counts <- read_count_matrix(untrimmed_counts_file)
trimmed_counts <- read_count_matrix(trimmed_counts_file)
sample_info_raw <- safe_read_csv(trimmed_sample_info_file)

pig_ensg_col <- first_existing_col(sig_raw, c("pig_ensg", "pig_gene_id", "pig_gene", "gene_id"), "signature pig Ensembl ID")
direction_col <- first_existing_col(sig_raw, c("signature_direction", "mouse_direction", "direction"), "signature direction")

signature_df <- data.frame(
  pig_ensg = clean_string(sig_raw[[pig_ensg_col]]),
  signature_direction = clean_string(sig_raw[[direction_col]]),
  stringsAsFactors = FALSE
)
signature_df <- signature_df[!is.na(signature_df$pig_ensg) & !is.na(signature_df$signature_direction), , drop = FALSE]
signature_df <- signature_df[signature_df$signature_direction %in% c("Up_in_ACLR", "Down_in_ACLR"), , drop = FALSE]
signature_df <- signature_df[!duplicated(signature_df$pig_ensg), , drop = FALSE]

sample_id_col <- first_existing_col(sample_info_raw, c("sample_id", "sample", "SampleID", "sample_name"), "sample info sample_id")
group_col <- first_existing_col(sample_info_raw, c("group", "condition", "time_group", "treatment_group", "group_label"), "sample info group")

sample_info <- data.frame(
  sample_id = clean_string(sample_info_raw[[sample_id_col]]),
  group_raw = clean_string(sample_info_raw[[group_col]]),
  stringsAsFactors = FALSE
)
sample_info$group <- normalize_group(sample_info$group_raw)
if (any(is.na(sample_info$group))) {
  write_csv(sample_info, file.path(out_dir, "stepH2w8_current78_sample_group_normalization_failed.csv"))
  stop("Some samples could not be normalized to Control/ACLT_t7/ACLT_t28. See diagnostic file.")
}
sample_info <- sample_info[!duplicated(sample_info$sample_id), , drop = FALSE]
sample_info$group <- factor(sample_info$group, levels = c("Control", "ACLT_t7", "ACLT_t28"))
sample_info <- sample_info[order(sample_info$group, sample_info$sample_id), , drop = FALSE]

cat("Signature gene count: ", nrow(signature_df), "\n", sep = "")
cat("Up signature genes: ", sum(signature_df$signature_direction == "Up_in_ACLR"), "\n", sep = "")
cat("Down signature genes: ", sum(signature_df$signature_direction == "Down_in_ACLR"), "\n", sep = "")
cat("Sample counts by group:\n")
print(table(sample_info$group))

if (nrow(signature_df) != 75 ||
    sum(signature_df$signature_direction == "Up_in_ACLR") != 65 ||
    sum(signature_df$signature_direction == "Down_in_ACLR") != 10) {
  stop("Current78 signature counts do not match expected 75/65/10. Please check signature input file.")
}

## =========================
## 4. Compute scores in both branches
## =========================

untrimmed_res <- compute_current78_scores(untrimmed_counts, signature_df, sample_info, "untrimmed")
trimmed_res <- compute_current78_scores(trimmed_counts, signature_df, sample_info, "trimmed")

if (length(untrimmed_res$missing_signature_genes) > 0) {
  write_csv(data.frame(pig_ensg = untrimmed_res$missing_signature_genes), file.path(out_dir, "stepH2w8_current78_missing_signature_genes_untrimmed.csv"))
  stop("Some current78 signature genes are missing in untrimmed count matrix.")
}
if (length(trimmed_res$missing_signature_genes) > 0) {
  write_csv(data.frame(pig_ensg = trimmed_res$missing_signature_genes), file.path(out_dir, "stepH2w8_current78_missing_signature_genes_trimmed.csv"))
  stop("Some current78 signature genes are missing in trimmed count matrix.")
}

score_untrimmed <- untrimmed_res$scores
score_trimmed <- trimmed_res$scores
score_all <- rbind(score_untrimmed, score_trimmed)

score_cols <- c("directional_score", "directional_score_sum_means", "up_score", "down_score_reoriented", "down_score_raw", "total_score_unoriented")

score_long <- do.call(rbind, lapply(score_cols, function(sc) {
  data.frame(
    branch = score_all$branch,
    sample_id = score_all$sample_id,
    group = as.character(score_all$group),
    score_name = sc,
    score_value = score_all[[sc]],
    stringsAsFactors = FALSE
  )
}))

score_wide <- do.call(rbind, lapply(score_cols, function(sc) {
  u <- score_untrimmed[, c("sample_id", "group", sc), drop = FALSE]
  t <- score_trimmed[, c("sample_id", sc), drop = FALSE]
  colnames(u)[colnames(u) == sc] <- "untrimmed_score"
  colnames(t)[colnames(t) == sc] <- "trimmed_score"
  m <- merge(u, t, by = "sample_id", all = FALSE)
  m$score_name <- sc
  m$delta_trimmed_minus_untrimmed <- m$trimmed_score - m$untrimmed_score
  m[, c("sample_id", "group", "score_name", "untrimmed_score", "trimmed_score", "delta_trimmed_minus_untrimmed")]
}))

score_stability <- do.call(rbind, lapply(score_cols, function(sc) {
  x <- score_wide[score_wide$score_name == sc, , drop = FALSE]
  sp <- suppressWarnings(cor(x$untrimmed_score, x$trimmed_score, method = "spearman", use = "complete.obs"))
  pe <- suppressWarnings(cor(x$untrimmed_score, x$trimmed_score, method = "pearson", use = "complete.obs"))
  maxad <- max(abs(x$delta_trimmed_minus_untrimmed), na.rm = TRUE)
  data.frame(
    score_name = sc,
    n_samples = nrow(x),
    spearman = sp,
    pearson = pe,
    mean_abs_delta = mean(abs(x$delta_trimmed_minus_untrimmed), na.rm = TRUE),
    median_abs_delta = median(abs(x$delta_trimmed_minus_untrimmed), na.rm = TRUE),
    max_abs_delta = maxad,
    stability_flag = flag_stability(sp, maxad),
    stringsAsFactors = FALSE
  )
}))

## =========================
## 5. Group comparisons by branch
## =========================

group_comparison <- do.call(rbind, lapply(c("untrimmed", "trimmed"), function(branch_name) {
  dfb <- if (branch_name == "untrimmed") score_untrimmed else score_trimmed
  do.call(rbind, lapply(score_cols, function(sc) {
    rbind(compare_one(dfb, sc, "ACLT_t7"), compare_one(dfb, sc, "ACLT_t28")) |>
      transform(branch = branch_name)
  }))
}))
group_comparison <- group_comparison[, c("branch", "comparison", "score",
                                         "n_control", "n_case",
                                         "median_control", "median_case",
                                         "mean_control", "mean_case",
                                         "median_difference_case_minus_control",
                                         "mean_difference_case_minus_control",
                                         "wilcox_p_value")]
group_comparison$BH_FDR_within_step <- p.adjust(group_comparison$wilcox_p_value, method = "BH")

## =========================
## 6. Audit official Step19 untrimmed scores if present
## =========================

official_audit <- data.frame(
  metric = c("official_step19_score_file_found", "max_abs_delta_directional_score_recomputed_vs_official", "max_abs_delta_up_score_recomputed_vs_official", "max_abs_delta_down_score_reoriented_recomputed_vs_official"),
  value = c(file.exists(official_untrimmed_score_file), NA, NA, NA),
  stringsAsFactors = FALSE
)

if (file.exists(official_untrimmed_score_file)) {
  official_df <- safe_read_csv(official_untrimmed_score_file)
  if (all(c("sample_id", "directional_score", "up_score", "down_score_reoriented") %in% colnames(official_df))) {
    m <- merge(score_untrimmed, official_df[, c("sample_id", "directional_score", "up_score", "down_score_reoriented")],
               by = "sample_id", suffixes = c("_recomputed", "_official"))
    official_audit$value[official_audit$metric == "max_abs_delta_directional_score_recomputed_vs_official"] <- max(abs(m$directional_score_recomputed - m$directional_score_official), na.rm = TRUE)
    official_audit$value[official_audit$metric == "max_abs_delta_up_score_recomputed_vs_official"] <- max(abs(m$up_score_recomputed - m$up_score_official), na.rm = TRUE)
    official_audit$value[official_audit$metric == "max_abs_delta_down_score_reoriented_recomputed_vs_official"] <- max(abs(m$down_score_reoriented_recomputed - m$down_score_reoriented_official), na.rm = TRUE)
  }
}

## =========================
## 7. Simple stability plot
## =========================

plot_df <- score_wide[score_wide$score_name %in% c("directional_score", "up_score", "down_score_reoriented"), , drop = FALSE]
plot_df$score_label <- factor(plot_df$score_name,
                              levels = c("directional_score", "up_score", "down_score_reoriented"),
                              labels = c("Directional score", "Up-score", "Down-score reoriented"))
p <- ggplot(plot_df, aes(x = untrimmed_score, y = trimmed_score)) +
  geom_point(aes(shape = group), size = 2.6, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ score_label, scales = "free", nrow = 1) +
  labs(
    title = "Current78 pig early signature score stability after fastp trimming",
    x = "Untrimmed score",
    y = "Trimmed score"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )

plot_png <- file.path(fig_dir, "stepH2w8_current78_signature_score_trimmed_vs_untrimmed_scatter.png")
plot_pdf <- file.path(fig_dir, "stepH2w8_current78_signature_score_trimmed_vs_untrimmed_scatter.pdf")
png(plot_png, width = 3000, height = 1000, res = 250)
print(p)
dev.off()
pdf(plot_pdf, width = 12, height = 4)
print(p)
dev.off()

## =========================
## 8. Save outputs
## =========================

signature_used_file <- file.path(out_dir, "stepH2w8_current78_signature_genes_used.csv")
scores_long_file <- file.path(out_dir, "stepH2w8_current78_signature_scores_long.csv")
scores_wide_file <- file.path(out_dir, "stepH2w8_current78_signature_score_trimmed_vs_untrimmed_wide.csv")
score_stability_file <- file.path(out_dir, "stepH2w8_current78_signature_score_stability_metrics.csv")
group_comparison_file <- file.path(out_dir, "stepH2w8_current78_signature_score_group_comparison_by_branch.csv")
official_audit_file <- file.path(out_dir, "stepH2w8_current78_official_step19_untrimmed_score_audit.csv")
untrimmed_logcpm_file <- file.path(out_dir, "stepH2w8_current78_untrimmed_signature_logCPM_matrix.csv")
trimmed_logcpm_file <- file.path(out_dir, "stepH2w8_current78_trimmed_signature_logCPM_matrix.csv")
untrimmed_z_file <- file.path(out_dir, "stepH2w8_current78_untrimmed_signature_zscore_matrix.csv")
trimmed_z_file <- file.path(out_dir, "stepH2w8_current78_trimmed_signature_zscore_matrix.csv")

write_csv(signature_df, signature_used_file)
write_csv(score_long, scores_long_file)
write_csv(score_wide, scores_wide_file)
write_csv(score_stability, score_stability_file)
write_csv(group_comparison, group_comparison_file)
write_csv(official_audit, official_audit_file)
write.csv(untrimmed_res$logcpm, untrimmed_logcpm_file, quote = TRUE)
write.csv(trimmed_res$logcpm, trimmed_logcpm_file, quote = TRUE)
write.csv(untrimmed_res$z, untrimmed_z_file, quote = TRUE)
write.csv(trimmed_res$z, trimmed_z_file, quote = TRUE)
write_csv(untrimmed_res$norm_factors, file.path(out_dir, "stepH2w8_current78_untrimmed_TMM_normalization_factors.csv"))
write_csv(trimmed_res$norm_factors, file.path(out_dir, "stepH2w8_current78_trimmed_TMM_normalization_factors.csv"))

primary_scores <- c("directional_score", "up_score", "down_score_reoriented")
primary_stability <- score_stability[score_stability$score_name %in% primary_scores, , drop = FALSE]
overall_stability_flag <- if (all(primary_stability$stability_flag == "highly_stable")) {
  "highly_stable"
} else if (all(primary_stability$stability_flag %in% c("highly_stable", "stable"))) {
  "stable"
} else if (all(primary_stability$stability_flag %in% c("highly_stable", "stable", "moderately_stable"))) {
  "moderately_stable"
} else {
  "needs_review"
}

summary_df <- data.frame(
  metric = c(
    "signature_file",
    "untrimmed_counts_file",
    "trimmed_counts_file",
    "trimmed_sample_info_file",
    "n_signature_genes_used",
    "n_up_signature_genes_used",
    "n_down_signature_genes_used",
    "n_detected_untrimmed",
    "n_detected_trimmed",
    "n_samples_compared",
    "primary_scores",
    "min_spearman_primary_scores",
    "max_abs_delta_primary_scores",
    "overall_stability_flag",
    "plot_png",
    "plot_pdf",
    "scores_wide_file",
    "score_stability_file",
    "group_comparison_file",
    "official_untrimmed_score_audit_file",
    "output_dir"
  ),
  value = c(
    signature_file,
    untrimmed_counts_file,
    trimmed_counts_file,
    trimmed_sample_info_file,
    nrow(signature_df),
    sum(signature_df$signature_direction == "Up_in_ACLR"),
    sum(signature_df$signature_direction == "Down_in_ACLR"),
    nrow(signature_df) - length(untrimmed_res$missing_signature_genes),
    nrow(signature_df) - length(trimmed_res$missing_signature_genes),
    nrow(sample_info),
    paste(primary_scores, collapse = "; "),
    round(min(primary_stability$spearman, na.rm = TRUE), 6),
    round(max(primary_stability$max_abs_delta, na.rm = TRUE), 6),
    overall_stability_flag,
    plot_png,
    plot_pdf,
    scores_wide_file,
    score_stability_file,
    group_comparison_file,
    official_audit_file,
    out_dir
  ),
  stringsAsFactors = FALSE
)

summary_file <- file.path(out_dir, "stepH2w8_current78_signature_score_stability_summary.csv")
write_csv(summary_df, summary_file)

saveRDS(
  list(
    signature_df = signature_df,
    sample_info = sample_info,
    score_untrimmed = score_untrimmed,
    score_trimmed = score_trimmed,
    score_wide = score_wide,
    score_stability = score_stability,
    group_comparison = group_comparison,
    official_audit = official_audit,
    summary_df = summary_df
  ),
  file.path(obj_dir, "stepH2w8_current78_signature_score_stability_workspace.rds")
)

cat("\n===== StepH2w8_current78 summary =====\n")
print(summary_df, row.names = FALSE)

cat("\n===== Score stability metrics =====\n")
print(score_stability, row.names = FALSE)

cat("\n===== Group comparison for manuscript-relevant scores =====\n")
print(group_comparison[group_comparison$score %in% primary_scores, ], row.names = FALSE)

cat("\n===== Official Step19 untrimmed score audit =====\n")
print(official_audit, row.names = FALSE)

version_df <- data.frame(
  item = c("R_version", "edgeR_version", "ggplot2_version", "score_formula", "logCPM_formula", "filtering_before_score", "signature_source", "branch_role"),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("ggplot2")),
    "directional_score = (sum z_up + sum -z_down)/(n_up+n_down)",
    "edgeR TMM-normalized logCPM: DGEList(full count matrix) -> calcNormFactors(method='TMM') -> cpm(log=TRUE, prior.count=1, normalized.lib.sizes=TRUE), matching current Step19",
    "No filterByExpr() filtering before signature scoring; all Step18-detected signature genes are retained",
    "fixed current78 Step18 pig signature gene table; no automatic signature search",
    "fastp trimming sensitivity analysis, not replacement for main untrimmed analysis"
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "stepH2w8_current78_versions_and_method_records.csv")
write_csv(version_df, version_file)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "stepH2w8_current78_sessionInfo.txt"))

## concise summary log
summary_con <- file(summary_log_file, open = "wt")
writeLines("===== StepH2w8_current78 signature score stability SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(summary_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Score stability metrics:", summary_con)
writeLines(capture.output(print(score_stability, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Group comparison for manuscript-relevant scores:", summary_con)
writeLines(capture.output(print(group_comparison[group_comparison$score %in% primary_scores, ], row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Official Step19 untrimmed score audit:", summary_con)
writeLines(capture.output(print(official_audit, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("Please verify n_signature_genes_used = 75, n_up_signature_genes_used = 65, n_down_signature_genes_used = 10.", summary_con)
writeLines("Please verify n_detected_trimmed = 75. If not, downstream trimmed score/core validation must use the detected subset.", summary_con)
writeLines("The main stability metric is the Spearman correlation between untrimmed and trimmed scores for directional_score, up_score, and down_score_reoriented.", summary_con)
close(summary_con)

sink()

cat("\nStepH2w8_current78 signature score stability completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
