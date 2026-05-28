# Chronic Step2 current78 TMM-aligned: 75-gene signature score in chronic pig synovium
# Purpose:
#   Recompute chronic pig signature scores using the fixed current78-derived pig early
#   75-gene signature in the GSE228848 chronic main comparison:
#     Control_52W vs ACLT_alone_52W.
#
# Key methodological update:
#   This version aligns the chronic signature-score input transformation with the
#   pig early signature-score method:
#     edgeR::DGEList()
#     edgeR::calcNormFactors(method = "TMM")
#     edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
#   No filterByExpr() is applied before scoring, to avoid losing predefined signature genes.
#   Salmon/tximport estimated counts are used as numeric values and are NOT rounded.
#
# Important interpretation:
#   - Chronic pig analysis is an extension validation, not a discovery/redefinition step.
#   - This step does NOT redefine the signature.
#   - The fixed signature is the pig early current78-derived 75-gene signature:
#       65 Up_in_ACLR genes + 10 Down_in_ACLR genes.
#   - The chronic main comparison is Control_52W vs ACLT_alone_52W, 12 vs 12.
#
# Outputs:
#   This script writes to a separate TMM-aligned output folder so the old non-TMM
#   Step2 results are not overwritten before review.

options(stringsAsFactors = FALSE)

status <- "FAILED"
method_version <- "2026-05-16_chronic_step2_current78_TMM_aligned_with_pig_early_Step19"

## =========================
## 0. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_early_dir <- file.path(rebuild_root, "02_pig_early")
pig_chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

early_step18_dir <- file.path(
  pig_early_dir,
  "tables",
  "step18_current78_pig_early_signature_remap"
)

chronic_tables_dir <- file.path(pig_chronic_dir, "tables")

out_dir <- file.path(
  chronic_tables_dir,
  "chronic_step2_current78_signature_score_TMM_aligned"
)
log_dir <- file.path(out_dir, "logs")
obj_dir <- file.path(out_dir, "objects")

chronic_figures_dir <- file.path(
  pig_chronic_dir,
  "figures",
  "Figure5A_current78_chronic_signature_score_TMM_aligned"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chronic_figures_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "chronic_step2_current78_TMM_aligned_signature_score_full_log.txt")
summary_log_file <- file.path(log_dir, "chronic_step2_current78_TMM_aligned_signature_score_summary_to_send_me.txt")

zz <- file(full_log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")

close_log <- function() {
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}

on.exit({
  cat("\n============================================================\n")
  cat("Chronic Step2 current78 TMM-aligned signature score finished with status:", status, "\n")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Summary log saved to:\n", summary_log_file, "\n", sep = "")
  cat("============================================================\n")
  close_log()
}, add = TRUE)

cat("===== CHRONIC STEP2 CURRENT78 TMM-ALIGNED SIGNATURE SCORE =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Method version: ", method_version, "\n", sep = "")
cat("Interpretation: chronic pig extension validation using fixed pig early current78 75-gene signature.\n")
cat("Normalization: edgeR TMM-normalized logCPM; no filterByExpr before scoring.\n")
cat("Primary score: directional_score = (sum z_up + sum -z_down)/(n_up+n_down).\n\n")

## =========================
## 1. Packages
## =========================

for (pkg in c("edgeR", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package not installed: ", pkg, call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
})

cat("Required packages loaded.\n")
cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n\n")

## =========================
## 2. Helper functions
## =========================

stop_if_missing <- function(file, label) {
  if (!file.exists(file)) stop("Missing ", label, ": ", file, call. = FALSE)
}

find_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_read_csv <- function(file) {
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(x, file) {
  write.csv(x, file, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
  cat("Saved: ", file, "\n", sep = "")
}

write_matrix_csv <- function(mat, file) {
  write.csv(mat, file, quote = TRUE, fileEncoding = "UTF-8")
  cat("Saved: ", file, "\n", sep = "")
}

clean_string <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}

find_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Could not find ", label, ". Checked: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

normalize_chronic_group <- function(x) {
  z <- tolower(clean_string(x))
  out <- rep(NA_character_, length(z))

  out[grepl("control", z) & grepl("52", z)] <- "Control_52W"
  out[grepl("ctrl", z) & grepl("52", z)] <- "Control_52W"

  out[grepl("aclt", z) & grepl("alone", z) & grepl("52", z)] <- "ACLT_alone_52W"
  out[grepl("acl transection", z) & grepl("alone", z) & grepl("52", z)] <- "ACLT_alone_52W"
  out[grepl("transection alone", z) & grepl("52", z)] <- "ACLT_alone_52W"

  exact <- clean_string(x)
  out[exact %in% c("Control_52W", "ACLT_alone_52W")] <- exact[exact %in% c("Control_52W", "ACLT_alone_52W")]

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

file_audit <- function(path, label) {
  info <- if (file.exists(path)) file.info(path) else data.frame(size = NA, mtime = NA)
  md5 <- NA_character_
  if (file.exists(path) && requireNamespace("tools", quietly = TRUE)) {
    md5 <- as.character(tools::md5sum(path))
  }
  data.frame(
    label = label,
    file = path,
    exists = file.exists(path),
    size_bytes = if (file.exists(path)) info$size else NA_real_,
    mtime = if (file.exists(path)) as.character(info$mtime) else NA_character_,
    md5 = md5,
    stringsAsFactors = FALSE
  )
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "p = NA",
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))
  )
}

compare_one <- function(df, score_col) {
  control_vals <- df[df$group == "Control_52W", score_col]
  case_vals <- df[df$group == "ACLT_alone_52W", score_col]
  wt <- suppressWarnings(
    wilcox.test(case_vals, control_vals, exact = FALSE, alternative = "two.sided")
  )
  data.frame(
    comparison = "ACLT_alone_52W_vs_Control_52W",
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

## =========================
## 3. Input files
## =========================

signature_file <- file.path(
  early_step18_dir,
  "step18_current78_pig_signature_gene_table.csv"
)

chronic_main_manifest_file <- find_existing_file(c(
  file.path(chronic_tables_dir, "step25v3_pig_chronic_main_comparison_manifest.csv"),
  file.path(chronic_tables_dir, "step25_pig_chronic_main_comparison_manifest.csv"),
  file.path(chronic_tables_dir, "pig_chronic_main_comparison_manifest.csv")
))

chronic_main_counts_file <- find_existing_file(c(
  file.path(chronic_tables_dir, "step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv"),
  file.path(chronic_tables_dir, "step25_pig_chronic_main_comparison_gene_level_counts_matrix.csv"),
  file.path(chronic_tables_dir, "pig_chronic_main_comparison_gene_level_counts_matrix.csv")
))

step1_summary_file <- file.path(
  chronic_tables_dir,
  "chronic_step1_current78_signature_core_detectability_audit",
  "chronic_step1_current78_signature_core_detectability_summary.csv"
)

stop_if_missing(signature_file, "current78 pig early 75-gene signature table")
if (is.na(chronic_main_manifest_file)) stop("Could not find chronic main comparison manifest under: ", chronic_tables_dir, call. = FALSE)
if (is.na(chronic_main_counts_file)) stop("Could not find chronic main comparison counts matrix under: ", chronic_tables_dir, call. = FALSE)

cat("Signature file: ", signature_file, "\n", sep = "")
cat("Chronic main manifest: ", chronic_main_manifest_file, "\n", sep = "")
cat("Chronic main count matrix: ", chronic_main_counts_file, "\n", sep = "")
cat("Step1 detectability summary: ", step1_summary_file, " exists=", file.exists(step1_summary_file), "\n\n", sep = "")

## =========================
## 4. Load fixed current78 signature
## =========================

sig_raw <- safe_read_csv(signature_file)

pig_ensg_col <- find_col(sig_raw, c("pig_ensg", "pig_gene_id", "pig_gene", "gene_id"), TRUE, "signature pig Ensembl ID")
pig_symbol_col <- find_col(sig_raw, c("pig_symbol", "pig_gene_symbol", "gene_symbol", "symbol", "gene_name"), FALSE, "signature pig symbol")
human_col <- find_col(sig_raw, c("human_gene", "ortholog_name", "Gene symbol", "gene_symbol_human"), FALSE, "signature human gene")
mouse_col <- find_col(sig_raw, c("mouse_symbol", "input", "mouse_gene", "mouse_gene_symbol"), FALSE, "signature mouse symbol")
direction_col <- find_col(sig_raw, c("signature_direction", "mouse_direction", "direction"), TRUE, "signature direction")

sig_df <- data.frame(
  pig_ensg = clean_string(sig_raw[[pig_ensg_col]]),
  pig_symbol = if (!is.na(pig_symbol_col)) clean_string(sig_raw[[pig_symbol_col]]) else NA_character_,
  human_gene = if (!is.na(human_col)) clean_string(sig_raw[[human_col]]) else NA_character_,
  mouse_symbol = if (!is.na(mouse_col)) clean_string(sig_raw[[mouse_col]]) else NA_character_,
  signature_direction = clean_string(sig_raw[[direction_col]]),
  stringsAsFactors = FALSE
)

sig_df <- sig_df[
  !is.na(sig_df$pig_ensg) &
    sig_df$pig_ensg != "" &
    sig_df$signature_direction %in% c("Up_in_ACLR", "Down_in_ACLR"),
  ,
  drop = FALSE
]
sig_df <- sig_df[!duplicated(sig_df$pig_ensg), , drop = FALSE]

up_genes <- sig_df$pig_ensg[sig_df$signature_direction == "Up_in_ACLR"]
down_genes <- sig_df$pig_ensg[sig_df$signature_direction == "Down_in_ACLR"]

if (nrow(sig_df) != 75 || length(up_genes) != 65 || length(down_genes) != 10) {
  stop(
    "Current signature table does not match expected current78 75/65/10 counts. Detected: ",
    nrow(sig_df), " total, ", length(up_genes), " up, ", length(down_genes), " down.",
    call. = FALSE
  )
}

cat("Fixed signature genes: ", nrow(sig_df), "\n", sep = "")
cat("Up genes: ", length(up_genes), "\n", sep = "")
cat("Down genes: ", length(down_genes), "\n\n", sep = "")

## =========================
## 5. Load chronic counts and metadata
## =========================

manifest <- safe_read_csv(chronic_main_manifest_file)
count_df <- safe_read_csv(chronic_main_counts_file)

sample_col <- find_col(manifest, c("sample_id", "sample", "SampleID", "geo_accession", "GSM", "run", "Run"), TRUE, "chronic manifest sample ID")
group_col <- find_col(manifest, c("core_group", "group", "treatment", "condition", "title", "sample_title"), TRUE, "chronic manifest group")

sample_info <- data.frame(
  sample_id = clean_string(manifest[[sample_col]]),
  group_raw = clean_string(manifest[[group_col]]),
  stringsAsFactors = FALSE
)

sample_info$group <- normalize_chronic_group(sample_info$group_raw)

if ("core_group" %in% colnames(manifest)) {
  cg <- clean_string(manifest$core_group)
  ok <- cg %in% c("Control_52W", "ACLT_alone_52W")
  sample_info$group[ok] <- cg[ok]
}

sample_info <- sample_info[
  !is.na(sample_info$sample_id) &
    sample_info$sample_id != "" &
    sample_info$group %in% c("Control_52W", "ACLT_alone_52W"),
  ,
  drop = FALSE
]
sample_info <- sample_info[!duplicated(sample_info$sample_id), , drop = FALSE]
sample_info$group <- factor(sample_info$group, levels = c("Control_52W", "ACLT_alone_52W"))
sample_info <- sample_info[order(sample_info$group, sample_info$sample_id), , drop = FALSE]

if (sum(sample_info$group == "Control_52W") != 12 ||
    sum(sample_info$group == "ACLT_alone_52W") != 12) {
  stop(
    "Expected 12 Control_52W and 12 ACLT_alone_52W samples. Observed: ",
    sum(sample_info$group == "Control_52W"), " Control_52W and ",
    sum(sample_info$group == "ACLT_alone_52W"), " ACLT_alone_52W.",
    call. = FALSE
  )
}

gene_col <- find_col(count_df, c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id"), TRUE, "chronic count matrix gene ID")
count_df[[gene_col]] <- clean_string(count_df[[gene_col]])

sample_cols <- setdiff(colnames(count_df), gene_col)
numeric_sample_cols <- sample_cols[vapply(count_df[sample_cols], function(x) {
  suppressWarnings(all(is.na(x) | !is.na(as.numeric(as.character(x)))))
}, logical(1))]

missing_samples <- setdiff(sample_info$sample_id, numeric_sample_cols)
if (length(missing_samples) > 0) {
  write_csv(data.frame(missing_sample_id = missing_samples), file.path(out_dir, "chronic_step2_TMM_aligned_missing_samples_in_count_matrix.csv"))
  stop("Some chronic samples are missing from count matrix. See diagnostic file.", call. = FALSE)
}

count_mat <- as.matrix(data.frame(
  lapply(count_df[, numeric_sample_cols, drop = FALSE], function(x) as.numeric(as.character(x))),
  check.names = FALSE
))
rownames(count_mat) <- count_df[[gene_col]]
valid_rows <- !is.na(rownames(count_mat)) & rownames(count_mat) != ""
count_mat <- count_mat[valid_rows, , drop = FALSE]

if (any(duplicated(rownames(count_mat)))) {
  count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
}

count_mat <- count_mat[, sample_info$sample_id, drop = FALSE]

if (anyNA(count_mat)) stop("NA values detected after numeric conversion of chronic count matrix.", call. = FALSE)
if (any(count_mat < 0, na.rm = TRUE)) stop("Negative values detected in chronic count matrix.", call. = FALSE)

cat("Chronic main count matrix dimensions used: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples\n", sep = "")
cat("Sample group counts:\n")
print(table(sample_info$group))
cat("\n")

## =========================
## 6. TMM-normalized logCPM without filterByExpr
## =========================

dge <- edgeR::DGEList(counts = count_mat, group = sample_info$group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm_all <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 1,
  normalized.lib.sizes = TRUE
)

normalization_summary <- data.frame(
  metric = c(
    "input_count_matrix_genes",
    "input_count_matrix_samples",
    "filterByExpr_applied_before_scoring",
    "edgeR_DGEList_used",
    "TMM_normalization",
    "logCPM_function",
    "prior_count",
    "normalized_lib_sizes",
    "estimated_counts_rounded"
  ),
  value = c(
    nrow(count_mat),
    ncol(count_mat),
    "FALSE",
    "TRUE",
    "edgeR::calcNormFactors(method = 'TMM')",
    "edgeR::cpm(log = TRUE)",
    "1",
    "TRUE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)

lib_summary <- data.frame(
  sample_id = colnames(count_mat),
  group = as.character(sample_info$group),
  raw_estimated_count_sum = colSums(count_mat),
  TMM_norm_factor = dge$samples$norm.factors,
  effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
  detected_genes_count_gt0 = colSums(count_mat > 0),
  stringsAsFactors = FALSE
)

cat("Normalization summary:\n")
print(normalization_summary)
cat("\nTMM library summary:\n")
print(lib_summary)
cat("\n")

## =========================
## 7. Match signature genes and compute scores
## =========================

missing_sig <- sig_df[!(sig_df$pig_ensg %in% rownames(logcpm_all)), , drop = FALSE]
if (nrow(missing_sig) > 0) {
  write_csv(missing_sig, file.path(out_dir, "chronic_step2_TMM_aligned_signature_genes_missing_in_chronic_count_matrix.csv"))
  stop("Some fixed signature genes are missing in chronic count matrix. See diagnostic file.", call. = FALSE)
}

sig_logcpm <- logcpm_all[sig_df$pig_ensg, , drop = FALSE]
z_mat <- row_z(sig_logcpm)

if (anyNA(z_mat)) {
  na_rows <- rownames(z_mat)[apply(is.na(z_mat), 1, any)]
  write_csv(data.frame(pig_ensg = na_rows), file.path(out_dir, "chronic_step2_TMM_aligned_signature_genes_with_NA_zscore.csv"))
  stop("NA z-scores detected for some signature genes, usually due to zero variance. See output file.", call. = FALSE)
}

up_idx <- sig_df$signature_direction == "Up_in_ACLR"
down_idx <- sig_df$signature_direction == "Down_in_ACLR"

up_score <- colMeans(z_mat[up_idx, , drop = FALSE])
down_score_raw <- colMeans(z_mat[down_idx, , drop = FALSE])
down_score_reoriented <- -down_score_raw
directional_score <- colSums(rbind(z_mat[up_idx, , drop = FALSE], -z_mat[down_idx, , drop = FALSE])) / nrow(z_mat)
directional_score_sum_means <- up_score + down_score_reoriented
total_score_unoriented <- colMeans(z_mat)

scores_df <- data.frame(
  sample_id = colnames(z_mat),
  group = as.character(sample_info$group),
  group_raw = sample_info$group_raw,
  total_score_unoriented = as.numeric(total_score_unoriented),
  up_score = as.numeric(up_score),
  down_score_raw = as.numeric(down_score_raw),
  down_score_reoriented = as.numeric(down_score_reoriented),
  directional_score = as.numeric(directional_score),
  directional_score_sum_means = as.numeric(directional_score_sum_means),
  stringsAsFactors = FALSE
)
scores_df$group <- factor(scores_df$group, levels = c("Control_52W", "ACLT_alone_52W"))
scores_df <- scores_df[order(scores_df$group, scores_df$sample_id), , drop = FALSE]

score_cols <- c(
  "directional_score",
  "directional_score_sum_means",
  "up_score",
  "down_score_reoriented",
  "down_score_raw",
  "total_score_unoriented"
)

stats_df <- do.call(rbind, lapply(score_cols, function(sc) compare_one(scores_df, sc)))
stats_df$BH_FDR_across_all_score_comparisons <- p.adjust(stats_df$wilcox_p_value, method = "BH")

scores_long <- do.call(rbind, lapply(score_cols, function(sc) {
  data.frame(
    sample_id = scores_df$sample_id,
    group = as.character(scores_df$group),
    score_name = sc,
    score_value = scores_df[[sc]],
    stringsAsFactors = FALSE
  )
}))

## =========================
## 8. Figure5A-style score plot
## =========================

score_plot_levels <- c("directional_score", "up_score", "down_score_reoriented")
score_plot_labels <- c(
  directional_score = "Directional score",
  up_score = "Up-score",
  down_score_reoriented = "Down-score reoriented"
)

plot_long <- scores_long[scores_long$score_name %in% score_plot_levels, , drop = FALSE]
plot_long$score_label <- factor(score_plot_labels[plot_long$score_name], levels = unname(score_plot_labels[score_plot_levels]))
plot_long$group <- factor(plot_long$group, levels = c("Control_52W", "ACLT_alone_52W"))

stats_plot <- stats_df[stats_df$score %in% score_plot_levels, , drop = FALSE]
stats_plot$score_label <- factor(score_plot_labels[stats_plot$score], levels = unname(score_plot_labels[score_plot_levels]))
stats_plot$p_label <- fmt_p(stats_plot$wilcox_p_value)

y_ranges <- aggregate(score_value ~ score_label, plot_long, function(x) diff(range(x, na.rm = TRUE)))
y_max <- aggregate(score_value ~ score_label, plot_long, max, na.rm = TRUE)
y_pos <- merge(y_max, y_ranges, by = "score_label")
colnames(y_pos) <- c("score_label", "y_max", "y_range")
y_pos$y_range[!is.finite(y_pos$y_range) | y_pos$y_range == 0] <- 1

stats_plot <- merge(stats_plot, y_pos, by = "score_label", all.x = TRUE)
stats_plot$xmid <- 1.5
stats_plot$xend <- 2
stats_plot$y_segment <- stats_plot$y_max + 0.07 * stats_plot$y_range
stats_plot$y_label <- stats_plot$y_max + 0.15 * stats_plot$y_range

point_fill <- c(Control_52W = "white", ACLT_alone_52W = "#D95F02")
point_edge <- c(Control_52W = "#7A7A7A", ACLT_alone_52W = "#7F3300")

p_combined <- ggplot(plot_long, aes(x = group, y = score_value)) +
  geom_point(
    aes(fill = group, color = group),
    shape = 21,
    stroke = 0.55,
    size = 3.1,
    position = position_jitter(width = 0.08, height = 0, seed = 1),
    alpha = 0.95
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.44,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_segment(
    data = stats_plot,
    aes(x = 1, xend = xend, y = y_segment, yend = y_segment),
    inherit.aes = FALSE,
    color = "#666666",
    linewidth = 0.35
  ) +
  geom_text(
    data = stats_plot,
    aes(x = xmid, y = y_label, label = p_label),
    inherit.aes = FALSE,
    size = 4.0
  ) +
  facet_grid(rows = vars(score_label), switch = "y", scales = "free_y") +
  scale_fill_manual(values = point_fill) +
  scale_color_manual(values = point_edge) +
  labs(
    title = "Chronic pig validation of current78-derived 75-gene signature",
    subtitle = "TMM-normalized logCPM; no filtering before signature scoring",
    x = NULL,
    y = "Mean row-z signature score"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    panel.grid.major = element_line(color = "#E6E6E6"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = "black"),
    strip.text = element_text(size = 13),
    strip.text.y.left = element_text(angle = -90, size = 13),
    axis.text.x = element_text(size = 12, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    legend.position = "none"
  )

plot_png <- file.path(chronic_figures_dir, "Figure5A_current78_chronic_signature_scores_TMM_aligned.png")
plot_pdf <- file.path(chronic_figures_dir, "Figure5A_current78_chronic_signature_scores_TMM_aligned.pdf")

ggsave(plot_png, p_combined, width = 7.0, height = 8.8, dpi = 320)
ggsave(plot_pdf, p_combined, width = 7.0, height = 8.8)

plot_files <- data.frame(
  plot = c("Figure5A_current78_chronic_signature_scores_TMM_aligned", "Figure5A_current78_chronic_signature_scores_TMM_aligned"),
  file = c(plot_png, plot_pdf),
  stringsAsFactors = FALSE
)

cat("Saved plot: ", plot_png, "\n", sep = "")
cat("Saved plot: ", plot_pdf, "\n\n", sep = "")

## =========================
## 9. Save outputs
## =========================

signature_metadata_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_signature_gene_metadata_used.csv")
all_logcpm_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_all_gene_TMM_logCPM_matrix.csv")
sig_logcpm_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_signature_TMM_logCPM_matrix.csv")
zscore_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_signature_zscore_matrix.csv")
scores_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_signature_scores_by_sample.csv")
scores_long_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_signature_scores_long.csv")
stats_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_signature_score_group_comparisons.csv")
sample_metadata_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_sample_metadata_used.csv")
plot_files_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_plot_files.csv")
normalization_summary_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_normalization_summary.csv")
library_summary_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_library_summary.csv")

sig_meta <- sig_df
sig_meta$used_for_chronic_step2_score <- TRUE

write_csv(sig_meta, signature_metadata_file)
write_matrix_csv(logcpm_all, all_logcpm_file)
write_matrix_csv(sig_logcpm, sig_logcpm_file)
write_matrix_csv(z_mat, zscore_file)
write_csv(scores_df, scores_file)
write_csv(scores_long, scores_long_file)
write_csv(stats_df, stats_file)
write_csv(sample_info, sample_metadata_file)
write_csv(plot_files, plot_files_file)
write_csv(normalization_summary, normalization_summary_file)
write_csv(lib_summary, library_summary_file)

input_audit <- rbind(
  file_audit(signature_file, "current78 pig early 75-gene signature table"),
  file_audit(chronic_main_manifest_file, "chronic main comparison manifest"),
  file_audit(chronic_main_counts_file, "chronic main comparison gene-level counts matrix"),
  file_audit(step1_summary_file, "chronic Step1 detectability summary")
)
input_audit_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_input_file_audit.csv")
write_csv(input_audit, input_audit_file)

run_summary <- data.frame(
  metric = c(
    "method_version",
    "signature_input_file",
    "chronic_main_manifest_file",
    "chronic_main_counts_file",
    "n_signature_genes_used",
    "n_up_signature_genes",
    "n_down_signature_genes",
    "n_samples",
    "n_Control_52W",
    "n_ACLT_alone_52W",
    "score_scale",
    "filterByExpr_before_scoring",
    "TMM_normalization",
    "estimated_counts_rounded",
    "primary_directional_score_definition",
    "audit_directional_score_sum_means_saved",
    "group_test",
    "directional_score_median_Control_52W",
    "directional_score_median_ACLT_alone_52W",
    "directional_score_median_difference_case_minus_control",
    "directional_score_wilcox_p_value",
    "directional_score_BH_FDR_across_all_score_comparisons",
    "main_score_table",
    "group_comparison_table",
    "figure_png",
    "figure_pdf",
    "output_dir"
  ),
  value = c(
    method_version,
    signature_file,
    chronic_main_manifest_file,
    chronic_main_counts_file,
    nrow(sig_df),
    length(up_genes),
    length(down_genes),
    nrow(scores_df),
    sum(scores_df$group == "Control_52W"),
    sum(scores_df$group == "ACLT_alone_52W"),
    "edgeR TMM-normalized logCPM; row-wise z-score across 24 chronic main-comparison samples",
    "FALSE",
    "edgeR::calcNormFactors(method = 'TMM'); edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)",
    "FALSE",
    "(sum up z + sum -down z)/(n_up+n_down)",
    "TRUE; directional_score_sum_means = up_score + down_score_reoriented",
    "two-sided Wilcoxon rank-sum test for ACLT_alone_52W vs Control_52W; exact = FALSE",
    stats_df$median_control[stats_df$score == "directional_score"],
    stats_df$median_case[stats_df$score == "directional_score"],
    stats_df$median_difference_case_minus_control[stats_df$score == "directional_score"],
    stats_df$wilcox_p_value[stats_df$score == "directional_score"],
    stats_df$BH_FDR_across_all_score_comparisons[stats_df$score == "directional_score"],
    scores_file,
    stats_file,
    plot_png,
    plot_pdf,
    out_dir
  ),
  stringsAsFactors = FALSE
)
run_summary_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_signature_score_run_summary.csv")
write_csv(run_summary, run_summary_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "edgeR_version",
    "ggplot2_version",
    "signature_definition",
    "score_input",
    "normalization",
    "filtering_before_scoring",
    "row_standardization",
    "primary_score_formula",
    "chronic_analysis_role"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("ggplot2")),
    "current78-derived pig early 75-gene signature; 65 Up_in_ACLR and 10 Down_in_ACLR",
    "GSE228848 chronic main-comparison tximport gene-level estimated counts",
    "edgeR::DGEList; edgeR::calcNormFactors(method = 'TMM'); edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)",
    "No filterByExpr before scoring, to avoid losing predefined signature genes",
    "row-wise z-score across 24 chronic main-comparison samples only",
    "directional_score = (sum z_up + sum -z_down)/(n_up+n_down)",
    "extension validation only; no chronic signature/core redefinition"
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "chronic_step2_current78_TMM_aligned_versions_and_method_records.csv")
write_csv(version_df, version_file)

saveRDS(
  list(
    sig_df = sig_df,
    sample_info = sample_info,
    count_mat = count_mat,
    dge = dge,
    logcpm_all = logcpm_all,
    sig_logcpm = sig_logcpm,
    z_mat = z_mat,
    scores_df = scores_df,
    scores_long = scores_long,
    stats_df = stats_df,
    run_summary = run_summary,
    version_df = version_df,
    normalization_summary = normalization_summary,
    lib_summary = lib_summary
  ),
  file.path(obj_dir, "chronic_step2_current78_TMM_aligned_signature_score_workspace.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "chronic_step2_current78_TMM_aligned_sessionInfo.txt"))

## =========================
## 10. Console and summary-to-send log
## =========================

cat("\n===== Chronic Step2 current78 TMM-aligned signature score summary =====\n")
print(run_summary, row.names = FALSE)

cat("\nSignature direction counts:\n")
print(table(sig_df$signature_direction))

cat("\nGroup comparison statistics:\n")
print(stats_df, row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

cat("\nChronic Step2 current78 TMM-aligned signature score completed successfully.\n")

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== Chronic Step2 current78 TMM-aligned signature score SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(run_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Normalization summary:", summary_con)
writeLines(capture.output(print(normalization_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Signature direction counts:", summary_con)
writeLines(capture.output(print(table(sig_df$signature_direction))), summary_con)
writeLines("", summary_con)
writeLines("Group comparison statistics:", summary_con)
writeLines(capture.output(print(stats_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This step uses the fixed current78-derived pig early 75-gene signature in the chronic pig main comparison.", summary_con)
writeLines("It does not redefine a chronic signature or chronic core set.", summary_con)
writeLines("It uses edgeR TMM-normalized logCPM and no filterByExpr before scoring, aligned with pig early Step19.", summary_con)
writeLines("Primary manuscript-relevant score is directional_score = (sum z_up + sum -z_down)/(n_up+n_down).", summary_con)
close(summary_con)

status <- "SUCCESS"

cat("\nChronic Step2 current78 TMM-aligned signature score completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
