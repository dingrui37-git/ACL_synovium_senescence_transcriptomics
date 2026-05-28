# Chronic Step3 current78 TMM-aligned: singscore sensitivity analysis for chronic pig signature score
# Purpose:
#   Evaluate whether the chronic current78-derived 75-gene signature result depends on
#   row-wise z-score standardization by performing a singscore-based sensitivity analysis.
#
# Important interpretation:
#   - Chronic pig analysis is an extension validation, not a discovery/redefinition step.
#   - This step does NOT redefine the signature.
#   - The fixed signature is the pig early current78-derived 75-gene signature:
#       65 Up_in_ACLR genes + 10 Down_in_ACLR genes.
#   - The chronic main comparison is Control_52W vs ACLT_alone_52W, 12 vs 12.
#   - This singscore analysis is a sensitivity analysis for Chronic Step2, not a replacement
#     of the primary z-score based signature score.
#
# Singscore implementation:
#   - Build genome-wide edgeR TMM-normalized logCPM matrix for the same 24 chronic main-comparison samples.
#   - Use singscore::rankGenes() to generate within-sample expression ranks.
#   - Compute:
#       1) directional_singscore using upSet = 65 up genes and downSet = 10 down genes,
#          centerScore = TRUE, knownDirection = TRUE.
#       2) up_singscore using up genes alone.
#       3) down_singscore_raw by treating down genes as an upSet.
#       4) down_singscore_reoriented = -down_singscore_raw.
#   - Compare ACLT_alone_52W vs Control_52W using two-sided Wilcoxon rank-sum test.
#
# This script also compares singscore outputs with Chronic Step2 z-score outputs.
# This script does NOT auto-archive or overwrite any manually saved R script.

options(stringsAsFactors = FALSE)

method_version <- "2026-05-16_chronic_step3_current78_TMM_aligned_singscore_sensitivity_v1"

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
chronic_figures_dir <- file.path(pig_chronic_dir, "figures", "Supplementary_chronic_current78_singscore_TMM_aligned")

out_dir <- file.path(chronic_tables_dir, "chronic_step3_current78_singscore_sensitivity_TMM_aligned")
log_dir <- file.path(out_dir, "logs")
obj_dir <- file.path(out_dir, "objects")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chronic_figures_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "chronic_step3_current78_TMM_aligned_singscore_sensitivity_full_log.txt")
summary_log_file <- file.path(log_dir, "chronic_step3_current78_TMM_aligned_singscore_sensitivity_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== CHRONIC STEP3 CURRENT78 TMM-ALIGNED SINGSCORE SENSITIVITY =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Interpretation: singscore sensitivity analysis using fixed pig early current78 75-gene signature.\n")
cat("Input transformation: edgeR TMM-normalized logCPM; no filterByExpr before ranking/scoring.\n\n")

## =========================
## 1. Packages
## =========================

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("Required package is not installed: edgeR.", call. = FALSE)
}
if (!requireNamespace("singscore", quietly = TRUE)) {
  stop(
    "Required package is not installed: singscore. Please install it with BiocManager::install('singscore') and rerun.",
    call. = FALSE
  )
}
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

suppressPackageStartupMessages({
  library(edgeR)
  library(singscore)
  library(ggplot2)
})

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

get_singscore_col <- function(score_df, label = "singscore result") {
  candidates <- c("TotalScore", "Score", "Total.Score", "score")
  hit <- candidates[candidates %in% colnames(score_df)]
  if (length(hit) > 0) return(hit[1])

  numeric_cols <- colnames(score_df)[vapply(score_df, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("Dispersion", "NumGenes", "nGenes", "Rank"))
  if (length(numeric_cols) == 0) {
    stop("Could not identify score column in ", label, ". Available columns: ", paste(colnames(score_df), collapse = ", "), call. = FALSE)
  }
  numeric_cols[1]
}

fmt_p <- function(p) {
  ifelse(is.na(p), "p = NA", ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3))))
}

compare_one <- function(df, score_col) {
  control_vals <- df[df$group == "Control_52W", score_col]
  case_vals <- df[df$group == "ACLT_alone_52W", score_col]
  wt <- suppressWarnings(wilcox.test(case_vals, control_vals, exact = FALSE, alternative = "two.sided"))
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

step2_scores_file <- find_existing_file(c(
  file.path(
    chronic_tables_dir,
    "chronic_step2_current78_signature_score_TMM_aligned",
    "chronic_step2_current78_TMM_aligned_signature_scores_by_sample.csv"
  ),
  file.path(
    chronic_tables_dir,
    "chronic_step2_current78_signature_score_TMM_aligned",
    "chronic_step2_current78_signature_scores_by_sample.csv"
  )
))

stop_if_missing(signature_file, "current78 pig early 75-gene signature table")
if (is.na(chronic_main_manifest_file)) stop("Could not find chronic main comparison manifest under: ", chronic_tables_dir, call. = FALSE)
if (is.na(chronic_main_counts_file)) stop("Could not find chronic main comparison counts matrix under: ", chronic_tables_dir, call. = FALSE)
if (is.na(step2_scores_file)) {
  stop("Could not find updated TMM-aligned Chronic Step2 z-score signature scores. Please run chronic_step2_current78_TMM_aligned_signature_score.R first.", call. = FALSE)
}

cat("Signature file: ", signature_file, "\n", sep = "")
cat("Chronic main manifest: ", chronic_main_manifest_file, "\n", sep = "")
cat("Chronic main count matrix: ", chronic_main_counts_file, "\n", sep = "")
cat("Updated TMM-aligned Step2 z-score scores: ", step2_scores_file, "\n\n", sep = "")

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
  write_csv(data.frame(missing_sample_id = missing_samples), file.path(out_dir, "chronic_step3_missing_samples_in_count_matrix.csv"))
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

missing_sig <- setdiff(sig_df$pig_ensg, rownames(count_mat))
if (length(missing_sig) > 0) {
  write_csv(sig_df[sig_df$pig_ensg %in% missing_sig, , drop = FALSE], file.path(out_dir, "chronic_step3_signature_genes_missing_in_count_matrix.csv"))
  stop("Some fixed signature genes are missing in chronic count matrix.", call. = FALSE)
}

## =========================
## 6. Genome-wide TMM-normalized logCPM and singscore
## =========================

dge <- edgeR::DGEList(counts = count_mat, group = sample_info$group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm_mat <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 1,
  normalized.lib.sizes = TRUE
)

normalization_summary <- data.frame(
  metric = c(
    "input_count_matrix_genes",
    "input_count_matrix_samples",
    "filterByExpr_applied_before_singscore",
    "edgeR_DGEList_used",
    "TMM_normalization",
    "logCPM_function",
    "prior_count",
    "normalized_lib_sizes",
    "estimated_counts_rounded",
    "n_signature_genes_available_for_scoring",
    "n_up_signature_genes",
    "n_down_signature_genes"
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
    "FALSE",
    nrow(sig_df),
    length(up_genes),
    length(down_genes)
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

# Singscore ranks genes within each sample.
ranked <- singscore::rankGenes(logcpm_mat)

directional_sc <- singscore::simpleScore(
  rankData = ranked,
  upSet = up_genes,
  downSet = down_genes,
  centerScore = TRUE,
  knownDirection = TRUE
)

up_sc <- singscore::simpleScore(
  rankData = ranked,
  upSet = up_genes,
  centerScore = TRUE,
  knownDirection = TRUE
)

down_raw_sc <- singscore::simpleScore(
  rankData = ranked,
  upSet = down_genes,
  centerScore = TRUE,
  knownDirection = TRUE
)

directional_col <- get_singscore_col(directional_sc, "directional singscore")
up_col <- get_singscore_col(up_sc, "up singscore")
down_col <- get_singscore_col(down_raw_sc, "down raw singscore")

# Most singscore outputs keep samples in rownames.
sample_ids_from_score <- rownames(directional_sc)
if (is.null(sample_ids_from_score) || any(is.na(sample_ids_from_score)) || any(sample_ids_from_score == "")) {
  if ("sample_id" %in% colnames(directional_sc)) {
    sample_ids_from_score <- as.character(directional_sc$sample_id)
  } else {
    stop("Could not identify sample IDs from singscore output.", call. = FALSE)
  }
}

scores_df <- data.frame(
  sample_id = sample_ids_from_score,
  directional_singscore = as.numeric(directional_sc[[directional_col]]),
  up_singscore = as.numeric(up_sc[[up_col]]),
  down_singscore_raw = as.numeric(down_raw_sc[[down_col]]),
  stringsAsFactors = FALSE
)
scores_df$down_singscore_reoriented <- -scores_df$down_singscore_raw

scores_df <- merge(
  data.frame(sample_id = sample_info$sample_id, group = as.character(sample_info$group), group_raw = sample_info$group_raw, stringsAsFactors = FALSE),
  scores_df,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(scores_df$directional_singscore)) {
  stop("NA directional singscore detected after merging sample metadata.", call. = FALSE)
}

scores_df$group <- factor(scores_df$group, levels = c("Control_52W", "ACLT_alone_52W"))
scores_df <- scores_df[order(scores_df$group, scores_df$sample_id), , drop = FALSE]

score_cols <- c("directional_singscore", "up_singscore", "down_singscore_reoriented", "down_singscore_raw")
stats_df <- do.call(rbind, lapply(score_cols, function(sc) compare_one(scores_df, sc)))
stats_df$BH_FDR_across_all_singscore_comparisons <- p.adjust(stats_df$wilcox_p_value, method = "BH")

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
## 7. Compare with Chronic Step2 z-score outputs
## =========================

step2_scores <- safe_read_csv(step2_scores_file)
if (!("sample_id" %in% colnames(step2_scores))) stop("TMM-aligned Step2 score file lacks sample_id column.", call. = FALSE)

merge_z_sing <- merge(
  step2_scores[, c("sample_id", "directional_score", "up_score", "down_score_reoriented"), drop = FALSE],
  scores_df[, c("sample_id", "directional_singscore", "up_singscore", "down_singscore_reoriented"), drop = FALSE],
  by = "sample_id",
  all = FALSE
)

cor_df <- data.frame(
  comparison = c("TMM_aligned_directional_zscore_vs_singscore", "TMM_aligned_up_zscore_vs_singscore", "TMM_aligned_down_reoriented_zscore_vs_singscore"),
  n_samples = c(nrow(merge_z_sing), nrow(merge_z_sing), nrow(merge_z_sing)),
  spearman = c(
    suppressWarnings(cor(merge_z_sing$directional_score, merge_z_sing$directional_singscore, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(merge_z_sing$up_score, merge_z_sing$up_singscore, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(merge_z_sing$down_score_reoriented, merge_z_sing$down_singscore_reoriented, method = "spearman", use = "complete.obs"))
  ),
  pearson = c(
    suppressWarnings(cor(merge_z_sing$directional_score, merge_z_sing$directional_singscore, method = "pearson", use = "complete.obs")),
    suppressWarnings(cor(merge_z_sing$up_score, merge_z_sing$up_singscore, method = "pearson", use = "complete.obs")),
    suppressWarnings(cor(merge_z_sing$down_score_reoriented, merge_z_sing$down_singscore_reoriented, method = "pearson", use = "complete.obs"))
  ),
  stringsAsFactors = FALSE
)

## =========================
## 8. Supplementary plot
## =========================

plot_files <- data.frame(plot = character(), file = character(), stringsAsFactors = FALSE)

score_plot_levels <- c("directional_singscore", "up_singscore", "down_singscore_reoriented")
score_plot_labels <- c(
  directional_singscore = "Directional singscore",
  up_singscore = "Up-signature singscore",
  down_singscore_reoriented = "Down-signature singscore reoriented"
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
    title = "Chronic pig current78 signature singscore sensitivity analysis",
    subtitle = "TMM-normalized logCPM; rank-based scoring",
    x = NULL,
    y = "Centered singscore"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 17),
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

plot_png <- file.path(chronic_figures_dir, "Supplementary_chronic_current78_TMM_aligned_singscore_scores_mouse_style.png")
plot_pdf <- file.path(chronic_figures_dir, "Supplementary_chronic_current78_TMM_aligned_singscore_scores_mouse_style.pdf")

ggsave(plot_png, p_combined, width = 7.0, height = 8.8, dpi = 320)
ggsave(plot_pdf, p_combined, width = 7.0, height = 8.8)

plot_files <- rbind(plot_files, data.frame(plot = "Supplementary_chronic_current78_TMM_aligned_singscore_scores_mouse_style", file = plot_png, stringsAsFactors = FALSE))
plot_files <- rbind(plot_files, data.frame(plot = "Supplementary_chronic_current78_TMM_aligned_singscore_scores_mouse_style", file = plot_pdf, stringsAsFactors = FALSE))

## =========================
## 9. Save outputs
## =========================

signature_metadata_file <- file.path(out_dir, "chronic_step3_current78_TMM_aligned_signature_gene_metadata_used.csv")
genome_logcpm_file <- file.path(out_dir, "chronic_step3_TMM_aligned_genomewide_logCPM_matrix_for_singscore.csv")
normalization_summary_file <- file.path(out_dir, "chronic_step3_TMM_aligned_normalization_summary.csv")
library_summary_file <- file.path(out_dir, "chronic_step3_TMM_aligned_library_summary.csv")
scores_file <- file.path(out_dir, "chronic_step3_current78_TMM_aligned_singscore_scores_by_sample.csv")
scores_long_file <- file.path(out_dir, "chronic_step3_current78_TMM_aligned_singscore_scores_long.csv")
stats_file <- file.path(out_dir, "chronic_step3_current78_TMM_aligned_singscore_group_comparisons.csv")
cor_file <- file.path(out_dir, "chronic_step3_current78_TMM_aligned_singscore_vs_zscore_correlation.csv")
sample_metadata_file <- file.path(out_dir, "chronic_step3_TMM_aligned_sample_metadata_used.csv")
plot_files_file <- file.path(out_dir, "chronic_step3_TMM_aligned_plot_files.csv")

sig_meta <- sig_df
sig_meta$used_for_chronic_step3_singscore <- TRUE

write_csv(sig_meta, signature_metadata_file)
write_matrix_csv(logcpm_mat, genome_logcpm_file)
write_csv(normalization_summary, normalization_summary_file)
write_csv(lib_summary, library_summary_file)
write_csv(scores_df, scores_file)
write_csv(scores_long, scores_long_file)
write_csv(stats_df, stats_file)
write_csv(cor_df, cor_file)
write_csv(sample_info, sample_metadata_file)
write_csv(plot_files, plot_files_file)

input_audit <- rbind(
  file_audit(signature_file, "current78 pig early 75-gene signature table"),
  file_audit(chronic_main_manifest_file, "chronic main comparison manifest"),
  file_audit(chronic_main_counts_file, "chronic main comparison gene-level counts matrix"),
  file_audit(step2_scores_file, "updated TMM-aligned chronic Step2 z-score signature scores")
)
input_audit_file <- file.path(out_dir, "chronic_step3_TMM_aligned_input_file_audit.csv")
write_csv(input_audit, input_audit_file)

run_summary <- data.frame(
  metric = c(
    "method_version",
    "signature_input_file",
    "chronic_main_manifest_file",
    "chronic_main_counts_file",
    "updated_TMM_aligned_step2_scores_file",
    "n_signature_genes_used",
    "n_up_signature_genes",
    "n_down_signature_genes",
    "n_genomewide_genes_ranked",
    "n_samples",
    "n_Control_52W",
    "n_ACLT_alone_52W",
    "singscore_scale",
    "filterByExpr_before_singscore",
    "TMM_normalization",
    "estimated_counts_rounded",
    "directional_singscore_definition",
    "down_singscore_reorientation",
    "group_test",
    "directional_singscore_median_Control_52W",
    "directional_singscore_median_ACLT_alone_52W",
    "directional_singscore_median_difference_case_minus_control",
    "directional_singscore_wilcox_p_value",
    "directional_singscore_BH_FDR_across_all_singscore_comparisons",
    "TMM_aligned_directional_zscore_vs_singscore_spearman",
    "main_singscore_table",
    "group_comparison_table",
    "zscore_singscore_correlation_table",
    "figure_png",
    "figure_pdf",
    "output_dir"
  ),
  value = c(
    method_version,
    signature_file,
    chronic_main_manifest_file,
    chronic_main_counts_file,
    step2_scores_file,
    nrow(sig_df),
    length(up_genes),
    length(down_genes),
    nrow(logcpm_mat),
    nrow(scores_df),
    sum(scores_df$group == "Control_52W"),
    sum(scores_df$group == "ACLT_alone_52W"),
    "genome-wide sample-wise ranks from TMM-normalized logCPM; centerScore=TRUE",
    "FALSE",
    "edgeR::calcNormFactors(method = 'TMM'); edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)",
    "FALSE",
    "simpleScore(rankData, upSet=65 genes, downSet=10 genes, centerScore=TRUE, knownDirection=TRUE)",
    "down genes scored as upSet to obtain raw score, then multiplied by -1",
    "two-sided Wilcoxon rank-sum test for ACLT_alone_52W vs Control_52W; exact = FALSE",
    stats_df$median_control[stats_df$score == "directional_singscore"],
    stats_df$median_case[stats_df$score == "directional_singscore"],
    stats_df$median_difference_case_minus_control[stats_df$score == "directional_singscore"],
    stats_df$wilcox_p_value[stats_df$score == "directional_singscore"],
    stats_df$BH_FDR_across_all_singscore_comparisons[stats_df$score == "directional_singscore"],
    cor_df$spearman[cor_df$comparison == "TMM_aligned_directional_zscore_vs_singscore"],
    scores_file,
    stats_file,
    cor_file,
    plot_png,
    plot_pdf,
    out_dir
  ),
  stringsAsFactors = FALSE
)
run_summary_file <- file.path(out_dir, "chronic_step3_current78_TMM_aligned_singscore_run_summary.csv")
write_csv(run_summary, run_summary_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "edgeR_version",
    "singscore_version",
    "ggplot2_version",
    "signature_definition",
    "score_input",
    "rank_input",
    "filtering_before_ranking",
    "centerScore",
    "knownDirection",
    "chronic_analysis_role"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("singscore")),
    as.character(utils::packageVersion("ggplot2")),
    "current78-derived pig early 75-gene signature; 65 Up_in_ACLR and 10 Down_in_ACLR",
    "GSE228848 chronic main-comparison tximport gene-level estimated counts",
    "genome-wide TMM-normalized logCPM matrix for 24 chronic main-comparison samples",
    "No filterByExpr before singscore ranking/scoring, to retain predefined signature genes",
    TRUE,
    TRUE,
    "sensitivity analysis only; no chronic signature/core redefinition"
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "chronic_step3_TMM_aligned_versions_and_method_records.csv")
write_csv(version_df, version_file)

saveRDS(
  list(
    sig_df = sig_df,
    sample_info = sample_info,
    logcpm_mat = logcpm_mat,
    ranked = ranked,
    directional_sc = directional_sc,
    up_sc = up_sc,
    down_raw_sc = down_raw_sc,
    scores_df = scores_df,
    scores_long = scores_long,
    stats_df = stats_df,
    merge_z_sing = merge_z_sing,
    cor_df = cor_df,
    run_summary = run_summary,
    version_df = version_df
  ),
  file.path(obj_dir, "chronic_step3_current78_TMM_aligned_singscore_sensitivity_workspace.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "chronic_step3_TMM_aligned_sessionInfo.txt"))

## =========================
## 10. Console and summary-to-send log
## =========================

cat("\n===== Chronic Step3 current78 TMM-aligned singscore summary =====\n")
print(run_summary, row.names = FALSE)

cat("\nNormalization summary:\n")
print(normalization_summary, row.names = FALSE)

cat("\nSignature direction counts:\n")
print(table(sig_df$signature_direction))

cat("\nSingscore group comparison statistics:\n")
print(stats_df, row.names = FALSE)

cat("\nTMM-aligned z-score vs singscore correlation:\n")
print(cor_df, row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

cat("\nChronic Step3 current78 TMM-aligned singscore sensitivity completed successfully.\n")

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== Chronic Step3 current78 TMM-aligned singscore sensitivity SUMMARY TO SEND ME =====", summary_con)
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
writeLines("Singscore group comparison statistics:", summary_con)
writeLines(capture.output(print(stats_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("TMM-aligned z-score vs singscore correlation:", summary_con)
writeLines(capture.output(print(cor_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This step uses the fixed current78-derived pig early 75-gene signature in the chronic pig main comparison.", summary_con)
writeLines("It is a sensitivity analysis for the updated TMM-aligned Chronic Step2 z-score scoring and does not redefine a chronic signature/core set.", summary_con)
writeLines("It uses genome-wide edgeR TMM-normalized logCPM and no filterByExpr before singscore ranking/scoring.", summary_con)
writeLines("Primary sensitivity readout is directional_singscore.", summary_con)
close(summary_con)

status <- "SUCCESS"

sink()

cat("\nChronic Step3 current78 TMM-aligned singscore sensitivity completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
