# ============================================================
# ACL × pig early fastp sensitivity branch
# Step H2w6: edgeR QLF differential expression on trimmed counts
#
# Purpose:
#   Use the gene-level count matrix generated from fastp-trimmed FASTQ
#   to run the same pig early edgeR QLF DE framework for:
#     1) ACLT_t7  vs Control_t0
#     2) ACLT_t28 vs Control_t0
#
# Notes:
#   - This script does NOT rerun fastp, alignment, or featureCounts.
#   - This is a sensitivity / robustness branch, not a replacement
#     for the main untrimmed analysis.
#   - Counts are featureCounts integer counts; no rounding is performed.
#   - Corrected locked-methodology version: glmQLFit(robust = TRUE),
#     matching the main untrimmed pig early DE analysis.
#   - Old Step H2w6 outputs are removed at the start to avoid mixing
#     previous robust = FALSE results with the corrected robust = TRUE run.
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n===== Step H2w6: trimmed branch edgeR QLF DE =====\n")

# ---------- Basic paths ----------
project_root <- "E:/R/ACLsenescence2"
branch_dir <- file.path(project_root, "rebuild_submission", "02_pig_early_fastp_sensitivity")
tables_dir <- file.path(branch_dir, "tables")
objects_dir <- file.path(branch_dir, "objects")
logs_dir <- file.path(branch_dir, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "stepH2w6_trimmed_edgeR_QLF_DE_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  txt <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste(..., collapse = " "))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

fail_stop <- function(...) {
  msg <- paste(..., collapse = " ")
  log_msg("ERROR:", msg)
  stop(msg, call. = FALSE)
}

# ---------- Remove previous Step H2w6 outputs ----------
# This corrected script keeps the same Step H2w6 naming but removes old H2w6
# outputs before rerunning, so previous robust = FALSE results cannot be
# accidentally mixed with the corrected robust = TRUE outputs.
old_h2w6_patterns <- c(
  file.path(tables_dir, "stepH2w6_trimmed_*.csv"),
  file.path(objects_dir, "stepH2w6_trimmed_*.rds"),
  file.path(logs_dir, "stepH2w6_trimmed_*.txt")
)
old_h2w6_files <- unique(unlist(lapply(old_h2w6_patterns, Sys.glob), use.names = FALSE))
if (length(old_h2w6_files) > 0) {
  unlink(old_h2w6_files, force = TRUE)
}

log_msg("Removed previous Step H2w6 output files:", length(old_h2w6_files))

# ---------- Input files from Step H2w5 ----------
counts_file <- file.path(
  tables_dir,
  "stepH2w5_v2_trimmed_gene_level_counts_matrix.csv"
)

sample_info_file <- file.path(
  tables_dir,
  "stepH2w5_v2_trimmed_sample_info_for_DE.csv"
)

h2w5_summary_file <- file.path(
  tables_dir,
  "stepH2w5_v2_trimmed_featureCounts_overall_summary.csv"
)

if (!file.exists(counts_file)) fail_stop("Counts file not found:", counts_file)
if (!file.exists(sample_info_file)) fail_stop("Sample info file not found:", sample_info_file)

log_msg("Project root:", project_root)
log_msg("Branch dir:", branch_dir)
log_msg("Counts file:", counts_file)
log_msg("Sample info file:", sample_info_file)

# ---------- Package checks ----------
required_pkgs <- c("edgeR")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  fail_stop(
    "Missing required package(s):",
    paste(missing_pkgs, collapse = ", "),
    ". Please install before rerunning this script."
  )
}

suppressPackageStartupMessages(library(edgeR))

edgeR_version <- as.character(utils::packageVersion("edgeR"))
log_msg("edgeR version:", edgeR_version)

# ---------- Read sample information ----------
sample_info <- read.csv(sample_info_file, check.names = FALSE, stringsAsFactors = FALSE)

required_sample_cols <- c("sample_id", "group")
missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_info))
if (length(missing_sample_cols) > 0) {
  fail_stop(
    "Sample info missing required column(s):",
    paste(missing_sample_cols, collapse = ", "),
    "\nCurrent columns:",
    paste(colnames(sample_info), collapse = ", ")
  )
}

sample_info$sample_id <- as.character(sample_info$sample_id)
sample_info$group <- as.character(sample_info$group)

expected_groups <- c("Control_t0", "ACLT_t7", "ACLT_t28")
unexpected_groups <- setdiff(unique(sample_info$group), expected_groups)
if (length(unexpected_groups) > 0) {
  fail_stop(
    "Unexpected group name(s) in sample_info:",
    paste(unexpected_groups, collapse = ", "),
    "\nExpected groups:",
    paste(expected_groups, collapse = ", ")
  )
}

sample_info <- sample_info[match(
  c(
    sample_info$sample_id[sample_info$group == "Control_t0"],
    sample_info$sample_id[sample_info$group == "ACLT_t7"],
    sample_info$sample_id[sample_info$group == "ACLT_t28"]
  ),
  sample_info$sample_id
), ]

if (any(is.na(sample_info$sample_id))) {
  fail_stop("Internal sample ordering failed. Please inspect sample_info_file:", sample_info_file)
}

# ---------- Read count matrix ----------
count_df <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)

if (ncol(count_df) < 3) {
  fail_stop("Count matrix has too few columns:", counts_file)
}

# Detect gene-id column.
possible_gene_cols <- c("gene_id", "Geneid", "GeneID", "geneid", "genes", "Gene")
gene_col <- intersect(possible_gene_cols, colnames(count_df))[1]

if (is.na(gene_col)) {
  # Fallback: use the first column if it is not a sample column.
  if (!(colnames(count_df)[1] %in% sample_info$sample_id)) {
    gene_col <- colnames(count_df)[1]
    log_msg("Gene-id column not explicitly named; using first column:", gene_col)
  } else {
    fail_stop(
      "Could not identify gene-id column in count matrix. Current columns:",
      paste(colnames(count_df), collapse = ", ")
    )
  }
} else {
  log_msg("Gene-id column detected:", gene_col)
}

gene_ids <- as.character(count_df[[gene_col]])
if (any(is.na(gene_ids)) || any(gene_ids == "")) {
  fail_stop("Gene-id column contains NA or empty values:", gene_col)
}
if (any(duplicated(gene_ids))) {
  dup_n <- sum(duplicated(gene_ids))
  fail_stop(
    "Duplicate gene IDs found in count matrix:",
    dup_n,
    ". This script stops rather than aggregating automatically."
  )
}

sample_cols_available <- intersect(sample_info$sample_id, colnames(count_df))
missing_count_samples <- setdiff(sample_info$sample_id, colnames(count_df))
if (length(missing_count_samples) > 0) {
  fail_stop(
    "Count matrix is missing sample column(s):",
    paste(missing_count_samples, collapse = ", ")
  )
}

count_mat_raw <- as.matrix(count_df[, sample_info$sample_id, drop = FALSE])
rownames(count_mat_raw) <- gene_ids

suppressWarnings(storage.mode(count_mat_raw) <- "numeric")
if (anyNA(count_mat_raw)) {
  fail_stop("Count matrix contains NA after numeric conversion.")
}
if (any(count_mat_raw < 0)) {
  fail_stop("Count matrix contains negative values.")
}

# Because featureCounts returns integer counts, this should be TRUE.
# We do NOT round. If non-integer values appear, the script stops.
non_integer_n <- sum(abs(count_mat_raw - round(count_mat_raw)) > 1e-8)
if (non_integer_n > 0) {
  fail_stop(
    "Non-integer count values detected:",
    non_integer_n,
    ". edgeR DE should use raw integer counts; no rounding is performed here."
  )
}

count_mat <- count_mat_raw
mode(count_mat) <- "integer"

sample_library_summary <- data.frame(
  sample_id = colnames(count_mat),
  group = sample_info$group,
  library_size = colSums(count_mat),
  detected_genes_count_gt0 = colSums(count_mat > 0),
  stringsAsFactors = FALSE
)

sample_library_file <- file.path(
  tables_dir,
  "stepH2w6_trimmed_count_library_summary_18samples.csv"
)
write.csv(sample_library_summary, sample_library_file, row.names = FALSE, fileEncoding = "UTF-8")

log_msg("Input genes:", nrow(count_mat))
log_msg("Input samples:", ncol(count_mat))
log_msg("Sample groups:", paste(table(sample_info$group), collapse = "; "))

# ---------- edgeR QLF DE function ----------
run_one_contrast <- function(treatment_group, control_group = "Control_t0") {
  contrast_name <- paste0(treatment_group, "_vs_", control_group)
  log_msg("Starting contrast:", contrast_name)

  keep_samples <- sample_info$group %in% c(control_group, treatment_group)
  si <- sample_info[keep_samples, , drop = FALSE]
  si$group <- factor(si$group, levels = c(control_group, treatment_group))
  si <- si[order(si$group, si$sample_id), , drop = FALSE]

  mat <- count_mat[, si$sample_id, drop = FALSE]

  y0 <- edgeR::DGEList(counts = mat, group = si$group)
  n_input_genes <- nrow(y0)

  # Pairwise filtering, matching the early pig comparison-specific logic.
  keep_gene <- edgeR::filterByExpr(y0, group = si$group)
  n_tested_genes <- sum(keep_gene)

  y <- y0[keep_gene, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")

  design <- model.matrix(~ si$group)
  colnames(design) <- c("Intercept", paste0(treatment_group, "_vs_", control_group))

  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)

  de_tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  de_tab$gene_id <- rownames(de_tab)
  de_tab <- de_tab[, c("gene_id", setdiff(colnames(de_tab), "gene_id")), drop = FALSE]

  de_tab$comparison <- contrast_name
  de_tab$direction <- ifelse(
    de_tab$logFC > 0, "up",
    ifelse(de_tab$logFC < 0, "down", "zero")
  )
  de_tab$FDR_lt_0_05 <- de_tab$FDR < 0.05
  de_tab$strict_DE_FDR_0_05_abslogFC_1 <- de_tab$FDR < 0.05 & abs(de_tab$logFC) > 1
  de_tab$strict_direction <- ifelse(
    de_tab$strict_DE_FDR_0_05_abslogFC_1 & de_tab$logFC > 1, "strict_up",
    ifelse(
      de_tab$strict_DE_FDR_0_05_abslogFC_1 & de_tab$logFC < -1, "strict_down",
      "not_strict"
    )
  )

  de_tab <- de_tab[order(de_tab$PValue, de_tab$FDR, -abs(de_tab$logFC)), , drop = FALSE]

  de_file <- file.path(
    tables_dir,
    paste0("stepH2w6_trimmed_DE_", contrast_name, "_QLF.csv")
  )
  write.csv(de_tab, de_file, row.names = FALSE, fileEncoding = "UTF-8")

  keep_file <- file.path(
    tables_dir,
    paste0("stepH2w6_trimmed_filterByExpr_keep_", contrast_name, ".csv")
  )
  keep_df <- data.frame(
    gene_id = rownames(y0$counts),
    keep_by_filterByExpr = as.logical(keep_gene),
    stringsAsFactors = FALSE
  )
  write.csv(keep_df, keep_file, row.names = FALSE, fileEncoding = "UTF-8")

  dge_rds <- file.path(
    objects_dir,
    paste0("stepH2w6_trimmed_edgeR_DGEList_", contrast_name, ".rds")
  )
  fit_rds <- file.path(
    objects_dir,
    paste0("stepH2w6_trimmed_edgeR_QLF_fit_test_", contrast_name, ".rds")
  )

  saveRDS(
    list(
      comparison = contrast_name,
      sample_info = si,
      design = design,
      dge = y,
      keep_gene = keep_gene
    ),
    dge_rds
  )

  saveRDS(
    list(
      comparison = contrast_name,
      fit = fit,
      qlf = qlf,
      de_table = de_tab
    ),
    fit_rds
  )

  summary_row <- data.frame(
    comparison = contrast_name,
    control_group = control_group,
    treatment_group = treatment_group,
    n_control_samples = sum(si$group == control_group),
    n_treatment_samples = sum(si$group == treatment_group),
    n_input_genes = n_input_genes,
    n_tested_genes_after_filterByExpr = n_tested_genes,
    n_FDR_lt_0_05 = sum(de_tab$FDR < 0.05, na.rm = TRUE),
    n_strict_FDR_0_05_abslogFC_1 = sum(de_tab$strict_DE_FDR_0_05_abslogFC_1, na.rm = TRUE),
    n_strict_up = sum(de_tab$strict_direction == "strict_up", na.rm = TRUE),
    n_strict_down = sum(de_tab$strict_direction == "strict_down", na.rm = TRUE),
    min_PValue = min(de_tab$PValue, na.rm = TRUE),
    min_FDR = min(de_tab$FDR, na.rm = TRUE),
    de_file = de_file,
    keep_file = keep_file,
    dge_rds = dge_rds,
    fit_rds = fit_rds,
    stringsAsFactors = FALSE
  )

  log_msg(
    "Finished contrast:", contrast_name,
    "| tested genes:", n_tested_genes,
    "| FDR<0.05:", summary_row$n_FDR_lt_0_05,
    "| strict:", summary_row$n_strict_FDR_0_05_abslogFC_1
  )

  return(summary_row)
}

# ---------- Run two early pig contrasts ----------
contrast_summaries <- do.call(
  rbind,
  list(
    run_one_contrast("ACLT_t7", "Control_t0"),
    run_one_contrast("ACLT_t28", "Control_t0")
  )
)

run_summary_file <- file.path(
  tables_dir,
  "stepH2w6_trimmed_edgeR_QLF_DE_run_summary.csv"
)
write.csv(contrast_summaries, run_summary_file, row.names = FALSE, fileEncoding = "UTF-8")

# ---------- Save method parameters ----------
method_params <- data.frame(
  parameter = c(
    "analysis_branch",
    "counts_source",
    "de_method",
    "comparison_strategy",
    "filtering",
    "normalization",
    "design",
    "glmQLFit_robust",
    "strict_DE_threshold",
    "edgeR_version"
  ),
  value = c(
    "pig early fastp sensitivity / trimmed FASTQ branch",
    counts_file,
    "edgeR glmQLFit + glmQLFTest",
    "pairwise comparisons: ACLT_t7 vs Control_t0; ACLT_t28 vs Control_t0",
    "edgeR::filterByExpr applied separately within each pairwise comparison",
    "TMM via edgeR::calcNormFactors",
    "~ group, treatment coefficient vs Control_t0 baseline",
    "TRUE",
    "FDR < 0.05 and |logFC| > 1",
    edgeR_version
  ),
  stringsAsFactors = FALSE
)

method_params_file <- file.path(
  tables_dir,
  "stepH2w6_trimmed_edgeR_QLF_method_parameters.csv"
)
write.csv(method_params, method_params_file, row.names = FALSE, fileEncoding = "UTF-8")

# ---------- Overall summary ----------
batch_ok_h2w6 <- all(
  file.exists(contrast_summaries$de_file),
  contrast_summaries$n_control_samples == 6,
  contrast_summaries$n_treatment_samples == 6,
  contrast_summaries$n_tested_genes_after_filterByExpr > 0
)

summary_file <- file.path(
  tables_dir,
  "stepH2w6_trimmed_edgeR_QLF_overall_summary.csv"
)

config_rds <- file.path(
  objects_dir,
  "stepH2w6_trimmed_edgeR_QLF_config.rds"
)

overall_summary <- data.frame(
  metric = c(
    "project_root",
    "branch_dir",
    "counts_file",
    "sample_info_file",
    "h2w5_summary_file",
    "edgeR_version",
    "n_samples",
    "n_input_genes",
    "n_contrasts",
    "n_t7_tested_genes",
    "n_t7_FDR_lt_0_05",
    "n_t7_strict_FDR_0_05_abslogFC_1",
    "n_t28_tested_genes",
    "n_t28_FDR_lt_0_05",
    "n_t28_strict_FDR_0_05_abslogFC_1",
    "batch_ok_h2w6",
    "sample_library_file",
    "run_summary_file",
    "method_params_file",
    "summary_file",
    "config_rds",
    "log_file",
    "status"
  ),
  value = c(
    project_root,
    branch_dir,
    counts_file,
    sample_info_file,
    h2w5_summary_file,
    edgeR_version,
    ncol(count_mat),
    nrow(count_mat),
    nrow(contrast_summaries),
    contrast_summaries$n_tested_genes_after_filterByExpr[contrast_summaries$comparison == "ACLT_t7_vs_Control_t0"],
    contrast_summaries$n_FDR_lt_0_05[contrast_summaries$comparison == "ACLT_t7_vs_Control_t0"],
    contrast_summaries$n_strict_FDR_0_05_abslogFC_1[contrast_summaries$comparison == "ACLT_t7_vs_Control_t0"],
    contrast_summaries$n_tested_genes_after_filterByExpr[contrast_summaries$comparison == "ACLT_t28_vs_Control_t0"],
    contrast_summaries$n_FDR_lt_0_05[contrast_summaries$comparison == "ACLT_t28_vs_Control_t0"],
    contrast_summaries$n_strict_FDR_0_05_abslogFC_1[contrast_summaries$comparison == "ACLT_t28_vs_Control_t0"],
    batch_ok_h2w6,
    sample_library_file,
    run_summary_file,
    method_params_file,
    summary_file,
    config_rds,
    log_file,
    if (isTRUE(batch_ok_h2w6)) {
      "Step H2w6 completed successfully"
    } else {
      "Step H2w6 completed with warnings"
    }
  ),
  stringsAsFactors = FALSE
)

write.csv(overall_summary, summary_file, row.names = FALSE, fileEncoding = "UTF-8")

saveRDS(
  list(
    project_root = project_root,
    branch_dir = branch_dir,
    counts_file = counts_file,
    sample_info_file = sample_info_file,
    sample_info = sample_info,
    count_dim = dim(count_mat),
    contrast_summaries = contrast_summaries,
    method_params = method_params,
    batch_ok_h2w6 = batch_ok_h2w6,
    created_at = Sys.time()
  ),
  config_rds
)

cat("\n===== Step H2w6 contrast summary =====\n")
print(contrast_summaries[, c(
  "comparison",
  "n_control_samples",
  "n_treatment_samples",
  "n_tested_genes_after_filterByExpr",
  "n_FDR_lt_0_05",
  "n_strict_FDR_0_05_abslogFC_1",
  "n_strict_up",
  "n_strict_down"
)], row.names = FALSE)

cat("\n===== Step H2w6 overall summary =====\n")
print(overall_summary, row.names = FALSE)

cat("\n输出文件：\n")
cat(run_summary_file, "\n")
cat(summary_file, "\n")
cat(method_params_file, "\n")
cat(sample_library_file, "\n")
cat(log_file, "\n")

if (isTRUE(batch_ok_h2w6)) {
  cat("\n结果：H2w6 成功。trimmed branch 已完成 edgeR QLF DE。下一步可以比较 untrimmed 与 trimmed 分支的 DE/score/GSEA 稳定性。\n")
} else {
  cat("\n结果：H2w6 完成但有警告。请优先查看 overall summary 和 run summary。\n")
}

cat("\n===== End of Step H2w6 =====\n")
