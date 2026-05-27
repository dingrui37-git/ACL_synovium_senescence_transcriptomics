# ============================================================
# Step H2w7: pig early fastp sensitivity
# Compare untrimmed vs trimmed edgeR QLF DE stability
# Project root: E:/R/ACLsenescence2
# Author: ChatGPT
# ============================================================

# 这一步是干什么：
# 1) 精确读取主分析 untrimmed Step17 pig early edgeR QLF DE 结果；
# 2) 精确读取 fastp sensitivity StepH2w6 trimmed branch edgeR QLF DE 结果；
# 3) 分别比较 ACLT_t7 vs Control_t0 和 ACLT_t28 vs Control_t0；
# 4) 输出 logFC 相关性、strict DEG 重叠、方向一致性等稳定性指标。
#
# Locked update:
# - 不再递归搜索或打分选择 DE 文件。
# - 只允许读取本脚本中明确定义的 4 个 DE 输入文件。

options(stringsAsFactors = FALSE)

# ---------- Basic paths ----------
project_root <- "E:/R/ACLsenescence2"
branch_dir <- file.path(project_root, "rebuild_submission", "02_pig_early_fastp_sensitivity")
untrimmed_tables_dir <- file.path(project_root, "rebuild_submission", "02_pig_early", "tables")
trimmed_tables_dir <- file.path(branch_dir, "tables")

tables_dir <- file.path(branch_dir, "tables")
logs_dir <- file.path(branch_dir, "logs")
objects_dir <- file.path(branch_dir, "objects")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "stepH2w7_trimmed_vs_untrimmed_DE_stability_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste(..., collapse = " ")
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

stop_with_log <- function(...) {
  msg <- paste(..., collapse = " ")
  log_msg("ERROR:", msg)
  stop(msg, call. = FALSE)
}

log_msg("===== Step H2w7 started =====")
log_msg("project_root:", project_root)
log_msg("untrimmed_tables_dir:", untrimmed_tables_dir)
log_msg("trimmed_tables_dir:", trimmed_tables_dir)

if (!dir.exists(untrimmed_tables_dir)) {
  stop_with_log("Untrimmed pig early tables dir not found:", untrimmed_tables_dir)
}
if (!dir.exists(trimmed_tables_dir)) {
  stop_with_log("Trimmed branch tables dir not found:", trimmed_tables_dir)
}

# ---------- Helpers ----------
find_col <- function(df, candidates) {
  nm <- names(df)
  nm_lower <- tolower(nm)
  cand_lower <- tolower(candidates)
  hit <- which(nm_lower %in% cand_lower)
  if (length(hit) > 0) return(nm[hit[1]])
  # allow loose matching for some common variants
  for (cand in cand_lower) {
    hit2 <- which(nm_lower == cand | grepl(paste0("^", cand, "$"), nm_lower))
    if (length(hit2) > 0) return(nm[hit2[1]])
  }
  NA_character_
}

is_probably_de_table <- function(file) {
  x <- tryCatch(read.csv(file, nrows = 5, check.names = FALSE), error = function(e) NULL)
  if (is.null(x) || ncol(x) < 5) return(FALSE)
  logfc_col <- find_col(x, c("logFC", "log2FC", "log2FoldChange"))
  p_col <- find_col(x, c("PValue", "P.Value", "P.Value.", "pvalue", "p_value", "P"))
  fdr_col <- find_col(x, c("FDR", "padj", "adj.P.Val", "adj.P.Val.", "qvalue", "q_value"))
  all(!is.na(c(logfc_col, p_col, fdr_col)))
}

score_candidate <- function(file, contrast_key, source_type) {
  b <- tolower(basename(file))
  p <- tolower(file)
  score <- 0

  # DE / QLF evidence
  if (grepl("de|deg|edger|qlf|glmqlf", b)) score <- score + 20
  if (grepl("qlf|glmqlf|edger", b)) score <- score + 20
  if (grepl("aclt|acl", b)) score <- score + 10
  if (grepl("control|ctrl|con", b)) score <- score + 5

  # contrast evidence
  if (contrast_key == "t7") {
    if (grepl("t7|day7|d7|_7_", b)) score <- score + 50
    if (grepl("t28|day28|d28|_28_", b)) score <- score - 50
  }
  if (contrast_key == "t28") {
    if (grepl("t28|day28|d28|_28_", b)) score <- score + 50
    if (grepl("t7|day7|d7|_7_", b)) score <- score - 50
  }

  # source-type evidence
  if (source_type == "trimmed") {
    if (grepl("h2w6|trimmed|fastp_sensitivity", p)) score <- score + 30
  } else {
    if (grepl("fastp_sensitivity|trimmed|h2w", p)) score <- score - 200
    if (grepl("pig_early", p)) score <- score + 10
  }

  # exclude non-DE summary tables
  if (grepl("summary|parameter|method|library|assignment|featurecounts|counts_matrix|manifest|qc|stat|sample_info|candidate|status|overall", b)) {
    score <- score - 80
  }
  score
}

list_de_candidates <- function(dir, contrast_key, source_type) {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) {
    return(data.frame())
  }

  out <- data.frame(
    source_type = source_type,
    contrast_key = contrast_key,
    file = files,
    basename = basename(files),
    score = vapply(files, score_candidate, numeric(1), contrast_key = contrast_key, source_type = source_type),
    valid_de_table = vapply(files, is_probably_de_table, logical(1)),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$valid_de_table, -out$score, out$basename), ]
  row.names(out) <- NULL
  out
}

choose_best_de_file <- function(dir, contrast_key, source_type) {
  cand <- list_de_candidates(dir, contrast_key, source_type)
  if (nrow(cand) == 0) {
    return(list(file = NA_character_, candidates = cand))
  }
  usable <- cand[cand$valid_de_table & cand$score > 0, , drop = FALSE]
  if (nrow(usable) == 0) {
    return(list(file = NA_character_, candidates = cand))
  }
  list(file = usable$file[1], candidates = cand)
}

infer_gene_col <- function(df) {
  gene_candidates <- c(
    "gene_id", "GeneID", "gene", "Gene", "genes", "Genes",
    "ensembl_gene_id", "ENSEMBL", "Ensembl", "id", "ID",
    "X", "...1", "rowname", "rownames"
  )
  hit <- find_col(df, gene_candidates)
  if (!is.na(hit)) return(hit)

  stat_cols <- tolower(c("logFC", "log2FC", "log2FoldChange", "logCPM", "F", "LR", "PValue", "P.Value", "pvalue", "FDR", "padj", "adj.P.Val"))
  for (nm in names(df)) {
    v <- df[[nm]]
    if (!tolower(nm) %in% stat_cols && !is.numeric(v)) {
      non_na <- sum(!is.na(v) & v != "")
      unique_n <- length(unique(v[!is.na(v) & v != ""]))
      if (non_na > 0 && unique_n / max(non_na, 1) > 0.8) return(nm)
    }
  }
  NA_character_
}

read_de_table <- function(file, label) {
  if (!file.exists(file)) stop_with_log(label, "file not found:", file)
  df <- read.csv(file, check.names = FALSE)
  gene_col <- infer_gene_col(df)
  logfc_col <- find_col(df, c("logFC", "log2FC", "log2FoldChange"))
  p_col <- find_col(df, c("PValue", "P.Value", "P.Value.", "pvalue", "p_value", "P"))
  fdr_col <- find_col(df, c("FDR", "padj", "adj.P.Val", "adj.P.Val.", "qvalue", "q_value"))

  if (any(is.na(c(gene_col, logfc_col, p_col, fdr_col)))) {
    stop_with_log(
      label, "does not contain required DE columns. File:", file,
      "| inferred gene_col=", gene_col,
      "logFC_col=", logfc_col,
      "P_col=", p_col,
      "FDR_col=", fdr_col,
      "| current columns:", paste(names(df), collapse = ", ")
    )
  }

  out <- data.frame(
    gene_id = as.character(df[[gene_col]]),
    logFC = suppressWarnings(as.numeric(df[[logfc_col]])),
    PValue = suppressWarnings(as.numeric(df[[p_col]])),
    FDR = suppressWarnings(as.numeric(df[[fdr_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene_id) & out$gene_id != "", , drop = FALSE]
  out <- out[is.finite(out$logFC) & is.finite(out$PValue) & is.finite(out$FDR), , drop = FALSE]
  out$PValue[out$PValue <= 0] <- .Machine$double.xmin
  out$FDR[out$FDR <= 0] <- .Machine$double.xmin

  # If duplicate genes exist, keep the row with the smallest FDR, then PValue.
  out <- out[order(out$gene_id, out$FDR, out$PValue), ]
  out <- out[!duplicated(out$gene_id), ]
  row.names(out) <- NULL

  attr(out, "gene_col") <- gene_col
  attr(out, "logfc_col") <- logfc_col
  attr(out, "p_col") <- p_col
  attr(out, "fdr_col") <- fdr_col
  out
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

pct <- function(x) {
  if (is.na(x)) return(NA_real_)
  round(100 * x, 3)
}

compare_one_contrast <- function(contrast_label, contrast_key, untrimmed_file, trimmed_file) {
  log_msg("Comparing contrast:", contrast_label)
  log_msg("  untrimmed_file:", untrimmed_file)
  log_msg("  trimmed_file:", trimmed_file)

  u <- read_de_table(untrimmed_file, paste0("untrimmed ", contrast_label))
  t <- read_de_table(trimmed_file, paste0("trimmed ", contrast_label))

  u$strict <- u$FDR < 0.05 & abs(u$logFC) > 1
  t$strict <- t$FDR < 0.05 & abs(t$logFC) > 1
  u$minuslog10P <- -log10(u$PValue)
  t$minuslog10P <- -log10(t$PValue)

  m <- merge(
    u, t,
    by = "gene_id",
    suffixes = c("_untrimmed", "_trimmed"),
    all = FALSE
  )
  m$logFC_diff_trimmed_minus_untrimmed <- m$logFC_trimmed - m$logFC_untrimmed
  m$abs_logFC_diff <- abs(m$logFC_diff_trimmed_minus_untrimmed)
  m$FDR_diff_trimmed_minus_untrimmed <- m$FDR_trimmed - m$FDR_untrimmed
  m$strict_status <- ifelse(
    m$strict_untrimmed & m$strict_trimmed, "strict_in_both",
    ifelse(m$strict_untrimmed & !m$strict_trimmed, "strict_untrimmed_only",
           ifelse(!m$strict_untrimmed & m$strict_trimmed, "strict_trimmed_only", "not_strict_in_either"))
  )
  m$direction_same <- sign(m$logFC_untrimmed) == sign(m$logFC_trimmed)

  strict_u <- u$gene_id[u$strict]
  strict_t <- t$gene_id[t$strict]
  strict_intersection <- intersect(strict_u, strict_t)
  strict_union <- union(strict_u, strict_t)

  top_overlap <- function(n) {
    top_u <- u$gene_id[order(u$FDR, u$PValue, -abs(u$logFC))][seq_len(min(n, nrow(u)))]
    top_t <- t$gene_id[order(t$FDR, t$PValue, -abs(t$logFC))][seq_len(min(n, nrow(t)))]
    length(intersect(top_u, top_t))
  }

  idx_common_nonzero <- sign(m$logFC_untrimmed) != 0 & sign(m$logFC_trimmed) != 0
  idx_strict_union <- m$gene_id %in% strict_union

  summary <- data.frame(
    contrast = contrast_label,
    untrimmed_file = untrimmed_file,
    trimmed_file = trimmed_file,
    n_untrimmed_genes = nrow(u),
    n_trimmed_genes = nrow(t),
    n_common_genes = nrow(m),
    pearson_logFC_common = round(safe_cor(m$logFC_untrimmed, m$logFC_trimmed, "pearson"), 6),
    spearman_logFC_common = round(safe_cor(m$logFC_untrimmed, m$logFC_trimmed, "spearman"), 6),
    pearson_minuslog10P_common = round(safe_cor(m$minuslog10P_untrimmed, m$minuslog10P_trimmed, "pearson"), 6),
    n_strict_untrimmed = length(strict_u),
    n_strict_trimmed = length(strict_t),
    n_strict_intersection = length(strict_intersection),
    n_strict_union = length(strict_union),
    strict_jaccard = round(length(strict_intersection) / max(length(strict_union), 1), 6),
    pct_untrimmed_strict_recovered_in_trimmed = pct(length(strict_intersection) / max(length(strict_u), 1)),
    pct_trimmed_strict_also_in_untrimmed = pct(length(strict_intersection) / max(length(strict_t), 1)),
    direction_agreement_pct_common = pct(mean(m$direction_same[idx_common_nonzero], na.rm = TRUE)),
    direction_agreement_pct_strict_union = pct(mean(m$direction_same[idx_strict_union], na.rm = TRUE)),
    median_abs_logFC_diff_common = round(median(m$abs_logFC_diff, na.rm = TRUE), 6),
    p95_abs_logFC_diff_common = round(as.numeric(quantile(m$abs_logFC_diff, 0.95, na.rm = TRUE)), 6),
    top100_FDR_overlap = top_overlap(100),
    top500_FDR_overlap = top_overlap(500),
    stringsAsFactors = FALSE
  )

  gene_comparison_file <- file.path(
    tables_dir,
    paste0("stepH2w7_DE_stability_", contrast_label, "_gene_level_comparison.csv")
  )
  strict_membership_file <- file.path(
    tables_dir,
    paste0("stepH2w7_DE_stability_", contrast_label, "_strict_membership.csv")
  )

  write.csv(m, gene_comparison_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(m[m$strict_status != "not_strict_in_either", ], strict_membership_file, row.names = FALSE, fileEncoding = "UTF-8")

  log_msg("  wrote gene comparison:", gene_comparison_file)
  log_msg("  wrote strict membership:", strict_membership_file)

  list(summary = summary, gene_file = gene_comparison_file, strict_file = strict_membership_file)
}

# ---------- Select DE files: locked exact paths, no recursive search ----------
# This corrected version intentionally does NOT use recursive/scored file discovery.
# It reads only the locked primary untrimmed Step17 DE files and the locked trimmed StepH2w6 DE files.
contrasts <- data.frame(
  contrast_label = c("ACLT_t7_vs_Control_t0", "ACLT_t28_vs_Control_t0"),
  contrast_key = c("t7", "t28"),
  stringsAsFactors = FALSE
)

selected <- data.frame(
  contrast = rep(contrasts$contrast_label, each = 2),
  source_type = rep(c("untrimmed", "trimmed"), times = nrow(contrasts)),
  selected_file = c(
    file.path(untrimmed_tables_dir, "step17_pig_early_DE_t7_vs_CON_QLF.csv"),
    file.path(trimmed_tables_dir, "stepH2w6_trimmed_DE_ACLT_t7_vs_Control_t0_QLF.csv"),
    file.path(untrimmed_tables_dir, "step17_pig_early_DE_t28_vs_CON_QLF.csv"),
    file.path(trimmed_tables_dir, "stepH2w6_trimmed_DE_ACLT_t28_vs_Control_t0_QLF.csv")
  ),
  selection_mode = "locked_exact_path_no_recursive_search",
  stringsAsFactors = FALSE
)

# Keep the historical output filenames so downstream checking code remains compatible.
# The candidate file now records only the locked exact paths rather than scored candidates.
candidate_file <- file.path(tables_dir, "stepH2w7_DE_file_selection_candidates.csv")
selected_file <- file.path(tables_dir, "stepH2w7_DE_file_selection_selected_files.csv")
write.csv(selected, candidate_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(selected, selected_file, row.names = FALSE, fileEncoding = "UTF-8")

log_msg("Locked exact DE input files written to:", selected_file)
print(selected)

missing_locked <- selected[!file.exists(selected$selected_file), , drop = FALSE]
if (nrow(missing_locked) > 0) {
  print(missing_locked)
  stop_with_log(
    "Locked exact DE input file(s) not found. No fallback search is allowed in this locked version. Inspect:",
    selected_file
  )
}

# ---------- Compare contrasts ----------
results <- list()
for (i in seq_len(nrow(contrasts))) {
  cl <- contrasts$contrast_label[i]
  ck <- contrasts$contrast_key[i]
  u_file <- selected$selected_file[selected$contrast == cl & selected$source_type == "untrimmed"]
  t_file <- selected$selected_file[selected$contrast == cl & selected$source_type == "trimmed"]
  results[[cl]] <- compare_one_contrast(cl, ck, u_file, t_file)
}

summary_df <- do.call(rbind, lapply(results, `[[`, "summary"))
row.names(summary_df) <- NULL

summary_file <- file.path(tables_dir, "stepH2w7_trimmed_vs_untrimmed_DE_stability_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE, fileEncoding = "UTF-8")

# Compact interpretation flags: conservative, not used as formal statistics.
interpretation_df <- summary_df
interpretation_df$stability_flag <- ifelse(
  interpretation_df$spearman_logFC_common >= 0.98 &
    interpretation_df$strict_jaccard >= 0.95 &
    interpretation_df$direction_agreement_pct_strict_union >= 99,
  "highly_stable",
  ifelse(
    interpretation_df$spearman_logFC_common >= 0.95 &
      interpretation_df$strict_jaccard >= 0.90,
    "stable_with_minor_differences",
    "needs_manual_review"
  )
)
interpretation_file <- file.path(tables_dir, "stepH2w7_trimmed_vs_untrimmed_DE_stability_interpretation_flags.csv")
write.csv(interpretation_df, interpretation_file, row.names = FALSE, fileEncoding = "UTF-8")

# Overall summary for the run.
batch_ok_h2w7 <- all(!is.na(summary_df$spearman_logFC_common)) && all(summary_df$n_common_genes > 10000)
overall_summary <- data.frame(
  metric = c(
    "project_root",
    "branch_dir",
    "untrimmed_tables_dir",
    "trimmed_tables_dir",
    "n_contrasts_compared",
    "min_common_genes",
    "min_spearman_logFC_common",
    "min_strict_jaccard",
    "min_direction_agreement_pct_strict_union",
    "batch_ok_h2w7",
    "candidate_file",
    "selected_file",
    "summary_file",
    "interpretation_file",
    "log_file",
    "status"
  ),
  value = c(
    project_root,
    branch_dir,
    untrimmed_tables_dir,
    trimmed_tables_dir,
    nrow(summary_df),
    min(summary_df$n_common_genes, na.rm = TRUE),
    min(summary_df$spearman_logFC_common, na.rm = TRUE),
    min(summary_df$strict_jaccard, na.rm = TRUE),
    min(summary_df$direction_agreement_pct_strict_union, na.rm = TRUE),
    batch_ok_h2w7,
    candidate_file,
    selected_file,
    summary_file,
    interpretation_file,
    log_file,
    ifelse(batch_ok_h2w7, "Step H2w7 completed successfully", "Step H2w7 completed but needs review")
  ),
  stringsAsFactors = FALSE
)

overall_summary_file <- file.path(tables_dir, "stepH2w7_trimmed_vs_untrimmed_DE_stability_overall_summary.csv")
write.csv(overall_summary, overall_summary_file, row.names = FALSE, fileEncoding = "UTF-8")

config_rds <- file.path(objects_dir, "stepH2w7_trimmed_vs_untrimmed_DE_stability_config.rds")
saveRDS(
  list(
    project_root = project_root,
    branch_dir = branch_dir,
    untrimmed_tables_dir = untrimmed_tables_dir,
    trimmed_tables_dir = trimmed_tables_dir,
    selected = selected,
    summary_df = summary_df,
    created_at = Sys.time()
  ),
  config_rds
)

cat("\n===== Step H2w7 DE stability summary =====\n")
print(summary_df, row.names = FALSE)

cat("\n===== Step H2w7 interpretation flags =====\n")
print(interpretation_df[, c("contrast", "stability_flag", "spearman_logFC_common", "strict_jaccard", "direction_agreement_pct_strict_union")], row.names = FALSE)

cat("\n===== Step H2w7 overall summary =====\n")
print(overall_summary, row.names = FALSE)

cat("\n输出文件：\n")
cat(summary_file, "\n")
cat(interpretation_file, "\n")
cat(candidate_file, "\n")
cat(selected_file, "\n")
cat(overall_summary_file, "\n")
cat(config_rds, "\n")
cat(log_file, "\n")

if (isTRUE(batch_ok_h2w7)) {
  cat("\n结果：H2w7 成功。已完成 trimmed 与 untrimmed 分支的 DE 稳定性比较。下一步可以进入 signature score / GSEA 稳定性比较。\n")
} else {
  cat("\n结果：H2w7 完成，但建议人工查看 summary 中的相关性和 strict DEG 重叠指标。\n")
}

log_msg("===== Step H2w7 finished =====")
cat("\n===== End of Step H2w7 =====\n")
