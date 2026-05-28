# Chronic Step5B current78: early-defined 24 core genes chronic DE audit
# Purpose:
#   Audit whether the fixed early-defined 24 core ortholog genes remain direction-consistent
#   and/or strict differentially expressed in the chronic pig main comparison:
#     ACLT_alone_52W vs Control_52W.
#
# Important interpretation:
#   - This is NOT a chronic core redefinition step.
#   - The 24 genes are fixed from the pig early primary-analysis workflow.
#   - Chronic DE audit uses the genome-wide Chronic Step4 edgeR QLF result.
#   - The audit quantifies:
#       1) detected/tested in chronic DE table
#       2) chronic direction consistency relative to early signature direction
#       3) chronic FDR < 0.05
#       4) chronic |logFC| > 1
#       5) chronic strict: FDR < 0.05 and |logFC| > 1
#       6) chronic strict + direction-consistent
#
# Direction consistency rule:
#   - early Up_in_ACLR core gene: chronic logFC > 0 is direction-consistent
#   - early Down_in_ACLR core gene: chronic logFC < 0 is direction-consistent
#
# This script does NOT auto-archive or overwrite any manually saved R script.
#
# Revision in this version:
#   - Optional Step5 heatmap summary input now points to the new supplementary TMM-aligned
#     pig-early-style compact heatmap log produced by script J.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_early_dir <- file.path(rebuild_root, "02_pig_early")
pig_chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

chronic_tables_dir <- file.path(pig_chronic_dir, "tables")

out_dir <- file.path(chronic_tables_dir, "chronic_step5B_current78_core_DE_audit_Step4_locked")
log_dir <- file.path(out_dir, "logs")
obj_dir <- file.path(out_dir, "objects")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "chronic_step5B_current78_core_DE_audit_Step4_locked_full_log.txt")
summary_log_file <- file.path(log_dir, "chronic_step5B_current78_core_DE_audit_Step4_locked_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== CHRONIC STEP5B CURRENT78 CORE DE AUDIT =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Interpretation: chronic DE audit of fixed early-defined 24 core genes; no chronic core redefinition.\n")
cat("DE source lock: final Chronic Step4 DE/GSEA output only; no fallback to old Step26.\n\n")

## =========================
## 1. Helper functions
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

direction_consistent_chronic <- function(logfc, direction) {
  out <- rep(NA, length(logfc))
  out[direction == "Up_in_ACLR" & !is.na(logfc)] <- logfc[direction == "Up_in_ACLR" & !is.na(logfc)] > 0
  out[direction == "Down_in_ACLR" & !is.na(logfc)] <- logfc[direction == "Down_in_ACLR" & !is.na(logfc)] < 0
  out
}

classify_retention <- function(direction_consistent, fdr_sig, abslogfc_gt1, strict) {
  out <- rep(NA_character_, length(direction_consistent))
  out[is.na(direction_consistent)] <- "not_tested"
  out[direction_consistent %in% FALSE] <- "opposite_or_reversed"
  out[direction_consistent %in% TRUE & strict %in% TRUE] <- "direction_consistent_strict"
  out[direction_consistent %in% TRUE & fdr_sig %in% TRUE & abslogfc_gt1 %in% FALSE] <- "direction_consistent_FDR_only"
  out[direction_consistent %in% TRUE & fdr_sig %in% FALSE & abslogfc_gt1 %in% TRUE] <- "direction_consistent_large_effect_only"
  out[direction_consistent %in% TRUE & fdr_sig %in% FALSE & abslogfc_gt1 %in% FALSE] <- "direction_consistent_but_attenuated"
  out
}

find_core_gene_file <- function() {
  candidates <- c(
    file.path(pig_early_dir, "tables", "step21_current78_Figure4C_core_heatmap", "step21_current78_pig_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_pig_core_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_core_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_audit", "step21_current78_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_pig_early_signature_validation", "step21_current78_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_pig_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_pig_core_ortholog_gene_table.csv")
  )
  hit <- find_existing_file(candidates)
  if (!is.na(hit)) return(hit)

  all_csv <- list.files(file.path(pig_early_dir, "tables"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(all_csv) == 0) return(NA_character_)

  bn <- tolower(basename(all_csv))
  score <- rep(0, length(all_csv))
  score[grepl("core", bn)] <- score[grepl("core", bn)] + 60
  score[grepl("ortholog", bn)] <- score[grepl("ortholog", bn)] + 20
  score[grepl("gene", bn)] <- score[grepl("gene", bn)] + 10
  score[grepl("table", bn)] <- score[grepl("table", bn)] + 5
  score[grepl("heatmap|zscore|matrix|summary|score|comparison|logcpm", bn)] <- score[grepl("heatmap|zscore|matrix|summary|score|comparison|logcpm", bn)] - 80

  cand <- data.frame(file = all_csv, score = score, stringsAsFactors = FALSE)
  cand <- cand[order(-cand$score, cand$file), , drop = FALSE]
  write_csv(cand, file.path(out_dir, "chronic_step5B_core_gene_file_candidates.csv"))

  if (nrow(cand) > 0 && cand$score[1] >= 70) return(cand$file[1])
  NA_character_
}

## =========================
## 2. Input files
## =========================

core_gene_file <- find_core_gene_file()

# Use ONLY the final locked chronic Step4 DE/GSEA workflow output.
# Do not fall back to the old standalone Step26 DE table, to avoid mixing deprecated
# intermediate DE results with the final submission workflow.
chronic_de_file <- file.path(
  chronic_tables_dir,
  "chronic_step4_current78_DE_GSEA",
  "chronic_step4_DE_ACLT_alone_52W_vs_Control_52W_edgeR_QLF_tximport_offset.csv"
)

step5_heatmap_summary_file <- file.path(
  chronic_tables_dir,
  "supplementary_chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact",
  "logs",
  "chronic_step5_TMM_aligned_pig_early_style_core_heatmap_summary_to_send_me.txt"
)

stop_if_missing(core_gene_file, "fixed early-defined 24 core gene table")
stop_if_missing(chronic_de_file, "final locked chronic Step4 edgeR QLF DE table")
# Heatmap summary is useful but not strictly required.
if (!file.exists(step5_heatmap_summary_file)) {
  warning("Chronic Step5 heatmap summary file not found; continuing DE audit without it.")
}

cat("Core gene file: ", core_gene_file, "\n", sep = "")
cat("Chronic DE file: ", chronic_de_file, "\n", sep = "")
cat("Step5 heatmap summary file: ", ifelse(file.exists(step5_heatmap_summary_file), step5_heatmap_summary_file, "not found"), "\n\n", sep = "")

## =========================
## 3. Read fixed early-defined core genes
## =========================

core_raw <- safe_read_csv(core_gene_file)

core_pig_col <- find_col(core_raw, c("pig_ensg", "pig_gene_id", "gene_id", "pig_gene"), TRUE, "core gene pig Ensembl ID")
core_symbol_col <- find_col(core_raw, c("pig_symbol", "gene_symbol", "symbol", "gene_name"), FALSE, "core gene pig symbol")
core_human_col <- find_col(core_raw, c("human_gene", "ortholog_name", "human_symbol", "Gene symbol"), FALSE, "core gene human symbol")
core_mouse_col <- find_col(core_raw, c("mouse_symbol", "input", "mouse_gene"), FALSE, "core gene mouse symbol")
core_direction_col <- find_col(core_raw, c("signature_direction", "mouse_direction", "direction"), TRUE, "core gene early signature direction")

core_df <- data.frame(
  pig_ensg = clean_string(core_raw[[core_pig_col]]),
  pig_symbol = if (!is.na(core_symbol_col)) clean_string(core_raw[[core_symbol_col]]) else NA_character_,
  human_gene = if (!is.na(core_human_col)) clean_string(core_raw[[core_human_col]]) else NA_character_,
  mouse_symbol = if (!is.na(core_mouse_col)) clean_string(core_raw[[core_mouse_col]]) else NA_character_,
  early_signature_direction = clean_string(core_raw[[core_direction_col]]),
  stringsAsFactors = FALSE
)

core_df <- core_df[!is.na(core_df$pig_ensg) & core_df$pig_ensg != "", , drop = FALSE]
core_df <- core_df[!duplicated(core_df$pig_ensg), , drop = FALSE]
core_df <- core_df[core_df$early_signature_direction %in% c("Up_in_ACLR", "Down_in_ACLR"), , drop = FALSE]

if (nrow(core_df) != 24) {
  stop("Expected 24 fixed early-defined core genes, but detected ", nrow(core_df), ".", call. = FALSE)
}

cat("Fixed early-defined core genes loaded: ", nrow(core_df), "\n", sep = "")
cat("Core direction counts:\n")
print(table(core_df$early_signature_direction))
cat("\n")

## =========================
## 4. Read chronic DE table and audit core genes
## =========================

de_raw <- safe_read_csv(chronic_de_file)

de_gene_col <- find_col(de_raw, c("gene_id", "pig_ensg", "GeneID", "Geneid", "ensembl_gene_id", "pig_gene_id"), TRUE, "chronic DE gene ID")
de_symbol_col <- find_col(de_raw, c("gene_symbol", "gene_name", "SYMBOL", "symbol", "pig_symbol", "external_gene_name"), FALSE, "chronic DE gene symbol")
de_logfc_col <- find_col(de_raw, c("logFC", "LogFC", "log2FC", "log2FoldChange", "log2_fold_change"), TRUE, "chronic DE logFC")
de_p_col <- find_col(de_raw, c("PValue", "P.Value", "pvalue", "P.value", "PVAL", "P", "p"), TRUE, "chronic DE P value")
de_fdr_col <- find_col(de_raw, c("FDR", "adj.P.Val", "padj", "qvalue"), TRUE, "chronic DE FDR")

de_df <- data.frame(
  pig_ensg = clean_string(de_raw[[de_gene_col]]),
  chronic_gene_symbol = if (!is.na(de_symbol_col)) clean_string(de_raw[[de_symbol_col]]) else NA_character_,
  chronic_logFC = suppressWarnings(as.numeric(de_raw[[de_logfc_col]])),
  chronic_PValue = suppressWarnings(as.numeric(de_raw[[de_p_col]])),
  chronic_FDR = suppressWarnings(as.numeric(de_raw[[de_fdr_col]])),
  stringsAsFactors = FALSE
)
de_df <- de_df[!is.na(de_df$pig_ensg) & de_df$pig_ensg != "", , drop = FALSE]
de_df <- de_df[!duplicated(de_df$pig_ensg), , drop = FALSE]

audit_df <- merge(core_df, de_df, by = "pig_ensg", all.x = TRUE, sort = FALSE)

audit_df$tested_in_chronic_DE <- !is.na(audit_df$chronic_logFC)
audit_df$chronic_direction_consistent_with_early <- direction_consistent_chronic(
  audit_df$chronic_logFC,
  audit_df$early_signature_direction
)
audit_df$chronic_FDR_lt_0.05 <- !is.na(audit_df$chronic_FDR) & audit_df$chronic_FDR < 0.05
audit_df$chronic_abslogFC_gt_1 <- !is.na(audit_df$chronic_logFC) & abs(audit_df$chronic_logFC) > 1
audit_df$chronic_strict_DEG <- audit_df$chronic_FDR_lt_0.05 & audit_df$chronic_abslogFC_gt_1
audit_df$chronic_strict_and_direction_consistent <- audit_df$chronic_strict_DEG & audit_df$chronic_direction_consistent_with_early
audit_df$chronic_retention_class <- classify_retention(
  audit_df$chronic_direction_consistent_with_early,
  audit_df$chronic_FDR_lt_0.05,
  audit_df$chronic_abslogFC_gt_1,
  audit_df$chronic_strict_DEG
)

audit_df$display_symbol <- ifelse(
  !is.na(audit_df$human_gene) & audit_df$human_gene != "",
  audit_df$human_gene,
  ifelse(!is.na(audit_df$pig_symbol) & audit_df$pig_symbol != "", audit_df$pig_symbol, audit_df$pig_ensg)
)

audit_df <- audit_df[order(
  audit_df$early_signature_direction != "Up_in_ACLR",
  !audit_df$chronic_strict_and_direction_consistent,
  -abs(audit_df$chronic_logFC),
  audit_df$display_symbol
), , drop = FALSE]

not_tested_df <- audit_df[!audit_df$tested_in_chronic_DE, , drop = FALSE]
direction_consistent_df <- audit_df[audit_df$chronic_direction_consistent_with_early %in% TRUE, , drop = FALSE]
strict_df <- audit_df[audit_df$chronic_strict_DEG %in% TRUE, , drop = FALSE]
strict_dc_df <- audit_df[audit_df$chronic_strict_and_direction_consistent %in% TRUE, , drop = FALSE]
opposite_df <- audit_df[audit_df$chronic_direction_consistent_with_early %in% FALSE, , drop = FALSE]

## =========================
## 5. Summaries
## =========================

summary_df <- data.frame(
  metric = c(
    "fixed_early_defined_core_genes",
    "core_genes_tested_in_chronic_DE",
    "core_genes_not_tested_in_chronic_DE",
    "direction_consistent_with_early",
    "opposite_or_reversed_direction",
    "chronic_FDR_lt_0.05",
    "chronic_abslogFC_gt_1",
    "chronic_strict_DEG_FDR0.05_abslogFC1",
    "chronic_strict_and_direction_consistent",
    "early_up_core_genes",
    "early_up_direction_consistent",
    "early_up_chronic_strict",
    "early_up_chronic_strict_and_direction_consistent",
    "early_down_core_genes",
    "early_down_direction_consistent",
    "early_down_chronic_strict",
    "early_down_chronic_strict_and_direction_consistent",
    "interpretation_policy",
    "output_dir"
  ),
  value = c(
    nrow(audit_df),
    sum(audit_df$tested_in_chronic_DE),
    sum(!audit_df$tested_in_chronic_DE),
    sum(audit_df$chronic_direction_consistent_with_early, na.rm = TRUE),
    sum(audit_df$chronic_direction_consistent_with_early %in% FALSE, na.rm = TRUE),
    sum(audit_df$chronic_FDR_lt_0.05, na.rm = TRUE),
    sum(audit_df$chronic_abslogFC_gt_1, na.rm = TRUE),
    sum(audit_df$chronic_strict_DEG, na.rm = TRUE),
    sum(audit_df$chronic_strict_and_direction_consistent, na.rm = TRUE),
    sum(audit_df$early_signature_direction == "Up_in_ACLR"),
    sum(audit_df$early_signature_direction == "Up_in_ACLR" & audit_df$chronic_direction_consistent_with_early, na.rm = TRUE),
    sum(audit_df$early_signature_direction == "Up_in_ACLR" & audit_df$chronic_strict_DEG, na.rm = TRUE),
    sum(audit_df$early_signature_direction == "Up_in_ACLR" & audit_df$chronic_strict_and_direction_consistent, na.rm = TRUE),
    sum(audit_df$early_signature_direction == "Down_in_ACLR"),
    sum(audit_df$early_signature_direction == "Down_in_ACLR" & audit_df$chronic_direction_consistent_with_early, na.rm = TRUE),
    sum(audit_df$early_signature_direction == "Down_in_ACLR" & audit_df$chronic_strict_DEG, na.rm = TRUE),
    sum(audit_df$early_signature_direction == "Down_in_ACLR" & audit_df$chronic_strict_and_direction_consistent, na.rm = TRUE),
    "This is a chronic DE audit of fixed early-defined core genes; it does not redefine chronic core genes.",
    out_dir
  ),
  stringsAsFactors = FALSE
)

direction_summary <- do.call(rbind, lapply(c("Up_in_ACLR", "Down_in_ACLR"), function(dir) {
  x <- audit_df[audit_df$early_signature_direction == dir, , drop = FALSE]
  data.frame(
    early_signature_direction = dir,
    n_core_genes = nrow(x),
    tested_in_chronic_DE = sum(x$tested_in_chronic_DE),
    direction_consistent_with_early = sum(x$chronic_direction_consistent_with_early, na.rm = TRUE),
    opposite_or_reversed_direction = sum(x$chronic_direction_consistent_with_early %in% FALSE, na.rm = TRUE),
    chronic_FDR_lt_0.05 = sum(x$chronic_FDR_lt_0.05, na.rm = TRUE),
    chronic_abslogFC_gt_1 = sum(x$chronic_abslogFC_gt_1, na.rm = TRUE),
    chronic_strict_DEG = sum(x$chronic_strict_DEG, na.rm = TRUE),
    chronic_strict_and_direction_consistent = sum(x$chronic_strict_and_direction_consistent, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

retention_summary <- as.data.frame(table(audit_df$chronic_retention_class), stringsAsFactors = FALSE)
colnames(retention_summary) <- c("chronic_retention_class", "n_genes")
retention_summary <- retention_summary[order(-retention_summary$n_genes, retention_summary$chronic_retention_class), , drop = FALSE]

## =========================
## 6. Write outputs
## =========================

audit_file <- file.path(out_dir, "chronic_step5B_early_defined_24_core_chronic_DE_audit_table.csv")
summary_file <- file.path(out_dir, "chronic_step5B_early_defined_24_core_chronic_DE_audit_summary.csv")
direction_summary_file <- file.path(out_dir, "chronic_step5B_core_chronic_DE_audit_by_early_direction.csv")
retention_summary_file <- file.path(out_dir, "chronic_step5B_core_retention_class_summary.csv")
strict_file <- file.path(out_dir, "chronic_step5B_core_genes_chronic_strict_DEG.csv")
strict_dc_file <- file.path(out_dir, "chronic_step5B_core_genes_chronic_strict_and_direction_consistent.csv")
direction_consistent_file <- file.path(out_dir, "chronic_step5B_core_genes_direction_consistent_with_early.csv")
opposite_file <- file.path(out_dir, "chronic_step5B_core_genes_opposite_or_reversed_in_chronic.csv")
not_tested_file <- file.path(out_dir, "chronic_step5B_core_genes_not_tested_in_chronic_DE.csv")

write_csv(audit_df, audit_file)
write_csv(summary_df, summary_file)
write_csv(direction_summary, direction_summary_file)
write_csv(retention_summary, retention_summary_file)
write_csv(strict_df, strict_file)
write_csv(strict_dc_df, strict_dc_file)
write_csv(direction_consistent_df, direction_consistent_file)
write_csv(opposite_df, opposite_file)
write_csv(not_tested_df, not_tested_file)

input_audit <- rbind(
  file_audit(core_gene_file, "fixed early-defined 24 core gene table"),
  file_audit(chronic_de_file, "final locked chronic Step4 edgeR QLF DE table"),
  file_audit(ifelse(file.exists(step5_heatmap_summary_file), step5_heatmap_summary_file, ""), "chronic Step5 heatmap summary")
)
input_audit_file <- file.path(out_dir, "chronic_step5B_input_file_audit.csv")
write_csv(input_audit, input_audit_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "method_version",
    "fixed_core_definition",
    "chronic_DE_source",
    "chronic_comparison",
    "direction_consistency_rule",
    "strict_DEG_rule",
    "analysis_role"
  ),
  value = c(
    as.character(getRversion()),
    "2026-05-27_Step5B_core_DE_audit_final_Step4_DE_source_locked_no_Step26_fallback_J_summary_path_fixed",
    "fixed early-defined 24 core ortholog genes from pig early primary-analysis workflow",
    chronic_de_file,
    "ACLT_alone_52W vs Control_52W; logFC > 0 means higher in ACLT_alone_52W",
    "Up_in_ACLR requires chronic logFC > 0; Down_in_ACLR requires chronic logFC < 0",
    "chronic FDR < 0.05 and |chronic logFC| > 1",
    "chronic DE audit only; final Step4 DE source locked; no chronic core redefinition"
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "chronic_step5B_versions_and_method_records.csv")
write_csv(version_df, version_file)

saveRDS(
  list(
    core_df = core_df,
    de_df = de_df,
    audit_df = audit_df,
    summary_df = summary_df,
    direction_summary = direction_summary,
    retention_summary = retention_summary,
    version_df = version_df
  ),
  file.path(obj_dir, "chronic_step5B_current78_core_DE_audit_workspace.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "chronic_step5B_sessionInfo.txt"))

## =========================
## 7. Console and summary-to-send log
## =========================

cat("\n===== Chronic Step5B current78 core DE audit summary =====\n")
print(summary_df, row.names = FALSE)

cat("\nDirection-level summary:\n")
print(direction_summary, row.names = FALSE)

cat("\nRetention class summary:\n")
print(retention_summary, row.names = FALSE)

cat("\nStrict and direction-consistent genes:\n")
if (nrow(strict_dc_df) == 0) {
  cat("None.\n")
} else {
  print(strict_dc_df[, c("display_symbol", "pig_ensg", "early_signature_direction", "chronic_logFC", "chronic_PValue", "chronic_FDR", "chronic_retention_class"), drop = FALSE], row.names = FALSE)
}

cat("\nOpposite/reversed genes:\n")
if (nrow(opposite_df) == 0) {
  cat("None.\n")
} else {
  print(opposite_df[, c("display_symbol", "pig_ensg", "early_signature_direction", "chronic_logFC", "chronic_PValue", "chronic_FDR", "chronic_retention_class"), drop = FALSE], row.names = FALSE)
}

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== Chronic Step5B current78 early-defined core chronic DE audit SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main summary:", summary_con)
writeLines(capture.output(print(summary_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Direction-level summary:", summary_con)
writeLines(capture.output(print(direction_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Retention class summary:", summary_con)
writeLines(capture.output(print(retention_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Strict and direction-consistent genes:", summary_con)
if (nrow(strict_dc_df) == 0) {
  writeLines("None.", summary_con)
} else {
  writeLines(capture.output(print(strict_dc_df[, c("display_symbol", "pig_ensg", "early_signature_direction", "chronic_logFC", "chronic_PValue", "chronic_FDR", "chronic_retention_class"), drop = FALSE], row.names = FALSE)), summary_con)
}
writeLines("", summary_con)
writeLines("Opposite/reversed genes:", summary_con)
if (nrow(opposite_df) == 0) {
  writeLines("None.", summary_con)
} else {
  writeLines(capture.output(print(opposite_df[, c("display_symbol", "pig_ensg", "early_signature_direction", "chronic_logFC", "chronic_PValue", "chronic_FDR", "chronic_retention_class"), drop = FALSE], row.names = FALSE)), summary_con)
}
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This is a chronic DE audit of the fixed early-defined 24 core ortholog genes.", summary_con)
writeLines("It does not redefine a chronic core set.", summary_con)
writeLines("Use chronic_strict_and_direction_consistent to describe strict retained genes, and use retention_class to describe partial retention/attenuation/reversal.", summary_con)
close(summary_con)

sink()

cat("\nChronic Step5B current78 core DE audit completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
