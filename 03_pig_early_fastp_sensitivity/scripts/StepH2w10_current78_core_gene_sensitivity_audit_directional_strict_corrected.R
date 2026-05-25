# StepH2w10_current78: core gene sensitivity audit in the fastp-trimmed pig early branch
# Purpose:
# This step audits gene-level validation robustness after fastp trimming.
#
# IMPORTANT INTERPRETATION:
#   - The fastp branch is a sensitivity-analysis workflow, not a replacement of the main untrimmed analysis.
#   - Main core ortholog genes remain defined by the primary untrimmed current78 analysis.
#   - This script does NOT redefine the official main core set.
#   - It tests whether the current78 signature-level and main-core gene-level validation structure is preserved
#     in the trimmed branch.
#
# v2 update for manuscript-method consistency:
#   - In the main analysis, strict_t7 / strict_t28 are directional strict flags:
#       FDR < 0.05, |logFC| > 1, and direction consistent with the mouse-derived signature direction.
#   - This script therefore keeps separate stat_strict_* flags for the statistical threshold alone,
#     and defines strict_t7_trimmed / strict_t28_trimmed as directional strict flags.
#
# The script outputs two key audit layers:
#   1) Trimmed branch validation summary for the current78-derived 75 pig signature genes.
#   2) Retention/preservation audit for the main-analysis 24 core ortholog genes in the trimmed branch.
#
# Inputs:
#   Main current78 Step18 signature table:
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
#   Main current78 Step18 validation table:
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_validation_table.csv
#   Trimmed count matrix:
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_gene_level_counts_matrix.csv
#   Trimmed QLF DE tables:
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t7_vs_Control_t0_QLF.csv
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t28_vs_Control_t0_QLF.csv
#
# This script does NOT auto-archive or overwrite any manually saved R script.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

main_pig_dir <- file.path(rebuild_root, "02_pig_early")
fastp_dir <- file.path(rebuild_root, "02_pig_early_fastp_sensitivity")

main_step18_dir <- file.path(main_pig_dir, "tables", "step18_current78_pig_early_signature_remap")
fastp_tables_dir <- file.path(fastp_dir, "tables")

out_dir <- file.path(fastp_tables_dir, "stepH2w10_current78_core_gene_sensitivity_audit")
log_dir <- file.path(out_dir, "logs")
obj_dir <- file.path(out_dir, "objects")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "stepH2w10_current78_core_gene_sensitivity_audit_full_log.txt")
summary_log_file <- file.path(log_dir, "stepH2w10_current78_core_gene_sensitivity_audit_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== StepH2w10_current78 core gene sensitivity audit v2 =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Interpretation: fastp branch is sensitivity audit only; main core set remains defined by untrimmed primary analysis.\n\n")

## =========================
## 1. Helper functions
## =========================

stop_if_missing <- function(file, label) {
  if (!file.exists(file)) stop("Missing ", label, ": ", file, call. = FALSE)
}

safe_read_csv <- function(file) {
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
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

find_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Could not find ", label, ". Checked: ", paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  NA_character_
}

as_logical_robust <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  y <- tolower(clean_string(x))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "yes", "y", "1")] <- TRUE
  out[y %in% c("false", "f", "no", "n", "0")] <- FALSE
  out
}

read_count_gene_ids <- function(file) {
  df <- safe_read_csv(file)
  gene_col <- find_col(df, c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id"), TRUE, "count matrix gene id")
  gene_ids <- clean_string(df[[gene_col]])
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  unique(gene_ids)
}

standardize_de <- function(file, contrast_label) {
  df <- safe_read_csv(file)

  gene_col <- find_col(df, c("gene_id", "pig_ensg", "GeneID", "Geneid", "feature_id", "ensembl_gene_id"), TRUE, paste0(contrast_label, " gene id"))
  symbol_col <- find_col(df, c("gene_name", "gene_symbol", "SYMBOL", "symbol", "pig_symbol", "external_gene_name"), FALSE, paste0(contrast_label, " gene symbol"))
  logfc_col <- find_col(df, c("logFC", "LogFC", "log2FC", "log2FoldChange", "log2_fold_change"), TRUE, paste0(contrast_label, " logFC"))
  p_col <- find_col(df, c("PValue", "P.Value", "pvalue", "P.value", "PVAL", "P", "p"), FALSE, paste0(contrast_label, " P value"))
  fdr_col <- find_col(df, c("FDR", "adj.P.Val", "padj", "qvalue"), TRUE, paste0(contrast_label, " FDR"))

  out <- data.frame(
    pig_ensg = clean_string(df[[gene_col]]),
    pig_symbol_from_trimmed_DE = if (!is.na(symbol_col)) clean_string(df[[symbol_col]]) else NA_character_,
    logFC = suppressWarnings(as.numeric(df[[logfc_col]])),
    PValue = if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_,
    FDR = suppressWarnings(as.numeric(df[[fdr_col]])),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$pig_ensg) & out$pig_ensg != "", , drop = FALSE]
  out <- out[!duplicated(out$pig_ensg), , drop = FALSE]

  attr(out, "gene_col") <- gene_col
  attr(out, "symbol_col") <- symbol_col
  attr(out, "logfc_col") <- logfc_col
  attr(out, "p_col") <- p_col
  attr(out, "fdr_col") <- fdr_col

  out
}

direction_consistent <- function(logfc, signature_direction) {
  out <- rep(NA, length(logfc))
  out[signature_direction == "Up_in_ACLR" & !is.na(logfc)] <- logfc[signature_direction == "Up_in_ACLR" & !is.na(logfc)] > 0
  out[signature_direction == "Down_in_ACLR" & !is.na(logfc)] <- logfc[signature_direction == "Down_in_ACLR" & !is.na(logfc)] < 0
  out
}

strict_flag <- function(logfc, fdr) {
  !is.na(logfc) & !is.na(fdr) & fdr < 0.05 & abs(logfc) > 1
}

## =========================
## 2. Input files
## =========================

signature_file <- file.path(main_step18_dir, "step18_current78_pig_signature_gene_table.csv")
main_validation_file <- file.path(main_step18_dir, "step18_current78_pig_signature_validation_table.csv")

trimmed_counts_file <- file.path(fastp_tables_dir, "stepH2w5_v2_trimmed_gene_level_counts_matrix.csv")
trimmed_t7_de_file <- file.path(fastp_tables_dir, "stepH2w6_trimmed_DE_ACLT_t7_vs_Control_t0_QLF.csv")
trimmed_t28_de_file <- file.path(fastp_tables_dir, "stepH2w6_trimmed_DE_ACLT_t28_vs_Control_t0_QLF.csv")

stop_if_missing(signature_file, "main current78 signature table")
stop_if_missing(main_validation_file, "main current78 validation table")
stop_if_missing(trimmed_counts_file, "trimmed gene count matrix")
stop_if_missing(trimmed_t7_de_file, "trimmed t7 QLF DE table")
stop_if_missing(trimmed_t28_de_file, "trimmed t28 QLF DE table")

cat("Signature file: ", signature_file, "\n", sep = "")
cat("Main validation file: ", main_validation_file, "\n", sep = "")
cat("Trimmed count matrix: ", trimmed_counts_file, "\n", sep = "")
cat("Trimmed t7 DE file: ", trimmed_t7_de_file, "\n", sep = "")
cat("Trimmed t28 DE file: ", trimmed_t28_de_file, "\n\n", sep = "")

## =========================
## 3. Read main signature and main core set
## =========================

sig_raw <- safe_read_csv(signature_file)
main_val_raw <- safe_read_csv(main_validation_file)

sig_pig_col <- find_col(sig_raw, c("pig_ensg", "pig_gene_id", "pig_gene", "gene_id"), TRUE, "signature pig Ensembl ID")
sig_symbol_col <- find_col(sig_raw, c("pig_symbol", "pig_gene_symbol", "gene_symbol", "symbol", "gene_name"), FALSE, "signature pig symbol")
sig_human_col <- find_col(sig_raw, c("human_gene", "ortholog_name", "Gene symbol", "gene_symbol_human"), FALSE, "signature human gene")
sig_mouse_col <- find_col(sig_raw, c("mouse_symbol", "input", "mouse_gene", "mouse_gene_symbol"), FALSE, "signature mouse symbol")
sig_dir_col <- find_col(sig_raw, c("signature_direction", "mouse_direction", "direction"), TRUE, "signature direction")

signature_df <- data.frame(
  pig_ensg = clean_string(sig_raw[[sig_pig_col]]),
  pig_symbol = if (!is.na(sig_symbol_col)) clean_string(sig_raw[[sig_symbol_col]]) else NA_character_,
  human_gene = if (!is.na(sig_human_col)) clean_string(sig_raw[[sig_human_col]]) else NA_character_,
  mouse_symbol = if (!is.na(sig_mouse_col)) clean_string(sig_raw[[sig_mouse_col]]) else NA_character_,
  signature_direction = clean_string(sig_raw[[sig_dir_col]]),
  stringsAsFactors = FALSE
)
signature_df <- signature_df[!is.na(signature_df$pig_ensg) & signature_df$pig_ensg != "" &
                               signature_df$signature_direction %in% c("Up_in_ACLR", "Down_in_ACLR"), , drop = FALSE]
signature_df <- signature_df[!duplicated(signature_df$pig_ensg), , drop = FALSE]

if (nrow(signature_df) != 75 ||
    sum(signature_df$signature_direction == "Up_in_ACLR") != 65 ||
    sum(signature_df$signature_direction == "Down_in_ACLR") != 10) {
  stop("Signature table does not match expected current78 75/65/10 counts. Please check input.", call. = FALSE)
}

main_pig_col <- find_col(main_val_raw, c("pig_ensg", "pig_gene_id", "pig_gene", "gene_id"), TRUE, "main validation pig Ensembl ID")
main_dir_col <- find_col(main_val_raw, c("signature_direction", "mouse_direction", "direction"), TRUE, "main validation direction")
main_dc_col <- find_col(main_val_raw, c("direction_consistent_both_timepoints", "direction_consistent_both", "direction_consistent"), TRUE, "main validation direction-consistent both")
main_strict_t7_col <- find_col(main_val_raw, c("strict_t7", "strict_at_t7"), TRUE, "main validation strict_t7")
main_strict_t28_col <- find_col(main_val_raw, c("strict_t28", "strict_at_t28"), TRUE, "main validation strict_t28")

main_core_df <- main_val_raw
main_core_df$pig_ensg_std <- clean_string(main_core_df[[main_pig_col]])
main_core_df$signature_direction_std <- clean_string(main_core_df[[main_dir_col]])
main_core_df$main_direction_consistent_both <- as_logical_robust(main_core_df[[main_dc_col]])
main_core_df$main_strict_t7 <- as_logical_robust(main_core_df[[main_strict_t7_col]])
main_core_df$main_strict_t28 <- as_logical_robust(main_core_df[[main_strict_t28_col]])
main_core_df$main_core <- main_core_df$main_direction_consistent_both & main_core_df$main_strict_t7 & main_core_df$main_strict_t28
main_core_ids <- unique(main_core_df$pig_ensg_std[main_core_df$main_core %in% TRUE])
main_core_ids <- main_core_ids[!is.na(main_core_ids) & main_core_ids != ""]

if (length(main_core_ids) != 24) {
  stop("Main-analysis core gene count is not 24 from the validation table. Detected: ", length(main_core_ids), call. = FALSE)
}

cat("Main current78 signature genes: ", nrow(signature_df), "\n", sep = "")
cat("Main current78 up/down: ", sum(signature_df$signature_direction == "Up_in_ACLR"), "/", sum(signature_df$signature_direction == "Down_in_ACLR"), "\n", sep = "")
cat("Main-analysis core genes: ", length(main_core_ids), "\n\n", sep = "")

## =========================
## 4. Read trimmed count/DE and compute audit flags
## =========================

trimmed_count_genes <- read_count_gene_ids(trimmed_counts_file)
trimmed_t7 <- standardize_de(trimmed_t7_de_file, "trimmed_t7")
trimmed_t28 <- standardize_de(trimmed_t28_de_file, "trimmed_t28")

cat("Trimmed count genes: ", length(trimmed_count_genes), "\n", sep = "")
cat("Trimmed t7 DE rows: ", nrow(trimmed_t7), "\n", sep = "")
cat("Trimmed t28 DE rows: ", nrow(trimmed_t28), "\n\n", sep = "")

trimmed_validation <- signature_df
trimmed_validation$detected_in_trimmed_counts <- trimmed_validation$pig_ensg %in% trimmed_count_genes

trimmed_validation <- merge(
  trimmed_validation,
  trimmed_t7[, c("pig_ensg", "pig_symbol_from_trimmed_DE", "logFC", "PValue", "FDR")],
  by = "pig_ensg",
  all.x = TRUE,
  sort = FALSE
)
colnames(trimmed_validation)[colnames(trimmed_validation) == "pig_symbol_from_trimmed_DE"] <- "pig_symbol_from_trimmed_DE_t7"
colnames(trimmed_validation)[colnames(trimmed_validation) == "logFC"] <- "logFC_t7_trimmed"
colnames(trimmed_validation)[colnames(trimmed_validation) == "PValue"] <- "PValue_t7_trimmed"
colnames(trimmed_validation)[colnames(trimmed_validation) == "FDR"] <- "FDR_t7_trimmed"

trimmed_validation <- merge(
  trimmed_validation,
  trimmed_t28[, c("pig_ensg", "pig_symbol_from_trimmed_DE", "logFC", "PValue", "FDR")],
  by = "pig_ensg",
  all.x = TRUE,
  sort = FALSE
)
colnames(trimmed_validation)[colnames(trimmed_validation) == "pig_symbol_from_trimmed_DE"] <- "pig_symbol_from_trimmed_DE_t28"
colnames(trimmed_validation)[colnames(trimmed_validation) == "logFC"] <- "logFC_t28_trimmed"
colnames(trimmed_validation)[colnames(trimmed_validation) == "PValue"] <- "PValue_t28_trimmed"
colnames(trimmed_validation)[colnames(trimmed_validation) == "FDR"] <- "FDR_t28_trimmed"

trimmed_validation$tested_in_trimmed_t7_DE <- !is.na(trimmed_validation$logFC_t7_trimmed)
trimmed_validation$tested_in_trimmed_t28_DE <- !is.na(trimmed_validation$logFC_t28_trimmed)

trimmed_validation$direction_consistent_t7_trimmed <- direction_consistent(
  trimmed_validation$logFC_t7_trimmed,
  trimmed_validation$signature_direction
)
trimmed_validation$direction_consistent_t28_trimmed <- direction_consistent(
  trimmed_validation$logFC_t28_trimmed,
  trimmed_validation$signature_direction
)
trimmed_validation$direction_consistent_both_trimmed <- trimmed_validation$direction_consistent_t7_trimmed & trimmed_validation$direction_consistent_t28_trimmed

# Statistical strict threshold only: FDR < 0.05 and |logFC| > 1.
# These columns are kept for auditing and should not be interpreted as the manuscript strict_t7/strict_t28 definition.
trimmed_validation$stat_strict_t7_trimmed <- strict_flag(trimmed_validation$logFC_t7_trimmed, trimmed_validation$FDR_t7_trimmed)
trimmed_validation$stat_strict_t28_trimmed <- strict_flag(trimmed_validation$logFC_t28_trimmed, trimmed_validation$FDR_t28_trimmed)
trimmed_validation$stat_strict_both_trimmed <- trimmed_validation$stat_strict_t7_trimmed & trimmed_validation$stat_strict_t28_trimmed

# Directional strict flags matching the primary untrimmed manuscript definition:
# strict_t7/strict_t28 = statistical strict threshold + direction consistency with the mouse-derived signature.
trimmed_validation$strict_t7_trimmed <- trimmed_validation$stat_strict_t7_trimmed & trimmed_validation$direction_consistent_t7_trimmed
trimmed_validation$strict_t28_trimmed <- trimmed_validation$stat_strict_t28_trimmed & trimmed_validation$direction_consistent_t28_trimmed
trimmed_validation$strict_both_trimmed <- trimmed_validation$strict_t7_trimmed & trimmed_validation$strict_t28_trimmed

trimmed_validation$core_like_trimmed <- trimmed_validation$direction_consistent_both_trimmed &
  trimmed_validation$strict_t7_trimmed & trimmed_validation$strict_t28_trimmed

trimmed_validation$main_core_gene <- trimmed_validation$pig_ensg %in% main_core_ids

# Add main-analysis validation flags for side-by-side comparison if available
main_flags <- data.frame(
  pig_ensg = clean_string(main_val_raw[[main_pig_col]]),
  main_direction_consistent_both = as_logical_robust(main_val_raw[[main_dc_col]]),
  main_strict_t7 = as_logical_robust(main_val_raw[[main_strict_t7_col]]),
  main_strict_t28 = as_logical_robust(main_val_raw[[main_strict_t28_col]]),
  stringsAsFactors = FALSE
)
main_flags$main_core <- main_flags$main_direction_consistent_both & main_flags$main_strict_t7 & main_flags$main_strict_t28
main_flags <- main_flags[!duplicated(main_flags$pig_ensg), , drop = FALSE]

trimmed_validation <- merge(trimmed_validation, main_flags, by = "pig_ensg", all.x = TRUE, sort = FALSE)

## =========================
## 5. Summaries
## =========================

main_core_trimmed <- trimmed_validation[trimmed_validation$main_core_gene, , drop = FALSE]

summary_df <- data.frame(
  metric = c(
    "current78_signature_genes",
    "current78_up_signature_genes",
    "current78_down_signature_genes",
    "trimmed_signature_genes_detected_in_counts",
    "trimmed_signature_genes_tested_in_t7_DE",
    "trimmed_signature_genes_tested_in_t28_DE",
    "trimmed_direction_consistent_t7",
    "trimmed_direction_consistent_t28",
    "trimmed_direction_consistent_both_timepoints",
    "trimmed_stat_strict_t7_FDR_abslogFC_only",
    "trimmed_stat_strict_t28_FDR_abslogFC_only",
    "trimmed_stat_strict_both_FDR_abslogFC_only",
    "trimmed_strict_t7_direction_consistent",
    "trimmed_strict_t28_direction_consistent",
    "trimmed_strict_both_direction_consistent",
    "trimmed_core_like_genes",
    "main_core_genes_from_untrimmed_primary_analysis",
    "main_core_genes_detected_in_trimmed_counts",
    "main_core_genes_tested_in_trimmed_t7_DE",
    "main_core_genes_tested_in_trimmed_t28_DE",
    "main_core_genes_direction_consistent_t7_trimmed",
    "main_core_genes_direction_consistent_t28_trimmed",
    "main_core_genes_direction_consistent_both_trimmed",
    "main_core_genes_stat_strict_t7_FDR_abslogFC_only_trimmed",
    "main_core_genes_stat_strict_t28_FDR_abslogFC_only_trimmed",
    "main_core_genes_stat_strict_both_FDR_abslogFC_only_trimmed",
    "main_core_genes_strict_t7_direction_consistent_trimmed",
    "main_core_genes_strict_t28_direction_consistent_trimmed",
    "main_core_genes_strict_both_direction_consistent_trimmed",
    "main_core_genes_preserved_as_core_like_in_trimmed",
    "official_interpretation",
    "output_dir"
  ),
  value = c(
    nrow(trimmed_validation),
    sum(trimmed_validation$signature_direction == "Up_in_ACLR"),
    sum(trimmed_validation$signature_direction == "Down_in_ACLR"),
    sum(trimmed_validation$detected_in_trimmed_counts),
    sum(trimmed_validation$tested_in_trimmed_t7_DE),
    sum(trimmed_validation$tested_in_trimmed_t28_DE),
    sum(trimmed_validation$direction_consistent_t7_trimmed, na.rm = TRUE),
    sum(trimmed_validation$direction_consistent_t28_trimmed, na.rm = TRUE),
    sum(trimmed_validation$direction_consistent_both_trimmed, na.rm = TRUE),
    sum(trimmed_validation$stat_strict_t7_trimmed, na.rm = TRUE),
    sum(trimmed_validation$stat_strict_t28_trimmed, na.rm = TRUE),
    sum(trimmed_validation$stat_strict_both_trimmed, na.rm = TRUE),
    sum(trimmed_validation$strict_t7_trimmed, na.rm = TRUE),
    sum(trimmed_validation$strict_t28_trimmed, na.rm = TRUE),
    sum(trimmed_validation$strict_both_trimmed, na.rm = TRUE),
    sum(trimmed_validation$core_like_trimmed, na.rm = TRUE),
    length(main_core_ids),
    sum(main_core_trimmed$detected_in_trimmed_counts),
    sum(main_core_trimmed$tested_in_trimmed_t7_DE),
    sum(main_core_trimmed$tested_in_trimmed_t28_DE),
    sum(main_core_trimmed$direction_consistent_t7_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$direction_consistent_t28_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$direction_consistent_both_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$stat_strict_t7_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$stat_strict_t28_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$stat_strict_both_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$strict_t7_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$strict_t28_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$strict_both_trimmed, na.rm = TRUE),
    sum(main_core_trimmed$core_like_trimmed, na.rm = TRUE),
    "fastp-trimmed branch is a sensitivity audit only; official core genes remain defined by the untrimmed primary analysis",
    out_dir
  ),
  stringsAsFactors = FALSE
)

direction_summary <- do.call(rbind, lapply(c("Up_in_ACLR", "Down_in_ACLR"), function(dir) {
  x <- trimmed_validation[trimmed_validation$signature_direction == dir, , drop = FALSE]
  data.frame(
    signature_direction = dir,
    n_signature_genes = nrow(x),
    detected_in_trimmed_counts = sum(x$detected_in_trimmed_counts),
    direction_consistent_t7_trimmed = sum(x$direction_consistent_t7_trimmed, na.rm = TRUE),
    direction_consistent_t28_trimmed = sum(x$direction_consistent_t28_trimmed, na.rm = TRUE),
    direction_consistent_both_trimmed = sum(x$direction_consistent_both_trimmed, na.rm = TRUE),
    stat_strict_t7_FDR_abslogFC_only_trimmed = sum(x$stat_strict_t7_trimmed, na.rm = TRUE),
    stat_strict_t28_FDR_abslogFC_only_trimmed = sum(x$stat_strict_t28_trimmed, na.rm = TRUE),
    stat_strict_both_FDR_abslogFC_only_trimmed = sum(x$stat_strict_both_trimmed, na.rm = TRUE),
    strict_t7_direction_consistent_trimmed = sum(x$strict_t7_trimmed, na.rm = TRUE),
    strict_t28_direction_consistent_trimmed = sum(x$strict_t28_trimmed, na.rm = TRUE),
    strict_both_direction_consistent_trimmed = sum(x$strict_both_trimmed, na.rm = TRUE),
    core_like_trimmed = sum(x$core_like_trimmed, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

main_core_preservation_table <- main_core_trimmed[
  order(!main_core_trimmed$core_like_trimmed, main_core_trimmed$signature_direction, main_core_trimmed$pig_symbol),
  c(
    "pig_ensg", "pig_symbol", "human_gene", "mouse_symbol", "signature_direction",
    "detected_in_trimmed_counts",
    "logFC_t7_trimmed", "FDR_t7_trimmed", "direction_consistent_t7_trimmed", "stat_strict_t7_trimmed", "strict_t7_trimmed",
    "logFC_t28_trimmed", "FDR_t28_trimmed", "direction_consistent_t28_trimmed", "stat_strict_t28_trimmed", "strict_t28_trimmed",
    "direction_consistent_both_trimmed", "stat_strict_both_trimmed", "strict_both_trimmed", "core_like_trimmed",
    "main_direction_consistent_both", "main_strict_t7", "main_strict_t28", "main_core"
  ),
  drop = FALSE
]

changed_main_core <- main_core_preservation_table[!main_core_preservation_table$core_like_trimmed, , drop = FALSE]

## =========================
## 6. Write outputs
## =========================

trimmed_validation_file <- file.path(out_dir, "stepH2w10_current78_trimmed_signature_gene_level_validation_table.csv")
summary_file <- file.path(out_dir, "stepH2w10_current78_core_gene_sensitivity_audit_summary.csv")
direction_summary_file <- file.path(out_dir, "stepH2w10_current78_trimmed_gene_level_summary_by_signature_direction.csv")
main_core_preservation_file <- file.path(out_dir, "stepH2w10_current78_main24_core_genes_trimmed_preservation_table.csv")
changed_main_core_file <- file.path(out_dir, "stepH2w10_current78_main24_core_genes_not_preserved_as_core_like_in_trimmed.csv")

write_csv(trimmed_validation, trimmed_validation_file)
write_csv(summary_df, summary_file)
write_csv(direction_summary, direction_summary_file)
write_csv(main_core_preservation_table, main_core_preservation_file)
write_csv(changed_main_core, changed_main_core_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "signature_source",
    "main_core_source",
    "trimmed_count_source",
    "trimmed_DE_source",
    "direction_consistency_rule",
    "strict_rule",
    "core_like_trimmed_rule",
    "interpretation"
  ),
  value = c(
    as.character(getRversion()),
    signature_file,
    main_validation_file,
    trimmed_counts_file,
    paste(trimmed_t7_de_file, trimmed_t28_de_file, sep = "; "),
    "Up_in_ACLR requires logFC > 0; Down_in_ACLR requires logFC < 0",
    "strict_t7/strict_t28 include FDR < 0.05, |logFC| > 1, and direction consistency with the mouse-derived signature; stat_strict_* records FDR/|logFC| only",
    "direction-consistent at both t7/t28 and directional strict at both t7/t28 in the trimmed branch",
    "This is a sensitivity audit only; it does not redefine the official main core gene set."
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "stepH2w10_current78_versions_and_method_records.csv")
write_csv(version_df, version_file)

saveRDS(
  list(
    signature_df = signature_df,
    main_core_ids = main_core_ids,
    trimmed_validation = trimmed_validation,
    main_core_preservation_table = main_core_preservation_table,
    changed_main_core = changed_main_core,
    summary_df = summary_df,
    direction_summary = direction_summary,
    version_df = version_df
  ),
  file.path(obj_dir, "stepH2w10_current78_core_gene_sensitivity_audit_workspace.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "stepH2w10_current78_sessionInfo.txt"))

cat("\n===== StepH2w10_current78 summary =====\n")
print(summary_df, row.names = FALSE)

cat("\n===== Trimmed gene-level summary by signature direction =====\n")
print(direction_summary, row.names = FALSE)

cat("\n===== Main 24 core genes not preserved as core-like in trimmed branch =====\n")
if (nrow(changed_main_core) == 0) {
  cat("All main 24 core genes were preserved as core-like in the trimmed branch.\n")
} else {
  print(changed_main_core[, c(
    "pig_ensg", "pig_symbol", "human_gene", "signature_direction",
    "direction_consistent_both_trimmed",
    "stat_strict_t7_trimmed", "strict_t7_trimmed",
    "stat_strict_t28_trimmed", "strict_t28_trimmed",
    "core_like_trimmed",
    "logFC_t7_trimmed", "FDR_t7_trimmed", "logFC_t28_trimmed", "FDR_t28_trimmed"
  ), drop = FALSE], row.names = FALSE)
}

cat("\n===== Version and method records =====\n")
print(version_df, row.names = FALSE)

## concise summary log
summary_con <- file(summary_log_file, open = "wt")
writeLines("===== StepH2w10_current78 core gene sensitivity audit SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(summary_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Trimmed gene-level summary by signature direction:", summary_con)
writeLines(capture.output(print(direction_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Main 24 core genes not preserved as core-like in trimmed branch:", summary_con)
if (nrow(changed_main_core) == 0) {
  writeLines("All main 24 core genes were preserved as core-like in the trimmed branch.", summary_con)
} else {
  writeLines(capture.output(print(changed_main_core[, c(
    "pig_ensg", "pig_symbol", "human_gene", "signature_direction",
    "direction_consistent_both_trimmed",
    "stat_strict_t7_trimmed", "strict_t7_trimmed",
    "stat_strict_t28_trimmed", "strict_t28_trimmed",
    "core_like_trimmed",
    "logFC_t7_trimmed", "FDR_t7_trimmed", "logFC_t28_trimmed", "FDR_t28_trimmed"
  ), drop = FALSE], row.names = FALSE)), summary_con)
}
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("The official main core remains the untrimmed primary-analysis 24-gene core set.", summary_con)
writeLines("The trimmed core-like count is only a sensitivity audit result, not a replacement core definition.", summary_con)
writeLines("In this corrected version, strict_t7_trimmed/strict_t28_trimmed include both statistical threshold and direction consistency, matching the primary-analysis manuscript definition.", summary_con)
writeLines("The stat_strict_* columns record FDR < 0.05 and |logFC| > 1 without direction consistency for auditing only.", summary_con)
writeLines("Please verify whether all or most main 24 core genes are preserved as direction-consistent and strict in the trimmed branch.", summary_con)
close(summary_con)

sink()

cat("\nStepH2w10_current78 core gene sensitivity audit v2 completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
