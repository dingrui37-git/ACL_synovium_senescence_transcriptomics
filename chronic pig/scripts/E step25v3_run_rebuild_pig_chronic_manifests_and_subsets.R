# ============================================================
# Step 25 v3: rebuild pig-chronic group labels from synovium sample_id
# Project root: E:/R/ACLsenescence2
#
# Purpose:
#   - reuse Step 25 v2 synovium-only manifest (48 samples)
#   - derive chronic group labels directly from sample_id pattern
#   - build main comparison manifest:
#       Control_52W vs ACLT_alone_52W  (12 vs 12)
#   - save subset matrices for downstream chronic DE / GSEA / Figure 5
# ============================================================

setwd("E:/R/ACLsenescence2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rebuild_root <- file.path(project_root, "rebuild_submission")
chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

scripts_dir <- file.path(chronic_dir, "scripts")
objects_dir <- file.path(chronic_dir, "objects")
tables_dir <- file.path(chronic_dir, "tables")
logs_dir <- file.path(chronic_dir, "logs")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

syn_manifest_file <- file.path(tables_dir, "step25v2_pig_chronic_synovium_manifest.csv")
counts_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_estimated_counts_matrix.csv")
abundance_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_abundance_tpm_matrix.csv")
length_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_length_matrix.csv")

needed <- c(syn_manifest_file, counts_file, abundance_file, length_file)
if (!all(file.exists(needed))) {
  stop("Required files are missing. Please make sure Step 24 and Step 25 v2 finished successfully.")
}

syn_manifest_df <- read.csv(syn_manifest_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
counts_df <- read.csv(counts_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)
abundance_df <- read.csv(abundance_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)
length_df <- read.csv(length_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)

if (!("sample_id" %in% colnames(syn_manifest_df))) stop("Synovium manifest must contain sample_id.")
if (!("gene_id" %in% colnames(counts_df))) stop("Counts matrix must contain gene_id.")

derive_core_group_from_sample_id <- function(x) {
  x_up <- toupper(x)
  out <- ifelse(grepl("_CON[0-9]+_", x_up), "Control_52W",
         ifelse(grepl("_ACLT[0-9]+_", x_up), "ACLT_alone_52W",
         ifelse(grepl("_RECON[0-9]+_", x_up), "Reconstruction_52W",
         ifelse(grepl("_REPAIR[0-9]+_", x_up), "Repair_52W", NA_character_))))
  out
}

derive_time_label_from_sample_id <- function(x) {
  rep("52W", length(x))
}

syn_manifest_df$core_group <- derive_core_group_from_sample_id(syn_manifest_df$sample_id)
syn_manifest_df$time_label <- derive_time_label_from_sample_id(syn_manifest_df$sample_id)

if (all(is.na(syn_manifest_df$core_group))) {
  stop("Failed to derive chronic group labels from sample_id.")
}

group_summary_df <- as.data.frame(table(syn_manifest_df$core_group), stringsAsFactors = FALSE)
colnames(group_summary_df) <- c("core_group", "n_samples")
group_summary_df <- group_summary_df[group_summary_df$core_group != "", , drop = FALSE]

main_manifest_df <- syn_manifest_df[syn_manifest_df$core_group %in% c("Control_52W", "ACLT_alone_52W"), , drop = FALSE]

subset_matrix_df <- function(mat_df, keep_samples) {
  keep_samples <- intersect(keep_samples, colnames(mat_df))
  mat_df[, c("gene_id", keep_samples), drop = FALSE]
}

syn_counts_df <- subset_matrix_df(counts_df, syn_manifest_df$sample_id)
syn_abundance_df <- subset_matrix_df(abundance_df, syn_manifest_df$sample_id)
syn_length_df <- subset_matrix_df(length_df, syn_manifest_df$sample_id)

main_counts_df <- subset_matrix_df(counts_df, main_manifest_df$sample_id)
main_abundance_df <- subset_matrix_df(abundance_df, main_manifest_df$sample_id)
main_length_df <- subset_matrix_df(length_df, main_manifest_df$sample_id)

summary_df <- data.frame(
  metric = c(
    "project_root",
    "source_synovium_manifest",
    "n_synovium_samples",
    "n_main_comparison_samples",
    "n_main_control_52w",
    "n_main_aclt_alone_52w",
    "synovium_counts_dim",
    "main_counts_dim",
    "status"
  ),
  value = c(
    project_root,
    "step25v2_pig_chronic_synovium_manifest.csv",
    nrow(syn_manifest_df),
    nrow(main_manifest_df),
    sum(main_manifest_df$core_group == "Control_52W", na.rm = TRUE),
    sum(main_manifest_df$core_group == "ACLT_alone_52W", na.rm = TRUE),
    paste(nrow(syn_counts_df), "x", ncol(syn_counts_df) - 1),
    paste(nrow(main_counts_df), "x", ncol(main_counts_df) - 1),
    "Step 25 v3 completed successfully"
  ),
  stringsAsFactors = FALSE
)

write.csv(syn_manifest_df,
          file.path(tables_dir, "step25v3_pig_chronic_synovium_manifest.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(group_summary_df,
          file.path(tables_dir, "step25v3_pig_chronic_group_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(main_manifest_df,
          file.path(tables_dir, "step25v3_pig_chronic_main_comparison_manifest.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(syn_counts_df,
          file.path(tables_dir, "step25v3_pig_chronic_synovium_gene_level_counts_matrix.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(syn_abundance_df,
          file.path(tables_dir, "step25v3_pig_chronic_synovium_gene_level_abundance_tpm_matrix.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(syn_length_df,
          file.path(tables_dir, "step25v3_pig_chronic_synovium_gene_level_length_matrix.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(main_counts_df,
          file.path(tables_dir, "step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(main_abundance_df,
          file.path(tables_dir, "step25v3_pig_chronic_main_comparison_gene_level_abundance_tpm_matrix.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(main_length_df,
          file.path(tables_dir, "step25v3_pig_chronic_main_comparison_gene_level_length_matrix.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(summary_df,
          file.path(tables_dir, "step25v3_pig_chronic_manifest_rebuild_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

save(
  syn_manifest_df, group_summary_df, main_manifest_df,
  syn_counts_df, syn_abundance_df, syn_length_df,
  main_counts_df, main_abundance_df, main_length_df,
  summary_df,
  file = file.path(objects_dir, "step25v3_pig_chronic_manifest_rebuild_workspace.RData")
)

writeLines(
  c(
    "Step 25 v3 completed successfully.",
    paste("Synovium samples:", nrow(syn_manifest_df)),
    paste("Main comparison samples:", nrow(main_manifest_df)),
    paste("Control_52W:", sum(main_manifest_df$core_group == 'Control_52W', na.rm = TRUE)),
    paste("ACLT_alone_52W:", sum(main_manifest_df$core_group == 'ACLT_alone_52W', na.rm = TRUE))
  ),
  file.path(logs_dir, "step25v3_pig_chronic_manifest_rebuild_log.txt"),
  useBytes = TRUE
)

message("Step 25 v3 finished successfully.")
