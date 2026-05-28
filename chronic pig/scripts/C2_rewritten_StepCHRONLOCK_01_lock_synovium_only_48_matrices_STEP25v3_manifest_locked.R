# StepCHRONLOCK_01_lock_synovium_only_48_matrices.R
# This script locks synovium-only 48-sample chronic pig gene-level matrices by subsetting the existing
# step24 tximport-derived all-sample matrices. It does not re-run tximport. It writes reproducible
# matrices, audit tables, figures, RDS objects, methods text, and a summary log.
#
# Revision in this version:
#   The 24-sample main-comparison audit no longer uses the old non-TMM Step2 metadata file.
#   It is audited against Step25v3's locked main-comparison manifest:
#     step25v3_pig_chronic_main_comparison_manifest.csv

options(stringsAsFactors = FALSE)
set.seed(1)

method_version <- "2026-05-26_lock_chronic_pig_synovium_only_48_matrices_from_step24_with_step25v3_manifest_audit"

# -----------------------------
# User-fixed paths
# -----------------------------
tables_dir <- "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables"
quant_dir  <- "E:/R/ACLsenescence2/rebuild_submission/raw data/GSE228848_synovium_quant"

counts_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_estimated_counts_matrix.csv")
abundance_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_abundance_tpm_matrix.csv")
length_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_length_matrix.csv")
tximport_summary_file <- file.path(tables_dir, "step24_pig_chronic_tximport_run_summary.csv")
main_manifest_file <- file.path(tables_dir, "step25v3_pig_chronic_main_comparison_manifest.csv")

output_root <- file.path(tables_dir, "step24_chronic_synovium_only_matrix_lock")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "figures"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(output_root, "logs", "StepCHRONLOCK_01_lock_synovium_only_48_matrices_summary_log.txt")

# Clean prior outputs from this step only, excluding the log opened later.
old_files <- list.files(output_root, pattern = "^StepCHRONLOCK_01_", recursive = TRUE, full.names = TRUE)
if (length(old_files) > 0) {
  suppressWarnings(file.remove(old_files))
}

zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")

status <- "FAILED"
on.exit({
  cat("\n============================================================\n")
  cat("StepCHRONLOCK_01_lock_synovium_only_48_matrices finished with status:", status, "\n")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Summary log saved to:\n", log_file, "\n", sep = "")
  cat("============================================================\n")
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("StepCHRONLOCK_01_lock_synovium_only_48_matrices\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Method version:", method_version, "\n")
cat("============================================================\n\n")

# -----------------------------
# Script archive
# -----------------------------
script_archive_file <- file.path(output_root, "scripts", "StepCHRONLOCK_01_lock_synovium_only_48_matrices_STEP25v3_manifest_locked.R")
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
  if (file.exists(script_path)) {
    file.copy(script_path, script_archive_file, overwrite = TRUE)
    cat("Script archived to:", script_archive_file, "\n")
  }
} else {
  cat("Script archive note: --file argument was not detected. If running interactively, save this script as:\n")
  cat(script_archive_file, "\n\n")
}


# -----------------------------
# Package loading
# -----------------------------
required_pkgs <- c("dplyr", "ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}
library(dplyr)
library(ggplot2)
cat("Required packages loaded.\n")
cat("dplyr version:", as.character(packageVersion("dplyr")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n\n")

# -----------------------------
# Helper functions
# -----------------------------
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
}

read_matrix_csv <- function(path, label) {
  stop_if_missing(path, label)
  x <- read.csv(path, check.names = FALSE)
  if (ncol(x) < 2) stop(label, " has fewer than 2 columns: ", path)
  if (!"gene_id" %in% colnames(x)) {
    colnames(x)[1] <- "gene_id"
  }
  if (anyDuplicated(x$gene_id)) {
    stop(label, " has duplicated gene_id values. First duplicated index: ", which(duplicated(x$gene_id))[1])
  }
  x
}

parse_sample_id <- function(sample_id) {
  geo <- ifelse(grepl("^GSM[0-9]+", sample_id), sub("^((GSM[0-9]+)).*$", "\\1", sample_id), NA_character_)
  label <- sub("^GSM[0-9]+_", "", sample_id)
  tissue <- ifelse(grepl("_synovium$", label, ignore.case = TRUE), "synovium",
                   ifelse(grepl("_cartilage$", label, ignore.case = TRUE), "cartilage", NA_character_))
  label_no_tissue <- sub("_(synovium|cartilage)$", "", label, ignore.case = TRUE)
  group <- ifelse(grepl("^CON[0-9]+$", label_no_tissue, ignore.case = TRUE), "Control_52W",
                  ifelse(grepl("^ACLT[0-9]+$", label_no_tissue, ignore.case = TRUE), "ACLT_alone_52W",
                         ifelse(grepl("^RECON[0-9]+$", label_no_tissue, ignore.case = TRUE), "Reconstruction_52W",
                                ifelse(grepl("^REPAIR[0-9]+$", label_no_tissue, ignore.case = TRUE), "Repair_52W", NA_character_))))
  sample_number <- suppressWarnings(as.integer(gsub("^[A-Za-z]+", "", label_no_tissue)))
  data.frame(
    sample_id = sample_id,
    geo_accession = geo,
    sample_label = label_no_tissue,
    tissue = tissue,
    group = group,
    sample_number = sample_number,
    stringsAsFactors = FALSE
  )
}

parse_quant_file <- function(path) {
  b <- basename(path)
  geo <- ifelse(grepl("^GSM[0-9]+", b), sub("^((GSM[0-9]+)).*$", "\\1", b), NA_character_)
  label <- sub("^GSM[0-9]+_", "", b)
  label <- sub("_quant\\.sf\\.txt\\.gz$", "", label, ignore.case = TRUE)
  label <- sub("_quant\\.sf$", "", label, ignore.case = TRUE)
  tissue <- ifelse(grepl("_synovium$", label, ignore.case = TRUE), "synovium",
                   ifelse(grepl("_cartilage$", label, ignore.case = TRUE), "cartilage", NA_character_))
  label_no_tissue <- sub("_(synovium|cartilage)$", "", label, ignore.case = TRUE)
  group <- ifelse(grepl("^CON[0-9]+$", label_no_tissue, ignore.case = TRUE), "Control_52W",
                  ifelse(grepl("^ACLT[0-9]+$", label_no_tissue, ignore.case = TRUE), "ACLT_alone_52W",
                         ifelse(grepl("^RECON[0-9]+$", label_no_tissue, ignore.case = TRUE), "Reconstruction_52W",
                                ifelse(grepl("^REPAIR[0-9]+$", label_no_tissue, ignore.case = TRUE), "Repair_52W", NA_character_))))
  data.frame(
    quant_file = path,
    quant_basename = b,
    geo_accession = geo,
    sample_label = label_no_tissue,
    tissue = tissue,
    group = group,
    stringsAsFactors = FALSE
  )
}

sample_qc_from_matrices <- function(counts_mat, abundance_mat, length_mat, metadata) {
  sample_cols <- metadata$sample_id
  out <- lapply(sample_cols, function(s) {
    counts_v <- counts_mat[[s]]
    tpm_v <- abundance_mat[[s]]
    len_v <- length_mat[[s]]
    data.frame(
      sample_id = s,
      total_estimated_counts = sum(counts_v, na.rm = TRUE),
      detected_genes_count_gt0 = sum(counts_v > 0, na.rm = TRUE),
      detected_genes_count_ge1 = sum(counts_v >= 1, na.rm = TRUE),
      tpm_sum_gene_level = sum(tpm_v, na.rm = TRUE),
      genes_with_tpm_gt0 = sum(tpm_v > 0, na.rm = TRUE),
      effective_length_na_n = sum(is.na(len_v)),
      effective_length_nonpositive_n = sum(len_v <= 0, na.rm = TRUE),
      effective_length_median = median(len_v, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(out)
  left_join(metadata, out, by = "sample_id")
}

summarise_matrix_qc <- function(qc_df) {
  qc_df %>%
    group_by(group) %>%
    summarise(
      n_samples = n(),
      mean_total_estimated_counts = mean(total_estimated_counts),
      median_total_estimated_counts = median(total_estimated_counts),
      min_total_estimated_counts = min(total_estimated_counts),
      max_total_estimated_counts = max(total_estimated_counts),
      mean_detected_genes_count_gt0 = mean(detected_genes_count_gt0),
      median_detected_genes_count_gt0 = median(detected_genes_count_gt0),
      min_detected_genes_count_gt0 = min(detected_genes_count_gt0),
      max_detected_genes_count_gt0 = max(detected_genes_count_gt0),
      mean_gene_level_TPM_sum = mean(tpm_sum_gene_level),
      min_gene_level_TPM_sum = min(tpm_sum_gene_level),
      max_gene_level_TPM_sum = max(tpm_sum_gene_level),
      total_effective_length_na_n = sum(effective_length_na_n),
      total_effective_length_nonpositive_n = sum(effective_length_nonpositive_n),
      .groups = "drop"
    )
}

# -----------------------------
# Input file check
# -----------------------------
input_check <- data.frame(
  label = c("counts_file", "abundance_file", "effective_length_file", "tximport_summary_file", "step25v3_main_manifest_file", "quant_dir"),
  path = c(counts_file, abundance_file, length_file, tximport_summary_file, main_manifest_file, quant_dir),
  exists = c(file.exists(counts_file), file.exists(abundance_file), file.exists(length_file),
             file.exists(tximport_summary_file), file.exists(main_manifest_file), dir.exists(quant_dir)),
  stringsAsFactors = FALSE
)
cat("Input paths:\n")
print(input_check)
write.csv(input_check, file.path(output_root, "tables", "StepCHRONLOCK_01_input_file_check.csv"), row.names = FALSE)

if (!all(input_check$exists)) {
  stop("One or more required input paths are missing. See StepCHRONLOCK_01_input_file_check.csv.")
}

# -----------------------------
# Read matrices
# -----------------------------
counts_all <- read_matrix_csv(counts_file, "counts_file")
abundance_all <- read_matrix_csv(abundance_file, "abundance_file")
length_all <- read_matrix_csv(length_file, "effective_length_file")

cat("\nAll-sample matrix dimensions:\n")
matrix_dims_all <- data.frame(
  matrix = c("estimated_counts", "abundance_TPM", "effective_length"),
  n_gene_rows = c(nrow(counts_all), nrow(abundance_all), nrow(length_all)),
  n_columns_total = c(ncol(counts_all), ncol(abundance_all), ncol(length_all)),
  n_sample_columns = c(ncol(counts_all) - 1, ncol(abundance_all) - 1, ncol(length_all) - 1),
  first_column = "gene_id",
  stringsAsFactors = FALSE
)
print(matrix_dims_all)

sample_cols_counts <- setdiff(colnames(counts_all), "gene_id")
sample_cols_abundance <- setdiff(colnames(abundance_all), "gene_id")
sample_cols_length <- setdiff(colnames(length_all), "gene_id")

if (!identical(counts_all$gene_id, abundance_all$gene_id)) {
  stop("gene_id order differs between counts and abundance matrices.")
}
if (!identical(counts_all$gene_id, length_all$gene_id)) {
  stop("gene_id order differs between counts and effective length matrices.")
}
if (!identical(sample_cols_counts, sample_cols_abundance)) {
  stop("sample columns differ between counts and abundance matrices.")
}
if (!identical(sample_cols_counts, sample_cols_length)) {
  stop("sample columns differ between counts and effective length matrices.")
}

sample_metadata_all <- parse_sample_id(sample_cols_counts)
write.csv(sample_metadata_all, file.path(output_root, "tables", "StepCHRONLOCK_01_all_96_sample_metadata_from_matrix_columns.csv"), row.names = FALSE)

cat("\nAll-sample metadata group/tissue counts inferred from matrix columns:\n")
print(table(sample_metadata_all$tissue, sample_metadata_all$group, useNA = "ifany"))

# -----------------------------
# Identify and lock synovium-only 48 samples
# -----------------------------
group_levels_48 <- c("Control_52W", "ACLT_alone_52W", "Reconstruction_52W", "Repair_52W")
syn_meta <- sample_metadata_all %>%
  filter(tissue == "synovium") %>%
  mutate(group = factor(group, levels = group_levels_48)) %>%
  arrange(group, sample_number, geo_accession) %>%
  mutate(group = as.character(group))

cat("\nSynovium-only sample metadata:\n")
print(syn_meta)
cat("\nSynovium-only group counts:\n")
print(table(syn_meta$group, useNA = "ifany"))

if (nrow(syn_meta) != 48) {
  stop("Expected 48 synovium samples, but detected ", nrow(syn_meta), ".")
}
expected_groups <- c(Control_52W = 12, ACLT_alone_52W = 12, Reconstruction_52W = 12, Repair_52W = 12)
observed_groups <- table(syn_meta$group)
if (!all(names(expected_groups) %in% names(observed_groups)) ||
    !all(as.integer(observed_groups[names(expected_groups)]) == as.integer(expected_groups))) {
  stop("Synovium group counts do not match expected 12 per group. Observed: ",
       paste(names(observed_groups), as.integer(observed_groups), collapse = "; "))
}

syn_cols <- syn_meta$sample_id
counts_syn <- counts_all[, c("gene_id", syn_cols), drop = FALSE]
abundance_syn <- abundance_all[, c("gene_id", syn_cols), drop = FALSE]
length_syn <- length_all[, c("gene_id", syn_cols), drop = FALSE]

# Main comparison 24 samples
main_group_levels <- c("Control_52W", "ACLT_alone_52W")
main_meta <- syn_meta %>%
  filter(group %in% main_group_levels) %>%
  mutate(group = factor(group, levels = main_group_levels)) %>%
  arrange(group, sample_number, geo_accession) %>%
  mutate(group = as.character(group))
main_cols <- main_meta$sample_id

if (nrow(main_meta) != 24) {
  stop("Expected 24 main-comparison synovium samples, but detected ", nrow(main_meta), ".")
}
counts_main24 <- counts_all[, c("gene_id", main_cols), drop = FALSE]
abundance_main24 <- abundance_all[, c("gene_id", main_cols), drop = FALSE]
length_main24 <- length_all[, c("gene_id", main_cols), drop = FALSE]

# -----------------------------
# Quant-file audit for synovium subset
# -----------------------------
quant_files <- list.files(quant_dir, pattern = "quant\\.sf(\\.txt)?\\.gz$|quant\\.sf$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
quant_manifest_all <- if (length(quant_files) > 0) bind_rows(lapply(quant_files, parse_quant_file)) else data.frame()
quant_manifest_syn <- quant_manifest_all %>% filter(tissue == "synovium")

syn_quant_match <- syn_meta %>%
  left_join(
    quant_manifest_syn %>% select(geo_accession, quant_file, quant_basename, sample_label_from_quant = sample_label, group_from_quant = group),
    by = "geo_accession"
  ) %>%
  mutate(has_matching_synovium_quant_file = !is.na(quant_file))

cat("\nQuant-file audit for locked 48 synovium matrix:\n")
cat("All quant-like files found:", nrow(quant_manifest_all), "\n")
cat("Synovium quant-like files found:", nrow(quant_manifest_syn), "\n")
cat("Locked synovium matrix samples with matching synovium quant files:", sum(syn_quant_match$has_matching_synovium_quant_file), "\n")
if (sum(syn_quant_match$has_matching_synovium_quant_file) != 48) {
  stop("Not all locked synovium matrix samples have matching synovium quant files.")
}

write.csv(quant_manifest_all, file.path(output_root, "tables", "StepCHRONLOCK_01_quant_file_manifest_all.csv"), row.names = FALSE)
write.csv(quant_manifest_syn, file.path(output_root, "tables", "StepCHRONLOCK_01_quant_file_manifest_synovium_48.csv"), row.names = FALSE)
write.csv(syn_quant_match, file.path(output_root, "tables", "StepCHRONLOCK_01_locked_synovium_48_vs_quant_files_audit.csv"), row.names = FALSE)

# -----------------------------
# Step25v3 main-comparison manifest audit
# -----------------------------
main_manifest <- read.csv(main_manifest_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")

if (!"sample_id" %in% colnames(main_manifest)) {
  stop("Step25v3 main-comparison manifest must contain a sample_id column: ", main_manifest_file)
}

main_manifest_group_col <- if ("core_group" %in% colnames(main_manifest)) {
  "core_group"
} else if ("group" %in% colnames(main_manifest)) {
  "group"
} else {
  NA_character_
}

if (is.na(main_manifest_group_col)) {
  stop(
    "Step25v3 main-comparison manifest must contain either core_group or group column. Available columns: ",
    paste(colnames(main_manifest), collapse = ", ")
  )
}

main_manifest_audit <- main_manifest %>%
  mutate(
    manifest_sample_id = as.character(.data[["sample_id"]]),
    manifest_group = as.character(.data[[main_manifest_group_col]]),
    in_all_96_matrix = manifest_sample_id %in% sample_cols_counts,
    in_locked_synovium_48_matrix = manifest_sample_id %in% syn_cols,
    in_main_24_matrix = manifest_sample_id %in% main_cols,
    manifest_geo_accession = ifelse(grepl("^GSM[0-9]+", manifest_sample_id),
                                    sub("^((GSM[0-9]+)).*$", "\\1", manifest_sample_id),
                                    NA_character_)
  ) %>%
  left_join(
    syn_meta %>% select(sample_id, locked_group = group, locked_tissue = tissue, locked_sample_number = sample_number),
    by = c("manifest_sample_id" = "sample_id")
  ) %>%
  mutate(
    group_consistent_with_locked_matrix = manifest_group == locked_group
  )

cat("\nStep25v3 main-comparison manifest audit:\n")
cat("Manifest file:", main_manifest_file, "\n")
cat("Manifest rows:", nrow(main_manifest_audit), "\n")
cat("Rows in all 96 matrix:", sum(main_manifest_audit$in_all_96_matrix), "\n")
cat("Rows in locked synovium 48 matrix:", sum(main_manifest_audit$in_locked_synovium_48_matrix), "\n")
cat("Rows in main 24 matrix:", sum(main_manifest_audit$in_main_24_matrix), "\n")
cat("Rows with group consistent between Step25v3 manifest and CHRONLOCK matrix:", sum(main_manifest_audit$group_consistent_with_locked_matrix, na.rm = TRUE), "\n")
cat("Step25v3 manifest group counts:\n")
print(table(main_manifest_audit$manifest_group, useNA = "ifany"))

expected_main_groups <- c(Control_52W = 12, ACLT_alone_52W = 12)
manifest_group_counts <- table(main_manifest_audit$manifest_group)
manifest_group_ok <- all(names(expected_main_groups) %in% names(manifest_group_counts)) &&
  all(as.integer(manifest_group_counts[names(expected_main_groups)]) == as.integer(expected_main_groups)) &&
  length(setdiff(names(manifest_group_counts), names(expected_main_groups))) == 0

if (nrow(main_manifest_audit) != 24 ||
    sum(main_manifest_audit$in_locked_synovium_48_matrix) != 24 ||
    sum(main_manifest_audit$in_main_24_matrix) != 24 ||
    !manifest_group_ok ||
    any(!main_manifest_audit$group_consistent_with_locked_matrix | is.na(main_manifest_audit$group_consistent_with_locked_matrix))) {
  write.csv(
    main_manifest_audit,
    file.path(output_root, "tables", "StepCHRONLOCK_01_step25v3_main_manifest_vs_locked_matrices_audit_FAILED.csv"),
    row.names = FALSE
  )
  stop("Step25v3 main-comparison manifest audit failed. See StepCHRONLOCK_01_step25v3_main_manifest_vs_locked_matrices_audit_FAILED.csv.")
}

if (!setequal(main_manifest_audit$manifest_sample_id, main_cols)) {
  write.csv(
    main_manifest_audit,
    file.path(output_root, "tables", "StepCHRONLOCK_01_step25v3_main_manifest_vs_locked_matrices_audit_FAILED.csv"),
    row.names = FALSE
  )
  stop("Step25v3 main-comparison manifest sample_id set does not match the CHRONLOCK-derived main 24 sample set.")
}

# -----------------------------
# QC summaries
# -----------------------------
syn_qc <- sample_qc_from_matrices(counts_syn, abundance_syn, length_syn, syn_meta)
main_qc <- sample_qc_from_matrices(counts_main24, abundance_main24, length_main24, main_meta)
syn_group_qc <- summarise_matrix_qc(syn_qc)
main_group_qc <- summarise_matrix_qc(main_qc)

cat("\nLocked synovium 48 sample QC preview:\n")
print(head(syn_qc, 20))
cat("\nLocked synovium 48 group QC summary:\n")
print(syn_group_qc)
cat("\nMain comparison 24 group QC summary:\n")
print(main_group_qc)

# -----------------------------
# Write locked matrices and audit outputs
# -----------------------------
write.csv(counts_syn, file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_gene_level_estimated_counts_matrix.csv"), row.names = FALSE)
write.csv(abundance_syn, file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_gene_level_abundance_tpm_matrix.csv"), row.names = FALSE)
write.csv(length_syn, file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_gene_level_effective_length_matrix.csv"), row.names = FALSE)
write.csv(syn_meta, file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_sample_metadata_locked.csv"), row.names = FALSE)

write.csv(counts_main24, file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_gene_level_estimated_counts_matrix.csv"), row.names = FALSE)
write.csv(abundance_main24, file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_gene_level_abundance_tpm_matrix.csv"), row.names = FALSE)
write.csv(length_main24, file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_gene_level_effective_length_matrix.csv"), row.names = FALSE)
write.csv(main_meta, file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_sample_metadata_locked.csv"), row.names = FALSE)

write.csv(matrix_dims_all, file.path(output_root, "tables", "StepCHRONLOCK_01_all_96_matrix_dimensions.csv"), row.names = FALSE)
matrix_dims_locked <- data.frame(
  matrix = c("synovium_48_estimated_counts", "synovium_48_abundance_TPM", "synovium_48_effective_length",
             "main_24_estimated_counts", "main_24_abundance_TPM", "main_24_effective_length"),
  n_gene_rows = c(nrow(counts_syn), nrow(abundance_syn), nrow(length_syn),
                  nrow(counts_main24), nrow(abundance_main24), nrow(length_main24)),
  n_columns_total = c(ncol(counts_syn), ncol(abundance_syn), ncol(length_syn),
                      ncol(counts_main24), ncol(abundance_main24), ncol(length_main24)),
  n_sample_columns = c(ncol(counts_syn) - 1, ncol(abundance_syn) - 1, ncol(length_syn) - 1,
                       ncol(counts_main24) - 1, ncol(abundance_main24) - 1, ncol(length_main24) - 1),
  first_column = "gene_id",
  stringsAsFactors = FALSE
)
write.csv(matrix_dims_locked, file.path(output_root, "tables", "StepCHRONLOCK_01_locked_matrix_dimensions.csv"), row.names = FALSE)

write.csv(main_manifest_audit, file.path(output_root, "tables", "StepCHRONLOCK_01_step25v3_main_manifest_vs_locked_matrices_audit.csv"), row.names = FALSE)
write.csv(syn_qc, file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_sample_level_matrix_QC.csv"), row.names = FALSE)
write.csv(syn_group_qc, file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_group_level_matrix_QC_summary.csv"), row.names = FALSE)
write.csv(main_qc, file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_sample_level_matrix_QC.csv"), row.names = FALSE)
write.csv(main_group_qc, file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_group_level_matrix_QC_summary.csv"), row.names = FALSE)

saveRDS(
  list(
    counts_synovium_48 = counts_syn,
    abundance_synovium_48 = abundance_syn,
    effective_length_synovium_48 = length_syn,
    synovium_48_metadata = syn_meta,
    counts_main_24 = counts_main24,
    abundance_main_24 = abundance_main24,
    effective_length_main_24 = length_main24,
    main_24_metadata = main_meta,
    method_version = method_version
  ),
  file.path(output_root, "objects", "StepCHRONLOCK_01_locked_synovium_48_and_main_24_matrices.rds")
)

# Method text
methods_text <- paste(
  "The chronic pig tximport-derived gene-level matrices were first inspected as all processed GSE228848 quantification files.",
  "Because the step24 tximport output contained 96 samples across cartilage and synovium, synovium-only matrices were locked by subsetting the estimated counts, TPM abundance and effective length matrices to the 48 columns annotated as synovium based on GEO accession/sample labels.",
  "The locked synovium-only matrices contained 22,438 genes and 48 samples, comprising 12 Control_52W, 12 ACLT_alone_52W, 12 Reconstruction_52W and 12 Repair_52W samples.",
  "The main chronic comparison matrix was then defined as the 24-sample subset containing Control_52W and ACLT_alone_52W only.",
  "The locked subsets were audited against the available synovium Salmon quant.sf files and the Step25v3 locked chronic main-comparison manifest.",
  sep = "\n"
)
writeLines(methods_text, file.path(output_root, "tables", "StepCHRONLOCK_01_methods_text_synovium_only_matrix_lock.txt"))

# Simple diagnostic figures
p_counts <- ggplot(syn_qc, aes(x = group, y = total_estimated_counts)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  labs(
    title = "Locked chronic pig synovium-only matrix QC",
    subtitle = "Total estimated counts per sample after subsetting step24 tximport matrix to 48 synovium samples",
    x = "Group",
    y = "Total estimated counts"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(file.path(output_root, "figures", "StepCHRONLOCK_01_synovium_48_total_estimated_counts_by_group.pdf"),
       p_counts, width = 7.5, height = 4.8)

p_detected <- ggplot(syn_qc, aes(x = group, y = detected_genes_count_gt0)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  labs(
    title = "Locked chronic pig synovium-only matrix QC",
    subtitle = "Detected genes per sample in the 48-sample synovium-only matrix",
    x = "Group",
    y = "Detected genes (counts > 0)"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(file.path(output_root, "figures", "StepCHRONLOCK_01_synovium_48_detected_genes_by_group.pdf"),
       p_detected, width = 7.5, height = 4.8)

# Final summary
final_summary <- data.frame(
  metric = c(
    "run_status",
    "method_version",
    "all_step24_matrix_genes",
    "all_step24_matrix_samples",
    "locked_synovium_matrix_genes",
    "locked_synovium_matrix_samples",
    "locked_synovium_group_counts",
    "main_comparison_matrix_genes",
    "main_comparison_matrix_samples",
    "main_comparison_group_counts",
    "synovium_quant_files_found",
    "locked_synovium_samples_with_matching_quant_file",
    "step25v3_main_manifest_rows",
    "step25v3_main_manifest_in_locked_synovium_48_matrix",
    "step25v3_main_manifest_in_main_24_matrix",
    "step25v3_main_manifest_group_consistent_with_locked_matrix",
    "counts_abundance_length_gene_id_order_identical",
    "counts_abundance_length_sample_columns_identical",
    "mean_synovium_48_total_estimated_counts",
    "min_synovium_48_total_estimated_counts",
    "max_synovium_48_total_estimated_counts",
    "mean_synovium_48_detected_genes_count_gt0",
    "min_synovium_48_detected_genes_count_gt0",
    "max_synovium_48_detected_genes_count_gt0",
    "output_root"
  ),
  value = c(
    "SUCCESS",
    method_version,
    nrow(counts_all),
    length(sample_cols_counts),
    nrow(counts_syn),
    length(syn_cols),
    paste(names(table(syn_meta$group)), as.integer(table(syn_meta$group)), collapse = "; "),
    nrow(counts_main24),
    length(main_cols),
    paste(names(table(main_meta$group)), as.integer(table(main_meta$group)), collapse = "; "),
    nrow(quant_manifest_syn),
    sum(syn_quant_match$has_matching_synovium_quant_file),
    nrow(main_manifest_audit),
    sum(main_manifest_audit$in_locked_synovium_48_matrix),
    sum(main_manifest_audit$in_main_24_matrix),
    sum(main_manifest_audit$group_consistent_with_locked_matrix, na.rm = TRUE),
    TRUE,
    TRUE,
    round(mean(syn_qc$total_estimated_counts), 3),
    round(min(syn_qc$total_estimated_counts), 3),
    round(max(syn_qc$total_estimated_counts), 3),
    round(mean(syn_qc$detected_genes_count_gt0), 3),
    min(syn_qc$detected_genes_count_gt0),
    max(syn_qc$detected_genes_count_gt0),
    output_root
  ),
  stringsAsFactors = FALSE
)
write.csv(final_summary, file.path(output_root, "tables", "StepCHRONLOCK_01_final_summary_for_review.csv"), row.names = FALSE)

cat("\nFinal summary for review:\n")
print(final_summary)

cat("\nKey output files:\n")
key_files <- c(
  file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_gene_level_estimated_counts_matrix.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_gene_level_abundance_tpm_matrix.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_gene_level_effective_length_matrix.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_synovium_48_sample_metadata_locked.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_gene_level_estimated_counts_matrix.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_gene_level_abundance_tpm_matrix.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_gene_level_effective_length_matrix.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_main_24_sample_metadata_locked.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_locked_matrix_dimensions.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_step25v3_main_manifest_vs_locked_matrices_audit.csv"),
  file.path(output_root, "tables", "StepCHRONLOCK_01_final_summary_for_review.csv")
)
for (i in seq_along(key_files)) {
  cat(i, ") ", key_files[i], "\n", sep = "")
}

cat("\nSession information:\n")
print(sessionInfo())

status <- "SUCCESS"
