# Step24B: Clean audit of chronic pig 96-sample tximport gene-level matrices
# Purpose:
#   This script audits the tximport-derived 96-sample gene-level matrices for GSE228848
#   after separate quant.sf-level QC has already been performed.
#   It intentionally does NOT repeat full per-file Salmon quant.sf integrity QC.
#   Instead, it checks:
#     1) Step24 tximport output files exist;
#     2) counts / TPM abundance / effective-length matrices have consistent dimensions,
#        gene_id order and sample-column order;
#     3) the 96 matrix sample columns are balanced as 48 synovium + 48 cartilage and
#        four treatment groups per tissue if sample IDs contain tissue/group labels;
#     4) sample-level gene-level QC metrics from tximport matrices;
#     5) optional matching of matrix samples to the Step23 quant-file inventory;
#     6) optional representative transcript-ID matching between one quant.sf file and tx2gene
#        to support the use of ignoreTxVersion = TRUE.
#
# Note:
#   Use this as the clean companion script for Methods section 12.3:
#   "tx2gene construction, 96-sample tximport summarization, and gene-level matrix consistency audit".
#   Use Step23 FIXED quant-level QC for Methods section 12.2.

method_version <- "2026-05-16_step24B_clean_tximport_gene_level_matrix_audit_v1"

# -----------------------------
# User-editable project paths
# -----------------------------
project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")
chronic_dir <- file.path(rebuild_root, "03_pig_chronic")
tables_dir <- file.path(chronic_dir, "tables")

output_root <- file.path(chronic_dir, "supplement", "Step24B_clean_tximport_gene_level_matrix_audit")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
for (d in c("logs", "tables", "source_data", "objects", "scripts")) {
  dir.create(file.path(output_root, d), recursive = TRUE, showWarnings = FALSE)
}

# Remove previous outputs from this clean audit step only
for (d in c("logs", "tables", "source_data", "objects", "scripts")) {
  old <- list.files(file.path(output_root, d), pattern = "^Step24B_", full.names = TRUE)
  if (length(old) > 0) suppressWarnings(file.remove(old))
}

# -----------------------------
# Logging
# -----------------------------
log_file <- file.path(output_root, "logs", "Step24B_clean_tximport_gene_level_matrix_audit_log.txt")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
run_status <- "FAILED"
on.exit({
  cat("\n============================================================\n")
  cat("Step24B clean tximport gene-level matrix audit finished with status: ", run_status, "\n", sep = "")
  cat("Finished at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("Summary log saved to:\n", log_file, "\n", sep = "")
  cat("============================================================\n")
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("Step24B_clean_tximport_gene_level_matrix_audit\n")
cat("Started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("Method version: ", method_version, "\n", sep = "")
cat("============================================================\n\n")

# Archive this script if launched by Rscript
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
  if (file.exists(script_path)) {
    file.copy(script_path, file.path(output_root, "scripts", basename(script_path)), overwrite = TRUE)
  }
} else {
  cat("Script archive note: --file argument was not detected. If running interactively, manually save this script in:\n",
      file.path(output_root, "scripts"), "\n\n", sep = "")
}

# -----------------------------
# Helper functions
# -----------------------------
read_csv_base <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

standardize_gene_matrix_df <- function(df, file_label) {
  if (ncol(df) < 2) stop(file_label, " must contain at least one gene_id column and one sample column.")
  first_col <- names(df)[1]
  if (!first_col %in% c("gene_id", "Geneid", "gene", "X", "")) {
    cat("Note: first column of ", file_label, " is named '", first_col,
        "'. Treating it as gene_id.\n", sep = "")
  }
  names(df)[1] <- "gene_id"
  if (any(is.na(df$gene_id) | df$gene_id == "")) {
    stop(file_label, " contains missing/empty gene_id values.")
  }
  df
}

to_numeric_matrix <- function(df, file_label) {
  mat <- as.matrix(df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "double")
  rownames(mat) <- df$gene_id
  if (any(!is.finite(mat) & !is.na(mat))) {
    warning(file_label, " contains non-finite numeric values.")
  }
  mat
}

extract_gsm <- function(x) {
  m <- regexpr("GSM[0-9]+", x)
  out <- regmatches(x, m)
  ifelse(m > 0, out, NA_character_)
}

strip_tx_version <- function(x) sub("\\.[0-9]+$", "", x)

clean_label <- function(x) {
  b <- basename(x)
  b <- sub("\\.gz$", "", b, ignore.case = TRUE)
  b <- sub("\\.txt$", "", b, ignore.case = TRUE)
  b <- sub("_quant\\.sf$", "", b, ignore.case = TRUE)
  b <- sub("\\.quant\\.sf$", "", b, ignore.case = TRUE)
  b <- sub("^GSM[0-9]+_", "", b)
  b
}

infer_tissue <- function(x) {
  bx <- tolower(basename(x))
  ifelse(grepl("synov", bx), "synovium",
         ifelse(grepl("cartilage", bx), "cartilage", NA_character_))
}

infer_group <- function(x) {
  z <- toupper(x)
  z <- sub("_(SYNOVIUM|CARTILAGE).*$", "", z)
  z <- sub("^GSM[0-9]+_", "", z)
  ifelse(grepl("^CON", z), "Control_52W",
         ifelse(grepl("^ACLT", z), "ACLT_alone_52W",
                ifelse(grepl("^RECON", z), "Reconstruction_52W",
                       ifelse(grepl("^REPAIR", z), "Repair_52W", NA_character_))))
}

find_optional_file <- function(dir, patterns) {
  if (!dir.exists(dir)) return(NA_character_)
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  if (length(files) == 0) return(NA_character_)
  keep <- rep(FALSE, length(files))
  for (p in patterns) keep <- keep | grepl(p, basename(files), ignore.case = TRUE)
  hits <- files[keep]
  if (length(hits) == 0) return(NA_character_)
  hits[order(file.info(hits)$mtime, decreasing = TRUE)][1]
}

find_col <- function(df, patterns) {
  nms <- names(df)
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

read_quant_name_column <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  q <- read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("Name" %in% names(q))) stop("Representative quant file does not contain Name column: ", path)
  q$Name
}

# -----------------------------
# Input files
# -----------------------------
counts_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_estimated_counts_matrix.csv")
abundance_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_abundance_tpm_matrix.csv")
length_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_length_matrix.csv")
gene_annot_file <- file.path(tables_dir, "step24_pig_chronic_gene_annotation_from_gtf.csv")
tx2gene_file <- file.path(tables_dir, "step24_pig_chronic_tx2gene.csv")
tximport_run_summary_file <- file.path(tables_dir, "step24_pig_chronic_tximport_run_summary.csv")

# Optional outputs from Step23 FIXED quant-level QC and Step25v3 subset locking
quant_inventory_file <- find_optional_file(
  tables_dir,
  c("^step23.*quant.*file.*inventory.*FIXED", "^step23.*quant.*inventory.*FIXED", "^step23.*quant.*file.*inventory")
)
step25v3_syn_manifest_file <- file.path(tables_dir, "step25v3_pig_chronic_synovium_manifest.csv")
step25v3_main_manifest_file <- file.path(tables_dir, "step25v3_pig_chronic_main_comparison_manifest.csv")

input_check <- data.frame(
  label = c(
    "counts_file",
    "abundance_file",
    "effective_length_file",
    "gene_annotation_file",
    "tx2gene_file",
    "tximport_run_summary_file",
    "optional_step23_quant_inventory",
    "optional_step25v3_synovium_manifest",
    "optional_step25v3_main_comparison_manifest"
  ),
  path = c(
    counts_file,
    abundance_file,
    length_file,
    gene_annot_file,
    tx2gene_file,
    tximport_run_summary_file,
    quant_inventory_file,
    step25v3_syn_manifest_file,
    step25v3_main_manifest_file
  ),
  exists = c(
    file.exists(counts_file),
    file.exists(abundance_file),
    file.exists(length_file),
    file.exists(gene_annot_file),
    file.exists(tx2gene_file),
    file.exists(tximport_run_summary_file),
    !is.na(quant_inventory_file) && file.exists(quant_inventory_file),
    file.exists(step25v3_syn_manifest_file),
    file.exists(step25v3_main_manifest_file)
  ),
  stringsAsFactors = FALSE
)

cat("Input file check:\n")
print(input_check)
write.csv(input_check, file.path(output_root, "tables", "Step24B_input_file_check.csv"), row.names = FALSE)

required <- c(counts_file, abundance_file, length_file)
if (!all(file.exists(required))) {
  stop("Required Step24 matrix files are missing. Please complete Step24 first.")
}

# -----------------------------
# Read Step24 matrices
# -----------------------------
counts_df <- standardize_gene_matrix_df(read_csv_base(counts_file), "counts_file")
abundance_df <- standardize_gene_matrix_df(read_csv_base(abundance_file), "abundance_file")
length_df <- standardize_gene_matrix_df(read_csv_base(length_file), "length_file")

count_mat <- to_numeric_matrix(counts_df, "counts_file")
abundance_mat <- to_numeric_matrix(abundance_df, "abundance_file")
length_mat <- to_numeric_matrix(length_df, "length_file")

sample_ids <- colnames(count_mat)

# -----------------------------
# Matrix dimension and consistency audit
# -----------------------------
matrix_dimensions <- data.frame(
  matrix = c("estimated_counts", "TPM_abundance", "effective_length"),
  file = c(counts_file, abundance_file, length_file),
  n_gene_rows = c(nrow(counts_df), nrow(abundance_df), nrow(length_df)),
  n_sample_columns = c(ncol(counts_df) - 1, ncol(abundance_df) - 1, ncol(length_df) - 1),
  first_column = c(names(counts_df)[1], names(abundance_df)[1], names(length_df)[1]),
  stringsAsFactors = FALSE
)

gene_id_consistency <- data.frame(
  comparison = c("counts_vs_abundance", "counts_vs_effective_length"),
  identical_gene_id_order = c(
    identical(counts_df$gene_id, abundance_df$gene_id),
    identical(counts_df$gene_id, length_df$gene_id)
  ),
  counts_n_genes = nrow(counts_df),
  other_n_genes = c(nrow(abundance_df), nrow(length_df)),
  stringsAsFactors = FALSE
)

sample_column_consistency <- data.frame(
  comparison = c("counts_vs_abundance", "counts_vs_effective_length"),
  identical_sample_columns = c(
    identical(colnames(count_mat), colnames(abundance_mat)),
    identical(colnames(count_mat), colnames(length_mat))
  ),
  counts_n_samples = ncol(count_mat),
  other_n_samples = c(ncol(abundance_mat), ncol(length_mat)),
  stringsAsFactors = FALSE
)

cat("\nMatrix dimensions:\n")
print(matrix_dimensions)
cat("\nGene ID consistency:\n")
print(gene_id_consistency)
cat("\nSample-column consistency:\n")
print(sample_column_consistency)

write.csv(matrix_dimensions, file.path(output_root, "tables", "Step24B_gene_level_matrix_dimensions.csv"), row.names = FALSE)
write.csv(gene_id_consistency, file.path(output_root, "tables", "Step24B_gene_id_order_consistency.csv"), row.names = FALSE)
write.csv(sample_column_consistency, file.path(output_root, "tables", "Step24B_sample_column_consistency.csv"), row.names = FALSE)

# -----------------------------
# Sample identity audit from matrix columns
# -----------------------------
sample_identity <- data.frame(
  sample_id = sample_ids,
  geo_accession = vapply(sample_ids, extract_gsm, character(1)),
  label_from_sample_id = vapply(sample_ids, clean_label, character(1)),
  tissue_from_sample_id = vapply(sample_ids, infer_tissue, character(1)),
  inferred_group_from_sample_id = vapply(vapply(sample_ids, clean_label, character(1)), infer_group, character(1)),
  stringsAsFactors = FALSE
)

sample_identity$file_order_index <- seq_len(nrow(sample_identity))

tissue_group_summary <- as.data.frame(
  xtabs(~ tissue_from_sample_id + inferred_group_from_sample_id, data = sample_identity),
  stringsAsFactors = FALSE
)
colnames(tissue_group_summary) <- c("tissue_inferred", "group_inferred", "n_samples")

cat("\nTissue/group summary inferred from Step24 matrix sample IDs:\n")
print(tissue_group_summary)

write.csv(sample_identity, file.path(output_root, "tables", "Step24B_matrix_sample_identity_inferred_from_columns.csv"), row.names = FALSE)
write.csv(tissue_group_summary, file.path(output_root, "tables", "Step24B_matrix_tissue_group_summary_from_columns.csv"), row.names = FALSE)

# -----------------------------
# Sample-level gene-level QC metrics
# -----------------------------
sample_matrix_qc <- data.frame(
  sample_id = sample_ids,
  geo_accession = vapply(sample_ids, extract_gsm, character(1)),
  tissue_from_sample_id = sample_identity$tissue_from_sample_id,
  inferred_group_from_sample_id = sample_identity$inferred_group_from_sample_id,
  total_estimated_counts = colSums(count_mat, na.rm = TRUE),
  detected_genes_count_gt0 = colSums(count_mat > 0, na.rm = TRUE),
  detected_genes_count_ge1 = colSums(count_mat >= 1, na.rm = TRUE),
  gene_level_TPM_sum = colSums(abundance_mat, na.rm = TRUE),
  genes_with_TPM_gt0 = colSums(abundance_mat > 0, na.rm = TRUE),
  effective_length_NA_n = colSums(is.na(length_mat)),
  effective_length_nonpositive_n = colSums(length_mat <= 0, na.rm = TRUE),
  effective_length_median = apply(length_mat, 2, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)

sample_matrix_qc_summary <- data.frame(
  metric = c(
    "n_samples",
    "mean_total_estimated_counts",
    "min_total_estimated_counts",
    "max_total_estimated_counts",
    "mean_detected_genes_count_gt0",
    "min_detected_genes_count_gt0",
    "max_detected_genes_count_gt0",
    "mean_gene_level_TPM_sum",
    "min_gene_level_TPM_sum",
    "max_gene_level_TPM_sum",
    "max_effective_length_NA_n",
    "max_effective_length_nonpositive_n"
  ),
  value = c(
    nrow(sample_matrix_qc),
    mean(sample_matrix_qc$total_estimated_counts, na.rm = TRUE),
    min(sample_matrix_qc$total_estimated_counts, na.rm = TRUE),
    max(sample_matrix_qc$total_estimated_counts, na.rm = TRUE),
    mean(sample_matrix_qc$detected_genes_count_gt0, na.rm = TRUE),
    min(sample_matrix_qc$detected_genes_count_gt0, na.rm = TRUE),
    max(sample_matrix_qc$detected_genes_count_gt0, na.rm = TRUE),
    mean(sample_matrix_qc$gene_level_TPM_sum, na.rm = TRUE),
    min(sample_matrix_qc$gene_level_TPM_sum, na.rm = TRUE),
    max(sample_matrix_qc$gene_level_TPM_sum, na.rm = TRUE),
    max(sample_matrix_qc$effective_length_NA_n, na.rm = TRUE),
    max(sample_matrix_qc$effective_length_nonpositive_n, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

cat("\nSample-level tximport matrix QC summary:\n")
print(sample_matrix_qc_summary)

write.csv(sample_matrix_qc, file.path(output_root, "tables", "Step24B_gene_level_matrix_sample_QC.csv"), row.names = FALSE)
write.csv(sample_matrix_qc_summary, file.path(output_root, "tables", "Step24B_gene_level_matrix_sample_QC_summary.csv"), row.names = FALSE)

# -----------------------------
# Optional: match matrix samples to Step23 quant-file inventory
# This does not repeat quant.sf integrity QC. It only verifies sample identity linkage.
# -----------------------------
matrix_vs_quant_inventory <- NULL
quant_inventory_summary <- data.frame(
  quant_inventory_file = quant_inventory_file,
  quant_inventory_exists = !is.na(quant_inventory_file) && file.exists(quant_inventory_file),
  n_inventory_rows = NA_integer_,
  matrix_samples_matched_by_GSM = NA_integer_,
  matrix_samples_unmatched_by_GSM = NA_integer_,
  stringsAsFactors = FALSE
)

if (!is.na(quant_inventory_file) && file.exists(quant_inventory_file)) {
  inv <- read_csv_base(quant_inventory_file)
  quant_inventory_summary$n_inventory_rows <- nrow(inv)
  gsm_col <- find_col(inv, c("^geo_accession$", "GSM", "accession"))
  path_col <- find_col(inv, c("^file$", "absolute_path", "quant_file", "path"))
  tissue_col <- find_col(inv, c("tissue"))
  group_col <- find_col(inv, c("core_group", "group", "treatment"))

  if (is.na(gsm_col)) {
    inv$geo_accession_inferred <- vapply(apply(inv, 1, paste, collapse = " "), extract_gsm, character(1))
    gsm_col <- "geo_accession_inferred"
  }

  inv_keep <- data.frame(
    inventory_geo_accession = as.character(inv[[gsm_col]]),
    inventory_quant_path = if (!is.na(path_col)) as.character(inv[[path_col]]) else NA_character_,
    inventory_tissue = if (!is.na(tissue_col)) as.character(inv[[tissue_col]]) else NA_character_,
    inventory_group = if (!is.na(group_col)) as.character(inv[[group_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  inv_keep <- inv_keep[!is.na(inv_keep$inventory_geo_accession) & inv_keep$inventory_geo_accession != "", , drop = FALSE]
  inv_keep <- inv_keep[!duplicated(inv_keep$inventory_geo_accession), , drop = FALSE]

  matrix_vs_quant_inventory <- merge(
    sample_identity,
    inv_keep,
    by.x = "geo_accession",
    by.y = "inventory_geo_accession",
    all.x = TRUE,
    sort = FALSE
  )
  matrix_vs_quant_inventory$matched_quant_inventory_by_GSM <- !is.na(matrix_vs_quant_inventory$inventory_quant_path) |
    matrix_vs_quant_inventory$geo_accession %in% inv_keep$inventory_geo_accession

  quant_inventory_summary$matrix_samples_matched_by_GSM <- sum(matrix_vs_quant_inventory$matched_quant_inventory_by_GSM, na.rm = TRUE)
  quant_inventory_summary$matrix_samples_unmatched_by_GSM <- sum(!matrix_vs_quant_inventory$matched_quant_inventory_by_GSM, na.rm = TRUE)

  write.csv(matrix_vs_quant_inventory, file.path(output_root, "tables", "Step24B_matrix_samples_vs_step23_quant_inventory.csv"), row.names = FALSE)
}

cat("\nMatrix vs Step23 quant-file inventory summary:\n")
print(quant_inventory_summary)
write.csv(quant_inventory_summary, file.path(output_root, "tables", "Step24B_matrix_vs_step23_quant_inventory_summary.csv"), row.names = FALSE)

# -----------------------------
# Optional: representative tx2gene transcript-ID matching audit
# This supports ignoreTxVersion = TRUE and does not replace Step23 full quant.sf QC.
# -----------------------------
tx2gene_match_audit <- data.frame(
  tx2gene_file = tx2gene_file,
  tx2gene_exists = file.exists(tx2gene_file),
  representative_quant_file = NA_character_,
  quant_n_transcripts = NA_integer_,
  tx2gene_n_transcripts = NA_integer_,
  matched_with_version = NA_integer_,
  matched_after_stripping_version = NA_integer_,
  match_rate_with_version_pct = NA_real_,
  match_rate_after_stripping_version_pct = NA_real_,
  status = "NOT_RUN",
  stringsAsFactors = FALSE
)

if (file.exists(tx2gene_file) && !is.na(quant_inventory_file) && file.exists(quant_inventory_file)) {
  inv <- read_csv_base(quant_inventory_file)
  path_col <- find_col(inv, c("^file$", "absolute_path", "quant_file", "path"))
  if (!is.na(path_col)) {
    candidate_paths <- as.character(inv[[path_col]])
    candidate_paths <- candidate_paths[file.exists(candidate_paths)]
    if (length(candidate_paths) > 0) {
      rep_quant <- candidate_paths[1]
      tx2gene <- read_csv_base(tx2gene_file)
      tx_col <- find_col(tx2gene, c("^transcript_id$", "transcript", "tx", "target", "Name"))

      if (!is.na(tx_col)) {
        q_tx <- tryCatch(read_quant_name_column(rep_quant), error = function(e) e)
        if (!inherits(q_tx, "error")) {
          q_tx <- unique(as.character(q_tx))
          tx_ref <- unique(as.character(tx2gene[[tx_col]]))
          tx2gene_match_audit$representative_quant_file <- rep_quant
          tx2gene_match_audit$quant_n_transcripts <- length(q_tx)
          tx2gene_match_audit$tx2gene_n_transcripts <- length(tx_ref)
          tx2gene_match_audit$matched_with_version <- sum(q_tx %in% tx_ref)
          tx2gene_match_audit$matched_after_stripping_version <- sum(strip_tx_version(q_tx) %in% strip_tx_version(tx_ref))
          tx2gene_match_audit$match_rate_with_version_pct <- round(
            100 * tx2gene_match_audit$matched_with_version / tx2gene_match_audit$quant_n_transcripts, 3
          )
          tx2gene_match_audit$match_rate_after_stripping_version_pct <- round(
            100 * tx2gene_match_audit$matched_after_stripping_version / tx2gene_match_audit$quant_n_transcripts, 3
          )
          tx2gene_match_audit$status <- "PASS"
        } else {
          tx2gene_match_audit$status <- paste("FAILED_TO_READ_REPRESENTATIVE_QUANT:", conditionMessage(q_tx))
        }
      } else {
        tx2gene_match_audit$status <- "NO_TRANSCRIPT_COLUMN_IN_TX2GENE"
      }
    } else {
      tx2gene_match_audit$status <- "NO_EXISTING_QUANT_PATH_IN_INVENTORY"
    }
  } else {
    tx2gene_match_audit$status <- "NO_QUANT_PATH_COLUMN_IN_INVENTORY"
  }
}

cat("\nRepresentative tx2gene transcript-ID matching audit:\n")
print(tx2gene_match_audit)
write.csv(tx2gene_match_audit, file.path(output_root, "tables", "Step24B_representative_tx2gene_transcript_matching_audit.csv"), row.names = FALSE)

# -----------------------------
# Optional: Step25v3 downstream subset availability audit
# This links 96-sample Step24 to 48/24 subset-locking without replacing Step25v3.
# -----------------------------
step25v3_subset_audit <- data.frame(
  file_label = c("step25v3_synovium_manifest", "step25v3_main_comparison_manifest"),
  file = c(step25v3_syn_manifest_file, step25v3_main_manifest_file),
  exists = c(file.exists(step25v3_syn_manifest_file), file.exists(step25v3_main_manifest_file)),
  n_rows = NA_integer_,
  n_control_52w = NA_integer_,
  n_aclt_alone_52w = NA_integer_,
  n_reconstruction_52w = NA_integer_,
  n_repair_52w = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(step25v3_subset_audit))) {
  f <- step25v3_subset_audit$file[i]
  if (file.exists(f)) {
    m <- read_csv_base(f)
    step25v3_subset_audit$n_rows[i] <- nrow(m)
    group_col <- find_col(m, c("^core_group$", "group", "treatment_clean", "treatment"))
    if (!is.na(group_col)) {
      step25v3_subset_audit$n_control_52w[i] <- sum(m[[group_col]] == "Control_52W", na.rm = TRUE)
      step25v3_subset_audit$n_aclt_alone_52w[i] <- sum(m[[group_col]] == "ACLT_alone_52W", na.rm = TRUE)
      step25v3_subset_audit$n_reconstruction_52w[i] <- sum(m[[group_col]] == "Reconstruction_52W", na.rm = TRUE)
      step25v3_subset_audit$n_repair_52w[i] <- sum(m[[group_col]] == "Repair_52W", na.rm = TRUE)
    }
  }
}

cat("\nStep25v3 downstream subset availability audit:\n")
print(step25v3_subset_audit)
write.csv(step25v3_subset_audit, file.path(output_root, "tables", "Step24B_step25v3_downstream_subset_availability_audit.csv"), row.names = FALSE)

# -----------------------------
# Overall pass/fail criteria for 12.3 matrix audit
# -----------------------------
expected_n_samples <- 96
expected_tissue_counts <- c(cartilage = 48, synovium = 48)
expected_per_tissue_group_n <- 12

tissue_counts <- table(sample_identity$tissue_from_sample_id, useNA = "ifany")
tissue_group_xtab <- xtabs(~ tissue_from_sample_id + inferred_group_from_sample_id, data = sample_identity)

check_df <- data.frame(
  check = c(
    "counts_matrix_has_96_samples",
    "abundance_matrix_has_96_samples",
    "length_matrix_has_96_samples",
    "counts_vs_abundance_gene_ids_identical",
    "counts_vs_length_gene_ids_identical",
    "counts_vs_abundance_sample_columns_identical",
    "counts_vs_length_sample_columns_identical",
    "no_duplicate_gene_ids_in_counts",
    "tissue_counts_are_48_cartilage_48_synovium",
    "each_tissue_group_has_12_samples",
    "max_effective_length_NA_is_zero",
    "max_effective_length_nonpositive_is_zero"
  ),
  pass = c(
    ncol(count_mat) == expected_n_samples,
    ncol(abundance_mat) == expected_n_samples,
    ncol(length_mat) == expected_n_samples,
    gene_id_consistency$identical_gene_id_order[gene_id_consistency$comparison == "counts_vs_abundance"],
    gene_id_consistency$identical_gene_id_order[gene_id_consistency$comparison == "counts_vs_effective_length"],
    sample_column_consistency$identical_sample_columns[sample_column_consistency$comparison == "counts_vs_abundance"],
    sample_column_consistency$identical_sample_columns[sample_column_consistency$comparison == "counts_vs_effective_length"],
    !any(duplicated(counts_df$gene_id)),
    all(names(expected_tissue_counts) %in% names(tissue_counts)) &&
      all(as.integer(tissue_counts[names(expected_tissue_counts)]) == expected_tissue_counts),
    all(as.integer(tissue_group_xtab[c("cartilage", "synovium"),
                                     c("ACLT_alone_52W", "Control_52W", "Reconstruction_52W", "Repair_52W")]) ==
          expected_per_tissue_group_n),
    max(sample_matrix_qc$effective_length_NA_n, na.rm = TRUE) == 0,
    max(sample_matrix_qc$effective_length_nonpositive_n, na.rm = TRUE) == 0
  ),
  stringsAsFactors = FALSE
)

overall_pass <- all(check_df$pass)

cat("\nFormal checks for Methods 12.3 matrix audit:\n")
print(check_df)
cat("\nOverall Step24B matrix-audit pass: ", overall_pass, "\n", sep = "")

write.csv(check_df, file.path(output_root, "tables", "Step24B_formal_matrix_audit_checks.csv"), row.names = FALSE)

# -----------------------------
# Final summary and methods text
# -----------------------------
final_summary <- data.frame(
  metric = c(
    "run_status",
    "method_version",
    "overall_matrix_audit_pass",
    "counts_matrix_genes",
    "counts_matrix_samples",
    "abundance_matrix_genes",
    "abundance_matrix_samples",
    "length_matrix_genes",
    "length_matrix_samples",
    "counts_vs_abundance_gene_ids_identical",
    "counts_vs_length_gene_ids_identical",
    "counts_vs_abundance_sample_columns_identical",
    "counts_vs_length_sample_columns_identical",
    "matrix_inferred_synovium_samples",
    "matrix_inferred_cartilage_samples",
    "mean_total_estimated_counts",
    "min_total_estimated_counts",
    "max_total_estimated_counts",
    "mean_detected_genes_count_gt0",
    "min_detected_genes_count_gt0",
    "max_detected_genes_count_gt0",
    "mean_gene_level_TPM_sum",
    "min_gene_level_TPM_sum",
    "max_gene_level_TPM_sum",
    "max_effective_length_NA_n",
    "max_effective_length_nonpositive_n",
    "matrix_samples_matched_to_step23_inventory_by_GSM",
    "tx2gene_match_rate_after_stripping_version_pct",
    "step25v3_synovium_manifest_rows",
    "step25v3_main_comparison_manifest_rows",
    "output_root"
  ),
  value = c(
    ifelse(overall_pass, "SUCCESS", "CHECK_WARNINGS"),
    method_version,
    overall_pass,
    nrow(count_mat),
    ncol(count_mat),
    nrow(abundance_mat),
    ncol(abundance_mat),
    nrow(length_mat),
    ncol(length_mat),
    gene_id_consistency$identical_gene_id_order[gene_id_consistency$comparison == "counts_vs_abundance"],
    gene_id_consistency$identical_gene_id_order[gene_id_consistency$comparison == "counts_vs_effective_length"],
    sample_column_consistency$identical_sample_columns[sample_column_consistency$comparison == "counts_vs_abundance"],
    sample_column_consistency$identical_sample_columns[sample_column_consistency$comparison == "counts_vs_effective_length"],
    if ("synovium" %in% names(tissue_counts)) as.integer(tissue_counts["synovium"]) else 0,
    if ("cartilage" %in% names(tissue_counts)) as.integer(tissue_counts["cartilage"]) else 0,
    round(mean(sample_matrix_qc$total_estimated_counts, na.rm = TRUE), 3),
    round(min(sample_matrix_qc$total_estimated_counts, na.rm = TRUE), 3),
    round(max(sample_matrix_qc$total_estimated_counts, na.rm = TRUE), 3),
    round(mean(sample_matrix_qc$detected_genes_count_gt0, na.rm = TRUE), 3),
    round(min(sample_matrix_qc$detected_genes_count_gt0, na.rm = TRUE), 3),
    round(max(sample_matrix_qc$detected_genes_count_gt0, na.rm = TRUE), 3),
    round(mean(sample_matrix_qc$gene_level_TPM_sum, na.rm = TRUE), 3),
    round(min(sample_matrix_qc$gene_level_TPM_sum, na.rm = TRUE), 3),
    round(max(sample_matrix_qc$gene_level_TPM_sum, na.rm = TRUE), 3),
    max(sample_matrix_qc$effective_length_NA_n, na.rm = TRUE),
    max(sample_matrix_qc$effective_length_nonpositive_n, na.rm = TRUE),
    quant_inventory_summary$matrix_samples_matched_by_GSM,
    tx2gene_match_audit$match_rate_after_stripping_version_pct,
    step25v3_subset_audit$n_rows[step25v3_subset_audit$file_label == "step25v3_synovium_manifest"],
    step25v3_subset_audit$n_rows[step25v3_subset_audit$file_label == "step25v3_main_comparison_manifest"],
    output_root
  ),
  stringsAsFactors = FALSE
)

write.csv(final_summary, file.path(output_root, "tables", "Step24B_final_summary_for_review.csv"), row.names = FALSE)

methods_text <- paste0(
  "After the 96 processed Salmon quant.sf files passed the separate quantification-level integrity QC, ",
  "tximport-derived gene-level estimated counts, TPM abundance and effective-length matrices were audited for downstream use. ",
  "The audit confirmed matrix dimensions, gene_id ordering, sample-column consistency across the three tximport-derived matrices, ",
  "sample-level total estimated counts, detected-gene counts, gene-level TPM totals and effective-length validity. ",
  "Matrix sample identifiers were additionally checked for tissue and treatment balance, and optional matching to the Step23 quant-file inventory was recorded. ",
  "This clean audit focuses on tximport gene-level matrix consistency and does not repeat the full per-file Salmon quant.sf integrity QC."
)

writeLines(methods_text, file.path(output_root, "tables", "Step24B_methods_text_for_12_3_matrix_audit.txt"))

summary_lines <- c(
  "===== Step24B clean tximport gene-level matrix audit summary =====",
  paste0("Run status: ", ifelse(overall_pass, "SUCCESS", "CHECK_WARNINGS")),
  paste0("Overall matrix audit pass: ", overall_pass),
  paste0("Counts matrix: ", nrow(count_mat), " genes x ", ncol(count_mat), " samples"),
  paste0("TPM abundance matrix: ", nrow(abundance_mat), " genes x ", ncol(abundance_mat), " samples"),
  paste0("Effective-length matrix: ", nrow(length_mat), " genes x ", ncol(length_mat), " samples"),
  paste0("Gene IDs identical counts vs abundance: ",
         gene_id_consistency$identical_gene_id_order[gene_id_consistency$comparison == "counts_vs_abundance"]),
  paste0("Gene IDs identical counts vs length: ",
         gene_id_consistency$identical_gene_id_order[gene_id_consistency$comparison == "counts_vs_effective_length"]),
  paste0("Sample columns identical counts vs abundance: ",
         sample_column_consistency$identical_sample_columns[sample_column_consistency$comparison == "counts_vs_abundance"]),
  paste0("Sample columns identical counts vs length: ",
         sample_column_consistency$identical_sample_columns[sample_column_consistency$comparison == "counts_vs_effective_length"]),
  paste0("Inferred synovium samples in 96-sample matrix: ",
         if ("synovium" %in% names(tissue_counts)) as.integer(tissue_counts["synovium"]) else 0),
  paste0("Inferred cartilage samples in 96-sample matrix: ",
         if ("cartilage" %in% names(tissue_counts)) as.integer(tissue_counts["cartilage"]) else 0),
  paste0("Max effective-length NA count per sample: ", max(sample_matrix_qc$effective_length_NA_n, na.rm = TRUE)),
  paste0("Max effective-length nonpositive count per sample: ", max(sample_matrix_qc$effective_length_nonpositive_n, na.rm = TRUE)),
  paste0("Matrix samples matched to Step23 inventory by GSM: ",
         quant_inventory_summary$matrix_samples_matched_by_GSM),
  paste0("Representative tx2gene match rate after stripping version (%): ",
         tx2gene_match_audit$match_rate_after_stripping_version_pct),
  "",
  "Key output files:",
  file.path(output_root, "tables", "Step24B_final_summary_for_review.csv"),
  file.path(output_root, "tables", "Step24B_formal_matrix_audit_checks.csv"),
  file.path(output_root, "tables", "Step24B_gene_level_matrix_dimensions.csv"),
  file.path(output_root, "tables", "Step24B_sample_column_consistency.csv"),
  file.path(output_root, "tables", "Step24B_gene_level_matrix_sample_QC_summary.csv"),
  file.path(output_root, "tables", "Step24B_methods_text_for_12_3_matrix_audit.txt")
)

writeLines(summary_lines, file.path(output_root, "logs", "Step24B_summary_to_send_me.txt"))
cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n")

save(
  input_check,
  matrix_dimensions,
  gene_id_consistency,
  sample_column_consistency,
  sample_identity,
  tissue_group_summary,
  sample_matrix_qc,
  sample_matrix_qc_summary,
  quant_inventory_summary,
  tx2gene_match_audit,
  step25v3_subset_audit,
  check_df,
  final_summary,
  file = file.path(output_root, "objects", "Step24B_clean_tximport_gene_level_matrix_audit_workspace.RData")
)

cat("\nSession information:\n")
print(sessionInfo())

run_status <- ifelse(overall_pass, "SUCCESS", "CHECK_WARNINGS")
