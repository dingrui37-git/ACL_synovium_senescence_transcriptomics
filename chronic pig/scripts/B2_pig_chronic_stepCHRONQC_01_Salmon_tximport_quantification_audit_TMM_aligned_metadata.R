# StepCHRONQC_01_chronic_pig_Salmon_tximport_quantification_audit.R
# Purpose:
#   Audit chronic pig GSE228848 processed Salmon quant.sf files and existing tximport
#   gene-level matrices. This is quantification-level QC, not alignment-level QC and
#   not featureCounts-level QC.
#
# Revision in this version:
#   - metadata_file now points to the TMM-aligned Chronic Step2 sample metadata.

method_version <- "2026-05-27_chronic_pig_Salmon_tximport_quantification_audit_v1_TMM_aligned_metadata_path"

quant_dir <- "E:/R/ACLsenescence2/rebuild_submission/raw data/GSE228848_synovium_quant"
tables_dir <- "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables"
metadata_file <- "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_sample_metadata_used.csv"
tximport_run_summary_file <- "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_tximport_run_summary.csv"
gtf_file <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf"

output_root <- "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
for (d in c("logs","tables","source_data","figures","objects","scripts")) {
  dir.create(file.path(output_root, d), recursive = TRUE, showWarnings = FALSE)
}

# Clean previous files from this step before opening log
for (d in c("logs","tables","source_data","figures","objects","scripts")) {
  old <- list.files(file.path(output_root, d), pattern = "^StepCHRONQC_01_", full.names = TRUE)
  if (length(old) > 0) suppressWarnings(file.remove(old))
}

log_file <- file.path(output_root, "logs", "StepCHRONQC_01_chronic_pig_Salmon_tximport_QC_summary_log.txt")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
run_status <- "FAILED"
on.exit({
  cat("\n============================================================\n")
  cat("StepCHRONQC_01_chronic_pig_Salmon_tximport_QC finished with status:", run_status, "\n")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Summary log saved to:\n", log_file, "\n", sep = "")
  cat("============================================================\n")
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("StepCHRONQC_01_chronic_pig_Salmon_tximport_QC\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Method version:", method_version, "\n")
cat("============================================================\n\n")

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

for (pkg in c("dplyr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package not installed: ", pkg)
}
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
cat("Required packages loaded.\n")
cat("dplyr version:", as.character(packageVersion("dplyr")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n\n")

read_csv_base <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
read_quant <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
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
  ifelse(grepl("^CON", z), "Control_52W",
         ifelse(grepl("^ACLT", z), "ACLT_alone_52W",
                ifelse(grepl("^RECON", z), "Reconstruction_52W",
                       ifelse(grepl("^REPAIR", z), "Repair_52W", NA_character_))))
}
find_col <- function(df, patterns) {
  nms <- names(df)
  for (p in patterns) {
    hit <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}
pick_file <- function(dir, patterns) {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
  keep <- rep(FALSE, length(files))
  for (p in patterns) keep <- keep | grepl(p, basename(files), ignore.case = TRUE)
  hits <- files[keep]
  if (length(hits) == 0) return(NA_character_)
  score <- ifelse(grepl("step24", basename(hits), ignore.case = TRUE), 100, 0) +
    ifelse(grepl("pig_chronic", basename(hits), ignore.case = TRUE), 50, 0) +
    ifelse(grepl("gene_level", basename(hits), ignore.case = TRUE), 20, 0)
  hits[order(score, decreasing = TRUE)][1]
}

cat("Input paths:\n")
input_paths <- data.frame(
  label = c("quant_dir", "tables_dir", "metadata_file", "tximport_run_summary_file", "gtf_file"),
  path = c(quant_dir, tables_dir, metadata_file, tximport_run_summary_file, gtf_file),
  exists = c(dir.exists(quant_dir), dir.exists(tables_dir), file.exists(metadata_file), file.exists(tximport_run_summary_file), file.exists(gtf_file))
)
print(input_paths)
write.csv(input_paths, file.path(output_root, "tables", "StepCHRONQC_01_input_path_check.csv"), row.names = FALSE)
if (!dir.exists(quant_dir)) stop("Quant directory not found: ", quant_dir)
if (!dir.exists(tables_dir)) stop("Tables directory not found: ", tables_dir)

# Quant.sf audit
quant_files <- list.files(quant_dir, pattern = "quant\\.sf", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
quant_manifest <- data.frame(
  quant_file = quant_files,
  basename = basename(quant_files),
  geo_accession = vapply(basename(quant_files), extract_gsm, character(1)),
  label_from_filename = vapply(quant_files, clean_label, character(1)),
  tissue_from_filename = vapply(quant_files, infer_tissue, character(1)),
  stringsAsFactors = FALSE
)
quant_manifest$inferred_group_from_filename <- vapply(quant_manifest$label_from_filename, infer_group, character(1))
quant_manifest$file_size_bytes <- file.info(quant_manifest$quant_file)$size
write.csv(quant_manifest, file.path(output_root, "tables", "StepCHRONQC_01_quant_file_manifest_all.csv"), row.names = FALSE)
cat("\nQuant-like files found:", nrow(quant_manifest), "\n")
cat("Tissue counts from filename:\n"); print(table(quant_manifest$tissue_from_filename, useNA = "ifany"))
cat("Group counts from filename:\n"); print(table(quant_manifest$inferred_group_from_filename, useNA = "ifany"))

syn_quant <- quant_manifest[quant_manifest$tissue_from_filename == "synovium", , drop = FALSE]
if (nrow(syn_quant) == 0) {
  cat("WARNING: no synovium files detected by filename; using all quant files for integrity audit.\n")
  syn_quant <- quant_manifest
}
write.csv(syn_quant, file.path(output_root, "tables", "StepCHRONQC_01_quant_file_manifest_synovium_candidate.csv"), row.names = FALSE)
cat("\nSynovium-candidate quant files:", nrow(syn_quant), "\n")
print(table(syn_quant$inferred_group_from_filename, useNA = "ifany"))

standard_cols <- c("Name","Length","EffectiveLength","TPM","NumReads")
qc_list <- vector("list", nrow(syn_quant))
for (i in seq_len(nrow(syn_quant))) {
  f <- syn_quant$quant_file[i]
  cat("Auditing quant file", i, "of", nrow(syn_quant), ":", basename(f), "\n")
  q <- tryCatch(read_quant(f), error = function(e) e)
  if (inherits(q, "error")) {
    qc_list[[i]] <- data.frame(
      geo_accession = syn_quant$geo_accession[i],
      label_from_filename = syn_quant$label_from_filename[i],
      inferred_group_from_filename = syn_quant$inferred_group_from_filename[i],
      quant_file = f,
      read_status = "FAIL",
      error_message = conditionMessage(q),
      n_transcripts = NA_integer_,
      has_standard_columns = FALSE,
      missing_standard_columns = paste(standard_cols, collapse = ";"),
      tpm_sum = NA_real_,
      numreads_sum = NA_real_,
      effective_length_na_n = NA_integer_,
      effective_length_nonpositive_n = NA_integer_,
      stringsAsFactors = FALSE
    )
  } else {
    miss <- setdiff(standard_cols, names(q))
    for (cn in intersect(c("Length","EffectiveLength","TPM","NumReads"), names(q))) {
      q[[cn]] <- suppressWarnings(as.numeric(q[[cn]]))
    }
    qc_list[[i]] <- data.frame(
      geo_accession = syn_quant$geo_accession[i],
      label_from_filename = syn_quant$label_from_filename[i],
      inferred_group_from_filename = syn_quant$inferred_group_from_filename[i],
      quant_file = f,
      read_status = "PASS",
      error_message = NA_character_,
      n_transcripts = nrow(q),
      has_standard_columns = length(miss) == 0,
      missing_standard_columns = ifelse(length(miss) == 0, "none", paste(miss, collapse = ";")),
      tpm_sum = if ("TPM" %in% names(q)) sum(q$TPM, na.rm = TRUE) else NA_real_,
      numreads_sum = if ("NumReads" %in% names(q)) sum(q$NumReads, na.rm = TRUE) else NA_real_,
      effective_length_na_n = if ("EffectiveLength" %in% names(q)) sum(is.na(q$EffectiveLength)) else NA_integer_,
      effective_length_nonpositive_n = if ("EffectiveLength" %in% names(q)) sum(q$EffectiveLength <= 0, na.rm = TRUE) else NA_integer_,
      stringsAsFactors = FALSE
    )
  }
}
quant_qc <- do.call(rbind, qc_list)
write.csv(quant_qc, file.path(output_root, "tables", "StepCHRONQC_01_quant_sf_integrity_QC_by_file.csv"), row.names = FALSE)

quant_summary <- data.frame(
  metric = c("quant_files_all_found", "synovium_candidate_quant_files", "quant_files_read_PASS",
             "quant_files_standard_columns_PASS", "min_n_transcripts", "max_n_transcripts",
             "mean_n_transcripts", "min_tpm_sum", "max_tpm_sum", "mean_tpm_sum",
             "min_numreads_sum", "max_numreads_sum", "mean_numreads_sum"),
  value = c(nrow(quant_manifest), nrow(syn_quant), sum(quant_qc$read_status == "PASS"),
            sum(quant_qc$has_standard_columns),
            min(quant_qc$n_transcripts, na.rm = TRUE),
            max(quant_qc$n_transcripts, na.rm = TRUE),
            round(mean(quant_qc$n_transcripts, na.rm = TRUE), 3),
            round(min(quant_qc$tpm_sum, na.rm = TRUE), 3),
            round(max(quant_qc$tpm_sum, na.rm = TRUE), 3),
            round(mean(quant_qc$tpm_sum, na.rm = TRUE), 3),
            round(min(quant_qc$numreads_sum, na.rm = TRUE), 3),
            round(max(quant_qc$numreads_sum, na.rm = TRUE), 3),
            round(mean(quant_qc$numreads_sum, na.rm = TRUE), 3))
)
cat("\nQuant.sf integrity summary:\n"); print(quant_summary)
write.csv(quant_summary, file.path(output_root, "tables", "StepCHRONQC_01_quant_sf_integrity_summary.csv"), row.names = FALSE)

# Detect existing tximport outputs
counts_file <- pick_file(tables_dir, c("estimated.*count.*matrix", "count.*matrix"))
abundance_file <- pick_file(tables_dir, c("abundance.*tpm.*matrix", "tpm.*matrix", "abundance.*matrix"))
length_file <- pick_file(tables_dir, c("effective.*length.*matrix", "length.*matrix"))
gene_anno_file <- pick_file(tables_dir, c("gene_annotation_from_gtf", "gene.*annotation"))
tx2gene_file <- pick_file(tables_dir, c("tx2gene", "transcript.*gene"))

detected_files <- data.frame(
  label = c("counts_file", "abundance_file", "effective_length_file", "gene_annotation_file", "tx2gene_file", "tximport_run_summary_file"),
  path = c(counts_file, abundance_file, length_file, gene_anno_file, tx2gene_file, tximport_run_summary_file),
  exists = file.exists(c(counts_file, abundance_file, length_file, gene_anno_file, tx2gene_file, tximport_run_summary_file))
)
cat("\nDetected tximport/gene-level files:\n"); print(detected_files)
write.csv(detected_files, file.path(output_root, "tables", "StepCHRONQC_01_detected_tximport_matrix_files.csv"), row.names = FALSE)
if (!file.exists(counts_file)) stop("Could not detect gene-level counts matrix in: ", tables_dir)

counts_df <- read_csv_base(counts_file)
counts_sample_cols <- names(counts_df)[-1]
matrix_dims <- data.frame(
  matrix = "counts",
  path = counts_file,
  n_genes_rows = nrow(counts_df),
  n_columns_total = ncol(counts_df),
  n_sample_columns = length(counts_sample_cols),
  first_column = names(counts_df)[1],
  stringsAsFactors = FALSE
)

abundance_df <- NULL
if (file.exists(abundance_file)) {
  abundance_df <- read_csv_base(abundance_file)
  matrix_dims <- rbind(matrix_dims, data.frame(
    matrix = "abundance_TPM", path = abundance_file, n_genes_rows = nrow(abundance_df),
    n_columns_total = ncol(abundance_df), n_sample_columns = ncol(abundance_df)-1,
    first_column = names(abundance_df)[1], stringsAsFactors = FALSE))
}
length_df <- NULL
if (file.exists(length_file)) {
  length_df <- read_csv_base(length_file)
  matrix_dims <- rbind(matrix_dims, data.frame(
    matrix = "effective_length", path = length_file, n_genes_rows = nrow(length_df),
    n_columns_total = ncol(length_df), n_sample_columns = ncol(length_df)-1,
    first_column = names(length_df)[1], stringsAsFactors = FALSE))
}
cat("\nGene-level matrix dimensions:\n"); print(matrix_dims)
write.csv(matrix_dims, file.path(output_root, "tables", "StepCHRONQC_01_gene_level_matrix_dimensions.csv"), row.names = FALSE)

sample_col_consistency <- data.frame(comparison=character(), identical_sample_columns=logical(), counts_n=integer(), other_n=integer())
if (!is.null(abundance_df)) {
  sample_col_consistency <- rbind(sample_col_consistency, data.frame(
    comparison="counts_vs_abundance",
    identical_sample_columns=identical(names(counts_df)[-1], names(abundance_df)[-1]),
    counts_n=length(names(counts_df)[-1]), other_n=length(names(abundance_df)[-1])))
}
if (!is.null(length_df)) {
  sample_col_consistency <- rbind(sample_col_consistency, data.frame(
    comparison="counts_vs_effective_length",
    identical_sample_columns=identical(names(counts_df)[-1], names(length_df)[-1]),
    counts_n=length(names(counts_df)[-1]), other_n=length(names(length_df)[-1])))
}
cat("\nMatrix sample-column consistency:\n"); print(sample_col_consistency)
write.csv(sample_col_consistency, file.path(output_root, "tables", "StepCHRONQC_01_matrix_sample_column_consistency.csv"), row.names = FALSE)

# Matrix sample QC
count_mat <- as.matrix(counts_df[, counts_sample_cols, drop = FALSE])
suppressWarnings(storage.mode(count_mat) <- "numeric")
sample_qc <- data.frame(
  sample_id = counts_sample_cols,
  geo_accession = vapply(counts_sample_cols, extract_gsm, character(1)),
  label_from_sample_id = vapply(counts_sample_cols, clean_label, character(1)),
  total_estimated_counts = colSums(count_mat, na.rm = TRUE),
  detected_genes_count_gt0 = colSums(count_mat > 0, na.rm = TRUE),
  detected_genes_count_ge1 = colSums(count_mat >= 1, na.rm = TRUE),
  stringsAsFactors = FALSE
)
sample_qc$tissue_from_sample_id <- vapply(sample_qc$sample_id, infer_tissue, character(1))
sample_qc$inferred_group_from_sample_id <- vapply(sample_qc$label_from_sample_id, infer_group, character(1))

if (!is.null(abundance_df)) {
  a_cols <- names(abundance_df)[-1]
  a_mat <- as.matrix(abundance_df[, a_cols, drop=FALSE])
  suppressWarnings(storage.mode(a_mat) <- "numeric")
  abundance_qc <- data.frame(sample_id=a_cols, tpm_sum_gene_level=colSums(a_mat, na.rm=TRUE), genes_with_tpm_gt0=colSums(a_mat > 0, na.rm=TRUE))
  sample_qc <- merge(sample_qc, abundance_qc, by="sample_id", all.x=TRUE)
}
if (!is.null(length_df)) {
  l_cols <- names(length_df)[-1]
  l_mat <- as.matrix(length_df[, l_cols, drop=FALSE])
  suppressWarnings(storage.mode(l_mat) <- "numeric")
  length_qc <- data.frame(sample_id=l_cols, effective_length_na_n=colSums(is.na(l_mat)),
                          effective_length_nonpositive_n=colSums(l_mat <= 0, na.rm=TRUE),
                          effective_length_median=apply(l_mat, 2, median, na.rm=TRUE))
  sample_qc <- merge(sample_qc, length_qc, by="sample_id", all.x=TRUE)
}
cat("\nGene-level matrix sample QC preview:\n"); print(head(sample_qc, 20))
write.csv(sample_qc, file.path(output_root, "tables", "StepCHRONQC_01_gene_level_matrix_sample_QC.csv"), row.names = FALSE)

matrix_group_counts <- as.data.frame(table(sample_qc$inferred_group_from_sample_id, useNA="ifany"))
names(matrix_group_counts) <- c("inferred_group_from_sample_id", "n_samples")
cat("\nGroup counts inferred from matrix sample IDs:\n"); print(matrix_group_counts)
write.csv(matrix_group_counts, file.path(output_root, "tables", "StepCHRONQC_01_matrix_inferred_group_counts.csv"), row.names = FALSE)

# Match matrix samples to quant files by GSM
matrix_quant_audit <- merge(
  data.frame(matrix_sample_id=counts_sample_cols, matrix_geo_accession=vapply(counts_sample_cols, extract_gsm, character(1))),
  syn_quant[, c("geo_accession","quant_file","basename","label_from_filename","tissue_from_filename","inferred_group_from_filename")],
  by.x="matrix_geo_accession", by.y="geo_accession", all.x=TRUE
)
matrix_quant_audit$has_matching_synovium_quant_file <- !is.na(matrix_quant_audit$quant_file)
cat("\nMatrix samples vs quant files audit preview:\n"); print(head(matrix_quant_audit, 20))
write.csv(matrix_quant_audit, file.path(output_root, "tables", "StepCHRONQC_01_matrix_samples_vs_quant_files_audit.csv"), row.names = FALSE)

# Metadata audit
metadata_audit <- data.frame(metadata_file=metadata_file, exists=file.exists(metadata_file), n_rows=NA_integer_, n_cols=NA_integer_, inferred_sample_col=NA_character_, inferred_group_col=NA_character_)
metadata_matrix_audit <- NULL
if (file.exists(metadata_file)) {
  meta <- read_csv_base(metadata_file)
  metadata_audit$n_rows <- nrow(meta)
  metadata_audit$n_cols <- ncol(meta)
  metadata_audit$inferred_sample_col <- find_col(meta, c("^sample_id$", "geo_accession", "GSM", "accession", "sample"))
  metadata_audit$inferred_group_col <- find_col(meta, c("^group$", "treatment", "condition", "phenotype", "title"))
  cat("\nMetadata file loaded. Columns:\n"); print(names(meta))
  print(metadata_audit)
  write.csv(meta, file.path(output_root, "source_data", "StepCHRONQC_01_metadata_file_user_copy.csv"), row.names = FALSE)
  if (!is.na(metadata_audit$inferred_sample_col)) {
    mv <- as.character(meta[[metadata_audit$inferred_sample_col]])
    metadata_matrix_audit <- data.frame(metadata_sample_value=mv, metadata_geo_accession=vapply(mv, extract_gsm, character(1)))
    metadata_matrix_audit$matched_by_exact_matrix_col <- metadata_matrix_audit$metadata_sample_value %in% counts_sample_cols
    metadata_matrix_audit$matched_by_geo_accession <- metadata_matrix_audit$metadata_geo_accession %in% vapply(counts_sample_cols, extract_gsm, character(1))
    if (!is.na(metadata_audit$inferred_group_col)) metadata_matrix_audit$metadata_group_value <- as.character(meta[[metadata_audit$inferred_group_col]])
    write.csv(metadata_matrix_audit, file.path(output_root, "tables", "StepCHRONQC_01_metadata_vs_matrix_samples_audit.csv"), row.names = FALSE)
    cat("\nMetadata vs matrix sample audit summary:\n")
    cat("Metadata rows:", nrow(meta), "\n")
    cat("Matched by exact matrix column:", sum(metadata_matrix_audit$matched_by_exact_matrix_col, na.rm=TRUE), "\n")
    cat("Matched by GEO accession:", sum(metadata_matrix_audit$matched_by_geo_accession, na.rm=TRUE), "\n")
  }
}
write.csv(metadata_audit, file.path(output_root, "tables", "StepCHRONQC_01_metadata_audit.csv"), row.names = FALSE)

# tx2gene audit if available
tx2gene_audit <- data.frame(tx2gene_file=tx2gene_file, tx2gene_exists=file.exists(tx2gene_file), quant_file_used=NA_character_, quant_n_transcripts=NA_integer_, tx2gene_n_transcripts=NA_integer_, matched_with_version=NA_integer_, matched_after_stripping_version=NA_integer_, match_rate_with_version_pct=NA_real_, match_rate_after_stripping_version_pct=NA_real_)
if (file.exists(tx2gene_file) && nrow(syn_quant) > 0) {
  tx2gene <- read_csv_base(tx2gene_file)
  tx_col <- find_col(tx2gene, c("transcript", "tx", "target", "Name"))
  if (!is.na(tx_col)) {
    q <- read_quant(syn_quant$quant_file[1])
    if ("Name" %in% names(q)) {
      q_tx <- unique(as.character(q$Name)); tx_ref <- unique(as.character(tx2gene[[tx_col]]))
      tx2gene_audit$quant_file_used <- syn_quant$quant_file[1]
      tx2gene_audit$quant_n_transcripts <- length(q_tx)
      tx2gene_audit$tx2gene_n_transcripts <- length(tx_ref)
      tx2gene_audit$matched_with_version <- sum(q_tx %in% tx_ref)
      tx2gene_audit$matched_after_stripping_version <- sum(strip_tx_version(q_tx) %in% strip_tx_version(tx_ref))
      tx2gene_audit$match_rate_with_version_pct <- round(100 * tx2gene_audit$matched_with_version / tx2gene_audit$quant_n_transcripts, 3)
      tx2gene_audit$match_rate_after_stripping_version_pct <- round(100 * tx2gene_audit$matched_after_stripping_version / tx2gene_audit$quant_n_transcripts, 3)
    }
  }
}
cat("\ntx2gene transcript matching audit:\n"); print(tx2gene_audit)
write.csv(tx2gene_audit, file.path(output_root, "tables", "StepCHRONQC_01_tx2gene_transcript_matching_audit.csv"), row.names = FALSE)

# tximport run summary copy
tximport_summary_audit <- data.frame(run_summary_file=tximport_run_summary_file, exists=file.exists(tximport_run_summary_file), n_rows=NA_integer_, n_cols=NA_integer_)
if (file.exists(tximport_run_summary_file)) {
  txsum <- read_csv_base(tximport_run_summary_file)
  tximport_summary_audit$n_rows <- nrow(txsum); tximport_summary_audit$n_cols <- ncol(txsum)
  write.csv(txsum, file.path(output_root, "source_data", "StepCHRONQC_01_tximport_run_summary_copy.csv"), row.names = FALSE)
  cat("\ntximport run summary preview:\n"); print(head(txsum, 20))
}
write.csv(tximport_summary_audit, file.path(output_root, "tables", "StepCHRONQC_01_tximport_run_summary_audit.csv"), row.names = FALSE)

# Figures
fig1 <- ggplot(sample_qc, aes(x=inferred_group_from_sample_id, y=total_estimated_counts)) +
  geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(width=0.12), size=2) +
  labs(title="Chronic pig tximport gene-level estimated counts", x="Inferred group", y="Total estimated counts") +
  theme_classic(base_size=12) + theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(output_root, "figures", "StepCHRONQC_01_gene_level_estimated_count_sum_by_group.pdf"), fig1, width=7, height=5)

fig2 <- ggplot(sample_qc, aes(x=inferred_group_from_sample_id, y=detected_genes_count_gt0)) +
  geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(width=0.12), size=2) +
  labs(title="Chronic pig detected genes after tximport", x="Inferred group", y="Detected genes with estimated count > 0") +
  theme_classic(base_size=12) + theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave(file.path(output_root, "figures", "StepCHRONQC_01_detected_genes_by_group.pdf"), fig2, width=7, height=5)

if (nrow(quant_qc) > 0) {
  fig3 <- ggplot(quant_qc, aes(x=inferred_group_from_filename, y=numreads_sum)) +
    geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(width=0.12), size=2) +
    labs(title="Chronic pig Salmon quant.sf NumReads sum", x="Inferred group", y="Sum of NumReads") +
    theme_classic(base_size=12) + theme(axis.text.x=element_text(angle=30, hjust=1))
  ggsave(file.path(output_root, "figures", "StepCHRONQC_01_quant_sf_NumReads_sum_by_group.pdf"), fig3, width=7, height=5)
}

methods_text <- paste0(
  "Chronic pig Salmon/tximport quantification-level QC was performed because GSE228848 provided processed per-sample Salmon quantification files rather than raw FASTQ or BAM files. ",
  "Quantification files were audited for file completeness, standard Salmon columns, transcript counts, TPM sums, NumReads sums and effective-length fields. ",
  "Existing tximport gene-level matrices were audited for dimensions, sample-column consistency across counts, abundance and effective-length matrices, and sample-level estimated count and detected-gene summaries. ",
  "This chronic-pig QC therefore focuses on Salmon quantification integrity and tximport gene-level import consistency, rather than alignment-level or featureCounts assignment-level metrics."
)
writeLines(methods_text, file.path(output_root, "tables", "StepCHRONQC_01_methods_text_Salmon_tximport_QC.txt"))

# Final summary
ident_abund <- if ("counts_vs_abundance" %in% sample_col_consistency$comparison) sample_col_consistency$identical_sample_columns[sample_col_consistency$comparison=="counts_vs_abundance"][1] else NA
ident_len <- if ("counts_vs_effective_length" %in% sample_col_consistency$comparison) sample_col_consistency$identical_sample_columns[sample_col_consistency$comparison=="counts_vs_effective_length"][1] else NA
metadata_match_exact <- if (!is.null(metadata_matrix_audit)) sum(metadata_matrix_audit$matched_by_exact_matrix_col, na.rm=TRUE) else NA
metadata_match_gsm <- if (!is.null(metadata_matrix_audit)) sum(metadata_matrix_audit$matched_by_geo_accession, na.rm=TRUE) else NA

final_summary <- data.frame(
  metric = c("run_status","method_version","quant_files_all_found","synovium_candidate_quant_files","quant_files_read_PASS","quant_files_standard_columns_PASS","counts_matrix_genes","counts_matrix_samples","abundance_matrix_detected","effective_length_matrix_detected","counts_vs_abundance_sample_columns_identical","counts_vs_effective_length_sample_columns_identical","mean_total_estimated_counts","min_total_estimated_counts","max_total_estimated_counts","mean_detected_genes_count_gt0","min_detected_genes_count_gt0","max_detected_genes_count_gt0","mean_gene_level_TPM_sum","min_gene_level_TPM_sum","max_gene_level_TPM_sum","metadata_file_rows","metadata_vs_matrix_matched_by_exact_sample_col","metadata_vs_matrix_matched_by_geo_accession","tx2gene_match_rate_after_stripping_version_pct","tximport_run_summary_exists","output_root"),
  value = c("SUCCESS", method_version, nrow(quant_manifest), nrow(syn_quant), sum(quant_qc$read_status=="PASS"), sum(quant_qc$has_standard_columns), nrow(counts_df), length(counts_sample_cols), file.exists(abundance_file), file.exists(length_file), ident_abund, ident_len, round(mean(sample_qc$total_estimated_counts, na.rm=TRUE),3), round(min(sample_qc$total_estimated_counts, na.rm=TRUE),3), round(max(sample_qc$total_estimated_counts, na.rm=TRUE),3), round(mean(sample_qc$detected_genes_count_gt0, na.rm=TRUE),3), round(min(sample_qc$detected_genes_count_gt0, na.rm=TRUE),3), round(max(sample_qc$detected_genes_count_gt0, na.rm=TRUE),3), if ("tpm_sum_gene_level" %in% names(sample_qc)) round(mean(sample_qc$tpm_sum_gene_level, na.rm=TRUE),3) else NA, if ("tpm_sum_gene_level" %in% names(sample_qc)) round(min(sample_qc$tpm_sum_gene_level, na.rm=TRUE),3) else NA, if ("tpm_sum_gene_level" %in% names(sample_qc)) round(max(sample_qc$tpm_sum_gene_level, na.rm=TRUE),3) else NA, metadata_audit$n_rows, metadata_match_exact, metadata_match_gsm, tx2gene_audit$match_rate_after_stripping_version_pct, file.exists(tximport_run_summary_file), output_root),
  stringsAsFactors = FALSE
)
cat("\nFinal summary for review:\n"); print(final_summary)
write.csv(final_summary, file.path(output_root, "tables", "StepCHRONQC_01_final_summary_for_review.csv"), row.names = FALSE)

cat("\nKey output files:\n")
key_files <- c(
  "tables/StepCHRONQC_01_quant_file_manifest_all.csv",
  "tables/StepCHRONQC_01_quant_file_manifest_synovium_candidate.csv",
  "tables/StepCHRONQC_01_quant_sf_integrity_QC_by_file.csv",
  "tables/StepCHRONQC_01_quant_sf_integrity_summary.csv",
  "tables/StepCHRONQC_01_gene_level_matrix_dimensions.csv",
  "tables/StepCHRONQC_01_matrix_sample_column_consistency.csv",
  "tables/StepCHRONQC_01_gene_level_matrix_sample_QC.csv",
  "tables/StepCHRONQC_01_matrix_samples_vs_quant_files_audit.csv",
  "tables/StepCHRONQC_01_metadata_vs_matrix_samples_audit.csv",
  "tables/StepCHRONQC_01_tx2gene_transcript_matching_audit.csv",
  "tables/StepCHRONQC_01_methods_text_Salmon_tximport_QC.txt",
  "tables/StepCHRONQC_01_final_summary_for_review.csv"
)
for (i in seq_along(key_files)) cat(i, ") ", file.path(output_root, key_files[i]), "\n", sep="")

cat("\nSession information:\n")
print(sessionInfo())

run_status <- "SUCCESS"
