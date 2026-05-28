# Step 23: GSE228848 processed Salmon quant-level QC for all 96 quant.sf files
# Purpose:
#   This script performs file-level and quantification-level QC for the GSE228848
#   processed Salmon quant.sf files before tximport. It is designed to support the
#   chronic pig 96-sample tximport route: all 96 processed quantification files are
#   checked first, then downstream steps can subset the tximport-derived matrices
#   to 48 synovium samples and finally to the 24-sample Control_52W vs ACLT_alone_52W
#   main comparison.
#
# Main checks:
#   1) Find all quant.sf-like files in the confirmed GSE228848 quantification folder.
#   2) Infer tissue and treatment group from filenames.
#   3) Read every quant.sf file robustly, including .gz and .txt.gz files.
#   4) Check Salmon standard columns: Name, Length, EffectiveLength, TPM, NumReads.
#   5) Record transcript count, TPM sum, NumReads sum, and read errors if any.
#   6) Save manuscript-ready QC tables and a compact summary log for review.

# -----------------------------
# 0) Project paths
# -----------------------------
project_root <- "E:/R/ACLsenescence2"
if (!dir.exists(project_root)) {
  stop("Project root does not exist: ", project_root)
}
setwd(project_root)

rebuild_root <- file.path(project_root, "rebuild_submission")
chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

scripts_dir <- file.path(chronic_dir, "scripts")
tables_dir  <- file.path(chronic_dir, "tables")
logs_dir    <- file.path(chronic_dir, "logs")
objects_dir <- file.path(chronic_dir, "objects")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)

# This is the folder you checked manually. It contains cartilage and synovium files.
quant_dir <- file.path(rebuild_root, "raw data", "GSE228848_synovium_quant")
if (!dir.exists(quant_dir)) {
  stop("Quantification directory does not exist: ", quant_dir)
}

# Archive this script for reproducibility.
script_out <- file.path(scripts_dir, "step23_pig_chronic_quant_level_QC_96_files_FIXED.R")
try({
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
  if (!is.na(this_file) && file.exists(this_file)) {
    file.copy(this_file, script_out, overwrite = TRUE)
  } else {
    writeLines(readLines(commandArgs(trailingOnly = FALSE)[1], warn = FALSE), script_out, useBytes = TRUE)
  }
}, silent = TRUE)

# -----------------------------
# 1) Locate quant.sf-like files
# -----------------------------
all_files <- list.files(
  path = quant_dir,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

quant_pattern <- "quant\\.sf(\\.txt)?(\\.gz)?$"
quant_files <- all_files[grepl(quant_pattern, basename(all_files), ignore.case = TRUE)]
quant_files <- sort(normalizePath(quant_files, winslash = "/", mustWork = TRUE))

if (length(quant_files) == 0) {
  stop("No quant.sf-like files found under: ", quant_dir)
}

# -----------------------------
# 2) Filename-based annotation helpers
# -----------------------------
infer_tissue <- function(x) {
  b <- basename(x)
  ifelse(grepl("synovium", b, ignore.case = TRUE), "synovium",
         ifelse(grepl("cartilage", b, ignore.case = TRUE), "cartilage", "unknown"))
}

infer_group <- function(x) {
  b <- basename(x)
  ifelse(grepl("_CON[0-9]+_", b, ignore.case = TRUE), "Control_52W",
         ifelse(grepl("_ACLT[0-9]+_", b, ignore.case = TRUE), "ACLT_alone_52W",
                ifelse(grepl("_RECON[0-9]+_", b, ignore.case = TRUE), "Reconstruction_52W",
                       ifelse(grepl("_REPAIR[0-9]+_", b, ignore.case = TRUE), "Repair_52W", "unknown"))))
}

extract_geo <- function(x) {
  b <- basename(x)
  out <- sub("^(GSM[0-9]+).*$", "\\1", b)
  ifelse(grepl("^GSM[0-9]+$", out), out, NA_character_)
}

extract_sample_id <- function(x) {
  b <- basename(x)
  b <- sub("\\.gz$", "", b, ignore.case = TRUE)
  b <- sub("\\.txt$", "", b, ignore.case = TRUE)
  b <- sub("_quant\\.sf$", "", b, ignore.case = TRUE)
  b <- sub("\\.sf$", "", b, ignore.case = TRUE)
  b
}

inventory_df <- data.frame(
  sample_id = extract_sample_id(quant_files),
  geo_accession = extract_geo(quant_files),
  tissue_inferred = infer_tissue(quant_files),
  group_inferred = infer_group(quant_files),
  filename = basename(quant_files),
  file = quant_files,
  file_exists = file.exists(quant_files),
  size_bytes = as.numeric(file.info(quant_files)$size),
  stringsAsFactors = FALSE
)
inventory_df$size_mb <- round(inventory_df$size_bytes / 1024^2, 3)

# -----------------------------
# 3) Robust reader for quant.sf files
# -----------------------------
required_cols <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")
expected_n_transcripts <- 46295L

to_numeric_safe <- function(z) suppressWarnings(as.numeric(z))

read_quant_with_method <- function(f, method) {
  if (method == "direct") {
    return(utils::read.delim(
      file = f,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (method == "gzfile") {
    con <- gzfile(f, open = "rt")
    on.exit(close(con), add = TRUE)
    return(utils::read.delim(
      file = con,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  stop("Unknown read method: ", method)
}

qc_one_file <- function(f) {
  methods <- if (grepl("\\.gz$", f, ignore.case = TRUE)) c("gzfile", "direct") else c("direct", "gzfile")
  q <- NULL
  used_method <- NA_character_
  err_msgs <- character(0)

  for (m in methods) {
    attempt <- try(read_quant_with_method(f, m), silent = TRUE)
    if (!inherits(attempt, "try-error")) {
      q <- attempt
      used_method <- m
      break
    } else {
      err_msgs <- c(err_msgs, paste0(m, ": ", as.character(attempt)))
    }
  }

  if (is.null(q)) {
    return(data.frame(
      readable = FALSE,
      read_method = NA_character_,
      read_error = paste(err_msgs, collapse = " || "),
      n_transcripts = NA_integer_,
      has_Name = FALSE,
      has_Length = FALSE,
      has_EffectiveLength = FALSE,
      has_TPM = FALSE,
      has_NumReads = FALSE,
      has_all_required_cols = FALSE,
      tpm_sum = NA_real_,
      numreads_sum = NA_real_,
      length_min = NA_real_,
      effective_length_min = NA_real_,
      tpm_sum_close_1e6 = FALSE,
      expected_transcript_count_46295 = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  cn <- colnames(q)
  has_vec <- required_cols %in% cn
  names(has_vec) <- required_cols
  has_all <- all(has_vec)

  tpm_sum <- NA_real_
  numreads_sum <- NA_real_
  length_min <- NA_real_
  effective_length_min <- NA_real_

  if ("TPM" %in% cn) {
    tpm_sum <- sum(to_numeric_safe(q$TPM), na.rm = TRUE)
  }
  if ("NumReads" %in% cn) {
    numreads_sum <- sum(to_numeric_safe(q$NumReads), na.rm = TRUE)
  }
  if ("Length" %in% cn) {
    length_min <- suppressWarnings(min(to_numeric_safe(q$Length), na.rm = TRUE))
  }
  if ("EffectiveLength" %in% cn) {
    effective_length_min <- suppressWarnings(min(to_numeric_safe(q$EffectiveLength), na.rm = TRUE))
  }

  data.frame(
    readable = TRUE,
    read_method = used_method,
    read_error = "",
    n_transcripts = nrow(q),
    has_Name = unname(has_vec["Name"]),
    has_Length = unname(has_vec["Length"]),
    has_EffectiveLength = unname(has_vec["EffectiveLength"]),
    has_TPM = unname(has_vec["TPM"]),
    has_NumReads = unname(has_vec["NumReads"]),
    has_all_required_cols = has_all,
    tpm_sum = tpm_sum,
    numreads_sum = numreads_sum,
    length_min = length_min,
    effective_length_min = effective_length_min,
    tpm_sum_close_1e6 = is.finite(tpm_sum) && abs(tpm_sum - 1000000) <= 1,
    expected_transcript_count_46295 = !is.na(nrow(q)) && nrow(q) == expected_n_transcripts,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 4) Run QC across all files
# -----------------------------
qc_list <- vector("list", length(quant_files))
for (i in seq_along(quant_files)) {
  message(sprintf("QC %03d/%03d: %s", i, length(quant_files), basename(quant_files[i])))
  qc_list[[i]] <- qc_one_file(quant_files[i])
}

qc_metrics_df <- do.call(rbind, qc_list)
qc_table <- cbind(inventory_df, qc_metrics_df)

# -----------------------------
# 5) Summaries
# -----------------------------
finite_tpm <- qc_table$tpm_sum[is.finite(qc_table$tpm_sum)]
finite_numreads <- qc_table$numreads_sum[is.finite(qc_table$numreads_sum)]
unique_n_tx <- sort(unique(qc_table$n_transcripts[!is.na(qc_table$n_transcripts)]))

overall_qc_pass <- all(qc_table$readable) &&
  all(qc_table$has_all_required_cols) &&
  all(qc_table$tpm_sum_close_1e6) &&
  all(qc_table$expected_transcript_count_46295)

summary_df <- data.frame(
  metric = c(
    "quant_dir",
    "n_quant_sf_like_files_found",
    "n_readable_files",
    "n_files_with_all_salmon_standard_columns",
    "n_files_with_expected_transcript_count_46295",
    "n_files_with_TPM_sum_close_to_1000000",
    "unique_transcript_count_values",
    "TPM_sum_min",
    "TPM_sum_max",
    "NumReads_sum_min",
    "NumReads_sum_max",
    "n_synovium_files_inferred",
    "n_cartilage_files_inferred",
    "n_unknown_tissue_files_inferred",
    "n_Control_52W_files_inferred",
    "n_ACLT_alone_52W_files_inferred",
    "n_Reconstruction_52W_files_inferred",
    "n_Repair_52W_files_inferred",
    "overall_qc_pass"
  ),
  value = c(
    quant_dir,
    length(quant_files),
    sum(qc_table$readable, na.rm = TRUE),
    sum(qc_table$has_all_required_cols, na.rm = TRUE),
    sum(qc_table$expected_transcript_count_46295, na.rm = TRUE),
    sum(qc_table$tpm_sum_close_1e6, na.rm = TRUE),
    paste(unique_n_tx, collapse = ";"),
    ifelse(length(finite_tpm) > 0, format(min(finite_tpm), digits = 12), NA_character_),
    ifelse(length(finite_tpm) > 0, format(max(finite_tpm), digits = 12), NA_character_),
    ifelse(length(finite_numreads) > 0, format(min(finite_numreads), digits = 12), NA_character_),
    ifelse(length(finite_numreads) > 0, format(max(finite_numreads), digits = 12), NA_character_),
    sum(qc_table$tissue_inferred == "synovium", na.rm = TRUE),
    sum(qc_table$tissue_inferred == "cartilage", na.rm = TRUE),
    sum(qc_table$tissue_inferred == "unknown", na.rm = TRUE),
    sum(qc_table$group_inferred == "Control_52W", na.rm = TRUE),
    sum(qc_table$group_inferred == "ACLT_alone_52W", na.rm = TRUE),
    sum(qc_table$group_inferred == "Reconstruction_52W", na.rm = TRUE),
    sum(qc_table$group_inferred == "Repair_52W", na.rm = TRUE),
    as.character(overall_qc_pass)
  ),
  stringsAsFactors = FALSE
)

tissue_group_summary <- as.data.frame.matrix(table(qc_table$tissue_inferred, qc_table$group_inferred))
tissue_group_summary$tissue_inferred <- rownames(tissue_group_summary)
tissue_group_summary <- tissue_group_summary[, c("tissue_inferred", setdiff(colnames(tissue_group_summary), "tissue_inferred")), drop = FALSE]
rownames(tissue_group_summary) <- NULL

read_error_df <- qc_table[!qc_table$readable | !qc_table$has_all_required_cols,
                          c("sample_id", "filename", "file", "readable", "read_method", "read_error",
                            "has_Name", "has_Length", "has_EffectiveLength", "has_TPM", "has_NumReads"),
                          drop = FALSE]

# -----------------------------
# 6) Save outputs
# -----------------------------
write.csv(inventory_df,
          file.path(tables_dir, "step23_pig_chronic_quant_file_inventory_FIXED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(qc_table,
          file.path(tables_dir, "step23_pig_chronic_quant_level_QC_table_FIXED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(summary_df,
          file.path(tables_dir, "step23_pig_chronic_quant_level_QC_summary_FIXED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(tissue_group_summary,
          file.path(tables_dir, "step23_pig_chronic_quant_level_QC_tissue_group_summary_FIXED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(read_error_df,
          file.path(tables_dir, "step23_pig_chronic_quant_level_QC_read_errors_FIXED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

save(qc_table, summary_df, tissue_group_summary, read_error_df,
     file = file.path(objects_dir, "step23_pig_chronic_quant_level_QC_96_files_FIXED_workspace.RData"))

summary_lines <- c(
  "===== Step 23 FIXED GSE228848 processed Salmon quant-level QC summary =====",
  paste("Quantification directory:", quant_dir),
  "",
  paste("Quant.sf-like files found:", length(quant_files)),
  paste("Readable files:", sum(qc_table$readable, na.rm = TRUE), "/", length(quant_files)),
  paste("Files with all Salmon standard columns [Name, Length, EffectiveLength, TPM, NumReads]:",
        sum(qc_table$has_all_required_cols, na.rm = TRUE), "/", length(quant_files)),
  paste("Files with expected transcript count 46295:",
        sum(qc_table$expected_transcript_count_46295, na.rm = TRUE), "/", length(quant_files)),
  paste("Files with TPM sum close to 1,000,000:",
        sum(qc_table$tpm_sum_close_1e6, na.rm = TRUE), "/", length(quant_files)),
  paste("Unique transcript count values observed:", paste(unique_n_tx, collapse = "; ")),
  paste("TPM sum range:",
        ifelse(length(finite_tpm) > 0, format(min(finite_tpm), digits = 12), "NA"), "to",
        ifelse(length(finite_tpm) > 0, format(max(finite_tpm), digits = 12), "NA")),
  paste("NumReads sum range:",
        ifelse(length(finite_numreads) > 0, format(min(finite_numreads), digits = 12), "NA"), "to",
        ifelse(length(finite_numreads) > 0, format(max(finite_numreads), digits = 12), "NA")),
  "",
  paste("Inferred synovium files:", sum(qc_table$tissue_inferred == "synovium", na.rm = TRUE)),
  paste("Inferred cartilage files:", sum(qc_table$tissue_inferred == "cartilage", na.rm = TRUE)),
  paste("Unknown tissue files:", sum(qc_table$tissue_inferred == "unknown", na.rm = TRUE)),
  "",
  "Tissue/group breakdown:",
  paste(capture.output(print(tissue_group_summary, row.names = FALSE)), collapse = "\n"),
  "",
  paste("Overall QC pass:", overall_qc_pass),
  "",
  "Key output files:",
  file.path(tables_dir, "step23_pig_chronic_quant_file_inventory_FIXED.csv"),
  file.path(tables_dir, "step23_pig_chronic_quant_level_QC_table_FIXED.csv"),
  file.path(tables_dir, "step23_pig_chronic_quant_level_QC_summary_FIXED.csv"),
  file.path(tables_dir, "step23_pig_chronic_quant_level_QC_tissue_group_summary_FIXED.csv"),
  file.path(tables_dir, "step23_pig_chronic_quant_level_QC_read_errors_FIXED.csv")
)

writeLines(summary_lines,
           file.path(tables_dir, "step23_pig_chronic_quant_level_QC_summary_to_send_me_FIXED.txt"),
           useBytes = TRUE)
writeLines(summary_lines,
           file.path(logs_dir, "step23_pig_chronic_quant_level_QC_FIXED_log.txt"),
           useBytes = TRUE)

cat(paste(summary_lines, collapse = "\n"), "\n")
message("Step 23 FIXED quant-level QC finished.")
