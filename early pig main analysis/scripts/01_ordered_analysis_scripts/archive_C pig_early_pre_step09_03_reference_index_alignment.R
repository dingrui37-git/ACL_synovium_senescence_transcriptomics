# ============================================================
# Pig early E-MTAB-6664 upstream Step 03
# Purpose:
#   Prepare the Sus scrofa Sscrofa11.1 / Ensembl release 115 reference,
#   build the Rsubread index, and align the verified original paired-end
#   gzFASTQ files for the 18 E-MTAB-6664 core samples with Rsubread::align().
#
# Manuscript-methods mapping:
#   Main analysis alignment branch based on integrity-checked original
#   gzFASTQ files. Main Rsubread::align() parameters:
#     input_format = "gzFASTQ"
#     output_format = "BAM"
#     type = "rna"
#     phredOffset = 33
# ============================================================

setwd("D:/R/ACL ∩ senescence2")

project_root <- getwd()
data_dir <- file.path(project_root, "data_raw", "E-MTAB-6664")
manifest_file <- file.path(data_dir, "E-MTAB-6664_core_validation_manifest_wide.csv")

if (!file.exists(manifest_file)) {
  stop("Cannot find core FASTQ manifest: ", manifest_file,
       "\nRun pig_early_pre_step09_01_metadata_core_manifest.R first.")
}

# -----------------------------
# 0) Install/load Rsubread
# -----------------------------
if (!requireNamespace("Rsubread", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Rsubread", ask = FALSE, update = FALSE)
}
library(Rsubread)

options(timeout = 7200)

# -----------------------------
# 1) Download/verify Sus scrofa reference FASTA and Ensembl release 115 GTF
# -----------------------------
ref_dir <- file.path(project_root, "reference", "Sus_scrofa_Ensembl115")
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

url_fasta <- "https://ftp.ensembl.org/pub/release-115/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
url_gtf   <- "https://ftp.ensembl.org/pub/release-115/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.115.gtf.gz"

fasta_file <- file.path(ref_dir, "Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz")
gtf_file   <- file.path(ref_dir, "Sus_scrofa.Sscrofa11.1.115.gtf.gz")

download_reference_if_needed <- function(url, destfile) {
  if (file.exists(destfile) && !is.na(file.info(destfile)$size) && file.info(destfile)$size > 0) {
    message("Reference file already exists: ", destfile)
    return(invisible(TRUE))
  }

  tmpfile <- paste0(destfile, ".partial")
  if (file.exists(tmpfile)) file.remove(tmpfile)

  message("Downloading reference file: ", basename(destfile))
  result <- try(
    download.file(
      url = url,
      destfile = tmpfile,
      mode = "wb",
      method = "libcurl",
      quiet = FALSE
    ),
    silent = TRUE
  )

  if (inherits(result, "try-error")) {
    if (file.exists(tmpfile)) file.remove(tmpfile)
    stop("Reference download failed for: ", url, "\n", as.character(result))
  }

  if (!file.exists(tmpfile) || file.info(tmpfile)$size <= 0) {
    if (file.exists(tmpfile)) file.remove(tmpfile)
    stop("Downloaded reference file is missing or zero size: ", destfile)
  }

  if (file.exists(destfile)) file.remove(destfile)
  file.rename(tmpfile, destfile)
  invisible(TRUE)
}

download_reference_if_needed(url_fasta, fasta_file)
download_reference_if_needed(url_gtf, gtf_file)

ref_download_summary <- data.frame(
  file = c(fasta_file, gtf_file),
  exists = file.exists(c(fasta_file, gtf_file)),
  size_bytes = file.info(c(fasta_file, gtf_file))$size,
  stringsAsFactors = FALSE
)

write.csv(
  ref_download_summary,
  file = file.path(ref_dir, "Sus_scrofa_Ensembl115_reference_file_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# 2) Build Rsubread gapped index
# -----------------------------
index_dir <- file.path(ref_dir, "Rsubread_index")
dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)

index_basename <- file.path(index_dir, "Sus_scrofa_Sscrofa11.1_gapped")
index_files_existing <- list.files(index_dir, pattern = "Sus_scrofa_Sscrofa11\\.1_gapped", full.names = TRUE)

if (length(index_files_existing) == 0) {
  message("Building Rsubread gapped index. This can take a long time.")
  buildindex(
    basename = index_basename,
    reference = fasta_file,
    gappedIndex = TRUE
  )
} else {
  message("Rsubread index files already exist. Skipping buildindex().")
}

index_files <- list.files(index_dir, pattern = "Sus_scrofa_Sscrofa11\\.1_gapped", full.names = TRUE)
if (length(index_files) == 0) {
  stop("No Rsubread index files were found after buildindex().")
}

index_summary <- data.frame(
  file = index_files,
  basename = basename(index_files),
  size_bytes = file.info(index_files)$size,
  stringsAsFactors = FALSE
)

write.csv(
  index_summary,
  file = file.path(ref_dir, "Sus_scrofa_Ensembl115_Rsubread_index_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# 3) Load manifest and check verified original gzFASTQ inputs
# -----------------------------
manifest_wide <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("sample_id", "core_group", "R1_file", "R2_file")
if (!all(required_cols %in% colnames(manifest_wide))) {
  stop("Manifest is missing required columns: ",
       paste(setdiff(required_cols, colnames(manifest_wide)), collapse = ", "))
}

manifest_wide$core_group <- factor(
  manifest_wide$core_group,
  levels = c("CON_t0", "ACLT_untreated_t7", "ACLT_untreated_t28")
)
manifest_wide <- manifest_wide[order(manifest_wide$core_group, manifest_wide$sample_id), ]
manifest_wide$core_group <- as.character(manifest_wide$core_group)

if (nrow(manifest_wide) != 18) {
  stop("Expected 18 core samples in manifest, found: ", nrow(manifest_wide))
}

fastq_input_check <- data.frame(
  sample_id = manifest_wide$sample_id,
  core_group = manifest_wide$core_group,
  R1_path = file.path(data_dir, manifest_wide$R1_file),
  R2_path = file.path(data_dir, manifest_wide$R2_file),
  stringsAsFactors = FALSE
)

fastq_input_check$R1_exists <- file.exists(fastq_input_check$R1_path)
fastq_input_check$R2_exists <- file.exists(fastq_input_check$R2_path)
fastq_input_check$R1_size_bytes <- ifelse(fastq_input_check$R1_exists, file.info(fastq_input_check$R1_path)$size, NA_real_)
fastq_input_check$R2_size_bytes <- ifelse(fastq_input_check$R2_exists, file.info(fastq_input_check$R2_path)$size, NA_real_)
fastq_input_check$pair_ok <- with(
  fastq_input_check,
  R1_exists & R2_exists & R1_size_bytes > 0 & R2_size_bytes > 0
)

write.csv(
  fastq_input_check,
  file = file.path(data_dir, "E-MTAB-6664_alignment_input_fastq_check.csv"),
  row.names = FALSE
)

if (!all(fastq_input_check$pair_ok)) {
  print(subset(fastq_input_check, !pair_ok))
  stop("Some paired FASTQ inputs are missing or zero size. Run pig_early_pre_step09_02_fastq_download_integrity.R first.")
}

# -----------------------------
# 4) Align all 18 core samples with Rsubread::align()
# -----------------------------
bam_dir <- file.path(project_root, "results", "E-MTAB-6664", "alignment_bam_batch")
log_dir <- file.path(project_root, "results", "E-MTAB-6664", "alignment_logs")
dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

status_df <- data.frame(
  sample_id = fastq_input_check$sample_id,
  core_group = fastq_input_check$core_group,
  R1_path = fastq_input_check$R1_path,
  R2_path = fastq_input_check$R2_path,
  bam_file = file.path(bam_dir, paste0(fastq_input_check$sample_id, ".BAM")),
  align_started = FALSE,
  align_success = FALSE,
  bam_exists_after = FALSE,
  bam_size_bytes = NA_real_,
  error_message = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(status_df))) {
  cat("\n==============================\n")
  cat("Processing sample:", status_df$sample_id[i], "\n")
  cat("Group:", status_df$core_group[i], "\n")

  if (file.exists(status_df$bam_file[i]) && file.info(status_df$bam_file[i])$size > 0) {
    cat("BAM already exists and is non-zero size; skipping alignment.\n")
    status_df$align_success[i] <- TRUE
    status_df$bam_exists_after[i] <- TRUE
    status_df$bam_size_bytes[i] <- file.info(status_df$bam_file[i])$size
    next
  }

  status_df$align_started[i] <- TRUE

  try_result <- try(
    align(
      index = index_basename,
      readfile1 = status_df$R1_path[i],
      readfile2 = status_df$R2_path[i],
      input_format = "gzFASTQ",
      output_file = status_df$bam_file[i],
      output_format = "BAM",
      type = "rna",
      nthreads = 1,
      phredOffset = 33
    ),
    silent = TRUE
  )

  bam_ok <- file.exists(status_df$bam_file[i]) && file.info(status_df$bam_file[i])$size > 0
  status_df$bam_exists_after[i] <- bam_ok
  status_df$bam_size_bytes[i] <- if (bam_ok) file.info(status_df$bam_file[i])$size else NA_real_
  status_df$align_success[i] <- bam_ok

  if (inherits(try_result, "try-error")) {
    status_df$error_message[i] <- as.character(try_result)
    cat("ERROR:", status_df$error_message[i], "\n")
  } else {
    cat("Finished sample:", status_df$sample_id[i], "\n")
  }

  write.csv(
    status_df,
    file = file.path(project_root, "results", "E-MTAB-6664", "E-MTAB-6664_batch_alignment_status.csv"),
    row.names = FALSE
  )
}

alignment_summary <- data.frame(
  metric = c(
    "n_samples_total",
    "n_align_success",
    "n_align_failed"
  ),
  value = c(
    nrow(status_df),
    sum(status_df$align_success, na.rm = TRUE),
    sum(!status_df$align_success, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  status_df,
  file = file.path(project_root, "results", "E-MTAB-6664", "E-MTAB-6664_batch_alignment_status.csv"),
  row.names = FALSE
)

write.csv(
  alignment_summary,
  file = file.path(project_root, "results", "E-MTAB-6664", "E-MTAB-6664_batch_alignment_summary.csv"),
  row.names = FALSE
)

cat("\n===== Final batch alignment summary =====\n")
print(alignment_summary)

if (!all(status_df$align_success)) {
  failed <- subset(status_df, !align_success)
  print(failed[, c("sample_id", "core_group", "bam_file", "error_message")])
  stop("Some samples failed alignment. See E-MTAB-6664_batch_alignment_status.csv.")
}

cat("\nAll 18 core samples were aligned or already had non-zero BAM files.\n")
