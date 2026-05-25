# ============================================================
# Pig early E-MTAB-6664 upstream Step 02
# Purpose:
#   Download the paired-end FASTQ files for the 18 frozen core samples
#   using E-MTAB-6664_core_validation_manifest_wide.csv, verify local
#   FASTQ integrity by comparing local file sizes with remote
#   Content-Length, and re-download any missing or incomplete FASTQ files
#   in one clean final script.
#
# Manuscript-methods mapping:
#   FASTQ download and file-level integrity verification before the
#   submission rebuild.
# ============================================================

setwd("D:/R/ACL ∩ senescence2")

project_root <- getwd()
data_dir <- file.path(project_root, "data_raw", "E-MTAB-6664")
manifest_file <- file.path(data_dir, "E-MTAB-6664_core_validation_manifest_wide.csv")

if (!file.exists(manifest_file)) {
  stop("Cannot find core FASTQ manifest: ", manifest_file,
       "\nRun pig_early_pre_step09_01_metadata_core_manifest.R first.")
}

if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
library(httr)

options(timeout = 7200)

# -----------------------------
# 1) Load the frozen 18-sample FASTQ manifest
# -----------------------------
manifest_wide_6664 <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("sample_id", "core_group", "R1_file", "R1_URI", "R2_file", "R2_URI")
if (!all(required_cols %in% colnames(manifest_wide_6664))) {
  stop("Manifest is missing required columns: ",
       paste(setdiff(required_cols, colnames(manifest_wide_6664)), collapse = ", "))
}

manifest_wide_6664$core_group <- factor(
  manifest_wide_6664$core_group,
  levels = c("CON_t0", "ACLT_untreated_t7", "ACLT_untreated_t28")
)
manifest_wide_6664 <- manifest_wide_6664[order(manifest_wide_6664$core_group, manifest_wide_6664$sample_id), ]
manifest_wide_6664$core_group <- as.character(manifest_wide_6664$core_group)

if (nrow(manifest_wide_6664) != 18) {
  stop("Expected 18 core samples in manifest, found: ", nrow(manifest_wide_6664))
}

# -----------------------------
# 2) Build a file-level table for all 36 FASTQ files
# -----------------------------
fastq_table_6664 <- rbind(
  data.frame(
    sample_id = manifest_wide_6664$sample_id,
    core_group = manifest_wide_6664$core_group,
    read_label = "R1",
    file_name = manifest_wide_6664$R1_file,
    file_uri = manifest_wide_6664$R1_URI,
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = manifest_wide_6664$sample_id,
    core_group = manifest_wide_6664$core_group,
    read_label = "R2",
    file_name = manifest_wide_6664$R2_file,
    file_uri = manifest_wide_6664$R2_URI,
    stringsAsFactors = FALSE
  )
)

fastq_table_6664$file_url_https <- sub("^ftp://", "https://", fastq_table_6664$file_uri)
fastq_table_6664$local_path <- file.path(data_dir, fastq_table_6664$file_name)

fastq_table_6664$core_group <- factor(
  fastq_table_6664$core_group,
  levels = c("CON_t0", "ACLT_untreated_t7", "ACLT_untreated_t28")
)
fastq_table_6664 <- fastq_table_6664[order(fastq_table_6664$core_group, fastq_table_6664$sample_id, fastq_table_6664$read_label), ]
fastq_table_6664$core_group <- as.character(fastq_table_6664$core_group)

if (nrow(fastq_table_6664) != 36) {
  stop("Expected 36 paired FASTQ files, found: ", nrow(fastq_table_6664))
}

# -----------------------------
# 3) Helper functions for remote size, local size, and download
# -----------------------------
remote_size_head_6664 <- function(url) {
  resp <- try(httr::HEAD(url, httr::timeout(120)), silent = TRUE)
  if (inherits(resp, "try-error")) return(NA_real_)
  cl <- httr::headers(resp)[["content-length"]]
  if (is.null(cl)) return(NA_real_)
  suppressWarnings(as.numeric(cl))
}

local_size_6664 <- function(path) {
  if (file.exists(path)) {
    as.numeric(file.info(path)$size)
  } else {
    NA_real_
  }
}

download_one_fastq_6664 <- function(url, destfile) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  tmpfile <- paste0(destfile, ".partial")
  if (file.exists(tmpfile)) file.remove(tmpfile)

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
    return(as.character(result))
  }

  if (!file.exists(tmpfile) || is.na(file.info(tmpfile)$size) || file.info(tmpfile)$size <= 0) {
    if (file.exists(tmpfile)) file.remove(tmpfile)
    return("Downloaded temporary file is missing or zero size.")
  }

  if (file.exists(destfile)) file.remove(destfile)
  file.rename(tmpfile, destfile)
  return(NA_character_)
}

# -----------------------------
# 4) Query remote sizes and decide which files need download/re-download
# -----------------------------
fastq_table_6664$remote_size_bytes <- vapply(
  fastq_table_6664$file_url_https,
  remote_size_head_6664,
  numeric(1)
)

fastq_table_6664$local_exists_before <- file.exists(fastq_table_6664$local_path)
fastq_table_6664$local_size_before <- vapply(fastq_table_6664$local_path, local_size_6664, numeric(1))

fastq_table_6664$size_match_before <- with(
  fastq_table_6664,
  !is.na(local_size_before) &
    !is.na(remote_size_bytes) &
    local_size_before == remote_size_bytes
)

fastq_table_6664$needs_download <- with(
  fastq_table_6664,
  !local_exists_before |
    is.na(local_size_before) |
    local_size_before <= 0 |
    (!is.na(remote_size_bytes) & local_size_before != remote_size_bytes)
)

write.csv(
  fastq_table_6664,
  file = file.path(data_dir, "E-MTAB-6664_core_fastq_integrity_before_download.csv"),
  row.names = FALSE
)

cat("===== FASTQ integrity before download/re-download =====\n")
print(table(fastq_table_6664$needs_download, useNA = "ifany"))

# -----------------------------
# 5) Download/re-download missing or incomplete FASTQ files
# -----------------------------
fastq_table_6664$download_attempted <- FALSE
fastq_table_6664$download_error <- NA_character_

if (any(fastq_table_6664$needs_download)) {
  for (i in which(fastq_table_6664$needs_download)) {
    cat("\n==============================\n")
    cat("Downloading/re-downloading:", fastq_table_6664$file_name[i], "\n")
    cat("Sample:", fastq_table_6664$sample_id[i], " Read:", fastq_table_6664$read_label[i], "\n")

    fastq_table_6664$download_attempted[i] <- TRUE

    err <- download_one_fastq_6664(
      url = fastq_table_6664$file_url_https[i],
      destfile = fastq_table_6664$local_path[i]
    )
    fastq_table_6664$download_error[i] <- err
  }
} else {
  cat("\nAll FASTQ files already present and size-matched where remote size was available.\n")
}

# -----------------------------
# 6) Final integrity check after download/re-download
# -----------------------------
fastq_table_6664$local_exists_after <- file.exists(fastq_table_6664$local_path)
fastq_table_6664$local_size_after <- vapply(fastq_table_6664$local_path, local_size_6664, numeric(1))

fastq_table_6664$size_match_after <- with(
  fastq_table_6664,
  !is.na(local_size_after) &
    !is.na(remote_size_bytes) &
    local_size_after == remote_size_bytes
)

fastq_table_6664$final_fastq_ok <- with(
  fastq_table_6664,
  local_exists_after &
    !is.na(local_size_after) &
    local_size_after > 0 &
    (is.na(remote_size_bytes) | local_size_after == remote_size_bytes)
)

fastq_summary_6664 <- data.frame(
  metric = c(
    "n_core_samples",
    "n_fastq_expected",
    "n_download_attempted",
    "n_fastq_final_ok",
    "n_fastq_final_not_ok"
  ),
  value = c(
    length(unique(fastq_table_6664$sample_id)),
    nrow(fastq_table_6664),
    sum(fastq_table_6664$download_attempted, na.rm = TRUE),
    sum(fastq_table_6664$final_fastq_ok, na.rm = TRUE),
    sum(!fastq_table_6664$final_fastq_ok, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  fastq_table_6664,
  file = file.path(data_dir, "E-MTAB-6664_core_fastq_download_integrity_final.csv"),
  row.names = FALSE
)

write.csv(
  fastq_summary_6664,
  file = file.path(data_dir, "E-MTAB-6664_core_fastq_download_integrity_summary.csv"),
  row.names = FALSE
)

cat("\n===== Final FASTQ download/integrity summary =====\n")
print(fastq_summary_6664)

cat("\n===== Final FASTQ status table =====\n")
print(table(fastq_table_6664$final_fastq_ok, useNA = "ifany"))

if (!all(fastq_table_6664$final_fastq_ok)) {
  failed <- subset(fastq_table_6664, !final_fastq_ok)
  print(failed[, c(
    "sample_id", "core_group", "read_label", "file_name",
    "local_size_after", "remote_size_bytes", "download_error"
  )])
  stop("Some FASTQ files are still missing or size-mismatched. See E-MTAB-6664_core_fastq_download_integrity_final.csv.")
}

cat("\nAll 18 core samples have complete paired FASTQ files available locally.\n")
