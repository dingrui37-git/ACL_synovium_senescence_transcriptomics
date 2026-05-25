# Step02: Download GEO supplementary data for GSE271903
# Purpose:
# This step downloads the official supplementary files from GEO (GSE271903),
# which contain the expression matrix used for downstream PCA/DE/GSEA analysis.
# The script ensures reproducibility by saving download logs, file records,
# and archiving itself.

## =========================
## 0. Load packages
## =========================

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  install.packages("GEOquery")
}
library(GEOquery)

## =========================
## 1. Define paths
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

raw_dir <- file.path(base_dir, "00_raw_data/GSE271903")
log_dir <- file.path(base_dir, "08_logs")
table_dir <- file.path(base_dir, "07_tables")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "step02_download_GSE271903_log.txt")

sink(log_file, split = TRUE)

cat("===== STEP02 DOWNLOAD GEO DATA =====\n")
cat("Time:", format(Sys.time()), "\n")

## =========================
## 2. Download supplementary files
## =========================

getGEOSuppFiles(
  GEO = "GSE271903",
  baseDir = raw_dir,
  makeDirectory = FALSE
)

cat("Download finished.\n")

## =========================
## 3. List downloaded files
## =========================

files <- list.files(raw_dir, recursive = TRUE, full.names = TRUE)

file_info <- data.frame(
  file_name = basename(files),
  full_path = files,
  size_MB = round(file.info(files)$size / 1024^2, 3),
  stringsAsFactors = FALSE
)

write.csv(
  file_info,
  file = file.path(table_dir, "step02_downloaded_files_GSE271903.csv"),
  row.names = FALSE
)

cat("File list saved.\n")

## =========================
## 4. Save session info
## =========================

session_file <- file.path(table_dir, "step02_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)

cat("Session info saved.\n")

## =========================
## 5. Archive script
## =========================

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }
  return(NA)
}

script_path <- get_script_path()
archive_path <- file.path(script_dir, "step02_download_GSE271903.R")

if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path), archive_path)
} else {
  writeLines("# Please manually save this script here.", archive_path)
}

cat("Script archived.\n")

## =========================
## 6. Console summary
## =========================

cat("\n===== STEP02 SUMMARY =====\n")
print(file_info)

sink()