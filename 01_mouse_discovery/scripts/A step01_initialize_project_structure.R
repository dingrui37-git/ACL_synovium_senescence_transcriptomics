# Step01: Initialize manuscript-ready project directory structure for mouse discovery
# Purpose:
# This step creates a fully reproducible, submission-ready directory structure
# for the mouse discovery analysis starting from raw GEO data.
# It ensures standardized organization of raw data, metadata, results, figures,
# scripts, and logs, and archives this script itself for reproducibility.

## =========================
## 0. Base directory
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

## =========================
## 1. Define directory structure
## =========================

dir_list <- c(
  "00_raw_data",
  "01_metadata",
  "02_expression_matrix",
  "03_DE_analysis",
  "04_GSEA",
  "05_CellAge",
  "06_figures/Figure1",
  "06_figures/Figure2",
  "06_figures/Figure3",
  "07_tables",
  "08_logs",
  "09_scripts",
  "10_submission"
)

full_paths <- file.path(base_dir, dir_list)

## =========================
## 2. Create directories
## =========================

for (p in full_paths) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

cat("All directories created successfully.\n")

## =========================
## 3. Save directory structure summary
## =========================

dir_df <- data.frame(
  relative_path = dir_list,
  full_path = full_paths,
  stringsAsFactors = FALSE
)

write.csv(
  dir_df,
  file = file.path(base_dir, "07_tables", "step01_directory_structure.csv"),
  row.names = FALSE
)

cat("Directory structure table saved.\n")

## =========================
## 4. Save session info (for submission)
## =========================

session_file <- file.path(base_dir, "07_tables", "step01_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)

cat("Session info saved.\n")

## =========================
## 5. Archive this script
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
archive_path <- file.path(base_dir, "09_scripts", "step01_initialize_project_structure.R")

if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path), archive_path)
} else {
  writeLines("# Please manually save this script here for reproducibility.", archive_path)
}

cat("Script archived.\n")

## =========================
## 6. Console summary
## =========================

cat("\n===== STEP01 SUMMARY =====\n")
print(dir_df)

cat("\nProject initialized successfully.\n")