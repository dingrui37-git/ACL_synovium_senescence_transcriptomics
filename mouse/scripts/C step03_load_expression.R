# Step03: Load and construct expression matrix and metadata from GEO supplementary files
# Purpose:
# This step reads the GEO count matrix and sample metadata, constructs a clean
# expression matrix and annotation table, and prepares standardized inputs for
# PCA and downstream analysis. It ensures reproducibility by saving outputs,
# logs, and archiving this script.

## =========================
## 0. Load packages
## =========================

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

## =========================
## 1. Define paths
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

raw_dir <- file.path(base_dir, "00_raw_data/GSE271903")
expr_dir <- file.path(base_dir, "02_expression_matrix")
meta_dir <- file.path(base_dir, "01_metadata")
table_dir <- file.path(base_dir, "07_tables")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "step03_load_expression_log.txt")
sink(log_file, split = TRUE)

cat("===== STEP03 LOAD EXPRESSION =====\n")

## =========================
## 2. Read count matrix
## =========================

count_file <- file.path(raw_dir, "GSE271903_count_matrix.xlsx")

count_df <- read_excel(count_file)

cat("Count matrix loaded.\n")
cat("Dimensions:\n")
print(dim(count_df))

## =========================
## 3. Construct expression matrix
## =========================

# 第一列是基因名
gene_col <- count_df[[1]]
expr_mat <- as.data.frame(count_df[, -1])
rownames(expr_mat) <- gene_col

# 转为 numeric
expr_mat <- as.data.frame(lapply(expr_mat, as.numeric))
rownames(expr_mat) <- gene_col

cat("Expression matrix constructed.\n")

## =========================
## 4. Read metadata
## =========================

meta_file <- file.path(raw_dir, "GSE271903_Mouse_Identifiers_and_Groups_2_NaN_.xlsx")

meta_df <- read_excel(meta_file)

cat("Metadata loaded.\n")
print(head(meta_df))

## =========================
## 5. Save cleaned data
## =========================

write.csv(
  expr_mat,
  file = file.path(expr_dir, "step03_expression_matrix_raw_counts.csv")
)

write.csv(
  meta_df,
  file = file.path(meta_dir, "step03_metadata_raw.csv"),
  row.names = FALSE
)

cat("Data saved.\n")

## =========================
## 6. Save session info
## =========================

session_file <- file.path(table_dir, "step03_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)

## =========================
## 7. Archive script
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
archive_path <- file.path(script_dir, "step03_load_expression.R")

if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path), archive_path)
} else {
  writeLines("# Please manually save this script here.", archive_path)
}

cat("Script archived.\n")

## =========================
## 8. Summary
## =========================

cat("\n===== STEP03 SUMMARY =====\n")
cat("Expression matrix dim:\n")
print(dim(expr_mat))

cat("\nMetadata dim:\n")
print(dim(meta_df))

sink()