# Step11A: Reconstruct CellAge clean table from raw CellAge file in project raw_data/cellAge folder.
# English note:
# This script reconstructs the cleaned CellAge human gene table used for Figure 3C.
# It only handles CellAge raw-table loading and cleaning. It does not perform gorth mapping
# and does not compute the mouse persistent ∩ CellAge overlap.
#
# Locked purpose:
# Raw CellAge table -> cleaned unique CellAge human symbols (cellage_clean_271903).
#
# Expected output:
# - step11A_CellAge_raw_detected.csv
# - step11A_CellAge_clean_271903.csv
# - step11A_CellAge_clean_symbol_list.csv
# - step11A_CellAge_cleaning_summary.csv
#
# The next step can intersect this cleaned CellAge table with the archived strict 1:1
# mouse-to-human ortholog mapping table.

options(stringsAsFactors = FALSE)

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"
raw_dir <- file.path(base_dir, "00_raw_data", "cellAge")

table_dir <- file.path(base_dir, "07_tables", "step11A_CellAge_raw_to_clean")
log_dir <- file.path(base_dir, "08_logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step11A_CellAge_raw_to_clean"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  log_file
)
unlink(old_outputs[file.exists(old_outputs)], recursive = TRUE, force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP11A CELLAGE RAW TO CLEAN =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("CellAge raw data directory: ", raw_dir, "\n\n", sep = "")


## Archive this script if running via Rscript.
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }
  return(NA_character_)
}
script_dir <- file.path(base_dir, "09_scripts")
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
script_path <- get_script_path()
archive_path <- file.path(script_dir, "step11A_CellAge_raw_to_clean_from_cellAge_folder.R")
if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path, warn = FALSE), archive_path, useBytes = TRUE)
  cat("Archived script: ", archive_path, "\n", sep = "")
} else {
  writeLines("# Script path not detected. Please manually save the executed script here.", archive_path)
  cat("Script path not detected; archive fallback saved: ", archive_path, "\n", sep = "")
}

safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (cc in colnames(x)) {
    if (is.list(x[[cc]])) {
      x[[cc]] <- vapply(x[[cc]], function(v) paste(as.character(v), collapse = ";"), character(1))
    }
  }
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "-", ".")] <- NA_character_
  x
}

normalize_colnames <- function(x) {
  y <- trimws(x)
  y <- gsub("\u00A0", " ", y, fixed = TRUE)
  y <- trimws(y)
  y
}

read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("csv")) {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (ext %in% c("tsv", "txt")) {
    ## Try tab-delimited first; if it produces one column, try comma.
    tab <- tryCatch(
      read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (!is.null(tab) && ncol(tab) > 1) return(tab)

    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (ext %in% c("xlsx", "xls")) {
    safe_library("readxl")
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }

  stop("Unsupported CellAge file extension: ", ext)
}

if (!dir.exists(raw_dir)) {
  stop("CellAge raw data directory does not exist: ", raw_dir, "\nExpected CellAge file such as cellage3.tsv under: ", raw_dir)
}

## Candidate CellAge raw files.
candidate_files <- list.files(
  raw_dir,
  pattern = "(CellAge|cellage|cells|senescence).*[.](csv|tsv|txt|xlsx|xls)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(candidate_files) == 0) {
  all_tables <- list.files(
    raw_dir,
    pattern = "[.](csv|tsv|txt|xlsx|xls)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  write_csv(
    data.frame(path = all_tables, stringsAsFactors = FALSE),
    file.path(table_dir, "step11A_all_table_files_seen_when_CellAge_not_found.csv")
  )
  stop("No CellAge-like raw file found in raw_dir. See all_table_files_seen file.")
}

file_info <- file.info(candidate_files)
candidate_summary <- data.frame(
  path = candidate_files,
  file = basename(candidate_files),
  size = file_info$size,
  ctime = as.character(file_info$ctime),
  mtime = as.character(file_info$mtime),
  md5 = unname(tools::md5sum(candidate_files)),
  stringsAsFactors = FALSE
)
write_csv(candidate_summary, file.path(table_dir, "step11A_CellAge_candidate_files.csv"))

cat("Candidate CellAge files:\n")
print(candidate_summary)

## Pick the most likely CellAge gene table:
## Prefer file names containing CellAge; then largest size.
candidate_summary$score <- 0
candidate_summary$score <- candidate_summary$score + ifelse(grepl("CellAge|cellage", candidate_summary$file, ignore.case = TRUE), 10, 0)
candidate_summary$score <- candidate_summary$score + ifelse(grepl("gene|genes|database|data", candidate_summary$file, ignore.case = TRUE), 3, 0)
candidate_summary <- candidate_summary[order(-candidate_summary$score, -candidate_summary$size, candidate_summary$file), ]

cellage_file <- candidate_summary$path[1]
cat("\nSelected CellAge raw file:\n", cellage_file, "\n", sep = "")
cat("Inferred CellAge version from filename: ", tools::file_path_sans_ext(basename(cellage_file)), "\n", sep = "")
cat("Selected CellAge file MD5: ", unname(tools::md5sum(cellage_file)), "\n\n", sep = "")

cellage_raw <- read_any_table(cellage_file)
colnames(cellage_raw) <- normalize_colnames(colnames(cellage_raw))

## Some CellAge CSVs may include an empty first column from row names.
empty_name_cols <- which(colnames(cellage_raw) == "" | grepl("^\\.\\.\\.", colnames(cellage_raw)))
if (length(empty_name_cols) > 0) {
  ## Drop empty columns only if they have no informative values.
  drop_cols <- empty_name_cols[vapply(cellage_raw[empty_name_cols], function(v) all(is.na(v) | v == ""), logical(1))]
  if (length(drop_cols) > 0) {
    cellage_raw <- cellage_raw[, -drop_cols, drop = FALSE]
  }
}

required_cols <- c(
  "Entrez ID",
  "Gene symbol",
  "Gene name",
  "Cancer Cell",
  "Type of senescence",
  "Senescence Effect",
  "Reference"
)

## Allow slightly different spelling for Gene symbol if needed.
if (!("Gene symbol" %in% colnames(cellage_raw))) {
  symbol_candidates <- colnames(cellage_raw)[
    grepl("^gene.*symbol$|symbol", colnames(cellage_raw), ignore.case = TRUE)
  ]
  if (length(symbol_candidates) > 0) {
    colnames(cellage_raw)[match(symbol_candidates[1], colnames(cellage_raw))] <- "Gene symbol"
  }
}

missing_required <- setdiff(required_cols, colnames(cellage_raw))
if (length(missing_required) > 0) {
  write_csv(
    data.frame(columns_detected = colnames(cellage_raw), stringsAsFactors = FALSE),
    file.path(table_dir, "step11A_CellAge_columns_detected_when_required_missing.csv")
  )
  stop(
    "CellAge raw file is missing required columns: ",
    paste(missing_required, collapse = ", "),
    "\nDetected columns were saved."
  )
}

## Keep standard CellAge columns in standard order.
cellage_raw_standard <- cellage_raw[, required_cols, drop = FALSE]

cellage_raw_standard[["Gene symbol"]] <- clean_symbol(cellage_raw_standard[["Gene symbol"]])
cellage_raw_standard[["Senescence Effect"]] <- trimws(as.character(cellage_raw_standard[["Senescence Effect"]]))
cellage_raw_standard[["Senescence Effect"]][is.na(cellage_raw_standard[["Senescence Effect"]]) |
                                             cellage_raw_standard[["Senescence Effect"]] == ""] <- "Unclear"

write_csv(cellage_raw_standard, file.path(table_dir, "step11A_CellAge_raw_detected.csv"))

## Clean table:
## - remove rows with missing Gene symbol
## - create CellAge_symbol
## - keep one row per CellAge_symbol
## If duplicate gene symbols exist, keep the first occurrence in original file order.
cellage_clean <- cellage_raw_standard
cellage_clean$CellAge_symbol <- clean_symbol(cellage_clean[["Gene symbol"]])
cellage_clean <- cellage_clean[!is.na(cellage_clean$CellAge_symbol) & cellage_clean$CellAge_symbol != "", , drop = FALSE]

duplicate_symbols <- cellage_clean$CellAge_symbol[duplicated(cellage_clean$CellAge_symbol)]
duplicate_rows <- cellage_clean[cellage_clean$CellAge_symbol %in% duplicate_symbols, , drop = FALSE]

if (nrow(duplicate_rows) > 0) {
  write_csv(duplicate_rows, file.path(table_dir, "step11A_CellAge_duplicate_symbol_rows_removed_by_first_occurrence_rule.csv"))
}

cellage_clean <- cellage_clean[!duplicated(cellage_clean$CellAge_symbol), , drop = FALSE]

## Stable ordering by CellAge symbol.
cellage_clean <- cellage_clean[order(cellage_clean$CellAge_symbol), , drop = FALSE]

write_csv(cellage_clean, file.path(table_dir, "step11A_CellAge_clean_271903.csv"))

symbol_list <- data.frame(CellAge_symbol = cellage_clean$CellAge_symbol, stringsAsFactors = FALSE)
write_csv(symbol_list, file.path(table_dir, "step11A_CellAge_clean_symbol_list.csv"))

effect_counts_raw <- as.data.frame(table(cellage_raw_standard[["Senescence Effect"]]), stringsAsFactors = FALSE)
colnames(effect_counts_raw) <- c("Senescence_Effect", "n_raw_rows")
effect_counts_clean <- as.data.frame(table(cellage_clean[["Senescence Effect"]]), stringsAsFactors = FALSE)
colnames(effect_counts_clean) <- c("Senescence_Effect", "n_clean_unique_symbols")

effect_counts <- merge(effect_counts_raw, effect_counts_clean, by = "Senescence_Effect", all = TRUE)
effect_counts$n_raw_rows[is.na(effect_counts$n_raw_rows)] <- 0
effect_counts$n_clean_unique_symbols[is.na(effect_counts$n_clean_unique_symbols)] <- 0
write_csv(effect_counts, file.path(table_dir, "step11A_CellAge_Senescence_Effect_counts.csv"))

selected_info <- file.info(cellage_file)

summary_df <- data.frame(
  metric = c(
    "selected_CellAge_file",
    "selected_file_size_bytes",
    "selected_file_ctime",
    "selected_file_mtime",
    "selected_file_md5",
    "inferred_CellAge_version_from_filename",
    "raw_rows",
    "raw_columns",
    "rows_with_nonmissing_Gene_symbol",
    "duplicate_Gene_symbol_rows",
    "clean_unique_CellAge_symbols",
    "cleaning_rule",
    "output_clean_table"
  ),
  value = c(
    cellage_file,
    selected_info$size,
    as.character(selected_info$ctime),
    as.character(selected_info$mtime),
    unname(tools::md5sum(cellage_file)),
    tools::file_path_sans_ext(basename(cellage_file)),
    nrow(cellage_raw_standard),
    ncol(cellage_raw_standard),
    sum(!is.na(cellage_raw_standard[["Gene symbol"]]) & cellage_raw_standard[["Gene symbol"]] != ""),
    nrow(duplicate_rows),
    nrow(cellage_clean),
    "remove missing Gene symbol; create CellAge_symbol; keep first occurrence for duplicated symbols; order by CellAge_symbol",
    file.path(table_dir, "step11A_CellAge_clean_271903.csv")
  ),
  stringsAsFactors = FALSE
)

write_csv(summary_df, file.path(table_dir, "step11A_CellAge_cleaning_summary.csv"))

versions <- data.frame(
  item = c("R"),
  version = c(as.character(getRversion())),
  stringsAsFactors = FALSE
)
if (requireNamespace("readxl", quietly = TRUE)) {
  versions <- rbind(
    versions,
    data.frame(item = "readxl", version = as.character(utils::packageVersion("readxl")), stringsAsFactors = FALSE)
  )
}
write_csv(versions, file.path(table_dir, "step11A_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step11A_sessionInfo.txt"))

cat("\n===== STEP11A SUMMARY =====\n")
print(summary_df)

cat("\nSenescence Effect counts:\n")
print(effect_counts)

cat("\nFirst rows of clean CellAge table:\n")
print(head(cellage_clean, 10))

cat("\nStep11A completed successfully.\n")
cat("Clean CellAge table:\n", file.path(table_dir, "step11A_CellAge_clean_271903.csv"), "\n", sep = "")

sink()

cat("\nStep11A CellAge raw-to-clean completed.\n")
cat("Clean CellAge table:\n", file.path(table_dir, "step11A_CellAge_clean_271903.csv"), "\n", sep = "")
