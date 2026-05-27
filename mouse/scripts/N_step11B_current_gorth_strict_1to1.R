# Step11B: Current gprofiler2::gorth mouse-to-human strict 1:1 ortholog remapping without frozen workspace
# English note:
# This script reruns the current online gprofiler2::gorth mapping for the locked Step07
# persistent direction-consistent mouse genes. It does NOT read or depend on the old
# Figure3C workspace. It applies a strict one-to-one mouse-human ortholog filtering rule,
# saves the mapping outputs, and records Cell/input provenance plus gprofiler2/g:Profiler
# version information for reproducibility.
#
# Locked strict 1:1 rule:
# Keep a mouse-human ortholog pair only if:
# 1) the mouse input gene maps to exactly one human ortholog in the gorth output; and
# 2) the human ortholog maps back to exactly one mouse input within the same gorth output.
#
# Main outputs:
# - current gorth raw output
# - all non-missing unique mouse-human pairs
# - strict 1:1 mouse-human ortholog mapping table
# - compact strict 1:1 mapping table
# - unmapped mouse input list
# - input file audit, gorth parameters, gprofiler2/g:Profiler version information

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

persistent_file <- file.path(
  base_dir, "07_tables", "step07_strict_DEG_upset_persistent",
  "step07_persistent_direction_consistent_genes.csv"
)

table_dir <- file.path(base_dir, "07_tables", "step11B_current_gorth_strict_1to1_remap")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step11B_current_gorth_strict_1to1_remap"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  log_file
)
unlink(old_outputs[file.exists(old_outputs)], recursive = TRUE, force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP11B CURRENT GORTH STRICT 1:1 REMAPPING =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Base directory: ", base_dir, "\n", sep = "")
cat("Persistent input file:\n", persistent_file, "\n\n", sep = "")
cat("Important: this script does NOT read the old Figure3C workspace.\n")
cat("It reruns current online gprofiler2::gorth directly.\n\n")

## =========================
## 1. Helper functions
## =========================

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
  x[x %in% c("", "NA", "NaN", "NULL", "null", "-", ".")] <- NA_character_
  x
}

pick_col <- function(df, candidates, required = TRUE, what = NULL) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Cannot find required column",
      if (!is.null(what)) paste0(" for ", what) else "",
      ". Checked: ", paste(candidates, collapse = ", ")
    )
  }
  NA_character_
}

file_audit <- function(path, label) {
  if (!file.exists(path)) {
    return(data.frame(
      label = label,
      path = path,
      exists = FALSE,
      size_bytes = NA_real_,
      modified_time = NA_character_,
      md5 = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  info <- file.info(path)
  data.frame(
    label = label,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    exists = TRUE,
    size_bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    md5 = unname(tools::md5sum(path)),
    stringsAsFactors = FALSE
  )
}

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

## =========================
## 2. Load packages and input
## =========================

safe_library("gprofiler2")

if (!file.exists(persistent_file)) {
  stop("Missing Step07 persistent table: ", persistent_file)
}

input_audit <- file_audit(persistent_file, "Step07 persistent direction-consistent genes")
write_csv(input_audit, file.path(table_dir, "step11B_input_file_audit.csv"))
cat("\nInput file audit:\n")
print(input_audit)

persistent <- read.csv(persistent_file, stringsAsFactors = FALSE, check.names = FALSE)

symbol_col <- pick_col(
  persistent,
  c("SYMBOL", "mouse_symbol", "symbol", "gene_symbol", "label_final"),
  what = "persistent mouse symbol"
)

persistent$mouse_symbol_for_gorth <- clean_symbol(persistent[[symbol_col]])
persistent_all_n <- nrow(persistent)

persistent_with_symbol <- persistent[
  !is.na(persistent$mouse_symbol_for_gorth) & persistent$mouse_symbol_for_gorth != "",
  ,
  drop = FALSE
]

persistent_with_symbol <- persistent_with_symbol[
  !duplicated(persistent_with_symbol$mouse_symbol_for_gorth),
  ,
  drop = FALSE
]

write_csv(persistent, file.path(table_dir, "step11B_input_persistent_mouse_genes_all.csv"))
write_csv(
  persistent_with_symbol,
  file.path(table_dir, "step11B_input_persistent_mouse_genes_with_SYMBOL_for_gorth.csv")
)

mouse_inputs <- persistent_with_symbol$mouse_symbol_for_gorth

if (length(mouse_inputs) == 0) {
  stop("No mouse symbols available for gorth after cleaning the persistent input table.")
}

cat("\nPersistent rows total: ", persistent_all_n, "\n", sep = "")
cat("Symbol column used: ", symbol_col, "\n", sep = "")
cat("Unique mouse symbols submitted to current gorth: ", length(mouse_inputs), "\n\n", sep = "")

## =========================
## 3. Version and mapping-parameter records
## =========================

gorth_parameters <- data.frame(
  parameter = c(
    "mapping_tool",
    "function",
    "source_organism",
    "target_organism",
    "mthreshold",
    "filter_na",
    "numeric_ns",
    "strict_1to1_rule"
  ),
  value = c(
    "gprofiler2 / g:Profiler",
    "gprofiler2::gorth",
    "mmusculus",
    "hsapiens",
    "Inf",
    "FALSE",
    "empty string",
    "n_human_orthologs == 1 and n_mouse_inputs == 1 within non-missing unique gorth pairs"
  ),
  stringsAsFactors = FALSE
)
write_csv(gorth_parameters, file.path(table_dir, "step11B_gorth_parameters.csv"))

base_url_value <- tryCatch(
  gprofiler2::get_base_url(),
  error = function(e) paste0("ERROR: ", conditionMessage(e))
)

basic_version_record <- data.frame(
  item = c("R", "platform", "gprofiler2_package", "gprofiler_base_url"),
  value = c(
    as.character(getRversion()),
    R.version$platform,
    as.character(utils::packageVersion("gprofiler2")),
    base_url_value
  ),
  stringsAsFactors = FALSE
)
write_csv(basic_version_record, file.path(table_dir, "step11B_current_gprofiler_basic_version_record.csv"))

version_info_mmusculus <- tryCatch(
  gprofiler2::get_version_info(organism = "mmusculus"),
  error = function(e) data.frame(error = conditionMessage(e), stringsAsFactors = FALSE)
)
version_info_hsapiens <- tryCatch(
  gprofiler2::get_version_info(organism = "hsapiens"),
  error = function(e) data.frame(error = conditionMessage(e), stringsAsFactors = FALSE)
)

write_csv(version_info_mmusculus, file.path(table_dir, "step11B_current_gprofiler_version_info_mmusculus.csv"))
write_csv(version_info_hsapiens, file.path(table_dir, "step11B_current_gprofiler_version_info_hsapiens.csv"))

cat("\nBasic version record:\n")
print(basic_version_record)
cat("\ngorth parameters:\n")
print(gorth_parameters)
cat("\ng:Profiler version info for mmusculus:\n")
print(version_info_mmusculus)
cat("\ng:Profiler version info for hsapiens:\n")
print(version_info_hsapiens)

## =========================
## 4. Rerun current online gorth
## =========================

cat("\nRunning current gprofiler2::gorth mouse -> human...\n")

gorth_raw <- gprofiler2::gorth(
  query = mouse_inputs,
  source_organism = "mmusculus",
  target_organism = "hsapiens",
  mthreshold = Inf,
  filter_na = FALSE,
  numeric_ns = ""
)

gorth_raw <- as.data.frame(gorth_raw, stringsAsFactors = FALSE)
write_csv(gorth_raw, file.path(table_dir, "step11B_current_gorth_mouse_to_human_raw.csv"))

required_gorth_cols <- c("input", "ortholog_name")
missing_gorth_cols <- setdiff(required_gorth_cols, colnames(gorth_raw))
if (length(missing_gorth_cols) > 0) {
  stop("Current gorth output is missing required columns: ", paste(missing_gorth_cols, collapse = ", "))
}

## =========================
## 5. Clean mapping and apply strict 1:1 filter
## =========================

mapping_all <- gorth_raw
mapping_all$input <- clean_symbol(mapping_all$input)
mapping_all$ortholog_name <- clean_symbol(mapping_all$ortholog_name)

unmapped_inputs <- setdiff(mouse_inputs, unique(mapping_all$input[!is.na(mapping_all$ortholog_name)]))
unmapped_df <- data.frame(mouse_symbol = unmapped_inputs, stringsAsFactors = FALSE)
write_csv(unmapped_df, file.path(table_dir, "step11B_current_gorth_unmapped_mouse_inputs.csv"))

mapping_all <- mapping_all[
  !is.na(mapping_all$input) & mapping_all$input != "" &
    !is.na(mapping_all$ortholog_name) & mapping_all$ortholog_name != "",
  ,
  drop = FALSE
]

mapping_all <- mapping_all[!duplicated(mapping_all[, c("input", "ortholog_name")]), , drop = FALSE]

n_human_per_mouse <- as.data.frame(table(mapping_all$input), stringsAsFactors = FALSE)
colnames(n_human_per_mouse) <- c("input", "n_human_orthologs")

n_mouse_per_human <- as.data.frame(table(mapping_all$ortholog_name), stringsAsFactors = FALSE)
colnames(n_mouse_per_human) <- c("ortholog_name", "n_mouse_inputs")

mapping_all <- merge(mapping_all, n_human_per_mouse, by = "input", all.x = TRUE, sort = FALSE)
mapping_all <- merge(mapping_all, n_mouse_per_human, by = "ortholog_name", all.x = TRUE, sort = FALSE)

mapping_all$n_human_orthologs <- as.integer(mapping_all$n_human_orthologs)
mapping_all$n_mouse_inputs <- as.integer(mapping_all$n_mouse_inputs)

mapping_all <- mapping_all[order(mapping_all$input, mapping_all$ortholog_name), , drop = FALSE]

mapping_strict <- mapping_all[
  mapping_all$n_human_orthologs == 1 & mapping_all$n_mouse_inputs == 1,
  ,
  drop = FALSE
]
mapping_strict <- mapping_strict[order(mapping_strict$input, mapping_strict$ortholog_name), , drop = FALSE]

compact_cols <- intersect(
  c("input", "ortholog_name", "input_ensg", "ortholog_ensg", "input_ensp", "ortholog_ensp"),
  colnames(mapping_strict)
)
if (length(compact_cols) < 2) compact_cols <- c("input", "ortholog_name")

mapping_strict_compact <- mapping_strict[, compact_cols, drop = FALSE]

write_csv(mapping_all, file.path(table_dir, "step11B_current_mouse_human_ortholog_mapping_all_unique_pairs.csv"))
write_csv(mapping_strict, file.path(table_dir, "step11B_current_mouse_human_ortholog_mapping_strict_1to1.csv"))
write_csv(mapping_strict_compact, file.path(table_dir, "step11B_current_mouse_human_ortholog_mapping_strict_1to1_compact.csv"))

## =========================
## 6. Summary, software versions, and script archive
## =========================

summary_df <- data.frame(
  metric = c(
    "analysis_role",
    "frozen_workspace_used",
    "persistent_input_file",
    "persistent_rows_total",
    "symbol_column_used",
    "unique_mouse_symbols_submitted_to_gorth",
    "current_gorth_raw_rows",
    "current_unique_mouse_human_pairs_nonmissing",
    "current_strict_1to1_pairs",
    "current_strict_unique_mouse_inputs",
    "current_strict_unique_human_symbols",
    "unmapped_mouse_inputs",
    "mapping_tool",
    "gprofiler2_package_version",
    "gprofiler_base_url"
  ),
  value = c(
    "current remapping; no frozen Figure3C workspace was used",
    "FALSE",
    normalizePath(persistent_file, winslash = "/", mustWork = FALSE),
    persistent_all_n,
    symbol_col,
    length(mouse_inputs),
    nrow(gorth_raw),
    nrow(mapping_all),
    nrow(mapping_strict),
    length(unique(mapping_strict$input)),
    length(unique(mapping_strict$ortholog_name)),
    length(unmapped_inputs),
    "gprofiler2::gorth",
    as.character(utils::packageVersion("gprofiler2")),
    base_url_value
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_df, file.path(table_dir, "step11B_current_gorth_strict_1to1_remap_summary.csv"))

software_versions <- data.frame(
  item = c("R", "platform", "gprofiler2"),
  version = c(
    as.character(getRversion()),
    R.version$platform,
    as.character(utils::packageVersion("gprofiler2"))
  ),
  stringsAsFactors = FALSE
)
write_csv(software_versions, file.path(table_dir, "step11B_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step11B_sessionInfo.txt"))

script_path <- get_script_path()
archive_path <- file.path(script_dir, "step11B_current_gorth_strict_1to1_remap_no_frozen_workspace.R")
if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path, warn = FALSE), archive_path)
} else {
  writeLines("# Please manually save this script here for reproducibility.", archive_path)
}
cat("Saved script archive placeholder/copy: ", archive_path, "\n", sep = "")

cat("\n===== STEP11B SUMMARY =====\n")
print(summary_df)

cat("\nFirst rows of strict 1:1 current mapping:\n")
print(head(mapping_strict, 20))

cat("\nUnmapped mouse inputs:\n")
print(unmapped_df)

cat("\nStep11B current gorth strict 1:1 remapping completed successfully.\n")

sink()

cat("\nStep11B current gorth strict 1:1 remapping completed.\n")
cat("Strict mapping:\n", file.path(table_dir, "step11B_current_mouse_human_ortholog_mapping_strict_1to1.csv"), "\n", sep = "")
cat("Summary:\n", file.path(table_dir, "step11B_current_gorth_strict_1to1_remap_summary.csv"), "\n", sep = "")
