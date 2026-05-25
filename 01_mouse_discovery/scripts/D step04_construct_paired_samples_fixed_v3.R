# Step04: Construct paired ACLR-Contra sample sets for mouse discovery.
# Purpose:
# This corrected Step04 script diagnoses and fixes sample-name mismatch between
# the GEO count matrix and the mouse metadata. It retains Rupture/ACLR mice at
# 1W and 4W, constructs expected paired samples (mouse_id + L/R), standardizes
# expression column names when needed, matches samples robustly, and saves clean
# expression matrices and annotation tables for downstream PCA and paired limma-voom.

## =========================
## 0. Define paths
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

meta_file <- file.path(base_dir, "01_metadata", "step03_metadata_raw.csv")
expr_file <- file.path(base_dir, "02_expression_matrix", "step03_expression_matrix_raw_counts.csv")

expr_dir <- file.path(base_dir, "02_expression_matrix")
meta_dir <- file.path(base_dir, "01_metadata")
table_dir <- file.path(base_dir, "07_tables")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step04_construct_paired_samples"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

## =========================
## 1. Clean old Step04 outputs
## =========================

old_files <- c(
  file.path(expr_dir, "step04_expr_1W.csv"),
  file.path(expr_dir, "step04_expr_4W.csv"),
  file.path(meta_dir, "step04_anno_1W.csv"),
  file.path(meta_dir, "step04_anno_4W.csv"),
  file.path(table_dir, "step04_sample_name_diagnosis.csv"),
  file.path(table_dir, "step04_unmatched_expression_columns.csv"),
  file.path(table_dir, "step04_sample_matching_summary.csv"),
  file.path(table_dir, "step04_group_summary.csv"),
  file.path(table_dir, "step04_sessionInfo.txt"),
  file.path(log_dir, "step04_construct_paired_samples_log.txt"),
  file.path(script_dir, "step04_construct_paired_samples.R")
)
unlink(old_files[file.exists(old_files)], force = TRUE)

## =========================
## 2. Start log
## =========================

sink(log_file, split = TRUE)

cat("===== STEP04 CONSTRUCT PAIRED MOUSE SAMPLES =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Base directory: ", base_dir, "\n", sep = "")
cat("Metadata file: ", meta_file, "\n", sep = "")
cat("Expression file: ", expr_file, "\n\n", sep = "")

## =========================
## 3. Helper functions
## =========================

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }
  return(NA_character_)
}

archive_current_script <- function(archive_path) {
  src <- get_script_path()
  if (!is.na(src) && file.exists(src)) {
    writeLines(readLines(src, warn = FALSE, encoding = "UTF-8"), archive_path, useBytes = TRUE)
    cat("Script archived to: ", archive_path, "\n", sep = "")
  } else {
    writeLines(
      c(
        "# Step04: Construct paired ACLR-Contra sample sets for mouse discovery.",
        "# Script archive fallback: R could not detect the executed script path."
      ),
      archive_path,
      useBytes = TRUE
    )
    cat("Script path not detected; archive placeholder saved to: ", archive_path, "\n", sep = "")
  }
}

find_col <- function(df, exact_candidates, regex_candidates = character(0)) {
  for (cn in exact_candidates) {
    if (cn %in% colnames(df)) return(cn)
  }
  for (pat in regex_candidates) {
    hit <- grep(pat, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

clean_sample_name <- function(x) {
  x0 <- as.character(x)
  x1 <- gsub("^X", "", x0)
  x1 <- gsub("\\.", "", x1)
  x1 <- gsub("_", "", x1)
  x1 <- gsub("-", "", x1)
  x1 <- toupper(x1)

  out <- rep(NA_character_, length(x1))

  for (i in seq_along(x1)) {
    m <- regexpr("M[0-9]+[LR]", x1[i], perl = TRUE)
    if (!is.na(m) && m > 0) {
      out[i] <- regmatches(x1[i], m)
    } else {
      fallback <- x1[i]
      if (is.na(fallback) || fallback == "") fallback <- paste0("COLUMN", i)
      out[i] <- paste0("UNMATCHED_", fallback)
    }
  }

  if (length(out) != length(x0)) {
    stop("Internal error: clean_sample_name returned wrong length.")
  }
  out
}

standardize_mouse_id <- function(x) {
  x <- as.character(x)
  x <- toupper(gsub("[^A-Za-z0-9]", "", x))
  out <- rep(NA_character_, length(x))

  for (i in seq_along(x)) {
    m <- regexpr("M[0-9]+", x[i], perl = TRUE)
    if (!is.na(m) && m > 0) {
      out[i] <- regmatches(x[i], m)
    } else {
      out[i] <- x[i]
    }
  }
  out
}

expand_pairs <- function(df, mouse_col, sex_col, timepoint) {
  out <- list()
  for (i in seq_len(nrow(df))) {
    mouse_id <- standardize_mouse_id(df[[mouse_col]][i])
    sex <- if (!is.na(sex_col)) as.character(df[[sex_col]][i]) else NA_character_

    out[[i]] <- data.frame(
      sample_id = c(paste0(mouse_id, "L"), paste0(mouse_id, "R")),
      mouse_id = mouse_id,
      sex = sex,
      timepoint = timepoint,
      treatment = c("Contra", "ACLR"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

match_expression <- function(anno, expr, timepoint) {
  expr_original_names <- colnames(expr)
  expr_clean_names <- clean_sample_name(expr_original_names)

  real_clean <- expr_clean_names[
    !is.na(expr_clean_names) &
      expr_clean_names != "" &
      !grepl("^UNMATCHED_", expr_clean_names)
  ]

  if (any(duplicated(real_clean))) {
    dup <- unique(real_clean[duplicated(real_clean)])
    stop(timepoint, ": duplicated cleaned expression sample names: ", paste(dup, collapse = ", "))
  }

  idx <- match(anno$sample_id, expr_clean_names)
  matched <- !is.na(idx)

  anno_keep <- anno[matched, , drop = FALSE]
  idx_keep <- idx[matched]

  expr_keep <- expr[, idx_keep, drop = FALSE]
  colnames(expr_keep) <- anno_keep$sample_id

  missing_samples <- anno$sample_id[!matched]
  matched_original_cols <- expr_original_names[idx_keep]

  list(
    anno = anno_keep,
    expr = expr_keep,
    missing_samples = missing_samples,
    matched_original_cols = matched_original_cols,
    expr_clean_names = expr_clean_names,
    expr_original_names = expr_original_names
  )
}

## =========================
## 4. Load Step03 data
## =========================

if (!file.exists(meta_file)) stop("Missing metadata file: ", meta_file)
if (!file.exists(expr_file)) stop("Missing expression matrix file: ", expr_file)

meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
expr <- read.csv(expr_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

cat("Loaded metadata dimensions: ", paste(dim(meta), collapse = " x "), "\n", sep = "")
cat("Loaded expression dimensions: ", paste(dim(expr), collapse = " x "), "\n", sep = "")

## =========================
## 5. Detect columns
## =========================

mouse_col <- find_col(
  meta,
  exact_candidates = c("Mouse.ID", "Mouse ID", "mouse_id", "Mouse", "...11"),
  regex_candidates = c("mouse", "identifier")
)
group_col <- find_col(
  meta,
  exact_candidates = c("Group", "group"),
  regex_candidates = c("group")
)
sex_col <- find_col(
  meta,
  exact_candidates = c("Sex", "sex"),
  regex_candidates = c("^sex$")
)

if (is.na(mouse_col)) stop("Could not detect mouse ID column.")
if (is.na(group_col)) stop("Could not detect group column.")

cat("Detected mouse column: ", mouse_col, "\n", sep = "")
cat("Detected group column: ", group_col, "\n", sep = "")
cat("Detected sex column: ", ifelse(is.na(sex_col), "NA", sex_col), "\n", sep = "")

## =========================
## 6. Keep Rupture 1W and Rupture 4W mice
## =========================

meta_rupture <- meta[grepl("Rupture", meta[[group_col]], ignore.case = TRUE), , drop = FALSE]

meta_1w <- meta_rupture[grepl("1W", meta_rupture[[group_col]], ignore.case = TRUE), , drop = FALSE]
meta_4w <- meta_rupture[grepl("4W", meta_rupture[[group_col]], ignore.case = TRUE), , drop = FALSE]

if (nrow(meta_1w) == 0) stop("No Rupture 1W mice found.")
if (nrow(meta_4w) == 0) stop("No Rupture 4W mice found.")

anno_1w_expected <- expand_pairs(meta_1w, mouse_col = mouse_col, sex_col = sex_col, timepoint = "1W")
anno_4w_expected <- expand_pairs(meta_4w, mouse_col = mouse_col, sex_col = sex_col, timepoint = "4W")

## =========================
## 7. Diagnose and match expression columns
## =========================

sample_diag <- data.frame(
  original_expr_colname = colnames(expr),
  cleaned_expr_colname = clean_sample_name(colnames(expr)),
  stringsAsFactors = FALSE
)
write.csv(sample_diag, file.path(table_dir, "step04_sample_name_diagnosis.csv"), row.names = FALSE)

unmatched_cols <- sample_diag[grepl("^UNMATCHED_", sample_diag$cleaned_expr_colname), , drop = FALSE]
write.csv(unmatched_cols, file.path(table_dir, "step04_unmatched_expression_columns.csv"), row.names = FALSE)

cat("\nFirst 30 expression column-name diagnosis rows:\n")
print(head(sample_diag, 30))

cat("\nNumber of expression columns not matching M[digits][L/R] pattern: ", nrow(unmatched_cols), "\n", sep = "")

cat("\nExpected 1W sample IDs:\n")
print(anno_1w_expected$sample_id)

cat("\nExpected 4W sample IDs:\n")
print(anno_4w_expected$sample_id)

res_1w <- match_expression(anno_1w_expected, expr, "1W")
res_4w <- match_expression(anno_4w_expected, expr, "4W")

if (ncol(res_1w$expr) == 0) {
  stop("No 1W paired samples matched expression matrix after sample-name cleaning. Check step04_sample_name_diagnosis.csv.")
}
if (ncol(res_4w$expr) == 0) {
  stop("No 4W paired samples matched expression matrix after sample-name cleaning. Check step04_sample_name_diagnosis.csv.")
}

## =========================
## 8. Save clean paired datasets
## =========================

write.csv(res_1w$expr, file.path(expr_dir, "step04_expr_1W.csv"), row.names = TRUE)
write.csv(res_4w$expr, file.path(expr_dir, "step04_expr_4W.csv"), row.names = TRUE)

write.csv(res_1w$anno, file.path(meta_dir, "step04_anno_1W.csv"), row.names = FALSE)
write.csv(res_4w$anno, file.path(meta_dir, "step04_anno_4W.csv"), row.names = FALSE)

matching_summary <- data.frame(
  timepoint = c("1W", "4W"),
  expected_mice = c(nrow(meta_1w), nrow(meta_4w)),
  expected_samples = c(nrow(anno_1w_expected), nrow(anno_4w_expected)),
  matched_samples = c(nrow(res_1w$anno), nrow(res_4w$anno)),
  matched_mice = c(length(unique(res_1w$anno$mouse_id)), length(unique(res_4w$anno$mouse_id))),
  missing_samples = c(
    paste(res_1w$missing_samples, collapse = ";"),
    paste(res_4w$missing_samples, collapse = ";")
  ),
  stringsAsFactors = FALSE
)

group_summary <- rbind(
  data.frame(timepoint = "1W", as.data.frame(table(res_1w$anno$treatment)), stringsAsFactors = FALSE),
  data.frame(timepoint = "4W", as.data.frame(table(res_4w$anno$treatment)), stringsAsFactors = FALSE)
)
colnames(group_summary) <- c("timepoint", "treatment", "n")

write.csv(matching_summary, file.path(table_dir, "step04_sample_matching_summary.csv"), row.names = FALSE)
write.csv(group_summary, file.path(table_dir, "step04_group_summary.csv"), row.names = FALSE)

## =========================
## 9. Save software/session info and archive script
## =========================

writeLines(capture.output(sessionInfo()), file.path(table_dir, "step04_sessionInfo.txt"))
archive_current_script(file.path(script_dir, "step04_construct_paired_samples.R"))

## =========================
## 10. Console summary
## =========================

cat("\n===== STEP04 SUMMARY =====\n")
cat("\nSample matching summary:\n")
print(matching_summary)

cat("\nGroup summary:\n")
print(group_summary)

cat("\n1W annotation preview:\n")
print(head(res_1w$anno, 12))

cat("\n4W annotation preview:\n")
print(head(res_4w$anno, 12))

cat("\nStep04 completed successfully.\n")

sink()
