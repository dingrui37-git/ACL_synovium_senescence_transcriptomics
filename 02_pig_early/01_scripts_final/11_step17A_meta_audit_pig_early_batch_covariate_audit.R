# ============================================================
# Step META-AUDIT: audit E-MTAB-6664 pig-early metadata for
# batch/covariate fields before edgeR differential expression
#
# Purpose:
#   This script audits the public E-MTAB-6664 IDF/SDRF metadata and
#   derived pig-early core metadata files to determine which metadata
#   fields are available for the 18 core validation samples and whether
#   any batch-like variables (batch, library preparation batch,
#   sequencing lane, processing date, flowcell, instrument, etc.) can be
#   consistently and stably included in the differential-expression
#   design matrix.
#
# Main inputs:
#   E:/R/ACLsenescence2/rebuild_submission/raw data/E-MTAB-6664
#
# Main outputs:
#   E:/R/ACLsenescence2/rebuild_submission/02_pig_early/early pig metadata batch covariate audit
#
# Notes:
#   - This audit does not change expression data or DE results.
#   - It produces transparent tables for reviewer-facing methodological
#     justification of the early-pig DE design.
#   - The 18 core samples are expected to be CON1-CON6, INJS1-INJS6,
#     and INJL1-INJL6.
# ============================================================

options(stringsAsFactors = FALSE)

project_root <- "E:/R/ACLsenescence2"
input_dir <- file.path(project_root, "rebuild_submission", "raw data", "E-MTAB-6664")
out_dir <- file.path(project_root, "rebuild_submission", "02_pig_early", "early pig metadata batch covariate audit")

tables_dir <- file.path(out_dir, "tables")
logs_dir <- file.path(out_dir, "logs")
objects_dir <- file.path(out_dir, "objects")
scripts_dir <- file.path(out_dir, "scripts")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step_meta_audit_pig_early_batch_covariate_audit"
log_file <- file.path(logs_dir, paste0(step_name, "_log.txt"))
script_archive <- file.path(scripts_dir, paste0(step_name, ".R"))

log_lines <- character(0)
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg))
}

write_log <- function() {
  writeLines(log_lines, con = log_file, useBytes = TRUE)
}

cleanup_outputs <- function() {
  files_to_remove <- c(
    file.path(tables_dir, paste0(step_name, "_input_file_status.csv")),
    file.path(tables_dir, paste0(step_name, "_idf_key_lines.csv")),
    file.path(tables_dir, paste0(step_name, "_metadata_source_column_audit.csv")),
    file.path(tables_dir, paste0(step_name, "_core18_metadata_merged.csv")),
    file.path(tables_dir, paste0(step_name, "_core18_group_summary.csv")),
    file.path(tables_dir, paste0(step_name, "_core18_available_covariates.csv")),
    file.path(tables_dir, paste0(step_name, "_batch_like_fields_audit.csv")),
    file.path(tables_dir, paste0(step_name, "_core18_covariate_by_group_summary.csv")),
    file.path(tables_dir, paste0(step_name, "_modeling_decision_summary.csv")),
    file.path(objects_dir, paste0(step_name, "_workspace.RData")),
    log_file,
    script_archive
  )
  existing <- files_to_remove[file.exists(files_to_remove)]
  if (length(existing) > 0) file.remove(existing)
}

cleanup_outputs()

safe_read_delim <- function(path, header = TRUE) {
  if (!file.exists(path)) return(NULL)
  read.delim(
    path,
    header = header,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    fill = TRUE,
    fileEncoding = "UTF-8"
  )
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

normalize_empty <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "na", "n/a", "null", "NULL", "None", "none", "not available", "not applicable")] <- NA_character_
  x
}

first_existing_col <- function(df, candidates, regex = FALSE) {
  if (is.null(df)) return(NA_character_)
  cn <- colnames(df)
  if (!regex) {
    hit <- candidates[candidates %in% cn]
    if (length(hit) > 0) return(hit[1])
    hit2 <- cn[tolower(cn) %in% tolower(candidates)]
    if (length(hit2) > 0) return(hit2[1])
  } else {
    for (p in candidates) {
      hit <- cn[grepl(p, cn, ignore.case = TRUE)]
      if (length(hit) > 0) return(hit[1])
    }
  }
  NA_character_
}

compact_values <- function(x, max_values = 12, max_chars = 400) {
  x <- normalize_empty(x)
  vals <- unique(x[!is.na(x)])
  vals <- vals[order(vals)]
  if (length(vals) == 0) return(NA_character_)
  out <- paste(head(vals, max_values), collapse = " | ")
  if (length(vals) > max_values) out <- paste0(out, " | ...")
  if (nchar(out) > max_chars) out <- paste0(substr(out, 1, max_chars), "...")
  out
}

summarise_columns <- function(df, source_name) {
  if (is.null(df)) {
    return(data.frame(source = character(0), column_name = character(0), n_rows = integer(0)))
  }
  do.call(rbind, lapply(colnames(df), function(col) {
    x <- normalize_empty(df[[col]])
    data.frame(
      source = source_name,
      column_name = col,
      n_rows = nrow(df),
      n_nonmissing = sum(!is.na(x)),
      n_missing = sum(is.na(x)),
      missing_rate_pct = round(100 * mean(is.na(x)), 2),
      n_unique_nonmissing = length(unique(x[!is.na(x)])),
      example_values = compact_values(x),
      stringsAsFactors = FALSE
    )
  }))
}

prefix_metadata <- function(df, prefix, sample_col) {
  if (is.null(df) || is.na(sample_col) || !(sample_col %in% colnames(df))) return(NULL)
  df <- df[!duplicated(df[[sample_col]]), , drop = FALSE]
  df$sample_id_merge <- normalize_empty(df[[sample_col]])
  keep <- !is.na(df$sample_id_merge)
  df <- df[keep, , drop = FALSE]
  colnames(df) <- ifelse(
    colnames(df) == "sample_id_merge",
    "sample_id_merge",
    paste0(prefix, colnames(df))
  )
  df
}

infer_core_group <- function(sample_id) {
  ifelse(
    grepl("^CON[1-6]$", sample_id), "CON_t0",
    ifelse(grepl("^INJS[1-6]$", sample_id), "ACLT_untreated_t7",
           ifelse(grepl("^INJL[1-6]$", sample_id), "ACLT_untreated_t28", NA_character_))
  )
}

# ------------------------------------------------------------
# 1. Locate input files
# ------------------------------------------------------------
idf_file <- file.path(input_dir, "E-MTAB-6664.idf.txt")
sdrf_file <- file.path(input_dir, "E-MTAB-6664.sdrf.txt")
analysis_ready_file <- file.path(input_dir, "E-MTAB-6664_analysis_ready_metadata.csv")
core_manifest_wide_file <- file.path(input_dir, "E-MTAB-6664_core_validation_manifest_wide.csv")
core_manifest_long_file <- file.path(input_dir, "E-MTAB-6664_core_validation_manifest_long.csv")
core_validation_metadata_file <- file.path(input_dir, "E-MTAB-6664_core_validation_metadata.csv")
sample_level_metadata_file <- file.path(input_dir, "E-MTAB-6664_sample_level_metadata.csv")

input_files <- data.frame(
  label = c(
    "idf", "sdrf", "analysis_ready_metadata", "core_manifest_wide",
    "core_manifest_long", "core_validation_metadata", "sample_level_metadata"
  ),
  path = c(
    idf_file, sdrf_file, analysis_ready_file, core_manifest_wide_file,
    core_manifest_long_file, core_validation_metadata_file, sample_level_metadata_file
  ),
  exists = file.exists(c(
    idf_file, sdrf_file, analysis_ready_file, core_manifest_wide_file,
    core_manifest_long_file, core_validation_metadata_file, sample_level_metadata_file
  )),
  size_bytes = ifelse(
    file.exists(c(
      idf_file, sdrf_file, analysis_ready_file, core_manifest_wide_file,
      core_manifest_long_file, core_validation_metadata_file, sample_level_metadata_file
    )),
    file.info(c(
      idf_file, sdrf_file, analysis_ready_file, core_manifest_wide_file,
      core_manifest_long_file, core_validation_metadata_file, sample_level_metadata_file
    ))$size,
    NA
  ),
  stringsAsFactors = FALSE
)

write.csv(input_files, file.path(tables_dir, paste0(step_name, "_input_file_status.csv")), row.names = FALSE, fileEncoding = "UTF-8")

if (!file.exists(sdrf_file)) stop("Missing required SDRF file: ", sdrf_file)
if (!file.exists(idf_file)) stop("Missing required IDF file: ", idf_file)

log_msg("Input directory: ", input_dir)
log_msg("Output directory: ", out_dir)

# ------------------------------------------------------------
# 2. Read metadata files
# ------------------------------------------------------------
idf_lines <- readLines(idf_file, warn = FALSE, encoding = "UTF-8")
sdrf_df <- safe_read_delim(sdrf_file, header = TRUE)
analysis_ready_df <- safe_read_csv(analysis_ready_file)
core_manifest_wide_df <- safe_read_csv(core_manifest_wide_file)
core_manifest_long_df <- safe_read_csv(core_manifest_long_file)
core_validation_metadata_df <- safe_read_csv(core_validation_metadata_file)
sample_level_metadata_df <- safe_read_csv(sample_level_metadata_file)

# IDF key-line audit for batch/protocol/sequencing clues
idf_key_pattern <- paste(
  c("batch", "lane", "library", "preparation", "prep", "date", "flowcell", "flow cell",
    "instrument", "platform", "sequenc", "protocol", "center", "centre", "run", "ENA", "ArrayExpress"),
  collapse = "|"
)
idf_key_lines <- data.frame(
  line_number = which(grepl(idf_key_pattern, idf_lines, ignore.case = TRUE)),
  line_text = idf_lines[grepl(idf_key_pattern, idf_lines, ignore.case = TRUE)],
  stringsAsFactors = FALSE
)
write.csv(idf_key_lines, file.path(tables_dir, paste0(step_name, "_idf_key_lines.csv")), row.names = FALSE, fileEncoding = "UTF-8")

# Column-level audit across sources
column_audit <- do.call(rbind, list(
  summarise_columns(sdrf_df, "sdrf"),
  summarise_columns(analysis_ready_df, "analysis_ready_metadata"),
  summarise_columns(core_manifest_wide_df, "core_manifest_wide"),
  summarise_columns(core_manifest_long_df, "core_manifest_long"),
  summarise_columns(core_validation_metadata_df, "core_validation_metadata"),
  summarise_columns(sample_level_metadata_df, "sample_level_metadata")
))
write.csv(column_audit, file.path(tables_dir, paste0(step_name, "_metadata_source_column_audit.csv")), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------------------------------------------
# 3. Build merged core-18 metadata table
# ------------------------------------------------------------
# Prefer analysis-ready metadata as the core sample anchor; fall back to manifest.
if (!is.null(analysis_ready_df) && "sample_id" %in% colnames(analysis_ready_df)) {
  ar <- analysis_ready_df
  ar$sample_id <- normalize_empty(ar$sample_id)
  if ("core_validation_use" %in% colnames(ar)) {
    core_flag <- tolower(as.character(ar$core_validation_use)) %in% c("true", "t", "1", "yes")
    core_base <- ar[core_flag, , drop = FALSE]
  } else if ("core_group" %in% colnames(ar)) {
    core_base <- ar[!is.na(normalize_empty(ar$core_group)), , drop = FALSE]
  } else {
    core_base <- ar[grepl("^(CON[1-6]|INJS[1-6]|INJL[1-6])$", ar$sample_id), , drop = FALSE]
  }
  core_base$sample_id_merge <- core_base$sample_id
} else if (!is.null(core_manifest_wide_df)) {
  sample_col <- first_existing_col(core_manifest_wide_df, c("sample_id", "sample", "Assay Name"), regex = FALSE)
  if (is.na(sample_col)) {
    sample_col <- first_existing_col(core_manifest_wide_df, c("sample", "assay"), regex = TRUE)
  }
  if (is.na(sample_col)) stop("Cannot identify sample column in core manifest wide file.")
  core_base <- core_manifest_wide_df
  core_base$sample_id <- normalize_empty(core_base[[sample_col]])
  core_base <- core_base[grepl("^(CON[1-6]|INJS[1-6]|INJL[1-6])$", core_base$sample_id), , drop = FALSE]
  core_base$core_group <- infer_core_group(core_base$sample_id)
  core_base$sample_id_merge <- core_base$sample_id
} else {
  stop("Cannot build core metadata: neither analysis-ready metadata nor core manifest is available.")
}

core_base <- core_base[!duplicated(core_base$sample_id_merge), , drop = FALSE]
core_base$core_group <- if ("core_group" %in% colnames(core_base)) normalize_empty(core_base$core_group) else NA_character_
core_base$core_group[is.na(core_base$core_group)] <- infer_core_group(core_base$sample_id_merge[is.na(core_base$core_group)])

# Restrict to expected core IDs.
core_base <- core_base[grepl("^(CON[1-6]|INJS[1-6]|INJL[1-6])$", core_base$sample_id_merge), , drop = FALSE]

if (nrow(core_base) != 18) {
  warning("Expected 18 core samples, but core metadata contains ", nrow(core_base), " rows. Please inspect outputs.")
}

# Merge additional sources with prefixes.
sdrf_sample_col <- first_existing_col(sdrf_df, c("Assay Name", "Source Name", "Sample Name", "Comment[ENA_SAMPLE]"), regex = FALSE)
analysis_sample_col <- if (!is.null(analysis_ready_df)) first_existing_col(analysis_ready_df, c("sample_id", "Assay Name", "source_name"), regex = FALSE) else NA_character_
coremeta_sample_col <- if (!is.null(core_validation_metadata_df)) first_existing_col(core_validation_metadata_df, c("sample_id", "Assay Name", "source_name"), regex = FALSE) else NA_character_
samplelevel_sample_col <- if (!is.null(sample_level_metadata_df)) first_existing_col(sample_level_metadata_df, c("sample_id", "Assay Name", "source_name"), regex = FALSE) else NA_character_
manifest_wide_sample_col <- if (!is.null(core_manifest_wide_df)) first_existing_col(core_manifest_wide_df, c("sample_id", "sample", "Assay Name"), regex = FALSE) else NA_character_

merge_sources <- list(
  prefix_metadata(sdrf_df, "sdrf__", sdrf_sample_col),
  prefix_metadata(analysis_ready_df, "analysis__", analysis_sample_col),
  prefix_metadata(core_validation_metadata_df, "coremeta__", coremeta_sample_col),
  prefix_metadata(sample_level_metadata_df, "samplelevel__", samplelevel_sample_col),
  prefix_metadata(core_manifest_wide_df, "manifest__", manifest_wide_sample_col)
)

core_merged <- core_base
for (src in merge_sources) {
  if (!is.null(src)) {
    core_merged <- merge(core_merged, src, by = "sample_id_merge", all.x = TRUE, sort = FALSE)
  }
}

# Ensure a concise front section.
front_cols <- unique(c("sample_id_merge", "sample_id", "core_group"))
front_cols <- front_cols[front_cols %in% colnames(core_merged)]
other_cols <- setdiff(colnames(core_merged), front_cols)
core_merged <- core_merged[, c(front_cols, other_cols), drop = FALSE]

write.csv(core_merged, file.path(tables_dir, paste0(step_name, "_core18_metadata_merged.csv")), row.names = FALSE, fileEncoding = "UTF-8")

core_group_summary <- as.data.frame(table(core_merged$core_group), stringsAsFactors = FALSE)
colnames(core_group_summary) <- c("core_group", "n_samples")
write.csv(core_group_summary, file.path(tables_dir, paste0(step_name, "_core18_group_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------------------------------------------
# 4. Core-18 covariate audit
# ------------------------------------------------------------
metadata_keywords <- paste(
  c("batch", "lane", "library", "prep", "preparation", "date", "flowcell", "flow cell",
    "instrument", "platform", "sequenc", "center", "centre", "facility", "run",
    "accession", "experiment", "assay", "protocol", "file", "uri", "submitted",
    "sample", "source", "organism", "site", "injury", "treatment", "time", "sex", "age", "group"),
  collapse = "|"
)

batch_like_keywords <- paste(
  c("batch", "lane", "library", "prep", "preparation", "date", "flowcell", "flow cell",
    "instrument", "platform", "sequenc", "center", "centre", "facility", "run",
    "accession", "experiment", "assay", "protocol", "file", "uri", "submitted"),
  collapse = "|"
)

unique_per_sample <- function(x) {
  x <- normalize_empty(x)
  vals <- x[!is.na(x)]
  length(vals) > 0 && length(unique(vals)) >= length(vals)
}

is_confounded_with_group <- function(x, group) {
  x <- normalize_empty(x)
  group <- normalize_empty(group)
  ok <- !is.na(x) & !is.na(group)
  if (sum(ok) == 0) return(NA)
  x <- x[ok]
  group <- group[ok]
  vals <- unique(x)
  if (length(vals) <= 1) return(FALSE)
  all(vapply(vals, function(v) length(unique(group[x == v])) == 1, logical(1)))
}

field_decision <- function(col, x, group) {
  x_norm <- normalize_empty(x)
  n <- length(x_norm)
  n_nonmissing <- sum(!is.na(x_norm))
  n_unique <- length(unique(x_norm[!is.na(x_norm)]))
  missing_rate <- ifelse(n == 0, NA_real_, mean(is.na(x_norm)))
  col_l <- tolower(col)
  is_batch_like <- grepl(batch_like_keywords, col_l, ignore.case = TRUE)
  is_primary_design <- grepl("core_group|group|injury|treatment|time", col_l, ignore.case = TRUE)
  is_id_or_file <- grepl("sample|source|assay|run|accession|file|uri|submitted", col_l, ignore.case = TRUE)
  unique_sample <- unique_per_sample(x_norm)
  confounded <- is_confounded_with_group(x_norm, group)
  complete <- n_nonmissing == n
  
  reason <- character(0)
  usable_for_batch_model <- FALSE
  
  if (!is_batch_like) reason <- c(reason, "not_batch_like")
  if (!complete) reason <- c(reason, "incomplete_or_missing")
  if (n_unique <= 1) reason <- c(reason, "constant_or_single_level")
  if (unique_sample) reason <- c(reason, "unique_per_sample_or_file_level")
  if (is_id_or_file) reason <- c(reason, "identifier_or_file_run_field")
  if (is_primary_design) reason <- c(reason, "primary_group_or_biological_design_field")
  if (isTRUE(confounded)) reason <- c(reason, "confounded_with_core_group")
  
  if (is_batch_like && complete && n_unique > 1 && !unique_sample && !is_id_or_file && !is_primary_design && !isTRUE(confounded)) {
    usable_for_batch_model <- TRUE
    reason <- "potentially_modelable_batch_like_field"
  }
  
  data.frame(
    column_name = col,
    n_core_samples = n,
    n_nonmissing = n_nonmissing,
    n_missing = n - n_nonmissing,
    missing_rate_pct = round(100 * missing_rate, 2),
    n_unique_nonmissing = n_unique,
    is_batch_like_field = is_batch_like,
    is_primary_design_field = is_primary_design,
    is_identifier_or_file_field = is_id_or_file,
    unique_per_sample_or_file_level = unique_sample,
    confounded_with_core_group = confounded,
    usable_for_batch_model = usable_for_batch_model,
    decision_reason = paste(unique(reason), collapse = "; "),
    example_values = compact_values(x_norm),
    stringsAsFactors = FALSE
  )
}

covariate_audit <- do.call(rbind, lapply(colnames(core_merged), function(col) {
  field_decision(col, core_merged[[col]], core_merged$core_group)
}))

covariate_audit <- covariate_audit[order(
  !covariate_audit$is_batch_like_field,
  !covariate_audit$usable_for_batch_model,
  covariate_audit$column_name
), , drop = FALSE]

write.csv(covariate_audit, file.path(tables_dir, paste0(step_name, "_core18_available_covariates.csv")), row.names = FALSE, fileEncoding = "UTF-8")

batch_like_audit <- subset(covariate_audit, is_batch_like_field)
write.csv(batch_like_audit, file.path(tables_dir, paste0(step_name, "_batch_like_fields_audit.csv")), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------------------------------------------
# 5. Candidate field distribution by group
# ------------------------------------------------------------
candidate_cols <- unique(c(
  covariate_audit$column_name[covariate_audit$is_batch_like_field],
  covariate_audit$column_name[grepl("sex|age|organism|site|source|time|treatment|injury|group", covariate_audit$column_name, ignore.case = TRUE)]
))
candidate_cols <- candidate_cols[candidate_cols %in% colnames(core_merged)]

by_group_list <- list()
idx <- 1
for (col in candidate_cols) {
  x <- normalize_empty(core_merged[[col]])
  group <- normalize_empty(core_merged$core_group)
  tmp <- as.data.frame(table(core_group = group, value = x, useNA = "ifany"), stringsAsFactors = FALSE)
  tmp <- tmp[tmp$Freq > 0, , drop = FALSE]
  if (nrow(tmp) > 0) {
    tmp$column_name <- col
    by_group_list[[idx]] <- tmp[, c("column_name", "core_group", "value", "Freq")]
    idx <- idx + 1
  }
}

covariate_by_group <- if (length(by_group_list) > 0) {
  do.call(rbind, by_group_list)
} else {
  data.frame(column_name = character(0), core_group = character(0), value = character(0), Freq = integer(0))
}
write.csv(covariate_by_group, file.path(tables_dir, paste0(step_name, "_core18_covariate_by_group_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------------------------------------------
# 6. Modeling decision summary
# ------------------------------------------------------------
modelable_batch_fields <- batch_like_audit$column_name[batch_like_audit$usable_for_batch_model]
complete_batch_like_fields <- batch_like_audit$column_name[batch_like_audit$n_missing == 0]
nonmodelable_reasons <- if (nrow(batch_like_audit) > 0) {
  paste(unique(batch_like_audit$decision_reason), collapse = " | ")
} else {
  "No batch-like fields detected in merged core-18 metadata columns."
}

modeling_decision_summary <- data.frame(
  metric = c(
    "n_core_samples_audited",
    "core_group_distribution",
    "n_metadata_sources_found",
    "n_total_merged_metadata_columns",
    "n_batch_like_fields_detected",
    "n_complete_batch_like_fields",
    "n_modelable_batch_like_fields",
    "modelable_batch_like_fields",
    "overall_decision",
    "reviewer_ready_method_sentence"
  ),
  value = c(
    nrow(core_merged),
    paste(paste(core_group_summary$core_group, core_group_summary$n_samples, sep = "="), collapse = "; "),
    sum(input_files$exists),
    ncol(core_merged),
    nrow(batch_like_audit),
    length(complete_batch_like_fields),
    length(modelable_batch_fields),
    ifelse(length(modelable_batch_fields) == 0, "NONE", paste(modelable_batch_fields, collapse = "; ")),
    ifelse(
      length(modelable_batch_fields) == 0,
      "No batch-like metadata field was consistently defined and suitable for stable inclusion in the core-18 group-based DE design matrix.",
      "At least one potentially modelable batch-like field was detected; inspect batch_like_fields_audit.csv before finalizing the DE design."
    ),
    ifelse(
      length(modelable_batch_fields) == 0,
      "Available E-MTAB-6664 metadata fields were audited before modeling; no explicit batch, library-preparation batch, sequencing-lane, flowcell, or processing-date variable was consistently defined and suitable for inclusion across the 18 core pig-early samples, so the primary edgeR QLF model focused on the predefined group effect.",
      "Available E-MTAB-6664 metadata fields were audited before modeling; at least one potentially modelable batch-like field was detected and should be evaluated before finalizing the primary DE design."
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(modeling_decision_summary, file.path(tables_dir, paste0(step_name, "_modeling_decision_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")

save(
  input_files,
  idf_key_lines,
  sdrf_df,
  analysis_ready_df,
  core_manifest_wide_df,
  core_manifest_long_df,
  core_validation_metadata_df,
  sample_level_metadata_df,
  column_audit,
  core_merged,
  core_group_summary,
  covariate_audit,
  batch_like_audit,
  covariate_by_group,
  modeling_decision_summary,
  file = file.path(objects_dir, paste0(step_name, "_workspace.RData"))
)

# Archive this script as executed, when possible.
this_file <- tryCatch(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
if (!is.na(this_file) && file.exists(this_file)) {
  file.copy(this_file, script_archive, overwrite = TRUE)
} else {
  writeLines(c(
    "# Script archive note:",
    "# This script was run from an interactive console or a non-file source.",
    "# Please save the distributed script as step_meta_audit_pig_early_batch_covariate_audit.R."
  ), con = script_archive, useBytes = TRUE)
}

log_msg("Core samples audited: ", nrow(core_merged))
log_msg("Core group distribution: ", paste(paste(core_group_summary$core_group, core_group_summary$n_samples, sep = "="), collapse = "; "))
log_msg("Batch-like fields detected: ", nrow(batch_like_audit))
log_msg("Modelable batch-like fields: ", ifelse(length(modelable_batch_fields) == 0, "NONE", paste(modelable_batch_fields, collapse = "; ")))
log_msg("Main decision: ", modeling_decision_summary$value[modeling_decision_summary$metric == "overall_decision"])
log_msg("Main output table: ", file.path(tables_dir, paste0(step_name, "_modeling_decision_summary.csv")))
write_log()

cat("\n===== Metadata batch/covariate audit completed =====\n")
cat("Output folder:\n", out_dir, "\n", sep = "")
cat("Key tables:\n")
cat("- ", file.path(tables_dir, paste0(step_name, "_modeling_decision_summary.csv")), "\n", sep = "")
cat("- ", file.path(tables_dir, paste0(step_name, "_batch_like_fields_audit.csv")), "\n", sep = "")
cat("- ", file.path(tables_dir, paste0(step_name, "_core18_available_covariates.csv")), "\n", sep = "")
cat("- ", file.path(tables_dir, paste0(step_name, "_core18_covariate_by_group_summary.csv")), "\n", sep = "")
cat("\nReviewer-ready sentence:\n")
cat(modeling_decision_summary$value[modeling_decision_summary$metric == "reviewer_ready_method_sentence"], "\n")
