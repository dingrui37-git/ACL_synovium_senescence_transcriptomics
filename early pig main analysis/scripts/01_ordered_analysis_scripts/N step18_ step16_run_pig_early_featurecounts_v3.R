# ============================================================
# Step 16 v3: pig-early gene-level featureCounts
# Project root: E:/R/ACLsenescence2
#
# Fix:
#   - no longer fails if fc$stat does not contain a uniquely parsed
#     "Assigned" row in the expected format
#   - always exports the raw stat table for inspection
#   - if Assigned cannot be identified uniquely, assigned-related
#     columns are set to NA, but the main count outputs are still saved
# ============================================================

setwd("E:/R/ACLsenescence2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_dir <- file.path(rebuild_root, "02_pig_early")
scripts_dir <- file.path(pig_dir, "scripts")
objects_dir <- file.path(pig_dir, "objects")
tables_dir <- file.path(pig_dir, "tables")
logs_dir <- file.path(pig_dir, "logs")
bam_dir <- file.path(pig_dir, "bam", "step15J_rsubread_align")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("Rsubread", quietly = TRUE)) {
  stop("Rsubread is not installed. Please install/load it before running Step 16 v3.")
}

normalize_slash <- function(x) {
  gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
}

manifest_file <- file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv")
canon_file <- file.path(tables_dir, "step14_pig_early_canonical_reference_manifest.csv")
status_candidates <- c(
  file.path(tables_dir, "step15K_pig_early_alignment_status_table.csv"),
  file.path(tables_dir, "step15J_pig_early_alignment_status_table.csv")
)

needed <- c(manifest_file, canon_file)
if (!all(file.exists(needed))) {
  stop("Required Step 09/14 files are missing. Please run Step 09 and Step 14 first.")
}

existing_status <- status_candidates[file.exists(status_candidates)]
if (length(existing_status) == 0) {
  stop("Cannot find Step 15 alignment status table. Please complete Step 15 first.")
}

manifest_df <- read.csv(manifest_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
canon_df <- read.csv(canon_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
align_status_df <- read.csv(existing_status[1], stringsAsFactors = FALSE, fileEncoding = "UTF-8")

gtf_rel <- canon_df$relative_path[canon_df$item == "annotation_gtf"]
if (length(gtf_rel) != 1 || is.na(gtf_rel)) {
  stop("Canonical GTF was not found in Step 14 manifest.")
}
gtf_abs <- normalize_slash(file.path(project_root, gtf_rel))
if (!file.exists(gtf_abs)) {
  stop("Canonical GTF does not exist on disk: ", gtf_abs)
}

bam_manifest <- data.frame(
  sample_id = manifest_df$sample_id,
  core_group = manifest_df$core_group,
  bam_filename = paste0(manifest_df$sample_id, ".bam"),
  bam_absolute = normalize_slash(file.path(bam_dir, paste0(manifest_df$sample_id, ".bam"))),
  stringsAsFactors = FALSE
)

bam_manifest$bam_exists <- file.exists(bam_manifest$bam_absolute)
bam_manifest$bam_size_mb <- ifelse(
  bam_manifest$bam_exists,
  round(file.info(bam_manifest$bam_absolute)$size / 1024^2, 3),
  NA_real_
)

status_match <- match(bam_manifest$sample_id, align_status_df$sample_id)
bam_manifest$align_completed_flag <- ifelse(
  !is.na(status_match),
  align_status_df$completed[status_match],
  FALSE
)

if (!all(bam_manifest$bam_exists)) {
  missing_samples <- bam_manifest$sample_id[!bam_manifest$bam_exists]
  stop("Some BAM files are missing in step15J_rsubread_align: ", paste(missing_samples, collapse = ", "))
}

if (!all(as.logical(bam_manifest$align_completed_flag))) {
  incomplete_samples <- bam_manifest$sample_id[!as.logical(bam_manifest$align_completed_flag)]
  stop("Some samples are not marked complete in Step 15 status: ", paste(incomplete_samples, collapse = ", "))
}

output_full_table <- file.path(tables_dir, "step16_pig_early_featurecounts_full_table.csv")
output_count_matrix <- file.path(tables_dir, "step16_pig_early_gene_count_matrix.csv")
output_annotation <- file.path(tables_dir, "step16_pig_early_featurecounts_annotation.csv")
output_stats <- file.path(tables_dir, "step16_pig_early_featurecounts_summary_by_sample.csv")
output_stat_table <- file.path(tables_dir, "step16_pig_early_featurecounts_stat_table.csv")
output_run_summary <- file.path(tables_dir, "step16_pig_early_featurecounts_run_summary.csv")
output_manifest <- file.path(tables_dir, "step16_pig_early_bam_manifest_used.csv")
output_workspace <- file.path(objects_dir, "step16_pig_early_featurecounts_workspace.RData")
output_log <- file.path(logs_dir, "step16_pig_early_featurecounts_log.txt")

step16_targets <- c(
  output_full_table,
  output_count_matrix,
  output_annotation,
  output_stats,
  output_stat_table,
  output_run_summary,
  output_manifest,
  output_workspace,
  output_log
)

cleanup_outputs <- function() {
  existing <- step16_targets[file.exists(step16_targets)]
  if (length(existing) > 0) unlink(existing, force = TRUE)
}

cleanup_outputs()

start_time <- Sys.time()

fc <- tryCatch(
  {
    Rsubread::featureCounts(
      files = bam_manifest$bam_absolute,
      annot.ext = gtf_abs,
      isGTFAnnotationFile = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      useMetaFeatures = TRUE,
      isPairedEnd = TRUE,
      requireBothEndsMapped = TRUE,
      countReadPairs = TRUE,
      strandSpecific = 0,
      nthreads = 1,
      allowMultiOverlap = FALSE,
      countMultiMappingReads = FALSE,
      primaryOnly = FALSE,
      checkFragLength = FALSE
    )
  },
  error = function(e) {
    cleanup_outputs()
    stop("featureCounts failed. Error message: ", conditionMessage(e))
  }
)

end_time <- Sys.time()
elapsed_min <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)

count_mat <- fc$counts
colnames(count_mat) <- bam_manifest$sample_id

annotation_df <- fc$annotation
if ("GeneID" %in% colnames(annotation_df)) {
  gene_id_col <- "GeneID"
} else if ("Geneid" %in% colnames(annotation_df)) {
  gene_id_col <- "Geneid"
} else {
  gene_id_col <- colnames(annotation_df)[1]
}

full_table <- cbind(annotation_df, as.data.frame(count_mat, stringsAsFactors = FALSE))

count_matrix_df <- data.frame(
  gene_id = annotation_df[[gene_id_col]],
  as.data.frame(count_mat, stringsAsFactors = FALSE),
  stringsAsFactors = FALSE
)

# ---- robust fc$stat handling ----
stat_obj <- fc$stat

to_utf8 <- function(x) {
  y <- suppressWarnings(iconv(as.character(x), from = "", to = "UTF-8", sub = "byte"))
  y[is.na(y)] <- ""
  y
}

if (is.matrix(stat_obj)) {
  stat_table_df <- data.frame(
    category = rownames(stat_obj),
    as.data.frame(stat_obj, stringsAsFactors = FALSE, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  stat_table_df <- as.data.frame(stat_obj, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("category" %in% colnames(stat_table_df))) {
    if ("Status" %in% colnames(stat_table_df)) {
      names(stat_table_df)[names(stat_table_df) == "Status"] <- "category"
    } else if ("status" %in% colnames(stat_table_df)) {
      names(stat_table_df)[names(stat_table_df) == "status"] <- "category"
    } else if (!is.null(rownames(stat_obj)) && all(nzchar(rownames(stat_obj)))) {
      stat_table_df <- cbind(category = rownames(stat_obj), stat_table_df, stringsAsFactors = FALSE)
    } else {
      stat_table_df <- cbind(category = seq_len(nrow(stat_table_df)), stat_table_df, stringsAsFactors = FALSE)
    }
  }
}

stat_table_df$category <- to_utf8(stat_table_df$category)

sample_cols_in_stat <- intersect(bam_manifest$bam_filename, colnames(stat_table_df))
if (length(sample_cols_in_stat) != nrow(bam_manifest)) {
  alt_cols <- intersect(bam_manifest$sample_id, colnames(stat_table_df))
  if (length(alt_cols) == nrow(bam_manifest)) {
    sample_cols_in_stat <- alt_cols
  }
}

if (length(sample_cols_in_stat) != nrow(bam_manifest)) {
  cleanup_outputs()
  stop("Could not identify all sample columns in fc$stat for sample-level summary.")
}

sample_stat_numeric <- stat_table_df[, sample_cols_in_stat, drop = FALSE]
for (j in seq_len(ncol(sample_stat_numeric))) {
  sample_stat_numeric[[j]] <- suppressWarnings(as.numeric(sample_stat_numeric[[j]]))
}
sample_stat_numeric_mat <- as.matrix(sample_stat_numeric)
colnames(sample_stat_numeric_mat) <- sample_cols_in_stat

normalize_category <- function(x) {
  x <- tolower(trimws(to_utf8(x)))
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

category_norm <- normalize_category(stat_table_df$category)
assigned_idx <- which(category_norm == "assigned")

assigned_identified <- length(assigned_idx) == 1

if (assigned_identified) {
  assigned_vec <- as.numeric(sample_stat_numeric_mat[assigned_idx, ])
} else {
  assigned_vec <- rep(NA_real_, ncol(sample_stat_numeric_mat))
}

total_counted_records_vec <- colSums(sample_stat_numeric_mat, na.rm = TRUE)

if (identical(sample_cols_in_stat, bam_manifest$bam_filename)) {
  assigned_vec <- assigned_vec[match(bam_manifest$bam_filename, sample_cols_in_stat)]
  total_counted_records_vec <- total_counted_records_vec[match(bam_manifest$bam_filename, sample_cols_in_stat)]
} else {
  assigned_vec <- assigned_vec[match(bam_manifest$sample_id, sample_cols_in_stat)]
  total_counted_records_vec <- total_counted_records_vec[match(bam_manifest$sample_id, sample_cols_in_stat)]
}

sample_summary_df <- data.frame(
  sample_id = bam_manifest$sample_id,
  core_group = bam_manifest$core_group,
  bam_filename = bam_manifest$bam_filename,
  bam_size_mb = bam_manifest$bam_size_mb,
  assigned = assigned_vec,
  total_counted_records = total_counted_records_vec,
  stringsAsFactors = FALSE
)

sample_summary_df$assignment_rate_pct <- ifelse(
  is.na(sample_summary_df$assigned) | is.na(sample_summary_df$total_counted_records) | sample_summary_df$total_counted_records == 0,
  NA_real_,
  round(100 * sample_summary_df$assigned / sample_summary_df$total_counted_records, 3)
)

run_summary_df <- data.frame(
  metric = c(
    "project_root",
    "bam_dir",
    "gtf_relative",
    "n_bam_input",
    "n_genes_counted",
    "nthreads",
    "isPairedEnd",
    "countReadPairs",
    "requireBothEndsMapped",
    "GTF_feature_type",
    "GTF_attr_type",
    "elapsed_minutes",
    "assigned_row_identified_uniquely",
    "assigned_row_match_count"
  ),
  value = c(
    project_root,
    bam_dir,
    gtf_rel,
    nrow(bam_manifest),
    nrow(count_mat),
    1,
    TRUE,
    TRUE,
    TRUE,
    "exon",
    "gene_id",
    elapsed_min,
    assigned_identified,
    length(assigned_idx)
  ),
  stringsAsFactors = FALSE
)

write.csv(full_table, output_full_table, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(count_matrix_df, output_count_matrix, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(annotation_df, output_annotation, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sample_summary_df, output_stats, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(stat_table_df, output_stat_table, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(run_summary_df, output_run_summary, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bam_manifest, output_manifest, row.names = FALSE, fileEncoding = "UTF-8")

save(
  fc,
  bam_manifest,
  annotation_df,
  count_mat,
  full_table,
  count_matrix_df,
  stat_table_df,
  sample_summary_df,
  run_summary_df,
  file = output_workspace
)

log_lines <- c(
  "Step 16 v3 completed successfully.",
  paste("Project root:", project_root),
  paste("BAM dir:", bam_dir),
  paste("GTF:", gtf_abs),
  paste("n BAM input:", nrow(bam_manifest)),
  paste("n genes counted:", nrow(count_mat)),
  paste("Elapsed minutes:", elapsed_min),
  paste("Assigned row identified uniquely:", assigned_identified),
  paste("Assigned row match count:", length(assigned_idx))
)

writeLines(log_lines, output_log, useBytes = TRUE)

writeLines(
  c(
    "# Step 16 v3 run script archived automatically.",
    "setwd(\"E:/R/ACLsenescence2\")"
  ),
  con = file.path(scripts_dir, "step16_run_pig_early_featurecounts_v3.R"),
  useBytes = TRUE
)

writeLines(
  c(
    "# Step 16 v3 check script archived automatically.",
    "setwd(\"E:/R/ACLsenescence2\")"
  ),
  con = file.path(scripts_dir, "step16_check_pig_early_featurecounts_v3.R"),
  useBytes = TRUE
)

message("Step 16 v3 finished successfully.")
