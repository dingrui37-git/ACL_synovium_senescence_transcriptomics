# StepSTRAND_01: Early pig strandSpecific featureCounts comparison (v3 auto-GTF)
# This script reruns featureCounts on the 18 core early pig BAM files using strandSpecific = 0, 1, and 2.
# Purpose: provide a clean, auditable strandness/quantification-mode check supporting the final strandSpecific setting.
# Outputs: assignment-rate tables, strandSpecific mode summary, recommended mode, count-matrix consistency audit, figures, and summary log.

options(stringsAsFactors = FALSE)

method_version <- "2026-05-10_clean_strandSpecific_0_1_2_featureCounts_comparison_v4_fixed_gtf_path"
step_name <- "StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison"

project_root <- "E:/R/ACLsenescence2"
bam_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align"
locked_count_matrix_file <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv"
fixed_gtf_file <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf"
output_root <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison"

# Candidate paths. The previous v2 script assumed the .gtf.gz path below; v3 also searches common project folders.
manual_gtf_candidates <- c(
  "E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf",
  "E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz",
  "E:/R/ACLsenescence2/rebuild_submission/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf",
  "E:/R/ACLsenescence2/rebuild_submission/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz",
  "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf",
  "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf.gz",
  "E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf",
  "E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf.gz"
)

subdirs <- c("logs", "tables", "source_data", "figures", "objects", "scripts", "featureCounts_tmp")
for (d in file.path(output_root, subdirs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Clean previous outputs after directories are ready and before opening the log connection.
old_files <- list.files(output_root, pattern = "^StepSTRAND_01_", recursive = TRUE, full.names = TRUE)
if (length(old_files) > 0) suppressWarnings(file.remove(old_files))

log_file <- file.path(output_root, "logs", "StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison_summary_log.txt")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")
status <- "FAILED"

safe_finish <- function(status_value = status) {
  cat("\n============================================================\n")
  cat(step_name, "finished with status:", status_value, "\n")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Summary log saved to:\n", log_file, "\n", sep = "")
  cat("============================================================\n")
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}

fail_stop <- function(msg) {
  cat("\nERROR:\n", msg, "\n", sep = "")
  safe_finish("FAILED")
  stop(msg, call. = FALSE)
}

cat("============================================================\n")
cat(step_name, "\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Method version:", method_version, "\n")
cat("============================================================\n\n")

required_pkgs <- c("Rsubread", "dplyr", "ggplot2", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) fail_stop(paste("Missing required packages:", paste(missing_pkgs, collapse = ", ")))

suppressPackageStartupMessages({
  library(Rsubread)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})
cat("Required packages loaded.\n")
cat("Rsubread version:", as.character(packageVersion("Rsubread")), "\n")
cat("dplyr version:", as.character(packageVersion("dplyr")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")
cat("tidyr version:", as.character(packageVersion("tidyr")), "\n\n")

# Archive running script if called via Rscript --file.
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) > 0) {
  src_script <- sub("^--file=", "", file_arg[1])
  if (file.exists(src_script)) {
    dst_script <- file.path(output_root, "scripts", basename(src_script))
    try(file.copy(src_script, dst_script, overwrite = TRUE), silent = TRUE)
    cat("Script archived to:\n", dst_script, "\n\n", sep = "")
  }
} else {
  cat("Script archive note: --file argument was not detected. If running interactively, manually save this script in:\n")
  cat(file.path(output_root, "scripts"), "\n\n")
}

cat("Input BAM directory:\n", bam_dir, "\n\n", sep = "")
cat("Locked count matrix file:\n", locked_count_matrix_file, "\n\n", sep = "")
cat("Output root:\n", output_root, "\n\n", sep = "")

if (!dir.exists(bam_dir)) fail_stop(paste("BAM directory does not exist:", bam_dir))
if (!file.exists(locked_count_matrix_file)) fail_stop(paste("Locked count matrix file does not exist:", locked_count_matrix_file))

sample_ids <- c(paste0("CON", 1:6), paste0("INJS", 1:6), paste0("INJL", 1:6))
sample_meta <- data.frame(
  sample_id = sample_ids,
  group = factor(c(rep("Control", 6), rep("ACLT_1W", 6), rep("ACLT_4W", 6)), levels = c("Control", "ACLT_1W", "ACLT_4W")),
  time_point = c(rep("t0", 6), rep("t7", 6), rep("t28", 6)),
  sample_number = rep(1:6, 3),
  stringsAsFactors = FALSE
)
cat("Expected 18 core samples:\n")
print(sample_meta, row.names = FALSE)
cat("\nGroup counts:\n")
print(table(sample_meta$group))
cat("\n")

bam_paths <- file.path(bam_dir, paste0(sample_ids, ".bam"))
bam_check <- data.frame(sample_id = sample_ids, bam_path = bam_paths, bam_exists = file.exists(bam_paths), stringsAsFactors = FALSE)
cat("BAM file check:\n")
print(bam_check, row.names = FALSE)
cat("\n")
if (any(!bam_check$bam_exists)) {
  fail_stop(paste("Missing BAM files for:", paste(bam_check$sample_id[!bam_check$bam_exists], collapse = ", ")))
}

# Fixed GTF path --------------------------------------------------------------
# User-confirmed Ensembl release 115 Sus scrofa GTF used for this strandSpecific comparison.
cat("Fixed GTF file supplied by user:\n", fixed_gtf_file, "\n\n", sep = "")
if (!file.exists(fixed_gtf_file)) {
  fail_stop(paste("Fixed GTF file does not exist:", fixed_gtf_file))
}

gtf_candidates_df <- data.frame(
  candidate_path = fixed_gtf_file,
  basename = basename(fixed_gtf_file),
  priority_score = 999,
  selection_reason = "User-confirmed fixed GTF path",
  stringsAsFactors = FALSE
)
write.csv(gtf_candidates_df, file.path(output_root, "tables", "StepSTRAND_01_gtf_candidates_detected.csv"), row.names = FALSE)
cat("Selected fixed GTF candidate:\n")
print(gtf_candidates_df, row.names = FALSE)
cat("\n")

gtf_original <- fixed_gtf_file

# Decompress .gtf.gz to a local reproducible copy if needed; otherwise copy the uncompressed GTF.
copy_or_gunzip_gtf <- function(src, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.gz$", src, ignore.case = TRUE)) {
    out_name <- sub("\\.gz$", "", basename(src), ignore.case = TRUE)
    dest <- file.path(out_dir, out_name)
    cat("Selected GTF is gzipped. Creating uncompressed working copy for featureCounts:\n", dest, "\n", sep = "")
    con_in <- gzfile(src, "rb")
    con_out <- file(dest, "wb")
    on.exit({
      try(close(con_in), silent = TRUE)
      try(close(con_out), silent = TRUE)
    }, add = TRUE)
    repeat {
      chunk <- readBin(con_in, what = "raw", n = 16 * 1024 * 1024)
      if (length(chunk) == 0) break
      writeBin(chunk, con_out)
    }
    return(dest)
  } else {
    dest <- file.path(out_dir, basename(src))
    if (normalizePath(src, winslash = "/", mustWork = FALSE) != normalizePath(dest, winslash = "/", mustWork = FALSE)) {
      ok <- file.copy(src, dest, overwrite = TRUE)
      if (!ok) stop("Failed to copy GTF working file: ", src, " -> ", dest)
    }
    return(dest)
  }
}

gtf_working <- copy_or_gunzip_gtf(gtf_original, file.path(output_root, "source_data", "gtf_working_copy"))
if (!file.exists(gtf_working)) fail_stop(paste("Working GTF file could not be created:", gtf_working))
cat("Working GTF used by featureCounts:\n", gtf_working, "\n\n", sep = "")

# Read locked count matrix ----------------------------------------------------
locked_counts_df <- read.csv(locked_count_matrix_file, check.names = FALSE)
if (!"gene_id" %in% colnames(locked_counts_df)) fail_stop("Locked count matrix does not contain a 'gene_id' column.")
missing_count_samples <- setdiff(sample_ids, colnames(locked_counts_df))
if (length(missing_count_samples) > 0) fail_stop(paste("Locked count matrix missing samples:", paste(missing_count_samples, collapse = ", ")))
locked_count_colsum <- colSums(as.matrix(locked_counts_df[, sample_ids, drop = FALSE]), na.rm = TRUE)
cat("Locked count matrix dimensions:", nrow(locked_counts_df), "rows x", ncol(locked_counts_df), "columns\n")
cat("Locked count matrix sample column sums:\n")
print(data.frame(sample_id = names(locked_count_colsum), locked_assigned_count = as.numeric(locked_count_colsum)), row.names = FALSE)
cat("\n")

# Helper to parse featureCounts result ---------------------------------------
parse_fc_stat <- function(fc, mode) {
  stat_df <- as.data.frame(fc$stat, check.names = FALSE)
  # Rsubread usually returns first column "Status" and one column per input file.
  status_col <- if ("Status" %in% names(stat_df)) "Status" else names(stat_df)[1]
  names(stat_df)[names(stat_df) == status_col] <- "status"
  sample_cols <- setdiff(names(stat_df), "status")
  # Convert file paths or basenames to sample IDs.
  sample_col_map <- data.frame(
    source_column = sample_cols,
    sample_id = sub("\\.bam$", "", basename(sample_cols), ignore.case = TRUE),
    stringsAsFactors = FALSE
  )
  long_list <- lapply(sample_cols, function(sc) {
    data.frame(
      strandSpecific = mode,
      sample_id = sample_col_map$sample_id[sample_col_map$source_column == sc][1],
      source_column = sc,
      status = as.character(stat_df$status),
      count = suppressWarnings(as.numeric(stat_df[[sc]])),
      stringsAsFactors = FALSE
    )
  })
  stat_long <- bind_rows(long_list)
  stat_wide <- tidyr::pivot_wider(stat_long, id_cols = c(strandSpecific, sample_id, source_column), names_from = status, values_from = count, values_fill = 0)
  if (!"Assigned" %in% names(stat_wide)) stat_wide$Assigned <- NA_real_
  status_cols <- setdiff(names(stat_wide), c("strandSpecific", "sample_id", "source_column"))
  stat_wide$total_fragments_considered <- rowSums(stat_wide[, status_cols, drop = FALSE], na.rm = TRUE)
  stat_wide$assigned_fragments <- stat_wide$Assigned
  stat_wide$unassigned_fragments <- stat_wide$total_fragments_considered - stat_wide$assigned_fragments
  stat_wide$assignment_rate_pct <- ifelse(stat_wide$total_fragments_considered > 0, 100 * stat_wide$assigned_fragments / stat_wide$total_fragments_considered, NA_real_)
  stat_wide
}

# Run featureCounts for each strand mode -------------------------------------
all_by_sample <- list()
all_counts_colsum <- list()
all_errors <- list()

for (mode in c(0, 1, 2)) {
  cat("\n------------------------------------------------------------\n")
  cat("Running Rsubread::featureCounts with strandSpecific =", mode, "\n")
  cat("------------------------------------------------------------\n")
  fc_tmp_dir <- file.path(output_root, "featureCounts_tmp", paste0("strandSpecific_", mode))
  dir.create(fc_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  fc <- tryCatch({
    Rsubread::featureCounts(
      files = bam_paths,
      annot.ext = gtf_working,
      isGTFAnnotationFile = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "gene_id",
      useMetaFeatures = TRUE,
      isPairedEnd = TRUE,
      countReadPairs = TRUE,
      requireBothEndsMapped = TRUE,
      checkFragLength = FALSE,
      strandSpecific = mode,
      nthreads = max(1, min(4, parallel::detectCores(logical = TRUE) - 1))
    )
  }, error = function(e) {
    all_errors[[as.character(mode)]] <<- conditionMessage(e)
    NULL
  })
  if (is.null(fc)) {
    cat("featureCounts failed for strandSpecific =", mode, "\n")
    cat("Error:", all_errors[[as.character(mode)]], "\n")
    next
  }
  saveRDS(fc, file.path(output_root, "objects", paste0("StepSTRAND_01_featureCounts_result_strandSpecific_", mode, ".rds")))
  by_sample <- parse_fc_stat(fc, mode)
  by_sample <- sample_meta %>% left_join(by_sample, by = "sample_id")
  by_sample$parse_status <- ifelse(!is.na(by_sample$assigned_fragments) & !is.na(by_sample$total_fragments_considered), "PASS", "FAIL")
  all_by_sample[[as.character(mode)]] <- by_sample

  fc_counts <- as.data.frame(fc$counts, check.names = FALSE)
  # First annotation columns can vary; sample columns are usually bam paths/basenames.
  numeric_cols <- names(fc_counts)[vapply(fc_counts, is.numeric, logical(1))]
  # Use the last 18 numeric columns as counts if annotation has numeric Start/End/Length columns before sample counts.
  count_cols <- tail(numeric_cols, 18)
  if (length(count_cols) != 18) {
    cat("Warning: could not confidently identify 18 count columns for mode", mode, "\n")
  } else {
    count_col_sample_ids <- sub("\\.bam$", "", basename(count_cols), ignore.case = TRUE)
    colsum_df <- data.frame(
      strandSpecific = mode,
      source_column = count_cols,
      sample_id = count_col_sample_ids,
      assigned_from_fc_counts_colsum = as.numeric(colSums(fc_counts[, count_cols, drop = FALSE], na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
    all_counts_colsum[[as.character(mode)]] <- colsum_df
  }
  cat("Completed strandSpecific =", mode, "\n")
  print(by_sample[, c("sample_id", "group", "assigned_fragments", "total_fragments_considered", "assignment_rate_pct", "parse_status")], row.names = FALSE)
}

if (length(all_by_sample) == 0) {
  error_df <- data.frame(strandSpecific = names(all_errors), error_message = unlist(all_errors), stringsAsFactors = FALSE)
  write.csv(error_df, file.path(output_root, "tables", "StepSTRAND_01_featureCounts_errors.csv"), row.names = FALSE)
  fail_stop("featureCounts failed for all strandSpecific modes. See StepSTRAND_01_featureCounts_errors.csv.")
}

by_sample_all <- bind_rows(all_by_sample)
write.csv(by_sample_all, file.path(output_root, "tables", "StepSTRAND_01_strandSpecific_assignment_rate_by_sample.csv"), row.names = FALSE)

mode_summary <- by_sample_all %>%
  group_by(strandSpecific) %>%
  summarise(
    n_samples_total = n(),
    n_parse_pass = sum(parse_status == "PASS", na.rm = TRUE),
    mean_assigned_fragments = mean(assigned_fragments, na.rm = TRUE),
    median_assigned_fragments = median(assigned_fragments, na.rm = TRUE),
    mean_total_fragments_considered = mean(total_fragments_considered, na.rm = TRUE),
    mean_assignment_rate_pct = mean(assignment_rate_pct, na.rm = TRUE),
    median_assignment_rate_pct = median(assignment_rate_pct, na.rm = TRUE),
    min_assignment_rate_pct = min(assignment_rate_pct, na.rm = TRUE),
    max_assignment_rate_pct = max(assignment_rate_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(desc(mean_assignment_rate_pct))
write.csv(mode_summary, file.path(output_root, "tables", "StepSTRAND_01_strandSpecific_mode_summary.csv"), row.names = FALSE)
cat("\nStrand-specific mode summary:\n")
print(mode_summary)

recommended_mode <- mode_summary$strandSpecific[1]
recommended_rule <- "Mode with the highest mean assignment_rate_pct across the 18 core samples."
recommended_df <- data.frame(
  recommended_strandSpecific = recommended_mode,
  recommendation_rule = recommended_rule,
  method_note = "All modes were run using the same BAM files, GTF annotation and featureCounts parameters except for strandSpecific.",
  stringsAsFactors = FALSE
)
write.csv(recommended_df, file.path(output_root, "tables", "StepSTRAND_01_recommended_strandSpecific_mode.csv"), row.names = FALSE)
cat("\nRecommended strandSpecific mode:\n")
print(recommended_df, row.names = FALSE)

# Mode 0 consistency with locked count matrix --------------------------------
fc_counts_colsum_all <- if (length(all_counts_colsum) > 0) bind_rows(all_counts_colsum) else data.frame()
if (nrow(fc_counts_colsum_all) > 0) {
  write.csv(fc_counts_colsum_all, file.path(output_root, "tables", "StepSTRAND_01_featureCounts_counts_colsum_by_mode.csv"), row.names = FALSE)
}

mode0_counts <- fc_counts_colsum_all %>% filter(strandSpecific == 0)
if (nrow(mode0_counts) > 0) {
  mode0_audit <- sample_meta %>%
    left_join(mode0_counts[, c("sample_id", "assigned_from_fc_counts_colsum")], by = "sample_id") %>%
    mutate(
      locked_count_matrix_colsum = as.numeric(locked_count_colsum[sample_id]),
      difference_mode0_rerun_minus_locked = assigned_from_fc_counts_colsum - locked_count_matrix_colsum,
      mode0_rerun_matches_locked_count_matrix = difference_mode0_rerun_minus_locked == 0
    )
} else {
  mode0_assigned <- by_sample_all %>% filter(strandSpecific == 0) %>% select(sample_id, assigned_fragments)
  mode0_audit <- sample_meta %>%
    left_join(mode0_assigned, by = "sample_id") %>%
    mutate(
      assigned_from_fc_counts_colsum = assigned_fragments,
      locked_count_matrix_colsum = as.numeric(locked_count_colsum[sample_id]),
      difference_mode0_rerun_minus_locked = assigned_from_fc_counts_colsum - locked_count_matrix_colsum,
      mode0_rerun_matches_locked_count_matrix = difference_mode0_rerun_minus_locked == 0
    )
}
write.csv(mode0_audit, file.path(output_root, "tables", "StepSTRAND_01_strandSpecific0_rerun_vs_locked_count_matrix_audit.csv"), row.names = FALSE)
cat("\nstrandSpecific = 0 rerun vs locked count matrix audit:\n")
print(mode0_audit, row.names = FALSE)

mode0_match <- all(mode0_audit$mode0_rerun_matches_locked_count_matrix, na.rm = FALSE)

# Figures ---------------------------------------------------------------------
p <- ggplot(by_sample_all, aes(x = factor(strandSpecific), y = assignment_rate_pct)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, height = 0, size = 2.4, alpha = 0.85) +
  labs(
    title = "Early pig featureCounts strandSpecific mode comparison",
    x = "featureCounts strandSpecific mode",
    y = "Assignment rate (%)"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(output_root, "figures", "StepSTRAND_01_strandSpecific_assignment_rate_comparison.pdf"), p, width = 6.5, height = 4.8, useDingbats = FALSE)
ggsave(file.path(output_root, "figures", "StepSTRAND_01_strandSpecific_assignment_rate_comparison.png"), p, width = 6.5, height = 4.8, dpi = 300)

# Method parameters and methods text -----------------------------------------
method_params <- data.frame(
  parameter = c(
    "method_version", "bam_dir", "gtf_original_selected", "gtf_working_used", "locked_count_matrix_file",
    "featureCounts_function", "isGTFAnnotationFile", "GTF.featureType", "GTF.attrType", "useMetaFeatures",
    "isPairedEnd", "countReadPairs", "requireBothEndsMapped", "checkFragLength", "strandSpecific_modes_compared",
    "n_core_samples", "recommendation_rule"
  ),
  value = c(
    method_version, bam_dir, gtf_original, gtf_working, locked_count_matrix_file,
    "Rsubread::featureCounts", "TRUE", "exon", "gene_id", "TRUE",
    "TRUE", "TRUE", "TRUE", "FALSE", "0, 1, 2",
    length(sample_ids), recommended_rule
  ),
  stringsAsFactors = FALSE
)
write.csv(method_params, file.path(output_root, "tables", "StepSTRAND_01_method_parameters.csv"), row.names = FALSE)

methods_text <- paste0(
  "To assess the appropriate featureCounts strand-specific counting mode for the early pig RNA-seq dataset, ",
  "the 18 core BAM files were re-counted using Rsubread::featureCounts under three otherwise identical settings ",
  "with strandSpecific = 0, 1, or 2. Gene-level counts were assigned to exon features from the Sus scrofa ",
  "Sscrofa11.1 Ensembl release 115 GTF annotation using gene_id as the meta-feature attribute. ",
  "Paired-end fragments were counted with countReadPairs = TRUE and requireBothEndsMapped = TRUE. ",
  "For each mode, assignment rates were calculated as Assigned fragments divided by the sum of featureCounts Assigned and Unassigned categories. ",
  "The mode with the highest mean assignment rate across the 18 core samples was selected for downstream gene-level quantification."
)
writeLines(methods_text, file.path(output_root, "tables", "StepSTRAND_01_methods_text_strandSpecific_comparison.txt"))

# Final summary ---------------------------------------------------------------
get_mode_metric <- function(mode, column) {
  val <- mode_summary[mode_summary$strandSpecific == mode, column, drop = TRUE]
  if (length(val) == 0) return(NA)
  as.character(round(as.numeric(val[1]), 3))
}
interpretation_flag <- if (recommended_mode == 0) {
  "strandSpecific_0_has_highest_mean_assignment_rate"
} else {
  paste0("strandSpecific_", recommended_mode, "_has_highest_mean_assignment_rate")
}

final_summary <- data.frame(
  metric = c(
    "run_status", "method_version", "bam_files_present", "gtf_original_selected", "gtf_working_used",
    "strandSpecific0_mean_assignment_rate_pct", "strandSpecific1_mean_assignment_rate_pct", "strandSpecific2_mean_assignment_rate_pct",
    "recommended_strandSpecific", "interpretation_flag", "strandSpecific0_rerun_assigned_counts_match_locked_count_matrix",
    "output_root"
  ),
  value = c(
    "SUCCESS", method_version, sum(bam_check$bam_exists), gtf_original, gtf_working,
    get_mode_metric(0, "mean_assignment_rate_pct"), get_mode_metric(1, "mean_assignment_rate_pct"), get_mode_metric(2, "mean_assignment_rate_pct"),
    as.character(recommended_mode), interpretation_flag, as.character(mode0_match), output_root
  ),
  stringsAsFactors = FALSE
)
write.csv(final_summary, file.path(output_root, "tables", "StepSTRAND_01_final_summary_for_review.csv"), row.names = FALSE)
cat("\nFinal summary for review:\n")
print(final_summary, row.names = FALSE)

cat("\nKey output files:\n")
cat("1) Assignment rate by sample:\n", file.path(output_root, "tables", "StepSTRAND_01_strandSpecific_assignment_rate_by_sample.csv"), "\n", sep = "")
cat("2) Mode-level summary:\n", file.path(output_root, "tables", "StepSTRAND_01_strandSpecific_mode_summary.csv"), "\n", sep = "")
cat("3) Recommended mode:\n", file.path(output_root, "tables", "StepSTRAND_01_recommended_strandSpecific_mode.csv"), "\n", sep = "")
cat("4) strandSpecific = 0 rerun vs locked count matrix audit:\n", file.path(output_root, "tables", "StepSTRAND_01_strandSpecific0_rerun_vs_locked_count_matrix_audit.csv"), "\n", sep = "")
cat("5) GTF candidates detected:\n", file.path(output_root, "tables", "StepSTRAND_01_gtf_candidates_detected.csv"), "\n", sep = "")
cat("6) Assignment-rate figure:\n", file.path(output_root, "figures", "StepSTRAND_01_strandSpecific_assignment_rate_comparison.pdf"), "\n", sep = "")
cat("7) Final summary:\n", file.path(output_root, "tables", "StepSTRAND_01_final_summary_for_review.csv"), "\n", sep = "")

cat("\nSession information:\n")
print(sessionInfo())

status <- "SUCCESS"
safe_finish(status)
