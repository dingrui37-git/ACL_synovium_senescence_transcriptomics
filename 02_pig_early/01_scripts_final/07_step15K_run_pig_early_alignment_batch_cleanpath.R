# ============================================================
# Step 15K: pig-early batch low-memory alignment (auto-continue)
# Project root: E:/R/ACLsenescence2
#
# Purpose:
#   - continue from Step 15J v2 success state
#   - automatically process ALL remaining pending samples in one run
#   - keep nthreads = 1 for lower memory usage
#   - save status/summary/log AFTER EACH SAMPLE
#   - can be safely re-run; completed BAMs are skipped automatically
#
# Notes:
#   - This uses the SAME BAM directory as Step 15J v2:
#       rebuild_submission/02_pig_early/bam/step15J_rsubread_align
#   - If one sample fails, the script records the error and CONTINUES
#     to later samples. Review the status table afterward.
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
dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("Rsubread", quietly = TRUE)) {
  stop("Rsubread is not installed. Please install/load it before running Step 15K.")
}

normalize_slash <- function(x) {
  gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
}

manifest_file <- file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv")
canon_file <- file.path(tables_dir, "step14_pig_early_canonical_reference_manifest.csv")
decision_file <- file.path(tables_dir, "step14_pig_early_index_reuse_decision_summary.csv")

needed <- c(manifest_file, canon_file, decision_file)
if (!all(file.exists(needed))) {
  stop("Required Step 09/14 files are missing. Please run Step 09 and Step 14 first.")
}

fastq_manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
canon_df <- read.csv(canon_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
decision_df <- read.csv(decision_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

reuse_flag <- decision_df$value[decision_df$metric == "reuse_existing_index"]
if (length(reuse_flag) != 1 || !isTRUE(as.logical(reuse_flag))) {
  stop("Step 14 did not approve reuse_existing_index. Step 15K expects the reuse_existing_index branch.")
}

index_rel <- canon_df$relative_path[canon_df$item == "subread_index_basename"]
fasta_rel <- canon_df$relative_path[canon_df$item == "genome_fasta"]
gtf_rel   <- canon_df$relative_path[canon_df$item == "annotation_gtf"]

if (length(index_rel) != 1 || is.na(index_rel)) stop("Canonical Subread index basename not found in Step 14 manifest.")
if (length(fasta_rel) != 1 || is.na(fasta_rel)) stop("Canonical FASTA not found in Step 14 manifest.")
if (length(gtf_rel) != 1 || is.na(gtf_rel)) stop("Canonical GTF not found in Step 14 manifest.")

index_abs <- normalize_slash(file.path(project_root, index_rel))
fasta_abs <- normalize_slash(file.path(project_root, fasta_rel))
gtf_abs   <- normalize_slash(file.path(project_root, gtf_rel))

if (!file.exists(fasta_abs)) stop("Canonical FASTA not found on disk: ", fasta_abs)
if (!file.exists(gtf_abs)) stop("Canonical GTF not found on disk: ", gtf_abs)

idx_dir <- dirname(index_abs)
idx_prefix <- basename(index_abs)

if (!dir.exists(idx_dir)) {
  stop("Canonical Subread index directory not found on disk: ", idx_dir)
}

bundle_files <- list.files(idx_dir, full.names = TRUE)
bundle_files <- normalize_slash(bundle_files)
bundle_files <- bundle_files[file.exists(bundle_files)]
bundle_files <- bundle_files[startsWith(basename(bundle_files), idx_prefix)]

if (length(bundle_files) == 0) {
  stop("No Subread index bundle files were found for basename prefix: ", index_abs)
}

if (!any(grepl("\\.[0-9]{2}\\.[A-Za-z]\\.array$", basename(bundle_files)))) {
  stop("Canonical Subread index bundle does not look complete: ", index_abs)
}

run_manifest <- data.frame(
  sample_id = fastq_manifest$sample_id,
  core_group = fastq_manifest$core_group,
  fastq_r1_relative = fastq_manifest$fastq_r1_relative,
  fastq_r2_relative = fastq_manifest$fastq_r2_relative,
  fastq_r1_absolute = normalize_slash(file.path(project_root, fastq_manifest$fastq_r1_relative)),
  fastq_r2_absolute = normalize_slash(file.path(project_root, fastq_manifest$fastq_r2_relative)),
  bam_filename = paste0(fastq_manifest$sample_id, ".bam"),
  bam_absolute = normalize_slash(file.path(bam_dir, paste0(fastq_manifest$sample_id, ".bam"))),
  stringsAsFactors = FALSE
)

run_manifest$fastq_r1_exists <- file.exists(run_manifest$fastq_r1_absolute)
run_manifest$fastq_r2_exists <- file.exists(run_manifest$fastq_r2_absolute)

if (!all(run_manifest$fastq_r1_exists) || !all(run_manifest$fastq_r2_exists)) {
  stop("One or more FASTQ files listed in Step 09 manifest are missing under the new project root.")
}

status_file <- file.path(tables_dir, "step15K_pig_early_alignment_status_table.csv")
summary_file <- file.path(tables_dir, "step15K_pig_early_alignment_status_summary.csv")
workspace_file <- file.path(objects_dir, "step15K_pig_early_alignment_workspace.RData")
log_file <- file.path(logs_dir, "step15K_pig_early_alignment_log.txt")
manifest_out_file <- file.path(tables_dir, "step15K_pig_early_alignment_manifest.csv")

old_status <- data.frame()
candidate_old_status <- c(
  file.path(tables_dir, "step15K_pig_early_alignment_status_table.csv"),
  file.path(tables_dir, "step15J_pig_early_alignment_status_table.csv")
)
existing_old <- candidate_old_status[file.exists(candidate_old_status)]
if (length(existing_old) > 0) {
  old_status <- read.csv(existing_old[1], stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}

status_df <- data.frame(
  sample_id = run_manifest$sample_id,
  core_group = run_manifest$core_group,
  bam_exists = FALSE,
  bam_size_mb = NA_real_,
  completed = FALSE,
  run_attempts = 0L,
  last_attempt_time = NA_character_,
  last_success_time = NA_character_,
  last_error_message = NA_character_,
  stringsAsFactors = FALSE
)

if (nrow(old_status) > 0) {
  common_cols <- intersect(colnames(status_df), colnames(old_status))
  match_idx <- match(status_df$sample_id, old_status$sample_id)
  matched <- !is.na(match_idx)
  status_df[matched, common_cols] <- old_status[match_idx[matched], common_cols]
}

status_df$bam_exists <- file.exists(run_manifest$bam_absolute)
status_df$bam_size_mb <- ifelse(
  status_df$bam_exists,
  round(file.info(run_manifest$bam_absolute)$size / 1024^2, 3),
  NA_real_
)
status_df$completed <- status_df$bam_exists & !is.na(status_df$bam_size_mb) & status_df$bam_size_mb > 0

nthreads_each <- 1L
phred_offset_used <- 33L
batch_start_time <- Sys.time()

write_outputs <- function(summary_df, log_lines) {
  write.csv(run_manifest, manifest_out_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(status_df, status_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(summary_df, summary_file, row.names = FALSE, fileEncoding = "UTF-8")
  save(
    run_manifest,
    status_df,
    summary_df,
    index_abs,
    bundle_files,
    fasta_abs,
    gtf_abs,
    bam_dir,
    file = workspace_file
  )
  writeLines(log_lines, log_file, useBytes = TRUE)
  writeLines(
    c(
      "# Step 15K run script archived automatically.",
      "setwd(\"E:/R/ACLsenescence2\")"
    ),
    con = file.path(scripts_dir, "step15K_run_pig_early_alignment_batch_cleanpath.R"),
    useBytes = TRUE
  )
  writeLines(
    c(
      "# Step 15K check script archived automatically.",
      "setwd(\"E:/R/ACLsenescence2\")"
    ),
    con = file.path(scripts_dir, "step15K_check_pig_early_alignment_batch_cleanpath.R"),
    useBytes = TRUE
  )
}

make_summary <- function(current_target_sample = NA_character_,
                         current_target_group = NA_character_,
                         current_run_success = NA,
                         current_run_elapsed_minutes = NA_real_,
                         batch_elapsed_minutes = NA_real_) {
  remaining_idx <- which(!status_df$completed)
  next_pending_sample <- if (length(remaining_idx) > 0) status_df$sample_id[remaining_idx[1]] else NA_character_
  data.frame(
    metric = c(
      "mode", "project_root", "n_samples_input", "n_completed", "n_remaining",
      "current_target_sample", "current_target_group", "current_run_success",
      "current_run_elapsed_minutes", "next_pending_sample", "bam_dir",
      "nthreads_each_sample", "phred_offset_used", "batch_elapsed_minutes"
    ),
    value = c(
      "batch_resume_cleanpath",
      project_root,
      nrow(run_manifest),
      sum(status_df$completed),
      sum(!status_df$completed),
      current_target_sample,
      current_target_group,
      current_run_success,
      current_run_elapsed_minutes,
      next_pending_sample,
      bam_dir,
      nthreads_each,
      phred_offset_used,
      batch_elapsed_minutes
    ),
    stringsAsFactors = FALSE
  )
}

initial_pending <- which(!status_df$completed)
if (length(initial_pending) == 0) {
  summary_df <- make_summary(
    batch_elapsed_minutes = round(as.numeric(difftime(Sys.time(), batch_start_time, units = "mins")), 3)
  )
  log_lines <- c(
    "All Step 15K BAMs are already complete.",
    paste("Project root:", project_root),
    paste("Completed:", sum(status_df$completed)),
    paste("Remaining:", 0),
    paste("BAM dir:", bam_dir)
  )
  write_outputs(summary_df, log_lines)
  message("All Step 15K BAMs are already complete.")
} else {
  run_log <- character(0)

  for (target_i in initial_pending) {

    if (isTRUE(status_df$completed[target_i])) next

    target_sample <- run_manifest$sample_id[target_i]
    target_group <- run_manifest$core_group[target_i]
    target_bam <- run_manifest$bam_absolute[target_i]

    if (file.exists(target_bam)) suppressWarnings(file.remove(target_bam))

    gc()
    start_time <- Sys.time()

    fit <- tryCatch(
      {
        Rsubread::align(
          index = index_abs,
          readfile1 = run_manifest$fastq_r1_absolute[target_i],
          readfile2 = run_manifest$fastq_r2_absolute[target_i],
          input_format = "gzFASTQ",
          output_file = target_bam,
          output_format = "BAM",
          phredOffset = phred_offset_used,
          nthreads = nthreads_each,
          unique = TRUE,
          nBestLocations = 1
        )
        list(success = TRUE, message = NA_character_)
      },
      error = function(e) {
        if (file.exists(target_bam)) suppressWarnings(file.remove(target_bam))
        list(success = FALSE, message = conditionMessage(e))
      }
    )

    end_time <- Sys.time()
    elapsed_min <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)

    bam_exists_i <- file.exists(target_bam)
    bam_size_mb_i <- if (bam_exists_i) round(file.info(target_bam)$size / 1024^2, 3) else NA_real_

    status_df$run_attempts[target_i] <- as.integer(status_df$run_attempts[target_i]) + 1L
    status_df$last_attempt_time[target_i] <- format(end_time, "%Y-%m-%d %H:%M:%S")
    status_df$last_error_message[target_i] <- fit$message
    status_df$bam_exists[target_i] <- bam_exists_i
    status_df$bam_size_mb[target_i] <- bam_size_mb_i
    status_df$completed[target_i] <- isTRUE(fit$success) && bam_exists_i && !is.na(bam_size_mb_i) && bam_size_mb_i > 0
    if (isTRUE(status_df$completed[target_i])) {
      status_df$last_success_time[target_i] <- format(end_time, "%Y-%m-%d %H:%M:%S")
      run_log <- c(run_log, paste0(target_sample, " | SUCCESS | ", elapsed_min, " min | BAM ", bam_size_mb_i, " MB"))
    } else {
      run_log <- c(run_log, paste0(target_sample, " | FAIL | ", elapsed_min, " min | ", fit$message))
    }

    summary_df <- make_summary(
      current_target_sample = target_sample,
      current_target_group = target_group,
      current_run_success = status_df$completed[target_i],
      current_run_elapsed_minutes = elapsed_min,
      batch_elapsed_minutes = round(as.numeric(difftime(end_time, batch_start_time, units = "mins")), 3)
    )

    log_lines <- c(
      "Step 15K batch auto-continue run.",
      paste("Project root:", project_root),
      paste("BAM dir:", bam_dir),
      paste("Completed total:", sum(status_df$completed)),
      paste("Remaining total:", sum(!status_df$completed)),
      "----- per-sample run log -----",
      run_log
    )

    write_outputs(summary_df, log_lines)
  }

  final_elapsed <- round(as.numeric(difftime(Sys.time(), batch_start_time, units = "mins")), 3)
  summary_df <- make_summary(batch_elapsed_minutes = final_elapsed)
  log_lines <- c(
    "Step 15K batch auto-continue run finished.",
    paste("Project root:", project_root),
    paste("BAM dir:", bam_dir),
    paste("Completed total:", sum(status_df$completed)),
    paste("Remaining total:", sum(!status_df$completed)),
    paste("Batch elapsed minutes:", final_elapsed),
    "----- per-sample run log -----",
    run_log
  )
  write_outputs(summary_df, log_lines)

  n_fail <- sum(!status_df$completed)
  if (n_fail > 0) {
    message(paste0(
      "Step 15K finished batch run. ",
      sum(status_df$completed), " completed, ",
      n_fail, " still incomplete. Check step15K_pig_early_alignment_status_table.csv."
    ))
  } else {
    message("Step 15K finished batch run successfully. All samples are complete.")
  }
}
