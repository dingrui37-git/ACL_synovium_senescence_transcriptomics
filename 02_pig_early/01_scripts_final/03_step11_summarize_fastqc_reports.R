# ============================================================
# Step 11: summarize 36 FastQC reports into structured QC tables
# Uses existing FastQC zip reports from Step 10J
# Includes automatic cleanup on failure
# ============================================================

setwd("D:/R/ACL ∩ senescence2")

project_root <- getwd()
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_dir <- file.path(rebuild_root, "02_pig_early")
scripts_dir <- file.path(pig_dir, "scripts")
objects_dir <- file.path(pig_dir, "objects")
tables_dir <- file.path(pig_dir, "tables")
logs_dir <- file.path(pig_dir, "logs")

ascii_batch_dir <- "D:/R/FastQC_ascii_batch_36"

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

step11_targets <- c(
  file.path(objects_dir, "step11_fastqc_summary_workspace.RData"),
  file.path(tables_dir, "step11_fastqc_basic_stats_by_fastq.csv"),
  file.path(tables_dir, "step11_fastqc_module_status_long.csv"),
  file.path(tables_dir, "step11_fastqc_module_status_wide.csv"),
  file.path(tables_dir, "step11_fastqc_module_summary_counts.csv"),
  file.path(tables_dir, "step11_fastqc_overall_qc_summary.csv"),
  file.path(logs_dir, "step11_fastqc_summary_log.txt"),
  file.path(scripts_dir, "step11_summarize_fastqc_reports.R"),
  file.path(scripts_dir, "step11_check_summarize_fastqc_reports.R")
)

cleanup_step11 <- function() {
  existing <- step11_targets[file.exists(step11_targets)]
  if (length(existing) > 0) file.remove(existing)
}

cleanup_step11()

step11_main <- function() {

  manifest_file <- file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv")
  if (!file.exists(manifest_file)) {
    stop("Cannot find Step 09 FASTQ manifest: ", manifest_file)
  }

  if (!dir.exists(ascii_batch_dir)) {
    stop("Cannot find ASCII FastQC batch directory: ", ascii_batch_dir)
  }

  zip_files <- list.files(ascii_batch_dir, pattern = "_fastqc\\.zip$", full.names = TRUE)
  if (length(zip_files) != 36) {
    stop("Expected 36 FastQC zip files, but found: ", length(zip_files))
  }

  manifest_df <- read.csv(manifest_file, stringsAsFactors = FALSE)
  manifest_long <- rbind(
    data.frame(
      sample_id = manifest_df$sample_id,
      core_group = manifest_df$core_group,
      read_label = "R1",
      fastq_relative = manifest_df$fastq_r1_relative,
      stringsAsFactors = FALSE
    ),
    data.frame(
      sample_id = manifest_df$sample_id,
      core_group = manifest_df$core_group,
      read_label = "R2",
      fastq_relative = manifest_df$fastq_r2_relative,
      stringsAsFactors = FALSE
    )
  )
  manifest_long$fastq_basename <- basename(manifest_long$fastq_relative)
  manifest_long$fastq_stub <- sub("\\.fastq\\.gz$", "", manifest_long$fastq_basename, ignore.case = TRUE)

  parse_one_fastqc_zip <- function(zip_path) {

    zip_listing <- utils::unzip(zip_path, list = TRUE)

    summary_member <- zip_listing$Name[grepl("/summary\\.txt$", zip_listing$Name)]
    data_member    <- zip_listing$Name[grepl("/fastqc_data\\.txt$", zip_listing$Name)]

    if (length(summary_member) != 1 || length(data_member) != 1) {
      stop("Could not uniquely locate summary.txt or fastqc_data.txt in: ", zip_path)
    }

    exdir <- tempfile("fastqc_parse_")
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

    utils::unzip(zip_path, files = c(summary_member, data_member), exdir = exdir)

    summary_path <- file.path(exdir, summary_member)
    data_path    <- file.path(exdir, data_member)

    summary_df <- read.delim(
      summary_path,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "",
      fill = TRUE
    )
    colnames(summary_df) <- c("status", "module_name", "reported_filename")

    raw_lines <- readLines(data_path, warn = FALSE, encoding = "UTF-8")

    start_idx <- grep("^>>Basic Statistics", raw_lines)
    if (length(start_idx) == 0) {
      stop("Basic Statistics section not found in: ", zip_path)
    }
    start_idx <- start_idx[1]

    end_rel <- grep("^>>END_MODULE", raw_lines[(start_idx + 1):length(raw_lines)])
    if (length(end_rel) == 0) {
      stop("Basic Statistics section end not found in: ", zip_path)
    }
    end_idx <- start_idx + end_rel[1]

    basic_lines <- raw_lines[(start_idx + 1):(end_idx - 1)]
    basic_lines <- basic_lines[!grepl("^#", basic_lines)]
    basic_lines <- basic_lines[nzchar(trimws(basic_lines))]

    basic_df <- read.delim(
      text = paste(basic_lines, collapse = "\n"),
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE,
      quote = "",
      fill = TRUE
    )
    colnames(basic_df) <- c("measure", "value")

    stats <- setNames(as.list(basic_df$value), basic_df$measure)

    fastq_zip_name <- basename(zip_path)
    fastq_stub <- sub("_fastqc\\.zip$", "", fastq_zip_name, ignore.case = TRUE)

    basic_row <- data.frame(
      fastq_zip = fastq_zip_name,
      fastq_stub = fastq_stub,
      total_sequences = if ("Total Sequences" %in% names(stats)) stats[["Total Sequences"]] else NA,
      poor_quality_sequences = if ("Sequences flagged as poor quality" %in% names(stats)) stats[["Sequences flagged as poor quality"]] else NA,
      sequence_length = if ("Sequence length" %in% names(stats)) stats[["Sequence length"]] else NA,
      gc_percent = if ("%GC" %in% names(stats)) stats[["%GC"]] else NA,
      stringsAsFactors = FALSE
    )

    module_long <- data.frame(
      fastq_zip = fastq_zip_name,
      fastq_stub = fastq_stub,
      module_name = summary_df$module_name,
      status = summary_df$status,
      stringsAsFactors = FALSE
    )

    list(
      basic_row = basic_row,
      module_long = module_long
    )
  }

  parsed_list <- lapply(zip_files, parse_one_fastqc_zip)

  basic_stats_df <- do.call(rbind, lapply(parsed_list, `[[`, "basic_row"))
  module_long_df <- do.call(rbind, lapply(parsed_list, `[[`, "module_long"))

  basic_stats_df <- merge(
    manifest_long[, c("sample_id", "core_group", "read_label", "fastq_stub")],
    basic_stats_df,
    by = "fastq_stub",
    all.y = TRUE,
    sort = FALSE
  )

  module_long_df <- merge(
    manifest_long[, c("sample_id", "core_group", "read_label", "fastq_stub")],
    module_long_df,
    by = "fastq_stub",
    all.y = TRUE,
    sort = FALSE
  )

  basic_stats_df$total_sequences <- suppressWarnings(as.numeric(basic_stats_df$total_sequences))
  basic_stats_df$poor_quality_sequences <- suppressWarnings(as.numeric(basic_stats_df$poor_quality_sequences))
  basic_stats_df$gc_percent <- suppressWarnings(as.numeric(basic_stats_df$gc_percent))

  module_wide_df <- reshape(
    module_long_df[, c("sample_id", "core_group", "read_label", "fastq_stub", "module_name", "status")],
    idvar = c("sample_id", "core_group", "read_label", "fastq_stub"),
    timevar = "module_name",
    direction = "wide"
  )

  colnames(module_wide_df) <- sub("^status\\.", "", colnames(module_wide_df))

  module_counts_df <- as.data.frame(table(module_long_df$module_name, module_long_df$status), stringsAsFactors = FALSE)
  colnames(module_counts_df) <- c("module_name", "status", "n_fastq")

  count_status <- function(module_name, wanted_statuses = c("WARN", "FAIL")) {
    sum(module_long_df$module_name == module_name & module_long_df$status %in% wanted_statuses, na.rm = TRUE)
  }

  overall_qc_summary_df <- data.frame(
    metric = c(
      "n_fastq_reports",
      "n_fastq_with_any_fail",
      "n_fastq_with_any_warn_or_fail",
      "n_adapter_content_warn_or_fail",
      "n_per_base_sequence_quality_warn_or_fail",
      "n_per_sequence_quality_scores_warn_or_fail",
      "n_overrepresented_sequences_warn_or_fail",
      "n_sequence_duplication_levels_warn_or_fail",
      "n_per_base_n_content_warn_or_fail",
      "median_total_sequences",
      "median_gc_percent"
    ),
    value = c(
      length(unique(basic_stats_df$fastq_stub)),
      length(unique(module_long_df$fastq_stub[module_long_df$status == "FAIL"])),
      length(unique(module_long_df$fastq_stub[module_long_df$status %in% c("WARN", "FAIL")])),
      count_status("Adapter Content"),
      count_status("Per base sequence quality"),
      count_status("Per sequence quality scores"),
      count_status("Overrepresented sequences"),
      count_status("Sequence Duplication Levels"),
      count_status("Per base N content"),
      median(basic_stats_df$total_sequences, na.rm = TRUE),
      median(basic_stats_df$gc_percent, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    basic_stats_df,
    file = file.path(tables_dir, "step11_fastqc_basic_stats_by_fastq.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    module_long_df,
    file = file.path(tables_dir, "step11_fastqc_module_status_long.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    module_wide_df,
    file = file.path(tables_dir, "step11_fastqc_module_status_wide.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    module_counts_df,
    file = file.path(tables_dir, "step11_fastqc_module_summary_counts.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    overall_qc_summary_df,
    file = file.path(tables_dir, "step11_fastqc_overall_qc_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  save(
    basic_stats_df,
    module_long_df,
    module_wide_df,
    module_counts_df,
    overall_qc_summary_df,
    file = file.path(objects_dir, "step11_fastqc_summary_workspace.RData")
  )

  log_lines <- c(
    "Step 11 completed: summarized 36 FastQC reports into structured QC tables.",
    paste("FastQC zip files parsed:", length(zip_files)),
    paste("Basic stats rows:", nrow(basic_stats_df)),
    paste("Module status rows:", nrow(module_long_df)),
    paste("FASTQ with any FAIL:", overall_qc_summary_df$value[overall_qc_summary_df$metric == "n_fastq_with_any_fail"]),
    paste("FASTQ with any WARN/FAIL:", overall_qc_summary_df$value[overall_qc_summary_df$metric == "n_fastq_with_any_warn_or_fail"])
  )

  writeLines(
    log_lines,
    con = file.path(logs_dir, "step11_fastqc_summary_log.txt"),
    useBytes = TRUE
  )

  writeLines(
    c(
      "# Step 11 operation script was executed from console and archived automatically.",
      "setwd(\"D:/R/ACL ∩ senescence2\")",
      "cat(\"Step 11 outputs are saved under rebuild_submission/02_pig_early/.\\n\")"
    ),
    con = file.path(scripts_dir, "step11_summarize_fastqc_reports.R"),
    useBytes = TRUE
  )

  writeLines(
    c(
      "setwd(\"D:/R/ACL ∩ senescence2\")",
      "pig_dir <- file.path(getwd(), \"rebuild_submission\", \"02_pig_early\")",
      "objects_dir <- file.path(pig_dir, \"objects\")",
      "tables_dir <- file.path(pig_dir, \"tables\")",
      "logs_dir <- file.path(pig_dir, \"logs\")",
      "",
      "cat(\"===== Step 11 key file check =====\\n\")",
      "key_files <- c(",
      "  file.path(objects_dir, \"step11_fastqc_summary_workspace.RData\"),",
      "  file.path(tables_dir, \"step11_fastqc_basic_stats_by_fastq.csv\"),",
      "  file.path(tables_dir, \"step11_fastqc_module_status_long.csv\"),",
      "  file.path(tables_dir, \"step11_fastqc_module_status_wide.csv\"),",
      "  file.path(tables_dir, \"step11_fastqc_module_summary_counts.csv\"),",
      "  file.path(tables_dir, \"step11_fastqc_overall_qc_summary.csv\"),",
      "  file.path(logs_dir, \"step11_fastqc_summary_log.txt\")",
      ")",
      "print(data.frame(file = key_files, exists = file.exists(key_files), size_bytes = ifelse(file.exists(key_files), file.info(key_files)$size, NA), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 11 overall QC summary =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step11_fastqc_overall_qc_summary.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 11 module summary counts (first 30 rows) =====\\n\")",
      "print(head(read.csv(file.path(tables_dir, \"step11_fastqc_module_summary_counts.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"), 30))",
      "",
      "cat(\"\\n===== Step 11 basic stats preview =====\\n\")",
      "print(head(read.csv(file.path(tables_dir, \"step11_fastqc_basic_stats_by_fastq.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"), 20))"
    ),
    con = file.path(scripts_dir, "step11_check_summarize_fastqc_reports.R"),
    useBytes = TRUE
  )
}

tryCatch(
  {
    step11_main()
    cat("Step 11 completed successfully.\n")
    cat("Main outputs:\n")
    print(c(
      file.path(objects_dir, "step11_fastqc_summary_workspace.RData"),
      file.path(tables_dir, "step11_fastqc_basic_stats_by_fastq.csv"),
      file.path(tables_dir, "step11_fastqc_module_status_long.csv"),
      file.path(tables_dir, "step11_fastqc_module_summary_counts.csv"),
      file.path(tables_dir, "step11_fastqc_overall_qc_summary.csv")
    ))
  },
  error = function(e) {
    cleanup_step11()
    message("Step 11 failed, and partial output files for this step were cleaned automatically.")
    stop(e)
  }
)