# ============================================================
# Step 10J: batch FastQC on all 36 pig-early FASTQ files
# Strategy:
#   - stage FASTQ files to ASCII-only working directory
#   - run run_fastqc.bat on each file with NO options
#   - FastQC outputs html/zip into the same staged directory
# Includes automatic cleanup on prerequisite failure
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
ascii_batch_log_dir <- file.path(ascii_batch_dir, "per_fastq_logs")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

step10j_targets <- c(
  file.path(objects_dir, "step10_fastqc_run_workspace.RData"),
  file.path(tables_dir, "step10_fastqc_tool_summary.csv"),
  file.path(tables_dir, "step10_fastqc_execution_table.csv"),
  file.path(tables_dir, "step10_fastqc_run_summary.csv"),
  file.path(logs_dir, "step10_fastqc_run_log.txt"),
  file.path(scripts_dir, "step10_run_fastqc_on_pig_early_core_fastq.R"),
  file.path(scripts_dir, "step10_check_run_fastqc_on_pig_early_core_fastq.R")
)

cleanup_step10j <- function() {
  existing <- step10j_targets[file.exists(step10j_targets)]
  if (length(existing) > 0) file.remove(existing)
  if (dir.exists(ascii_batch_dir)) unlink(ascii_batch_dir, recursive = TRUE, force = TRUE)
}

cleanup_step10j()

step10j_main <- function() {

  manifest_file <- file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv")
  if (!file.exists(manifest_file)) {
    stop("Cannot find Step 09 FASTQ manifest: ", manifest_file)
  }

  manifest_df <- read.csv(manifest_file, stringsAsFactors = FALSE)
  manifest_df <- as.data.frame(manifest_df, stringsAsFactors = FALSE)

  required_cols <- c("sample_id", "core_group", "fastq_r1_relative", "fastq_r2_relative", "pair_complete")
  if (!all(required_cols %in% colnames(manifest_df))) {
    stop("Step 09 manifest is missing required columns.")
  }
  if (!all(manifest_df$pair_complete)) {
    stop("Step 09 manifest contains incomplete FASTQ pairs.")
  }

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

  manifest_long$source_fastq <- file.path(project_root, manifest_long$fastq_relative)
  manifest_long$source_exists <- file.exists(manifest_long$source_fastq)

  if (!all(manifest_long$source_exists)) {
    stop("One or more FASTQ files from Step 09 manifest do not exist.")
  }

  fastqc_exec <- "D:/R/FastQC/fastqc_v0.12.1/FastQC/run_fastqc.bat"
  java_exec   <- "D:/丁锐/Eclipse Temrin/bin/java.exe"
  fastqc_workdir <- dirname(fastqc_exec)
  java_home_dir  <- dirname(dirname(java_exec))

  if (!file.exists(fastqc_exec)) stop("FastQC launcher not found: ", fastqc_exec)
  if (!file.exists(java_exec)) stop("Java executable not found: ", java_exec)

  dir.create(ascii_batch_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ascii_batch_log_dir, recursive = TRUE, showWarnings = FALSE)

  # ----------------------------------------------------------
  # 1. Stage 36 FASTQ files into ASCII directory
  # ----------------------------------------------------------
  manifest_long$staged_fastq <- file.path(ascii_batch_dir, basename(manifest_long$source_fastq))

  copy_ok <- mapply(
    FUN = function(src, dst) file.copy(src, dst, overwrite = TRUE),
    src = manifest_long$source_fastq,
    dst = manifest_long$staged_fastq
  )

  manifest_long$copy_ok <- as.logical(copy_ok)
  manifest_long$staged_exists <- file.exists(manifest_long$staged_fastq)

  if (!all(manifest_long$copy_ok & manifest_long$staged_exists)) {
    stop("Failed to stage one or more FASTQ files into ASCII batch directory.")
  }

  # ----------------------------------------------------------
  # 2. Tool summary
  # ----------------------------------------------------------
  tool_summary_df <- data.frame(
    tool = c("FastQC", "Java"),
    executable = c(fastqc_exec, java_exec),
    version = c(NA, NA),
    stringsAsFactors = FALSE
  )

  write.csv(
    tool_summary_df,
    file = file.path(tables_dir, "step10_fastqc_tool_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  # ----------------------------------------------------------
  # 3. Run FastQC one file at a time, NO options
  # ----------------------------------------------------------
  run_fastqc_one <- function(staged_fastq, sample_id, read_label) {

    fq_base <- basename(staged_fastq)
    fq_stub <- sub("\\.fastq\\.gz$", "", fq_base, ignore.case = TRUE)

    stdout_file <- file.path(ascii_batch_log_dir, paste0(fq_stub, "_stdout.txt"))
    stderr_file <- file.path(ascii_batch_log_dir, paste0(fq_stub, "_stderr.txt"))

    cmd_string <- paste0(
      "cd /d ", shQuote(fastqc_workdir),
      " && set \"JAVA_HOME=", java_home_dir, "\"",
      " && set \"PATH=", dirname(java_exec), ";%PATH%\"",
      " && call ", shQuote(fastqc_exec), " ", shQuote(staged_fastq)
    )

    status <- tryCatch(
      system2(
        "cmd.exe",
        args = c("/c", cmd_string),
        stdout = stdout_file,
        stderr = stderr_file
      ),
      error = function(e) structure(NA_integer_, class = "cmd_error", msg = conditionMessage(e))
    )

    if (inherits(status, "cmd_error")) {
      exit_code <- NA_integer_
    } else {
      exit_code <- ifelse(is.null(status), 0L, as.integer(status))
    }

    html_out <- file.path(ascii_batch_dir, paste0(fq_stub, "_fastqc.html"))
    zip_out  <- file.path(ascii_batch_dir, paste0(fq_stub, "_fastqc.zip"))

    data.frame(
      sample_id = sample_id,
      read_label = read_label,
      staged_fastq = staged_fastq,
      stdout_file = stdout_file,
      stderr_file = stderr_file,
      exit_code = exit_code,
      html_exists = file.exists(html_out),
      zip_exists = file.exists(zip_out),
      success = !is.na(exit_code) & exit_code == 0L & file.exists(html_out) & file.exists(zip_out),
      stringsAsFactors = FALSE
    )
  }

  execution_list <- vector("list", nrow(manifest_long))

  for (i in seq_len(nrow(manifest_long))) {
    execution_list[[i]] <- cbind(
      manifest_long[i, c("sample_id", "core_group", "read_label", "fastq_relative", "source_fastq", "staged_fastq")],
      run_fastqc_one(
        staged_fastq = manifest_long$staged_fastq[i],
        sample_id = manifest_long$sample_id[i],
        read_label = manifest_long$read_label[i]
      )[, c("stdout_file", "stderr_file", "exit_code", "html_exists", "zip_exists", "success")],
      stringsAsFactors = FALSE
    )
  }

  execution_df <- do.call(rbind, execution_list)

  write.csv(
    execution_df,
    file = file.path(tables_dir, "step10_fastqc_execution_table.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  # ----------------------------------------------------------
  # 4. Run summary
  # ----------------------------------------------------------
  run_summary_df <- data.frame(
    metric = c(
      "n_fastq_input",
      "n_fastq_staged",
      "n_fastqc_success",
      "n_fastqc_failed",
      "n_fastqc_html_found",
      "n_fastqc_zip_found"
    ),
    value = c(
      nrow(manifest_long),
      sum(manifest_long$staged_exists, na.rm = TRUE),
      sum(execution_df$success, na.rm = TRUE),
      sum(!execution_df$success, na.rm = TRUE),
      sum(execution_df$html_exists, na.rm = TRUE),
      sum(execution_df$zip_exists, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    run_summary_df,
    file = file.path(tables_dir, "step10_fastqc_run_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  save(
    manifest_long,
    tool_summary_df,
    execution_df,
    run_summary_df,
    ascii_batch_dir,
    ascii_batch_log_dir,
    file = file.path(objects_dir, "step10_fastqc_run_workspace.RData")
  )

  log_lines <- c(
    "Step 10J completed: batch FastQC on 36 pig-early FASTQ files using ASCII staging directory and launcher with NO options.",
    paste("ASCII batch directory:", ascii_batch_dir),
    paste("Per-fastq log directory:", ascii_batch_log_dir),
    paste("Input FASTQ files:", nrow(manifest_long)),
    paste("Staged FASTQ files:", sum(manifest_long$staged_exists, na.rm = TRUE)),
    paste("FastQC successes:", sum(execution_df$success, na.rm = TRUE)),
    paste("FastQC failures:", sum(!execution_df$success, na.rm = TRUE))
  )

  writeLines(
    log_lines,
    con = file.path(logs_dir, "step10_fastqc_run_log.txt"),
    useBytes = TRUE
  )

  writeLines(
    c(
      "# Step 10J operation script was executed from console and archived automatically.",
      "setwd(\"D:/R/ACL ∩ senescence2\")",
      "cat(\"Step 10 outputs are saved under rebuild_submission/02_pig_early/ plus ASCII FastQC working directory.\\n\")"
    ),
    con = file.path(scripts_dir, "step10_run_fastqc_on_pig_early_core_fastq.R"),
    useBytes = TRUE
  )

  writeLines(
    c(
      "setwd(\"D:/R/ACL ∩ senescence2\")",
      "pig_dir <- file.path(getwd(), \"rebuild_submission\", \"02_pig_early\")",
      "objects_dir <- file.path(pig_dir, \"objects\")",
      "tables_dir <- file.path(pig_dir, \"tables\")",
      "logs_dir <- file.path(pig_dir, \"logs\")",
      "ascii_batch_dir <- \"D:/R/FastQC_ascii_batch_36\"",
      "",
      "cat(\"===== Step 10 key file check =====\\n\")",
      "key_files <- c(",
      "  file.path(objects_dir, \"step10_fastqc_run_workspace.RData\"),",
      "  file.path(tables_dir, \"step10_fastqc_tool_summary.csv\"),",
      "  file.path(tables_dir, \"step10_fastqc_execution_table.csv\"),",
      "  file.path(tables_dir, \"step10_fastqc_run_summary.csv\"),",
      "  file.path(logs_dir, \"step10_fastqc_run_log.txt\")",
      ")",
      "print(data.frame(file = key_files, exists = file.exists(key_files), size_bytes = ifelse(file.exists(key_files), file.info(key_files)$size, NA), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 10 tool summary =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step10_fastqc_tool_summary.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 10 run summary =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step10_fastqc_run_summary.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 10 execution preview =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step10_fastqc_execution_table.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\")[, c(\"sample_id\",\"read_label\",\"success\",\"html_exists\",\"zip_exists\",\"exit_code\")])",
      "",
      "cat(\"\\n===== Step 10 FastQC output count =====\\n\")",
      "html_n <- length(list.files(ascii_batch_dir, pattern = \"_fastqc\\\\.html$\", full.names = TRUE))",
      "zip_n  <- length(list.files(ascii_batch_dir, pattern = \"_fastqc\\\\.zip$\", full.names = TRUE))",
      "print(data.frame(file_type = c(\"html\",\"zip\"), n = c(html_n, zip_n), stringsAsFactors = FALSE))"
    ),
    con = file.path(scripts_dir, "step10_check_run_fastqc_on_pig_early_core_fastq.R"),
    useBytes = TRUE
  )
}

tryCatch(
  {
    step10j_main()
    cat("Step 10J completed successfully.\n")
    cat("Main outputs:\n")
    print(c(
      file.path(objects_dir, "step10_fastqc_run_workspace.RData"),
      file.path(tables_dir, "step10_fastqc_tool_summary.csv"),
      file.path(tables_dir, "step10_fastqc_execution_table.csv"),
      file.path(tables_dir, "step10_fastqc_run_summary.csv"),
      ascii_batch_dir
    ))
  },
  error = function(e) {
    cleanup_step10j()
    message("Step 10J failed, and partial output files for this step were cleaned automatically.")
    stop(e)
  }
)