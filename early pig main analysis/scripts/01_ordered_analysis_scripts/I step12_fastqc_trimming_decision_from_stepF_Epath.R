# English note: This E-path standardized Step12 script formalizes the pig early trimming decision from the StepE/StepF FastQC-MultiQC branch. It reads FastQC PASS/WARN/FAIL summaries generated from rebuild_submission/02_pig_early/qc/fastqc_raw, checks the key trimming-related FastQC modules, writes manuscript-ready trimming-decision outputs using the original Step12 output names, and preserves reproducibility logs under E:/R/ACLsenescence2/rebuild_submission/02_pig_early.

# ============================================================
# Step 12: formalize FastQC-based trimming decision from StepE/StepF
# Input branch:
#   StepE: run FastQC + MultiQC, store outputs in qc/fastqc_raw and qc/multiqc_raw
#   StepF: extract FastQC PASS/WARN/FAIL summaries from qc/fastqc_raw
# Purpose:
#   Use StepF module-level PASS/WARN/FAIL summaries to decide whether an
#   additional trimming branch is required for the 18 pig early core samples.
# ============================================================

options(stringsAsFactors = FALSE)

setwd("E:/R/ACLsenescence2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_dir <- file.path(rebuild_root, "02_pig_early")
scripts_dir <- file.path(pig_dir, "scripts")
objects_dir <- file.path(pig_dir, "objects")
tables_dir <- file.path(pig_dir, "tables")
logs_dir <- file.path(pig_dir, "logs")
qc_dir <- file.path(pig_dir, "qc")
fastqc_raw_dir <- file.path(qc_dir, "fastqc_raw")
multiqc_raw_dir <- file.path(qc_dir, "multiqc_raw")

for (d in c(scripts_dir, objects_dir, tables_dir, logs_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

step12_targets <- c(
  file.path(objects_dir, "step12_fastqc_trimming_decision_workspace.RData"),
  file.path(tables_dir, "step12_fastqc_trimming_decision_summary.csv"),
  file.path(tables_dir, "step12_fastqc_key_module_evidence.csv"),
  file.path(tables_dir, "step12_fastqc_input_file_audit_from_stepF.csv"),
  file.path(logs_dir, "step12_fastqc_methods_ready_note.txt"),
  file.path(logs_dir, "step12_fastqc_response_ready_note.txt"),
  file.path(logs_dir, "step12_fastqc_trimming_decision_log.txt"),
  file.path(scripts_dir, "step12_fastqc_trimming_decision.R"),
  file.path(scripts_dir, "step12_check_fastqc_trimming_decision.R")
)

cleanup_step12 <- function() {
  existing <- step12_targets[file.exists(step12_targets)]
  if (length(existing) > 0) file.remove(existing)
}
cleanup_step12()

safe_read_csv <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

as_bool <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}

get_metric <- function(df, metric_name, default = NA_character_) {
  if (is.null(df) || !all(c("metric", "value") %in% colnames(df))) return(default)
  hit <- df$value[df$metric == metric_name]
  if (length(hit) == 0) return(default)
  hit[1]
}

find_col <- function(df, candidates, label) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop(
      "Could not find ", label, ". Checked: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }
  hit[1]
}

get_count <- function(module_df, module_name, status_name) {
  module_col <- find_col(module_df, c("module", "module_name", "Module"), "FastQC module column")
  status_col <- find_col(module_df, c(status_name, tolower(status_name), tools::toTitleCase(tolower(status_name))), paste0(status_name, " count column"))
  hit <- module_df[[status_col]][module_df[[module_col]] == module_name]
  if (length(hit) == 0) return(0)
  suppressWarnings(as.numeric(hit[1]))
}

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  NA_character_
}

archive_current_script <- function(path) {
  src <- get_script_path()
  if (!is.na(src) && file.exists(src)) {
    writeLines(readLines(src, warn = FALSE, encoding = "UTF-8"), path, useBytes = TRUE)
  } else {
    writeLines(c(
      "# Step12 FastQC trimming decision from StepF branch.",
      "# Running script path was not detected because the script may have been sourced interactively.",
      "# Use the external saved script file step12_fastqc_trimming_decision_from_stepF_Epath.R as the reproducible source."
    ), path, useBytes = TRUE)
  }
}

step12_main <- function() {
  module_summary_file <- file.path(tables_dir, "stepF_fastqc_module_status_summary.csv")
  file_summary_file <- file.path(tables_dir, "stepF_fastqc_file_status_summary.csv")
  stepF_summary_file <- file.path(tables_dir, "stepF_fastqc_multiqc_pass_warn_fail_summary.csv")
  stepE_summary_file <- file.path(tables_dir, "stepE_pig_early_fastqc_multiqc_summary.csv")
  stepF_long_file <- file.path(tables_dir, "stepF_fastqc_module_status_long.csv")
  stepF_wide_file <- file.path(tables_dir, "stepF_fastqc_module_status_wide.csv")

  input_files <- c(
    module_summary_file,
    file_summary_file,
    stepF_summary_file,
    stepE_summary_file,
    stepF_long_file,
    stepF_wide_file
  )
  input_audit <- data.frame(
    file_role = c(
      "StepF module PASS/WARN/FAIL summary",
      "StepF file-level PASS/WARN/FAIL summary",
      "StepF run summary",
      "StepE FastQC/MultiQC summary",
      "StepF module status long table",
      "StepF module status wide table"
    ),
    file = input_files,
    exists = file.exists(input_files),
    size_bytes = ifelse(file.exists(input_files), file.info(input_files)$size, NA_real_),
    stringsAsFactors = FALSE
  )
  write.csv(input_audit, file.path(tables_dir, "step12_fastqc_input_file_audit_from_stepF.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  required <- c(module_summary_file, file_summary_file, stepF_summary_file)
  missing_required <- required[!file.exists(required)]
  if (length(missing_required) > 0) {
    stop("Required StepF outputs are missing:\n", paste(missing_required, collapse = "\n"), call. = FALSE)
  }

  module_df <- safe_read_csv(module_summary_file)
  file_df <- safe_read_csv(file_summary_file)
  stepF_summary <- safe_read_csv(stepF_summary_file)
  stepE_summary <- if (file.exists(stepE_summary_file)) safe_read_csv(stepE_summary_file) else NULL

  module_col <- find_col(module_df, c("module", "module_name", "Module"), "FastQC module column")
  required_modules <- c(
    "Adapter Content",
    "Per base sequence quality",
    "Per sequence quality scores"
  )
  missing_modules <- setdiff(required_modules, module_df[[module_col]])
  if (length(missing_modules) > 0) {
    stop("StepF module summary is missing key trimming-related FastQC modules: ", paste(missing_modules, collapse = ", "), call. = FALSE)
  }

  n_fastq_from_stepF <- suppressWarnings(as.numeric(get_metric(stepF_summary, "n_fastqc_zip")))
  if (is.na(n_fastq_from_stepF)) {
    n_fastq_from_stepF <- suppressWarnings(as.numeric(max(module_df$n_total, na.rm = TRUE)))
  }
  if (is.na(n_fastq_from_stepF) || n_fastq_from_stepF <= 0) {
    stop("Could not determine the number of FastQC reports from StepF outputs.", call. = FALSE)
  }
  if (n_fastq_from_stepF != 36) {
    stop("Expected 36 FastQC zip reports for 18 paired-end core samples, but StepF reports: ", n_fastq_from_stepF, call. = FALSE)
  }

  multiqc_report_exists <- as_bool(get_metric(stepF_summary, "multiqc_report_exists", FALSE))
  multiqc_data_exists <- as_bool(get_metric(stepF_summary, "multiqc_data_dir_exists", FALSE))
  if (!multiqc_report_exists && !is.null(stepE_summary)) {
    multiqc_report_exists <- as_bool(get_metric(stepE_summary, "multiqc_report_exists", FALSE))
  }

  key_modules <- data.frame(
    module_name = c(
      "Adapter Content",
      "Per base sequence quality",
      "Per sequence quality scores",
      "Per base N content",
      "Overrepresented sequences",
      "Sequence Duplication Levels",
      "Per base sequence content",
      "Per tile sequence quality",
      "Per sequence GC content"
    ),
    stringsAsFactors = FALSE
  )
  key_modules$pass_n <- vapply(key_modules$module_name, function(m) get_count(module_df, m, "PASS"), numeric(1))
  key_modules$warn_n <- vapply(key_modules$module_name, function(m) get_count(module_df, m, "WARN"), numeric(1))
  key_modules$fail_n <- vapply(key_modules$module_name, function(m) get_count(module_df, m, "FAIL"), numeric(1))
  key_modules$warn_or_fail_n <- key_modules$warn_n + key_modules$fail_n
  key_modules$warn_or_fail_pct <- round(100 * key_modules$warn_or_fail_n / n_fastq_from_stepF, 1)

  adapter_problem <- key_modules$warn_or_fail_n[key_modules$module_name == "Adapter Content"] > 0
  base_quality_problem <- key_modules$warn_or_fail_n[key_modules$module_name == "Per base sequence quality"] > 0
  seq_quality_problem <- key_modules$warn_or_fail_n[key_modules$module_name == "Per sequence quality scores"] > 0

  if (!adapter_problem && !base_quality_problem && !seq_quality_problem) {
    decision <- "no_trimming"
    decision_reason <- paste(
      "No substantial adapter or low-quality-tail evidence was detected by FastQC.",
      "Adapter Content, Per base sequence quality, and Per sequence quality scores were PASS for all 36 FASTQ files."
    )
  } else {
    decision <- "consider_trimming"
    decision_reason <- paste(
      "FastQC detected adapter and/or read-quality issues in one or more key trimming-related modules.",
      "Additional trimming should be considered before final alignment, or a trimming-sensitivity branch should be documented."
    )
  }

  decision_summary_df <- data.frame(
    metric = c(
      "input_branch",
      "project_root",
      "fastqc_raw_dir",
      "multiqc_raw_dir",
      "n_fastq_reports",
      "multiqc_report_exists",
      "multiqc_data_dir_exists",
      "adapter_content_warn_or_fail_n",
      "per_base_quality_warn_or_fail_n",
      "per_sequence_quality_warn_or_fail_n",
      "overrepresented_sequences_warn_or_fail_n",
      "sequence_duplication_warn_or_fail_n",
      "per_base_sequence_content_warn_or_fail_n",
      "recommended_branch",
      "decision_reason"
    ),
    value = c(
      "StepE/F FastQC-MultiQC branch using qc/fastqc_raw and qc/multiqc_raw",
      project_root,
      fastqc_raw_dir,
      multiqc_raw_dir,
      n_fastq_from_stepF,
      multiqc_report_exists,
      multiqc_data_exists,
      key_modules$warn_or_fail_n[key_modules$module_name == "Adapter Content"],
      key_modules$warn_or_fail_n[key_modules$module_name == "Per base sequence quality"],
      key_modules$warn_or_fail_n[key_modules$module_name == "Per sequence quality scores"],
      key_modules$warn_or_fail_n[key_modules$module_name == "Overrepresented sequences"],
      key_modules$warn_or_fail_n[key_modules$module_name == "Sequence Duplication Levels"],
      key_modules$warn_or_fail_n[key_modules$module_name == "Per base sequence content"],
      decision,
      decision_reason
    ),
    stringsAsFactors = FALSE
  )

  if (decision == "no_trimming") {
    methods_note <- c(
      "Pig early raw-read quality control note (methods-ready draft):",
      "",
      paste0(
        "Raw FASTQ quality was assessed using FastQC on all 36 paired-end files from the 18 retained pig early core samples. ",
        "FastQC outputs were collected under qc/fastqc_raw and summarized with MultiQC under qc/multiqc_raw. "
      ),
      paste0(
        "The key trimming-related modules were uniformly acceptable: Adapter Content, Per base sequence quality, ",
        "and Per sequence quality scores were PASS in all 36 FASTQ files. "
      ),
      paste0(
        "Warnings/failures in other modules, if present, were interpreted as library-composition or sequencing-context signals rather ",
        "than direct evidence requiring routine adapter/quality trimming. Accordingly, the main pig early analysis proceeded without an additional trimming step."
      )
    )
    response_note <- c(
      "Pig early raw-read QC interpretation note (response-ready draft):",
      "",
      paste0(
        "We performed FastQC for all 36 FASTQ files corresponding to the 18 retained pig early core samples and summarized the results using the StepE/F FastQC-MultiQC branch. "
      ),
      paste0(
        "The key trimming-related FastQC modules showed uniformly acceptable results: Adapter Content = PASS in 36/36 files, ",
        "Per base sequence quality = PASS in 36/36 files, and Per sequence quality scores = PASS in 36/36 files. "
      ),
      paste0(
        "Therefore, these QC results did not support mandatory adapter/quality trimming before the primary alignment/counting workflow."
      )
    )
  } else {
    methods_note <- c(
      "Pig early raw-read quality control note (methods-ready draft):",
      "",
      paste0(
        "Raw FASTQ quality was assessed using FastQC on all 36 paired-end files from the 18 retained pig early core samples. ",
        "The StepE/F summaries identified WARN/FAIL calls in one or more trimming-related modules. "
      ),
      paste0(
        "Because adapter or read-quality issues were detected in key trimming-related modules, an additional trimming or trimming-sensitivity branch should be considered before finalizing the primary alignment/counting workflow."
      )
    )
    response_note <- c(
      "Pig early raw-read QC interpretation note (response-ready draft):",
      "",
      paste0(
        "FastQC summaries from the StepE/F branch showed WARN/FAIL calls in one or more key trimming-related modules. ",
        "This supports considering an additional trimming or trimming-sensitivity analysis before finalizing the main raw-read processing branch."
      )
    )
  }

  write.csv(decision_summary_df, file.path(tables_dir, "step12_fastqc_trimming_decision_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(key_modules, file.path(tables_dir, "step12_fastqc_key_module_evidence.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(methods_note, con = file.path(logs_dir, "step12_fastqc_methods_ready_note.txt"), useBytes = TRUE)
  writeLines(response_note, con = file.path(logs_dir, "step12_fastqc_response_ready_note.txt"), useBytes = TRUE)

  save(
    module_df,
    file_df,
    stepF_summary,
    stepE_summary,
    key_modules,
    decision_summary_df,
    methods_note,
    response_note,
    input_audit,
    file = file.path(objects_dir, "step12_fastqc_trimming_decision_workspace.RData")
  )

  log_lines <- c(
    "Step 12 completed: formalized FastQC-based trimming decision from StepE/F outputs.",
    paste("Project root:", project_root),
    paste("Recommended branch:", decision),
    paste("FastQC reports:", n_fastq_from_stepF),
    paste("MultiQC report exists:", multiqc_report_exists),
    paste("Adapter warn/fail:", key_modules$warn_or_fail_n[key_modules$module_name == "Adapter Content"]),
    paste("Per base quality warn/fail:", key_modules$warn_or_fail_n[key_modules$module_name == "Per base sequence quality"]),
    paste("Per sequence quality warn/fail:", key_modules$warn_or_fail_n[key_modules$module_name == "Per sequence quality scores"]),
    paste("Decision reason:", decision_reason)
  )
  writeLines(log_lines, con = file.path(logs_dir, "step12_fastqc_trimming_decision_log.txt"), useBytes = TRUE)

  archive_current_script(file.path(scripts_dir, "step12_fastqc_trimming_decision.R"))

  writeLines(
    c(
      "setwd(\"E:/R/ACLsenescence2\")",
      "pig_dir <- file.path(getwd(), \"rebuild_submission\", \"02_pig_early\")",
      "objects_dir <- file.path(pig_dir, \"objects\")",
      "tables_dir <- file.path(pig_dir, \"tables\")",
      "logs_dir <- file.path(pig_dir, \"logs\")",
      "",
      "cat(\"===== Step 12 key file check =====\\n\")",
      "key_files <- c(",
      "  file.path(objects_dir, \"step12_fastqc_trimming_decision_workspace.RData\"),",
      "  file.path(tables_dir, \"step12_fastqc_trimming_decision_summary.csv\"),",
      "  file.path(tables_dir, \"step12_fastqc_key_module_evidence.csv\"),",
      "  file.path(tables_dir, \"step12_fastqc_input_file_audit_from_stepF.csv\"),",
      "  file.path(logs_dir, \"step12_fastqc_methods_ready_note.txt\"),",
      "  file.path(logs_dir, \"step12_fastqc_response_ready_note.txt\"),",
      "  file.path(logs_dir, \"step12_fastqc_trimming_decision_log.txt\")",
      ")",
      "print(data.frame(file = key_files, exists = file.exists(key_files), size_bytes = ifelse(file.exists(key_files), file.info(key_files)$size, NA), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 12 trimming decision summary =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step12_fastqc_trimming_decision_summary.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 12 key module evidence =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step12_fastqc_key_module_evidence.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 12 methods-ready note =====\\n\")",
      "cat(paste(readLines(file.path(logs_dir, \"step12_fastqc_methods_ready_note.txt\"), warn = FALSE), collapse = \"\\n\"))",
      "cat(\"\\n\")"
    ),
    con = file.path(scripts_dir, "step12_check_fastqc_trimming_decision.R"),
    useBytes = TRUE
  )

  invisible(list(decision_summary = decision_summary_df, key_modules = key_modules))
}

tryCatch(
  {
    out <- step12_main()
    cat("Step 12 completed successfully.\n")
    cat("Main outputs:\n")
    print(c(
      file.path(objects_dir, "step12_fastqc_trimming_decision_workspace.RData"),
      file.path(tables_dir, "step12_fastqc_trimming_decision_summary.csv"),
      file.path(tables_dir, "step12_fastqc_key_module_evidence.csv"),
      file.path(logs_dir, "step12_fastqc_methods_ready_note.txt"),
      file.path(logs_dir, "step12_fastqc_response_ready_note.txt")
    ))
    cat("\nDecision summary:\n")
    print(out$decision_summary, row.names = FALSE)
  },
  error = function(e) {
    cleanup_step12()
    message("Step 12 failed, and partial output files for this step were cleaned automatically.")
    stop(e)
  }
)
