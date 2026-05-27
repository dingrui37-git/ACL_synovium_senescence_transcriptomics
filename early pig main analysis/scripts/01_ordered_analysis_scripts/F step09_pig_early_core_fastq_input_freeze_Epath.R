# English note: E-path standardized version. This script uses E:/R/ACLsenescence2 as the locked project root and preserves the original Step logic and output file names.
# ============================================================
# Step 09: freeze pig early core sample FASTQ inputs for rebuild
# Uses existing local raw FASTQ files only; no redownload
# Includes automatic cleanup on failure
# ============================================================

setwd("E:/R/ACLsenescence2")

project_root <- getwd()
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_dir <- file.path(rebuild_root, "02_pig_early")
scripts_dir <- file.path(pig_dir, "scripts")
objects_dir <- file.path(pig_dir, "objects")
tables_dir <- file.path(pig_dir, "tables")
logs_dir <- file.path(pig_dir, "logs")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

step09_targets <- c(
  file.path(objects_dir, "step09_pig_early_core_fastq_input_freeze.RData"),
  file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv"),
  file.path(tables_dir, "step09_pig_early_core_group_summary.csv"),
  file.path(tables_dir, "step09_pig_early_fastq_pair_check_summary.csv"),
  file.path(tables_dir, "step09_pig_early_missing_or_incomplete_pairs.csv"),
  file.path(logs_dir, "step09_pig_early_core_fastq_input_freeze_log.txt"),
  file.path(scripts_dir, "step09_pig_early_core_fastq_input_freeze.R"),
  file.path(scripts_dir, "step09_check_pig_early_core_fastq_input_freeze.R")
)

cleanup_step09 <- function() {
  existing <- step09_targets[file.exists(step09_targets)]
  if (length(existing) > 0) file.remove(existing)
}

cleanup_step09()

step09_main <- function() {

  data_dir <- file.path(project_root, "data_raw", "E-MTAB-6664")
  manifest_file <- file.path(data_dir, "E-MTAB-6664_core_validation_manifest_wide.csv")

  if (!file.exists(manifest_file)) {
    stop("Cannot find core validation manifest: ", manifest_file)
  }

  manifest_raw <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  manifest_raw <- as.data.frame(manifest_raw, stringsAsFactors = FALSE)

  # ----------------------------------------------------------
  # 1. Identify sample_id column robustly
  # ----------------------------------------------------------
  colnames_raw <- colnames(manifest_raw)

  sample_col_by_name <- colnames_raw[
    grepl("^sample(_id)?$", colnames_raw, ignore.case = TRUE)
  ]

  if (length(sample_col_by_name) > 0) {
    sample_col <- sample_col_by_name[1]
  } else {
    match_counts <- vapply(
      manifest_raw,
      function(x) {
        vals <- trimws(as.character(x))
        sum(grepl("^(CON[1-6]|INJS[1-6]|INJL[1-6])$", vals, ignore.case = FALSE), na.rm = TRUE)
      },
      integer(1)
    )
    if (max(match_counts) == 0) {
      stop("Failed to identify sample ID column in manifest.")
    }
    sample_col <- names(which.max(match_counts))
  }

  sample_ids <- trimws(as.character(manifest_raw[[sample_col]]))
  keep_idx <- grepl("^(CON[1-6]|INJS[1-6]|INJL[1-6])$", sample_ids)
  core_manifest <- manifest_raw[keep_idx, , drop = FALSE]
  core_manifest$sample_id <- trimws(as.character(core_manifest[[sample_col]]))

  if (nrow(core_manifest) != 18) {
    stop("Expected 18 core pig early samples, but found: ", nrow(core_manifest))
  }

  # ----------------------------------------------------------
  # 2. Derive core group labels if not clearly present
  # ----------------------------------------------------------
  core_manifest$core_group <- ifelse(
    grepl("^CON", core_manifest$sample_id), "CON_t0",
    ifelse(
      grepl("^INJS", core_manifest$sample_id), "ACLT_untreated_t7",
      ifelse(grepl("^INJL", core_manifest$sample_id), "ACLT_untreated_t28", NA)
    )
  )

  if (any(is.na(core_manifest$core_group))) {
    stop("Failed to derive core_group for one or more samples.")
  }

  # ----------------------------------------------------------
  # 3. Match existing local FASTQ files
  # ----------------------------------------------------------
  all_fastq <- list.files(
    data_dir,
    pattern = "\\.fastq\\.gz$",
    full.names = TRUE
  )

  if (length(all_fastq) == 0) {
    stop("No FASTQ files found in: ", data_dir)
  }

  rel_path <- function(x) {
    sub(
      paste0("^", gsub("\\\\", "/", normalizePath(project_root, winslash = "/")), "/?"),
      "",
      gsub("\\\\", "/", normalizePath(x, winslash = "/"))
    )
  }

  find_fastq_for_sample <- function(sample_id, read_no) {
    bn <- basename(all_fastq)
    pat <- paste0("^", sample_id, "_R", read_no, "\\.fastq\\.gz$")
    hit <- all_fastq[grepl(pat, bn, ignore.case = FALSE)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }

  core_manifest$fastq_r1 <- vapply(core_manifest$sample_id, find_fastq_for_sample, character(1), read_no = 1)
  core_manifest$fastq_r2 <- vapply(core_manifest$sample_id, find_fastq_for_sample, character(1), read_no = 2)

  core_manifest$fastq_r1_exists <- file.exists(core_manifest$fastq_r1)
  core_manifest$fastq_r2_exists <- file.exists(core_manifest$fastq_r2)

  core_manifest$fastq_r1_size_bytes <- ifelse(
    core_manifest$fastq_r1_exists,
    file.info(core_manifest$fastq_r1)$size,
    NA
  )
  core_manifest$fastq_r2_size_bytes <- ifelse(
    core_manifest$fastq_r2_exists,
    file.info(core_manifest$fastq_r2)$size,
    NA
  )

  core_manifest$fastq_r1_size_mb <- round(core_manifest$fastq_r1_size_bytes / 1024^2, 3)
  core_manifest$fastq_r2_size_mb <- round(core_manifest$fastq_r2_size_bytes / 1024^2, 3)

  core_manifest$pair_complete <- core_manifest$fastq_r1_exists & core_manifest$fastq_r2_exists
  core_manifest$total_pair_size_gb <- round(
    (core_manifest$fastq_r1_size_bytes + core_manifest$fastq_r2_size_bytes) / 1024^3,
    3
  )

  core_manifest$fastq_r1_relative <- ifelse(core_manifest$fastq_r1_exists, rel_path(core_manifest$fastq_r1), NA)
  core_manifest$fastq_r2_relative <- ifelse(core_manifest$fastq_r2_exists, rel_path(core_manifest$fastq_r2), NA)

  core_manifest_out <- core_manifest[, c(
    "sample_id", "core_group",
    "fastq_r1_relative", "fastq_r2_relative",
    "fastq_r1_exists", "fastq_r2_exists",
    "fastq_r1_size_mb", "fastq_r2_size_mb",
    "pair_complete", "total_pair_size_gb"
  )]

  # ----------------------------------------------------------
  # 4. Summaries
  # ----------------------------------------------------------
  group_summary <- aggregate(
    sample_id ~ core_group,
    data = core_manifest_out,
    FUN = length
  )
  colnames(group_summary)[2] <- "n_samples"

  pair_check_summary <- data.frame(
    metric = c(
      "n_core_samples",
      "n_r1_found",
      "n_r2_found",
      "n_complete_pairs",
      "n_incomplete_pairs"
    ),
    value = c(
      nrow(core_manifest_out),
      sum(core_manifest_out$fastq_r1_exists),
      sum(core_manifest_out$fastq_r2_exists),
      sum(core_manifest_out$pair_complete),
      sum(!core_manifest_out$pair_complete)
    ),
    stringsAsFactors = FALSE
  )

  incomplete_df <- subset(
    core_manifest_out,
    !pair_complete,
    select = c(
      "sample_id", "core_group",
      "fastq_r1_relative", "fastq_r2_relative",
      "fastq_r1_exists", "fastq_r2_exists"
    )
  )

  # ----------------------------------------------------------
  # 5. Save
  # ----------------------------------------------------------
  write.csv(
    core_manifest_out,
    file = file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv"),
    row.names = FALSE
  )

  write.csv(
    group_summary,
    file = file.path(tables_dir, "step09_pig_early_core_group_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    pair_check_summary,
    file = file.path(tables_dir, "step09_pig_early_fastq_pair_check_summary.csv"),
    row.names = FALSE
  )

  write.csv(
    incomplete_df,
    file = file.path(tables_dir, "step09_pig_early_missing_or_incomplete_pairs.csv"),
    row.names = FALSE
  )

  save(
    core_manifest_out,
    group_summary,
    pair_check_summary,
    incomplete_df,
    file = file.path(objects_dir, "step09_pig_early_core_fastq_input_freeze.RData")
  )

  log_lines <- c(
    "Step 09 completed: froze pig early core FASTQ inputs for rebuild.",
    paste("Manifest file:", manifest_file),
    paste("Detected sample_id column:", sample_col),
    paste("Core samples:", nrow(core_manifest_out)),
    paste("Complete FASTQ pairs:", sum(core_manifest_out$pair_complete)),
    paste("Incomplete FASTQ pairs:", sum(!core_manifest_out$pair_complete))
  )

  writeLines(
    log_lines,
    con = file.path(logs_dir, "step09_pig_early_core_fastq_input_freeze_log.txt")
  )

  writeLines(
    c(
      "# Step 09 operation script was executed from console and archived automatically.",
      "setwd(\"E:/R/ACLsenescence2\")",
      "cat(\"Step 09 outputs are saved under rebuild_submission/02_pig_early/.\\n\")"
    ),
    con = file.path(scripts_dir, "step09_pig_early_core_fastq_input_freeze.R")
  )

  writeLines(
    c(
      "setwd(\"E:/R/ACLsenescence2\")",
      "pig_dir <- file.path(getwd(), \"rebuild_submission\", \"02_pig_early\")",
      "objects_dir <- file.path(pig_dir, \"objects\")",
      "tables_dir <- file.path(pig_dir, \"tables\")",
      "logs_dir <- file.path(pig_dir, \"logs\")",
      "",
      "cat(\"===== Step 09 key file check =====\\n\")",
      "key_files <- c(",
      "  file.path(objects_dir, \"step09_pig_early_core_fastq_input_freeze.RData\"),",
      "  file.path(tables_dir, \"step09_pig_early_core_sample_fastq_manifest.csv\"),",
      "  file.path(tables_dir, \"step09_pig_early_core_group_summary.csv\"),",
      "  file.path(tables_dir, \"step09_pig_early_fastq_pair_check_summary.csv\"),",
      "  file.path(tables_dir, \"step09_pig_early_missing_or_incomplete_pairs.csv\"),",
      "  file.path(logs_dir, \"step09_pig_early_core_fastq_input_freeze_log.txt\")",
      ")",
      "print(data.frame(file = key_files, exists = file.exists(key_files), size_bytes = ifelse(file.exists(key_files), file.info(key_files)$size, NA), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 09 group summary =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step09_pig_early_core_group_summary.csv\"), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 09 FASTQ pair check =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step09_pig_early_fastq_pair_check_summary.csv\"), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 09 manifest preview =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step09_pig_early_core_sample_fastq_manifest.csv\"), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 09 incomplete pairs =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step09_pig_early_missing_or_incomplete_pairs.csv\"), stringsAsFactors = FALSE))"
    ),
    con = file.path(scripts_dir, "step09_check_pig_early_core_fastq_input_freeze.R")
  )
}

tryCatch(
  {
    step09_main()
    cat("Step 09 completed successfully.\n")
    cat("Main outputs:\n")
    print(c(
      file.path(objects_dir, "step09_pig_early_core_fastq_input_freeze.RData"),
      file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv"),
      file.path(tables_dir, "step09_pig_early_fastq_pair_check_summary.csv")
    ))
  },
  error = function(e) {
    cleanup_step09()
    message("Step 09 failed, and partial output files for this step were cleaned automatically.")
    stop(e)
  }
)