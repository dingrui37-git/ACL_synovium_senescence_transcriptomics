# ============================================================
# Step 14: canonicalize pig-early reference choice and
# decide whether existing Subread index can be reused
# Fixed version: avoid regex escaping for index bundle detection
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

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

step14_targets <- c(
  file.path(objects_dir, "step14_pig_early_canonical_reference_workspace.RData"),
  file.path(tables_dir, "step14_pig_early_index_bundle_files.csv"),
  file.path(tables_dir, "step14_pig_early_canonical_reference_manifest.csv"),
  file.path(tables_dir, "step14_pig_early_index_reuse_decision_summary.csv"),
  file.path(logs_dir, "step14_pig_early_canonical_reference_log.txt"),
  file.path(scripts_dir, "step14_canonicalize_pig_early_reference_and_index.R"),
  file.path(scripts_dir, "step14_check_canonicalize_pig_early_reference_and_index.R")
)

cleanup_step14 <- function() {
  existing <- step14_targets[file.exists(step14_targets)]
  if (length(existing) > 0) file.remove(existing)
}

cleanup_step14()

step14_main <- function() {

  normalize_slash <- function(x) {
    gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
  }

  selection_file <- file.path(tables_dir, "step13_pig_early_reference_selection_summary.csv")
  fasta_candidates_file <- file.path(tables_dir, "step13_pig_early_fasta_candidates.csv")
  gtf_candidates_file <- file.path(tables_dir, "step13_pig_early_gtf_candidates.csv")
  idx_candidates_file <- file.path(tables_dir, "step13_pig_early_subread_index_candidates.csv")

  needed <- c(selection_file, fasta_candidates_file, gtf_candidates_file, idx_candidates_file)
  if (!all(file.exists(needed))) {
    stop("Step 13 outputs are missing. Please rerun Step 13 first.")
  }

  selection_df <- read.csv(selection_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  fasta_df <- read.csv(fasta_candidates_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  gtf_df <- read.csv(gtf_candidates_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  idx_df <- read.csv(idx_candidates_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

  fasta_ref_hits <- fasta_df[grepl("^reference/", fasta_df$relative_path), , drop = FALSE]
  if (nrow(fasta_ref_hits) > 0) {
    fasta_pick <- fasta_ref_hits[1, , drop = FALSE]
  } else if (nrow(fasta_df) > 0) {
    fasta_pick <- fasta_df[1, , drop = FALSE]
  } else {
    stop("No FASTA candidate found.")
  }

  gtf_ref_hits <- gtf_df[grepl("^reference/", gtf_df$relative_path), , drop = FALSE]
  if (nrow(gtf_ref_hits) > 0) {
    gtf_pick <- gtf_ref_hits[1, , drop = FALSE]
  } else if (nrow(gtf_df) > 0) {
    gtf_pick <- gtf_df[1, , drop = FALSE]
  } else {
    stop("No GTF candidate found.")
  }

  idx_ref_hits <- idx_df[grepl("^reference/", idx_df$index_basename_relative), , drop = FALSE]
  if (nrow(idx_ref_hits) > 0) {
    idx_pick <- idx_ref_hits[1, , drop = FALSE]
  } else if (nrow(idx_df) > 0) {
    idx_pick <- idx_df[1, , drop = FALSE]
  } else {
    idx_pick <- NULL
  }

  fasta_abs <- normalize_slash(fasta_pick$absolute_path[1])
  fasta_rel <- fasta_pick$relative_path[1]

  gtf_abs <- normalize_slash(gtf_pick$absolute_path[1])
  gtf_rel <- gtf_pick$relative_path[1]

  fasta_exists <- file.exists(fasta_abs)
  gtf_exists <- file.exists(gtf_abs)

  if (!fasta_exists) stop("Chosen FASTA does not exist: ", fasta_abs)
  if (!gtf_exists) stop("Chosen GTF does not exist: ", gtf_abs)

  if (!is.null(idx_pick)) {
    idx_base_abs <- normalize_slash(idx_pick$index_basename_absolute[1])
    idx_base_rel <- idx_pick$index_basename_relative[1]

    idx_dir <- dirname(idx_base_abs)
    idx_prefix <- basename(idx_base_abs)

    dir_files <- list.files(idx_dir, full.names = TRUE)
    dir_files <- normalize_slash(dir_files)
    dir_files <- dir_files[file.exists(dir_files)]

    bundle_files_abs <- dir_files[startsWith(basename(dir_files), idx_prefix)]

    bundle_df <- data.frame(
      absolute_path = bundle_files_abs,
      relative_path = sub(paste0("^", normalize_slash(project_root), "/?"), "", bundle_files_abs),
      filename = basename(bundle_files_abs),
      size_bytes = if (length(bundle_files_abs) > 0) file.info(bundle_files_abs)$size else numeric(0),
      stringsAsFactors = FALSE
    )

    has_array_part <- any(grepl("\\.[0-9]{2}\\.[A-Za-z]\\.array$", bundle_df$filename))
    reuse_index <- nrow(bundle_df) > 0 && has_array_part
  } else {
    idx_base_abs <- NA_character_
    idx_base_rel <- NA_character_
    bundle_df <- data.frame(
      absolute_path = character(0),
      relative_path = character(0),
      filename = character(0),
      size_bytes = numeric(0),
      stringsAsFactors = FALSE
    )
    reuse_index <- FALSE
  }

  canonical_manifest_df <- data.frame(
    item = c("genome_fasta", "annotation_gtf", "subread_index_basename"),
    absolute_path = c(fasta_abs, gtf_abs, idx_base_abs),
    relative_path = c(fasta_rel, gtf_rel, idx_base_rel),
    exists = c(file.exists(fasta_abs), file.exists(gtf_abs), ifelse(is.na(idx_base_abs), FALSE, TRUE)),
    stringsAsFactors = FALSE
  )

  decision_summary_df <- data.frame(
    metric = c(
      "canonical_fasta_relative",
      "canonical_gtf_relative",
      "canonical_index_basename_relative",
      "canonical_fasta_exists",
      "canonical_gtf_exists",
      "index_bundle_file_n",
      "index_has_array_part",
      "reuse_existing_index",
      "next_step_branch"
    ),
    value = c(
      fasta_rel,
      gtf_rel,
      idx_base_rel,
      file.exists(fasta_abs),
      file.exists(gtf_abs),
      nrow(bundle_df),
      ifelse(nrow(bundle_df) > 0, any(grepl("\\.[0-9]{2}\\.[A-Za-z]\\.array$", bundle_df$filename)), FALSE),
      reuse_index,
      ifelse(reuse_index, "reuse_existing_index", "build_new_index")
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    bundle_df,
    file = file.path(tables_dir, "step14_pig_early_index_bundle_files.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    canonical_manifest_df,
    file = file.path(tables_dir, "step14_pig_early_canonical_reference_manifest.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    decision_summary_df,
    file = file.path(tables_dir, "step14_pig_early_index_reuse_decision_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  save(
    fasta_pick,
    gtf_pick,
    idx_pick,
    bundle_df,
    canonical_manifest_df,
    decision_summary_df,
    file = file.path(objects_dir, "step14_pig_early_canonical_reference_workspace.RData")
  )

  log_lines <- c(
    "Step 14 completed: canonicalized pig-early reference choice and index reuse decision.",
    paste("Canonical FASTA:", fasta_rel),
    paste("Canonical GTF:", gtf_rel),
    paste("Canonical index basename:", ifelse(is.na(idx_base_rel), "NOT_FOUND", idx_base_rel)),
    paste("Index bundle file count:", nrow(bundle_df)),
    paste("Reuse existing index:", reuse_index)
  )

  writeLines(
    log_lines,
    con = file.path(logs_dir, "step14_pig_early_canonical_reference_log.txt"),
    useBytes = TRUE
  )

  writeLines(
    c(
      "# Step 14 operation script was executed from console and archived automatically.",
      "setwd(\"D:/R/ACL ∩ senescence2\")",
      "cat(\"Step 14 outputs are saved under rebuild_submission/02_pig_early/.\\n\")"
    ),
    con = file.path(scripts_dir, "step14_canonicalize_pig_early_reference_and_index.R"),
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
      "cat(\"===== Step 14 key file check =====\\n\")",
      "key_files <- c(",
      "  file.path(objects_dir, \"step14_pig_early_canonical_reference_workspace.RData\"),",
      "  file.path(tables_dir, \"step14_pig_early_index_bundle_files.csv\"),",
      "  file.path(tables_dir, \"step14_pig_early_canonical_reference_manifest.csv\"),",
      "  file.path(tables_dir, \"step14_pig_early_index_reuse_decision_summary.csv\"),",
      "  file.path(logs_dir, \"step14_pig_early_canonical_reference_log.txt\")",
      ")",
      "print(data.frame(file = key_files, exists = file.exists(key_files), size_bytes = ifelse(file.exists(key_files), file.info(key_files)$size, NA), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 14 canonical reference manifest =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step14_pig_early_canonical_reference_manifest.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 14 index reuse decision =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step14_pig_early_index_reuse_decision_summary.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 14 index bundle files (top 20) =====\\n\")",
      "print(head(read.csv(file.path(tables_dir, \"step14_pig_early_index_bundle_files.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"), 20))"
    ),
    con = file.path(scripts_dir, "step14_check_canonicalize_pig_early_reference_and_index.R"),
    useBytes = TRUE
  )
}

tryCatch(
  {
    step14_main()
    cat("Step 14 completed successfully.\n")
    cat("Main outputs:\n")
    print(c(
      file.path(objects_dir, "step14_pig_early_canonical_reference_workspace.RData"),
      file.path(tables_dir, "step14_pig_early_index_bundle_files.csv"),
      file.path(tables_dir, "step14_pig_early_canonical_reference_manifest.csv"),
      file.path(tables_dir, "step14_pig_early_index_reuse_decision_summary.csv")
    ))
  },
  error = function(e) {
    cleanup_step14()
    message("Step 14 failed, and partial output files for this step were cleaned automatically.")
    stop(e)
  }
)