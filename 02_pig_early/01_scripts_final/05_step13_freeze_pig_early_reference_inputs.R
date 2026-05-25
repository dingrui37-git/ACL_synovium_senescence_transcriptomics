# ============================================================
# Step 13: freeze pig-early alignment/counting reference inputs
# Goal:
#   - locate candidate genome FASTA
#   - locate candidate annotation GTF
#   - locate any existing Subread index basename
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

step13_targets <- c(
  file.path(objects_dir, "step13_pig_early_reference_inputs_workspace.RData"),
  file.path(tables_dir, "step13_pig_early_fasta_candidates.csv"),
  file.path(tables_dir, "step13_pig_early_gtf_candidates.csv"),
  file.path(tables_dir, "step13_pig_early_subread_index_candidates.csv"),
  file.path(tables_dir, "step13_pig_early_reference_selection_summary.csv"),
  file.path(logs_dir, "step13_pig_early_reference_inputs_log.txt"),
  file.path(scripts_dir, "step13_freeze_pig_early_reference_inputs.R"),
  file.path(scripts_dir, "step13_check_freeze_pig_early_reference_inputs.R")
)

cleanup_step13 <- function() {
  existing <- step13_targets[file.exists(step13_targets)]
  if (length(existing) > 0) file.remove(existing)
}

cleanup_step13()

step13_main <- function() {

  normalize_slash <- function(x) {
    gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
  }

  make_relative <- function(x, root) {
    x_norm <- normalize_slash(x)
    root_norm <- normalize_slash(root)
    prefix <- paste0(root_norm, "/")
    out <- ifelse(startsWith(x_norm, prefix), substring(x_norm, nchar(prefix) + 1), x_norm)
    out
  }

  rank_fasta <- function(path) {
    p <- tolower(path)
    score <- 0
    if (grepl("sus|scrofa|sscrofa|sus_scrofa|pig", p)) score <- score + 100
    if (grepl("genome|dna|toplevel|primary_assembly|assembly", p)) score <- score + 40
    if (grepl("\\.fa$|\\.fasta$|\\.fna$", p)) score <- score + 20
    if (grepl("\\.gz$", p)) score <- score - 5
    if (grepl("fastqc|rebuild_submission|ascii_batch|launcher_onefile|onefile", p)) score <- score - 200
    score
  }

  rank_gtf <- function(path) {
    p <- tolower(path)
    score <- 0
    if (grepl("sus|scrofa|sscrofa|sus_scrofa|pig", p)) score <- score + 100
    if (grepl("gtf|annotation|genes|ensembl|refseq", p)) score <- score + 40
    if (grepl("\\.gtf$", p)) score <- score + 20
    if (grepl("\\.gz$", p)) score <- score - 5
    if (grepl("fastqc|rebuild_submission|ascii_batch|launcher_onefile|onefile", p)) score <- score - 200
    score
  }

  rank_index <- function(path, n_parts, total_size_bytes) {
    p <- tolower(path)
    score <- 0
    if (grepl("sus|scrofa|sscrofa|sus_scrofa|pig", p)) score <- score + 100
    score <- score + n_parts
    score <- score + ifelse(total_size_bytes > 1e8, 20, 0)
    if (grepl("rebuild_submission|fastqc|ascii_batch|launcher_onefile|onefile", p)) score <- score - 200
    score
  }

  all_files <- list.files(
    project_root,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )
  all_files <- all_files[file.exists(all_files)]

  if (length(all_files) == 0) {
    stop("No files found under project root: ", project_root)
  }

  file_df <- data.frame(
    absolute_path = normalize_slash(all_files),
    relative_path = make_relative(all_files, project_root),
    filename = basename(all_files),
    size_bytes = file.info(all_files)$size,
    stringsAsFactors = FALSE
  )

  fasta_keep <- grepl("\\.(fa|fasta|fna)(\\.gz)?$", tolower(file_df$filename))
  fasta_df <- file_df[fasta_keep, c("absolute_path", "relative_path", "filename", "size_bytes")]
  if (nrow(fasta_df) > 0) {
    fasta_df$rank_score <- vapply(fasta_df$absolute_path, rank_fasta, numeric(1))
    fasta_df <- fasta_df[order(-fasta_df$rank_score, -fasta_df$size_bytes, fasta_df$absolute_path), ]
  } else {
    fasta_df <- data.frame(
      absolute_path = character(0),
      relative_path = character(0),
      filename = character(0),
      size_bytes = numeric(0),
      rank_score = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  gtf_keep <- grepl("\\.gtf(\\.gz)?$", tolower(file_df$filename))
  gtf_df <- file_df[gtf_keep, c("absolute_path", "relative_path", "filename", "size_bytes")]
  if (nrow(gtf_df) > 0) {
    gtf_df$rank_score <- vapply(gtf_df$absolute_path, rank_gtf, numeric(1))
    gtf_df <- gtf_df[order(-gtf_df$rank_score, -gtf_df$size_bytes, gtf_df$absolute_path), ]
  } else {
    gtf_df <- data.frame(
      absolute_path = character(0),
      relative_path = character(0),
      filename = character(0),
      size_bytes = numeric(0),
      rank_score = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  idx_keep <- grepl("\\.[0-9]{2}\\.[A-Za-z]\\.array$", file_df$filename)
  idx_parts_df <- file_df[idx_keep, c("absolute_path", "relative_path", "filename", "size_bytes")]

  if (nrow(idx_parts_df) > 0) {
    idx_parts_df$index_basename_absolute <- sub("\\.[0-9]{2}\\.[A-Za-z]\\.array$", "", idx_parts_df$absolute_path)
    idx_parts_df$index_basename_relative <- sub("\\.[0-9]{2}\\.[A-Za-z]\\.array$", "", idx_parts_df$relative_path)

    idx_summary <- aggregate(
      size_bytes ~ index_basename_absolute + index_basename_relative,
      data = idx_parts_df,
      FUN = sum
    )
    part_counts <- aggregate(
      filename ~ index_basename_absolute + index_basename_relative,
      data = idx_parts_df,
      FUN = length
    )
    colnames(part_counts)[colnames(part_counts) == "filename"] <- "n_parts"

    idx_df <- merge(
      idx_summary,
      part_counts,
      by = c("index_basename_absolute", "index_basename_relative"),
      all = TRUE,
      sort = FALSE
    )

    idx_df$rank_score <- mapply(
      FUN = rank_index,
      path = idx_df$index_basename_absolute,
      n_parts = idx_df$n_parts,
      total_size_bytes = idx_df$size_bytes
    )

    idx_df <- idx_df[order(-idx_df$rank_score, -idx_df$n_parts, -idx_df$size_bytes, idx_df$index_basename_absolute), ]
  } else {
    idx_df <- data.frame(
      index_basename_absolute = character(0),
      index_basename_relative = character(0),
      size_bytes = numeric(0),
      n_parts = numeric(0),
      rank_score = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  selected_fasta_abs <- if (nrow(fasta_df) > 0) fasta_df$absolute_path[1] else NA_character_
  selected_fasta_rel <- if (nrow(fasta_df) > 0) fasta_df$relative_path[1] else NA_character_

  selected_gtf_abs <- if (nrow(gtf_df) > 0) gtf_df$absolute_path[1] else NA_character_
  selected_gtf_rel <- if (nrow(gtf_df) > 0) gtf_df$relative_path[1] else NA_character_

  selected_idx_abs <- if (nrow(idx_df) > 0) idx_df$index_basename_absolute[1] else NA_character_
  selected_idx_rel <- if (nrow(idx_df) > 0) idx_df$index_basename_relative[1] else NA_character_
  selected_idx_parts <- if (nrow(idx_df) > 0) idx_df$n_parts[1] else NA

  selection_summary_df <- data.frame(
    metric = c(
      "n_fasta_candidates_found",
      "n_gtf_candidates_found",
      "n_subread_index_candidates_found",
      "selected_fasta_relative",
      "selected_gtf_relative",
      "selected_subread_index_basename_relative",
      "selected_subread_index_n_parts",
      "selected_fasta_exists",
      "selected_gtf_exists",
      "selected_subread_index_detected"
    ),
    value = c(
      nrow(fasta_df),
      nrow(gtf_df),
      nrow(idx_df),
      selected_fasta_rel,
      selected_gtf_rel,
      selected_idx_rel,
      selected_idx_parts,
      !is.na(selected_fasta_abs) && file.exists(selected_fasta_abs),
      !is.na(selected_gtf_abs) && file.exists(selected_gtf_abs),
      !is.na(selected_idx_abs)
    ),
    stringsAsFactors = FALSE
  )

  write.csv(
    fasta_df,
    file = file.path(tables_dir, "step13_pig_early_fasta_candidates.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    gtf_df,
    file = file.path(tables_dir, "step13_pig_early_gtf_candidates.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    idx_df,
    file = file.path(tables_dir, "step13_pig_early_subread_index_candidates.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    selection_summary_df,
    file = file.path(tables_dir, "step13_pig_early_reference_selection_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  save(
    fasta_df,
    gtf_df,
    idx_df,
    selection_summary_df,
    file = file.path(objects_dir, "step13_pig_early_reference_inputs_workspace.RData")
  )

  log_lines <- c(
    "Step 13 completed: froze pig-early alignment/counting reference inputs.",
    paste("FASTA candidates found:", nrow(fasta_df)),
    paste("GTF candidates found:", nrow(gtf_df)),
    paste("Subread index candidates found:", nrow(idx_df)),
    paste("Selected FASTA:", ifelse(is.na(selected_fasta_rel), "NOT_FOUND", selected_fasta_rel)),
    paste("Selected GTF:", ifelse(is.na(selected_gtf_rel), "NOT_FOUND", selected_gtf_rel)),
    paste("Selected Subread index basename:", ifelse(is.na(selected_idx_rel), "NOT_FOUND", selected_idx_rel))
  )

  writeLines(
    log_lines,
    con = file.path(logs_dir, "step13_pig_early_reference_inputs_log.txt"),
    useBytes = TRUE
  )

  writeLines(
    c(
      "# Step 13 operation script was executed from console and archived automatically.",
      "setwd(\"D:/R/ACL ∩ senescence2\")",
      "cat(\"Step 13 outputs are saved under rebuild_submission/02_pig_early/.\\n\")"
    ),
    con = file.path(scripts_dir, "step13_freeze_pig_early_reference_inputs.R"),
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
      "cat(\"===== Step 13 key file check =====\\n\")",
      "key_files <- c(",
      "  file.path(objects_dir, \"step13_pig_early_reference_inputs_workspace.RData\"),",
      "  file.path(tables_dir, \"step13_pig_early_fasta_candidates.csv\"),",
      "  file.path(tables_dir, \"step13_pig_early_gtf_candidates.csv\"),",
      "  file.path(tables_dir, \"step13_pig_early_subread_index_candidates.csv\"),",
      "  file.path(tables_dir, \"step13_pig_early_reference_selection_summary.csv\"),",
      "  file.path(logs_dir, \"step13_pig_early_reference_inputs_log.txt\")",
      ")",
      "print(data.frame(file = key_files, exists = file.exists(key_files), size_bytes = ifelse(file.exists(key_files), file.info(key_files)$size, NA), stringsAsFactors = FALSE))",
      "",
      "cat(\"\\n===== Step 13 selection summary =====\\n\")",
      "print(read.csv(file.path(tables_dir, \"step13_pig_early_reference_selection_summary.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"))",
      "",
      "cat(\"\\n===== Step 13 FASTA candidates (top 10) =====\\n\")",
      "print(head(read.csv(file.path(tables_dir, \"step13_pig_early_fasta_candidates.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"), 10))",
      "",
      "cat(\"\\n===== Step 13 GTF candidates (top 10) =====\\n\")",
      "print(head(read.csv(file.path(tables_dir, \"step13_pig_early_gtf_candidates.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"), 10))",
      "",
      "cat(\"\\n===== Step 13 Subread index candidates (top 10) =====\\n\")",
      "print(head(read.csv(file.path(tables_dir, \"step13_pig_early_subread_index_candidates.csv\"), stringsAsFactors = FALSE, fileEncoding = \"UTF-8\"), 10))"
    ),
    con = file.path(scripts_dir, "step13_check_freeze_pig_early_reference_inputs.R"),
    useBytes = TRUE
  )
}

tryCatch(
  {
    step13_main()
    cat("Step 13 completed successfully.\n")
    cat("Main outputs:\n")
    print(c(
      file.path(objects_dir, "step13_pig_early_reference_inputs_workspace.RData"),
      file.path(tables_dir, "step13_pig_early_fasta_candidates.csv"),
      file.path(tables_dir, "step13_pig_early_gtf_candidates.csv"),
      file.path(tables_dir, "step13_pig_early_subread_index_candidates.csv"),
      file.path(tables_dir, "step13_pig_early_reference_selection_summary.csv")
    ))
  },
  error = function(e) {
    cleanup_step13()
    message("Step 13 failed, and partial output files for this step were cleaned automatically.")
    stop(e)
  }
)
