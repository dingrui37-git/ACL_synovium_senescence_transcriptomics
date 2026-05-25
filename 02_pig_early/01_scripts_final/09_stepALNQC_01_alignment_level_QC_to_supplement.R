# StepALNQC_01: Early pig alignment-level QC from Rsubread BAM and BAM.summary files.
# This script summarizes alignment-level QC for the 18 core early pig synovium samples.
# It uses BAM file availability/size and Rsubread *.bam.summary files only; it does not use
# featureCounts assignment summaries or downstream differential-expression outputs.

method_version <- "2026-05-10_output_to_02_pig_early_supplement_alignment_level_QC"
step_name <- "StepALNQC_01_pig_early_alignment_level_QC"

input_bam_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align"
output_root <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/alignment_level_QC"

subdirs <- c("tables", "source_data", "logs", "objects", "scripts")
for (d in file.path(output_root, subdirs)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(output_root, "logs", "StepALNQC_01_pig_early_alignment_level_QC_summary_log.txt")

sink(log_file, split = TRUE)
status <- "FAILED"
start_time <- Sys.time()

cat("============================================================\n")
cat(step_name, "\n")
cat("Started at:", format(start_time), "\n")
cat("Method version:", method_version, "\n")
cat("============================================================\n\n")

cleanup_previous_outputs <- function() {
  cat("Cleaning previous StepALNQC_01 outputs from output_root subfolders...\n")
  removed <- 0L
  for (sd in subdirs) {
    dd <- file.path(output_root, sd)
    if (dir.exists(dd)) {
      old <- list.files(dd, pattern = "^StepALNQC_01_", full.names = TRUE)
      if (length(old) > 0) {
        ok <- file.remove(old)
        removed <- removed + sum(ok, na.rm = TRUE)
      }
    }
  }
  cat("Number of previous StepALNQC_01 files removed:", removed, "\n\n")
}

cleanup_previous_outputs()

tryCatch({
  suppressPackageStartupMessages({
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required but not installed.")
  })

  cat("Required packages loaded.\n")
  cat("dplyr version:", as.character(packageVersion("dplyr")), "\n\n")

  cat("Input BAM directory:\n")
  cat(input_bam_dir, "\n\n")
  cat("Output root:\n")
  cat(output_root, "\n\n")

  if (!dir.exists(input_bam_dir)) {
    stop("Input BAM directory does not exist: ", input_bam_dir)
  }

  expected_samples <- c(paste0("CON", 1:6), paste0("INJS", 1:6), paste0("INJL", 1:6))

  infer_group <- function(sample_id) {
    if (grepl("^CON", sample_id)) return("Control")
    if (grepl("^INJS", sample_id)) return("ACLT_1W")
    if (grepl("^INJL", sample_id)) return("ACLT_4W")
    return(NA_character_)
  }

  infer_time_point <- function(sample_id) {
    if (grepl("^CON", sample_id)) return("t0")
    if (grepl("^INJS", sample_id)) return("t7")
    if (grepl("^INJL", sample_id)) return("t28")
    return(NA_character_)
  }

  sample_number <- function(sample_id) {
    as.integer(sub("^[A-Za-z]+", "", sample_id))
  }

  sample_metadata <- data.frame(
    sample_id = expected_samples,
    group = factor(vapply(expected_samples, infer_group, character(1)),
                   levels = c("Control", "ACLT_1W", "ACLT_4W")),
    time_point = vapply(expected_samples, infer_time_point, character(1)),
    sample_number = vapply(expected_samples, sample_number, integer(1)),
    stringsAsFactors = FALSE
  )

  cat("Expected 18 core samples:\n")
  print(sample_metadata)
  cat("\nGroup counts:\n")
  print(table(sample_metadata$group))
  cat("\n")

  bam_files_all <- list.files(input_bam_dir, pattern = "\\.bam$", full.names = TRUE, ignore.case = TRUE)
  bam_summary_files_all <- list.files(input_bam_dir, pattern = "\\.bam\\.summary$", full.names = TRUE, ignore.case = TRUE)
  indel_vcf_files_all <- list.files(input_bam_dir, pattern = "\\.bam\\.indel\\.vcf$", full.names = TRUE, ignore.case = TRUE)

  bam_sample_ids_all <- sub("\\.bam$", "", basename(bam_files_all), ignore.case = TRUE)
  summary_sample_ids_all <- sub("\\.bam\\.summary$", "", basename(bam_summary_files_all), ignore.case = TRUE)

  cat("File discovery summary:\n")
  cat("BAM files found:", length(bam_files_all), "\n")
  cat("BAM.summary files found:", length(bam_summary_files_all), "\n")
  cat("BAM.indel.vcf files found but not used:", length(indel_vcf_files_all), "\n\n")

  extra_bam <- setdiff(bam_sample_ids_all, expected_samples)
  extra_summary <- setdiff(summary_sample_ids_all, expected_samples)
  if (length(extra_bam) > 0) {
    cat("Extra BAM sample IDs not in 18 core samples:\n")
    print(extra_bam)
  } else {
    cat("No extra BAM sample IDs outside the 18 core samples were detected.\n")
  }
  if (length(extra_summary) > 0) {
    cat("Extra BAM.summary sample IDs not in 18 core samples:\n")
    print(extra_summary)
  } else {
    cat("No extra BAM.summary sample IDs outside the 18 core samples were detected.\n")
  }
  cat("\n")

  first_count_like_number <- function(line) {
    line2 <- gsub(",", "", line)
    nums <- regmatches(line2, gregexpr("[0-9]+(?:\\.[0-9]+)?", line2, perl = TRUE))[[1]]
    if (length(nums) == 0) return(NA_real_)
    values <- suppressWarnings(as.numeric(nums))
    values <- values[!is.na(values)]
    if (length(values) == 0) return(NA_real_)
    large_values <- values[values >= 100]
    if (length(large_values) > 0) return(large_values[1])
    return(values[1])
  }

  extract_candidates <- function(lines, include_patterns, exclude_patterns = character(0)) {
    if (length(lines) == 0) return(numeric(0))
    lc <- tolower(lines)
    keep <- rep(TRUE, length(lines))
    for (p in include_patterns) keep <- keep & grepl(p, lc, perl = TRUE)
    if (length(exclude_patterns) > 0) {
      for (p in exclude_patterns) keep <- keep & !grepl(p, lc, perl = TRUE)
    }
    vals <- vapply(lines[keep], first_count_like_number, numeric(1))
    vals <- vals[!is.na(vals)]
    vals
  }

  parse_bam_summary <- function(summary_file) {
    if (!file.exists(summary_file)) {
      return(list(
        total_fragments = NA_real_, mapped_fragments = NA_real_, unmapped_fragments = NA_real_,
        mapping_rate_pct = NA_real_, parse_status = "MISSING_SUMMARY",
        parse_note = "Summary file does not exist"
      ))
    }

    lines <- readLines(summary_file, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]

    total_candidates <- c(
      extract_candidates(lines,
                         include_patterns = c("total|input", "read|fragment|pair"),
                         exclude_patterns = c("mapped|unmapped|not mapped|not-mapped|not aligned|assigned|alignment rate|bases|base")),
      extract_candidates(lines,
                         include_patterns = c("number of input"),
                         exclude_patterns = c("bases|base"))
    )

    mapped_total_candidates <- extract_candidates(lines,
                                                  include_patterns = c("mapped|aligned|alignment"),
                                                  exclude_patterns = c("unmapped|not mapped|not-mapped|not aligned|not-aligned|unique|uniquely|multi|multiple|multi-mapping|multimapping|rate|%"))

    successful_mapped_candidates <- extract_candidates(lines,
                                                       include_patterns = c("success"),
                                                       exclude_patterns = c("rate|%"))

    unique_candidates <- extract_candidates(lines,
                                            include_patterns = c("unique|uniquely", "mapped|aligned"),
                                            exclude_patterns = c("rate|%"))

    multi_candidates <- extract_candidates(lines,
                                           include_patterns = c("multi|multiple", "mapped|aligned|mapping"),
                                           exclude_patterns = c("rate|%"))

    unmapped_candidates <- extract_candidates(lines,
                                              include_patterns = c("unmapped|not mapped|not-mapped|not aligned|not-aligned|failed"),
                                              exclude_patterns = c("rate|%"))

    # Use the largest count candidate to avoid accidentally using small numbers from sample IDs or percentages.
    total <- if (length(total_candidates) > 0) max(total_candidates, na.rm = TRUE) else NA_real_

    mapped <- NA_real_
    if (length(mapped_total_candidates) > 0) {
      mapped <- max(mapped_total_candidates, na.rm = TRUE)
    } else if (length(successful_mapped_candidates) > 0) {
      mapped <- max(successful_mapped_candidates, na.rm = TRUE)
    } else if (length(unique_candidates) > 0 || length(multi_candidates) > 0) {
      mapped <- sum(c(unique_candidates, multi_candidates), na.rm = TRUE)
    }

    unmapped <- if (length(unmapped_candidates) > 0) max(unmapped_candidates, na.rm = TRUE) else NA_real_

    if (is.na(mapped) && !is.na(total) && !is.na(unmapped)) mapped <- total - unmapped
    if (is.na(unmapped) && !is.na(total) && !is.na(mapped)) unmapped <- total - mapped
    if (is.na(total) && !is.na(mapped) && !is.na(unmapped)) total <- mapped + unmapped

    # Guard against invalid negative values caused by unusual summary formatting.
    if (!is.na(mapped) && mapped < 0) mapped <- NA_real_
    if (!is.na(unmapped) && unmapped < 0) unmapped <- NA_real_
    if (!is.na(total) && total < 0) total <- NA_real_

    mapping_rate <- if (!is.na(total) && total > 0 && !is.na(mapped)) 100 * mapped / total else NA_real_

    parse_status <- if (!is.na(total) && !is.na(mapped) && !is.na(unmapped) && !is.na(mapping_rate)) "PASS" else "CHECK"
    parse_note <- if (parse_status == "PASS") {
      "Parsed total/mapped/unmapped counts from Rsubread BAM.summary file"
    } else {
      "Could not confidently parse one or more of total/mapped/unmapped counts; inspect raw summary lines"
    }

    list(
      total_fragments = total,
      mapped_fragments = mapped,
      unmapped_fragments = unmapped,
      mapping_rate_pct = mapping_rate,
      parse_status = parse_status,
      parse_note = parse_note
    )
  }

  bam_path_for_sample <- function(sample_id) {
    f <- file.path(input_bam_dir, paste0(sample_id, ".bam"))
    if (file.exists(f)) return(f)
    candidates <- bam_files_all[toupper(bam_sample_ids_all) == toupper(sample_id)]
    if (length(candidates) > 0) return(candidates[1])
    return(NA_character_)
  }

  summary_path_for_sample <- function(sample_id) {
    f <- file.path(input_bam_dir, paste0(sample_id, ".bam.summary"))
    if (file.exists(f)) return(f)
    candidates <- bam_summary_files_all[toupper(summary_sample_ids_all) == toupper(sample_id)]
    if (length(candidates) > 0) return(candidates[1])
    return(NA_character_)
  }

  qc_rows <- lapply(expected_samples, function(sid) {
    bam_file <- bam_path_for_sample(sid)
    summary_file <- summary_path_for_sample(sid)

    bam_exists <- !is.na(bam_file) && file.exists(bam_file)
    summary_exists <- !is.na(summary_file) && file.exists(summary_file)

    bam_size_bytes <- if (bam_exists) file.info(bam_file)$size else NA_real_
    bam_size_mb <- if (!is.na(bam_size_bytes)) bam_size_bytes / 1024^2 else NA_real_

    parsed <- if (summary_exists) parse_bam_summary(summary_file) else parse_bam_summary("__missing__")

    data.frame(
      sample_id = sid,
      group = infer_group(sid),
      time_point = infer_time_point(sid),
      sample_number = sample_number(sid),
      bam_exists = bam_exists,
      bam_file = ifelse(bam_exists, bam_file, NA_character_),
      bam_size_bytes = bam_size_bytes,
      bam_size_mb = round(bam_size_mb, 3),
      bam_summary_exists = summary_exists,
      bam_summary_file = ifelse(summary_exists, summary_file, NA_character_),
      total_fragments = parsed$total_fragments,
      mapped_fragments = parsed$mapped_fragments,
      unmapped_fragments = parsed$unmapped_fragments,
      mapping_rate_pct = round(parsed$mapping_rate_pct, 3),
      parse_status = parsed$parse_status,
      parse_note = parsed$parse_note,
      stringsAsFactors = FALSE
    )
  })

  qc_table <- do.call(rbind, qc_rows)
  qc_table$group <- factor(qc_table$group, levels = c("Control", "ACLT_1W", "ACLT_4W"))
  qc_table <- qc_table[order(qc_table$group, qc_table$sample_number), ]

  cat("Alignment-level QC table preview:\n")
  print(qc_table[, c("sample_id", "group", "time_point", "bam_exists", "bam_size_mb", "bam_summary_exists", "total_fragments", "mapped_fragments", "unmapped_fragments", "mapping_rate_pct", "parse_status")], row.names = FALSE)
  cat("\n")

  # Raw summary lines for audit.
  raw_lines_list <- lapply(expected_samples, function(sid) {
    summary_file <- summary_path_for_sample(sid)
    if (!is.na(summary_file) && file.exists(summary_file)) {
      lines <- readLines(summary_file, warn = FALSE)
      data.frame(
        sample_id = sid,
        bam_summary_file = summary_file,
        line_number = seq_along(lines),
        raw_line = lines,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        sample_id = sid,
        bam_summary_file = NA_character_,
        line_number = NA_integer_,
        raw_line = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  })
  raw_summary_lines <- do.call(rbind, raw_lines_list)

  group_summary <- qc_table |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      n_samples = dplyr::n(),
      bam_exists_n = sum(bam_exists, na.rm = TRUE),
      bam_summary_exists_n = sum(bam_summary_exists, na.rm = TRUE),
      parsed_pass_n = sum(parse_status == "PASS", na.rm = TRUE),
      mean_bam_size_mb = round(mean(bam_size_mb, na.rm = TRUE), 3),
      min_bam_size_mb = round(min(bam_size_mb, na.rm = TRUE), 3),
      max_bam_size_mb = round(max(bam_size_mb, na.rm = TRUE), 3),
      mean_total_fragments = round(mean(total_fragments, na.rm = TRUE), 0),
      mean_mapped_fragments = round(mean(mapped_fragments, na.rm = TRUE), 0),
      mean_unmapped_fragments = round(mean(unmapped_fragments, na.rm = TRUE), 0),
      mean_mapping_rate_pct = round(mean(mapping_rate_pct, na.rm = TRUE), 3),
      min_mapping_rate_pct = round(min(mapping_rate_pct, na.rm = TRUE), 3),
      max_mapping_rate_pct = round(max(mapping_rate_pct, na.rm = TRUE), 3),
      .groups = "drop"
    )

  overall_summary <- data.frame(
    metric = c(
      "method_version",
      "input_bam_dir",
      "output_root",
      "expected_core_samples",
      "bam_files_found_all",
      "bam_summary_files_found_all",
      "core_bam_files_present",
      "core_bam_summary_files_present",
      "core_bam_summary_files_parsed_PASS",
      "mean_bam_size_mb",
      "min_bam_size_mb",
      "max_bam_size_mb",
      "mean_total_fragments",
      "mean_mapped_fragments",
      "mean_unmapped_fragments",
      "mean_mapping_rate_pct",
      "min_mapping_rate_pct",
      "max_mapping_rate_pct",
      "failed_or_incomplete_samples"
    ),
    value = c(
      method_version,
      input_bam_dir,
      output_root,
      length(expected_samples),
      length(bam_files_all),
      length(bam_summary_files_all),
      sum(qc_table$bam_exists, na.rm = TRUE),
      sum(qc_table$bam_summary_exists, na.rm = TRUE),
      sum(qc_table$parse_status == "PASS", na.rm = TRUE),
      round(mean(qc_table$bam_size_mb, na.rm = TRUE), 3),
      round(min(qc_table$bam_size_mb, na.rm = TRUE), 3),
      round(max(qc_table$bam_size_mb, na.rm = TRUE), 3),
      round(mean(qc_table$total_fragments, na.rm = TRUE), 0),
      round(mean(qc_table$mapped_fragments, na.rm = TRUE), 0),
      round(mean(qc_table$unmapped_fragments, na.rm = TRUE), 0),
      round(mean(qc_table$mapping_rate_pct, na.rm = TRUE), 3),
      round(min(qc_table$mapping_rate_pct, na.rm = TRUE), 3),
      round(max(qc_table$mapping_rate_pct, na.rm = TRUE), 3),
      paste(qc_table$sample_id[!(qc_table$bam_exists & qc_table$bam_summary_exists & qc_table$parse_status == "PASS")], collapse = "; ")
    ),
    stringsAsFactors = FALSE
  )
  overall_summary$value[overall_summary$metric == "failed_or_incomplete_samples" & overall_summary$value == ""] <- "none"

  method_parameters <- data.frame(
    parameter = c(
      "step_name",
      "method_version",
      "input_bam_dir",
      "output_root",
      "expected_samples",
      "files_used",
      "files_not_used",
      "sample_group_rule",
      "mapping_rate_formula"
    ),
    value = c(
      step_name,
      method_version,
      input_bam_dir,
      output_root,
      paste(expected_samples, collapse = ", "),
      "*.bam for file availability/size; *.bam.summary for Rsubread alignment fragment/read summary",
      "*.bam.indel.vcf files; featureCounts assignment summary files; downstream DE/signature/core-gene tables",
      "CON=Control/t0; INJS=ACLT_1W/t7; INJL=ACLT_4W/t28",
      "mapped_fragments / total_fragments * 100"
    ),
    stringsAsFactors = FALSE
  )

  methods_text <- paste(
    "Alignment-level QC for the early pig validation dataset was summarized from the Rsubread alignment outputs of the 18 core synovium samples.",
    "For each sample, BAM file availability and BAM file size were recorded, and the corresponding Rsubread .bam.summary file was parsed to extract total fragments, mapped fragments, unmapped fragments, and mapping rate.",
    "The 18 samples consisted of six Control samples, six ACLT untreated 1-week samples, and six ACLT untreated 4-week samples.",
    "This table reports alignment-level QC only; featureCounts assignment counts and assignment rates were not included because the retained .bam.summary files reflect the alignment summary rather than gene-level feature assignment.",
    sep = "\n"
  )

  # Save outputs.
  qc_table_file <- file.path(output_root, "tables", "StepALNQC_01_pig_early_alignment_level_QC_table.csv")
  supp_table_file <- file.path(output_root, "tables", "StepALNQC_01_Supplementary_Table_S3_alignment_level_QC.csv")
  group_summary_file <- file.path(output_root, "tables", "StepALNQC_01_group_alignment_level_QC_summary.csv")
  overall_summary_file <- file.path(output_root, "tables", "StepALNQC_01_final_summary_for_review.csv")
  method_parameters_file <- file.path(output_root, "tables", "StepALNQC_01_method_parameters.csv")
  methods_text_file <- file.path(output_root, "tables", "StepALNQC_01_methods_text_alignment_level_QC.txt")
  raw_lines_file <- file.path(output_root, "source_data", "StepALNQC_01_raw_bam_summary_lines_for_audit.csv")
  rds_file <- file.path(output_root, "objects", "StepALNQC_01_alignment_level_QC_objects.rds")

  write.csv(qc_table, qc_table_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(qc_table[, c("sample_id", "group", "time_point", "bam_exists", "bam_size_mb", "total_fragments", "mapped_fragments", "unmapped_fragments", "mapping_rate_pct", "parse_status", "parse_note")],
            supp_table_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(group_summary, group_summary_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(overall_summary, overall_summary_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(method_parameters, method_parameters_file, row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(methods_text, methods_text_file, useBytes = TRUE)
  write.csv(raw_summary_lines, raw_lines_file, row.names = FALSE, fileEncoding = "UTF-8")
  saveRDS(list(qc_table = qc_table, group_summary = group_summary, overall_summary = overall_summary, method_parameters = method_parameters, raw_summary_lines = raw_summary_lines), rds_file)

  cat("Group-level alignment QC summary:\n")
  print(group_summary)
  cat("\n")

  cat("Final summary for review:\n")
  print(overall_summary, row.names = FALSE)
  cat("\n")

  cat("Key output files:\n")
  cat("1) Supplementary Table S3-style QC table:\n", supp_table_file, "\n", sep = "")
  cat("2) Full QC table with file paths:\n", qc_table_file, "\n", sep = "")
  cat("3) Group-level QC summary:\n", group_summary_file, "\n", sep = "")
  cat("4) Raw BAM.summary lines for audit:\n", raw_lines_file, "\n", sep = "")
  cat("5) Method parameters:\n", method_parameters_file, "\n", sep = "")
  cat("6) Methods text:\n", methods_text_file, "\n", sep = "")
  cat("7) Final summary table:\n", overall_summary_file, "\n\n", sep = "")

  cat("Session information:\n")
  print(sessionInfo())

  status <- "SUCCESS"
}, error = function(e) {
  cat("\nERROR encountered during ", step_name, ":\n", sep = "")
  cat(conditionMessage(e), "\n")
  status <<- "FAILED"
}, finally = {
  end_time <- Sys.time()
  cat("\n============================================================\n")
  cat(step_name, "finished with status:", status, "\n")
  cat("Finished at:", format(end_time), "\n")
  cat("Summary log saved to:\n")
  cat(log_file, "\n")
  cat("============================================================\n")
  sink()
})

if (!identical(status, "SUCCESS")) {
  stop(step_name, " failed. Please inspect the summary log: ", log_file)
}
