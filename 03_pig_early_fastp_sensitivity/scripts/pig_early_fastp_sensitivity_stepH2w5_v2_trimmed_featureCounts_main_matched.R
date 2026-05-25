# ============================================================
# Step H2w5 v2
# pig early fastp sensitivity branch:
# Run featureCounts on trimmed BAM files with a writable temp/work directory
#
# 这一步做什么：
# 1) 读取 H2w4 v2 生成的 trimmed BAM 状态表；
# 2) 检查 18 个 trimmed BAM 是否存在；
# 3) 使用 non-stranded paired-end 参数重新运行 featureCounts；
# 4) 生成 trimmed 分支的 gene-level count matrix；
# 5) 修复上一版报错：ERROR: temporary directory is not writable: '.'
# 6) 修正 v2 参数记录不完整的问题：显式采用与主分析一致的 featureCounts 参数
#    useMetaFeatures=TRUE, countReadPairs=TRUE, requireBothEndsMapped=TRUE, checkFragLength=FALSE
# ============================================================

options(stringsAsFactors = FALSE)

# ---------- Basic paths ----------
project_root <- "E:/R/ACLsenescence2"
branch_dir <- file.path(project_root, "rebuild_submission", "02_pig_early_fastp_sensitivity")

tables_dir <- file.path(branch_dir, "tables")
logs_dir <- file.path(branch_dir, "logs")
objects_dir <- file.path(branch_dir, "objects")
ref_cache_dir <- file.path(branch_dir, "reference_cache")
tmp_dir <- file.path(branch_dir, "featureCounts_tmp")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(logs_dir, "stepH2w5_v2_trimmed_featureCounts_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste(...)
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

stop_with_log <- function(...) {
  msg <- paste(...)
  log_msg("ERROR:", msg)
  stop(msg, call. = FALSE)
}

# ---------- Clean old Step H2w5 outputs from earlier parameter-incomplete runs ----------
# This corrected script keeps the same Step H2w5 v2 naming, but removes previous H2w5
# outputs before rerunning so that old featureCounts results cannot be mixed with the
# corrected main-analysis-matched parameter set. Upstream H2w4 BAM files are not removed.
old_h2w5_outputs <- c(
  file.path(tables_dir, "stepH2w5_v2_trimmed_bam_input_check_18samples.csv"),
  file.path(tables_dir, "stepH2w5_v2_featureCounts_arguments.csv"),
  file.path(tables_dir, "stepH2w5_v2_trimmed_gene_level_counts_matrix.csv"),
  file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_annotation.csv"),
  file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_stat.csv"),
  file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_sample_assignment_summary.csv"),
  file.path(tables_dir, "stepH2w5_v2_trimmed_sample_info_for_DE.csv"),
  file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_overall_summary.csv"),
  file.path(objects_dir, "stepH2w5_v2_trimmed_featureCounts_raw_object.rds"),
  file.path(objects_dir, "stepH2w5_v2_trimmed_featureCounts_config.rds")
)
old_h2w5_outputs <- old_h2w5_outputs[file.exists(old_h2w5_outputs)]
if (length(old_h2w5_outputs) > 0) {
  unlink(old_h2w5_outputs, force = TRUE)
}

log_msg("===== Step H2w5 v2 started =====")
log_msg("Corrected featureCounts parameter set: explicitly matches the primary pig early featureCounts settings.")
if (length(old_h2w5_outputs) > 0) {
  log_msg("Removed old H2w5 output files before rerun:", paste(old_h2w5_outputs, collapse = "; "))
}

log_msg("Project root:", project_root)
log_msg("Branch dir:", branch_dir)

# ---------- Package check ----------
if (!requireNamespace("Rsubread", quietly = TRUE)) {
  stop_with_log("Rsubread package is not installed or not available in this R session.")
}

# ---------- Make sure temporary/work directory is writable ----------
# featureCounts 有时会在当前工作目录 '.' 写临时文件。
# 所以这里显式切换到一个可写 tmp_dir，并设置 TMP/TEMP/TMPDIR。
test_file <- file.path(tmp_dir, paste0("write_test_", format(Sys.time(), "%Y%m%d%H%M%S"), ".tmp"))
write_ok <- tryCatch({
  cat("test\n", file = test_file)
  file.exists(test_file)
}, error = function(e) FALSE)
if (!write_ok) {
  stop_with_log("Temporary directory is not writable:", tmp_dir)
}
unlink(test_file, force = TRUE)

Sys.setenv(TMPDIR = tmp_dir, TMP = tmp_dir, TEMP = tmp_dir)
log_msg("Writable featureCounts temp/work dir:", tmp_dir)

# ---------- Input files ----------
h2w4_status_file <- file.path(tables_dir, "stepH2w4_v2_trimmed_alignment_status_18samples.csv")
if (!file.exists(h2w4_status_file)) {
  stop_with_log("Cannot find H2w4 v2 status file:", h2w4_status_file)
}

status_df <- read.csv(h2w4_status_file, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("sample_id", "group", "sample_ok_h2w4")
missing_required <- setdiff(required_cols, names(status_df))
if (length(missing_required) > 0) {
  stop_with_log(
    "H2w4 status file lacks required columns:",
    paste(missing_required, collapse = ", "),
    "\nCurrent columns:", paste(names(status_df), collapse = ", ")
  )
}

if (nrow(status_df) != 18) {
  log_msg("WARNING: H2w4 status file has", nrow(status_df), "rows, expected 18.")
}

ok_h2w4 <- as.logical(status_df$sample_ok_h2w4)
if (!all(ok_h2w4)) {
  bad <- status_df$sample_id[!ok_h2w4]
  stop_with_log("Some H2w4 samples are not OK:", paste(bad, collapse = ", "))
}

# ---------- Locate BAM files ----------
bam_dir <- file.path(branch_dir, "bam_trimmed")
if (!dir.exists(bam_dir)) {
  stop_with_log("Cannot find trimmed BAM directory:", bam_dir)
}

# 优先使用 status 表里已有的 BAM 路径列；如果没有，则从 bam_trimmed 文件夹里自动匹配 sample_id。
possible_bam_cols <- c(
  "bam_file", "bam_path", "bam", "output_bam", "bam_output",
  "out_bam", "BAM", "bam_file_after"
)
bam_col <- intersect(possible_bam_cols, names(status_df))
bam_col <- if (length(bam_col) > 0) bam_col[1] else NA_character_

all_bams <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE, recursive = FALSE)

escape_regex <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

find_bam_for_sample <- function(sample_id) {
  if (!is.na(bam_col)) {
    candidate <- status_df[status_df$sample_id == sample_id, bam_col][1]
    if (!is.na(candidate) && nzchar(candidate) && file.exists(candidate)) return(candidate)
  }
  base_names <- basename(all_bams)
  sid <- escape_regex(sample_id)
  idx1 <- which(grepl(paste0("^", sid, "([_.-]|$)"), base_names))
  if (length(idx1) == 1) return(all_bams[idx1])
  if (length(idx1) > 1) return(all_bams[idx1[1]])
  idx2 <- which(grepl(sid, base_names))
  if (length(idx2) == 1) return(all_bams[idx2])
  if (length(idx2) > 1) return(all_bams[idx2[1]])
  NA_character_
}

status_df$bam_file_h2w5 <- vapply(status_df$sample_id, find_bam_for_sample, character(1))
status_df$bam_exists_h2w5 <- file.exists(status_df$bam_file_h2w5)
status_df$bam_size_h2w5 <- ifelse(
  status_df$bam_exists_h2w5,
  file.info(status_df$bam_file_h2w5)$size,
  NA_real_
)

bam_check_file <- file.path(tables_dir, "stepH2w5_v2_trimmed_bam_input_check_18samples.csv")
write.csv(status_df, bam_check_file, row.names = FALSE, fileEncoding = "UTF-8")

if (!all(status_df$bam_exists_h2w5)) {
  bad <- status_df$sample_id[!status_df$bam_exists_h2w5]
  stop_with_log(
    "Some BAM files cannot be located:",
    paste(bad, collapse = ", "),
    "\nPlease check:", bam_check_file
  )
}

if (any(is.na(status_df$bam_size_h2w5) | status_df$bam_size_h2w5 <= 1000000)) {
  bad <- status_df$sample_id[is.na(status_df$bam_size_h2w5) | status_df$bam_size_h2w5 <= 1000000]
  stop_with_log(
    "Some BAM files are too small or invalid:",
    paste(bad, collapse = ", "),
    "\nPlease check:", bam_check_file
  )
}

bam_files <- normalizePath(status_df$bam_file_h2w5, winslash = "/", mustWork = TRUE)
sample_ids <- status_df$sample_id
groups <- status_df$group

log_msg("BAM files found:", length(bam_files))
log_msg("First BAM:", bam_files[1])

# ---------- GTF ----------
gtf_gz <- file.path(project_root, "reference", "Sus_scrofa_Ensembl115", "Sus_scrofa.Sscrofa11.1.115.gtf.gz")
gtf_plain <- file.path(ref_cache_dir, "Sus_scrofa.Sscrofa11.1.115.gtf")

if (!file.exists(gtf_plain)) {
  if (!file.exists(gtf_gz)) {
    stop_with_log("Cannot find GTF .gz:", gtf_gz)
  }
  log_msg("Decompressing GTF:", gtf_gz)
  gz_con <- gzfile(gtf_gz, open = "rb")
  out_con <- file(gtf_plain, open = "wb")
  ok <- tryCatch({
    repeat {
      x <- readBin(gz_con, what = raw(), n = 1024 * 1024)
      if (length(x) == 0) break
      writeBin(x, out_con)
    }
    TRUE
  }, error = function(e) {
    log_msg("GTF decompression error:", conditionMessage(e))
    FALSE
  })
  close(gz_con)
  close(out_con)
  if (!ok || !file.exists(gtf_plain)) {
    stop_with_log("Failed to decompress GTF.")
  }
} else {
  log_msg("Using cached GTF:", gtf_plain)
}

gtf_for_featureCounts <- normalizePath(gtf_plain, winslash = "/", mustWork = TRUE)
log_msg("GTF for featureCounts:", gtf_for_featureCounts)

# ---------- Run featureCounts ----------
threads <- 4

featurecounts_args <- list(
  files = bam_files,
  annot.ext = gtf_for_featureCounts,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  isPairedEnd = TRUE,
  countReadPairs = TRUE,
  requireBothEndsMapped = TRUE,
  checkFragLength = FALSE,
  strandSpecific = 0,
  nthreads = threads
)

fc_formals <- names(formals(Rsubread::featureCounts))
requested_featurecounts_arg_names <- names(featurecounts_args)
featurecounts_args <- featurecounts_args[names(featurecounts_args) %in% fc_formals]
missing_from_this_Rsubread <- setdiff(requested_featurecounts_arg_names, names(featurecounts_args))
if (length(missing_from_this_Rsubread) > 0) {
  log_msg(
    "WARNING: These requested featureCounts arguments are not available in this Rsubread version and were dropped:",
    paste(missing_from_this_Rsubread, collapse = ", ")
  )
}
if ("autosort" %in% fc_formals) featurecounts_args$autosort <- TRUE

required_main_args <- c("useMetaFeatures", "countReadPairs", "requireBothEndsMapped", "checkFragLength")
missing_required_main_args <- setdiff(required_main_args, names(featurecounts_args))
if (length(missing_required_main_args) > 0) {
  stop_with_log(
    "The current Rsubread::featureCounts() does not expose required main-analysis parameters:",
    paste(missing_required_main_args, collapse = ", "),
    ". Please use the same Rsubread version/environment as the main pig early analysis."
  )
}

featurecounts_args_file <- file.path(tables_dir, "stepH2w5_v2_featureCounts_arguments.csv")
write.csv(
  data.frame(
    argument = names(featurecounts_args),
    value = vapply(featurecounts_args, function(x) paste(x, collapse = "; "), character(1))
  ),
  featurecounts_args_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

log_msg("Running featureCounts in writable work dir. This may take a while ...")

old_wd <- getwd()
setwd(tmp_dir)
fc <- tryCatch({
  do.call(Rsubread::featureCounts, featurecounts_args)
}, error = function(e) {
  setwd(old_wd)
  stop_with_log("featureCounts failed:", conditionMessage(e))
})
setwd(old_wd)

log_msg("featureCounts finished.")

# ---------- Save raw RDS ----------
fc_rds <- file.path(objects_dir, "stepH2w5_v2_trimmed_featureCounts_raw_object.rds")
saveRDS(fc, fc_rds)

# ---------- Count matrix ----------
if (is.null(fc$counts)) {
  stop_with_log("featureCounts returned NULL counts.")
}

count_mat <- fc$counts
colnames(count_mat) <- sample_ids

count_df <- data.frame(
  gene_id = rownames(count_mat),
  count_mat,
  check.names = FALSE
)

counts_file <- file.path(tables_dir, "stepH2w5_v2_trimmed_gene_level_counts_matrix.csv")
write.csv(count_df, counts_file, row.names = FALSE, fileEncoding = "UTF-8")

# Annotation table from featureCounts
annotation_file <- file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_annotation.csv")
if (!is.null(fc$annotation)) {
  anno_df <- as.data.frame(fc$annotation, check.names = FALSE, stringsAsFactors = FALSE)
  anno_df$gene_id <- rownames(anno_df)
  anno_df <- anno_df[, c("gene_id", setdiff(names(anno_df), "gene_id")), drop = FALSE]
  write.csv(anno_df, annotation_file, row.names = FALSE, fileEncoding = "UTF-8")
}

# ---------- featureCounts stat ----------
stat_file <- file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_stat.csv")
sample_assignment_file <- file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_sample_assignment_summary.csv")

assignment_summary <- data.frame(
  sample_id = sample_ids,
  group = groups,
  assigned = NA_real_,
  total_status_reads = NA_real_,
  assignment_rate_pct = NA_real_,
  stringsAsFactors = FALSE
)

if (!is.null(fc$stat)) {
  stat_df <- as.data.frame(fc$stat, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(stat_df) >= length(sample_ids) + 1) {
    names(stat_df)[1] <- "Status"
    names(stat_df)[2:(length(sample_ids) + 1)] <- sample_ids
  }
  write.csv(stat_df, stat_file, row.names = FALSE, fileEncoding = "UTF-8")

  if ("Status" %in% names(stat_df) && all(sample_ids %in% names(stat_df))) {
    assigned_row <- which(stat_df$Status == "Assigned")
    if (length(assigned_row) == 1) {
      assigned <- as.numeric(stat_df[assigned_row, sample_ids])
      totals <- colSums(sapply(stat_df[, sample_ids, drop = FALSE], as.numeric), na.rm = TRUE)
      assignment_summary$assigned <- assigned
      assignment_summary$total_status_reads <- totals
      assignment_summary$assignment_rate_pct <- round(100 * assigned / totals, 3)
    }
  }
} else {
  stat_df <- NULL
}

write.csv(assignment_summary, sample_assignment_file, row.names = FALSE, fileEncoding = "UTF-8")

# ---------- Compact sample manifest for downstream DE ----------
sample_info <- data.frame(
  sample_id = sample_ids,
  group = groups,
  bam_file = bam_files,
  stringsAsFactors = FALSE
)

sample_info_file <- file.path(tables_dir, "stepH2w5_v2_trimmed_sample_info_for_DE.csv")
write.csv(sample_info, sample_info_file, row.names = FALSE, fileEncoding = "UTF-8")

# ---------- Overall summary ----------
batch_ok_h2w5 <- file.exists(counts_file) &&
  nrow(count_df) > 0 &&
  ncol(count_df) == length(sample_ids) + 1 &&
  all(!is.na(assignment_summary$assigned))

assignment_rates <- assignment_summary$assignment_rate_pct
arg_value <- function(nm) {
  if (nm %in% names(featurecounts_args)) {
    return(paste(featurecounts_args[[nm]], collapse = "; "))
  }
  "not_available_in_this_Rsubread_version"
}
overall_summary <- data.frame(
  metric = c(
    "project_root",
    "branch_dir",
    "h2w4_status_file",
    "bam_check_file",
    "gtf_for_featureCounts",
    "threads",
    "strandSpecific",
    "isPairedEnd",
    "useMetaFeatures",
    "countReadPairs",
    "requireBothEndsMapped",
    "checkFragLength",
    "n_samples",
    "n_genes_in_count_matrix",
    "n_bam_found",
    "min_assignment_rate_pct",
    "mean_assignment_rate_pct",
    "max_assignment_rate_pct",
    "batch_ok_h2w5",
    "counts_file",
    "sample_assignment_file",
    "featureCounts_stat_file",
    "sample_info_file",
    "featureCounts_args_file",
    "fc_rds",
    "tmp_dir_used",
    "log_file",
    "status"
  ),
  value = c(
    project_root,
    branch_dir,
    h2w4_status_file,
    bam_check_file,
    gtf_for_featureCounts,
    as.character(threads),
    arg_value("strandSpecific"),
    arg_value("isPairedEnd"),
    arg_value("useMetaFeatures"),
    arg_value("countReadPairs"),
    arg_value("requireBothEndsMapped"),
    arg_value("checkFragLength"),
    as.character(length(sample_ids)),
    as.character(nrow(count_df)),
    as.character(sum(status_df$bam_exists_h2w5)),
    as.character(round(min(assignment_rates, na.rm = TRUE), 3)),
    as.character(round(mean(assignment_rates, na.rm = TRUE), 3)),
    as.character(round(max(assignment_rates, na.rm = TRUE), 3)),
    as.character(batch_ok_h2w5),
    counts_file,
    sample_assignment_file,
    stat_file,
    sample_info_file,
    featurecounts_args_file,
    fc_rds,
    tmp_dir,
    log_file,
    if (batch_ok_h2w5) "Step H2w5 v2 completed successfully" else "Step H2w5 v2 completed with warnings"
  ),
  stringsAsFactors = FALSE
)

summary_file <- file.path(tables_dir, "stepH2w5_v2_trimmed_featureCounts_overall_summary.csv")
write.csv(overall_summary, summary_file, row.names = FALSE, fileEncoding = "UTF-8")

config_rds <- file.path(objects_dir, "stepH2w5_v2_trimmed_featureCounts_config.rds")
saveRDS(
  list(
    project_root = project_root,
    branch_dir = branch_dir,
    h2w4_status_file = h2w4_status_file,
    bam_files = bam_files,
    sample_ids = sample_ids,
    groups = groups,
    gtf_for_featureCounts = gtf_for_featureCounts,
    featurecounts_args = featurecounts_args,
    tmp_dir = tmp_dir,
    output_files = list(
      counts_file = counts_file,
      sample_assignment_file = sample_assignment_file,
      stat_file = stat_file,
      sample_info_file = sample_info_file,
      summary_file = summary_file,
      fc_rds = fc_rds
    )
  ),
  config_rds
)

# ---------- Print ----------
cat("\n===== Step H2w5 v2 overall summary =====\n")
print(overall_summary, row.names = FALSE)

cat("\n===== Step H2w5 v2 sample assignment summary =====\n")
print(assignment_summary, row.names = FALSE)

cat("\n输出文件：\n")
cat(counts_file, "\n")
cat(sample_assignment_file, "\n")
cat(stat_file, "\n")
cat(sample_info_file, "\n")
cat(summary_file, "\n")
cat(config_rds, "\n")
cat(log_file, "\n")

if (isTRUE(batch_ok_h2w5)) {
  cat("\n结果：H2w5 v2 成功。trimmed BAM 已完成 featureCounts，并生成 gene-level count matrix。下一步可以进入 trimmed branch 的 edgeR QLF DE。\n")
} else {
  cat("\n结果：H2w5 v2 完成但有警告。请优先查看 summary 和 assignment summary。\n")
}

cat("\n===== End of Step H2w5 v2 =====\n")
