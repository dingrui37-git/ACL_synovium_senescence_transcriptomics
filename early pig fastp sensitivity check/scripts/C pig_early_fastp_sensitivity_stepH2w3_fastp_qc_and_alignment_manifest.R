# ======================================================================
# Step H2w3: Summarize fastp trimming outputs and prepare trimmed FASTQ manifest
# Project: ACLsenescence2 / pig early fastp sensitivity branch
#
# 功能：
# 1) 检查 Step H2w2 v3 生成的 18 个样本 trimmed FASTQ / html / json 是否完整
# 2) 解析 fastp JSON，提取 before/after reads、bases、Q30、GC、adapter trimming 等指标
# 3) 输出 trimmed FASTQ 后续比对用 manifest
#
# 运行方式：
# source("E:/R/ACLsenescence2/pig_early_fastp_sensitivity_stepH2w3_fastp_qc_and_alignment_manifest.R")
# ======================================================================

options(stringsAsFactors = FALSE)

cat("\n===== Step H2w3: fastp QC summary + trimmed alignment manifest =====\n")

# -----------------------------
# 0. 固定项目路径
# -----------------------------
project_root <- "E:/R/ACLsenescence2"

branch_dir <- file.path(
  project_root,
  "rebuild_submission",
  "02_pig_early_fastp_sensitivity"
)

tables_dir <- file.path(branch_dir, "tables")
objects_dir <- file.path(branch_dir, "objects")
logs_dir <- file.path(branch_dir, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

command_sheet_file <- file.path(
  tables_dir,
  "stepH1_fastp_command_sheet_18samples.csv"
)

h2w2_status_file <- file.path(
  tables_dir,
  "stepH2w2_v3_fastp_batch_status_18samples.csv"
)

qc_summary_file <- file.path(
  tables_dir,
  "stepH2w3_fastp_trimmed_output_qc_summary_18samples.csv"
)

alignment_manifest_file <- file.path(
  tables_dir,
  "stepH2w3_trimmed_fastq_manifest_for_alignment_18samples.csv"
)

overall_summary_file <- file.path(
  tables_dir,
  "stepH2w3_fastp_qc_overall_summary.csv"
)

config_rds <- file.path(
  objects_dir,
  "stepH2w3_fastp_qc_and_alignment_manifest_config.rds"
)

log_file <- file.path(
  logs_dir,
  "stepH2w3_fastp_qc_and_alignment_manifest_log.txt"
)

sink(log_file, split = TRUE)
on.exit({
  try(sink(), silent = TRUE)
}, add = TRUE)

cat("project_root:", project_root, "\n")
cat("branch_dir:", branch_dir, "\n")
cat("command_sheet_file:", command_sheet_file, "\n")
cat("h2w2_status_file:", h2w2_status_file, "\n\n")

# -----------------------------
# 1. 依赖检查
# -----------------------------
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop(
    paste0(
      "缺少 R 包 jsonlite，无法解析 fastp JSON。\n",
      "请先在 R 里运行：install.packages('jsonlite')\n",
      "然后重新运行本脚本。"
    )
  )
}

if (!file.exists(command_sheet_file)) {
  stop("没找到 Step H1 command sheet：", command_sheet_file)
}

# -----------------------------
# 2. 辅助函数
# -----------------------------
read_csv_safely <- function(path) {
  x <- tryCatch(
    read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e) {
      read.csv(path, check.names = FALSE)
    }
  )
  x
}

file_mb <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NA_real_)
  }
  round(as.numeric(file.info(path)$size) / 1024^2, 3)
}

safe_num <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x[[1]]))
}

get_nested <- function(x, path, default = NA) {
  z <- x
  for (nm in path) {
    if (is.null(z)) {
      return(default)
    }
    if (is.list(z)) {
      if (!nm %in% names(z)) {
        return(default)
      }
      z <- z[[nm]]
    } else {
      return(default)
    }
  }
  if (is.null(z) || length(z) == 0) {
    return(default)
  }
  z
}

rate_pct <- function(numer, denom) {
  numer <- suppressWarnings(as.numeric(numer))
  denom <- suppressWarnings(as.numeric(denom))
  if (is.na(numer) || is.na(denom) || denom == 0) {
    return(NA_real_)
  }
  round(100 * numer / denom, 3)
}

round3 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA_real_, round(x, 3))
}

round6 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA_real_, round(x, 6))
}

# -----------------------------
# 3. 读取 H1 command sheet 和 H2w2 状态表
# -----------------------------
cmd <- read_csv_safely(command_sheet_file)

required_cols <- c(
  "sample_id",
  "group",
  "trimmed_R1",
  "trimmed_R2",
  "fastp_html",
  "fastp_json"
)

missing_cols <- setdiff(required_cols, colnames(cmd))
if (length(missing_cols) > 0) {
  stop(
    "Step H1 command sheet 缺少必要列：",
    paste(missing_cols, collapse = ", ")
  )
}

h2w2_status <- NULL
if (file.exists(h2w2_status_file)) {
  h2w2_status <- read_csv_safely(h2w2_status_file)
  if (!"sample_id" %in% colnames(h2w2_status)) {
    h2w2_status <- NULL
  }
}

cat("n_samples_in_command_sheet:", nrow(cmd), "\n")

# -----------------------------
# 4. 逐样本解析 fastp JSON
# -----------------------------
qc_rows <- vector("list", nrow(cmd))

for (i in seq_len(nrow(cmd))) {
  sample_id <- cmd$sample_id[i]
  group <- cmd$group[i]

  raw_R1 <- if ("absolute_path.R1" %in% colnames(cmd)) cmd[["absolute_path.R1"]][i] else NA_character_
  raw_R2 <- if ("absolute_path.R2" %in% colnames(cmd)) cmd[["absolute_path.R2"]][i] else NA_character_

  trimmed_R1 <- cmd$trimmed_R1[i]
  trimmed_R2 <- cmd$trimmed_R2[i]
  fastp_html <- cmd$fastp_html[i]
  fastp_json <- cmd$fastp_json[i]

  raw_R1_exists <- !is.na(raw_R1) && nzchar(raw_R1) && file.exists(raw_R1)
  raw_R2_exists <- !is.na(raw_R2) && nzchar(raw_R2) && file.exists(raw_R2)
  trimmed_R1_exists <- !is.na(trimmed_R1) && nzchar(trimmed_R1) && file.exists(trimmed_R1)
  trimmed_R2_exists <- !is.na(trimmed_R2) && nzchar(trimmed_R2) && file.exists(trimmed_R2)
  fastp_html_exists <- !is.na(fastp_html) && nzchar(fastp_html) && file.exists(fastp_html)
  fastp_json_exists <- !is.na(fastp_json) && nzchar(fastp_json) && file.exists(fastp_json)

  parsed_ok <- FALSE
  parse_error <- NA_character_

  before_reads <- NA_real_
  after_reads <- NA_real_
  before_bases <- NA_real_
  after_bases <- NA_real_
  before_q30_rate <- NA_real_
  after_q30_rate <- NA_real_
  before_gc <- NA_real_
  after_gc <- NA_real_
  read1_mean_length_before <- NA_real_
  read2_mean_length_before <- NA_real_
  read1_mean_length_after <- NA_real_
  read2_mean_length_after <- NA_real_
  passed_filter_reads <- NA_real_
  low_quality_reads <- NA_real_
  too_many_N_reads <- NA_real_
  too_short_reads <- NA_real_
  adapter_trimmed_reads <- NA_real_
  adapter_trimmed_bases <- NA_real_
  duplication_rate <- NA_real_

  if (fastp_json_exists) {
    j <- tryCatch(
      jsonlite::fromJSON(fastp_json, simplifyVector = FALSE),
      error = function(e) e
    )

    if (inherits(j, "error")) {
      parse_error <- conditionMessage(j)
    } else {
      parsed_ok <- TRUE

      before_reads <- safe_num(get_nested(j, c("summary", "before_filtering", "total_reads")))
      after_reads <- safe_num(get_nested(j, c("summary", "after_filtering", "total_reads")))

      before_bases <- safe_num(get_nested(j, c("summary", "before_filtering", "total_bases")))
      after_bases <- safe_num(get_nested(j, c("summary", "after_filtering", "total_bases")))

      before_q30_rate <- safe_num(get_nested(j, c("summary", "before_filtering", "q30_rate")))
      after_q30_rate <- safe_num(get_nested(j, c("summary", "after_filtering", "q30_rate")))

      before_gc <- safe_num(get_nested(j, c("summary", "before_filtering", "gc_content")))
      after_gc <- safe_num(get_nested(j, c("summary", "after_filtering", "gc_content")))

      read1_mean_length_before <- safe_num(get_nested(j, c("summary", "before_filtering", "read1_mean_length")))
      read2_mean_length_before <- safe_num(get_nested(j, c("summary", "before_filtering", "read2_mean_length")))
      read1_mean_length_after <- safe_num(get_nested(j, c("summary", "after_filtering", "read1_mean_length")))
      read2_mean_length_after <- safe_num(get_nested(j, c("summary", "after_filtering", "read2_mean_length")))

      passed_filter_reads <- safe_num(get_nested(j, c("filtering_result", "passed_filter_reads")))
      low_quality_reads <- safe_num(get_nested(j, c("filtering_result", "low_quality_reads")))
      too_many_N_reads <- safe_num(get_nested(j, c("filtering_result", "too_many_N_reads")))
      too_short_reads <- safe_num(get_nested(j, c("filtering_result", "too_short_reads")))

      adapter_trimmed_reads <- safe_num(get_nested(j, c("adapter_cutting", "adapter_trimmed_reads")))
      adapter_trimmed_bases <- safe_num(get_nested(j, c("adapter_cutting", "adapter_trimmed_bases")))

      duplication_rate <- safe_num(get_nested(j, c("duplication", "rate")))
    }
  }

  h2w2_exit_status <- NA
  h2w2_sample_ok <- NA

  if (!is.null(h2w2_status)) {
    m <- match(sample_id, h2w2_status$sample_id)
    if (!is.na(m)) {
      if ("exit_status" %in% colnames(h2w2_status)) {
        h2w2_exit_status <- h2w2_status$exit_status[m]
      }
      if ("sample_ok" %in% colnames(h2w2_status)) {
        h2w2_sample_ok <- h2w2_status$sample_ok[m]
      }
    }
  }

  retained_read_pct <- rate_pct(after_reads, before_reads)
  retained_base_pct <- rate_pct(after_bases, before_bases)
  adapter_trimmed_read_pct <- rate_pct(adapter_trimmed_reads, before_reads)
  low_quality_read_pct <- rate_pct(low_quality_reads, before_reads)
  too_many_N_read_pct <- rate_pct(too_many_N_reads, before_reads)
  too_short_read_pct <- rate_pct(too_short_reads, before_reads)

  sample_ok_h2w3 <- all(c(
    trimmed_R1_exists,
    trimmed_R2_exists,
    fastp_html_exists,
    fastp_json_exists,
    parsed_ok,
    !is.na(after_reads),
    after_reads > 0
  ))

  qc_rows[[i]] <- data.frame(
    sample_id = sample_id,
    group = group,

    raw_R1_exists = raw_R1_exists,
    raw_R2_exists = raw_R2_exists,
    trimmed_R1_exists = trimmed_R1_exists,
    trimmed_R2_exists = trimmed_R2_exists,
    fastp_html_exists = fastp_html_exists,
    fastp_json_exists = fastp_json_exists,

    raw_R1_size_mb = file_mb(raw_R1),
    raw_R2_size_mb = file_mb(raw_R2),
    trimmed_R1_size_mb = file_mb(trimmed_R1),
    trimmed_R2_size_mb = file_mb(trimmed_R2),

    h2w2_exit_status = h2w2_exit_status,
    h2w2_sample_ok = h2w2_sample_ok,

    fastp_json_parsed_ok = parsed_ok,
    fastp_json_parse_error = parse_error,

    before_total_reads = before_reads,
    after_total_reads = after_reads,
    retained_read_pct = retained_read_pct,

    before_total_bases = before_bases,
    after_total_bases = after_bases,
    retained_base_pct = retained_base_pct,

    before_q30_rate = round6(before_q30_rate),
    after_q30_rate = round6(after_q30_rate),
    before_q30_pct = round3(before_q30_rate * 100),
    after_q30_pct = round3(after_q30_rate * 100),

    before_gc_content = round6(before_gc),
    after_gc_content = round6(after_gc),
    before_gc_pct = round3(before_gc * 100),
    after_gc_pct = round3(after_gc * 100),

    read1_mean_length_before = round3(read1_mean_length_before),
    read2_mean_length_before = round3(read2_mean_length_before),
    read1_mean_length_after = round3(read1_mean_length_after),
    read2_mean_length_after = round3(read2_mean_length_after),

    passed_filter_reads = passed_filter_reads,
    low_quality_reads = low_quality_reads,
    too_many_N_reads = too_many_N_reads,
    too_short_reads = too_short_reads,

    low_quality_read_pct = low_quality_read_pct,
    too_many_N_read_pct = too_many_N_read_pct,
    too_short_read_pct = too_short_read_pct,

    adapter_trimmed_reads = adapter_trimmed_reads,
    adapter_trimmed_bases = adapter_trimmed_bases,
    adapter_trimmed_read_pct = adapter_trimmed_read_pct,

    duplication_rate = round6(duplication_rate),
    duplication_pct = round3(duplication_rate * 100),

    sample_ok_h2w3 = sample_ok_h2w3,

    raw_R1 = raw_R1,
    raw_R2 = raw_R2,
    trimmed_R1 = trimmed_R1,
    trimmed_R2 = trimmed_R2,
    fastp_html = fastp_html,
    fastp_json = fastp_json,

    stringsAsFactors = FALSE
  )

  cat(
    "checked", i, "/", nrow(cmd), ": ",
    sample_id,
    " | ok=", sample_ok_h2w3,
    " | retained_reads=", retained_read_pct, "%",
    " | after_Q30=", round3(after_q30_rate * 100), "%\n",
    sep = ""
  )
}

qc_df <- do.call(rbind, qc_rows)

# -----------------------------
# 5. 输出 trimmed alignment manifest
# -----------------------------
alignment_manifest <- qc_df[, c(
  "sample_id",
  "group",
  "trimmed_R1",
  "trimmed_R2",
  "fastp_json",
  "fastp_html",
  "before_total_reads",
  "after_total_reads",
  "retained_read_pct",
  "after_q30_pct",
  "sample_ok_h2w3"
)]

colnames(alignment_manifest)[colnames(alignment_manifest) == "trimmed_R1"] <- "R1_fastq"
colnames(alignment_manifest)[colnames(alignment_manifest) == "trimmed_R2"] <- "R2_fastq"

# -----------------------------
# 6. 总体 summary
# -----------------------------
batch_ok_h2w3 <- all(qc_df$sample_ok_h2w3)

overall_summary <- data.frame(
  metric = c(
    "project_root",
    "branch_dir",
    "command_sheet_file",
    "h2w2_status_file",
    "n_samples",
    "n_trimmed_R1_exists",
    "n_trimmed_R2_exists",
    "n_fastp_html_exists",
    "n_fastp_json_exists",
    "n_fastp_json_parsed_ok",
    "n_sample_ok_h2w3",
    "n_sample_failed_or_incomplete_h2w3",
    "min_retained_read_pct",
    "mean_retained_read_pct",
    "max_retained_read_pct",
    "min_after_q30_pct",
    "mean_after_q30_pct",
    "max_after_q30_pct",
    "mean_adapter_trimmed_read_pct",
    "batch_ok_h2w3",
    "qc_summary_file",
    "alignment_manifest_file",
    "overall_summary_file",
    "config_rds",
    "log_file",
    "status"
  ),
  value = as.character(c(
    project_root,
    branch_dir,
    command_sheet_file,
    h2w2_status_file,
    nrow(qc_df),
    sum(qc_df$trimmed_R1_exists),
    sum(qc_df$trimmed_R2_exists),
    sum(qc_df$fastp_html_exists),
    sum(qc_df$fastp_json_exists),
    sum(qc_df$fastp_json_parsed_ok),
    sum(qc_df$sample_ok_h2w3),
    sum(!qc_df$sample_ok_h2w3),
    round(min(qc_df$retained_read_pct, na.rm = TRUE), 3),
    round(mean(qc_df$retained_read_pct, na.rm = TRUE), 3),
    round(max(qc_df$retained_read_pct, na.rm = TRUE), 3),
    round(min(qc_df$after_q30_pct, na.rm = TRUE), 3),
    round(mean(qc_df$after_q30_pct, na.rm = TRUE), 3),
    round(max(qc_df$after_q30_pct, na.rm = TRUE), 3),
    round(mean(qc_df$adapter_trimmed_read_pct, na.rm = TRUE), 3),
    batch_ok_h2w3,
    qc_summary_file,
    alignment_manifest_file,
    overall_summary_file,
    config_rds,
    log_file,
    if (batch_ok_h2w3) {
      "Step H2w3 completed successfully"
    } else {
      "Step H2w3 completed with failed or incomplete samples"
    }
  )),
  stringsAsFactors = FALSE
)

config <- list(
  project_root = project_root,
  branch_dir = branch_dir,
  command_sheet_file = command_sheet_file,
  h2w2_status_file = h2w2_status_file,
  qc_summary_file = qc_summary_file,
  alignment_manifest_file = alignment_manifest_file,
  overall_summary_file = overall_summary_file,
  config_rds = config_rds,
  log_file = log_file,
  created_at = Sys.time()
)

# -----------------------------
# 7. 保存文件
# -----------------------------
write.csv(qc_df, qc_summary_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(alignment_manifest, alignment_manifest_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(overall_summary, overall_summary_file, row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(config, config_rds)

cat("\n===== Step H2w3 overall summary =====\n")
print(overall_summary, row.names = FALSE)

cat("\n===== Step H2w3 per-sample compact QC =====\n")
print(
  qc_df[, c(
    "sample_id",
    "group",
    "sample_ok_h2w3",
    "before_total_reads",
    "after_total_reads",
    "retained_read_pct",
    "after_q30_pct",
    "adapter_trimmed_read_pct"
  )],
  row.names = FALSE
)

cat("\n输出文件：\n")
cat(qc_summary_file, "\n")
cat(alignment_manifest_file, "\n")
cat(overall_summary_file, "\n")
cat(config_rds, "\n")
cat(log_file, "\n")

if (isTRUE(batch_ok_h2w3)) {
  cat("\n结果：H2w3 成功。trimmed FASTQ 输出完整，fastp JSON 可解析，后续可以进入 trimmed FASTQ 比对步骤。\n")
} else {
  cat("\n结果：H2w3 发现部分样本不完整。请优先查看 stepH2w3_fastp_trimmed_output_qc_summary_18samples.csv 中 sample_ok_h2w3=FALSE 的行。\n")
}

cat("\n===== End of Step H2w3 =====\n")
