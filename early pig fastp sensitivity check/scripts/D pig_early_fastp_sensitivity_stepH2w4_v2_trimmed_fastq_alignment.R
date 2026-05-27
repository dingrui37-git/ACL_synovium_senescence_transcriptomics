# ============================================================
# Step H2w4 v2
# pig early fastp sensitivity:
# 使用 fastp trimmed FASTQ 重新进行 Rsubread 比对
#
# 作用：
# 1. 读取 Step H2w3 生成的 trimmed FASTQ alignment manifest
# 2. 自动寻找项目中的 pig Rsubread genome index
# 3. 对 18 个样本逐个运行 Rsubread::align
# 4. 输出 BAM、每样本状态表、总 summary 和日志
#
# 注意：
# - 这一步只做 trimmed FASTQ 的比对
# - 不做 featureCounts
# - 不做 DE / GSEA
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n===== Step H2w4 v2: trimmed FASTQ alignment with Rsubread =====\n\n")

# -----------------------------
# 0. 固定项目路径
# -----------------------------
project_root <- "E:/R/ACLsenescence2"
branch_dir <- file.path(project_root, "rebuild_submission", "02_pig_early_fastp_sensitivity")

manifest_file <- file.path(
  branch_dir,
  "tables",
  "stepH2w3_trimmed_fastq_manifest_for_alignment_18samples.csv"
)

tables_dir  <- file.path(branch_dir, "tables")
logs_dir    <- file.path(branch_dir, "logs")
objects_dir <- file.path(branch_dir, "objects")
bam_dir     <- file.path(branch_dir, "bam_trimmed")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)

status_file <- file.path(tables_dir, "stepH2w4_v2_trimmed_alignment_status_18samples.csv")
summary_file <- file.path(tables_dir, "stepH2w4_v2_trimmed_alignment_summary.csv")
index_candidates_file <- file.path(tables_dir, "stepH2w4_v2_rsubread_index_candidates.csv")
main_log_file <- file.path(logs_dir, "stepH2w4_v2_trimmed_alignment_main_log.txt")
config_rds <- file.path(objects_dir, "stepH2w4_v2_trimmed_alignment_config.rds")

# 比对线程数：可以根据电脑情况调整；为了稳妥先用 4
threads <- 4

# 如果重复运行脚本，已经完整存在的 BAM 默认跳过
skip_existing_complete <- TRUE

# BAM 文件至少大于这个大小才视为非空有效；这里只是防止 0 字节假文件
min_bam_size_bytes <- 1024

# -----------------------------
# 1. 基础检查
# -----------------------------
if (!dir.exists(project_root)) {
  stop("找不到项目根目录：", project_root)
}

if (!file.exists(manifest_file)) {
  stop("找不到 H2w3 alignment manifest：", manifest_file)
}

if (!requireNamespace("Rsubread", quietly = TRUE)) {
  stop(
    "当前 R 环境没有安装 Rsubread。\n",
    "请先确认原 pig early 主分析使用的 R 环境是否已经安装 Rsubread。"
  )
}

# -----------------------------
# 2. 读取 trimmed FASTQ manifest
# -----------------------------
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)

# H2w3 manifest 在当前实际输出中使用的是 R1_fastq / R2_fastq；
# 早先 H2w4 脚本预期 trimmed_R1 / trimmed_R2，导致还没开始比对就报错。
# v2 在这里同时兼容两种列名。
base_cols <- c("sample_id", "group")
missing_base_cols <- setdiff(base_cols, colnames(manifest))
if (length(missing_base_cols) > 0) {
  stop(
    "H2w3 manifest 缺少必要基础列：",
    paste(missing_base_cols, collapse = ", "),
    "\n当前列名为：",
    paste(colnames(manifest), collapse = ", ")
  )
}

if (all(c("trimmed_R1", "trimmed_R2") %in% colnames(manifest))) {
  r1_col <- "trimmed_R1"
  r2_col <- "trimmed_R2"
} else if (all(c("R1_fastq", "R2_fastq") %in% colnames(manifest))) {
  r1_col <- "R1_fastq"
  r2_col <- "R2_fastq"
} else {
  stop(
    "H2w3 manifest 中没有找到 trimmed FASTQ 路径列。\n",
    "需要 trimmed_R1/trimmed_R2 或 R1_fastq/R2_fastq。\n",
    "当前列名为：",
    paste(colnames(manifest), collapse = ", ")
  )
}

if ("sample_ok_h2w3" %in% colnames(manifest)) {
  bad_h2w3 <- manifest[!isTRUE(all(manifest$sample_ok_h2w3 %in% c(TRUE, "TRUE", "true", 1, "1"))), , drop = FALSE]
  if (nrow(bad_h2w3) > 0) {
    write.csv(bad_h2w3, file.path(tables_dir, "stepH2w4_v2_manifest_rows_not_ok_h2w3.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    stop(
      "H2w3 manifest 中存在 sample_ok_h2w3 不是 TRUE 的行。请查看：",
      file.path(tables_dir, "stepH2w4_v2_manifest_rows_not_ok_h2w3.csv")
    )
  }
}

manifest <- data.frame(
  sample_id = as.character(manifest$sample_id),
  group = as.character(manifest$group),
  trimmed_R1 = as.character(manifest[[r1_col]]),
  trimmed_R2 = as.character(manifest[[r2_col]]),
  manifest_R1_column_used = r1_col,
  manifest_R2_column_used = r2_col,
  stringsAsFactors = FALSE
)

cat("H2w4 v2 读取 trimmed FASTQ 路径列：", r1_col, ", ", r2_col, "\n", sep = "")

# 检查 trimmed FASTQ 是否存在
manifest$trimmed_R1_exists <- file.exists(manifest$trimmed_R1)
manifest$trimmed_R2_exists <- file.exists(manifest$trimmed_R2)

if (any(!manifest$trimmed_R1_exists | !manifest$trimmed_R2_exists)) {
  bad <- manifest[!manifest$trimmed_R1_exists | !manifest$trimmed_R2_exists, ]
  write.csv(bad, file.path(tables_dir, "stepH2w4_missing_trimmed_fastq_files.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  stop(
    "发现 trimmed FASTQ 文件缺失。请查看：",
    file.path(tables_dir, "stepH2w4_missing_trimmed_fastq_files.csv")
  )
}

# -----------------------------
# 3. 自动寻找 Rsubread genome index
# -----------------------------
cat("正在项目目录中寻找 Rsubread index...\n")

index_array_files <- list.files(
  project_root,
  pattern = "\\.00\\.b\\.array$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(index_array_files) == 0) {
  stop(
    "没有在项目目录中找到 Rsubread index 文件（*.00.b.array）。\n",
    "请确认原 pig early 比对使用的 Rsubread index 是否仍在 E:/R/ACLsenescence2 下面。"
  )
}

index_bases <- sub("\\.00\\.b\\.array$", "", index_array_files, ignore.case = TRUE)

index_candidates <- data.frame(
  index_base = normalizePath(index_bases, winslash = "/", mustWork = FALSE),
  array_file = normalizePath(index_array_files, winslash = "/", mustWork = FALSE),
  stringsAsFactors = FALSE
)

# 给候选 index 打分：优先选择 Sus scrofa / Sscrofa / pig 相关 index
cand_lower <- tolower(index_candidates$index_base)
index_candidates$score <- 0
index_candidates$score <- index_candidates$score + ifelse(grepl("ssc|sscrofa|sus_scrofa|sus-scrofa|pig|porcine", cand_lower), 10, 0)
index_candidates$score <- index_candidates$score + ifelse(grepl("mouse|mus_musculus|mm10|grcm", cand_lower), -20, 0)
index_candidates$score <- index_candidates$score + ifelse(grepl("02_pig_early|pig_early|reference|genome|index", cand_lower), 2, 0)
index_candidates$score <- index_candidates$score + ifelse(grepl("fastp_sensitivity", cand_lower), -1, 0)

index_candidates <- index_candidates[order(-index_candidates$score, index_candidates$index_base), ]
write.csv(index_candidates, index_candidates_file, row.names = FALSE, fileEncoding = "UTF-8")

chosen_index <- index_candidates$index_base[1]
chosen_score <- index_candidates$score[1]

if (chosen_score <= 0 && nrow(index_candidates) > 1) {
  cat("\n警告：找到多个 index，但没有明显 pig/Sscrofa 标记。将使用候选列表第一项：\n")
  cat(chosen_index, "\n")
  cat("完整候选表已输出：", index_candidates_file, "\n\n")
}

cat("使用 Rsubread index:\n")
cat(chosen_index, "\n\n")

# -----------------------------
# 4. 逐样本比对
# -----------------------------
status_rows <- list()

cat("开始 trimmed FASTQ 比对，共 ", nrow(manifest), " 个样本。\n", sep = "")
cat("主日志文件：", main_log_file, "\n\n")

log_con <- file(main_log_file, open = "wt", encoding = "UTF-8")
writeLines(c(
  "===== Step H2w4 v2 trimmed alignment main log =====",
  paste("project_root:", project_root),
  paste("branch_dir:", branch_dir),
  paste("manifest_file:", manifest_file),
  paste("chosen_index:", chosen_index),
  paste("threads:", threads),
  paste("start_time:", as.character(Sys.time())),
  ""
), con = log_con)

for (i in seq_len(nrow(manifest))) {
  sid <- manifest$sample_id[i]
  grp <- manifest$group[i]
  r1 <- manifest$trimmed_R1[i]
  r2 <- manifest$trimmed_R2[i]

  bam_file <- file.path(bam_dir, paste0(sid, ".trimmed.Rsubread.bam"))
  sample_log_file <- file.path(logs_dir, paste0("stepH2w4_v2_align_", sid, ".log"))
  expected_summary_file <- paste0(bam_file, ".summary")

  bam_exists_before <- file.exists(bam_file)
  summary_exists_before <- file.exists(expected_summary_file)
  bam_size_before <- if (bam_exists_before) file.info(bam_file)$size else NA_real_

  already_complete <- isTRUE(skip_existing_complete) &&
    file.exists(bam_file) &&
    !is.na(file.info(bam_file)$size) &&
    file.info(bam_file)$size > min_bam_size_bytes &&
    file.exists(expected_summary_file)

  cat("开始样本 ", i, "/", nrow(manifest), ": ", sid, "\n", sep = "")
  writeLines(paste0("\n----- Sample ", i, "/", nrow(manifest), ": ", sid, " -----"), con = log_con)
  writeLines(paste("group:", grp), con = log_con)
  writeLines(paste("R1:", r1), con = log_con)
  writeLines(paste("R2:", r2), con = log_con)
  writeLines(paste("BAM:", bam_file), con = log_con)

  action <- if (already_complete) "skip_existing_complete" else "run_align"
  exit_status <- NA_integer_
  error_message <- NA_character_

  if (already_complete) {
    cat("样本 ", sid, " 已有完整 BAM，跳过。\n", sep = "")
    writeLines("Action: skip_existing_complete", con = log_con)
    exit_status <- 0L
  } else {
    # 如果存在旧的不完整 BAM / summary，先删除，避免混淆
    if (file.exists(bam_file)) unlink(bam_file)
    if (file.exists(expected_summary_file)) unlink(expected_summary_file)

    align_output <- NULL

    result <- tryCatch({
      align_output <- capture.output({
        Rsubread::align(
          index = chosen_index,
          readfile1 = r1,
          readfile2 = r2,
          type = "rna",
          input_format = "gzFASTQ",
          output_format = "BAM",
          output_file = bam_file,
          nthreads = threads,
          phredOffset = 33
        )
      }, type = "output")

      writeLines(align_output, con = sample_log_file, useBytes = TRUE)
      writeLines(align_output, con = log_con, useBytes = TRUE)
      TRUE
    }, error = function(e) {
      error_message <<- conditionMessage(e)
      writeLines(paste("ERROR:", error_message), con = sample_log_file, useBytes = TRUE)
      writeLines(paste("ERROR:", error_message), con = log_con, useBytes = TRUE)
      FALSE
    })

    exit_status <- if (isTRUE(result)) 0L else 1L
    cat("样本 ", sid, " align exit status: ", exit_status, "\n", sep = "")
  }

  bam_exists_after <- file.exists(bam_file)
  summary_exists_after <- file.exists(expected_summary_file)
  bam_size_after <- if (bam_exists_after) file.info(bam_file)$size else NA_real_

  sample_ok <- exit_status == 0L &&
    isTRUE(bam_exists_after) &&
    !is.na(bam_size_after) &&
    bam_size_after > min_bam_size_bytes &&
    isTRUE(summary_exists_after)

  status_rows[[i]] <- data.frame(
    sample_id = sid,
    group = grp,
    action = action,
    exit_status = exit_status,
    sample_ok_h2w4 = sample_ok,
    trimmed_R1 = r1,
    trimmed_R2 = r2,
    bam_file = bam_file,
    bam_exists_before = bam_exists_before,
    bam_size_before = bam_size_before,
    bam_exists_after = bam_exists_after,
    bam_size_after = bam_size_after,
    rsubread_summary_file = expected_summary_file,
    rsubread_summary_exists = summary_exists_after,
    sample_log_file = sample_log_file,
    error_message = error_message,
    stringsAsFactors = FALSE
  )

  status_df_tmp <- do.call(rbind, status_rows[seq_len(i)])
  write.csv(status_df_tmp, status_file, row.names = FALSE, fileEncoding = "UTF-8")

  writeLines(paste("sample_ok_h2w4:", sample_ok), con = log_con)
  flush(log_con)
}

writeLines(c("", paste("end_time:", as.character(Sys.time()))), con = log_con)
close(log_con)

status_df <- do.call(rbind, status_rows)
write.csv(status_df, status_file, row.names = FALSE, fileEncoding = "UTF-8")

# -----------------------------
# 5. 总结
# -----------------------------
n_samples <- nrow(status_df)
n_samples_ok <- sum(status_df$sample_ok_h2w4, na.rm = TRUE)
n_samples_failed_or_incomplete <- n_samples - n_samples_ok
n_run_align <- sum(status_df$action == "run_align", na.rm = TRUE)
n_skipped_existing_complete <- sum(status_df$action == "skip_existing_complete", na.rm = TRUE)
batch_ok_h2w4 <- (n_samples == 18 && n_samples_ok == 18)

summary_df <- data.frame(
  metric = c(
    "project_root",
    "branch_dir",
    "manifest_file",
    "manifest_R1_column_used",
    "manifest_R2_column_used",
    "chosen_rsubread_index",
    "threads",
    "n_samples_in_manifest",
    "n_samples_ok",
    "n_samples_failed_or_incomplete",
    "n_run_align",
    "n_skipped_existing_complete",
    "batch_ok_h2w4",
    "bam_dir",
    "status_file",
    "summary_file",
    "index_candidates_file",
    "main_log_file",
    "config_rds",
    "status"
  ),
  value = c(
    project_root,
    branch_dir,
    manifest_file,
    r1_col,
    r2_col,
    chosen_index,
    as.character(threads),
    as.character(n_samples),
    as.character(n_samples_ok),
    as.character(n_samples_failed_or_incomplete),
    as.character(n_run_align),
    as.character(n_skipped_existing_complete),
    as.character(batch_ok_h2w4),
    bam_dir,
    status_file,
    summary_file,
    index_candidates_file,
    main_log_file,
    config_rds,
    if (batch_ok_h2w4) {
      "Step H2w4 v2 completed successfully"
    } else {
      "Step H2w4 v2 completed with failed or incomplete samples"
    }
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_file, row.names = FALSE, fileEncoding = "UTF-8")

config <- list(
  project_root = project_root,
  branch_dir = branch_dir,
  manifest_file = manifest_file,
  chosen_rsubread_index = chosen_index,
  threads = threads,
  skip_existing_complete = skip_existing_complete,
  min_bam_size_bytes = min_bam_size_bytes,
  bam_dir = bam_dir,
  status_file = status_file,
  summary_file = summary_file,
  index_candidates_file = index_candidates_file,
  main_log_file = main_log_file,
  created_at = Sys.time()
)
saveRDS(config, config_rds)

cat("\n===== Step H2w4 v2 summary =====\n")
print(summary_df, row.names = FALSE)

cat("\n===== Step H2w4 v2 per-sample status =====\n")
print(status_df[, c("sample_id", "group", "action", "exit_status", "sample_ok_h2w4", "bam_size_after")],
      row.names = FALSE)

cat("\n输出文件：\n")
cat(status_file, "\n")
cat(summary_file, "\n")
cat(index_candidates_file, "\n")
cat(main_log_file, "\n")
cat(config_rds, "\n")

if (isTRUE(batch_ok_h2w4)) {
  cat("\n结果：H2w4 v2 成功。18 个 trimmed FASTQ 样本已经完成 Rsubread 比对，下一步可以进入 featureCounts。\n")
} else {
  cat("\n结果：H2w4 v2 有样本失败或 BAM 输出不完整。请优先查看 stepH2w4_v2_trimmed_alignment_status_18samples.csv 中 sample_ok_h2w4=FALSE 的行。\n")
}

cat("\n===== End of Step H2w4 v2 =====\n")
