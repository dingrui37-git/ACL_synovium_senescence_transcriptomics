# =========================================================
# pig_early_fastqc_multiqc_stepE_batch_move_and_multiqc.R
# 目的：
# 1) 对 pig early 18 个样本（36 个 FASTQ）逐个运行 FastQC
# 2) 由于当前 Windows/FastQC 环境下 -o/--outdir 参数不能稳定解析，
#    因此让 FastQC 按默认行为把结果写到原始 FASTQ 所在目录
# 3) 每跑完一个 FASTQ，立刻把 *_fastqc.html / *_fastqc.zip 移动到
#    rebuild_submission/02_pig_early/qc/fastqc_raw
# 4) 全部成功后，用 conda run -n multiqc_env multiqc 汇总生成 MultiQC
#
# 说明：
# - 这是对当前 Windows 本地环境的兼容性修复版，不改变主分析数据本身
# - 如果某个 FASTQ 的 html/zip 未生成，脚本会立即停止
# - 运行结束后，请再运行配套 quick check 脚本，把控制台输出发回 ChatGPT
# =========================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
})

# -----------------------------
# 0. 手工确认路径
# -----------------------------
FASTQC_DIR <- "D:/software/fastq/fastqc_v0.12.1/FastQC"
FASTQC_BAT <- file.path(FASTQC_DIR, "run_fastqc.bat")
CONDA_EXE <- "D:/software/Conda/1/Scripts/conda.exe"
MULTIQC_ENV_NAME <- "multiqc_env"

# -----------------------------
# 1. 基础路径
# -----------------------------
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_root, "data_raw", "E-MTAB-6664")
pig_dir <- file.path(project_root, "rebuild_submission", "02_pig_early")
qc_dir <- file.path(pig_dir, "qc")
fastqc_out_dir <- file.path(qc_dir, "fastqc_raw")
multiqc_out_dir <- file.path(qc_dir, "multiqc_raw")
tables_dir <- file.path(pig_dir, "tables")
logs_dir <- file.path(pig_dir, "logs")
objects_dir <- file.path(pig_dir, "objects")

for (d in c(qc_dir, fastqc_out_dir, multiqc_out_dir, tables_dir, logs_dir, objects_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(logs_dir, "stepE_batch_fastqc_multiqc_log.txt")
if (file.exists(log_file)) file.remove(log_file)

append_log <- function(...) {
  msg <- paste0(...)
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

append_log("Step E started.")
append_log("project_root = ", project_root)
append_log("raw_dir = ", raw_dir)
append_log("FASTQC_BAT exists = ", file.exists(FASTQC_BAT))
append_log("CONDA_EXE exists = ", file.exists(CONDA_EXE))

if (!dir.exists(raw_dir)) stop("找不到 raw_dir：", raw_dir)
if (!dir.exists(FASTQC_DIR)) stop("找不到 FASTQC_DIR：", FASTQC_DIR)
if (!file.exists(FASTQC_BAT)) stop("找不到 FASTQC_BAT：", FASTQC_BAT)
if (!file.exists(CONDA_EXE)) stop("找不到 CONDA_EXE：", CONDA_EXE)

# -----------------------------
# 2. 36 个 FASTQ 清单
# -----------------------------
expected_samples <- c(
  paste0("CON", 1:6),
  paste0("INJS", 1:6),
  paste0("INJL", 1:6)
)

all_fastq <- list.files(
  raw_dir,
  pattern = "\\.(fastq|fq)\\.gz$",
  recursive = FALSE,
  full.names = TRUE
)

if (length(all_fastq) == 0) stop("在 raw_dir 中没有找到 FASTQ：", raw_dir)

inventory <- data.frame(
  absolute_path = normalizePath(all_fastq, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE
)
inventory$filename <- basename(inventory$absolute_path)
inventory$sample_id <- sub("^(.*?)_(R1|R2)\\.(fastq|fq)\\.gz$", "\\1", inventory$filename, perl = TRUE)
inventory$read <- sub("^.*?_(R1|R2)\\.(fastq|fq)\\.gz$", "\\1", inventory$filename, perl = TRUE)
inventory <- inventory[inventory$sample_id %in% expected_samples, , drop = FALSE]
inventory <- inventory[order(match(inventory$sample_id, expected_samples), inventory$read), , drop = FALSE]

pair_check <- inventory %>%
  group_by(sample_id) %>%
  summarise(
    n_files = n(),
    has_R1 = any(read == "R1"),
    has_R2 = any(read == "R2"),
    n_R1 = sum(read == "R1"),
    n_R2 = sum(read == "R2"),
    paired_ok = has_R1 & has_R2 & n_R1 == 1 & n_R2 == 1,
    .groups = "drop"
  ) %>%
  arrange(match(sample_id, expected_samples))

if (!all(expected_samples %in% pair_check$sample_id)) {
  stop("以下样本未检测到 FASTQ：", paste(setdiff(expected_samples, pair_check$sample_id), collapse = ", "))
}
if (!all(pair_check$paired_ok)) {
  stop("以下样本 R1/R2 配对异常：", paste(pair_check$sample_id[!pair_check$paired_ok], collapse = ", "))
}

inventory_file <- file.path(tables_dir, "stepE_pig_early_fastq_inventory.csv")
pair_file <- file.path(tables_dir, "stepE_pig_early_fastq_pair_check.csv")
write.csv(inventory, inventory_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pair_check, pair_file, row.names = FALSE, fileEncoding = "UTF-8")
append_log("FASTQ inventory saved: ", inventory_file)
append_log("Pair check saved: ", pair_file)

# -----------------------------
# 3. 先清理旧的目标输出，避免混淆
# -----------------------------
fastq_base <- sub("\\.(fastq|fq)\\.gz$", "", inventory$filename, ignore.case = TRUE)
old_html_raw <- file.path(raw_dir, paste0(fastq_base, "_fastqc.html"))
old_zip_raw  <- file.path(raw_dir, paste0(fastq_base, "_fastqc.zip"))
old_html_out <- file.path(fastqc_out_dir, paste0(fastq_base, "_fastqc.html"))
old_zip_out  <- file.path(fastqc_out_dir, paste0(fastq_base, "_fastqc.zip"))
old_targets <- c(old_html_raw, old_zip_raw, old_html_out, old_zip_out)
old_targets <- old_targets[file.exists(old_targets)]
if (length(old_targets) > 0) {
  unlink(old_targets, force = TRUE)
  append_log("Removed old FastQC outputs: ", length(old_targets))
}

# 同样清理旧的 MultiQC 输出
if (dir.exists(multiqc_out_dir)) {
  old_multiqc <- list.files(multiqc_out_dir, full.names = TRUE, recursive = TRUE)
  if (length(old_multiqc) > 0) unlink(old_multiqc, recursive = TRUE, force = TRUE)
}

# -----------------------------
# 4. 逐个运行 FastQC（不传任何选项）
# -----------------------------
fastqc_dir_win <- normalizePath(FASTQC_DIR, winslash = "\\", mustWork = TRUE)
fastqc_detail_list <- vector("list", nrow(inventory))

move_one_file <- function(src, dst) {
  if (!file.exists(src)) return(FALSE)
  if (file.exists(dst)) file.remove(dst)
  ok <- file.rename(src, dst)
  if (!ok) {
    ok <- file.copy(src, dst, overwrite = TRUE)
    if (ok && file.exists(src)) file.remove(src)
  }
  ok && file.exists(dst)
}

for (i in seq_len(nrow(inventory))) {
  fastq <- inventory$absolute_path[i]
  fastq_file <- inventory$filename[i]
  base_name <- sub("\\.(fastq|fq)\\.gz$", "", fastq_file, ignore.case = TRUE)
  input_win <- normalizePath(fastq, winslash = "\\", mustWork = TRUE)

  expected_html_raw <- file.path(raw_dir, paste0(base_name, "_fastqc.html"))
  expected_zip_raw  <- file.path(raw_dir, paste0(base_name, "_fastqc.zip"))
  expected_html_out <- file.path(fastqc_out_dir, paste0(base_name, "_fastqc.html"))
  expected_zip_out  <- file.path(fastqc_out_dir, paste0(base_name, "_fastqc.zip"))

  if (file.exists(expected_html_raw)) file.remove(expected_html_raw)
  if (file.exists(expected_zip_raw)) file.remove(expected_zip_raw)
  if (file.exists(expected_html_out)) file.remove(expected_html_out)
  if (file.exists(expected_zip_out)) file.remove(expected_zip_out)

  cmd <- paste0(
    'cmd.exe /c cd /d ', fastqc_dir_win,
    ' && run_fastqc.bat "', input_win, '"'
  )

  append_log("FastQC [", i, "/", nrow(inventory), "] : ", fastq_file)

  res <- tryCatch(
    system(cmd, intern = TRUE, ignore.stderr = FALSE),
    warning = function(w) w,
    error = function(e) e
  )

  exit_status <- 0L
  out_lines <- character(0)
  warning_or_error <- ""

  if (inherits(res, "warning")) {
    warning_or_error <- conditionMessage(res)
    res2 <- tryCatch(system(cmd, intern = TRUE, ignore.stderr = FALSE), error = function(e) e)
    if (inherits(res2, "error")) {
      out_lines <- res2$message
      exit_status <- 999L
    } else {
      out_lines <- res2
      exit_status <- ifelse(is.null(attr(res2, "status")), 0L, as.integer(attr(res2, "status")))
    }
  } else if (inherits(res, "error")) {
    warning_or_error <- res$message
    out_lines <- res$message
    exit_status <- 999L
  } else {
    out_lines <- res
    exit_status <- ifelse(is.null(attr(res, "status")), 0L, as.integer(attr(res, "status")))
  }

  stdout_file_i <- file.path(logs_dir, paste0("stepE_fastqc_", base_name, "_stdout_stderr.txt"))
  writeLines(out_lines, con = stdout_file_i, useBytes = TRUE)

  html_in_raw <- file.exists(expected_html_raw)
  zip_in_raw  <- file.exists(expected_zip_raw)

  moved_html <- FALSE
  moved_zip <- FALSE
  if (html_in_raw) moved_html <- move_one_file(expected_html_raw, expected_html_out)
  if (zip_in_raw)  moved_zip  <- move_one_file(expected_zip_raw, expected_zip_out)

  ok <- (exit_status == 0L) && moved_html && moved_zip && file.exists(expected_html_out) && file.exists(expected_zip_out)

  fastqc_detail_list[[i]] <- data.frame(
    fastq_file = fastq_file,
    sample_id = inventory$sample_id[i],
    read = inventory$read[i],
    exit_status = exit_status,
    html_created_in_raw = html_in_raw,
    zip_created_in_raw = zip_in_raw,
    html_moved_to_outdir = moved_html,
    zip_moved_to_outdir = moved_zip,
    html_exists_final = file.exists(expected_html_out),
    zip_exists_final = file.exists(expected_zip_out),
    ok = ok,
    warning_or_error = warning_or_error,
    stdout_file = stdout_file_i,
    stringsAsFactors = FALSE
  )

  if (!ok) {
    detail_df_partial <- bind_rows(fastqc_detail_list[seq_len(i)])
    detail_file_partial <- file.path(tables_dir, "stepE_pig_early_fastqc_detail.csv")
    write.csv(detail_df_partial, detail_file_partial, row.names = FALSE, fileEncoding = "UTF-8")
    append_log("FastQC failed at ", fastq_file, ". Detail saved: ", detail_file_partial)
    stop("FastQC 未真正完成：请先检查 stepE_pig_early_fastqc_detail.csv 和对应 stdout 日志。")
  }
}

fastqc_detail <- bind_rows(fastqc_detail_list)
fastqc_detail_file <- file.path(tables_dir, "stepE_pig_early_fastqc_detail.csv")
write.csv(fastqc_detail, fastqc_detail_file, row.names = FALSE, fileEncoding = "UTF-8")
append_log("FastQC detail saved: ", fastqc_detail_file)

fastqc_html <- list.files(fastqc_out_dir, pattern = "_fastqc\\.html$", full.names = FALSE)
fastqc_zip  <- list.files(fastqc_out_dir, pattern = "_fastqc\\.zip$", full.names = FALSE)
fastqc_output_inventory <- data.frame(
  type = c(rep("html", length(fastqc_html)), rep("zip", length(fastqc_zip))),
  file = c(fastqc_html, fastqc_zip),
  stringsAsFactors = FALSE
)
fastqc_output_file <- file.path(tables_dir, "stepE_pig_early_fastqc_output_inventory.csv")
write.csv(fastqc_output_inventory, fastqc_output_file, row.names = FALSE, fileEncoding = "UTF-8")
append_log("FastQC output inventory saved: ", fastqc_output_file)

if (length(fastqc_html) != nrow(inventory) || length(fastqc_zip) != nrow(inventory)) {
  stop("FastQC 产物数量异常：html=", length(fastqc_html), ", zip=", length(fastqc_zip), ", expected=", nrow(inventory))
}

# -----------------------------
# 5. 运行 MultiQC（用 conda run）
# -----------------------------
append_log("Running MultiQC by conda run...")
multiqc_cmd <- c(
  "run", "-n", MULTIQC_ENV_NAME,
  "multiqc",
  normalizePath(fastqc_out_dir, winslash = "/", mustWork = TRUE),
  "--outdir", normalizePath(multiqc_out_dir, winslash = "/", mustWork = TRUE),
  "--filename", "multiqc_report"
)

multiqc_cmd_file <- file.path(logs_dir, "stepE_multiqc_command.txt")
writeLines(paste(shQuote(CONDA_EXE), paste(shQuote(multiqc_cmd), collapse = " ")), con = multiqc_cmd_file, useBytes = TRUE)

multiqc_res <- tryCatch(
  system2(CONDA_EXE, args = multiqc_cmd, stdout = TRUE, stderr = TRUE),
  error = function(e) e
)

multiqc_stdout_file <- file.path(logs_dir, "stepE_multiqc_stdout_stderr.txt")
if (inherits(multiqc_res, "error")) {
  stop("MultiQC 运行失败：", multiqc_res$message)
}
writeLines(multiqc_res, con = multiqc_stdout_file, useBytes = TRUE)

multiqc_report_html <- file.path(multiqc_out_dir, "multiqc_report.html")
multiqc_files <- list.files(multiqc_out_dir, recursive = TRUE, full.names = FALSE)
multiqc_output_inventory <- data.frame(file = multiqc_files, stringsAsFactors = FALSE)
multiqc_output_file <- file.path(tables_dir, "stepE_pig_early_multiqc_output_inventory.csv")
write.csv(multiqc_output_inventory, multiqc_output_file, row.names = FALSE, fileEncoding = "UTF-8")
append_log("MultiQC output inventory saved: ", multiqc_output_file)

if (!file.exists(multiqc_report_html)) {
  stop("MultiQC 未生成 multiqc_report.html，请检查 stepE_multiqc_stdout_stderr.txt")
}

# -----------------------------
# 6. summary
# -----------------------------
summary_df <- data.frame(
  metric = c(
    "project_root",
    "raw_dir",
    "n_input_fastq_files",
    "n_unique_samples",
    "fastqc_bat",
    "conda_exe",
    "multiqc_env_name",
    "n_fastqc_ok",
    "n_fastqc_html_final",
    "n_fastqc_zip_final",
    "multiqc_report_exists",
    "status"
  ),
  value = c(
    project_root,
    raw_dir,
    nrow(inventory),
    length(unique(inventory$sample_id)),
    FASTQC_BAT,
    CONDA_EXE,
    MULTIQC_ENV_NAME,
    sum(fastqc_detail$ok, na.rm = TRUE),
    length(fastqc_html),
    length(fastqc_zip),
    file.exists(multiqc_report_html),
    "Step E completed successfully"
  ),
  stringsAsFactors = FALSE
)
summary_file <- file.path(tables_dir, "stepE_pig_early_fastqc_multiqc_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE, fileEncoding = "UTF-8")
append_log("Summary saved: ", summary_file)

save(
  inventory,
  pair_check,
  fastqc_detail,
  fastqc_output_inventory,
  multiqc_output_inventory,
  summary_df,
  file = file.path(objects_dir, "stepE_pig_early_fastqc_multiqc_workspace.RData")
)
append_log("Workspace saved.")
append_log("Step E finished successfully.")

cat("===== Step E summary =====\n")
print(summary_df, row.names = FALSE)
cat("\n输出文件：\n")
cat(summary_file, "\n")
cat(fastqc_detail_file, "\n")
cat(fastqc_output_file, "\n")
cat(multiqc_output_file, "\n")
cat(log_file, "\n")
