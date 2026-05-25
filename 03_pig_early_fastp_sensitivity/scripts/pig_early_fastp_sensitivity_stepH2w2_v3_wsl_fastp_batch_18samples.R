# ============================================================
# Step H2w2 v3: pig early fastp sensitivity - batch trimming via WSL fastp
# 项目：ACL × pig early/chronic × senescence
#
# 这一步是干什么：
#   在 Step H2w1 已经成功的基础上，从 Windows R 调用 WSL，
#   激活 WSL 中的 conda 环境 fastp_env，对 pig early 的 18 个核心样本
#   批量运行 fastp trimming。
#
# v3 关键修改：
#   - 不再把一整串 bash 命令塞进 wsl.exe bash -lc "..."。
#   - 改为：R 先写出 .sh 文件，再让 WSL 执行这个 .sh 文件。
#   - 这样可以避开 Windows R -> wsl.exe -> bash -lc 的引号/转义问题。
#   - 支持断点续跑：如果某个样本四个输出文件已经完整存在，则自动跳过。
#
# 使用方法：
#   source("E:/R/ACLsenescence2/pig_early_fastp_sensitivity_stepH2w2_v3_wsl_fastp_batch_18samples.R")
# ============================================================

cat("\n===== Step H2w2 v3: batch fastp trimming for 18 pig early samples =====\n\n")

# ----------------------------
# 0. 固定参数
# ----------------------------
project_root <- "E:/R/ACLsenescence2"
fastp_env_name <- "fastp_env"
threads <- 4
skip_if_complete <- TRUE

branch_dir <- file.path(
  project_root,
  "rebuild_submission",
  "02_pig_early_fastp_sensitivity"
)

tables_dir  <- file.path(branch_dir, "tables")
logs_dir    <- file.path(branch_dir, "logs")
objects_dir <- file.path(branch_dir, "objects")
script_dir  <- file.path(branch_dir, "wsl_scripts")

trimmed_default_dir <- file.path(branch_dir, "trimmed_fastq")
reports_default_dir <- file.path(branch_dir, "fastp_reports")

for (d in c(tables_dir, logs_dir, objects_dir, script_dir, trimmed_default_dir, reports_default_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

command_sheet <- file.path(
  tables_dir,
  "stepH1_fastp_command_sheet_18samples.csv"
)

status_file <- file.path(
  tables_dir,
  "stepH2w2_v3_fastp_batch_status_18samples.csv"
)

summary_file <- file.path(
  tables_dir,
  "stepH2w2_v3_fastp_batch_summary.csv"
)

command_file <- file.path(
  logs_dir,
  "stepH2w2_v3_fastp_batch_commands.txt"
)

main_log_file <- file.path(
  logs_dir,
  "stepH2w2_v3_fastp_batch_main_log.txt"
)

config_rds <- file.path(
  objects_dir,
  "stepH2w2_v3_wsl_fastp_batch_runtime_config.rds"
)

preflight_sh <- file.path(script_dir, "stepH2w2_v3_preflight_fastp_check.sh")
preflight_log <- file.path(logs_dir, "stepH2w2_v3_preflight_wsl_fastp_check.log")

# ----------------------------
# 1. 小工具函数
# ----------------------------

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n", sep = "", file = main_log_file, append = TRUE)
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("找不到 %s：%s", label, path), call. = FALSE)
  }
}

find_col_required <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop(
      sprintf(
        "在 command sheet 中找不到必需列：%s。当前列名为：%s",
        label,
        paste(names(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  hit[1]
}

find_col_optional <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_file_size <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  as.numeric(file.info(path)$size)
}

file_complete <- function(path) {
  file.exists(path) && isTRUE(safe_file_size(path) > 0)
}

sanitize_id <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

# Windows 路径转换成 WSL 路径，例如 E:/R/a.fastq.gz -> /mnt/e/R/a.fastq.gz
windows_to_wsl <- function(path) {
  p <- gsub("\\\\", "/", path)

  if (grepl("^/mnt/[A-Za-z]/", p)) {
    return(p)
  }

  p2 <- suppressWarnings(normalizePath(p, winslash = "/", mustWork = FALSE))
  p2 <- gsub("\\\\", "/", p2)

  if (grepl("^[A-Za-z]:/", p2)) {
    drive <- tolower(substr(p2, 1, 1))
    rest  <- substring(p2, 3)
    return(paste0("/mnt/", drive, rest))
  }

  stop(sprintf("无法转换为 WSL 路径：%s", path), call. = FALSE)
}

# WSL bash 里的安全引号，用于路径和 sample_id
sh_quote_wsl <- function(x) {
  paste0("'", gsub("'", "'\"'\"'", x, fixed = TRUE), "'")
}

# 写出 LF 换行的 bash 脚本，避免 Windows CRLF 带来的隐患
write_text_lf <- function(lines, path) {
  txt <- paste(lines, collapse = "\n")
  txt <- paste0(txt, "\n")
  con <- base::file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(txt), con = con)
}

# 运行一个已经写好的 WSL shell 脚本文件
run_wsl_script <- function(script_file_win, log_file, wsl_bin) {
  script_file_wsl <- windows_to_wsl(script_file_win)

  out <- tryCatch(
    {
      suppressWarnings(
        system2(
          command = wsl_bin,
          args = c("bash", script_file_wsl),
          stdout = TRUE,
          stderr = TRUE
        )
      )
    },
    error = function(e) {
      structure(
        paste0("R/system2 error: ", conditionMessage(e)),
        status = 999
      )
    }
  )

  status <- attr(out, "status")
  if (is.null(status)) status <- 0

  writeLines(as.character(out), con = log_file, useBytes = TRUE)
  list(status = as.integer(status), output = as.character(out))
}

# WSL 中寻找 conda.sh：优先 miniforge3，其次 miniconda3/anaconda3
conda_block_lines <- c(
  'if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then',
  '  source "$HOME/miniforge3/etc/profile.d/conda.sh"',
  'elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then',
  '  source "$HOME/miniconda3/etc/profile.d/conda.sh"',
  'elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then',
  '  source "$HOME/anaconda3/etc/profile.d/conda.sh"',
  'else',
  '  echo "ERROR: conda.sh not found under ~/miniforge3, ~/miniconda3, or ~/anaconda3"',
  '  exit 91',
  'fi'
)

# ----------------------------
# 2. 读取 Step H1 command sheet
# ----------------------------

if (file.exists(main_log_file)) file.remove(main_log_file)

log_msg("读取 Step H1 command sheet：", command_sheet)
stop_if_missing(command_sheet, "Step H1 command sheet")

cmd_df <- read.csv(command_sheet, stringsAsFactors = FALSE, check.names = FALSE)

sample_col <- find_col_required(
  cmd_df,
  candidates = c("sample_id", "sample", "Sample", "SampleID"),
  label = "sample_id"
)

group_col <- find_col_optional(
  cmd_df,
  candidates = c("group", "core_group", "condition", "Group")
)

r1_col <- find_col_required(
  cmd_df,
  candidates = c("absolute_path.R1", "absolute_path_R1", "R1", "r1", "fastq_R1", "fq1"),
  label = "R1 FASTQ path"
)

r2_col <- find_col_required(
  cmd_df,
  candidates = c("absolute_path.R2", "absolute_path_R2", "R2", "r2", "fastq_R2", "fq2"),
  label = "R2 FASTQ path"
)

trim_r1_col <- find_col_optional(cmd_df, c("trimmed_R1", "trimmed.R1", "out_R1", "output_R1"))
trim_r2_col <- find_col_optional(cmd_df, c("trimmed_R2", "trimmed.R2", "out_R2", "output_R2"))
html_col    <- find_col_optional(cmd_df, c("fastp_html", "html", "fastp_report_html"))
json_col    <- find_col_optional(cmd_df, c("fastp_json", "json", "fastp_report_json"))

if (anyDuplicated(cmd_df[[sample_col]]) > 0) {
  dup_ids <- unique(cmd_df[[sample_col]][duplicated(cmd_df[[sample_col]])])
  stop(sprintf("command sheet 中存在重复 sample_id：%s", paste(dup_ids, collapse = ", ")), call. = FALSE)
}

n_samples <- nrow(cmd_df)
log_msg("读取到样本数：", n_samples)
if (!identical(n_samples, 18L)) {
  log_msg("提醒：当前 command sheet 样本数不是 18，请确认是否符合预期。脚本仍会按表内样本继续运行。")
}

# ----------------------------
# 3. 检查 wsl.exe 和 fastp_env
# ----------------------------

wsl_bin <- Sys.which("wsl.exe")
if (identical(unname(wsl_bin), "")) {
  wsl_bin <- Sys.which("wsl")
}
if (identical(unname(wsl_bin), "")) {
  stop("Windows R 找不到 wsl.exe。请先确认 PowerShell/CMD 中可以运行 wsl。", call. = FALSE)
}
wsl_bin <- unname(wsl_bin)
log_msg("检测到 WSL 调用程序：", wsl_bin)

preflight_lines <- c(
  '#!/usr/bin/env bash',
  'set -euo pipefail',
  'echo "===== Step H2w2 v3 WSL / fastp preflight ====="',
  'echo "WSL user: $(whoami)"',
  'echo "WSL home: $HOME"',
  conda_block_lines,
  sprintf('conda activate %s', fastp_env_name),
  'echo "Conda env: $CONDA_DEFAULT_ENV"',
  'echo "fastp path:"',
  'which fastp',
  'echo "fastp version:"',
  'fastp --version'
)

write_text_lf(preflight_lines, preflight_sh)

log_msg("开始 preflight：通过 .sh 文件检查 WSL 中的 fastp_env / fastp。")
preflight_res <- run_wsl_script(preflight_sh, preflight_log, wsl_bin)

if (!identical(preflight_res$status, 0L)) {
  cat("\npreflight 失败，下面是日志内容：\n")
  cat(paste(preflight_res$output, collapse = "\n"), "\n")
  stop(
    sprintf("Step H2w2 v3 preflight 失败，exit status = %s。请查看日志：%s", preflight_res$status, preflight_log),
    call. = FALSE
  )
}
log_msg("preflight 成功。日志：", preflight_log)

# ----------------------------
# 4. 构建每个样本的任务表
# ----------------------------

make_default_outputs <- function(sample_id) {
  list(
    trim_r1 = file.path(trimmed_default_dir, paste0(sample_id, ".R1.trimmed.fastq.gz")),
    trim_r2 = file.path(trimmed_default_dir, paste0(sample_id, ".R2.trimmed.fastq.gz")),
    html    = file.path(reports_default_dir, paste0(sample_id, ".fastp.html")),
    json    = file.path(reports_default_dir, paste0(sample_id, ".fastp.json"))
  )
}

nonempty_cell <- function(df, col, i) {
  if (is.na(col)) return(FALSE)
  value <- df[[col]][i]
  !is.na(value) && nzchar(value)
}

task_list <- vector("list", n_samples)

for (i in seq_len(n_samples)) {
  sample_id <- cmd_df[[sample_col]][i]
  default_out <- make_default_outputs(sample_id)

  trim_r1 <- if (nonempty_cell(cmd_df, trim_r1_col, i)) cmd_df[[trim_r1_col]][i] else default_out$trim_r1
  trim_r2 <- if (nonempty_cell(cmd_df, trim_r2_col, i)) cmd_df[[trim_r2_col]][i] else default_out$trim_r2
  html    <- if (nonempty_cell(cmd_df, html_col, i))    cmd_df[[html_col]][i]    else default_out$html
  json    <- if (nonempty_cell(cmd_df, json_col, i))    cmd_df[[json_col]][i]    else default_out$json

  task_list[[i]] <- data.frame(
    sample_id = sample_id,
    group = if (!is.na(group_col)) cmd_df[[group_col]][i] else NA_character_,
    r1 = cmd_df[[r1_col]][i],
    r2 = cmd_df[[r2_col]][i],
    trimmed_r1 = trim_r1,
    trimmed_r2 = trim_r2,
    fastp_html = html,
    fastp_json = json,
    log_file = file.path(logs_dir, paste0("stepH2w2_v3_fastp_", sanitize_id(sample_id), "_stdout_stderr.txt")),
    sh_file = file.path(script_dir, paste0("stepH2w2_v3_fastp_", sanitize_id(sample_id), ".sh")),
    stringsAsFactors = FALSE
  )
}

tasks <- do.call(rbind, task_list)

for (p in unique(c(dirname(tasks$trimmed_r1), dirname(tasks$trimmed_r2), dirname(tasks$fastp_html), dirname(tasks$fastp_json)))) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

saveRDS(
  list(
    timestamp = Sys.time(),
    project_root = project_root,
    branch_dir = branch_dir,
    command_sheet = command_sheet,
    fastp_env_name = fastp_env_name,
    threads = threads,
    skip_if_complete = skip_if_complete,
    wsl_bin = wsl_bin,
    preflight_sh = preflight_sh,
    preflight_log = preflight_log,
    tasks = tasks
  ),
  config_rds
)

# ----------------------------
# 5. 批量运行 fastp
# ----------------------------

status_rows <- vector("list", nrow(tasks))
all_command_lines <- character(0)

for (i in seq_len(nrow(tasks))) {
  one <- tasks[i, , drop = FALSE]
  sample_id <- one$sample_id

  log_msg("开始样本 ", i, "/", nrow(tasks), "：", sample_id)

  input_r1_ok <- file_complete(one$r1)
  input_r2_ok <- file_complete(one$r2)

  out_complete_before <- all(c(
    file_complete(one$trimmed_r1),
    file_complete(one$trimmed_r2),
    file_complete(one$fastp_html),
    file_complete(one$fastp_json)
  ))

  action <- NA_character_
  exit_status <- NA_integer_
  started_at <- Sys.time()
  ended_at <- as.POSIXct(NA)

  if (!input_r1_ok || !input_r2_ok) {
    action <- "input_missing_skip"
    exit_status <- 998L
    log_msg("样本 ", sample_id, " 输入 FASTQ 不完整，跳过。")
  } else if (isTRUE(skip_if_complete) && isTRUE(out_complete_before)) {
    action <- "skipped_existing_complete"
    exit_status <- 0L
    log_msg("样本 ", sample_id, " 已有完整输出，跳过。")
  } else {
    action <- "run_fastp"

    # 删除残留的部分输出，避免误判
    unlink(c(one$trimmed_r1, one$trimmed_r2, one$fastp_html, one$fastp_json), force = TRUE)

    r1_wsl <- windows_to_wsl(one$r1)
    r2_wsl <- windows_to_wsl(one$r2)
    trim_r1_wsl <- windows_to_wsl(one$trimmed_r1)
    trim_r2_wsl <- windows_to_wsl(one$trimmed_r2)
    html_wsl <- windows_to_wsl(one$fastp_html)
    json_wsl <- windows_to_wsl(one$fastp_json)

    sample_lines <- c(
      '#!/usr/bin/env bash',
      'set -euo pipefail',
      'echo "===== fastp sample start ====="',
      sprintf('echo "sample_id: %s"', sample_id),
      conda_block_lines,
      sprintf('conda activate %s', fastp_env_name),
      'echo "Conda env: $CONDA_DEFAULT_ENV"',
      'echo "fastp path: $(which fastp)"',
      sprintf(
        'mkdir -p %s %s %s %s',
        sh_quote_wsl(dirname(trim_r1_wsl)),
        sh_quote_wsl(dirname(trim_r2_wsl)),
        sh_quote_wsl(dirname(html_wsl)),
        sh_quote_wsl(dirname(json_wsl))
      ),
      sprintf(
        paste(
          'fastp',
          '-i %s',
          '-I %s',
          '-o %s',
          '-O %s',
          '-h %s',
          '-j %s',
          '--detect_adapter_for_pe',
          '--thread %d'
        ),
        sh_quote_wsl(r1_wsl),
        sh_quote_wsl(r2_wsl),
        sh_quote_wsl(trim_r1_wsl),
        sh_quote_wsl(trim_r2_wsl),
        sh_quote_wsl(html_wsl),
        sh_quote_wsl(json_wsl),
        threads
      ),
      'echo "===== fastp sample finished ====="'
    )

    write_text_lf(sample_lines, one$sh_file)

    all_command_lines <- c(
      all_command_lines,
      paste0(
        "\n# ===== ", sample_id, " =====\n",
        paste(sample_lines, collapse = "\n"),
        "\n"
      )
    )

    res <- run_wsl_script(one$sh_file, one$log_file, wsl_bin)
    exit_status <- res$status
    ended_at <- Sys.time()

    log_msg("样本 ", sample_id, " fastp exit status：", exit_status)
  }

  if (is.na(ended_at)) ended_at <- Sys.time()

  output_r1_ok <- file_complete(one$trimmed_r1)
  output_r2_ok <- file_complete(one$trimmed_r2)
  output_html_ok <- file_complete(one$fastp_html)
  output_json_ok <- file_complete(one$fastp_json)
  sample_ok <- identical(exit_status, 0L) && output_r1_ok && output_r2_ok && output_html_ok && output_json_ok

  status_rows[[i]] <- data.frame(
    sample_id = sample_id,
    group = one$group,
    action = action,
    exit_status = exit_status,
    sample_ok = sample_ok,
    input_r1_exists = file.exists(one$r1),
    input_r2_exists = file.exists(one$r2),
    input_r1_size_bytes = safe_file_size(one$r1),
    input_r2_size_bytes = safe_file_size(one$r2),
    trimmed_r1_exists = file.exists(one$trimmed_r1),
    trimmed_r2_exists = file.exists(one$trimmed_r2),
    fastp_html_exists = file.exists(one$fastp_html),
    fastp_json_exists = file.exists(one$fastp_json),
    trimmed_r1_size_bytes = safe_file_size(one$trimmed_r1),
    trimmed_r2_size_bytes = safe_file_size(one$trimmed_r2),
    fastp_html_size_bytes = safe_file_size(one$fastp_html),
    fastp_json_size_bytes = safe_file_size(one$fastp_json),
    r1 = one$r1,
    r2 = one$r2,
    trimmed_r1 = one$trimmed_r1,
    trimmed_r2 = one$trimmed_r2,
    fastp_html = one$fastp_html,
    fastp_json = one$fastp_json,
    log_file = one$log_file,
    sh_file = one$sh_file,
    started_at = format(started_at, "%Y-%m-%d %H:%M:%S"),
    ended_at = format(ended_at, "%Y-%m-%d %H:%M:%S"),
    elapsed_min = round(as.numeric(difftime(ended_at, started_at, units = "mins")), 3),
    stringsAsFactors = FALSE
  )

  status_df_tmp <- do.call(rbind, status_rows[seq_len(i)])
  write.csv(status_df_tmp, status_file, row.names = FALSE, fileEncoding = "UTF-8")
}

status_df <- do.call(rbind, status_rows)
write.csv(status_df, status_file, row.names = FALSE, fileEncoding = "UTF-8")
writeLines(all_command_lines, command_file, useBytes = TRUE)

# ----------------------------
# 6. 输出 summary
# ----------------------------

batch_ok <- all(status_df$sample_ok)

summary_df <- data.frame(
  metric = c(
    "project_root",
    "branch_dir",
    "command_sheet_file",
    "fastp_env_name",
    "threads",
    "n_samples_in_command_sheet",
    "n_samples_ok",
    "n_samples_failed_or_incomplete",
    "n_run_fastp",
    "n_skipped_existing_complete",
    "batch_ok",
    "preflight_sh",
    "preflight_log",
    "status_file",
    "summary_file",
    "command_file",
    "main_log_file",
    "config_rds",
    "status"
  ),
  value = c(
    project_root,
    branch_dir,
    command_sheet,
    fastp_env_name,
    as.character(threads),
    as.character(nrow(status_df)),
    as.character(sum(status_df$sample_ok)),
    as.character(sum(!status_df$sample_ok)),
    as.character(sum(status_df$action == "run_fastp", na.rm = TRUE)),
    as.character(sum(status_df$action == "skipped_existing_complete", na.rm = TRUE)),
    as.character(batch_ok),
    preflight_sh,
    preflight_log,
    status_file,
    summary_file,
    command_file,
    main_log_file,
    config_rds,
    if (isTRUE(batch_ok)) "Step H2w2 v3 completed successfully" else "Step H2w2 v3 completed with failed/incomplete samples"
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_file, row.names = FALSE, fileEncoding = "UTF-8")

cat("\n===== Step H2w2 v3 summary =====\n")
print(summary_df, row.names = FALSE)

cat("\n===== Step H2w2 v3 per-sample status =====\n")
print(status_df[, c("sample_id", "group", "action", "exit_status", "sample_ok")], row.names = FALSE)

cat("\n输出文件：\n")
cat(status_file, "\n")
cat(summary_file, "\n")
cat(command_file, "\n")
cat(main_log_file, "\n")
cat(config_rds, "\n")

if (isTRUE(batch_ok)) {
  cat("\n结果：H2w2 v3 成功。18 个样本的 fastp trimming 输出完整。下一步可以进入 trimmed FASTQ 的比对准备。\n")
} else {
  cat("\n结果：H2w2 v3 仍有样本失败或输出不完整。请优先查看 stepH2w2_v3_fastp_batch_status_18samples.csv 中 sample_ok=FALSE 的行。\n")
}

cat("\n===== End of Step H2w2 v3 =====\n")
