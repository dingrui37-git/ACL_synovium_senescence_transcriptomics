# ============================================================
# Step H2w1: pig early fastp sensitivity - WSL fastp smoke test
# 项目：ACL × pig early/chronic × senescence
# 目的：
#   1) 从 Windows 的 R 调用 WSL
#   2) 激活 WSL 里的 conda 环境 fastp_env
#   3) 检查 fastp 是否可用
#   4) 只对 CON1 做一个小规模 smoke test，确认调用链和输出路径都正常
#
# 使用方法：
#   1) 把本脚本放到：E:/R/ACLsenescence2
#   2) 在 R/RStudio 中运行：
#        source("E:/R/ACLsenescence2/pig_early_fastp_sensitivity_stepH2w1_wsl_fastp_smoke_test.R")
#
# 注意：
#   - 这个脚本不会覆盖正式 H2w2 的 trimmed FASTQ。
#   - smoke test 只处理 CON1 的前 100000 reads，用于验证流程是否跑通。
# ============================================================

cat("\n===== Step H2w1: WSL fastp smoke test for CON1 =====\n\n")

# ----------------------------
# 0. 固定项目路径
# ----------------------------
project_root <- "E:/R/ACLsenescence2"

branch_dir <- file.path(
  project_root,
  "rebuild_submission",
  "02_pig_early_fastp_sensitivity"
)

tables_dir <- file.path(branch_dir, "tables")
logs_dir   <- file.path(branch_dir, "logs")
smoke_dir  <- file.path(branch_dir, "smoke_test_CON1")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(smoke_dir, recursive = TRUE, showWarnings = FALSE)

command_sheet <- file.path(
  tables_dir,
  "stepH1_fastp_command_sheet_18samples.csv"
)

summary_csv <- file.path(
  tables_dir,
  "stepH2w1_wsl_fastp_smoke_test_summary.csv"
)

file_check_csv <- file.path(
  tables_dir,
  "stepH2w1_wsl_fastp_smoke_test_file_check.csv"
)

version_log <- file.path(
  logs_dir,
  "stepH2w1_wsl_fastp_version_check.log"
)

smoke_log <- file.path(
  logs_dir,
  "stepH2w1_CON1_smoke_fastp.log"
)

# ----------------------------
# 1. 小工具函数
# ----------------------------

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      sprintf("找不到 %s：%s", label, path),
      call. = FALSE
    )
  }
}

find_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    stop(
      sprintf(
        "在 command sheet 中找不到列：%s。当前列名为：%s",
        label,
        paste(names(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  hit[1]
}

# 把 Windows 路径转换成 WSL 可识别的 /mnt/e/... 路径
windows_to_wsl <- function(path) {
  p <- gsub("\\\\", "/", path)

  # 如果已经是 WSL 路径，则直接返回
  if (grepl("^/mnt/[A-Za-z]/", p)) {
    return(p)
  }

  # 尽量标准化 Windows 路径
  p2 <- suppressWarnings(normalizePath(p, winslash = "/", mustWork = FALSE))
  p2 <- gsub("\\\\", "/", p2)

  if (grepl("^[A-Za-z]:/", p2)) {
    drive <- tolower(substr(p2, 1, 1))
    rest  <- substring(p2, 3)
    return(paste0("/mnt/", drive, rest))
  }

  stop(sprintf("无法转换为 WSL 路径：%s", path), call. = FALSE)
}

# WSL bash 里的安全引号
sh_quote_wsl <- function(x) {
  paste0("'", gsub("'", "'\"'\"'", x, fixed = TRUE), "'")
}

# 在 Windows R 中运行 WSL bash 命令，并把 stdout/stderr 写入 log
run_wsl_bash <- function(command, log_file, wsl_bin) {
  out <- tryCatch(
    {
      suppressWarnings(
        system2(
          command = wsl_bin,
          args = c("bash", "-lc", shQuote(command)),
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

safe_file_size <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  as.numeric(file.info(path)$size)
}

# WSL 中寻找 conda.sh 的代码块：
# 优先 miniforge3；若用户之后换成 miniconda3 / anaconda3，也能自动尝试。
conda_block <- paste(
  'if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then',
  '  source "$HOME/miniforge3/etc/profile.d/conda.sh";',
  'elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then',
  '  source "$HOME/miniconda3/etc/profile.d/conda.sh";',
  'elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then',
  '  source "$HOME/anaconda3/etc/profile.d/conda.sh";',
  'else',
  '  echo "ERROR: conda.sh not found under ~/miniforge3, ~/miniconda3, or ~/anaconda3";',
  '  exit 91;',
  'fi',
  sep = " "
)

# ----------------------------
# 2. 读取 H1 command sheet，并锁定 CON1
# ----------------------------

stop_if_missing(command_sheet, "Step H1 command sheet")

cmd_df <- read.csv(command_sheet, stringsAsFactors = FALSE, check.names = FALSE)

sample_col <- find_col(
  cmd_df,
  candidates = c("sample_id", "sample", "Sample", "SampleID"),
  label = "sample_id"
)

r1_col <- find_col(
  cmd_df,
  candidates = c("absolute_path.R1", "absolute_path_R1", "R1", "r1", "fastq_R1", "fq1"),
  label = "R1 FASTQ path"
)

r2_col <- find_col(
  cmd_df,
  candidates = c("absolute_path.R2", "absolute_path_R2", "R2", "r2", "fastq_R2", "fq2"),
  label = "R2 FASTQ path"
)

con1_idx <- which(cmd_df[[sample_col]] == "CON1")
if (length(con1_idx) == 0) {
  stop("command sheet 中找不到 sample_id == CON1。请先检查 stepH1_fastp_command_sheet_18samples.csv。", call. = FALSE)
}
con1_row <- cmd_df[con1_idx[1], , drop = FALSE]

sample_id <- con1_row[[sample_col]]
r1_win <- con1_row[[r1_col]]
r2_win <- con1_row[[r2_col]]

stop_if_missing(r1_win, "CON1 R1 FASTQ")
stop_if_missing(r2_win, "CON1 R2 FASTQ")

# smoke test 输出文件：单独放在 smoke_test_CON1，不覆盖后续正式 trimmed 输出
out_r1_win <- file.path(smoke_dir, paste0(sample_id, ".smoke.R1.trimmed.fastq.gz"))
out_r2_win <- file.path(smoke_dir, paste0(sample_id, ".smoke.R2.trimmed.fastq.gz"))
html_win   <- file.path(smoke_dir, paste0(sample_id, ".smoke.fastp.html"))
json_win   <- file.path(smoke_dir, paste0(sample_id, ".smoke.fastp.json"))

r1_wsl      <- windows_to_wsl(r1_win)
r2_wsl      <- windows_to_wsl(r2_win)
out_r1_wsl  <- windows_to_wsl(out_r1_win)
out_r2_wsl  <- windows_to_wsl(out_r2_win)
html_wsl    <- windows_to_wsl(html_win)
json_wsl    <- windows_to_wsl(json_win)
smoke_log_wsl <- windows_to_wsl(smoke_log)

cat("已读取 H1 command sheet。\n")
cat("本次 smoke test 样本：", sample_id, "\n")
cat("R1: ", r1_win, "\n", sep = "")
cat("R2: ", r2_win, "\n\n", sep = "")

# ----------------------------
# 3. 检查 wsl.exe 是否存在
# ----------------------------

wsl_bin <- Sys.which("wsl.exe")
if (identical(unname(wsl_bin), "")) {
  wsl_bin <- Sys.which("wsl")
}
if (identical(unname(wsl_bin), "")) {
  stop("Windows R 找不到 wsl.exe。请先确认 PowerShell/CMD 中可以运行 wsl。", call. = FALSE)
}
wsl_bin <- unname(wsl_bin)

cat("检测到 WSL 调用程序：", wsl_bin, "\n\n", sep = "")

# ----------------------------
# 4. WSL 内检查 conda fastp_env 和 fastp
# ----------------------------

version_cmd <- paste(
  'echo "===== WSL / fastp preflight check =====";',
  'echo "WSL user: $(whoami)";',
  'echo "WSL home: $HOME";',
  conda_block,
  'conda activate fastp_env;',
  'echo "Conda env: $CONDA_DEFAULT_ENV";',
  'echo "fastp path:";',
  'which fastp;',
  'echo "fastp version:";',
  'fastp --version',
  sep = " "
)

cat("开始检查 WSL 中的 fastp_env / fastp ...\n")
version_res <- run_wsl_bash(version_cmd, version_log, wsl_bin)
which_fastp_exit_status <- version_res$status
fastp_check_ok <- identical(which_fastp_exit_status, 0L)

cat("fastp 检查 exit status：", which_fastp_exit_status, "\n", sep = "")
cat("fastp 检查日志：", version_log, "\n\n", sep = "")

# ----------------------------
# 5. 如果 fastp 可用，则对 CON1 跑 smoke test
# ----------------------------

smoke_exit_status <- NA_integer_
smoke_ok <- FALSE

if (fastp_check_ok) {
  smoke_cmd <- paste(
    'echo "===== CON1 fastp smoke test =====";',
    conda_block,
    'conda activate fastp_env;',
    sprintf('mkdir -p %s;', sh_quote_wsl(dirname(out_r1_wsl))),
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
        '--thread 4',
        '--reads_to_process 100000'
      ),
      sh_quote_wsl(r1_wsl),
      sh_quote_wsl(r2_wsl),
      sh_quote_wsl(out_r1_wsl),
      sh_quote_wsl(out_r2_wsl),
      sh_quote_wsl(html_wsl),
      sh_quote_wsl(json_wsl)
    ),
    sep = " "
  )

  cat("开始运行 CON1 smoke test（只处理前 100000 reads）...\n")
  smoke_res <- run_wsl_bash(smoke_cmd, smoke_log, wsl_bin)
  smoke_exit_status <- smoke_res$status

  file_check_tmp <- data.frame(
    file_role = c("smoke_trimmed_R1", "smoke_trimmed_R2", "smoke_fastp_html", "smoke_fastp_json", "smoke_log"),
    path = c(out_r1_win, out_r2_win, html_win, json_win, smoke_log),
    exists = c(
      file.exists(out_r1_win),
      file.exists(out_r2_win),
      file.exists(html_win),
      file.exists(json_win),
      file.exists(smoke_log)
    ),
    size_bytes = c(
      safe_file_size(out_r1_win),
      safe_file_size(out_r2_win),
      safe_file_size(html_win),
      safe_file_size(json_win),
      safe_file_size(smoke_log)
    ),
    stringsAsFactors = FALSE
  )

  smoke_ok <- identical(smoke_exit_status, 0L) &&
    all(file_check_tmp$exists[1:4]) &&
    all(file_check_tmp$size_bytes[1:4] > 0, na.rm = FALSE)

  cat("smoke test exit status：", smoke_exit_status, "\n", sep = "")
  cat("smoke test 是否成功：", smoke_ok, "\n\n", sep = "")

} else {
  cat("fastp preflight 未通过，因此跳过 CON1 smoke test。\n\n")

  file_check_tmp <- data.frame(
    file_role = c("smoke_trimmed_R1", "smoke_trimmed_R2", "smoke_fastp_html", "smoke_fastp_json", "smoke_log"),
    path = c(out_r1_win, out_r2_win, html_win, json_win, smoke_log),
    exists = c(
      file.exists(out_r1_win),
      file.exists(out_r2_win),
      file.exists(html_win),
      file.exists(json_win),
      file.exists(smoke_log)
    ),
    size_bytes = c(
      safe_file_size(out_r1_win),
      safe_file_size(out_r2_win),
      safe_file_size(html_win),
      safe_file_size(json_win),
      safe_file_size(smoke_log)
    ),
    stringsAsFactors = FALSE
  )
}

# ----------------------------
# 6. 输出 summary 和 file check
# ----------------------------

summary_df <- data.frame(
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  project_root = project_root,
  branch_dir = branch_dir,
  command_sheet = command_sheet,
  sample_id = sample_id,
  r1_windows = r1_win,
  r2_windows = r2_win,
  r1_wsl = r1_wsl,
  r2_wsl = r2_wsl,
  wsl_bin = wsl_bin,
  which_fastp_exit_status = which_fastp_exit_status,
  fastp_check_ok = fastp_check_ok,
  smoke_exit_status = smoke_exit_status,
  smoke_ok = smoke_ok,
  smoke_reads_to_process = 100000,
  out_r1_windows = out_r1_win,
  out_r2_windows = out_r2_win,
  fastp_html_windows = html_win,
  fastp_json_windows = json_win,
  version_log = version_log,
  smoke_log = smoke_log,
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(file_check_tmp, file_check_csv, row.names = FALSE, fileEncoding = "UTF-8")

cat("===== Step H2w1 输出完成 =====\n")
cat("summary:    ", summary_csv, "\n", sep = "")
cat("file check: ", file_check_csv, "\n", sep = "")
cat("version log:", version_log, "\n", sep = "")
cat("smoke log:  ", smoke_log, "\n\n", sep = "")

cat("关键判据：\n")
cat("  which_fastp_exit_status = ", which_fastp_exit_status, "\n", sep = "")
cat("  smoke_ok = ", smoke_ok, "\n\n", sep = "")

if (isTRUE(smoke_ok)) {
  cat("结果：H2w1 成功。下一步可以继续 H2w2：18 个样本批量 fastp trimming。\n")
} else {
  cat("结果：H2w1 未完全成功。请把 stepH2w1_wsl_fastp_smoke_test_summary.csv 和 smoke/version log 的报错内容发给我，我再只针对这个错误修。\n")
}

cat("\n===== End of Step H2w1 =====\n")
