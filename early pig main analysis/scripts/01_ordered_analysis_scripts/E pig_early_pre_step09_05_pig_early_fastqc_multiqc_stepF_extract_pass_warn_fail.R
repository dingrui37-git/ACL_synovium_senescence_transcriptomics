# =========================================================
# pig_early_fastqc_multiqc_stepF_extract_pass_warn_fail.R
# 这一步是干什么：
# 1) 从 pig early 已完成的 FastQC 结果中，提取每个样本、每个模块的 PASS / WARN / FAIL
# 2) 汇总生成模块级统计表（每个模块有多少 PASS / WARN / FAIL）
# 3) 生成每个文件的 FAIL/WARN 数量表，便于快速判断哪些样本问题更多
# 4) 检查 MultiQC 报告和 multiqc_data 是否存在，并记录到 summary
#
# 说明：
# - 本脚本优先解析 fastqc_raw 目录下的 *_fastqc.zip（最稳）
# - 不依赖额外 R 包，直接可跑
# - 运行后把控制台输出复制给我即可
# =========================================================

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
pig_dir <- file.path(project_root, "rebuild_submission", "02_pig_early")
qc_dir <- file.path(pig_dir, "qc")
fastqc_dir <- file.path(qc_dir, "fastqc_raw")
multiqc_dir <- file.path(qc_dir, "multiqc_raw")
tables_dir <- file.path(pig_dir, "tables")
logs_dir <- file.path(pig_dir, "logs")
objects_dir <- file.path(pig_dir, "objects")

for (d in c(tables_dir, logs_dir, objects_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(logs_dir, "stepF_extract_fastqc_multiqc_pass_warn_fail_log.txt")
if (file.exists(log_file)) file.remove(log_file)

append_log <- function(...) {
  msg <- paste0(...)
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

append_log("Step F started.")
append_log("project_root = ", project_root)
append_log("fastqc_dir = ", fastqc_dir)
append_log("multiqc_dir = ", multiqc_dir)

if (!dir.exists(fastqc_dir)) {
  stop("找不到 fastqc_dir：", fastqc_dir)
}

fastqc_zip <- list.files(fastqc_dir, pattern = "_fastqc\\.zip$", full.names = TRUE)
fastqc_html <- list.files(fastqc_dir, pattern = "_fastqc\\.html$", full.names = TRUE)

if (length(fastqc_zip) == 0) {
  stop("fastqc_raw 下没有找到 *_fastqc.zip：", fastqc_dir)
}

append_log("n_fastqc_zip = ", length(fastqc_zip))
append_log("n_fastqc_html = ", length(fastqc_html))

extract_summary_from_zip <- function(zip_path) {
  lst <- utils::unzip(zip_path, list = TRUE)
  target <- lst$Name[grepl("(^|/)summary\\.txt$", lst$Name)]
  if (length(target) == 0) {
    return(data.frame(
      zip_file = basename(zip_path),
      sample_id = sub("_fastqc\\.zip$", "", basename(zip_path), perl = TRUE),
      module = NA_character_,
      status = NA_character_,
      source_filename = NA_character_,
      parse_ok = FALSE,
      parse_note = "summary.txt not found in zip",
      stringsAsFactors = FALSE
    ))
  }

  target <- target[1]
  tmpdir <- tempfile(pattern = "fastqc_unzip_")
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  utils::unzip(zip_path, files = target, exdir = tmpdir)
  extracted_path <- file.path(tmpdir, target)

  if (!file.exists(extracted_path)) {
    return(data.frame(
      zip_file = basename(zip_path),
      sample_id = sub("_fastqc\\.zip$", "", basename(zip_path), perl = TRUE),
      module = NA_character_,
      status = NA_character_,
      source_filename = NA_character_,
      parse_ok = FALSE,
      parse_note = "summary.txt failed to extract",
      stringsAsFactors = FALSE
    ))
  }

  x <- read.delim(extracted_path, header = FALSE, sep = "\t", quote = "", comment.char = "", fill = TRUE)
  if (ncol(x) < 3) {
    return(data.frame(
      zip_file = basename(zip_path),
      sample_id = sub("_fastqc\\.zip$", "", basename(zip_path), perl = TRUE),
      module = NA_character_,
      status = NA_character_,
      source_filename = NA_character_,
      parse_ok = FALSE,
      parse_note = "summary.txt has <3 columns",
      stringsAsFactors = FALSE
    ))
  }

  colnames(x)[1:3] <- c("status", "module", "source_filename")
  x <- x[, c("status", "module", "source_filename"), drop = FALSE]
  x$zip_file <- basename(zip_path)
  x$sample_id <- sub("_fastqc\\.zip$", "", basename(zip_path), perl = TRUE)
  x$parse_ok <- TRUE
  x$parse_note <- "OK"
  x[, c("zip_file", "sample_id", "module", "status", "source_filename", "parse_ok", "parse_note")]
}

summary_list <- lapply(fastqc_zip, extract_summary_from_zip)
fastqc_long <- do.call(rbind, summary_list)

# 排序
ord_sample <- order(fastqc_long$sample_id, fastqc_long$module)
fastqc_long <- fastqc_long[ord_sample, , drop = FALSE]
rownames(fastqc_long) <- NULL

# 只保留解析成功且模块非 NA 的记录做统计
fastqc_long_ok <- fastqc_long[fastqc_long$parse_ok %in% TRUE & !is.na(fastqc_long$module), , drop = FALSE]

if (nrow(fastqc_long_ok) == 0) {
  stop("FastQC zip 找到了，但 summary.txt 解析后没有有效模块记录。")
}

# 模块级 PASS/WARN/FAIL 统计
status_levels <- c("PASS", "WARN", "FAIL")
module_status_tab <- as.data.frame.matrix(table(fastqc_long_ok$module, factor(fastqc_long_ok$status, levels = status_levels)))
for (nm in status_levels) {
  if (!nm %in% colnames(module_status_tab)) module_status_tab[[nm]] <- 0L
}
module_status_tab$module <- rownames(module_status_tab)
rownames(module_status_tab) <- NULL
module_status_tab <- module_status_tab[, c("module", status_levels), drop = FALSE]
module_status_tab$n_total <- rowSums(module_status_tab[, status_levels, drop = FALSE])
module_status_tab$fail_rate_pct <- round(ifelse(module_status_tab$n_total > 0, module_status_tab$FAIL / module_status_tab$n_total * 100, NA_real_), 2)
module_status_tab$warn_rate_pct <- round(ifelse(module_status_tab$n_total > 0, module_status_tab$WARN / module_status_tab$n_total * 100, NA_real_), 2)
module_status_tab <- module_status_tab[order(-module_status_tab$FAIL, -module_status_tab$WARN, module_status_tab$module), , drop = FALSE]

# 文件级别统计
file_status_tab <- as.data.frame.matrix(table(fastqc_long_ok$sample_id, factor(fastqc_long_ok$status, levels = status_levels)))
for (nm in status_levels) {
  if (!nm %in% colnames(file_status_tab)) file_status_tab[[nm]] <- 0L
}
file_status_tab$sample_id <- rownames(file_status_tab)
rownames(file_status_tab) <- NULL
file_status_tab <- file_status_tab[, c("sample_id", status_levels), drop = FALSE]
file_status_tab$n_total_modules <- rowSums(file_status_tab[, status_levels, drop = FALSE])
file_status_tab <- file_status_tab[order(-file_status_tab$FAIL, -file_status_tab$WARN, file_status_tab$sample_id), , drop = FALSE]

# 宽表：每个 sample_id 一行，每个 module 一列
modules_all <- sort(unique(fastqc_long_ok$module))
samples_all <- sort(unique(fastqc_long_ok$sample_id))
fastqc_wide <- data.frame(sample_id = samples_all, stringsAsFactors = FALSE)
for (m in modules_all) {
  tmp <- fastqc_long_ok[fastqc_long_ok$module == m, c("sample_id", "status"), drop = FALSE]
  names(tmp)[2] <- m
  fastqc_wide <- merge(fastqc_wide, tmp, by = "sample_id", all.x = TRUE, sort = FALSE)
}
fastqc_wide <- fastqc_wide[match(samples_all, fastqc_wide$sample_id), , drop = FALSE]

# MultiQC 检查
multiqc_report_html <- file.path(multiqc_dir, "multiqc_report.html")
multiqc_data_dir <- file.path(multiqc_dir, "multiqc_data")
multiqc_files <- if (dir.exists(multiqc_dir)) list.files(multiqc_dir, recursive = TRUE, full.names = FALSE) else character(0)
multiqc_data_files <- if (dir.exists(multiqc_data_dir)) list.files(multiqc_data_dir, recursive = TRUE, full.names = FALSE) else character(0)

summary_df <- data.frame(
  metric = c(
    "project_root",
    "fastqc_dir",
    "multiqc_dir",
    "n_fastqc_zip",
    "n_fastqc_html",
    "n_fastqc_samples_parsed",
    "n_fastqc_module_records",
    "n_unique_modules",
    "multiqc_report_exists",
    "multiqc_data_dir_exists",
    "n_files_under_multiqc_dir",
    "n_files_under_multiqc_data",
    "status"
  ),
  value = c(
    project_root,
    fastqc_dir,
    multiqc_dir,
    length(fastqc_zip),
    length(fastqc_html),
    length(unique(fastqc_long_ok$sample_id)),
    nrow(fastqc_long_ok),
    length(unique(fastqc_long_ok$module)),
    file.exists(multiqc_report_html),
    dir.exists(multiqc_data_dir),
    length(multiqc_files),
    length(multiqc_data_files),
    "Step F completed successfully"
  ),
  stringsAsFactors = FALSE
)

# 输出文件
long_file <- file.path(tables_dir, "stepF_fastqc_module_status_long.csv")
wide_file <- file.path(tables_dir, "stepF_fastqc_module_status_wide.csv")
module_file <- file.path(tables_dir, "stepF_fastqc_module_status_summary.csv")
file_file <- file.path(tables_dir, "stepF_fastqc_file_status_summary.csv")
summary_file <- file.path(tables_dir, "stepF_fastqc_multiqc_pass_warn_fail_summary.csv")

write.csv(fastqc_long, long_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fastqc_wide, wide_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(module_status_tab, module_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(file_status_tab, file_file, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(summary_df, summary_file, row.names = FALSE, fileEncoding = "UTF-8")

save(
  fastqc_long,
  fastqc_long_ok,
  fastqc_wide,
  module_status_tab,
  file_status_tab,
  summary_df,
  file = file.path(objects_dir, "stepF_fastqc_multiqc_pass_warn_fail_workspace.RData")
)

append_log("Long table saved: ", long_file)
append_log("Wide table saved: ", wide_file)
append_log("Module summary saved: ", module_file)
append_log("File summary saved: ", file_file)
append_log("Summary saved: ", summary_file)
append_log("Step F finished successfully.")

# 控制台打印
cat("===== Step F summary =====\n")
print(summary_df, row.names = FALSE)

cat("\n===== FastQC module PASS/WARN/FAIL summary =====\n")
print(module_status_tab, row.names = FALSE)

cat("\n===== Top samples by FAIL/WARN count =====\n")
print(utils::head(file_status_tab, 20), row.names = FALSE)

cat("\n输出文件：\n")
cat(summary_file, "\n")
cat(module_file, "\n")
cat(file_file, "\n")
cat(long_file, "\n")
cat(wide_file, "\n")
cat(log_file, "\n")
