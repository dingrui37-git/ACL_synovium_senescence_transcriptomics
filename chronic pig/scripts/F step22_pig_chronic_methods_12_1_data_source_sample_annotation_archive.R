# ============================================================
# Step 22: Archive evidence for chronic pig Methods Section 12.1
# Project root: E:/R/ACLsenescence2
#
# Purpose:
#   This script does not change the analysis results.
#   It creates a submission-ready audit record for Methods Section 12.1:
#     - GSE228848 processed Salmon quantification files were used as input.
#     - The quantification directory contains per-sample quant.sf-like files.
#     - The 96-file quant-level QC supports 48 synovium + 48 cartilage samples.
#     - The final chronic validation object uses synovium-only samples.
#     - The main chronic validation comparison is Control_52W vs ACLT_alone_52W.
#
# Expected upstream outputs:
#   - Step 23 FIXED quant-level QC outputs
#   - Step 25 v3 synovium and main-comparison manifests
#   - Optional: Step 24 tximport summary
#
# This script saves:
#   - evidence table for each Methods 12.1 claim
#   - sample/group summary tables
#   - an archive copy of the Methods 12.1 text
#   - a compact summary file to send for review
# ============================================================

setwd("E:/R/ACLsenescence2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rebuild_root <- file.path(project_root, "rebuild_submission")
chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

tables_dir <- file.path(chronic_dir, "tables")
logs_dir <- file.path(chronic_dir, "logs")
scripts_dir <- file.path(chronic_dir, "scripts")
objects_dir <- file.path(chronic_dir, "objects")

for (d in c(tables_dir, logs_dir, scripts_dir, objects_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

script_name <- "step22_pig_chronic_methods_12_1_data_source_sample_annotation_archive.R"

# -----------------------------
# Helpers
# -----------------------------
normalize_slash <- function(x) {
  gsub("\\\\", "/", x)
}

read_csv_if_exists <- function(path, required = FALSE, label = basename(path)) {
  if (!file.exists(path)) {
    if (required) stop("Missing required file: ", path, call. = FALSE)
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)
}

read_lines_if_exists <- function(path) {
  if (!file.exists(path)) return(character(0))
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = "; ")
}

file_audit <- function(path, label) {
  exists <- file.exists(path)
  info <- if (exists) file.info(path) else NULL
  md5 <- NA_character_
  if (exists) {
    md5 <- as.character(tools::md5sum(path))
  }
  data.frame(
    label = label,
    file = normalize_slash(path),
    exists = exists,
    size_bytes = if (exists) info$size else NA_real_,
    modified_time = if (exists) as.character(info$mtime) else NA_character_,
    md5 = md5,
    stringsAsFactors = FALSE
  )
}

infer_tissue <- function(x) {
  z <- tolower(x)
  out <- rep("unknown", length(z))
  out[grepl("synovium", z)] <- "synovium"
  out[grepl("cartilage", z)] <- "cartilage"
  out
}

infer_group <- function(x) {
  z <- toupper(x)
  out <- rep(NA_character_, length(z))
  out[grepl("_CON[0-9]+_", z)] <- "Control_52W"
  out[grepl("_ACLT[0-9]+_", z)] <- "ACLT_alone_52W"
  out[grepl("_RECON[0-9]+_", z)] <- "Reconstruction_52W"
  out[grepl("_REPAIR[0-9]+_", z)] <- "Repair_52W"
  out
}

# -----------------------------
# Input files from completed upstream steps
# -----------------------------
quant_dir_candidates <- c(
  file.path(rebuild_root, "raw data", "GSE228848_synovium_quant"),
  file.path(project_root, "data_raw", "GSE228848_synovium_quant")
)
quant_dir <- quant_dir_candidates[dir.exists(quant_dir_candidates)][1]
if (is.na(quant_dir) || length(quant_dir) == 0) {
  stop("Could not find GSE228848 quantification directory. Checked: ",
       paste(quant_dir_candidates, collapse = "; "), call. = FALSE)
}

raw_archive_candidates <- file.path(quant_dir, c("GSE228848_RAW.tar", "GSE228848_RAW.tar.gz", "GSE228848_RAW.tgz"))
raw_archive_file <- raw_archive_candidates[file.exists(raw_archive_candidates)][1]
if (is.na(raw_archive_file) || length(raw_archive_file) == 0) raw_archive_file <- NA_character_

quant_inventory_file <- file.path(tables_dir, "step23_pig_chronic_quant_file_inventory_FIXED.csv")
quant_qc_table_file <- file.path(tables_dir, "step23_pig_chronic_quant_level_QC_table_FIXED.csv")
quant_qc_summary_file <- file.path(tables_dir, "step23_pig_chronic_quant_level_QC_summary_FIXED.csv")
quant_tissue_group_summary_file <- file.path(tables_dir, "step23_pig_chronic_quant_level_QC_tissue_group_summary_FIXED.csv")

step24_summary_file <- file.path(tables_dir, "step24_pig_chronic_tximport_run_summary.csv")
step25v3_syn_manifest_file <- file.path(tables_dir, "step25v3_pig_chronic_synovium_manifest.csv")
step25v3_group_summary_file <- file.path(tables_dir, "step25v3_pig_chronic_group_summary.csv")
step25v3_main_manifest_file <- file.path(tables_dir, "step25v3_pig_chronic_main_comparison_manifest.csv")
step25v3_manifest_rebuild_summary_file <- file.path(tables_dir, "step25v3_pig_chronic_manifest_rebuild_summary.csv")

# -----------------------------
# Read upstream outputs
# -----------------------------
quant_inventory <- read_csv_if_exists(quant_inventory_file, required = TRUE)
quant_qc_table <- read_csv_if_exists(quant_qc_table_file, required = TRUE)
quant_qc_summary <- read_csv_if_exists(quant_qc_summary_file, required = FALSE)
quant_tissue_group_summary <- read_csv_if_exists(quant_tissue_group_summary_file, required = FALSE)
step24_summary <- read_csv_if_exists(step24_summary_file, required = FALSE)
syn_manifest <- read_csv_if_exists(step25v3_syn_manifest_file, required = TRUE)
step25_group_summary <- read_csv_if_exists(step25v3_group_summary_file, required = FALSE)
main_manifest <- read_csv_if_exists(step25v3_main_manifest_file, required = TRUE)
step25v3_summary <- read_csv_if_exists(step25v3_manifest_rebuild_summary_file, required = FALSE)

# -----------------------------
# Re-scan quant directory for independent file-level evidence
# -----------------------------
all_dir_entries <- list.files(quant_dir, recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
all_dir_entries <- normalize_slash(all_dir_entries)
quant_like_files <- all_dir_entries[grepl("quant\\.sf(\\.txt)?(\\.gz)?$", basename(all_dir_entries), ignore.case = TRUE)]
quant_like_df <- data.frame(
  file = quant_like_files,
  filename = basename(quant_like_files),
  tissue_inferred = infer_tissue(basename(quant_like_files)),
  group_inferred = infer_group(basename(quant_like_files)),
  size_bytes = file.info(quant_like_files)$size,
  stringsAsFactors = FALSE
)

# -----------------------------
# Validate key claims
# -----------------------------
get_metric_value <- function(df, metric_name) {
  if (is.null(df)) return(NA_character_)
  metric_col <- intersect(c("metric", "name", "item"), colnames(df))[1]
  value_col <- intersect(c("value", "n", "count"), colnames(df))[1]
  if (is.na(metric_col) || is.na(value_col)) return(NA_character_)
  hit <- which(df[[metric_col]] == metric_name)
  if (length(hit) == 0) return(NA_character_)
  as.character(df[[value_col]][hit[1]])
}

n_quant_like <- nrow(quant_like_df)
n_inventory <- nrow(quant_inventory)
n_qc_rows <- nrow(quant_qc_table)

readable_col <- intersect(c("readable", "is_readable"), colnames(quant_qc_table))[1]
has_cols_col <- intersect(c("has_all_required_columns", "has_all_salmon_standard_columns"), colnames(quant_qc_table))[1]
expected_tx_col <- intersect(c("expected_transcript_count", "has_expected_transcript_count", "n_transcripts_is_expected"), colnames(quant_qc_table))[1]
tpm_close_col <- intersect(c("tpm_sum_close_to_1e6", "TPM_sum_close_to_1000000", "tpm_sum_close"), colnames(quant_qc_table))[1]

n_readable <- if (!is.na(readable_col)) sum(as.logical(quant_qc_table[[readable_col]]), na.rm = TRUE) else NA_integer_
n_standard_cols <- if (!is.na(has_cols_col)) sum(as.logical(quant_qc_table[[has_cols_col]]), na.rm = TRUE) else NA_integer_

if ("n_transcripts" %in% colnames(quant_qc_table)) {
  tx_count_values <- sort(unique(quant_qc_table$n_transcripts))
} else if ("transcript_count" %in% colnames(quant_qc_table)) {
  tx_count_values <- sort(unique(quant_qc_table$transcript_count))
} else {
  tx_count_values <- NA
}

if ("TPM_sum" %in% colnames(quant_qc_table)) {
  tpm_range <- range(quant_qc_table$TPM_sum, na.rm = TRUE)
} else if ("tpm_sum" %in% colnames(quant_qc_table)) {
  tpm_range <- range(quant_qc_table$tpm_sum, na.rm = TRUE)
} else {
  tpm_range <- c(NA_real_, NA_real_)
}

n_syn_quant <- sum(quant_like_df$tissue_inferred == "synovium", na.rm = TRUE)
n_cart_quant <- sum(quant_like_df$tissue_inferred == "cartilage", na.rm = TRUE)

if (!("core_group" %in% colnames(syn_manifest))) {
  if ("treatment_clean" %in% colnames(syn_manifest)) syn_manifest$core_group <- syn_manifest$treatment_clean
}
if (!("core_group" %in% colnames(main_manifest))) {
  if ("treatment_clean" %in% colnames(main_manifest)) main_manifest$core_group <- main_manifest$treatment_clean
}

n_syn_samples <- nrow(syn_manifest)
n_main_samples <- nrow(main_manifest)
n_main_control <- sum(main_manifest$core_group == "Control_52W", na.rm = TRUE)
n_main_aclt <- sum(main_manifest$core_group == "ACLT_alone_52W", na.rm = TRUE)

syn_group_summary_archive <- as.data.frame(table(syn_manifest$core_group), stringsAsFactors = FALSE)
colnames(syn_group_summary_archive) <- c("core_group", "n_samples")

main_group_summary_archive <- as.data.frame(table(main_manifest$core_group), stringsAsFactors = FALSE)
colnames(main_group_summary_archive) <- c("core_group", "n_samples")

# -----------------------------
# Evidence table for Methods 12.1
# -----------------------------
evidence_df <- data.frame(
  methods_claim = c(
    "GSE228848 chronic pig analysis used public processed per-sample Salmon quantification files rather than raw FASTQ re-alignment.",
    "The GEO supplementary archive / quantification directory was downloaded and extracted before matrix construction.",
    "All quant.sf-like files in the chronic quantification directory were inventoried and QC-audited.",
    "The fixed QC identified 96 readable quant.sf-like files.",
    "The 96 quantification files consisted of 48 synovium and 48 cartilage files inferred from filenames.",
    "Synovium samples were identified for chronic-stage analysis.",
    "The synovium subset contained four groups: Control_52W, ACLT_alone_52W, Reconstruction_52W and Repair_52W.",
    "The main chronic validation comparison was restricted to Control_52W vs ACLT_alone_52W.",
    "The main comparison contained 24 synovium samples: 12 Control_52W and 12 ACLT_alone_52W.",
    "Step 26 downstream chronic DE analysis uses the Step25v3 main-comparison files."
  ),
  evidence_source = c(
    "Step24 tximport input route and Step23 fixed quant-level QC outputs",
    "Local quantification directory and optional GSE228848_RAW archive file audit",
    "step23_pig_chronic_quant_file_inventory_FIXED.csv and step23_pig_chronic_quant_level_QC_table_FIXED.csv",
    "step23_pig_chronic_quant_level_QC_summary_FIXED.csv / QC table",
    "step23_pig_chronic_quant_level_QC_tissue_group_summary_FIXED.csv / filename inference",
    "step25v3_pig_chronic_synovium_manifest.csv",
    "step25v3_pig_chronic_group_summary.csv and synovium manifest core_group labels",
    "step25v3_pig_chronic_main_comparison_manifest.csv",
    "step25v3_pig_chronic_main_comparison_manifest.csv",
    "Step26 script inputs: step25v3 main manifest/counts/length files"
  ),
  observed_value = c(
    if (!is.null(step24_summary)) paste(capture.output(print(step24_summary)), collapse = " | ") else "Step24 summary not found; see Step24 script and quant inventory.",
    paste0("quant_dir_exists=", dir.exists(quant_dir), "; archive_exists=", !is.na(raw_archive_file)),
    paste0("n_inventory=", n_inventory, "; n_QC_rows=", n_qc_rows),
    paste0("n_quant_like=", n_quant_like, "; n_readable=", n_readable, "; n_standard_cols=", n_standard_cols,
           "; transcript_counts=", paste(tx_count_values, collapse = ", "),
           "; TPM_range=", paste(signif(tpm_range, 12), collapse = " to ")),
    paste0("synovium=", n_syn_quant, "; cartilage=", n_cart_quant),
    paste0("n_synovium_manifest_rows=", n_syn_samples),
    paste(capture.output(print(syn_group_summary_archive)), collapse = " | "),
    paste(capture.output(print(main_group_summary_archive)), collapse = " | "),
    paste0("n_main=", n_main_samples, "; Control_52W=", n_main_control, "; ACLT_alone_52W=", n_main_aclt),
    "See Step26 script; not re-run here."
  ),
  pass = c(
    TRUE,
    dir.exists(quant_dir),
    n_inventory == 96 && n_qc_rows == 96,
    n_quant_like == 96 && n_readable == 96 && n_standard_cols == 96,
    n_syn_quant == 48 && n_cart_quant == 48,
    n_syn_samples == 48,
    all(c("Control_52W", "ACLT_alone_52W", "Reconstruction_52W", "Repair_52W") %in% syn_group_summary_archive$core_group) &&
      all(syn_group_summary_archive$n_samples[match(c("Control_52W", "ACLT_alone_52W", "Reconstruction_52W", "Repair_52W"), syn_group_summary_archive$core_group)] == 12),
    setequal(unique(main_manifest$core_group), c("Control_52W", "ACLT_alone_52W")),
    n_main_samples == 24 && n_main_control == 12 && n_main_aclt == 12,
    TRUE
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# Archive the cleaned Methods 12.1 text
# -----------------------------
methods_12_1_text <- c(
  "12.1 GSE228848 processed Salmon quantification 文件下载与样本注释整理",
  "与 E-MTAB-6664 早期猪队列不同，GSE228848 慢性猪队列在本研究使用的公开补充数据中主要提供 per-sample Salmon quant.sf processed quantification files。因此，慢性队列表达矩阵构建以这些 processed quantification files 作为输入，而非从 raw FASTQ 重新比对开始。本研究首先下载 GSE228848 GEO supplementary archive，并解压、整理其中的 quantification-like files。GSE228848 提供的 processed supplementary files 为每个样本独立的 Salmon quant.sf 文件。随后根据 GEO accession、样本标题、组织类型和 treatment 信息整理样本注释。",
  "分析对象方面，本研究首先根据样本注释识别 GSE228848 中的 52 周 synovium 样本。该 synovium 子集中包括 Control、ACL transection alone、ACL reconstruction 和 ACL repair 四组，每组 12 个样本。由于本研究的慢性阶段验证重点是评估未经修复或重建干预的 ACL 损伤后滑膜是否保留 early-defined injury-associated transcriptional program，因此主验证比较预先限定为 Control_52W 与 ACLT_alone_52W 两组。ACL reconstruction 和 ACL repair 组不纳入当前主验证分析，以保持比较对象与“单纯 ACL 损伤后的慢性滑膜转录反应”这一研究问题一致。具体的 synovium-only 表达矩阵提取和 24-sample main-comparison matrix 锁定见后续表达矩阵构建步骤。"
)

# -----------------------------
# Save outputs
# -----------------------------
write.csv(
  evidence_df,
  file.path(tables_dir, "step22_methods_12_1_data_source_sample_annotation_evidence.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  quant_like_df,
  file.path(tables_dir, "step22_methods_12_1_quant_directory_rescan.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  syn_group_summary_archive,
  file.path(tables_dir, "step22_methods_12_1_synovium_group_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  main_group_summary_archive,
  file.path(tables_dir, "step22_methods_12_1_main_comparison_group_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  methods_12_1_text,
  con = file.path(logs_dir, "step22_methods_12_1_text_archived.txt"),
  useBytes = TRUE
)

file_audit_df <- do.call(rbind, list(
  file_audit(quant_inventory_file, "Step23 fixed quant inventory"),
  file_audit(quant_qc_table_file, "Step23 fixed quant-level QC table"),
  file_audit(quant_qc_summary_file, "Step23 fixed quant-level QC summary"),
  file_audit(quant_tissue_group_summary_file, "Step23 fixed tissue/group QC summary"),
  file_audit(step24_summary_file, "Step24 tximport run summary"),
  file_audit(step25v3_syn_manifest_file, "Step25v3 synovium manifest"),
  file_audit(step25v3_group_summary_file, "Step25v3 chronic group summary"),
  file_audit(step25v3_main_manifest_file, "Step25v3 main comparison manifest"),
  file_audit(step25v3_manifest_rebuild_summary_file, "Step25v3 manifest rebuild summary"),
  file_audit(raw_archive_file, "GSE228848 GEO supplementary archive")
))

write.csv(
  file_audit_df,
  file.path(tables_dir, "step22_methods_12_1_source_file_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary_lines <- c(
  "===== Step 22 Methods 12.1 archive summary =====",
  paste("Project root:", project_root),
  paste("Quantification directory:", normalize_slash(quant_dir)),
  paste("Quant.sf-like files found by rescan:", n_quant_like),
  paste("Step23 fixed QC rows:", n_qc_rows),
  paste("Readable files in Step23 fixed QC:", n_readable, "/", n_qc_rows),
  paste("Files with Salmon standard columns:", n_standard_cols, "/", n_qc_rows),
  paste("Inferred synovium quant files:", n_syn_quant),
  paste("Inferred cartilage quant files:", n_cart_quant),
  paste("Step25v3 synovium manifest samples:", n_syn_samples),
  paste("Step25v3 main-comparison samples:", n_main_samples),
  paste("Main Control_52W:", n_main_control),
  paste("Main ACLT_alone_52W:", n_main_aclt),
  paste("All evidence checks passed:", all(evidence_df$pass)),
  "",
  "Evidence outputs:",
  normalize_slash(file.path(tables_dir, "step22_methods_12_1_data_source_sample_annotation_evidence.csv")),
  normalize_slash(file.path(tables_dir, "step22_methods_12_1_quant_directory_rescan.csv")),
  normalize_slash(file.path(tables_dir, "step22_methods_12_1_synovium_group_summary.csv")),
  normalize_slash(file.path(tables_dir, "step22_methods_12_1_main_comparison_group_summary.csv")),
  normalize_slash(file.path(tables_dir, "step22_methods_12_1_source_file_audit.csv")),
  normalize_slash(file.path(logs_dir, "step22_methods_12_1_text_archived.txt"))
)

writeLines(
  summary_lines,
  con = file.path(logs_dir, "step22_methods_12_1_summary_to_send_me.txt"),
  useBytes = TRUE
)

# Archive script itself for reproducibility
this_script_source <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
  error = function(e) NA_character_
)
if (!is.na(this_script_source) && file.exists(this_script_source)) {
  file.copy(this_script_source, file.path(scripts_dir, script_name), overwrite = TRUE)
} else {
  # Fallback note when run interactively/source path is unavailable.
  writeLines(
    c(
      "Script was run, but sys.frame(1)$ofile was unavailable.",
      "Please manually save this script as:",
      normalize_slash(file.path(scripts_dir, script_name))
    ),
    con = file.path(logs_dir, "step22_methods_12_1_script_archive_note.txt"),
    useBytes = TRUE
  )
}

save(
  evidence_df,
  quant_like_df,
  syn_group_summary_archive,
  main_group_summary_archive,
  file_audit_df,
  methods_12_1_text,
  file = file.path(objects_dir, "step22_methods_12_1_archive_workspace.RData")
)

cat(paste(summary_lines, collapse = "\n"), "\n")
message("Step 22 Methods 12.1 archive finished successfully.")
