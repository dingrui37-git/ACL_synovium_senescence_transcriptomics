# ============================================================
# Step 25 v2: rebuild pig-chronic manifests and subset matrices
# Project root: E:/R/ACLsenescence2
#
# What this fixes:
#   - Step 25 v1 incorrectly inferred "synovium" from free-text rows
#     in GEO metadata, which also mention synovium in protocol text.
#   - That caused cartilage rows to be mis-labeled as synovium.
#
# Strategy in v2:
#   1) Prefer the existing local 48-row synovium manifest files that
#      were already produced in the earlier successful chronic workflow.
#   2) Fall back to parsed sample-info only if needed.
#   3) Rebuild:
#        - 48-sample synovium manifest + matrices
#        - 24-sample main comparison manifest
#          (Control_52W vs ACLT_alone_52W, synovium only)
#   4) Auto-save all files for submission/reproducibility.
# ============================================================

setwd("E:/R/ACLsenescence2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rebuild_root <- file.path(project_root, "rebuild_submission")
chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

scripts_dir <- file.path(chronic_dir, "scripts")
objects_dir <- file.path(chronic_dir, "objects")
tables_dir <- file.path(chronic_dir, "tables")
logs_dir <- file.path(chronic_dir, "logs")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8"),
    error = function(e1) {
      tryCatch(
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e2) NULL
      )
    }
  )
}

extract_gsm <- function(x) {
  m <- regexpr("GSM[0-9]+", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  has_match <- m > 0
  out[has_match] <- regmatches(x, m)[has_match]
  out
}

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

assign_core_group <- function(x) {
  y <- tolower(normalize_text(x))
  out <- rep(NA_character_, length(y))
  out[grepl("control", y) & grepl("52", y)] <- "Control_52W"
  out[grepl("aclt|acl transection alone", y) & grepl("52", y)] <- "ACLT_alone_52W"
  out[grepl("reconstruction|recon", y) & grepl("52", y)] <- "Reconstruction_52W"
  out[grepl("repair", y) & grepl("52", y)] <- "Repair_52W"
  out
}

pick_first_existing <- function(paths) {
  idx <- which(file.exists(paths))
  if (length(idx) == 0) return(NA_character_)
  paths[idx[1]]
}

quant_inventory_file <- file.path(tables_dir, "step24_pig_chronic_quant_inventory.csv")
counts_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_estimated_counts_matrix.csv")
abundance_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_abundance_tpm_matrix.csv")
length_file <- file.path(tables_dir, "step24_pig_chronic_gene_level_length_matrix.csv")

needed <- c(quant_inventory_file, counts_file, abundance_file, length_file)
if (!all(file.exists(needed))) {
  stop("Required Step 24 files are missing. Please complete Step 24 first.")
}

quant_inventory_df <- read.csv(quant_inventory_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
counts_df <- read.csv(counts_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)
abundance_df <- read.csv(abundance_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)
length_df <- read.csv(length_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)

quant_inventory_df$sample_id <- normalize_text(quant_inventory_df$sample_id)
quant_inventory_df$gsm_accession <- extract_gsm(quant_inventory_df$sample_id)

# ------------------------------------------------------------
# Prefer previously generated local synovium manifest files
# ------------------------------------------------------------
preferred_manifest_candidates <- c(
  file.path(project_root, "results", "GSE228848_chronic_validation", "GSE228848_synovium_chronic_manifest_repaired.csv"),
  file.path(project_root, "results", "GSE228848_chronic_validation", "GSE228848_synovium_chronic_manifest.csv"),
  file.path(project_root, "results", "GSE228848_chronic_validation", "GSE228848_chronic_manifest_final.csv")
)

preferred_manifest_file <- pick_first_existing(preferred_manifest_candidates)

if (is.na(preferred_manifest_file)) {
  stop(
    "No previously generated chronic synovium manifest file was found. ",
    "Expected one of:\n",
    paste(preferred_manifest_candidates, collapse = "\n")
  )
}

manifest_raw <- safe_read_csv(preferred_manifest_file)
if (is.null(manifest_raw) || nrow(manifest_raw) == 0) {
  stop("Failed to read preferred chronic manifest file: ", preferred_manifest_file)
}

manifest_raw[] <- lapply(manifest_raw, normalize_text)
colnames(manifest_raw) <- normalize_text(colnames(manifest_raw))
cn_lower <- tolower(colnames(manifest_raw))

find_col <- function(patterns) {
  for (pat in patterns) {
    hit <- which(grepl(pat, cn_lower, perl = TRUE))
    if (length(hit) >= 1) return(colnames(manifest_raw)[hit[1]])
  }
  NA_character_
}

sample_col <- find_col(c("^sample_id$", "sample[_ ]?id", "sample", "title"))
gsm_col <- find_col(c("^geo_accession$", "^gsm$", "geo[_ ]?accession", "gsm"))
tissue_col <- find_col(c("^tissue$", "tissue"))
group_col <- find_col(c("^core_group$", "core[_ ]?group", "group", "treatment"))
time_col <- find_col(c("^time$", "time"))
source_name_col <- find_col(c("source[_ ]?name"))

if (is.na(sample_col) && is.na(gsm_col)) {
  stop("Could not identify sample ID / GSM column in: ", preferred_manifest_file)
}

manifest_df <- manifest_raw

manifest_df$gsm_accession <- if (!is.na(gsm_col)) extract_gsm(manifest_df[[gsm_col]]) else NA_character_
if (all(is.na(manifest_df$gsm_accession)) && !is.na(sample_col)) {
  manifest_df$gsm_accession <- extract_gsm(manifest_df[[sample_col]])
}

if (all(is.na(manifest_df$gsm_accession))) {
  stop("Could not recover GSM accession from preferred manifest file: ", preferred_manifest_file)
}

manifest_df$sample_id <- quant_inventory_df$sample_id[match(manifest_df$gsm_accession, quant_inventory_df$gsm_accession)]
if (any(is.na(manifest_df$sample_id))) {
  missing_gsm <- unique(manifest_df$gsm_accession[is.na(manifest_df$sample_id)])
  stop(
    "Some manifest GSM accessions could not be matched to Step 24 quant files: ",
    paste(missing_gsm, collapse = ", ")
  )
}

manifest_df$tissue_label <- if (!is.na(tissue_col)) manifest_df[[tissue_col]] else ""
if (all(manifest_df$tissue_label == "") && !is.na(source_name_col)) {
  manifest_df$tissue_label <- manifest_df[[source_name_col]]
}
if (all(manifest_df$tissue_label == "") && !is.na(sample_col)) {
  manifest_df$tissue_label <- manifest_df[[sample_col]]
}

manifest_df$core_group <- if (!is.na(group_col)) manifest_df[[group_col]] else ""
if (all(manifest_df$core_group == "") && !is.na(sample_col)) {
  manifest_df$core_group <- manifest_df[[sample_col]]
}
manifest_df$core_group <- assign_core_group(manifest_df$core_group)

manifest_df$time_label <- if (!is.na(time_col)) manifest_df[[time_col]] else ""
manifest_df$tissue_synovium <- grepl("synov", tolower(manifest_df$tissue_label))
manifest_df$tissue_cartilage <- grepl("cartilage", tolower(manifest_df$tissue_label))

# Keep only true synovium rows from the preferred manifest
synovium_manifest_df <- manifest_df[manifest_df$tissue_synovium & !manifest_df$tissue_cartilage, , drop = FALSE]

# If tissue column is absent or unhelpful, recover by sample_id suffix
if (nrow(synovium_manifest_df) == 0) {
  synovium_manifest_df <- manifest_df[grepl("_synovium$", manifest_df$sample_id, ignore.case = TRUE), , drop = FALSE]
  synovium_manifest_df$tissue_synovium <- TRUE
  synovium_manifest_df$tissue_cartilage <- FALSE
}

if (nrow(synovium_manifest_df) == 0) {
  stop("Failed to recover synovium-only rows from preferred chronic manifest.")
}

synovium_manifest_df <- synovium_manifest_df[!duplicated(synovium_manifest_df$sample_id), , drop = FALSE]

# Rebuild a clean output manifest with the user-facing fields we need
keep_cols <- c(
  "sample_id", "gsm_accession", "core_group", "tissue_label", "time_label",
  "tissue_synovium", "tissue_cartilage"
)
keep_cols <- keep_cols[keep_cols %in% colnames(synovium_manifest_df)]

synovium_manifest_df <- synovium_manifest_df[, keep_cols, drop = FALSE]
synovium_manifest_df$quant_relative <- quant_inventory_df$relative_path[match(synovium_manifest_df$sample_id, quant_inventory_df$sample_id)]
synovium_manifest_df$quant_exists <- !is.na(synovium_manifest_df$quant_relative)
synovium_manifest_df$source_manifest_file <- sub(paste0("^", gsub("\\\\", "/", project_root), "/?"), "", gsub("\\\\", "/", preferred_manifest_file))

synovium_manifest_df <- synovium_manifest_df[order(synovium_manifest_df$core_group, synovium_manifest_df$sample_id), , drop = FALSE]

# Expected: 48 synovium samples
# Main comparison: synovium only, Control_52W vs ACLT_alone_52W
main_manifest_df <- synovium_manifest_df[synovium_manifest_df$core_group %in% c("Control_52W", "ACLT_alone_52W"), , drop = FALSE]
main_manifest_df <- main_manifest_df[order(main_manifest_df$core_group, main_manifest_df$sample_id), , drop = FALSE]

subset_matrix_df <- function(mat_df, keep_samples) {
  keep_samples <- keep_samples[keep_samples %in% colnames(mat_df)]
  mat_df[, c("gene_id", keep_samples), drop = FALSE]
}

synovium_counts_df <- subset_matrix_df(counts_df, synovium_manifest_df$sample_id)
synovium_abundance_df <- subset_matrix_df(abundance_df, synovium_manifest_df$sample_id)
synovium_length_df <- subset_matrix_df(length_df, synovium_manifest_df$sample_id)

main_counts_df <- subset_matrix_df(counts_df, main_manifest_df$sample_id)
main_abundance_df <- subset_matrix_df(abundance_df, main_manifest_df$sample_id)
main_length_df <- subset_matrix_df(length_df, main_manifest_df$sample_id)

summary_df <- data.frame(
  metric = c(
    "project_root",
    "preferred_manifest_file",
    "n_total_quant_samples",
    "n_synovium_samples",
    "n_main_comparison_samples",
    "n_main_control_52w",
    "n_main_aclt_alone_52w",
    "synovium_counts_dim",
    "main_counts_dim",
    "status"
  ),
  value = c(
    project_root,
    sub(paste0("^", gsub("\\\\", "/", project_root), "/?"), "", gsub("\\\\", "/", preferred_manifest_file)),
    nrow(quant_inventory_df),
    nrow(synovium_manifest_df),
    nrow(main_manifest_df),
    sum(main_manifest_df$core_group == "Control_52W", na.rm = TRUE),
    sum(main_manifest_df$core_group == "ACLT_alone_52W", na.rm = TRUE),
    paste(nrow(synovium_counts_df), "x", ncol(synovium_counts_df) - 1),
    paste(nrow(main_counts_df), "x", ncol(main_counts_df) - 1),
    "Step 25 v2 completed successfully"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  synovium_manifest_df,
  file.path(tables_dir, "step25v2_pig_chronic_synovium_manifest.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  main_manifest_df,
  file.path(tables_dir, "step25v2_pig_chronic_main_comparison_manifest.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  synovium_counts_df,
  file.path(tables_dir, "step25v2_pig_chronic_synovium_gene_level_counts_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  synovium_abundance_df,
  file.path(tables_dir, "step25v2_pig_chronic_synovium_gene_level_abundance_tpm_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  synovium_length_df,
  file.path(tables_dir, "step25v2_pig_chronic_synovium_gene_level_length_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  main_counts_df,
  file.path(tables_dir, "step25v2_pig_chronic_main_comparison_gene_level_counts_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  main_abundance_df,
  file.path(tables_dir, "step25v2_pig_chronic_main_comparison_gene_level_abundance_tpm_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  main_length_df,
  file.path(tables_dir, "step25v2_pig_chronic_main_comparison_gene_level_length_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  summary_df,
  file.path(tables_dir, "step25v2_pig_chronic_manifest_rebuild_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

save(
  synovium_manifest_df, main_manifest_df,
  synovium_counts_df, synovium_abundance_df, synovium_length_df,
  main_counts_df, main_abundance_df, main_length_df,
  summary_df,
  file = file.path(objects_dir, "step25v2_pig_chronic_manifest_rebuild_workspace.RData")
)

writeLines(
  c(
    "Step 25 v2 completed successfully.",
    paste("Preferred manifest file:", preferred_manifest_file),
    paste("Synovium samples:", nrow(synovium_manifest_df)),
    paste("Main comparison samples:", nrow(main_manifest_df)),
    paste("Main Control_52W:", sum(main_manifest_df$core_group == "Control_52W", na.rm = TRUE)),
    paste("Main ACLT_alone_52W:", sum(main_manifest_df$core_group == "ACLT_alone_52W", na.rm = TRUE)),
    paste("Synovium counts dim:", paste(nrow(synovium_counts_df), "x", ncol(synovium_counts_df) - 1)),
    paste("Main counts dim:", paste(nrow(main_counts_df), "x", ncol(main_counts_df) - 1))
  ),
  file.path(logs_dir, "step25v2_pig_chronic_manifest_rebuild_log.txt"),
  useBytes = TRUE
)

writeLines(
  c(
    "# Step 25 v2 run script archived automatically.",
    "setwd(\"E:/R/ACLsenescence2\")"
  ),
  con = file.path(scripts_dir, "step25v2_run_rebuild_pig_chronic_manifests_and_subsets.R"),
  useBytes = TRUE
)

message("Step 25 v2 finished successfully.")
