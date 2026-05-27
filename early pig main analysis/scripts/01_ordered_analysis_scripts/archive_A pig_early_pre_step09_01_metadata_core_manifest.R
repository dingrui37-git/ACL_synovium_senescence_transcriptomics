# ============================================================
# Pig early E-MTAB-6664 upstream Step 01
# Purpose:
#   Download/locate E-MTAB-6664 IDF/SDRF metadata, read IDF/SDRF,
#   build analysis-ready sample metadata, define the 18 core synovium
#   validation samples, and generate the core FASTQ manifest used by
#   the later rebuild Step09.
#
# Manuscript-methods mapping:
#   Metadata processing and sample grouping for E-MTAB-6664.
#   Core groups retained:
#     - Control: CON_t0
#     - ACL transection untreated 1 week: ACLT_untreated_t7
#     - ACL transection untreated 4 weeks: ACLT_untreated_t28
#   Reconstruction and BEAR/repair samples are not retained in the
#   primary early validation manifest.
# ============================================================

setwd("D:/R/ACL ∩ senescence2")

project_root <- getwd()
data_dir <- file.path(project_root, "data_raw", "E-MTAB-6664")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 0) Helpers
# -----------------------------
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

safe_col <- function(df, col, default = NA_character_) {
  if (col %in% colnames(df)) {
    as.character(df[[col]])
  } else {
    rep(default, nrow(df))
  }
}

# -----------------------------
# 1) Obtain E-MTAB-6664 metadata files
# -----------------------------
# Original upstream record used ArrayExpress::getAE() to retrieve the
# E-MTAB-6664 metadata into data_raw/E-MTAB-6664.
idf_file <- file.path(data_dir, "E-MTAB-6664.idf.txt")
sdrf_file <- file.path(data_dir, "E-MTAB-6664.sdrf.txt")

if (!file.exists(idf_file) || !file.exists(sdrf_file)) {
  install_if_missing("ArrayExpress", bioc = TRUE)
  library(ArrayExpress)

  message("IDF/SDRF not found locally. Attempting to download E-MTAB-6664 metadata with ArrayExpress::getAE().")
  ae_6664 <- ArrayExpress::getAE(
    accession = "E-MTAB-6664",
    type = "processed",
    path = data_dir
  )
} else {
  message("IDF/SDRF already found locally. Skipping metadata download.")
}

if (!file.exists(idf_file) || !file.exists(sdrf_file)) {
  stop(
    "IDF/SDRF files are still missing after download attempt:\n",
    idf_file, "\n", sdrf_file,
    "\nPlease check ArrayExpress/EBI download status or place the files manually."
  )
}

# Save a file inventory for provenance
files_now <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
write.csv(
  data.frame(file = files_now, basename = basename(files_now), stringsAsFactors = FALSE),
  file = file.path(data_dir, "E-MTAB-6664_local_file_inventory_after_metadata_download.csv"),
  row.names = FALSE
)

# -----------------------------
# 2) Read IDF/SDRF and inspect key columns
# -----------------------------
idf_6664 <- read.delim(
  idf_file,
  header = FALSE,
  sep = "\t",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE
)

sdrf_6664 <- read.delim(
  sdrf_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

key_col_idx <- grepl(
  "sample|source|characteristics|factor|organism|tissue|disease|time|treatment|assay|comment|file|run|library|ENA|FASTQ",
  colnames(sdrf_6664),
  ignore.case = TRUE
)
key_cols <- colnames(sdrf_6664)[key_col_idx]

if (length(key_cols) > 0) {
  sdrf_key_6664 <- sdrf_6664[, key_cols, drop = FALSE]
  write.csv(
    sdrf_key_6664,
    file = file.path(data_dir, "E-MTAB-6664_SDRF_key_columns_preview.csv"),
    row.names = FALSE
  )
}

# -----------------------------
# 3) Build one-row-per-sample metadata
# -----------------------------
if (!"Assay Name" %in% colnames(sdrf_6664)) {
  stop("Required SDRF column missing: Assay Name")
}

sample_meta_6664 <- sdrf_6664[!duplicated(sdrf_6664$`Assay Name`), , drop = FALSE]

key_cols_order <- c(
  "Source Name",
  "Assay Name",
  "Characteristics[organism part]",
  "Characteristics[sampling site]",
  "Characteristics[injury]",
  "Characteristics[treatment]",
  "Factor Value[injury]",
  "Factor Value[treatment]",
  "Factor Value[time]",
  "Unit[time unit]",
  "Characteristics[age]",
  "Characteristics[sex]",
  "Comment[ENA_SAMPLE]",
  "Comment[ENA_EXPERIMENT]",
  "Comment[ENA_RUN]"
)
key_cols_order <- key_cols_order[key_cols_order %in% colnames(sample_meta_6664)]
sample_meta_key_6664 <- sample_meta_6664[, key_cols_order, drop = FALSE]

write.csv(
  sample_meta_key_6664,
  file = file.path(data_dir, "E-MTAB-6664_one_sample_per_row_key_metadata.csv"),
  row.names = FALSE
)

# -----------------------------
# 4) Build analysis-ready metadata and define core validation groups
# -----------------------------
meta_6664_clean <- sample_meta_6664

meta_6664_clean$sample_id      <- safe_col(meta_6664_clean, "Assay Name")
meta_6664_clean$source_name    <- safe_col(meta_6664_clean, "Source Name")
meta_6664_clean$injury_raw     <- safe_col(meta_6664_clean, "Factor Value[injury]")
meta_6664_clean$treatment_raw  <- safe_col(meta_6664_clean, "Factor Value[treatment]")
meta_6664_clean$time_raw       <- safe_col(meta_6664_clean, "Factor Value[time]")
meta_6664_clean$time_unit_raw  <- safe_col(meta_6664_clean, "Unit[time unit]")
meta_6664_clean$organism_part  <- safe_col(meta_6664_clean, "Characteristics[organism part]")
meta_6664_clean$sampling_site  <- safe_col(meta_6664_clean, "Characteristics[sampling site]")
meta_6664_clean$sex            <- safe_col(meta_6664_clean, "Characteristics[sex]")
meta_6664_clean$age            <- safe_col(meta_6664_clean, "Characteristics[age]")

meta_6664_clean$treatment_label <- ifelse(
  meta_6664_clean$treatment_raw == "none", "untreated",
  ifelse(
    meta_6664_clean$treatment_raw == "ligament reconstruction surgery", "reconstruction",
    ifelse(
      meta_6664_clean$treatment_raw == "ligament repair surgery (BEAR)", "repair",
      ifelse(meta_6664_clean$treatment_raw == "not applicable", "control", "other")
    )
  )
)

meta_6664_clean$injury_label <- ifelse(
  meta_6664_clean$injury_raw == "none", "control",
  ifelse(
    grepl("anterior cruciate ligament transection", meta_6664_clean$injury_raw, ignore.case = TRUE),
    "ACLT",
    "other"
  )
)

meta_6664_clean$study_group_raw <- paste(
  meta_6664_clean$injury_label,
  meta_6664_clean$treatment_label,
  paste0("t", meta_6664_clean$time_raw),
  sep = "__"
)

meta_6664_clean$core_group <- NA_character_
meta_6664_clean$core_group[meta_6664_clean$injury_label == "control"] <- "CON_t0"
meta_6664_clean$core_group[
  meta_6664_clean$injury_label == "ACLT" &
    meta_6664_clean$treatment_label == "untreated" &
    as.character(meta_6664_clean$time_raw) == "7"
] <- "ACLT_untreated_t7"
meta_6664_clean$core_group[
  meta_6664_clean$injury_label == "ACLT" &
    meta_6664_clean$treatment_label == "untreated" &
    as.character(meta_6664_clean$time_raw) == "28"
] <- "ACLT_untreated_t28"

# This flag is deliberately restricted to the three manuscript core groups.
meta_6664_clean$core_validation_use <- !is.na(meta_6664_clean$core_group)

meta_6664_analysis <- meta_6664_clean[, c(
  "sample_id",
  "source_name",
  "injury_raw",
  "treatment_raw",
  "time_raw",
  "time_unit_raw",
  "injury_label",
  "treatment_label",
  "study_group_raw",
  "core_validation_use",
  "core_group",
  "organism_part",
  "sampling_site",
  "sex",
  "age"
)]

write.csv(
  meta_6664_analysis,
  file = file.path(data_dir, "E-MTAB-6664_analysis_ready_metadata.csv"),
  row.names = FALSE
)

# -----------------------------
# 5) Freeze the exact 18 core samples and the two early comparisons
# -----------------------------
core_meta_6664 <- subset(meta_6664_analysis, core_validation_use)

core_meta_6664 <- core_meta_6664[, c(
  "sample_id",
  "source_name",
  "core_group",
  "injury_label",
  "treatment_label",
  "time_raw",
  "time_unit_raw",
  "sex",
  "age"
)]

ena_cols <- c("Assay Name", "Comment[ENA_SAMPLE]", "Comment[ENA_EXPERIMENT]", "Comment[ENA_RUN]")
if (!all(ena_cols %in% colnames(sample_meta_6664))) {
  stop("One or more ENA columns are missing from SDRF metadata.")
}

ena_map_6664 <- sample_meta_6664[!duplicated(sample_meta_6664$`Assay Name`), ena_cols]
colnames(ena_map_6664) <- c("sample_id", "ENA_SAMPLE", "ENA_EXPERIMENT", "ENA_RUN")

core_meta_6664 <- merge(core_meta_6664, ena_map_6664, by = "sample_id", all.x = TRUE, sort = FALSE)

core_meta_6664$core_group <- factor(
  core_meta_6664$core_group,
  levels = c("CON_t0", "ACLT_untreated_t7", "ACLT_untreated_t28")
)
core_meta_6664 <- core_meta_6664[order(core_meta_6664$core_group, core_meta_6664$sample_id), ]
core_meta_6664$core_group <- as.character(core_meta_6664$core_group)

core_group_counts <- table(core_meta_6664$core_group)
expected_counts <- c(CON_t0 = 6, ACLT_untreated_t7 = 6, ACLT_untreated_t28 = 6)
if (!all(names(expected_counts) %in% names(core_group_counts)) ||
    !all(as.integer(core_group_counts[names(expected_counts)]) == as.integer(expected_counts))) {
  print(core_group_counts)
  stop("Core group counts are not 6/6/6. Please inspect E-MTAB-6664_analysis_ready_metadata.csv.")
}

meta_6664_t7 <- subset(core_meta_6664, core_group %in% c("CON_t0", "ACLT_untreated_t7"))
meta_6664_t28 <- subset(core_meta_6664, core_group %in% c("CON_t0", "ACLT_untreated_t28"))

write.csv(core_meta_6664, file = file.path(data_dir, "E-MTAB-6664_core_validation_metadata_18samples.csv"), row.names = FALSE)
write.csv(meta_6664_t7, file = file.path(data_dir, "E-MTAB-6664_core_validation_metadata_t7_vs_control.csv"), row.names = FALSE)
write.csv(meta_6664_t28, file = file.path(data_dir, "E-MTAB-6664_core_validation_metadata_t28_vs_control.csv"), row.names = FALSE)

# -----------------------------
# 6) Build long and wide FASTQ manifest for the 18 core samples
# -----------------------------
required_manifest_cols <- c(
  "Assay Name",
  "Source Name",
  "Comment[ENA_SAMPLE]",
  "Comment[ENA_EXPERIMENT]",
  "Comment[ENA_RUN]",
  "Comment[SUBMITTED_FILE_NAME]",
  "Comment[FASTQ_URI]"
)
if (!all(required_manifest_cols %in% colnames(sdrf_6664))) {
  missing_cols <- setdiff(required_manifest_cols, colnames(sdrf_6664))
  stop("Required SDRF manifest columns missing: ", paste(missing_cols, collapse = ", "))
}

sdrf_core_6664 <- subset(sdrf_6664, `Assay Name` %in% core_meta_6664$sample_id)

manifest_long_6664 <- sdrf_core_6664[, required_manifest_cols]
colnames(manifest_long_6664) <- c(
  "sample_id",
  "source_name",
  "ENA_SAMPLE",
  "ENA_EXPERIMENT",
  "ENA_RUN",
  "submitted_file_name",
  "FASTQ_URI"
)

manifest_long_6664$read_label <- ifelse(
  grepl("_R1\\.fastq\\.gz$|_1\\.fastq\\.gz$", manifest_long_6664$submitted_file_name, ignore.case = TRUE),
  "R1",
  ifelse(
    grepl("_R2\\.fastq\\.gz$|_2\\.fastq\\.gz$", manifest_long_6664$submitted_file_name, ignore.case = TRUE),
    "R2",
    NA_character_
  )
)

manifest_long_6664 <- manifest_long_6664[order(manifest_long_6664$sample_id, manifest_long_6664$read_label), ]

r1_df_6664 <- subset(
  manifest_long_6664,
  read_label == "R1",
  select = c("sample_id", "ENA_RUN", "submitted_file_name", "FASTQ_URI")
)
colnames(r1_df_6664) <- c("sample_id", "ENA_RUN", "R1_file", "R1_URI")

r2_df_6664 <- subset(
  manifest_long_6664,
  read_label == "R2",
  select = c("sample_id", "submitted_file_name", "FASTQ_URI")
)
colnames(r2_df_6664) <- c("sample_id", "R2_file", "R2_URI")

manifest_wide_6664 <- merge(r1_df_6664, r2_df_6664, by = "sample_id", all = TRUE, sort = FALSE)

manifest_wide_6664 <- merge(
  core_meta_6664[, c("sample_id", "core_group", "sex", "age")],
  manifest_wide_6664,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

manifest_wide_6664$core_group <- factor(
  manifest_wide_6664$core_group,
  levels = c("CON_t0", "ACLT_untreated_t7", "ACLT_untreated_t28")
)
manifest_wide_6664 <- manifest_wide_6664[order(manifest_wide_6664$core_group, manifest_wide_6664$sample_id), ]
manifest_wide_6664$core_group <- as.character(manifest_wide_6664$core_group)

if (nrow(manifest_wide_6664) != 18) {
  stop("Expected 18 rows in the wide core FASTQ manifest, but found: ", nrow(manifest_wide_6664))
}
if (any(is.na(manifest_wide_6664[, c("R1_file", "R1_URI", "R2_file", "R2_URI")]))) {
  print(colSums(is.na(manifest_wide_6664[, c("R1_file", "R1_URI", "R2_file", "R2_URI")])))
  stop("Missing R1/R2 file names or URIs in the core FASTQ manifest.")
}

write.csv(
  manifest_long_6664,
  file = file.path(data_dir, "E-MTAB-6664_core_validation_manifest_long.csv"),
  row.names = FALSE
)

write.csv(
  manifest_wide_6664,
  file = file.path(data_dir, "E-MTAB-6664_core_validation_manifest_wide.csv"),
  row.names = FALSE
)

summary_6664 <- data.frame(
  metric = c(
    "n_sdrf_rows",
    "n_unique_samples",
    "n_core_samples",
    "n_manifest_long_rows",
    "n_manifest_wide_rows"
  ),
  value = c(
    nrow(sdrf_6664),
    nrow(sample_meta_6664),
    nrow(core_meta_6664),
    nrow(manifest_long_6664),
    nrow(manifest_wide_6664)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_6664,
  file = file.path(data_dir, "E-MTAB-6664_metadata_core_manifest_summary.csv"),
  row.names = FALSE
)

cat("===== E-MTAB-6664 metadata and core manifest summary =====\n")
print(summary_6664)
cat("\n===== Core group counts =====\n")
print(table(core_meta_6664$core_group, useNA = "ifany"))
cat("\nSaved wide manifest:\n")
print(file.path(data_dir, "E-MTAB-6664_core_validation_manifest_wide.csv"))
