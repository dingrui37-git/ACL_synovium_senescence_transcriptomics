# StepH2w13_current78: fastp-trimmed sensitivity heatmap for the primary main-analysis 24 core genes
# Purpose:
# Generate a Supplementary Figure-style heatmap showing the expression of the primary
# untrimmed-analysis 24 core ortholog genes in the fastp-trimmed sensitivity branch.
#
# IMPORTANT INTERPRETATION:
#   - This heatmap is a sensitivity-audit figure.
#   - It does NOT redefine or replace the official main-analysis core ortholog gene set.
#   - The official 24 core genes remain defined by the primary untrimmed current78 analysis.
#   - This figure only visualizes whether those same primary 24 genes show a comparable
#     expression pattern in the trimmed count matrix.
#
# Style:
#   - Same style as the main Figure4D heatmap:
#       rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red
#       z-score clipped to [-2, 2]
#       pheatmap with Group, CellAgeEffect and MouseDirection annotations
#   - Expanded canvas to avoid clipping
#   - Generates exactly ONE PNG and ONE PDF
#
# Inputs:
#   StepH2w10 main-core preservation audit table:
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w10_current78_core_gene_sensitivity_audit/stepH2w10_current78_main24_core_genes_trimmed_preservation_table.csv
#   Trimmed count matrix:
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_gene_level_counts_matrix.csv
#   Trimmed sample info:
#     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_sample_info_for_DE.csv
#   Step11C current mouse persistent ∩ CellAge overlap table for CellAgeEffect annotation

options(stringsAsFactors = FALSE)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
}
install_if_missing(c("readr", "dplyr", "edgeR", "pheatmap", "RColorBrewer"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(edgeR)
  library(pheatmap)
  library(RColorBrewer)
})

## =========================
## 1. Helper functions
## =========================

find_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

stop_if_missing <- function(file, label) {
  if (!file.exists(file)) stop("Missing ", label, ": ", file, call. = FALSE)
}

find_col <- function(df, candidates, label) {
  idx <- match(candidates, colnames(df))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    stop(sprintf("Could not find column for %s. Checked: %s\nAvailable columns: %s",
                 label, paste(candidates, collapse = ", "), paste(colnames(df), collapse = ", ")),
         call. = FALSE)
  }
  colnames(df)[idx[1]]
}

find_col_optional <- function(df, candidates) {
  idx <- match(candidates, colnames(df))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NA_character_)
  colnames(df)[idx[1]]
}

clean_char <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}

as_logical_robust <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  y <- tolower(clean_char(x))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "yes", "y", "1")] <- TRUE
  out[y %in% c("false", "f", "no", "n", "0")] <- FALSE
  out
}

row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}

normalize_group <- function(x) {
  lx <- tolower(as.character(x))
  out <- rep(NA_character_, length(lx))
  out[grepl("control|ctrl|con|t0", lx)] <- "Control"
  out[grepl("t7|7d|day7|1w|week1", lx) & grepl("acl|aclt|untreated|inj|transection", lx)] <- "ACLT-untreated-1W"
  out[grepl("t28|28d|day28|4w|week4", lx) & grepl("acl|aclt|untreated|inj|transection", lx)] <- "ACLT-untreated-4W"
  out[is.na(out) & grepl("t7|7d|day7|1w|week1", lx)] <- "ACLT-untreated-1W"
  out[is.na(out) & grepl("t28|28d|day28|4w|week4", lx)] <- "ACLT-untreated-4W"
  out
}

## =========================
## 2. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

fastp_dir <- file.path(rebuild_root, "02_pig_early_fastp_sensitivity")
main_mouse_dir <- file.path(rebuild_root, "02_mouse_discovery")

stepH2w10_dir <- file.path(fastp_dir, "tables", "stepH2w10_current78_core_gene_sensitivity_audit")
out_dir <- file.path(fastp_dir, "tables", "stepH2w13_current78_trimmed_main24_core_heatmap_sensitivity")
fig_dir <- file.path(fastp_dir, "figures", "stepH2w13_current78_trimmed_main24_core_heatmap_sensitivity")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "stepH2w13_current78_trimmed_main24_core_heatmap_sensitivity_full_log.txt")
summary_log_file <- file.path(log_dir, "stepH2w13_current78_trimmed_main24_core_heatmap_sensitivity_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== StepH2w13_current78 trimmed main-24 core heatmap sensitivity figure =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Interpretation: visualizes primary untrimmed-analysis 24 core genes in the trimmed sensitivity branch; does not redefine core genes.\n\n")

main24_preservation_file <- file.path(stepH2w10_dir, "stepH2w10_current78_main24_core_genes_trimmed_preservation_table.csv")
trimmed_count_file <- file.path(fastp_dir, "tables", "stepH2w5_v2_trimmed_gene_level_counts_matrix.csv")
trimmed_sample_info_file <- file.path(fastp_dir, "tables", "stepH2w5_v2_trimmed_sample_info_for_DE.csv")

step11C_overlap_file <- find_existing_file(c(
  file.path(main_mouse_dir, "07_tables", "step11C_intersect_current_gorth_mapping_with_CellAge_clean", "step11C_persistent_CellAge_overlap_78_genes.csv"),
  file.path(main_mouse_dir, "07_tables", "step11C_intersect_current_gorth_mapping_with_CellAge_clean", "step11C_persistent_CellAge_overlap_78_genes_compact.csv")
))

stop_if_missing(main24_preservation_file, "StepH2w10 main 24 core preservation table")
stop_if_missing(trimmed_count_file, "trimmed count matrix")
stop_if_missing(trimmed_sample_info_file, "trimmed sample info")
if (is.na(step11C_overlap_file)) stop("Missing Step11C overlap file for CellAge effect annotation.", call. = FALSE)

cat("Main 24 preservation file: ", main24_preservation_file, "\n", sep = "")
cat("Trimmed count matrix: ", trimmed_count_file, "\n", sep = "")
cat("Trimmed sample info: ", trimmed_sample_info_file, "\n", sep = "")
cat("Step11C overlap file: ", step11C_overlap_file, "\n\n", sep = "")

## =========================
## 3. Read inputs
## =========================

core_df_raw <- read_csv(main24_preservation_file, show_col_types = FALSE)
count_df <- read_csv(trimmed_count_file, show_col_types = FALSE)
sample_info_raw <- read_csv(trimmed_sample_info_file, show_col_types = FALSE)
overlap_df <- read_csv(step11C_overlap_file, show_col_types = FALSE)

pig_ensg_col <- find_col(core_df_raw, c("pig_ensg", "pig_gene_id", "pig_gene", "gene_id", "ensembl_gene_id"), "pig Ensembl ID")
pig_symbol_col <- find_col_optional(core_df_raw, c("pig_symbol", "pig_gene_symbol", "gene_symbol", "symbol", "gene_name"))
human_gene_col <- find_col_optional(core_df_raw, c("human_gene", "ortholog_name", "Gene symbol", "gene_symbol_human"))
mouse_symbol_col <- find_col_optional(core_df_raw, c("mouse_symbol", "input", "mouse_gene", "mouse_gene_symbol"))
sig_dir_col <- find_col(core_df_raw, c("signature_direction", "mouse_direction", "direction"), "signature direction")

preserved_col <- find_col_optional(core_df_raw, c("core_like_trimmed", "main_core_genes_preserved_as_core_like_in_trimmed"))
strict_t7_col <- find_col_optional(core_df_raw, c("strict_t7_trimmed", "main_core_genes_strict_t7_trimmed"))
strict_t28_col <- find_col_optional(core_df_raw, c("strict_t28_trimmed", "main_core_genes_strict_t28_trimmed"))
dc_both_col <- find_col_optional(core_df_raw, c("direction_consistent_both_trimmed", "main_core_genes_direction_consistent_both_trimmed"))

core_df <- core_df_raw %>%
  mutate(
    pig_ensg = clean_char(.data[[pig_ensg_col]]),
    pig_symbol = if (!is.na(pig_symbol_col)) clean_char(.data[[pig_symbol_col]]) else NA_character_,
    human_gene = if (!is.na(human_gene_col)) clean_char(.data[[human_gene_col]]) else NA_character_,
    mouse_symbol = if (!is.na(mouse_symbol_col)) clean_char(.data[[mouse_symbol_col]]) else NA_character_,
    signature_direction = clean_char(.data[[sig_dir_col]])
  ) %>%
  filter(!is.na(pig_ensg), pig_ensg != "") %>%
  distinct(pig_ensg, .keep_all = TRUE)

if (nrow(core_df) != 24) {
  stop("Expected 24 primary main-analysis core genes from StepH2w10 preservation table, but found: ", nrow(core_df), call. = FALSE)
}

if (!is.na(preserved_col)) {
  preserved_flag <- as_logical_robust(core_df_raw[[preserved_col]])
  if (sum(preserved_flag, na.rm = TRUE) != 24) {
    warning("Not all 24 genes are preserved as core-like in trimmed branch according to StepH2w10. Please inspect the preservation table.")
  }
}

## =========================
## 4. Add CellAge effect annotation
## =========================

ov_human_col <- find_col(overlap_df, c("ortholog_name", "Gene symbol", "human_gene", "gene_symbol"), "Step11C human gene")
possible_effect_cols <- intersect(c("Senescence Effect", "Senescence_Effect", "senescence_effect", "CellAgeEffect"), colnames(overlap_df))
if (length(possible_effect_cols) == 0) stop("Could not find CellAge effect column in Step11C overlap file.", call. = FALSE)
effect_col <- possible_effect_cols[1]

effect_df <- overlap_df %>%
  transmute(human_gene = clean_char(.data[[ov_human_col]]), CellAgeEffect = clean_char(.data[[effect_col]])) %>%
  distinct(human_gene, .keep_all = TRUE)

core_df <- core_df %>% left_join(effect_df, by = "human_gene")
core_df$CellAgeEffect[is.na(core_df$CellAgeEffect)] <- "Unclear"
core_df$CellAgeEffect <- ifelse(
  grepl("Induce", core_df$CellAgeEffect, ignore.case = TRUE), "Induces",
  ifelse(grepl("Inhibit", core_df$CellAgeEffect, ignore.case = TRUE), "Inhibits", "Unclear")
)

core_df$MouseDirection <- ifelse(core_df$signature_direction == "Down_in_ACLR", "mouse_persistent_down", "mouse_persistent_up")

## =========================
## 5. Trimmed expression matrix and sample annotation
## =========================

count_gene_col <- find_col(count_df, c("gene_id", "Geneid", "ensembl_gene_id", "pig_ensg"), "trimmed count matrix gene ID")
non_count_cols <- c(count_gene_col, intersect(c("gene_symbol", "symbol", "gene_name", "gene", "Length", "length", "Chr", "Start", "End", "Strand"), colnames(count_df)))
sample_cols <- setdiff(colnames(count_df), non_count_cols)

count_mat <- as.matrix(count_df[, sample_cols, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- clean_char(count_df[[count_gene_col]])

valid_rows <- !is.na(rownames(count_mat)) & rownames(count_mat) != ""
count_mat <- count_mat[valid_rows, , drop = FALSE]
if (any(duplicated(rownames(count_mat)))) {
  count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
}

missing_core <- setdiff(core_df$pig_ensg, rownames(count_mat))
if (length(missing_core) > 0) {
  write_csv(data.frame(missing_pig_ensg = missing_core), file.path(out_dir, "stepH2w13_missing_main24_core_genes_in_trimmed_count_matrix.csv"))
  stop("Some primary 24 core genes are missing in the trimmed count matrix.", call. = FALSE)
}

sample_id_col <- find_col(sample_info_raw, c("sample_id", "sample", "SampleID", "sample_name"), "sample info sample ID")
group_col <- find_col(sample_info_raw, c("group", "condition", "time_group", "treatment_group", "group_label"), "sample info group")

annotation_df <- sample_info_raw %>%
  transmute(
    sample_id = clean_char(.data[[sample_id_col]]),
    group_raw = clean_char(.data[[group_col]])
  ) %>%
  filter(!is.na(sample_id), sample_id %in% colnames(count_mat)) %>%
  distinct(sample_id, .keep_all = TRUE)

annotation_df$group <- normalize_group(annotation_df$group_raw)
if (any(is.na(annotation_df$group))) {
  write_csv(annotation_df, file.path(out_dir, "stepH2w13_sample_group_normalization_failed.csv"))
  stop("Some sample groups could not be normalized to Control/ACLT-untreated-1W/ACLT-untreated-4W.", call. = FALSE)
}

annotation_df$group <- factor(annotation_df$group, levels = c("Control", "ACLT-untreated-1W", "ACLT-untreated-4W"))
annotation_df <- annotation_df %>% arrange(group, sample_id)

if (nrow(annotation_df) != 18) {
  warning("Expected 18 pig early core samples, but found: ", nrow(annotation_df))
}
if (!all(annotation_df$sample_id %in% colnames(count_mat))) {
  missing_samples <- setdiff(annotation_df$sample_id, colnames(count_mat))
  write_csv(data.frame(missing_sample_id = missing_samples), file.path(out_dir, "stepH2w13_missing_samples_in_trimmed_count_matrix.csv"))
  stop("Some samples are missing in the trimmed count matrix.", call. = FALSE)
}

core_counts <- count_mat[core_df$pig_ensg, annotation_df$sample_id, drop = FALSE]

# Match the main heatmap style: TMM-normalized logCPM, then row-wise z-score, clipped to [-2, 2].
dge <- DGEList(counts = core_counts)
dge <- calcNormFactors(dge, method = "TMM")
logcpm_mat <- cpm(dge, log = TRUE, prior.count = 1)
z_mat <- row_zscore(logcpm_mat)

core_df$display_label <- ifelse(!is.na(core_df$pig_symbol) & core_df$pig_symbol != "", core_df$pig_symbol, core_df$human_gene)

# Preserve the same biological ordering approach as main Figure4D before row clustering.
core_df$dir_order <- ifelse(core_df$MouseDirection == "mouse_persistent_down", 1, 2)
core_df$effect_order <- ifelse(core_df$CellAgeEffect == "Induces", 1, ifelse(core_df$CellAgeEffect == "Inhibits", 2, 3))
core_df <- core_df %>% arrange(dir_order, effect_order, display_label)

rownames(z_mat) <- ifelse(!is.na(core_df$display_label[match(rownames(z_mat), core_df$pig_ensg)]),
                          core_df$display_label[match(rownames(z_mat), core_df$pig_ensg)],
                          rownames(z_mat))
z_mat <- z_mat[core_df$display_label, annotation_df$sample_id, drop = FALSE]

annotation_col <- data.frame(Group = annotation_df$group)
rownames(annotation_col) <- annotation_df$sample_id

annotation_row <- data.frame(
  CellAgeEffect = factor(core_df$CellAgeEffect, levels = c("Induces", "Inhibits", "Unclear")),
  MouseDirection = factor(core_df$MouseDirection, levels = c("mouse_persistent_down", "mouse_persistent_up"))
)
rownames(annotation_row) <- core_df$display_label

## =========================
## 6. Heatmap style from main Figure4D
## =========================

heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(101)
heat_breaks <- seq(-2, 2, length.out = length(heat_colors) + 1)

ann_colors <- list(
  Group = c(
    Control = "#66CC00",
    "ACLT-untreated-1W" = "#F564E3",
    "ACLT-untreated-4W" = "#00BFC4"
  ),
  CellAgeEffect = c(
    Induces = "#D9B300",
    Inhibits = "#00C1A2",
    Unclear = "#BDBDBD"
  ),
  MouseDirection = c(
    mouse_persistent_down = "#C77CFF",
    mouse_persistent_up = "#F8766D"
  )
)

main_title <- sprintf(
  "Sensitivity audit: fastp-trimmed expression of primary %d core genes",
  nrow(core_df)
)

pdf_file <- file.path(fig_dir, "Supplementary_fastp_trimmed_primary24_core_gene_heatmap_step09_style.pdf")
png_file <- file.path(fig_dir, "Supplementary_fastp_trimmed_primary24_core_gene_heatmap_step09_style.png")

make_heatmap <- function(filename, width, height) {
  pheatmap::pheatmap(
    mat = z_mat,
    color = heat_colors,
    breaks = heat_breaks,
    scale = "none",

    cluster_rows = TRUE,
    cluster_cols = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",

    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = ann_colors,

    show_rownames = TRUE,
    show_colnames = TRUE,
    labels_row = rownames(z_mat),

    annotation_names_row = FALSE,
    annotation_names_col = TRUE,

    cellwidth = 22,
    cellheight = 18,
    fontsize = 11,
    fontsize_row = 9,
    fontsize_col = 9,
    angle_col = 270,

    border_color = "#A8A8A8",
    treeheight_row = 55,
    treeheight_col = 0,
    legend = TRUE,
    main = main_title,

    filename = filename,
    width = width,
    height = height
  )
}

# Generate exactly one PDF and one PNG. No duplicate copies.
make_heatmap(filename = pdf_file, width = 13.5, height = 8.2)
make_heatmap(filename = png_file, width = 13.5, height = 8.2)

## =========================
## 7. Save source data and summary
## =========================

core_table_file <- file.path(out_dir, "stepH2w13_primary24_core_genes_used_for_trimmed_heatmap.csv")
z_matrix_file <- file.path(out_dir, "stepH2w13_trimmed_primary24_core_gene_heatmap_zscore_matrix.csv")
logcpm_file <- file.path(out_dir, "stepH2w13_trimmed_primary24_core_gene_logCPM_matrix.csv")
sample_annotation_file <- file.path(out_dir, "stepH2w13_trimmed_core_sample_annotation.csv")

write_csv(
  core_df %>%
    select(human_gene, mouse_symbol, pig_ensg, pig_symbol, display_label, CellAgeEffect, MouseDirection, signature_direction),
  core_table_file
)
write_csv(as.data.frame(z_mat) %>% mutate(display_label = rownames(z_mat), .before = 1), z_matrix_file)
write_csv(as.data.frame(logcpm_mat) %>% mutate(pig_ensg = rownames(logcpm_mat), .before = 1), logcpm_file)
write_csv(annotation_df, sample_annotation_file)

preservation_summary <- data.frame(
  metric = c(
    "primary_main_analysis_core_genes_visualized",
    "trimmed_samples_visualized",
    "control_n",
    "aclt_1w_n",
    "aclt_4w_n",
    "heatmap_color_scale",
    "zscore_clip",
    "figure_role",
    "official_core_interpretation",
    "png_file",
    "pdf_file",
    "core_table_file",
    "z_matrix_file",
    "sample_annotation_file",
    "output_dir"
  ),
  value = c(
    nrow(core_df),
    nrow(annotation_df),
    sum(annotation_df$group == "Control"),
    sum(annotation_df$group == "ACLT-untreated-1W"),
    sum(annotation_df$group == "ACLT-untreated-4W"),
    "rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red",
    "[-2, 2]",
    "supplementary fastp-trimming sensitivity heatmap",
    "Official core genes remain defined by the primary untrimmed analysis; this figure visualizes those same 24 genes in the trimmed branch.",
    png_file,
    pdf_file,
    core_table_file,
    z_matrix_file,
    sample_annotation_file,
    out_dir
  ),
  stringsAsFactors = FALSE
)

summary_file <- file.path(out_dir, "stepH2w13_current78_trimmed_main24_core_heatmap_sensitivity_summary.csv")
write_csv(preservation_summary, summary_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "edgeR_version",
    "pheatmap_version",
    "RColorBrewer_version",
    "input_core_definition",
    "trimmed_expression_source",
    "normalization",
    "row_standardization",
    "interpretation"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("pheatmap")),
    as.character(utils::packageVersion("RColorBrewer")),
    main24_preservation_file,
    trimmed_count_file,
    "TMM-normalized logCPM",
    "row-wise z-score clipped to [-2, 2]",
    "sensitivity visualization only; does not redefine the official core gene set"
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "stepH2w13_versions_and_method_records.csv")
write_csv(version_df, version_file)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "stepH2w13_sessionInfo.txt"))

cat("\n===== StepH2w13_current78 trimmed main-24 core heatmap summary =====\n")
print(preservation_summary, row.names = FALSE)

cat("\nCore genes preview:\n")
print(head(core_df %>% select(display_label, human_gene, pig_ensg, CellAgeEffect, MouseDirection, signature_direction), 20), row.names = FALSE)

cat("\nSample annotation:\n")
print(annotation_df, row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

cat("\nStepH2w13_current78 trimmed main-24 core heatmap completed successfully.\n")

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== StepH2w13_current78 trimmed main-24 core heatmap SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(preservation_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Core genes preview:", summary_con)
writeLines(capture.output(print(head(core_df %>% select(display_label, human_gene, pig_ensg, CellAgeEffect, MouseDirection, signature_direction), 20), row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Sample annotation:", summary_con)
writeLines(capture.output(print(annotation_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This is a sensitivity-audit heatmap of the primary main-analysis 24 core genes in the fastp-trimmed count matrix.", summary_con)
writeLines("It does not redefine the official core gene set and should not replace the main Figure4D heatmap.", summary_con)
close(summary_con)

sink()

cat("\nStepH2w13_current78 trimmed main-24 core heatmap completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
