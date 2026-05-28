# StepCHRONMDS_01: Chronic pig main-comparison PCA and MDS
# Purpose: Perform sample-level PCA and MDS for the locked 24-sample chronic pig synovium main comparison
#          (Control_52W vs ACLT_alone_52W) using gene-level estimated counts, edgeR::filterByExpr,
#          TMM normalization, PCA by prcomp, and MDS by limma::plotMDS. This script does not round
#          Salmon/tximport estimated counts.

options(stringsAsFactors = FALSE)

status <- "FAILED"
method_version <- "2026-05-10_chronic_pig_main24_PCA_MDS_unified_with_mouse_and_early_pig"

project_root <- "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic"
input_dir <- file.path(project_root, "tables", "step24_chronic_synovium_only_matrix_lock", "tables")
counts_file <- file.path(input_dir, "StepCHRONLOCK_01_main_24_gene_level_estimated_counts_matrix.csv")
metadata_file <- file.path(input_dir, "StepCHRONLOCK_01_main_24_sample_metadata_locked.csv")
output_root <- file.path(project_root, "supplement", "MDS_PCA")

subdirs <- c("logs", "scripts", "tables", "source_data", "figures", "objects")
for (d in subdirs) dir.create(file.path(output_root, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(output_root, "logs", "StepCHRONMDS_01_chronic_pig_main24_PCA_MDS_summary_log.txt")
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")

close_log <- function() {
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}

fail <- function(msg) {
  cat("\nERROR:\n", msg, "\n", sep = "")
  close_log()
  stop(msg, call. = FALSE)
}

on.exit({
  cat("\n============================================================\n")
  cat("StepCHRONMDS_01_chronic_pig_main24_PCA_MDS finished with status:", status, "\n")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Summary log saved to:\n", log_file, "\n", sep = "")
  cat("============================================================\n")
  close_log()
}, add = TRUE)

cat("============================================================\n")
cat("StepCHRONMDS_01_chronic_pig_main24_PCA_MDS\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Method version:", method_version, "\n")
cat("============================================================\n\n")

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
  if (file.exists(script_path)) {
    archive_path <- file.path(output_root, "scripts", basename(script_path))
    ok <- file.copy(script_path, archive_path, overwrite = TRUE)
    cat("Script archived to:", archive_path, " copy_status=", ok, "\n", sep = "")
  }
} else {
  cat("Script archive note: --file argument was not detected. If running interactively, manually save this script in:\n")
  cat(file.path(output_root, "scripts"), "\n\n")
}

required_pkgs <- c("edgeR", "limma", "ggplot2", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) fail(paste0("Required package not installed: ", pkg))
}
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
})
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

cat("Required packages loaded.\n")
cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
cat("limma version:", as.character(packageVersion("limma")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")
cat("ggrepel available for label plots:", has_ggrepel, "\n\n")

cat("Input files:\n")
cat("Counts matrix:", counts_file, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Output root:", output_root, "\n\n")

if (!file.exists(counts_file)) fail(paste0("Missing counts matrix: ", counts_file))
if (!file.exists(metadata_file)) fail(paste0("Missing metadata file: ", metadata_file))

counts_df <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("Counts table dimension:", nrow(counts_df), "rows x", ncol(counts_df), "columns\n")
cat("First 10 count table columns:\n")
print(head(colnames(counts_df), 10))
cat("\nMetadata dimension:", nrow(meta), "rows x", ncol(meta), "columns\n")
cat("Metadata columns:\n")
print(colnames(meta))

# Identify gene id column.
gene_col_candidates <- c("gene_id", "GeneID", "gene", "Gene", "ensembl_gene_id", "pig_ensg")
gene_col <- gene_col_candidates[gene_col_candidates %in% colnames(counts_df)][1]
if (is.na(gene_col) || length(gene_col) == 0) gene_col <- colnames(counts_df)[1]
cat("Inferred gene ID column:", gene_col, "\n")

# Identify sample and group columns in metadata.
sample_col_candidates <- c("sample_id", "sample", "SampleID", "sample_name", "geo_accession", "GSM", "matrix_sample_id")
group_col_candidates <- c("group", "group_raw", "condition", "treatment", "Group")
sample_col <- sample_col_candidates[sample_col_candidates %in% colnames(meta)][1]
group_col <- group_col_candidates[group_col_candidates %in% colnames(meta)][1]
if (is.na(sample_col) || length(sample_col) == 0) fail("Could not infer sample column from metadata.")
if (is.na(group_col) || length(group_col) == 0) fail("Could not infer group column from metadata.")
cat("Inferred metadata sample column:", sample_col, "\n")
cat("Inferred metadata group column:", group_col, "\n\n")

meta[[sample_col]] <- as.character(meta[[sample_col]])
meta[[group_col]] <- as.character(meta[[group_col]])

# Restrict to main comparison groups if metadata contains extra groups by accident.
main_groups <- c("Control_52W", "ACLT_alone_52W")
if (all(main_groups %in% unique(meta[[group_col]]))) {
  meta <- meta[meta[[group_col]] %in% main_groups, , drop = FALSE]
}

# Standardize metadata column names for downstream work.
meta_std <- meta
meta_std$sample_id <- meta[[sample_col]]
meta_std$group <- meta[[group_col]]
meta_std$group <- factor(meta_std$group, levels = c("Control_52W", "ACLT_alone_52W"))
if (any(is.na(meta_std$group))) {
  fail("Metadata group contains values outside Control_52W and ACLT_alone_52W after standardization.")
}

# Ensure sample columns are present in counts matrix.
count_sample_cols <- setdiff(colnames(counts_df), gene_col)
missing_in_counts <- setdiff(meta_std$sample_id, count_sample_cols)
extra_in_counts <- setdiff(count_sample_cols, meta_std$sample_id)
cat("Metadata sample count:", nrow(meta_std), "\n")
cat("Count matrix sample columns:", length(count_sample_cols), "\n")
cat("Missing metadata samples in count matrix:", ifelse(length(missing_in_counts) == 0, "none", paste(missing_in_counts, collapse = "; ")), "\n")
cat("Extra count matrix samples not in metadata:", ifelse(length(extra_in_counts) == 0, "none", paste(extra_in_counts, collapse = "; ")), "\n\n")
if (length(missing_in_counts) > 0) fail("Some metadata samples are missing from the count matrix.")

# Reorder count matrix columns to metadata order.
counts_raw <- counts_df[, meta_std$sample_id, drop = FALSE]
rownames(counts_raw) <- counts_df[[gene_col]]

# Convert to numeric without rounding Salmon/tximport estimated counts.
counts_mat <- as.matrix(counts_raw)
suppressWarnings(storage.mode(counts_mat) <- "numeric")
if (anyNA(counts_mat)) fail("Counts matrix contains NA after numeric conversion.")
if (any(counts_mat < 0, na.rm = TRUE)) fail("Counts matrix contains negative values.")
if (anyDuplicated(rownames(counts_mat))) fail("Duplicate gene IDs detected in counts matrix.")

cat("Validated analysis count matrix dimension:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")
cat("Group counts:\n")
print(table(meta_std$group))
cat("\n")

# Save standardized metadata and raw analysis matrix for traceability.
write.csv(meta_std, file.path(output_root, "tables", "StepCHRONMDS_01_main24_metadata_used.csv"), row.names = FALSE)
raw_counts_out <- data.frame(gene_id = rownames(counts_mat), counts_mat, check.names = FALSE)
write.csv(raw_counts_out, file.path(output_root, "source_data", "StepCHRONMDS_01_main24_gene_level_estimated_counts_matrix_used.csv"), row.names = FALSE)

# edgeR filtering and TMM normalization.
dge <- DGEList(counts = counts_mat, group = meta_std$group)
keep <- filterByExpr(dge, group = meta_std$group)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")
logcpm <- cpm(dge_filtered, log = TRUE, prior.count = 1)

filter_norm_summary <- data.frame(
  metric = c(
    "method_version",
    "input_count_matrix",
    "input_metadata",
    "genes_before_filterByExpr",
    "genes_after_filterByExpr",
    "genes_removed_by_filterByExpr",
    "samples",
    "groups",
    "normalization_method",
    "logCPM_prior_count",
    "PCA_method",
    "MDS_method",
    "MDS_input",
    "MDS_top_genes",
    "MDS_gene_selection",
    "Salmon_tximport_counts_note"
  ),
  value = c(
    method_version,
    counts_file,
    metadata_file,
    nrow(counts_mat),
    nrow(logcpm),
    nrow(counts_mat) - nrow(logcpm),
    ncol(counts_mat),
    paste(names(table(meta_std$group)), as.integer(table(meta_std$group)), collapse = "; "),
    "edgeR_TMM",
    1,
    "prcomp(t(filtered_TMM_logCPM), center=TRUE, scale.=FALSE)",
    "limma::plotMDS(filtered_TMM_logCPM, top=500, gene.selection='common')",
    "filtered TMM-normalized logCPM matrix",
    500,
    "common",
    "Estimated counts from Salmon/tximport were used as numeric values and were not rounded."
  ),
  stringsAsFactors = FALSE
)
write.csv(filter_norm_summary, file.path(output_root, "tables", "StepCHRONMDS_01_method_parameters.csv"), row.names = FALSE)

cat("Filtering and normalization summary:\n")
print(filter_norm_summary)
cat("\n")

# Library summary.
library_summary <- data.frame(
  sample_id = colnames(counts_mat),
  group = as.character(meta_std$group),
  raw_estimated_count_sum = colSums(counts_mat),
  filtered_estimated_count_sum = colSums(dge_filtered$counts),
  norm_factor_TMM = dge_filtered$samples$norm.factors,
  effective_library_size = dge_filtered$samples$lib.size * dge_filtered$samples$norm.factors,
  detected_genes_count_gt0_raw = colSums(counts_mat > 0),
  detected_genes_count_gt0_filtered = colSums(dge_filtered$counts > 0),
  stringsAsFactors = FALSE
)
write.csv(library_summary, file.path(output_root, "tables", "StepCHRONMDS_01_main24_library_summary_after_TMM.csv"), row.names = FALSE)

cat("Library summary after TMM:\n")
print(library_summary)
cat("\n")

# Save filtered logCPM.
logcpm_out <- data.frame(gene_id = rownames(logcpm), logcpm, check.names = FALSE)
write.csv(logcpm_out, file.path(output_root, "source_data", "StepCHRONMDS_01_main24_filtered_TMM_logCPM_matrix.csv"), row.names = FALSE)
saveRDS(list(dge_filtered = dge_filtered, logcpm = logcpm, metadata = meta_std), file.path(output_root, "objects", "StepCHRONMDS_01_main24_edgeR_logCPM_objects.rds"))

# PCA.
pca_obj <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
pca_var <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
pca_var_df <- data.frame(
  PC = paste0("PC", seq_along(pca_var)),
  variance_fraction = pca_var,
  variance_percent = pca_var * 100,
  stringsAsFactors = FALSE
)
pca_coord <- data.frame(
  sample_id = rownames(pca_obj$x),
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  PC3 = if (ncol(pca_obj$x) >= 3) pca_obj$x[, 3] else NA_real_,
  group = as.character(meta_std$group),
  stringsAsFactors = FALSE
)
write.csv(pca_var_df, file.path(output_root, "tables", "StepCHRONMDS_01_main24_PCA_variance_explained.csv"), row.names = FALSE)
write.csv(pca_coord, file.path(output_root, "source_data", "StepCHRONMDS_01_main24_PCA_coordinates.csv"), row.names = FALSE)

cat("PCA variance explained by first five PCs:\n")
print(head(pca_var_df, 5))
cat("\n")

# MDS using limma::plotMDS on filtered TMM logCPM matrix.
mds_obj <- limma::plotMDS(logcpm, top = 500, gene.selection = "common", plot = FALSE)
mds_coord <- data.frame(
  sample_id = colnames(logcpm),
  MDS1 = as.numeric(mds_obj$x),
  MDS2 = as.numeric(mds_obj$y),
  group = as.character(meta_std$group),
  stringsAsFactors = FALSE
)
write.csv(mds_coord, file.path(output_root, "source_data", "StepCHRONMDS_01_main24_MDS_coordinates.csv"), row.names = FALSE)

cat("MDS object diagnostic:\n")
cat("MDS function used: limma::plotMDS\n")
cat("MDS input object: filtered TMM-normalized logCPM matrix\n")
cat("length(mds_obj$x):", length(mds_obj$x), "\n")
cat("length(mds_obj$y):", length(mds_obj$y), "\n")
cat("length(colnames(logcpm)):", length(colnames(logcpm)), "\n")
cat("MDS coordinate table preview:\n")
print(mds_coord)
cat("\n")

# Plot helpers.
plot_base_theme <- function() {
  theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}

shape_values <- c("Control_52W" = 16, "ACLT_alone_52W" = 17)

p_pca <- ggplot(pca_coord, aes(x = PC1, y = PC2, shape = group)) +
  geom_point(size = 3.4) +
  scale_shape_manual(values = shape_values, drop = FALSE) +
  labs(
    title = "Chronic pig synovium PCA",
    subtitle = "Main comparison: Control_52W vs ACLT_alone_52W; filtered TMM-normalized logCPM",
    x = paste0("PC1 (", sprintf("%.1f", pca_var_df$variance_percent[1]), "%)"),
    y = paste0("PC2 (", sprintf("%.1f", pca_var_df$variance_percent[2]), "%)"),
    shape = "Group"
  ) +
  plot_base_theme()

p_mds <- ggplot(mds_coord, aes(x = MDS1, y = MDS2, shape = group)) +
  geom_point(size = 3.4) +
  scale_shape_manual(values = shape_values, drop = FALSE) +
  labs(
    title = "Chronic pig synovium MDS",
    subtitle = "Main comparison; limma::plotMDS on filtered TMM-normalized logCPM; top 500 genes; common gene selection",
    x = "MDS dimension 1",
    y = "MDS dimension 2",
    shape = "Group"
  ) +
  plot_base_theme()

p_pca_label <- p_pca
p_mds_label <- p_mds
if (has_ggrepel) {
  p_pca_label <- p_pca + ggrepel::geom_text_repel(aes(label = sample_id), size = 3.0, max.overlaps = Inf, show.legend = FALSE)
  p_mds_label <- p_mds + ggrepel::geom_text_repel(aes(label = sample_id), size = 3.0, max.overlaps = Inf, show.legend = FALSE)
} else {
  p_pca_label <- p_pca + geom_text(aes(label = sample_id), size = 3.0, vjust = -0.8, show.legend = FALSE)
  p_mds_label <- p_mds + geom_text(aes(label = sample_id), size = 3.0, vjust = -0.8, show.legend = FALSE)
}

fig_files <- c(
  PCA_unlabeled_pdf = file.path(output_root, "figures", "StepCHRONMDS_01_main24_PCA_unlabeled.pdf"),
  PCA_labeled_pdf   = file.path(output_root, "figures", "StepCHRONMDS_01_main24_PCA_labeled.pdf"),
  MDS_unlabeled_pdf = file.path(output_root, "figures", "StepCHRONMDS_01_main24_MDS_unlabeled.pdf"),
  MDS_labeled_pdf   = file.path(output_root, "figures", "StepCHRONMDS_01_main24_MDS_labeled.pdf"),
  PCA_unlabeled_png = file.path(output_root, "figures", "StepCHRONMDS_01_main24_PCA_unlabeled.png"),
  MDS_unlabeled_png = file.path(output_root, "figures", "StepCHRONMDS_01_main24_MDS_unlabeled.png")
)

ggsave(fig_files["PCA_unlabeled_pdf"], p_pca, width = 6.6, height = 5.0, units = "in", device = cairo_pdf)
ggsave(fig_files["PCA_labeled_pdf"], p_pca_label, width = 7.2, height = 5.4, units = "in", device = cairo_pdf)
ggsave(fig_files["MDS_unlabeled_pdf"], p_mds, width = 6.6, height = 5.0, units = "in", device = cairo_pdf)
ggsave(fig_files["MDS_labeled_pdf"], p_mds_label, width = 7.2, height = 5.4, units = "in", device = cairo_pdf)
ggsave(fig_files["PCA_unlabeled_png"], p_pca, width = 6.6, height = 5.0, units = "in", dpi = 300)
ggsave(fig_files["MDS_unlabeled_png"], p_mds, width = 6.6, height = 5.0, units = "in", dpi = 300)

# Optional simple sample-distance matrix for audit.
sample_distance <- as.matrix(dist(t(logcpm), method = "euclidean"))
sample_distance_out <- data.frame(sample_id = rownames(sample_distance), sample_distance, check.names = FALSE)
write.csv(sample_distance_out, file.path(output_root, "source_data", "StepCHRONMDS_01_main24_sample_distance_matrix_from_filtered_logCPM.csv"), row.names = FALSE)

final_summary <- data.frame(
  metric = c(
    "run_status",
    "method_version",
    "input_counts_file",
    "input_metadata_file",
    "output_root",
    "genes_before_filterByExpr",
    "genes_after_filterByExpr",
    "samples",
    "groups",
    "normalization_method",
    "PCA_PC1_percent",
    "PCA_PC2_percent",
    "MDS_function",
    "MDS_input",
    "MDS1_min",
    "MDS1_max",
    "MDS2_min",
    "MDS2_max",
    "estimated_counts_rounded",
    "group_Control_52W_n",
    "group_ACLT_alone_52W_n"
  ),
  value = c(
    "SUCCESS",
    method_version,
    counts_file,
    metadata_file,
    output_root,
    nrow(counts_mat),
    nrow(logcpm),
    ncol(counts_mat),
    paste(names(table(meta_std$group)), as.integer(table(meta_std$group)), collapse = "; "),
    "edgeR_TMM",
    sprintf("%.6f", pca_var_df$variance_percent[1]),
    sprintf("%.6f", pca_var_df$variance_percent[2]),
    "limma::plotMDS",
    "filtered TMM-normalized logCPM matrix",
    sprintf("%.6f", min(mds_coord$MDS1)),
    sprintf("%.6f", max(mds_coord$MDS1)),
    sprintf("%.6f", min(mds_coord$MDS2)),
    sprintf("%.6f", max(mds_coord$MDS2)),
    "FALSE",
    as.integer(sum(meta_std$group == "Control_52W")),
    as.integer(sum(meta_std$group == "ACLT_alone_52W"))
  ),
  stringsAsFactors = FALSE
)
write.csv(final_summary, file.path(output_root, "tables", "StepCHRONMDS_01_final_summary_for_review.csv"), row.names = FALSE)

methods_text <- paste(
  "Sample-level PCA and MDS for the chronic pig main comparison were performed using the locked 24-sample synovium gene-level estimated counts matrix, including 12 Control_52W and 12 ACLT_alone_52W samples.",
  "Lowly expressed genes were filtered using edgeR::filterByExpr, followed by TMM normalization with edgeR::calcNormFactors.",
  "Filtered TMM-normalized logCPM values were generated using edgeR::cpm with prior.count = 1.",
  "PCA was performed using prcomp on the transposed logCPM matrix, and MDS was performed using limma::plotMDS with the top 500 genes and common gene selection.",
  "Salmon/tximport estimated counts were used as numeric estimated counts and were not rounded.",
  sep = "\n"
)
writeLines(methods_text, file.path(output_root, "tables", "StepCHRONMDS_01_methods_text_main24_PCA_MDS.txt"))

cat("Final summary for review:\n")
print(final_summary)
cat("\nKey output files:\n")
key_files <- c(
  file.path(output_root, "source_data", "StepCHRONMDS_01_main24_filtered_TMM_logCPM_matrix.csv"),
  file.path(output_root, "source_data", "StepCHRONMDS_01_main24_PCA_coordinates.csv"),
  file.path(output_root, "source_data", "StepCHRONMDS_01_main24_MDS_coordinates.csv"),
  file.path(output_root, "tables", "StepCHRONMDS_01_main24_PCA_variance_explained.csv"),
  file.path(output_root, "tables", "StepCHRONMDS_01_main24_library_summary_after_TMM.csv"),
  fig_files["PCA_unlabeled_pdf"],
  fig_files["PCA_labeled_pdf"],
  fig_files["MDS_unlabeled_pdf"],
  fig_files["MDS_labeled_pdf"],
  file.path(output_root, "tables", "StepCHRONMDS_01_final_summary_for_review.csv")
)
for (i in seq_along(key_files)) cat(i, ") ", key_files[i], "\n", sep = "")

cat("\nSession information:\n")
print(sessionInfo())

status <- "SUCCESS"
