# Step09: Redraw Figure 2C persistent-gene heatmap in pheatmap style with expanded right margin.
# Purpose:
# This corrected Step09 does not rerun DE. It uses locked Step07 persistent genes and
# locked Step05 second voom logCPM matrices, then redraws Figure 2C in the requested
# pheatmap style. This version increases output canvas width and reduces annotation
# text size to prevent the right-side legends and row labels from being clipped.

options(stringsAsFactors = FALSE)

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

persistent_file <- file.path(
  base_dir, "07_tables", "step07_strict_DEG_upset_persistent",
  "step07_persistent_direction_consistent_genes.csv"
)

voom_1w_file <- file.path(base_dir, "03_DE_analysis", "step05_voom_logCPM_1W.csv")
voom_4w_file <- file.path(base_dir, "03_DE_analysis", "step05_voom_logCPM_4W.csv")

anno_1w_file <- file.path(base_dir, "01_metadata", "step04_anno_1W.csv")
anno_4w_file <- file.path(base_dir, "01_metadata", "step04_anno_4W.csv")

figure_dir <- file.path(base_dir, "06_figures", "Figure2")
table_dir <- file.path(base_dir, "07_tables", "step09_persistent_heatmap")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step09_persistent_heatmap"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

top_up_n <- 25
top_down_n <- 25

old_files <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  file.path(figure_dir, "Figure2C_mouse_top_persistent_genes_heatmap.png"),
  file.path(figure_dir, "Figure2C_mouse_top_persistent_genes_heatmap.pdf"),
  file.path(log_dir, paste0(step_name, "_log.txt")),
  file.path(script_dir, "step09_persistent_heatmap.R")
)
unlink(old_files[file.exists(old_files)], force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP09 PHEATMAP STYLE HEATMAP WIDTH-FIXED =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Fix: expanded output width and reduced font sizes to avoid right-side clipping.\n\n")

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  NA_character_
}

archive_current_script <- function(path) {
  src <- get_script_path()
  if (!is.na(src) && file.exists(src)) {
    writeLines(readLines(src, warn = FALSE, encoding = "UTF-8"), path, useBytes = TRUE)
  } else {
    writeLines("# Step09 archive fallback.", path, useBytes = TRUE)
  }
}

safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path, row_names = FALSE) {
  write.csv(x, path, row.names = row_names, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

make_unique_labels <- function(symbol, gene_id) {
  label <- ifelse(!is.na(symbol) & symbol != "", symbol, gene_id)
  duplicated_label <- duplicated(label) | duplicated(label, fromLast = TRUE)
  if (any(duplicated_label)) {
    label[duplicated_label] <- paste0(label[duplicated_label], "_", gene_id[duplicated_label])
  }
  label
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

safe_library("pheatmap")
safe_library("RColorBrewer")

archive_current_script(file.path(script_dir, "step09_persistent_heatmap.R"))

if (!file.exists(persistent_file)) stop("Missing Step07 persistent table: ", persistent_file)
if (!file.exists(voom_1w_file)) stop("Missing Step05 1W voom matrix: ", voom_1w_file)
if (!file.exists(voom_4w_file)) stop("Missing Step05 4W voom matrix: ", voom_4w_file)

persistent <- read.csv(persistent_file, stringsAsFactors = FALSE, check.names = FALSE)

voom_1w <- read.csv(voom_1w_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
voom_4w <- read.csv(voom_4w_file, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)

voom_1w <- as.matrix(voom_1w)
voom_4w <- as.matrix(voom_4w)
storage.mode(voom_1w) <- "numeric"
storage.mode(voom_4w) <- "numeric"

anno_1w <- read.csv(anno_1w_file, stringsAsFactors = FALSE)
anno_4w <- read.csv(anno_4w_file, stringsAsFactors = FALSE)

anno_1w$Timepoint <- "1W"
anno_4w$Timepoint <- "4W"

anno <- rbind(anno_1w, anno_4w)
anno$Treatment <- anno$treatment
anno$Timepoint <- factor(anno$Timepoint, levels = c("1W", "4W"))
anno$Treatment <- factor(anno$Treatment, levels = c("Contra", "ACLR"))

persistent$combined_score <- -log10(pmax(persistent$FDR_1W, .Machine$double.xmin)) * abs(persistent$logFC_1W) +
  -log10(pmax(persistent$FDR_4W, .Machine$double.xmin)) * abs(persistent$logFC_4W)

persistent_up <- persistent[persistent$direction_persistent == "Up", , drop = FALSE]
persistent_down <- persistent[persistent$direction_persistent == "Down", , drop = FALSE]

persistent_up <- persistent_up[order(-persistent_up$combined_score), , drop = FALSE]
persistent_down <- persistent_down[order(-persistent_down$combined_score), , drop = FALSE]

selected <- rbind(head(persistent_up, top_up_n), head(persistent_down, top_down_n))
selected <- selected[!duplicated(selected$gene_id), , drop = FALSE]
selected$Direction <- ifelse(selected$direction_persistent == "Up", "Persistent_up", "Persistent_down")
selected$Direction <- factor(selected$Direction, levels = c("Persistent_up", "Persistent_down"))
selected$heatmap_label <- make_unique_labels(selected$SYMBOL, selected$gene_id)

write_csv(selected, file.path(table_dir, "step09_selected_top_persistent_genes_pheatmap_style.csv"))

common_genes <- intersect(selected$gene_id, intersect(rownames(voom_1w), rownames(voom_4w)))
if (length(common_genes) == 0) stop("No selected persistent genes found in both Step05 voom matrices.")

selected <- selected[match(common_genes, selected$gene_id), , drop = FALSE]

expr_mat <- cbind(
  voom_1w[common_genes, , drop = FALSE],
  voom_4w[common_genes, , drop = FALSE]
)

sample_order_df <- anno[anno$sample_id %in% colnames(expr_mat), , drop = FALSE]
sample_order_df$Treatment <- factor(sample_order_df$Treatment, levels = c("Contra", "ACLR"))
sample_order_df$Timepoint <- factor(sample_order_df$Timepoint, levels = c("1W", "4W"))
sample_order_df <- sample_order_df[order(sample_order_df$Timepoint, sample_order_df$Treatment, sample_order_df$mouse_id, sample_order_df$sample_id), , drop = FALSE]

sample_order <- sample_order_df$sample_id
expr_mat <- expr_mat[, sample_order, drop = FALSE]
rownames(expr_mat) <- selected$heatmap_label

z_mat <- row_zscore(expr_mat)
z_mat[z_mat > 2] <- 2
z_mat[z_mat < -2] <- -2

annotation_col <- data.frame(
  Treatment = as.character(sample_order_df$Treatment),
  Timepoint = as.character(sample_order_df$Timepoint),
  row.names = sample_order_df$sample_id,
  stringsAsFactors = FALSE
)

annotation_row <- data.frame(
  Direction = as.character(selected$Direction),
  row.names = selected$heatmap_label,
  stringsAsFactors = FALSE
)

write_csv(as.data.frame(expr_mat), file.path(table_dir, "step09_heatmap_voom_logCPM_matrix_pheatmap_style.csv"), row_names = TRUE)
write_csv(as.data.frame(z_mat), file.path(table_dir, "step09_heatmap_row_zscore_matrix_pheatmap_style.csv"), row_names = TRUE)
write_csv(annotation_col, file.path(table_dir, "step09_heatmap_column_annotation_pheatmap_style.csv"), row_names = TRUE)
write_csv(annotation_row, file.path(table_dir, "step09_heatmap_row_annotation_pheatmap_style.csv"), row_names = TRUE)
write_csv(sample_order_df, file.path(table_dir, "step09_heatmap_sample_order_pheatmap_style.csv"))

ann_colors <- list(
  Treatment = c(
    Contra = "#00BA38",
    ACLR   = "#619CFF"
  ),
  Timepoint = c(
    "1W" = "#D9B300",
    "4W" = "#00BFC4"
  ),
  Direction = c(
    Persistent_down = "#D98CFF",
    Persistent_up   = "#F8766D"
  )
)

heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(101)
heat_breaks <- seq(-2, 2, length.out = length(heat_colors) + 1)

main_title <- "Top strict persistent genes across 1W and 4W"

pdf_file <- file.path(figure_dir, "Figure2C_mouse_top_persistent_genes_heatmap.pdf")
png_file <- file.path(figure_dir, "Figure2C_mouse_top_persistent_genes_heatmap.png")

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
    show_colnames = FALSE,
    labels_row = rownames(z_mat),

    annotation_names_row = FALSE,
    annotation_names_col = TRUE,

    cellwidth = 12,
    cellheight = 10,
    fontsize = 8,
    fontsize_row = 6,
    fontsize_col = 8,
    angle_col = 45,

    border_color = NA,
    treeheight_row = 55,
    treeheight_col = 0,
    legend = TRUE,
    main = main_title,

    filename = filename,
    width = width,
    height = height
  )
}

# Wider canvas prevents right-side row labels and legends from being clipped.
make_heatmap(filename = pdf_file, width = 11.5, height = 9.2)
make_heatmap(filename = png_file, width = 11.5, height = 9.2)

summary_df <- data.frame(
  metric = c(
    "persistent_total_input",
    "persistent_up_input",
    "persistent_down_input",
    "selected_top_up",
    "selected_top_down",
    "selected_total",
    "genes_found_in_both_voom_matrices",
    "samples_in_heatmap",
    "column_order",
    "top_annotations",
    "row_annotation",
    "row_dendrogram",
    "column_clustering",
    "heatmap_color_scale",
    "zscore_clip",
    "figure_width_in",
    "figure_height_in",
    "heatmap_input"
  ),
  value = c(
    nrow(persistent),
    nrow(persistent_up),
    nrow(persistent_down),
    sum(selected$direction_persistent == "Up"),
    sum(selected$direction_persistent == "Down"),
    nrow(selected),
    length(common_genes),
    ncol(z_mat),
    "1W Contra -> 1W ACLR -> 4W Contra -> 4W ACLR",
    "Treatment and Timepoint only",
    "Direction shown as annotation legend; annotation name hidden",
    "TRUE, euclidean complete clustering",
    "FALSE",
    "rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red",
    "[-2, 2]",
    "11.5",
    "9.2",
    "Step05 second voom logCPM, row z-score"
  ),
  stringsAsFactors = FALSE
)

write_csv(summary_df, file.path(table_dir, "step09_persistent_heatmap_summary.csv"))

software_versions <- data.frame(
  item = c("R", "pheatmap", "RColorBrewer"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("pheatmap")),
    as.character(utils::packageVersion("RColorBrewer"))
  ),
  stringsAsFactors = FALSE
)
write_csv(software_versions, file.path(table_dir, "step09_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step09_sessionInfo.txt"))

cat("\n===== STEP09 SUMMARY =====\n")
print(summary_df)

cat("\nSelected gene preview:\n")
print(selected[, c("gene_id", "ENTREZID", "SYMBOL", "heatmap_label", "Direction", "logFC_1W", "logFC_4W", "combined_score")])

cat("\nColumn order preview:\n")
print(sample_order_df[, c("sample_id", "Timepoint", "Treatment", "mouse_id", "sex")])

cat("\nFigure2C saved to:\n")
cat(pdf_file, "\n")
cat(png_file, "\n")

cat("\nStep09 completed successfully.\n")

sink()
