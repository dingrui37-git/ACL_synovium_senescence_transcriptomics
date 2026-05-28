# Chronic Supplementary Step5 current78 TMM-aligned: early-defined 24 core genes heatmap in pig early Figure4D style
# Purpose:
#   Visualize the expression pattern of the fixed early-defined 24 core ortholog genes as a supplementary figure
#   in the chronic pig main comparison (Control_52W vs ACLT_alone_52W).
#
# Methodological scale:
#   - Use chronic main-comparison gene-level estimated counts.
#   - Use edgeR TMM-normalized logCPM from the full chronic main-comparison matrix:
#       edgeR::DGEList()
#       edgeR::calcNormFactors(method = "TMM")
#       edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
#   - No filterByExpr() before heatmap extraction, to avoid losing fixed predefined core genes.
#   - Salmon/tximport estimated counts are used as numeric values and are NOT rounded.
#   - Extract the fixed early-defined 24 core genes.
#   - Compute row-wise z-scores across the 24 chronic main-comparison samples.
#
# Visual style update:
#   This version follows the earlier pig Figure4D Step22 heatmap style:
#   - rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red
#   - row clustering enabled with euclidean distance and complete linkage
#   - columns are not clustered and remain Control_52W first, ACLT_alone_52W second
#   - row dendrogram shown on the left
#   - border_color = "#A8A8A8"
#   - expanded output canvas
#   - one PDF and one PNG output
#
# Interpretation:
#   - This is chronic extension visualization only.
#   - This step does NOT redefine a chronic core set.
#   - The displayed genes are the fixed 24 core genes from the pig early primary-analysis workflow.
#   - The plotted heatmap uses z-scores clipped to [-2, 2] for readability;
#     unclipped z-scores are saved as source data.
#
# Output:
#   Written to a separate supplementary output directory. It is not named Figure5C to avoid conflict with the main Figure5C lollipop/dot audit plot.

options(stringsAsFactors = FALSE)

status <- "FAILED"
method_version <- "2026-05-26_supplementary_chronic_step5_TMM_aligned_core_heatmap_pig_early_style_compact_v3"

## =========================
## 0. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_early_dir <- file.path(rebuild_root, "02_pig_early")
pig_chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

chronic_tables_dir <- file.path(pig_chronic_dir, "tables")
chronic_figures_dir <- file.path(
  pig_chronic_dir,
  "figures",
  "Supplementary_chronic_core_heatmap_TMM_aligned_pig_early_style_compact"
)

out_dir <- file.path(
  chronic_tables_dir,
  "supplementary_chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact"
)
log_dir <- file.path(out_dir, "logs")
obj_dir <- file.path(out_dir, "objects")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chronic_figures_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "chronic_step5_TMM_aligned_pig_early_style_core_heatmap_full_log.txt")
summary_log_file <- file.path(log_dir, "chronic_step5_TMM_aligned_pig_early_style_core_heatmap_summary_to_send_me.txt")

zz <- file(full_log_file, open = "wt")
sink(zz, split = TRUE)
sink(zz, type = "message")

close_log <- function() {
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(zz), silent = TRUE)
}

on.exit({
  cat("\n============================================================\n")
  cat("Chronic Step5 TMM-aligned pig-early-style core heatmap finished with status:", status, "\n")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Summary log saved to:\n", summary_log_file, "\n", sep = "")
  cat("============================================================\n")
  close_log()
}, add = TRUE)

cat("===== CHRONIC SUPPLEMENTARY STEP5 CURRENT78 TMM-ALIGNED CORE HEATMAP, PIG EARLY FIGURE4D STYLE =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Method version: ", method_version, "\n", sep = "")
cat("Interpretation: supplementary chronic visualization of fixed early-defined 24 core genes.\n")
cat("This step does NOT redefine a chronic core set.\n\n")

## =========================
## 1. Packages
## =========================

required_pkgs <- c("edgeR", "pheatmap", "RColorBrewer")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("pheatmap", "RColorBrewer")) {
      install.packages(pkg)
    } else {
      stop("Required package not installed: ", pkg, call. = FALSE)
    }
  }
}
suppressPackageStartupMessages({
  library(edgeR)
  library(pheatmap)
  library(RColorBrewer)
})

cat("Required packages loaded.\n")
cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
cat("pheatmap version:", as.character(packageVersion("pheatmap")), "\n")
cat("RColorBrewer version:", as.character(packageVersion("RColorBrewer")), "\n\n")

## =========================
## 2. Helper functions
## =========================

stop_if_missing <- function(file, label) {
  if (is.na(file) || !file.exists(file)) stop("Missing ", label, ": ", file, call. = FALSE)
}

find_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_read_csv <- function(file) {
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(x, file) {
  write.csv(x, file, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
  cat("Saved: ", file, "\n", sep = "")
}

write_matrix_csv <- function(mat, file) {
  write.csv(mat, file, quote = TRUE, fileEncoding = "UTF-8")
  cat("Saved: ", file, "\n", sep = "")
}

clean_string <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}

find_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Could not find ", label, ". Checked: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

normalize_chronic_group <- function(x) {
  z <- tolower(clean_string(x))
  out <- rep(NA_character_, length(z))

  out[grepl("control", z) & grepl("52", z)] <- "Control_52W"
  out[grepl("ctrl", z) & grepl("52", z)] <- "Control_52W"

  out[grepl("aclt", z) & grepl("alone", z) & grepl("52", z)] <- "ACLT_alone_52W"
  out[grepl("acl transection", z) & grepl("alone", z) & grepl("52", z)] <- "ACLT_alone_52W"
  out[grepl("transection alone", z) & grepl("52", z)] <- "ACLT_alone_52W"

  exact <- clean_string(x)
  out[exact %in% c("Control_52W", "ACLT_alone_52W")] <- exact[exact %in% c("Control_52W", "ACLT_alone_52W")]

  out
}

row_zscore_clip <- function(mat) {
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}

row_zscore_unclipped <- function(mat) {
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

file_audit <- function(path, label) {
  info <- if (file.exists(path)) file.info(path) else data.frame(size = NA, mtime = NA)
  md5 <- NA_character_
  if (file.exists(path) && requireNamespace("tools", quietly = TRUE)) {
    md5 <- as.character(tools::md5sum(path))
  }
  data.frame(
    label = label,
    file = path,
    exists = file.exists(path),
    size_bytes = if (file.exists(path)) info$size else NA_real_,
    mtime = if (file.exists(path)) as.character(info$mtime) else NA_character_,
    md5 = md5,
    stringsAsFactors = FALSE
  )
}

find_core_gene_file <- function() {
  candidates <- c(
    file.path(pig_early_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_pig_core_genes.csv"),
    file.path(pig_early_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_pig_core_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_core_genes.csv"),
    file.path(pig_early_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_core_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_audit", "step21_current78_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_audit", "step21_current78_pig_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_pig_early_signature_validation", "step21_current78_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_pig_early_signature_validation", "step21_current78_pig_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_pig_core_ortholog_gene_table.csv"),
    file.path(pig_early_dir, "tables", "step21_current78_pig_core_ortholog_gene_table.csv")
  )
  hit <- find_existing_file(candidates)
  if (!is.na(hit)) return(hit)

  all_csv <- list.files(file.path(pig_early_dir, "tables"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(all_csv) == 0) return(NA_character_)

  bn <- tolower(basename(all_csv))
  score <- rep(0, length(all_csv))
  score[grepl("core", bn)] <- score[grepl("core", bn)] + 60
  score[grepl("ortholog", bn)] <- score[grepl("ortholog", bn)] + 20
  score[grepl("gene", bn)] <- score[grepl("gene", bn)] + 10
  score[grepl("table", bn)] <- score[grepl("table", bn)] + 5
  score[grepl("heatmap|zscore|matrix|summary|score|comparison|logcpm", bn)] <- score[grepl("heatmap|zscore|matrix|summary|score|comparison|logcpm", bn)] - 80

  cand <- data.frame(file = all_csv, score = score, stringsAsFactors = FALSE)
  cand <- cand[order(-cand$score, cand$file), , drop = FALSE]
  write_csv(cand, file.path(out_dir, "chronic_step5_pig_early_style_core_gene_file_candidates.csv"))

  if (nrow(cand) > 0 && cand$score[1] >= 70) return(cand$file[1])
  NA_character_
}

## =========================
## 3. Input files
## =========================

core_gene_file <- find_core_gene_file()

manifest_file <- find_existing_file(c(
  file.path(chronic_tables_dir, "step25v3_pig_chronic_main_comparison_manifest.csv"),
  file.path(chronic_tables_dir, "step25_pig_chronic_main_comparison_manifest.csv"),
  file.path(chronic_tables_dir, "pig_chronic_main_comparison_manifest.csv")
))

count_file <- find_existing_file(c(
  file.path(chronic_tables_dir, "step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv"),
  file.path(chronic_tables_dir, "step25_pig_chronic_main_comparison_gene_level_counts_matrix.csv"),
  file.path(chronic_tables_dir, "pig_chronic_main_comparison_gene_level_counts_matrix.csv")
))

stop_if_missing(core_gene_file, "early-defined 24 core gene table")
stop_if_missing(manifest_file, "chronic main comparison manifest")
stop_if_missing(count_file, "chronic main comparison gene-level count matrix")

cat("Core gene file: ", core_gene_file, "\n", sep = "")
cat("Manifest file: ", manifest_file, "\n", sep = "")
cat("Count matrix file: ", count_file, "\n\n", sep = "")

## =========================
## 4. Read fixed early-defined 24 core genes
## =========================

core_raw <- safe_read_csv(core_gene_file)

pig_ensg_col <- find_col(core_raw, c("pig_ensg", "pig_gene_id", "gene_id", "pig_gene"), TRUE, "core gene pig Ensembl ID")
pig_symbol_col <- find_col(core_raw, c("pig_symbol", "gene_symbol", "symbol", "gene_name"), FALSE, "core gene pig symbol")
human_col <- find_col(core_raw, c("human_gene", "ortholog_name", "human_symbol", "Gene symbol"), FALSE, "core gene human symbol")
mouse_col <- find_col(core_raw, c("mouse_symbol", "input", "mouse_gene"), FALSE, "core gene mouse symbol")
direction_col <- find_col(core_raw, c("signature_direction", "mouse_direction", "direction"), FALSE, "core gene direction")

core_df <- data.frame(
  pig_ensg = clean_string(core_raw[[pig_ensg_col]]),
  pig_symbol = if (!is.na(pig_symbol_col)) clean_string(core_raw[[pig_symbol_col]]) else NA_character_,
  human_gene = if (!is.na(human_col)) clean_string(core_raw[[human_col]]) else NA_character_,
  mouse_symbol = if (!is.na(mouse_col)) clean_string(core_raw[[mouse_col]]) else NA_character_,
  signature_direction = if (!is.na(direction_col)) clean_string(core_raw[[direction_col]]) else NA_character_,
  stringsAsFactors = FALSE
)

core_df <- core_df[!is.na(core_df$pig_ensg) & core_df$pig_ensg != "", , drop = FALSE]
core_df <- core_df[!duplicated(core_df$pig_ensg), , drop = FALSE]

if (nrow(core_df) != 24) {
  stop("Expected 24 early-defined core genes, but detected ", nrow(core_df), " genes in the selected core gene table: ", core_gene_file, call. = FALSE)
}

cat("Fixed early-defined core genes loaded: ", nrow(core_df), "\n\n", sep = "")

## =========================
## 5. Read chronic manifest and count matrix
## =========================

manifest <- safe_read_csv(manifest_file)
count_df <- safe_read_csv(count_file)

sample_col <- find_col(manifest, c("sample_id", "sample", "SampleID", "geo_accession", "GSM", "run", "Run"), TRUE, "chronic manifest sample ID")
group_col <- find_col(manifest, c("core_group", "group", "treatment", "condition", "title", "sample_title"), TRUE, "chronic manifest group")

sample_info <- data.frame(
  sample_id = clean_string(manifest[[sample_col]]),
  group_raw = clean_string(manifest[[group_col]]),
  stringsAsFactors = FALSE
)
sample_info$group <- normalize_chronic_group(sample_info$group_raw)

if ("core_group" %in% colnames(manifest)) {
  cg <- clean_string(manifest$core_group)
  ok <- cg %in% c("Control_52W", "ACLT_alone_52W")
  sample_info$group[ok] <- cg[ok]
}

sample_info <- sample_info[
  !is.na(sample_info$sample_id) &
    sample_info$sample_id != "" &
    sample_info$group %in% c("Control_52W", "ACLT_alone_52W"),
  ,
  drop = FALSE
]
sample_info <- sample_info[!duplicated(sample_info$sample_id), , drop = FALSE]
sample_info$group <- factor(sample_info$group, levels = c("Control_52W", "ACLT_alone_52W"))
sample_info <- sample_info[order(sample_info$group, sample_info$sample_id), , drop = FALSE]

if (sum(sample_info$group == "Control_52W") != 12 ||
    sum(sample_info$group == "ACLT_alone_52W") != 12) {
  stop(
    "Expected 12 Control_52W and 12 ACLT_alone_52W samples. Observed: ",
    sum(sample_info$group == "Control_52W"), " Control_52W and ",
    sum(sample_info$group == "ACLT_alone_52W"), " ACLT_alone_52W.",
    call. = FALSE
  )
}

gene_col <- find_col(count_df, c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id"), TRUE, "count matrix gene ID")
count_df[[gene_col]] <- clean_string(count_df[[gene_col]])

sample_cols <- setdiff(colnames(count_df), gene_col)
numeric_sample_cols <- sample_cols[vapply(count_df[sample_cols], function(x) {
  suppressWarnings(all(is.na(x) | !is.na(as.numeric(as.character(x)))))
}, logical(1))]

missing_samples <- setdiff(sample_info$sample_id, numeric_sample_cols)
if (length(missing_samples) > 0) {
  write_csv(data.frame(missing_sample_id = missing_samples), file.path(out_dir, "chronic_step5_pig_early_style_missing_samples_in_count_matrix.csv"))
  stop("Some chronic samples are missing from count matrix. See diagnostic file.", call. = FALSE)
}

count_mat <- as.matrix(data.frame(
  lapply(count_df[, numeric_sample_cols, drop = FALSE], function(x) as.numeric(as.character(x))),
  check.names = FALSE
))
colnames(count_mat) <- numeric_sample_cols
rownames(count_mat) <- count_df[[gene_col]]

keep_rows <- !is.na(rownames(count_mat)) & rownames(count_mat) != ""
count_mat <- count_mat[keep_rows, , drop = FALSE]

if (any(duplicated(rownames(count_mat)))) {
  count_mat <- rowsum(count_mat, group = rownames(count_mat), reorder = FALSE)
}

count_mat <- count_mat[, sample_info$sample_id, drop = FALSE]

if (anyNA(count_mat)) stop("NA values detected after numeric conversion of chronic count matrix.", call. = FALSE)
if (any(count_mat < 0, na.rm = TRUE)) stop("Negative values detected in chronic count matrix.", call. = FALSE)

## =========================
## 6. Build TMM-normalized logCPM and core heatmap matrices
## =========================

missing_core <- core_df[!(core_df$pig_ensg %in% rownames(count_mat)), , drop = FALSE]
if (nrow(missing_core) > 0) {
  write_csv(missing_core, file.path(out_dir, "chronic_step5_pig_early_style_missing_core_genes_in_count_matrix.csv"))
  stop("Some early-defined core genes are missing in chronic count matrix.", call. = FALSE)
}

dge_full <- edgeR::DGEList(counts = count_mat)
dge_full <- edgeR::calcNormFactors(dge_full, method = "TMM")
full_logcpm_mat <- edgeR::cpm(
  dge_full,
  log = TRUE,
  prior.count = 1,
  normalized.lib.sizes = TRUE
)

normalization_summary <- data.frame(
  metric = c(
    "input_count_matrix_genes",
    "input_count_matrix_samples",
    "edgeR_DGEList_used",
    "TMM_normalization",
    "logCPM_function",
    "prior_count",
    "normalized_lib_sizes",
    "filterByExpr_before_heatmap",
    "estimated_counts_rounded",
    "n_fixed_core_genes",
    "n_fixed_core_genes_detected"
  ),
  value = c(
    nrow(count_mat),
    ncol(count_mat),
    "TRUE",
    "edgeR::calcNormFactors(method = 'TMM')",
    "edgeR::cpm(log = TRUE)",
    "1",
    "TRUE",
    "FALSE",
    "FALSE",
    nrow(core_df),
    nrow(core_df)
  ),
  stringsAsFactors = FALSE
)

lib_summary <- data.frame(
  sample_id = colnames(count_mat),
  group = as.character(sample_info$group),
  library_size = dge_full$samples$lib.size,
  norm_factor = dge_full$samples$norm.factors,
  normalized_library_size = dge_full$samples$lib.size * dge_full$samples$norm.factors,
  detected_genes_count_gt0 = colSums(count_mat > 0),
  stringsAsFactors = FALSE
)

core_logcpm_mat <- full_logcpm_mat[core_df$pig_ensg, sample_info$sample_id, drop = FALSE]
z_mat_unclipped <- row_zscore_unclipped(core_logcpm_mat)
z_mat <- row_zscore_clip(core_logcpm_mat)

# Display label style:
# The previous pig early Figure4D script prioritized pig_symbol, then human_gene.
# This chronic style-aligned script follows that visual convention.
core_df$display_label <- ifelse(
  !is.na(core_df$pig_symbol) & core_df$pig_symbol != "",
  core_df$pig_symbol,
  ifelse(!is.na(core_df$human_gene) & core_df$human_gene != "", core_df$human_gene, core_df$pig_ensg)
)

if (any(duplicated(core_df$display_label))) {
  core_df$display_label <- make.unique(core_df$display_label)
}

rownames(core_logcpm_mat) <- core_df$display_label
rownames(z_mat_unclipped) <- core_df$display_label
rownames(z_mat) <- core_df$display_label
colnames(z_mat) <- sample_info$sample_id

core_df$MouseDirection <- ifelse(
  core_df$signature_direction == "Down_in_ACLR",
  "mouse_persistent_down",
  "mouse_persistent_up"
)
core_df$dir_order <- ifelse(core_df$MouseDirection == "mouse_persistent_down", 1, 2)

# Match previous early-pig visual logic:
# establish a stable biological pre-order, while pheatmap row clustering determines the final dendrogram order.
core_df <- core_df[order(core_df$dir_order, core_df$display_label), , drop = FALSE]
z_mat <- z_mat[core_df$display_label, sample_info$sample_id, drop = FALSE]
z_mat_unclipped <- z_mat_unclipped[core_df$display_label, sample_info$sample_id, drop = FALSE]
core_logcpm_mat <- core_logcpm_mat[core_df$display_label, sample_info$sample_id, drop = FALSE]

## =========================
## 7. Annotations and pig-early Figure4D style plotting
## =========================

annotation_col <- data.frame(
  Group = factor(as.character(sample_info$group), levels = c("Control_52W", "ACLT_alone_52W"))
)
rownames(annotation_col) <- sample_info$sample_id

annotation_row <- data.frame(
  MouseDirection = factor(
    core_df$MouseDirection,
    levels = c("mouse_persistent_down", "mouse_persistent_up")
  )
)
rownames(annotation_row) <- core_df$display_label

# Compact display labels for columns.
# Raw GSM/sample IDs are kept in the source-data files; these short labels are only for plotting.
sample_info$plot_label <- ave(
  as.character(sample_info$group),
  sample_info$group,
  FUN = function(x) {
    prefix <- ifelse(x[1] == "Control_52W", "CON", "ACLT")
    paste0(prefix, seq_along(x))
  }
)
labels_col_compact <- sample_info$plot_label
names(labels_col_compact) <- sample_info$sample_id

heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(101)
heat_breaks <- seq(-2, 2, length.out = length(heat_colors) + 1)

ann_colors <- list(
  Group = c(
    Control_52W = "#66CC00",
    ACLT_alone_52W = "#F564E3"
  ),
  MouseDirection = c(
    mouse_persistent_down = "#C77CFF",
    mouse_persistent_up = "#F8766D"
  )
)

main_title <- sprintf("Supplementary figure: chronic expression of early-defined core genes (%d genes)", nrow(core_df))

pdf_file <- file.path(chronic_figures_dir, "Supplementary_chronic_core_heatmap_TMM_aligned_pig_early_style_compact.pdf")
png_file <- file.path(chronic_figures_dir, "Supplementary_chronic_core_heatmap_TMM_aligned_pig_early_style_compact.png")

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
    labels_col = labels_col_compact[colnames(z_mat)],

    annotation_names_row = FALSE,
    annotation_names_col = TRUE,

    # Compact layout to avoid exceeding figure boundaries.
    cellwidth = 13,
    cellheight = 13,
    fontsize = 8.5,
    fontsize_row = 7.5,
    fontsize_col = 7.2,
    angle_col = 90,

    border_color = "#A8A8A8",
    treeheight_row = 42,
    treeheight_col = 0,
    legend = TRUE,
    main = main_title,

    filename = filename,
    width = width,
    height = height
  )
}

# More compact canvas than the first pig-early-style version.
make_heatmap(filename = pdf_file, width = 9.8, height = 6.8)
make_heatmap(filename = png_file, width = 9.8, height = 6.8)

cat("Saved plot: ", pdf_file, "\n", sep = "")
cat("Saved plot: ", png_file, "\n\n", sep = "")

## =========================
## 8. Save source data and summary
## =========================

core_gene_metadata_file <- file.path(out_dir, "chronic_step5_pig_early_style_fixed_early_defined_24_core_gene_metadata.csv")
all_logcpm_file <- file.path(out_dir, "chronic_step5_pig_early_style_all_gene_TMM_logCPM_matrix.csv")
core_logcpm_file <- file.path(out_dir, "chronic_step5_pig_early_style_core_TMM_logCPM_matrix.csv")
zscore_unclipped_file <- file.path(out_dir, "chronic_step5_pig_early_style_core_heatmap_zscore_matrix_unclipped.csv")
zscore_plot_file <- file.path(out_dir, "chronic_step5_pig_early_style_core_heatmap_zscore_matrix_clipped_for_plot.csv")
sample_metadata_file <- file.path(out_dir, "chronic_step5_pig_early_style_sample_metadata_used.csv")
plot_files_file <- file.path(out_dir, "chronic_step5_pig_early_style_plot_files.csv")
normalization_summary_file <- file.path(out_dir, "chronic_step5_pig_early_style_normalization_summary.csv")
library_summary_file <- file.path(out_dir, "chronic_step5_pig_early_style_TMM_normalization_factors.csv")

core_meta_out <- core_df
core_meta_out$display_label_priority <- "pig_symbol > human_gene > pig_ensg, matching pig early Figure4D visual convention"
core_meta_out$used_for_chronic_step5_TMM_aligned_pig_early_style <- TRUE

write_csv(core_meta_out, core_gene_metadata_file)
write_matrix_csv(full_logcpm_mat, all_logcpm_file)
write_matrix_csv(core_logcpm_mat, core_logcpm_file)
write_matrix_csv(z_mat_unclipped, zscore_unclipped_file)
write_matrix_csv(z_mat, zscore_plot_file)
write_csv(sample_info, sample_metadata_file)
write_csv(normalization_summary, normalization_summary_file)
write_csv(lib_summary, library_summary_file)
write_csv(
  data.frame(
    plot = "Supplementary_chronic_core_heatmap_TMM_aligned_pig_early_style_compact",
    pdf = pdf_file,
    png = png_file,
    stringsAsFactors = FALSE
  ),
  plot_files_file
)

input_audit <- rbind(
  file_audit(core_gene_file, "fixed early-defined 24 core gene table"),
  file_audit(manifest_file, "chronic main comparison manifest"),
  file_audit(count_file, "chronic main comparison gene-level count matrix")
)
input_audit_file <- file.path(out_dir, "chronic_step5_pig_early_style_input_file_audit.csv")
write_csv(input_audit, input_audit_file)

run_summary <- data.frame(
  metric = c(
    "method_version",
    "core_gene_input_file",
    "manifest_file",
    "count_file",
    "n_fixed_early_defined_core_genes",
    "n_core_genes_detected_in_chronic_count_matrix",
    "n_samples",
    "n_Control_52W",
    "n_ACLT_alone_52W",
    "normalization_scope",
    "normalization_method",
    "filterByExpr_before_heatmap",
    "estimated_counts_rounded",
    "heatmap_color_scale",
    "heatmap_value_source_data",
    "heatmap_value_plotted",
    "zscore_clip_for_plot",
    "display_label_priority",
    "row_clustering",
    "column_clustering",
    "row_dendrogram_treeheight",
    "column_order_rule",
    "pdf_file",
    "png_file",
    "output_dir",
    "note"
  ),
  value = c(
    method_version,
    core_gene_file,
    manifest_file,
    count_file,
    nrow(core_df),
    nrow(core_logcpm_mat),
    ncol(core_logcpm_mat),
    sum(sample_info$group == "Control_52W"),
    sum(sample_info$group == "ACLT_alone_52W"),
    "Full chronic main-comparison gene-level estimated count matrix across the 24 plotted samples; 24 core genes are extracted after TMM-normalized logCPM calculation",
    "edgeR::DGEList(full count matrix) -> edgeR::calcNormFactors(method = 'TMM') -> edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)",
    "FALSE",
    "FALSE",
    "rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red",
    "unclipped row-wise z-score across 24 chronic samples",
    "row-wise z-score clipped to [-2, 2] for pheatmap plotting",
    "[-2, 2]",
    "pig_symbol > human_gene > pig_ensg, matching pig early Figure4D visual convention",
    "TRUE; euclidean distance, complete linkage",
    "FALSE",
    "55",
    "Control_52W samples first, then ACLT_alone_52W samples",
    pdf_file,
    png_file,
    out_dir,
    "Supplementary style-aligned chronic Step5 heatmap: TMM normalization is estimated from the full chronic main-comparison count matrix before extracting fixed 24 core genes; pig early Figure4D heatmap style is used. This is not main Figure5C."
  ),
  stringsAsFactors = FALSE
)
run_summary_file <- file.path(out_dir, "chronic_step5_pig_early_style_core_heatmap_run_summary.csv")
write_csv(run_summary, run_summary_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "edgeR_version",
    "pheatmap_version",
    "RColorBrewer_version",
    "fixed_core_definition",
    "heatmap_scale",
    "filtering_before_heatmap",
    "heatmap_value",
    "display_label_priority",
    "visual_style",
    "compact_plot_labels",
    "analysis_role"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("pheatmap")),
    as.character(utils::packageVersion("RColorBrewer")),
    "fixed early-defined 24 core ortholog genes from pig early primary-analysis workflow",
    "edgeR TMM-normalized logCPM",
    "No filterByExpr before heatmap extraction, to retain predefined core genes",
    "row-wise z-score across 24 chronic main-comparison samples; clipped to [-2, 2] for plotting only",
    "pig_symbol > human_gene > pig_ensg",
    "pig early Figure4D Step22 style: rev(RdYlBu), row dendrogram, no column clustering; compact layout",
    "Column display labels shortened to CON1-CON12 and ACLT1-ACLT12 for plotting; raw sample IDs retained in source data",
    "chronic extension visualization only; no chronic core redefinition"
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "chronic_step5_pig_early_style_versions_and_method_records.csv")
write_csv(version_df, version_file)

saveRDS(
  list(
    core_df = core_df,
    sample_info = sample_info,
    count_mat = count_mat,
    dge_full = dge_full,
    full_logcpm_mat = full_logcpm_mat,
    core_logcpm_mat = core_logcpm_mat,
    z_mat_unclipped = z_mat_unclipped,
    z_mat = z_mat,
    run_summary = run_summary,
    version_df = version_df,
    normalization_summary = normalization_summary,
    lib_summary = lib_summary
  ),
  file.path(obj_dir, "chronic_step5_TMM_aligned_pig_early_style_core_heatmap_workspace.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "chronic_step5_pig_early_style_sessionInfo.txt"))

## =========================
## 9. Console and summary-to-send log
## =========================

cat("\n===== Chronic Step5 current78 TMM-aligned pig-early-style core heatmap summary =====\n")
print(run_summary, row.names = FALSE)

cat("\nNormalization summary:\n")
print(normalization_summary, row.names = FALSE)

cat("\nCore genes preview:\n")
print(head(core_meta_out[, c("display_label", "human_gene", "pig_ensg", "pig_symbol", "MouseDirection", "signature_direction")], 24), row.names = FALSE)

cat("\nSample annotation:\n")
print(sample_info, row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== Chronic Step5 current78 TMM-aligned pig-early-style core heatmap SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Run summary:", summary_con)
writeLines(capture.output(print(run_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Normalization summary:", summary_con)
writeLines(capture.output(print(normalization_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Core genes preview:", summary_con)
writeLines(capture.output(print(head(core_meta_out[, c("display_label", "human_gene", "pig_ensg", "pig_symbol", "MouseDirection", "signature_direction")], 24), row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Sample annotation:", summary_con)
writeLines(capture.output(print(sample_info, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This step visualizes the fixed early-defined 24 core genes in chronic samples.", summary_con)
writeLines("It does not redefine a chronic core gene set.", summary_con)
writeLines("Heatmap input is edgeR TMM-normalized logCPM from the full chronic main-comparison count matrix.", summary_con)
writeLines("Heatmap source values are row-wise z-scores computed across 24 chronic samples.", summary_con)
writeLines("The plotted heatmap uses z-scores clipped to [-2, 2] for readability; unclipped z-scores are saved as source data.", summary_con)
writeLines("Visual style follows pig early Figure4D Step22: rev(RdYlBu), row dendrogram, no column clustering, compact canvas; column labels shortened for plotting.", summary_con)
close(summary_con)

status <- "SUCCESS"

cat("\nChronic Step5 current78 TMM-aligned pig-early-style core heatmap completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
