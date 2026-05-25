# Step22_current78: Figure4D core ortholog heatmap in mouse Step09 heatmap style
# Purpose:
# Generate a current78-based pig early core ortholog heatmap using the 24 strict core genes.
# Style update:
#   - use rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red
#   - expand output canvas to prevent row labels / legends from being clipped
#   - generate exactly ONE PDF and ONE PNG, no duplicate canonical copies
#   - corrected normalization: estimate TMM factors from the full Step16 count matrix first,
#     then extract the 24 core genes for heatmap visualization
#   - delete old Step22 outputs before rerun to avoid mixing the previous 24-gene-submatrix TMM output
#
# Inputs:
#   Step18_current78 validation table
#   Step19_current78 score table for sample group annotation
#   Step16 pig early full gene-level count matrix
#   Optional Step11C/current mouse persistent ∩ CellAge overlap table (78-gene current version preferred) for CellAgeEffect annotation
#   If this optional file is not found, the heatmap is still generated with MouseDirection row annotation only.

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

find_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_file_recursive <- function(root_dirs, patterns) {
  root_dirs <- root_dirs[dir.exists(root_dirs)]
  if (length(root_dirs) == 0) return(NA_character_)
  for (root in root_dirs) {
    all_csv <- list.files(root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(all_csv) == 0) next
    for (pat in patterns) {
      hit <- all_csv[grepl(pat, basename(all_csv), ignore.case = TRUE)]
      if (length(hit) > 0) return(hit[1])
    }
  }
  NA_character_
}

find_col <- function(df, candidates, label) {
  idx <- match(candidates, colnames(df))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    stop(sprintf("Could not find column for %s. Checked: %s", label, paste(candidates, collapse = ", ")))
  }
  colnames(df)[idx[1]]
}

clean_char <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

row_zscore <- function(mat) {
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}

# =========================
# 1. Paths
# =========================
base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early"
step18_dir <- file.path(base_dir, "tables", "step18_current78_pig_early_signature_remap")
out_dir <- file.path(base_dir, "tables", "step22_current78_Figure4D_core_heatmap_step09_style")
fig_dir <- file.path(base_dir, "figures", "Figure4")
log_dir <- file.path(out_dir, "logs")

# Clean old Step22 outputs from the previous implementation before regenerating corrected outputs.
# This keeps the Step22 name but prevents mixing the previous 24-core-gene-submatrix TMM output
# with the corrected full-matrix TMM output.
if (dir.exists(out_dir)) {
  unlink(out_dir, recursive = TRUE, force = TRUE)
}
old_fig_files <- file.path(fig_dir, c(
  "Figure4D_current78_core_ortholog_heatmap_step09_style.pdf",
  "Figure4D_current78_core_ortholog_heatmap_step09_style.png"
))
unlink(old_fig_files[file.exists(old_fig_files)], force = TRUE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "step22_current78_Figure4D_core_heatmap_step09_style_full_log.txt")
summary_log_file <- file.path(log_dir, "step22_current78_Figure4D_core_heatmap_step09_style_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)
cat("===== STEP22 CURRENT78 FIGURE4D CORE HEATMAP STEP09 STYLE =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Style: rev(RColorBrewer::RdYlBu), z-score clipped to [-2, 2], expanded canvas.\n")
cat("Normalization correction: TMM factors are estimated from the full Step16 count matrix before extracting the 24 core genes.\n\n")

validation_file <- file.path(step18_dir, "step18_current78_pig_signature_validation_table.csv")
step18_summary_file <- file.path(step18_dir, "step18_current78_pig_early_signature_remap_summary.csv")
step19_score_file <- file.path(base_dir, "tables", "step19_current78_pig_signature_score", "step19_current78_pig_signature_scores_by_sample.csv")

count_file <- find_existing_file(c(
  "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv",
  "E:/R/ACLsenescence2 LD/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv"
))

step11C_overlap_file <- find_existing_file(c(
  # Current locked Step11C file after renaming to the 78-gene version.
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_78_genes.csv",
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_78_genes_compact.csv",
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_78_genes.csv",
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_78_genes_compact.csv",
  # Backward-compatible fallback for older archived Step11C file names.
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_79_genes.csv",
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_79_genes_compact.csv",
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_79_genes.csv",
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_79_genes_compact.csv"
))
if (is.na(step11C_overlap_file)) {
  step11C_overlap_file <- find_file_recursive(
    root_dirs = c(
      "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery",
      "E:/R/ACLsenescence2/rebuild_submission",
      "E:/R/ACLsenescence2 LD/rebuild_submission/02_mouse_discovery"
    ),
    patterns = c("step11C.*CellAge.*overlap.*\\.csv$", "persistent.*CellAge.*overlap.*\\.csv$")
  )
}

if (!file.exists(validation_file)) stop("Missing Step18 validation file: ", validation_file)
if (!file.exists(step18_summary_file)) stop("Missing Step18 summary file: ", step18_summary_file)
if (!file.exists(step19_score_file)) stop("Missing Step19 score file: ", step19_score_file)
if (is.na(count_file)) stop("Missing Step16 pig early gene count matrix.")

use_cellage_effect_annotation <- !is.na(step11C_overlap_file) && file.exists(step11C_overlap_file)

cat("Validation file: ", validation_file, "\n", sep = "")
cat("Step18 summary file: ", step18_summary_file, "\n", sep = "")
cat("Step19 score file: ", step19_score_file, "\n", sep = "")
cat("Count file: ", count_file, "\n", sep = "")
if (use_cellage_effect_annotation) {
  cat("Optional Step11C overlap file found for CellAgeEffect annotation: ", step11C_overlap_file, "\n\n", sep = "")
} else {
  cat("Optional Step11C overlap file was not found. The heatmap will be generated without CellAgeEffect row annotation.\n\n")
}

# =========================
# 2. Read inputs and identify current78 core genes
# =========================
val_df <- read_csv(validation_file, show_col_types = FALSE)
step18_summary <- read_csv(step18_summary_file, show_col_types = FALSE)
score_df <- read_csv(step19_score_file, show_col_types = FALSE)
overlap_df <- if (use_cellage_effect_annotation) read_csv(step11C_overlap_file, show_col_types = FALSE) else NULL
count_df <- read_csv(count_file, show_col_types = FALSE)

pig_ensg_col <- find_col(val_df, c("pig_ensg", "pig_gene_id", "pig_gene", "ensembl_gene_id"), "pig ensg")
pig_symbol_col <- find_col(val_df, c("pig_symbol", "pig_gene", "gene_symbol", "symbol"), "pig symbol")
human_gene_col <- find_col(val_df, c("human_gene", "ortholog_name", "Gene symbol", "gene_symbol"), "human gene")
mouse_symbol_col <- find_col(val_df, c("mouse_symbol", "input", "mouse_gene", "mouse_gene_symbol"), "mouse symbol")
sig_dir_col <- find_col(val_df, c("signature_direction", "direction", "mouse_direction"), "signature direction")
dir_both_col <- find_col(val_df, c("direction_consistent_both_timepoints", "direction_consistent_both", "direction_consistent"), "direction consistent both")
strict_t7_col <- find_col(val_df, c("strict_t7", "strict_at_t7"), "strict_t7")
strict_t28_col <- find_col(val_df, c("strict_t28", "strict_at_t28"), "strict_t28")
logFC_t7_col <- find_col(val_df, c("logFC_t7", "pig_logFC_t7", "t7_logFC"), "logFC_t7")
FDR_t7_col <- find_col(val_df, c("FDR_t7", "pig_FDR_t7", "t7_FDR"), "FDR_t7")
logFC_t28_col <- find_col(val_df, c("logFC_t28", "pig_logFC_t28", "t28_logFC"), "logFC_t28")
FDR_t28_col <- find_col(val_df, c("FDR_t28", "pig_FDR_t28", "t28_FDR"), "FDR_t28")

core_df <- val_df %>%
  mutate(
    pig_ensg = .data[[pig_ensg_col]],
    pig_symbol = clean_char(.data[[pig_symbol_col]]),
    human_gene = clean_char(.data[[human_gene_col]]),
    mouse_symbol = clean_char(.data[[mouse_symbol_col]]),
    signature_direction = clean_char(.data[[sig_dir_col]]),
    direction_consistent_both_timepoints = as.logical(.data[[dir_both_col]]),
    strict_t7 = as.logical(.data[[strict_t7_col]]),
    strict_t28 = as.logical(.data[[strict_t28_col]]),
    logFC_t7 = as.numeric(.data[[logFC_t7_col]]),
    FDR_t7 = as.numeric(.data[[FDR_t7_col]]),
    logFC_t28 = as.numeric(.data[[logFC_t28_col]]),
    FDR_t28 = as.numeric(.data[[FDR_t28_col]])
  ) %>%
  filter(direction_consistent_both_timepoints, strict_t7, strict_t28) %>%
  distinct(pig_ensg, .keep_all = TRUE)

if (nrow(core_df) == 0) stop("No core genes identified.")
cat("Core genes identified: ", nrow(core_df), "\n", sep = "")

# Optional CellAge effect annotation from Step11C overlap table.
# This annotation is useful but not required for Figure4D generation.
# If Step11C is unavailable, continue with MouseDirection-only row annotation.
if (use_cellage_effect_annotation) {
  ov_human_col <- find_col(overlap_df, c("ortholog_name", "Gene symbol", "human_gene", "gene_symbol"), "Step11C human gene")
  possible_effect_cols <- intersect(c("Senescence Effect", "Senescence_Effect", "senescence_effect", "CellAgeEffect"), colnames(overlap_df))
  if (length(possible_effect_cols) == 0) {
    warning("Step11C file was found, but no CellAge effect column was found. Continuing without CellAgeEffect row annotation.")
    use_cellage_effect_annotation <- FALSE
  } else {
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
  }
}
if (!use_cellage_effect_annotation) {
  core_df$CellAgeEffect <- "Not_annotated"
}

# =========================
# 3. Expression matrix and sample annotation
# =========================
count_gene_col <- find_col(count_df, c("gene_id", "Geneid", "ensembl_gene_id"), "count gene id")
non_count_cols <- c(count_gene_col, intersect(c("gene_symbol", "symbol", "gene_name", "gene", "Length", "length", "Chr", "Start", "End", "Strand"), colnames(count_df)))
sample_cols <- setdiff(colnames(count_df), non_count_cols)

count_mat <- as.matrix(count_df[, sample_cols, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- count_df[[count_gene_col]]

missing_core <- setdiff(core_df$pig_ensg, rownames(count_mat))
if (length(missing_core) > 0) {
  write_csv(data.frame(missing_pig_ensg = missing_core), file.path(out_dir, "step22_current78_missing_core_genes_in_count_matrix.csv"))
  stop("Some core genes are missing in the count matrix.")
}

score_sample_col <- find_col(score_df, c("sample_id", "sample", "SampleID", "sample_name"), "Step19 sample id")
score_group_col <- find_col(score_df, c("group", "condition", "time_group", "treatment_group", "group_label"), "Step19 group")

annotation_df <- score_df %>%
  transmute(sample_id = .data[[score_sample_col]], group = clean_char(.data[[score_group_col]])) %>%
  distinct()

normalize_group <- function(x) {
  lx <- tolower(as.character(x))
  if (grepl("control|ctrl|con|t0", lx)) return("Control")
  if (grepl("t7|1w|week1|day7", lx)) return("ACLT-untreated-1W")
  if (grepl("t28|4w|week4|day28", lx)) return("ACLT-untreated-4W")
  as.character(x)
}
annotation_df$group <- vapply(annotation_df$group, normalize_group, character(1))
annotation_df$group <- factor(annotation_df$group, levels = c("Control", "ACLT-untreated-1W", "ACLT-untreated-4W"))
annotation_df <- annotation_df %>% arrange(group, sample_id)

if (!all(annotation_df$sample_id %in% colnames(count_mat))) {
  missing_samples <- setdiff(annotation_df$sample_id, colnames(count_mat))
  write_csv(data.frame(missing_sample_id = missing_samples), file.path(out_dir, "step22_current78_missing_samples_in_count_matrix.csv"))
  stop("Some Step19 samples are missing in the count matrix.")
}

# Corrected normalization workflow:
#   1) use the full Step16 gene-level count matrix for all 18 plotted samples;
#   2) estimate TMM normalization factors from the full matrix;
#   3) compute TMM-normalized logCPM genome-wide;
#   4) extract the 24 core genes for heatmap z-score visualization.
full_counts_for_tmm <- count_mat[, annotation_df$sample_id, drop = FALSE]
if (any(!is.finite(full_counts_for_tmm))) stop("Non-finite values found in full count matrix.")
if (any(full_counts_for_tmm < 0, na.rm = TRUE)) stop("Negative values found in full count matrix.")

dge_full <- DGEList(counts = full_counts_for_tmm)
dge_full <- calcNormFactors(dge_full, method = "TMM")
full_logcpm_mat <- cpm(dge_full, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

core_logcpm_mat <- full_logcpm_mat[core_df$pig_ensg, annotation_df$sample_id, drop = FALSE]
z_mat <- row_zscore(core_logcpm_mat)

core_df$display_label <- ifelse(!is.na(core_df$pig_symbol) & core_df$pig_symbol != "", core_df$pig_symbol, core_df$human_gene)
rownames(core_logcpm_mat) <- core_df$display_label
rownames(z_mat) <- core_df$display_label
colnames(z_mat) <- annotation_df$sample_id

core_df$MouseDirection <- ifelse(core_df$signature_direction == "Down_in_ACLR", "mouse_persistent_down", "mouse_persistent_up")
core_df$dir_order <- ifelse(core_df$MouseDirection == "mouse_persistent_down", 1, 2)
core_df$effect_order <- ifelse(core_df$CellAgeEffect == "Induces", 1, ifelse(core_df$CellAgeEffect == "Inhibits", 2, 3))

# Keep a stable biological order: down first, then up; if CellAgeEffect is available, Induces before Inhibits.
# Row clustering is still enabled in pheatmap, so final visual order follows dendrogram, like the reference.
core_df <- core_df %>% arrange(dir_order, effect_order, display_label)
z_mat <- z_mat[core_df$display_label, annotation_df$sample_id, drop = FALSE]

annotation_col <- data.frame(Group = annotation_df$group)
rownames(annotation_col) <- annotation_df$sample_id

if (use_cellage_effect_annotation) {
  annotation_row <- data.frame(
    CellAgeEffect = factor(core_df$CellAgeEffect, levels = c("Induces", "Inhibits", "Unclear")),
    MouseDirection = factor(core_df$MouseDirection, levels = c("mouse_persistent_down", "mouse_persistent_up"))
  )
} else {
  annotation_row <- data.frame(
    MouseDirection = factor(core_df$MouseDirection, levels = c("mouse_persistent_down", "mouse_persistent_up"))
  )
}
rownames(annotation_row) <- core_df$display_label

# =========================
# 4. Heatmap style from mouse Step09
# =========================
heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(101)
heat_breaks <- seq(-2, 2, length.out = length(heat_colors) + 1)

ann_colors <- list(
  Group = c(
    Control = "#66CC00",
    "ACLT-untreated-1W" = "#F564E3",
    "ACLT-untreated-4W" = "#00BFC4"
  ),
  MouseDirection = c(
    mouse_persistent_down = "#C77CFF",
    mouse_persistent_up = "#F8766D"
  )
)
if (use_cellage_effect_annotation) {
  ann_colors$CellAgeEffect <- c(
    Induces = "#D9B300",
    Inhibits = "#00C1A2",
    Unclear = "#BDBDBD"
  )
}

main_title <- sprintf("Figure 4D: core ortholog heatmap (current78 %d genes)", nrow(core_df))

pdf_file <- file.path(fig_dir, "Figure4D_current78_core_ortholog_heatmap_step09_style.pdf")
png_file <- file.path(fig_dir, "Figure4D_current78_core_ortholog_heatmap_step09_style.png")

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

# Wide canvas prevents right-side row labels and legends from being clipped.
# Only one PDF and one PNG are generated.
make_heatmap(filename = pdf_file, width = 13.5, height = 8.2)
make_heatmap(filename = png_file, width = 13.5, height = 8.2)

# =========================
# 5. Save source data and summary
# =========================
write_csv(
  core_df %>%
    select(human_gene, mouse_symbol, pig_ensg, pig_symbol, display_label, CellAgeEffect, MouseDirection,
           signature_direction, logFC_t7, FDR_t7, logFC_t28, FDR_t28),
  file.path(out_dir, "step22_current78_Figure4D_core_ortholog_gene_table.csv")
)

write_csv(
  as.data.frame(z_mat) %>% mutate(display_label = rownames(z_mat), .before = 1),
  file.path(out_dir, "step22_current78_Figure4D_core_heatmap_zscore_matrix.csv")
)

write_csv(
  as.data.frame(core_logcpm_mat) %>% mutate(display_label = rownames(core_logcpm_mat), .before = 1),
  file.path(out_dir, "step22_current78_Figure4D_core_TMM_logCPM_matrix.csv")
)

write_csv(
  data.frame(
    sample_id = colnames(dge_full$counts),
    library_size = dge_full$samples$lib.size,
    norm_factor = dge_full$samples$norm.factors,
    normalized_library_size = dge_full$samples$lib.size * dge_full$samples$norm.factors,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "step22_current78_Figure4D_full_matrix_TMM_normalization_factors.csv")
)

write_csv(annotation_df, file.path(out_dir, "step22_current78_Figure4D_core_sample_annotation.csv"))

summary_df <- data.frame(
  metric = c(
    "current78_core_genes",
    "control_n",
    "aclt_1w_n",
    "aclt_4w_n",
    "normalization_scope",
    "normalization_method",
    "CellAgeEffect_annotation",
    "heatmap_color_scale",
    "zscore_clip",
    "pdf_file",
    "png_file",
    "output_dir",
    "note"
  ),
  value = c(
    nrow(core_df),
    sum(annotation_df$group == "Control"),
    sum(annotation_df$group == "ACLT-untreated-1W"),
    sum(annotation_df$group == "ACLT-untreated-4W"),
    "Full Step16 gene-level count matrix across the 18 plotted samples; 24 core genes are extracted after TMM-normalized logCPM calculation",
    "edgeR::DGEList(full count matrix) -> edgeR::calcNormFactors(method = 'TMM') -> edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)",
    ifelse(use_cellage_effect_annotation, paste0("Included from optional Step11C file: ", step11C_overlap_file), "Not included because no Step11C/CellAge overlap file was found; heatmap generated with MouseDirection annotation only"),
    "rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red",
    "[-2, 2]",
    pdf_file,
    png_file,
    out_dir,
    "Corrected Step22: TMM normalization is estimated from the full Step16 count matrix before extracting 24 core genes; one PDF and one PNG are generated."
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_df, file.path(out_dir, "step22_current78_Figure4D_core_heatmap_step09_style_summary.csv"))

cat("\n===== STEP22 CURRENT78 FIGURE4D SUMMARY =====\n")
print(summary_df)

cat("\nCore genes preview:\n")
print(head(core_df %>% select(display_label, human_gene, pig_ensg, CellAgeEffect, MouseDirection, logFC_t7, FDR_t7, logFC_t28, FDR_t28), 20))

cat("\nSample annotation:\n")
print(annotation_df)

cat("\nNormalization check: full Step16 count matrix was used to estimate TMM factors before extracting core genes.\n")
cat("Step22 current78 Figure4D core heatmap completed successfully.\n")

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== STEP22 CURRENT78 FIGURE4D SUMMARY TO SEND ME =====", summary_con)
writeLines("Correction: TMM normalization was estimated from the full Step16 count matrix before extracting 24 core genes.", summary_con)
writeLines(ifelse(use_cellage_effect_annotation, paste0("CellAgeEffect annotation included from: ", step11C_overlap_file), "CellAgeEffect annotation was not included because the optional Step11C/CellAge overlap file was not found."), summary_con)
writeLines("", summary_con)
writeLines(capture.output(print(summary_df)), summary_con)
writeLines("", summary_con)
writeLines("Core genes preview:", summary_con)
writeLines(capture.output(print(head(core_df %>% select(display_label, human_gene, pig_ensg, CellAgeEffect, MouseDirection, logFC_t7, FDR_t7, logFC_t28, FDR_t28), 20))), summary_con)
writeLines("", summary_con)
writeLines("Sample annotation:", summary_con)
writeLines(capture.output(print(annotation_df)), summary_con)
close(summary_con)

sink()
cat("\nStep22 current78 Figure4D core heatmap completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
