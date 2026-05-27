# Step05: Paired mouse differential expression analysis using limma-voom with duplicateCorrelation.
# Purpose:
# This locked-method Step05 strictly follows the manuscript methods:
# for each mouse timepoint (1W and 4W), it constructs a DGEList from the paired count
# matrix, removes genes containing missing values, applies filterByExpr(), performs
# TMM normalization, runs a first voom(), estimates within-mouse correlation using
# duplicateCorrelation(block = mouse_id), runs a second voom() with the estimated
# correlation, fits the model using lmFit(block = mouse_id, correlation = consensus),
# applies eBayes(), and exports topTable(coef = "treatmentACLR"). logFC > 0 is
# defined as ACLR higher than Contra. This step also generates Figure 1D/1E volcano
# plots with official mouse SYMBOL labels by mapping st_gene_id -> Entrez ID ->
# SYMBOL using clusterProfiler::bitr() and org.Mm.eg.db.

## =========================
## 0. Define paths and locked parameters
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

expr_1w_file <- file.path(base_dir, "02_expression_matrix", "step04_expr_1W.csv")
expr_4w_file <- file.path(base_dir, "02_expression_matrix", "step04_expr_4W.csv")
anno_1w_file <- file.path(base_dir, "01_metadata", "step04_anno_1W.csv")
anno_4w_file <- file.path(base_dir, "01_metadata", "step04_anno_4W.csv")
raw_xlsx <- file.path(base_dir, "00_raw_data", "GSE271903", "GSE271903_count_matrix.xlsx")

de_dir <- file.path(base_dir, "03_DE_analysis")
figure_dir <- file.path(base_dir, "06_figures", "Figure1")
table_dir <- file.path(base_dir, "07_tables", "step05_mouse_DE_volcano")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step05_mouse_DE_volcano"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

strict_FDR <- 0.05
strict_abs_logFC <- 1
top_n_per_direction <- 6

## =========================
## 1. Delete old Step05 outputs
## =========================

old_files <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  file.path(de_dir, "step05_DE_1W_ACLR_vs_Contra_limma_voom.csv"),
  file.path(de_dir, "step05_DE_4W_ACLR_vs_Contra_limma_voom.csv"),
  file.path(de_dir, "step05_voom_logCPM_1W.csv"),
  file.path(de_dir, "step05_voom_logCPM_4W.csv"),
  file.path(de_dir, "step05_design_matrix_1W.csv"),
  file.path(de_dir, "step05_design_matrix_4W.csv"),
  file.path(figure_dir, "Figure1D_mouse_1W_volcano.png"),
  file.path(figure_dir, "Figure1D_mouse_1W_volcano.pdf"),
  file.path(figure_dir, "Figure1E_mouse_4W_volcano.png"),
  file.path(figure_dir, "Figure1E_mouse_4W_volcano.pdf"),
  file.path(figure_dir, "Figure1D_mouse_1W_volcano_final_SYMBOL.png"),
  file.path(figure_dir, "Figure1D_mouse_1W_volcano_final_SYMBOL.pdf"),
  file.path(figure_dir, "Figure1E_mouse_4W_volcano_final_SYMBOL.png"),
  file.path(figure_dir, "Figure1E_mouse_4W_volcano_final_SYMBOL.pdf"),
  file.path(log_dir, paste0(step_name, "_log.txt")),
  file.path(script_dir, "step05_mouse_DE_volcano.R")
)

# Also remove legacy files from the previous numbering (old Step06) to avoid
# mixing old and newly reproduced outputs in the same folder.
legacy_step06_files <- c(
  file.path(de_dir, "step06_DE_1W_ACLR_vs_Contra_limma_voom.csv"),
  file.path(de_dir, "step06_DE_4W_ACLR_vs_Contra_limma_voom.csv"),
  file.path(de_dir, "step06_DE_1W_ACLR_vs_Contra_limma_voom_duplicateCorrelation.csv"),
  file.path(de_dir, "step06_DE_4W_ACLR_vs_Contra_limma_voom_duplicateCorrelation.csv"),
  file.path(de_dir, "step06_voom_logCPM_1W.csv"),
  file.path(de_dir, "step06_voom_logCPM_4W.csv"),
  file.path(de_dir, "step06_design_matrix_1W.csv"),
  file.path(de_dir, "step06_design_matrix_4W.csv"),
  file.path(script_dir, "step06_mouse_DE_volcano.R")
)

unlink(c(old_files, legacy_step06_files)[file.exists(c(old_files, legacy_step06_files))], force = TRUE)

## =========================
## 2. Start log and helper functions
## =========================

sink(log_file, split = TRUE)

cat("===== STEP05 LOCKED METHOD: DUPLICATECORRELATION LIMMA-VOOM =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Base directory: ", base_dir, "\n", sep = "")
cat("Locked method: DGEList -> filterByExpr -> calcNormFactors -> voom1 -> duplicateCorrelation(block=mouse_id) -> voom2 -> lmFit(block=mouse_id, correlation) -> eBayes -> topTable(coef=treatmentACLR)\n\n")

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
    cat("Script archived to: ", path, "\n", sep = "")
  } else {
    writeLines("# Step05 archive fallback.", path, useBytes = TRUE)
    cat("Script path not detected; archive fallback saved.\n")
  }
}

safe_bioc_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

safe_cran_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path, row_names = FALSE) {
  write.csv(x, path, row.names = row_names, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

standardize_sex <- function(x) {
  x <- as.character(x)
  out <- x
  out[toupper(x) == "M"] <- "Male"
  out[toupper(x) == "F"] <- "Female"
  factor(out, levels = c("Female", "Male"))
}

## =========================
## 3. Load packages and archive script
## =========================

safe_bioc_library("edgeR")
safe_bioc_library("limma")
safe_bioc_library("clusterProfiler")
safe_bioc_library("org.Mm.eg.db")
safe_cran_library("readxl")
safe_cran_library("ggplot2")
safe_cran_library("ggrepel")

archive_current_script(file.path(script_dir, "step05_mouse_DE_volcano.R"))

## =========================
## 4. Read count matrices and remove genes with missing values
## =========================

read_count_matrix_locked <- function(path, timepoint) {
  if (!file.exists(path)) stop("Missing count matrix: ", path)

  raw <- read.csv(path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  numeric_df <- as.data.frame(lapply(raw, function(z) suppressWarnings(as.numeric(as.character(z)))))
  rownames(numeric_df) <- rownames(raw)
  mat <- as.matrix(numeric_df)
  storage.mode(mat) <- "numeric"

  na_rows <- rowSums(is.na(mat)) > 0

  na_detail <- data.frame()
  if (any(na_rows)) {
    idx <- which(is.na(mat), arr.ind = TRUE)
    na_detail <- data.frame(
      timepoint = timepoint,
      gene_id = rownames(mat)[idx[, "row"]],
      sample_id = colnames(mat)[idx[, "col"]],
      raw_value = mapply(function(r, c) as.character(raw[r, c]), idx[, "row"], idx[, "col"]),
      stringsAsFactors = FALSE
    )
  }

  mat_clean <- mat[!na_rows, , drop = FALSE]

  if (any(mat_clean < 0, na.rm = TRUE)) stop(timepoint, ": negative counts detected.")

  qc <- data.frame(
    timepoint = timepoint,
    original_genes = nrow(mat),
    samples = ncol(mat),
    genes_removed_due_to_any_NA = sum(na_rows),
    cells_with_NA = sum(is.na(mat)),
    genes_after_NA_removal = nrow(mat_clean),
    non_integer_like_cells = sum(abs(mat_clean - round(mat_clean)) > 1e-6, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(counts = mat_clean, qc = qc, na_detail = na_detail)
}

qc1 <- read_count_matrix_locked(expr_1w_file, "1W")
qc4 <- read_count_matrix_locked(expr_4w_file, "4W")

counts_1w <- qc1$counts
counts_4w <- qc4$counts

anno_1w <- read.csv(anno_1w_file, stringsAsFactors = FALSE)
anno_4w <- read.csv(anno_4w_file, stringsAsFactors = FALSE)

write_csv(rbind(qc1$qc, qc4$qc), file.path(table_dir, "step05_count_matrix_QC_summary.csv"))
if (nrow(qc1$na_detail) > 0) write_csv(qc1$na_detail, file.path(table_dir, "step05_1W_NA_count_cells_removed.csv"))
if (nrow(qc4$na_detail) > 0) write_csv(qc4$na_detail, file.path(table_dir, "step05_4W_NA_count_cells_removed.csv"))

## =========================
## 5. Build st_gene_id -> Entrez -> SYMBOL mapping
## =========================

raw_anno <- readxl::read_excel(raw_xlsx)
raw_anno <- as.data.frame(raw_anno, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("st_gene_id", "gene_id") %in% colnames(raw_anno))) {
  stop("Expected columns st_gene_id and gene_id not found in raw count matrix.")
}

st_to_entrez <- data.frame(
  st_gene_id = clean_text(raw_anno$st_gene_id),
  ENTREZID = clean_text(raw_anno$gene_id),
  stringsAsFactors = FALSE
)
st_to_entrez <- st_to_entrez[st_to_entrez$st_gene_id != "" & st_to_entrez$ENTREZID != "", , drop = FALSE]
st_to_entrez <- st_to_entrez[!duplicated(st_to_entrez$st_gene_id), , drop = FALSE]

bitr_map <- clusterProfiler::bitr(
  unique(st_to_entrez$ENTREZID),
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Mm.eg.db
)
bitr_map <- as.data.frame(bitr_map, stringsAsFactors = FALSE)
bitr_map$ENTREZID <- as.character(bitr_map$ENTREZID)
bitr_map$SYMBOL <- as.character(bitr_map$SYMBOL)
bitr_map <- bitr_map[!duplicated(bitr_map$ENTREZID), , drop = FALSE]

gene_map <- merge(st_to_entrez, bitr_map, by = "ENTREZID", all.x = TRUE)
gene_map$SYMBOL[is.na(gene_map$SYMBOL)] <- ""
gene_map <- gene_map[, c("st_gene_id", "ENTREZID", "SYMBOL")]

write_csv(gene_map, file.path(table_dir, "step05_st_gene_id_entrez_symbol_map.csv"))

## =========================
## 6. Locked duplicateCorrelation DE analysis
## =========================

prepare_annotation <- function(anno, counts, timepoint) {
  anno$sex <- standardize_sex(anno$sex)
  anno$treatment <- factor(anno$treatment, levels = c("Contra", "ACLR"))
  anno$mouse_id <- factor(anno$mouse_id)

  anno <- anno[match(colnames(counts), anno$sample_id), , drop = FALSE]

  if (any(is.na(anno$sample_id))) stop(timepoint, ": annotation does not match count columns.")
  if (!all(c("Contra", "ACLR") %in% anno$treatment)) stop(timepoint, ": missing Contra or ACLR.")
  if (length(unique(anno$mouse_id)) * 2 != nrow(anno)) stop(timepoint, ": paired mouse structure incomplete.")

  anno
}

add_gene_annotation <- function(df) {
  idx <- match(as.character(df$gene_id), gene_map$st_gene_id)
  df$ENTREZID <- gene_map$ENTREZID[idx]
  df$SYMBOL <- gene_map$SYMBOL[idx]
  df$ENTREZID[is.na(df$ENTREZID)] <- ""
  df$SYMBOL[is.na(df$SYMBOL)] <- ""
  df
}

run_locked_de <- function(counts, anno, timepoint, panel, title_text, png_path, pdf_path) {
  anno <- prepare_annotation(anno, counts, timepoint)

  design <- model.matrix(~ treatment, data = anno)
  rownames(design) <- anno$sample_id
  write_csv(as.data.frame(design), file.path(de_dir, paste0("step05_design_matrix_", timepoint, ".csv")), row_names = TRUE)

  dge <- edgeR::DGEList(counts = counts, samples = anno)
  keep <- edgeR::filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  v1 <- limma::voom(dge, design = design, plot = FALSE)
  corfit <- limma::duplicateCorrelation(v1, design = design, block = anno$mouse_id)

  v2 <- limma::voom(dge, design = design, plot = FALSE, block = anno$mouse_id, correlation = corfit$consensus)
  fit <- limma::lmFit(v2, design = design, block = anno$mouse_id, correlation = corfit$consensus)
  fit <- limma::eBayes(fit)

  coef_name <- "treatmentACLR"
  if (!coef_name %in% colnames(fit$coefficients)) {
    stop(timepoint, ": coefficient treatmentACLR not found. Available: ", paste(colnames(fit$coefficients), collapse = ", "))
  }

  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = "BH")
  tt$gene_id <- rownames(tt)
  tt <- tt[, c("gene_id", setdiff(colnames(tt), "gene_id"))]
  colnames(tt)[colnames(tt) == "adj.P.Val"] <- "FDR"

  tt <- add_gene_annotation(tt)
  tt$strict_DEG <- tt$FDR < strict_FDR & abs(tt$logFC) > strict_abs_logFC
  tt$direction <- ifelse(tt$strict_DEG & tt$logFC > 0, "Up",
                         ifelse(tt$strict_DEG & tt$logFC < 0, "Down", "Not_strict"))

  tt$label_final <- ifelse(tt$SYMBOL != "", tt$SYMBOL, tt$gene_id)
  tt$label_score <- -log10(pmax(tt$FDR, .Machine$double.xmin)) * abs(tt$logFC)

  de_path <- file.path(de_dir, paste0("step05_DE_", timepoint, "_ACLR_vs_Contra_limma_voom_duplicateCorrelation.csv"))
  write_csv(tt, de_path)

  # Also write standard filename for downstream steps.
  standard_de_path <- file.path(de_dir, paste0("step05_DE_", timepoint, "_ACLR_vs_Contra_limma_voom.csv"))
  write_csv(tt, standard_de_path)

  voom_path <- file.path(de_dir, paste0("step05_voom_logCPM_", timepoint, ".csv"))
  write.csv(v2$E, voom_path, row.names = TRUE)

  volcano_df <- data.frame(
    panel = panel,
    timepoint = timepoint,
    gene_id = tt$gene_id,
    ENTREZID = tt$ENTREZID,
    SYMBOL = tt$SYMBOL,
    label_final = tt$label_final,
    logFC = tt$logFC,
    AveExpr = tt$AveExpr,
    t = tt$t,
    P.Value = tt$P.Value,
    FDR = tt$FDR,
    B = tt$B,
    strict_DEG = tt$strict_DEG,
    direction = tt$direction,
    label_score = tt$label_score,
    neg_log10_FDR = -log10(pmax(tt$FDR, .Machine$double.xmin)),
    stringsAsFactors = FALSE
  )

  volcano_df$DEG_class <- "Not strict"
  volcano_df$DEG_class[volcano_df$strict_DEG & volcano_df$logFC > 0] <- "Up strict"
  volcano_df$DEG_class[volcano_df$strict_DEG & volcano_df$logFC < 0] <- "Down strict"
  volcano_df$DEG_class <- factor(volcano_df$DEG_class, levels = c("Up strict", "Down strict", "Not strict"))

  volcano_source <- file.path(table_dir, paste0("step05_", panel, "_", timepoint, "_volcano_source_data.csv"))
  write_csv(volcano_df, volcano_source)

  up <- volcano_df[volcano_df$strict_DEG & volcano_df$logFC > 0, , drop = FALSE]
  down <- volcano_df[volcano_df$strict_DEG & volcano_df$logFC < 0, , drop = FALSE]
  up <- up[order(up$SYMBOL == "", -up$label_score), , drop = FALSE]
  down <- down[order(down$SYMBOL == "", -down$label_score), , drop = FALSE]
  labels <- rbind(head(up, top_n_per_direction), head(down, top_n_per_direction))
  labels$label_rule <- paste0("top_", top_n_per_direction, "_per_direction_by_SYMBOL_priority_and_score")

  label_path <- file.path(table_dir, paste0("step05_", panel, "_", timepoint, "_selected_SYMBOL_labels.csv"))
  write_csv(labels, label_path)

  p <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = logFC, y = neg_log10_FDR)) +
    ggplot2::geom_point(ggplot2::aes(color = DEG_class), size = 1.1, alpha = 0.65) +
    ggplot2::geom_vline(xintercept = c(-strict_abs_logFC, strict_abs_logFC), linetype = "dashed", linewidth = 0.3, color = "grey45") +
    ggplot2::geom_hline(yintercept = -log10(strict_FDR), linetype = "dashed", linewidth = 0.3, color = "grey45") +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = label_final),
      size = 3.0,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("Up strict" = "#D55E00", "Down strict" = "#0072B2", "Not strict" = "grey70"),
      drop = FALSE
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = title_text,
      x = "log2 fold change (ACLR vs Contra)",
      y = expression(-log[10]("FDR")),
      color = "DEG class"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(png_path, p, width = 6.5, height = 5.3, dpi = 300)
  ggplot2::ggsave(pdf_path, p, width = 6.5, height = 5.3)

  run_summary <- data.frame(
    panel = panel,
    timepoint = timepoint,
    raw_genes_after_NA_QC = nrow(counts),
    samples = ncol(counts),
    genes_after_filterByExpr = nrow(v2$E),
    design_formula = "~ treatment",
    block_factor = "mouse_id",
    duplicateCorrelation_consensus = corfit$consensus,
    normalization = "edgeR TMM",
    first_voom = TRUE,
    second_voom_with_block_correlation = TRUE,
    lmFit_block = "mouse_id",
    lmFit_correlation = corfit$consensus,
    eBayes_robust = FALSE,
    coefficient = coef_name,
    contrast_direction = "logFC > 0 means ACLR higher than Contra",
    FDR_threshold = strict_FDR,
    abs_logFC_threshold = strict_abs_logFC,
    strict_DEG_total = sum(tt$strict_DEG),
    strict_DEG_up = sum(tt$strict_DEG & tt$logFC > 0),
    strict_DEG_down = sum(tt$strict_DEG & tt$logFC < 0),
    de_table = de_path,
    volcano_source_data = volcano_source,
    selected_label_table = label_path,
    output_png = png_path,
    output_pdf = pdf_path,
    stringsAsFactors = FALSE
  )

  list(de = tt, volcano = volcano_df, labels = labels, run_summary = run_summary)
}

res1 <- run_locked_de(
  counts = counts_1w,
  anno = anno_1w,
  timepoint = "1W",
  panel = "Figure1D",
  title_text = "GSE271903_1W ACLR vs Contra",
  png_path = file.path(figure_dir, "Figure1D_mouse_1W_volcano.png"),
  pdf_path = file.path(figure_dir, "Figure1D_mouse_1W_volcano.pdf")
)

res4 <- run_locked_de(
  counts = counts_4w,
  anno = anno_4w,
  timepoint = "4W",
  panel = "Figure1E",
  title_text = "GSE271903_4W ACLR vs Contra",
  png_path = file.path(figure_dir, "Figure1E_mouse_4W_volcano.png"),
  pdf_path = file.path(figure_dir, "Figure1E_mouse_4W_volcano.pdf")
)

run_summary <- rbind(res1$run_summary, res4$run_summary)
write_csv(run_summary, file.path(table_dir, "step05_DE_run_summary.csv"))

strict_summary <- data.frame(
  timepoint = c("1W", "4W"),
  strict_DEG_total = c(res1$run_summary$strict_DEG_total, res4$run_summary$strict_DEG_total),
  strict_DEG_up = c(res1$run_summary$strict_DEG_up, res4$run_summary$strict_DEG_up),
  strict_DEG_down = c(res1$run_summary$strict_DEG_down, res4$run_summary$strict_DEG_down),
  stringsAsFactors = FALSE
)
write_csv(strict_summary, file.path(table_dir, "step05_strict_DEG_summary.csv"))

software_versions <- data.frame(
  item = c("R", "edgeR", "limma", "readxl", "clusterProfiler", "org.Mm.eg.db", "ggplot2", "ggrepel"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("limma")),
    as.character(utils::packageVersion("readxl")),
    as.character(utils::packageVersion("clusterProfiler")),
    as.character(utils::packageVersion("org.Mm.eg.db")),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("ggrepel"))
  ),
  stringsAsFactors = FALSE
)
write_csv(software_versions, file.path(table_dir, "step05_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step05_sessionInfo.txt"))

cat("\n===== STEP05 SUMMARY =====\n")
cat("\nCount matrix QC summary:\n")
print(rbind(qc1$qc, qc4$qc))
cat("\nDE run summary:\n")
print(run_summary)
cat("\nStrict DEG summary:\n")
print(strict_summary)
cat("\nSelected labels 1W:\n")
print(res1$labels[, c("gene_id", "ENTREZID", "SYMBOL", "label_final", "logFC", "FDR", "direction", "label_score")])
cat("\nSelected labels 4W:\n")
print(res4$labels[, c("gene_id", "ENTREZID", "SYMBOL", "label_final", "logFC", "FDR", "direction", "label_score")])
cat("\nStep05 locked-method completed successfully.\n")

sink()
