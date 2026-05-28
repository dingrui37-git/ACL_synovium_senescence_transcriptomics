# Chronic Step4 current78: edgeR QLF differential expression and full Hallmark GSEA
# Purpose:
#   Re-run chronic pig Control_52W vs ACLT_alone_52W differential expression using
#   edgeR QL with tximport effective-length offset, then perform full Hallmark GSEA.
#
# Important interpretation:
#   - Chronic pig analysis is an extension validation.
#   - This step does NOT redefine the current78-derived signature.
#   - This step does NOT redefine core ortholog genes.
#   - DE/GSEA is genome-wide and independent from the fixed 75-gene signature score.
#
# Current locked GSEA口径:
#   - full Sus scrofa MSigDB Hallmark collection
#   - rank statistic = sign(logFC) * -log10(PValue)
#   - pig gene symbols as ranking identifiers
#   - duplicate gene symbols resolved by largest absolute rank statistic
#   - fgsea parameters: minSize = 10, maxSize = 500, eps = 0, nproc = 1, seed = 1
#   - target pathways are extracted after full Hallmark-level analysis
#
# This script does NOT auto-archive or overwrite any manually saved R script.
# Revision in this version:
#   Step4 QC MDS is drawn inside this same Step4 script as a colored-point ggplot
#   instead of the old overlapping text-label limma::plotMDS graphic.

options(stringsAsFactors = FALSE)

method_version <- "2026-05-27_chronic_step4_current78_DE_GSEA_with_integrated_colored_point_MDS_v1"

## =========================
## 0. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_chronic_dir <- file.path(rebuild_root, "03_pig_chronic")
chronic_tables_dir <- file.path(pig_chronic_dir, "tables")
chronic_figures_dir <- file.path(pig_chronic_dir, "figures", "Figure5B_current78_chronic_Hallmark_GSEA")

out_dir <- file.path(chronic_tables_dir, "chronic_step4_current78_DE_GSEA")
log_dir <- file.path(out_dir, "logs")
obj_dir <- file.path(out_dir, "objects")
qc_fig_dir <- file.path(pig_chronic_dir, "figures", "chronic_step4_current78_QC")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chronic_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_fig_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "chronic_step4_current78_DE_GSEA_full_log.txt")
summary_log_file <- file.path(log_dir, "chronic_step4_current78_DE_GSEA_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== CHRONIC STEP4 CURRENT78 DE + FULL HALLMARK GSEA =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Comparison: ACLT_alone_52W vs Control_52W.\n")
cat("Method: edgeR QL with tximport effective-length offset; full Sus scrofa Hallmark fgsea.\n")
cat("Method version: ", method_version, "\n\n", sep = "")

## =========================
## 1. Packages
## =========================

required_pkgs <- c("edgeR", "limma", "fgsea", "msigdbr", "ggplot2", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
})

## =========================
## 2. Helper functions
## =========================

stop_if_missing <- function(file, label) {
  if (!file.exists(file)) stop("Missing ", label, ": ", file, call. = FALSE)
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

read_matrix_with_gene_ids <- function(file, label) {
  df <- safe_read_csv(file)
  gene_col <- find_col(df, c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id"), TRUE, paste0(label, " gene ID"))

  gene_ids <- clean_string(df[[gene_col]])
  sample_cols <- setdiff(colnames(df), gene_col)
  numeric_cols <- sample_cols[vapply(df[sample_cols], function(x) {
    suppressWarnings(all(is.na(x) | !is.na(as.numeric(as.character(x)))))
  }, logical(1))]

  if (length(numeric_cols) == 0) stop("No numeric sample columns found in ", label, ": ", file, call. = FALSE)

  mat <- as.matrix(data.frame(
    lapply(df[, numeric_cols, drop = FALSE], function(x) as.numeric(as.character(x))),
    check.names = FALSE
  ))
  colnames(mat) <- numeric_cols
  rownames(mat) <- gene_ids

  keep <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep, , drop = FALSE]

  if (any(duplicated(rownames(mat)))) {
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }

  list(df = df, gene_col = gene_col, mat = mat)
}

find_length_matrix_file <- function() {
  candidates <- c(
    file.path(chronic_tables_dir, "step25v3_pig_chronic_main_comparison_gene_level_length_matrix.csv"),
    file.path(chronic_tables_dir, "step25v3_pig_chronic_main_comparison_length_matrix.csv"),
    file.path(chronic_tables_dir, "step25_pig_chronic_main_comparison_gene_level_length_matrix.csv"),
    file.path(chronic_tables_dir, "step25_pig_chronic_main_comparison_length_matrix.csv"),
    file.path(chronic_tables_dir, "pig_chronic_main_comparison_gene_level_length_matrix.csv"),
    file.path(chronic_tables_dir, "pig_chronic_main_comparison_length_matrix.csv"),
    file.path(chronic_tables_dir, "step25v3_pig_chronic_gene_level_length_matrix.csv"),
    file.path(chronic_tables_dir, "step25_pig_chronic_gene_level_length_matrix.csv"),
    file.path(chronic_tables_dir, "pig_chronic_gene_level_length_matrix.csv")
  )

  hit <- find_existing_file(candidates)
  if (!is.na(hit)) return(hit)

  all_csv <- list.files(chronic_tables_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(all_csv) == 0) return(NA_character_)

  bn <- tolower(basename(all_csv))
  score <- rep(0, length(all_csv))
  score[grepl("length", bn)] <- score[grepl("length", bn)] + 50
  score[grepl("matrix", bn)] <- score[grepl("matrix", bn)] + 20
  score[grepl("gene", bn)] <- score[grepl("gene", bn)] + 10
  score[grepl("main|comparison|52", bn)] <- score[grepl("main|comparison|52", bn)] + 10
  score[grepl("count|abundance|tpm|summary|manifest|de|gsea|score|audit", bn)] <- score[grepl("count|abundance|tpm|summary|manifest|de|gsea|score|audit", bn)] - 100

  cand <- data.frame(file = normalizePath(all_csv, winslash = "/", mustWork = FALSE), score = score, stringsAsFactors = FALSE)
  cand <- cand[order(-cand$score, cand$file), , drop = FALSE]
  write_csv(cand, file.path(out_dir, "chronic_step4_length_matrix_file_candidates.csv"))

  if (nrow(cand) > 0 && cand$score[1] >= 60) return(cand$file[1])
  NA_character_
}

build_symbol_map <- function(count_df, count_gene_col) {
  # First try gene symbol column inside count matrix.
  symbol_col <- find_col(count_df, c("gene_name", "gene_symbol", "SYMBOL", "symbol", "pig_symbol", "external_gene_name"), FALSE, "count matrix gene symbol")
  if (!is.na(symbol_col)) {
    tmp <- data.frame(
      gene_id = clean_string(count_df[[count_gene_col]]),
      gene_symbol = clean_string(count_df[[symbol_col]]),
      source = paste0("count_matrix_column:", symbol_col),
      stringsAsFactors = FALSE
    )
    tmp <- tmp[!is.na(tmp$gene_id) & tmp$gene_id != "" & !is.na(tmp$gene_symbol) & tmp$gene_symbol != "", , drop = FALSE]
    tmp <- tmp[!duplicated(tmp$gene_id), , drop = FALSE]
    if (nrow(tmp) > 1000) return(tmp)
  }

  # Then search likely tx2gene / annotation files.
  all_csv <- list.files(chronic_tables_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  candidates <- all_csv[grepl("tx2gene|annotation|annot|gtf|gene.*name|gene.*symbol", tolower(basename(all_csv)))]
  if (length(candidates) > 0) {
    for (f in candidates) {
      df <- tryCatch(safe_read_csv(f), error = function(e) NULL)
      if (is.null(df)) next
      gid <- find_col(df, c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id"), FALSE, "annotation gene ID")
      sym <- find_col(df, c("gene_name", "gene_symbol", "SYMBOL", "symbol", "pig_symbol", "external_gene_name"), FALSE, "annotation gene symbol")
      if (!is.na(gid) && !is.na(sym)) {
        tmp <- data.frame(
          gene_id = clean_string(df[[gid]]),
          gene_symbol = clean_string(df[[sym]]),
          source = paste0("annotation_file:", f),
          stringsAsFactors = FALSE
        )
        tmp <- tmp[!is.na(tmp$gene_id) & tmp$gene_id != "" & !is.na(tmp$gene_symbol) & tmp$gene_symbol != "", , drop = FALSE]
        tmp <- tmp[!duplicated(tmp$gene_id), , drop = FALSE]
        if (nrow(tmp) > 1000) return(tmp)
      }
    }
  }

  # Last resort: use an old DE output only for gene_id -> gene_name annotation, not for statistics.
  old_de_file <- file.path(chronic_tables_dir, "step26_pig_chronic_DE_ACLT_alone_52W_vs_Control_52W_QLF_tximport_offset.csv")
  if (file.exists(old_de_file)) {
    df <- safe_read_csv(old_de_file)
    gid <- find_col(df, c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id"), FALSE, "old DE gene ID")
    sym <- find_col(df, c("gene_name", "gene_symbol", "SYMBOL", "symbol", "pig_symbol", "external_gene_name"), FALSE, "old DE gene symbol")
    if (!is.na(gid) && !is.na(sym)) {
      tmp <- data.frame(
        gene_id = clean_string(df[[gid]]),
        gene_symbol = clean_string(df[[sym]]),
        source = paste0("old_DE_annotation_only:", old_de_file),
        stringsAsFactors = FALSE
      )
      tmp <- tmp[!is.na(tmp$gene_id) & tmp$gene_id != "" & !is.na(tmp$gene_symbol) & tmp$gene_symbol != "", , drop = FALSE]
      tmp <- tmp[!duplicated(tmp$gene_id), , drop = FALSE]
      if (nrow(tmp) > 1000) return(tmp)
    }
  }

  data.frame(gene_id = character(0), gene_symbol = character(0), source = character(0), stringsAsFactors = FALSE)
}

flatten_fgsea <- function(res, contrast_label) {
  out <- as.data.frame(res, stringsAsFactors = FALSE)
  if (nrow(out) == 0) return(out)
  if ("leadingEdge" %in% colnames(out)) {
    out$leadingEdge_size <- vapply(out$leadingEdge, length, integer(1))
    out$leadingEdge_genes <- vapply(out$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    out$leadingEdge <- NULL
  }
  out$contrast <- contrast_label
  out <- out[, c("contrast", setdiff(colnames(out), "contrast")), drop = FALSE]
  out
}

enrichment_curve_df <- function(stats, pathway_genes, pathway_name, contrast_label, nes, padj, pval) {
  stats <- sort(stats, decreasing = TRUE)
  hits <- names(stats) %in% pathway_genes
  N <- length(stats)
  Nh <- sum(hits)
  Nm <- N - Nh
  if (Nh == 0 || Nm == 0) return(data.frame())

  weights <- abs(stats)
  Phit <- ifelse(hits, weights / sum(weights[hits]), 0)
  Pmiss <- ifelse(!hits, 1 / Nm, 0)
  running <- cumsum(Phit - Pmiss)

  data.frame(
    contrast = contrast_label,
    pathway = pathway_name,
    rank = seq_along(stats),
    running_ES = running,
    is_hit = hits,
    gene_symbol = names(stats),
    NES = nes,
    padj = padj,
    pval = pval,
    stringsAsFactors = FALSE
  )
}

pretty_pathway <- function(x) {
  x <- as.character(x)
  x[x == "HALLMARK_INFLAMMATORY_RESPONSE"] <- "Inflammatory response"
  x[x == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT / ECM remodeling"
  x
}

format_fdr <- function(x) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "e", digits = 2))
}


axis_label_from_mds <- function(mds_obj, dim_index) {
  pct <- NA_real_

  if (!is.null(mds_obj$var.explained)) {
    ve <- suppressWarnings(as.numeric(mds_obj$var.explained))
    if (length(ve) >= dim_index) {
      pct <- ve[dim_index]
      if (is.finite(pct) && pct <= 1) pct <- pct * 100
    }
  }

  if (is.finite(pct)) {
    return(paste0("Leading logFC dim ", dim_index, " (", round(pct), "%)"))
  }

  paste0("Leading logFC dim ", dim_index)
}

## =========================
## 3. Input files
## =========================

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

length_file <- find_length_matrix_file()

if (is.na(manifest_file)) stop("Could not find chronic main comparison manifest under: ", chronic_tables_dir, call. = FALSE)
if (is.na(count_file)) stop("Could not find chronic main comparison counts matrix under: ", chronic_tables_dir, call. = FALSE)
if (is.na(length_file)) {
  stop(
    "Could not find a tximport effective-length matrix. This step requires the length matrix for the locked chronic DE method with effective-length offset. ",
    "Candidate list, if generated, is saved in output directory.",
    call. = FALSE
  )
}

cat("Manifest file: ", manifest_file, "\n", sep = "")
cat("Count matrix file: ", count_file, "\n", sep = "")
cat("Effective length matrix file: ", length_file, "\n\n", sep = "")

## =========================
## 4. Load manifest, counts, and lengths
## =========================

manifest <- safe_read_csv(manifest_file)
count_obj <- read_matrix_with_gene_ids(count_file, "chronic count matrix")
length_obj <- read_matrix_with_gene_ids(length_file, "chronic effective length matrix")

count_df_raw <- count_obj$df
count_gene_col <- count_obj$gene_col
count_mat_all <- count_obj$mat
length_mat_all <- length_obj$mat

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

missing_count_samples <- setdiff(sample_info$sample_id, colnames(count_mat_all))
missing_length_samples <- setdiff(sample_info$sample_id, colnames(length_mat_all))
if (length(missing_count_samples) > 0) {
  write_csv(data.frame(missing_sample_id = missing_count_samples), file.path(out_dir, "chronic_step4_missing_samples_in_count_matrix.csv"))
  stop("Some main-comparison samples are missing from count matrix.", call. = FALSE)
}
if (length(missing_length_samples) > 0) {
  write_csv(data.frame(missing_sample_id = missing_length_samples), file.path(out_dir, "chronic_step4_missing_samples_in_length_matrix.csv"))
  stop("Some main-comparison samples are missing from length matrix.", call. = FALSE)
}

common_gene_ids <- intersect(rownames(count_mat_all), rownames(length_mat_all))
if (length(common_gene_ids) == 0) stop("No overlapping gene IDs between count and length matrices.", call. = FALSE)

count_mat <- count_mat_all[common_gene_ids, sample_info$sample_id, drop = FALSE]
length_mat <- length_mat_all[common_gene_ids, sample_info$sample_id, drop = FALSE]

# Preserve count-file gene order
count_mat <- count_mat[rownames(count_mat_all)[rownames(count_mat_all) %in% common_gene_ids], , drop = FALSE]
length_mat <- length_mat[rownames(count_mat), , drop = FALSE]

cat("Main comparison samples:\n")
print(table(sample_info$group))
cat("\n")
cat("Count/length common genes before filtering: ", nrow(count_mat), "\n\n", sep = "")

## =========================
## 5. Gene symbol annotation
## =========================

symbol_map <- build_symbol_map(count_df_raw, count_gene_col)
if (nrow(symbol_map) == 0) {
  stop("Could not obtain gene_id -> gene_symbol mapping for GSEA. Please provide a tx2gene/annotation table or count matrix with gene_name/gene_symbol.", call. = FALSE)
}

symbol_map_file <- file.path(out_dir, "chronic_step4_gene_id_to_symbol_map_used_for_GSEA.csv")
write_csv(symbol_map, symbol_map_file)

cat("Gene symbol map rows: ", nrow(symbol_map), "\n", sep = "")
cat("Gene symbol map source preview:\n")
print(head(unique(symbol_map$source), 5))
cat("\n")

## =========================
## 6. edgeR QL DE with effective-length offset
## =========================

group <- sample_info$group
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("^group", "", colnames(design))
contrast_vec <- makeContrasts(ACLT_alone_52W - Control_52W, levels = design)

dge0 <- DGEList(counts = count_mat, group = group)
keep_expr <- filterByExpr(dge0, design = design)
dge <- dge0[keep_expr, , keep.lib.sizes = FALSE]
length_f <- length_mat[rownames(dge), , drop = FALSE]

length_ok <- apply(length_f, 1, function(x) all(is.finite(x) & x > 0))
if (sum(!length_ok) > 0) {
  write_csv(
    data.frame(gene_id = rownames(length_f)[!length_ok]),
    file.path(out_dir, "chronic_step4_genes_removed_due_to_invalid_effective_length.csv")
  )
}
dge <- dge[length_ok, , keep.lib.sizes = FALSE]
length_f <- length_f[length_ok, , drop = FALSE]

normMat <- length_f / exp(rowMeans(log(length_f)))
normCounts <- dge$counts / normMat
eff_lib <- calcNormFactors(normCounts, method = "TMM") * colSums(normCounts)

dge$offset <- log(normMat) + matrix(log(eff_lib), nrow = nrow(normMat), ncol = ncol(normMat), byrow = TRUE)
dge$samples$norm.factors <- 1

dge <- estimateDisp(dge, design = design, robust = TRUE)
fit <- glmQLFit(dge, design = design, robust = TRUE)
qlf <- glmQLFTest(fit, contrast = contrast_vec)

de_table <- topTags(qlf, n = Inf, sort.by = "none")$table
de_table$gene_id <- rownames(de_table)

de_table <- merge(
  de_table,
  symbol_map[, c("gene_id", "gene_symbol")],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
de_table <- de_table[, c("gene_id", "gene_symbol", setdiff(colnames(de_table), c("gene_id", "gene_symbol"))), drop = FALSE]

# Sort for user readability, while retaining all rows.
de_table_sorted <- de_table[order(de_table$PValue, de_table$gene_id), , drop = FALSE]

strict_de <- de_table[!is.na(de_table$FDR) & de_table$FDR < 0.05 & abs(de_table$logFC) > 1, , drop = FALSE]
strict_up <- strict_de[strict_de$logFC > 1, , drop = FALSE]
strict_down <- strict_de[strict_de$logFC < -1, , drop = FALSE]

de_file <- file.path(out_dir, "chronic_step4_DE_ACLT_alone_52W_vs_Control_52W_edgeR_QLF_tximport_offset.csv")
de_sorted_file <- file.path(out_dir, "chronic_step4_DE_ACLT_alone_52W_vs_Control_52W_edgeR_QLF_tximport_offset_sorted_by_PValue.csv")
strict_de_file <- file.path(out_dir, "chronic_step4_strict_DEGs_FDR0.05_abslogFC1.csv")
write_csv(de_table, de_file)
write_csv(de_table_sorted, de_sorted_file)
write_csv(strict_de, strict_de_file)

# Save exploratory QC plots based on filtered genes and offset-normalized logCPM.
logcpm_qc <- cpm(dge, log = TRUE, prior.count = 1)
pca <- prcomp(t(logcpm_qc), center = TRUE, scale. = FALSE)
pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = group,
  stringsAsFactors = FALSE
)
pc_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = group, color = group)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.7, alpha = 0.95) +
  labs(
    title = "Chronic pig PCA before DE",
    x = sprintf("PC1 (%.1f%%)", pc_var[1]),
    y = sprintf("PC2 (%.1f%%)", pc_var[2])
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pca_png <- file.path(qc_fig_dir, "chronic_step4_PCA_Control52W_vs_ACLTalone52W.png")
pca_pdf <- file.path(qc_fig_dir, "chronic_step4_PCA_Control52W_vs_ACLTalone52W.pdf")
ggsave(pca_png, p_pca, width = 5.8, height = 4.8, dpi = 320)
ggsave(pca_pdf, p_pca, width = 5.8, height = 4.8)

# Step4 QC MDS: use the same DGEList object, but replace overlapping text labels
# with colored points for readability. This changes only the visualization layer;
# DE, offsets, dispersion estimates and GSEA are unchanged.
mds_obj <- limma::plotMDS(
  dge,
  labels = colnames(dge$counts),
  col = as.integer(group),
  main = "Chronic pig MDS before DE",
  plot = FALSE
)

mds_df <- data.frame(
  sample_id = colnames(dge$counts),
  group = as.character(group),
  MDS1 = as.numeric(mds_obj$x),
  MDS2 = as.numeric(mds_obj$y),
  stringsAsFactors = FALSE
)

mds_coord_file <- file.path(out_dir, "chronic_step4_MDS_Control52W_vs_ACLTalone52W_colored_points_coordinates.csv")
write_csv(mds_df, mds_coord_file)

p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", color = "grey75") +
  stat_ellipse(linewidth = 0.5, linetype = "dashed", alpha = 0.6, show.legend = FALSE) +
  geom_point(size = 3.8, alpha = 0.95) +
  scale_color_manual(
    values = c("Control_52W" = "black", "ACLT_alone_52W" = "#E64B5D"),
    breaks = c("Control_52W", "ACLT_alone_52W")
  ) +
  labs(
    title = "Chronic pig MDS before DE",
    subtitle = "Control_52W vs ACLT_alone_52W; colored points replace overlapping text labels",
    x = axis_label_from_mds(mds_obj, 1),
    y = axis_label_from_mds(mds_obj, 2),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.text = element_text(size = 10)
  )

mds_png <- file.path(qc_fig_dir, "chronic_step4_MDS_Control52W_vs_ACLTalone52W.png")
mds_pdf <- file.path(qc_fig_dir, "chronic_step4_MDS_Control52W_vs_ACLTalone52W.pdf")
ggsave(mds_png, p_mds, width = 6.2, height = 5.2, dpi = 300)
ggsave(mds_pdf, p_mds, width = 6.2, height = 5.2, device = cairo_pdf)

## =========================
## 7. Ranked list for GSEA
## =========================

rank_df <- de_table
rank_df$gene_symbol_clean <- clean_string(rank_df$gene_symbol)
rank_df$PValue_rank <- suppressWarnings(as.numeric(rank_df$PValue))
rank_df$PValue_rank[!is.finite(rank_df$PValue_rank) | rank_df$PValue_rank <= 0] <- .Machine$double.xmin
rank_df$rank_statistic <- sign(rank_df$logFC) * (-log10(rank_df$PValue_rank))
rank_df$abs_rank_statistic <- abs(rank_df$rank_statistic)

rank_keep <- !is.na(rank_df$gene_symbol_clean) & rank_df$gene_symbol_clean != "" &
  !is.na(rank_df$logFC) &
  !is.na(rank_df$PValue_rank) &
  is.finite(rank_df$rank_statistic)

rank_df <- rank_df[rank_keep, , drop = FALSE]
rank_df <- rank_df[order(-rank_df$abs_rank_statistic, rank_df$gene_symbol_clean), , drop = FALSE]
duplicated_symbol_rows_removed <- sum(duplicated(rank_df$gene_symbol_clean))
rank_df_unique <- rank_df[!duplicated(rank_df$gene_symbol_clean), , drop = FALSE]
rank_df_unique <- rank_df_unique[order(rank_df_unique$rank_statistic, decreasing = TRUE), , drop = FALSE]

rank_vector <- rank_df_unique$rank_statistic
names(rank_vector) <- rank_df_unique$gene_symbol_clean
rank_vector <- sort(rank_vector, decreasing = TRUE)

ranked_file <- file.path(out_dir, "chronic_step4_ranked_genes_ACLT_alone_52W_vs_Control_52W.csv")
write_csv(
  rank_df_unique[, c("gene_id", "gene_symbol", "logFC", "PValue", "FDR", "rank_statistic"), drop = FALSE],
  ranked_file
)

## =========================
## 8. Full Sus scrofa Hallmark GSEA
## =========================

msig_df <- tryCatch(
  msigdbr::msigdbr(species = "Sus scrofa", collection = "H"),
  error = function(e) msigdbr::msigdbr(species = "Sus scrofa", category = "H")
)

pathway_col <- find_col(msig_df, c("gs_name", "gs_exact_source"), TRUE, "MSigDB pathway column")
gene_col <- find_col(msig_df, c("gene_symbol", "db_gene_symbol"), TRUE, "MSigDB gene symbol column")

hallmark_members <- data.frame(
  pathway = as.character(msig_df[[pathway_col]]),
  gene_symbol = as.character(msig_df[[gene_col]]),
  stringsAsFactors = FALSE
)
hallmark_members <- hallmark_members[!is.na(hallmark_members$pathway) & !is.na(hallmark_members$gene_symbol) &
                                       hallmark_members$pathway != "" & hallmark_members$gene_symbol != "", , drop = FALSE]
hallmark_members <- unique(hallmark_members)

pathways <- split(hallmark_members$gene_symbol, hallmark_members$pathway)
pathways <- lapply(pathways, unique)
pathways <- pathways[order(names(pathways))]

target_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

missing_targets <- setdiff(target_pathways, names(pathways))
if (length(missing_targets) > 0) {
  stop("Target Hallmark pathways missing from Sus scrofa Hallmark collection: ", paste(missing_targets, collapse = ", "), call. = FALSE)
}

fgsea_minSize <- 10
fgsea_maxSize <- 500
fgsea_eps <- 0
fgsea_nproc <- 1
fgsea_seed <- 1

coverage_df <- data.frame(
  pathway = names(pathways),
  n_genes_in_pathway = vapply(pathways, length, integer(1)),
  n_genes_present_in_rank = vapply(pathways, function(gs) sum(gs %in% names(rank_vector)), integer(1)),
  stringsAsFactors = FALSE
)

n_pathways_passing <- sum(coverage_df$n_genes_present_in_rank >= fgsea_minSize &
                            coverage_df$n_genes_present_in_rank <= fgsea_maxSize)
if (n_pathways_passing == 0) {
  write_csv(coverage_df, file.path(out_dir, "chronic_step4_hallmark_coverage_failed_zero_pathways.csv"))
  stop("No Hallmark pathways passed fgsea size filtering. Check gene symbols.", call. = FALSE)
}

set.seed(fgsea_seed)
fg <- fgsea::fgsea(
  pathways = pathways,
  stats = rank_vector,
  minSize = fgsea_minSize,
  maxSize = fgsea_maxSize,
  eps = fgsea_eps,
  nproc = fgsea_nproc
)

fg_out <- flatten_fgsea(fg, "ACLT_alone_52W_vs_Control_52W")
fg_out <- fg_out[order(fg_out$padj, fg_out$pathway), , drop = FALSE]

target_results <- fg_out[fg_out$pathway %in% target_pathways, , drop = FALSE]
target_results <- target_results[match(target_pathways, target_results$pathway), , drop = FALSE]

hallmark_members_file <- file.path(out_dir, "chronic_step4_sus_scrofa_hallmark_membership.csv")
coverage_file <- file.path(out_dir, "chronic_step4_hallmark_gene_set_coverage.csv")
fgsea_all_file <- file.path(out_dir, "chronic_step4_fgsea_full_hallmark_ACLT_alone_52W_vs_Control_52W.csv")
target_results_file <- file.path(out_dir, "chronic_step4_target_pathway_results_from_full_hallmark.csv")
write_csv(hallmark_members, hallmark_members_file)
write_csv(coverage_df, coverage_file)
write_csv(fg_out, fgsea_all_file)
write_csv(target_results, target_results_file)

## =========================
## 9. Target running curve plot
## =========================

running_list <- list()
for (pw in target_pathways) {
  rr <- target_results[target_results$pathway == pw, , drop = FALSE]
  if (nrow(rr) != 1) stop("Target pathway result not found exactly once: ", pw, call. = FALSE)
  running_list[[pw]] <- enrichment_curve_df(
    stats = rank_vector,
    pathway_genes = pathways[[pw]],
    pathway_name = pw,
    contrast_label = "ACLT_alone_52W_vs_Control_52W",
    nes = rr$NES[1],
    padj = rr$padj[1],
    pval = rr$pval[1]
  )
}

running_df <- bind_rows(running_list)
running_df$pathway_pretty <- pretty_pathway(running_df$pathway)

label_df <- target_results
label_df$pathway_pretty <- pretty_pathway(label_df$pathway)
label_df$panel_label <- paste0(
  label_df$pathway_pretty,
  "\nNES = ",
  sprintf("%.3f", as.numeric(label_df$NES)),
  "   |   FDR = ",
  format_fdr(label_df$padj)
)

label_key <- label_df[, c("pathway", "panel_label"), drop = FALSE]
running_df <- merge(running_df, label_key, by = "pathway", all.x = TRUE)
running_df$panel_label <- factor(running_df$panel_label, levels = label_df$panel_label[match(target_pathways, label_df$pathway)])

peak_df <- aggregate(running_ES ~ panel_label, data = running_df, FUN = max)
colnames(peak_df)[2] <- "peak_ES"

rug_df <- running_df[running_df$is_hit, , drop = FALSE]
rug_df$rug_ymin <- -0.045
rug_df$rug_ymax <- 0.045

p_gsea <- ggplot(running_df, aes(x = rank, y = running_ES)) +
  geom_hline(
    data = peak_df,
    aes(yintercept = peak_ES),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_segment(
    data = rug_df,
    aes(x = rank, xend = rank, y = rug_ymin, yend = rug_ymax),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.35,
    linewidth = 0.18
  ) +
  geom_line(color = "#00CD00", linewidth = 0.75) +
  facet_wrap(~ panel_label, ncol = 2, scales = "free_y") +
  labs(x = NULL, y = "Running enrichment score") +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13.5, face = "plain", lineheight = 1.15),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = element_line(color = "#F0F0F0", linewidth = 0.25),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
  )

gsea_png <- file.path(chronic_figures_dir, "Figure5B_current78_chronic_Hallmark_GSEA_running_curves_mouse_style.png")
gsea_pdf <- file.path(chronic_figures_dir, "Figure5B_current78_chronic_Hallmark_GSEA_running_curves_mouse_style.pdf")
ggsave(gsea_png, p_gsea, width = 11.8, height = 4.2, dpi = 320)
ggsave(gsea_pdf, p_gsea, width = 11.8, height = 4.2)

running_file <- file.path(out_dir, "chronic_step4_target_pathway_running_curve_source_data.csv")
write_csv(running_df, running_file)

## =========================
## 10. Save method records and summaries
## =========================

de_summary <- data.frame(
  metric = c(
    "comparison",
    "input_genes_common_count_length",
    "genes_after_filterByExpr",
    "genes_removed_due_to_invalid_effective_length",
    "genes_tested_after_length_filter",
    "strict_DEG_FDR0.05_abslogFC1_total",
    "strict_DEG_up",
    "strict_DEG_down",
    "design",
    "contrast",
    "edgeR_method",
    "effective_length_offset_used"
  ),
  value = c(
    "ACLT_alone_52W_vs_Control_52W",
    nrow(count_mat),
    sum(keep_expr),
    sum(!length_ok),
    nrow(de_table),
    nrow(strict_de),
    nrow(strict_up),
    nrow(strict_down),
    paste(colnames(design), collapse = "; "),
    "ACLT_alone_52W - Control_52W",
    "DGEList + filterByExpr + tximport effective-length offset + estimateDisp + glmQLFit(robust=TRUE) + glmQLFTest",
    TRUE
  ),
  stringsAsFactors = FALSE
)
de_summary_file <- file.path(out_dir, "chronic_step4_DE_run_summary.csv")
write_csv(de_summary, de_summary_file)

gsea_summary <- data.frame(
  metric = c(
    "ranked_genes_after_symbol_dedup",
    "duplicated_symbol_rows_removed",
    "n_full_hallmark_gene_sets_loaded",
    "n_full_hallmark_gene_sets_tested",
    "fgsea_minSize",
    "fgsea_maxSize",
    "fgsea_eps",
    "fgsea_nproc",
    "fgsea_seed",
    "rank_statistic",
    "gene_identifier_for_GSEA",
    "full_collection_FDR_used",
    "target_pathways",
    "inflammatory_response_NES",
    "inflammatory_response_FDR",
    "EMT_NES",
    "EMT_FDR",
    "gsea_png",
    "gsea_pdf"
  ),
  value = c(
    length(rank_vector),
    duplicated_symbol_rows_removed,
    length(pathways),
    nrow(fg_out),
    fgsea_minSize,
    fgsea_maxSize,
    fgsea_eps,
    fgsea_nproc,
    fgsea_seed,
    "sign(logFC) * -log10(PValue)",
    "pig gene symbols",
    TRUE,
    paste(target_pathways, collapse = "; "),
    target_results$NES[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"],
    target_results$padj[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"],
    target_results$NES[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
    target_results$padj[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
    gsea_png,
    gsea_pdf
  ),
  stringsAsFactors = FALSE
)
gsea_summary_file <- file.path(out_dir, "chronic_step4_GSEA_run_summary.csv")
write_csv(gsea_summary, gsea_summary_file)

version_df <- data.frame(
  item = c(
    "method_version",
    "R_version",
    "edgeR_version",
    "limma_version",
    "fgsea_version",
    "msigdbr_version",
    "ggplot2_version",
    "DE_input_counts",
    "DE_input_lengths",
    "DE_input_manifest",
    "GSEA_collection",
    "GSEA_parameters",
    "QC_MDS_visualization",
    "analysis_role"
  ),
  value = c(
    method_version,
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("limma")),
    as.character(utils::packageVersion("fgsea")),
    as.character(utils::packageVersion("msigdbr")),
    as.character(utils::packageVersion("ggplot2")),
    count_file,
    length_file,
    manifest_file,
    "full Sus scrofa MSigDB Hallmark collection H",
    paste0("minSize=", fgsea_minSize, "; maxSize=", fgsea_maxSize, "; eps=", fgsea_eps, "; nproc=", fgsea_nproc, "; seed=", fgsea_seed),
    "limma::plotMDS coordinates from the same Step4 DGEList, displayed as colored points without text labels",
    "chronic extension DE/GSEA; no signature/core redefinition"
  ),
  stringsAsFactors = FALSE
)
version_file <- file.path(out_dir, "chronic_step4_versions_and_method_records.csv")
write_csv(version_df, version_file)

fgsea_reproducibility_settings <- data.frame(
  contrast = "ACLT_alone_52W_vs_Control_52W",
  seed = fgsea_seed,
  nproc = fgsea_nproc,
  minSize = fgsea_minSize,
  maxSize = fgsea_maxSize,
  eps = fgsea_eps,
  fgsea_implementation = "fgsea default multilevel implementation; no fixed nperm specified",
  full_collection_FDR_used = TRUE,
  stringsAsFactors = FALSE
)
fgsea_reproducibility_file <- file.path(out_dir, "chronic_step4_fgsea_reproducibility_settings.csv")
write_csv(fgsea_reproducibility_settings, fgsea_reproducibility_file)

input_audit <- rbind(
  file_audit(manifest_file, "chronic main comparison manifest"),
  file_audit(count_file, "chronic main comparison gene-level estimated count matrix"),
  file_audit(length_file, "chronic main comparison effective-length matrix")
)
input_audit_file <- file.path(out_dir, "chronic_step4_input_file_audit.csv")
write_csv(input_audit, input_audit_file)

sample_file <- file.path(out_dir, "chronic_step4_sample_metadata_used.csv")
write_csv(sample_info, sample_file)

plot_files <- data.frame(
  plot = c("PCA", "MDS_colored_points", "Figure5B_current78_chronic_Hallmark_GSEA_running_curves_mouse_style"),
  png = c(pca_png, mds_png, gsea_png),
  pdf = c(pca_pdf, mds_pdf, gsea_pdf),
  stringsAsFactors = FALSE
)
plot_files_file <- file.path(out_dir, "chronic_step4_plot_files.csv")
write_csv(plot_files, plot_files_file)

overall_summary <- data.frame(
  metric = c(
    "manifest_file",
    "count_file",
    "length_file",
    "n_Control_52W",
    "n_ACLT_alone_52W",
    "genes_tested_DE",
    "strict_DEG_total",
    "strict_DEG_up",
    "strict_DEG_down",
    "ranked_genes_GSEA",
    "hallmark_pathways_tested",
    "fgsea_seed",
    "fgsea_nproc",
    "inflammatory_response_NES",
    "inflammatory_response_FDR",
    "EMT_NES",
    "EMT_FDR",
    "DE_output_file",
    "GSEA_full_output_file",
    "target_pathway_output_file",
    "Figure5B_png",
    "Figure5B_pdf",
    "Step4_QC_MDS_png",
    "Step4_QC_MDS_pdf",
    "Step4_QC_MDS_coordinate_file",
    "output_dir"
  ),
  value = c(
    manifest_file,
    count_file,
    length_file,
    sum(sample_info$group == "Control_52W"),
    sum(sample_info$group == "ACLT_alone_52W"),
    nrow(de_table),
    nrow(strict_de),
    nrow(strict_up),
    nrow(strict_down),
    length(rank_vector),
    nrow(fg_out),
    fgsea_seed,
    fgsea_nproc,
    target_results$NES[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"],
    target_results$padj[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"],
    target_results$NES[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
    target_results$padj[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
    de_file,
    fgsea_all_file,
    target_results_file,
    gsea_png,
    gsea_pdf,
    mds_png,
    mds_pdf,
    mds_coord_file,
    out_dir
  ),
  stringsAsFactors = FALSE
)
overall_summary_file <- file.path(out_dir, "chronic_step4_current78_DE_GSEA_overall_summary.csv")
write_csv(overall_summary, overall_summary_file)

saveRDS(
  list(
    sample_info = sample_info,
    design = design,
    dge = dge,
    fit = fit,
    qlf = qlf,
    de_table = de_table,
    de_summary = de_summary,
    pca_df = pca_df,
    mds_df = mds_df,
    rank_df_unique = rank_df_unique,
    rank_vector = rank_vector,
    pathways = pathways,
    hallmark_members = hallmark_members,
    coverage_df = coverage_df,
    fg_out = fg_out,
    target_results = target_results,
    running_df = running_df,
    gsea_summary = gsea_summary,
    overall_summary = overall_summary,
    version_df = version_df
  ),
  file.path(obj_dir, "chronic_step4_current78_DE_GSEA_workspace.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "chronic_step4_sessionInfo.txt"))

## =========================
## 11. Console and summary-to-send log
## =========================

cat("\n===== Chronic Step4 overall summary =====\n")
print(overall_summary, row.names = FALSE)

cat("\nDE summary:\n")
print(de_summary, row.names = FALSE)

cat("\nTarget pathway GSEA results:\n")
print(target_results[, c("pathway", "NES", "pval", "padj", "size", "leadingEdge_size"), drop = FALSE], row.names = FALSE)

cat("\nGSEA summary:\n")
print(gsea_summary, row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== Chronic Step4 current78 DE + full Hallmark GSEA SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Overall summary:", summary_con)
writeLines(capture.output(print(overall_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("DE summary:", summary_con)
writeLines(capture.output(print(de_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Target pathway GSEA results:", summary_con)
writeLines(capture.output(print(target_results[, c("pathway", "NES", "pval", "padj", "size", "leadingEdge_size"), drop = FALSE], row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("GSEA summary:", summary_con)
writeLines(capture.output(print(gsea_summary, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("fgsea reproducibility settings:", summary_con)
writeLines(capture.output(print(fgsea_reproducibility_settings, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Step4 QC MDS visualization:", summary_con)
writeLines("MDS coordinates were generated with limma::plotMDS from the same Step4 DGEList object; the plotted figure uses colored points without overlapping group text labels.", summary_con)
writeLines(paste0("MDS coordinate file: ", mds_coord_file), summary_con)
writeLines(paste0("MDS PNG: ", mds_png), summary_con)
writeLines(paste0("MDS PDF: ", mds_pdf), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("DE/GSEA is genome-wide and independent of signature/core definition.", summary_con)
writeLines(paste0("GSEA is run on the full Sus scrofa Hallmark collection with minSize=", fgsea_minSize, ", maxSize=", fgsea_maxSize, ", eps=", fgsea_eps, ", nproc=", fgsea_nproc, ", seed=", fgsea_seed, "."), summary_con)
writeLines("Target pathway FDR values are retained from the full Hallmark-level analysis.", summary_con)
close(summary_con)

sink()

cat("\nChronic Step4 current78 DE + full Hallmark GSEA completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
