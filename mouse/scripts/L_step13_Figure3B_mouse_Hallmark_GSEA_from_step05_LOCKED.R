# Step13: Figure 3B mouse Hallmark GSEA running enrichment curves from locked Step05 DE tables.
# English note:
# This locked revision removes the previous fuzzy DE-file search logic.
# It reads only the exact Step05 paired limma-voom duplicateCorrelation DE output files:
#   - step05_DE_1W_ACLR_vs_Contra_limma_voom.csv
#   - step05_DE_4W_ACLR_vs_Contra_limma_voom.csv
#
# Locked method:
# 1. Do NOT use persistent genes or CellAge-overlapping genes as the GSEA input.
# 2. Use the full gene-level ranked list from each mouse paired limma-voom DE result.
# 3. Rank statistic = sign(logFC) * -log10(P value used for ranking), where logFC > 0 means ACLR > Contra.
#    Non-positive or non-finite P values are replaced by .Machine$double.xmin before ranking.
# 4. Use Mus musculus MSigDB Hallmark gene sets from msigdbr.
# 5. Run fgsea on the full Hallmark collection with minSize = 10, maxSize = 500, eps = 0.
# 6. Adjusted P values are retained from the full Hallmark collection analysis.
# 7. Figure 3B focuses only on the two pre-specified targeted pathway readouts:
#    - HALLMARK_INFLAMMATORY_RESPONSE
#    - HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
# 8. HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION is interpreted biologically as
#    an EMT / ECM remodeling-related transcriptional module, but the actual gene set is the
#    standard Hallmark EMT set, not a custom ECM gene set.
#
# This script does not rerun differential expression analysis.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths and locked inputs
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"
de_dir <- file.path(base_dir, "03_DE_analysis")

de_1w_file <- file.path(de_dir, "step05_DE_1W_ACLR_vs_Contra_limma_voom.csv")
de_4w_file <- file.path(de_dir, "step05_DE_4W_ACLR_vs_Contra_limma_voom.csv")

table_dir <- file.path(base_dir, "07_tables", "step13_Figure3B_mouse_Hallmark_GSEA")
figure_dir <- file.path(base_dir, "06_figures", "Figure3")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")
tmp_dir <- file.path(base_dir, "00_temp", "step13_fgsea_tmp")

for (d in c(table_dir, figure_dir, log_dir, script_dir, tmp_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

Sys.setenv(TMPDIR = tmp_dir, TEMP = tmp_dir, TMP = tmp_dir)

step_name <- "step13_Figure3B_mouse_Hallmark_GSEA"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))
archive_file <- file.path(script_dir, "step13_Figure3B_mouse_Hallmark_GSEA.R")

target_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

fgsea_min_size <- 10
fgsea_max_size <- 500
fgsea_eps <- 0
fgsea_nproc <- 1
fgsea_seed_1w <- 12345
fgsea_seed_4w <- 12346

png_file <- file.path(figure_dir, "Figure3B_mouse_Hallmark_GSEA_running_curves.png")
pdf_file <- file.path(figure_dir, "Figure3B_mouse_Hallmark_GSEA_running_curves.pdf")

## =========================
## 1. Clean only Step13 outputs
## =========================

old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  png_file,
  pdf_file,
  log_file,
  archive_file
)
unlink(old_outputs[file.exists(old_outputs)], recursive = TRUE, force = TRUE)

## =========================
## 2. Start log and helper functions
## =========================

sink(log_file, split = TRUE)
on.exit({
  sink()
}, add = TRUE)

cat("===== STEP13 FIGURE3B MOUSE HALLMARK GSEA FROM LOCKED STEP05 DE =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Base directory: ", base_dir, "\n", sep = "")
cat("Locked 1W DE input: ", de_1w_file, "\n", sep = "")
cat("Locked 4W DE input: ", de_4w_file, "\n", sep = "")
cat("Rank statistic: sign(logFC) * -log10(P.Value_for_rank)\n")
cat("P-value handling: P.Value <= 0 or non-finite P.Value is replaced with .Machine$double.xmin before rank calculation.\n")
cat("fgsea parameters: minSize = ", fgsea_min_size, ", maxSize = ", fgsea_max_size, ", eps = ", fgsea_eps, ", nproc = ", fgsea_nproc, "\n", sep = "")
cat("fgsea fixed seeds: 1W = ", fgsea_seed_1w, ", 4W = ", fgsea_seed_4w, "\n\n", sep = "")
cat("Important: no fuzzy DE-file search is used in this locked revision.\n\n")

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

archive_current_script <- function(path) {
  src <- get_script_path()
  if (!is.na(src) && file.exists(src)) {
    writeLines(readLines(src, warn = FALSE, encoding = "UTF-8"), path, useBytes = TRUE)
    cat("Script archived to: ", path, "\n", sep = "")
  } else {
    writeLines(
      c(
        "# Step13 archive fallback.",
        "# R could not detect the executed script path.",
        "# Please keep the externally executed script as the authoritative copy."
      ),
      path,
      useBytes = TRUE
    )
    cat("Script path not detected; archive fallback saved to: ", path, "\n", sep = "")
  }
}

safe_library <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (cc in colnames(x)) {
    if (is.list(x[[cc]])) {
      x[[cc]] <- vapply(
        x[[cc]],
        function(v) paste(as.character(v), collapse = ";"),
        character(1)
      )
    }
  }
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

pick_col <- function(df, candidates, required = TRUE, what = NULL) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Cannot find required column",
      if (!is.null(what)) paste0(" for ", what) else "",
      ". Checked: ", paste(candidates, collapse = ", ")
    )
  }
  NA_character_
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing required ", label, ": ", path)
  }
  invisible(TRUE)
}

make_rank_table <- function(de_file, timepoint) {
  require_file(de_file, paste0(timepoint, " locked Step05 DE file"))
  de <- read.csv(de_file, stringsAsFactors = FALSE, check.names = FALSE)

  logfc_col <- pick_col(de, c("logFC", "log2FC", "log2FoldChange"), what = paste0(timepoint, " logFC"))
  p_col <- pick_col(de, c("P.Value", "PValue", "P.Value.", "pvalue", "P"), what = paste0(timepoint, " P value"))
  symbol_col <- pick_col(de, c("SYMBOL", "gene_symbol", "mouse_symbol", "GeneSymbol", "symbol", "Symbol"), what = paste0(timepoint, " mouse gene symbol"))
  fdr_col <- pick_col(de, c("FDR", "adj.P.Val", "padj", "adj_p", "adj.P.Val.", "qvalue", "q_value"), required = FALSE)

  rank_df <- data.frame(
    timepoint = timepoint,
    gene_symbol = as.character(de[[symbol_col]]),
    logFC = as.numeric(de[[logfc_col]]),
    P.Value = as.numeric(de[[p_col]]),
    FDR = if (!is.na(fdr_col)) as.numeric(de[[fdr_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )

  rank_df$gene_symbol <- trimws(rank_df$gene_symbol)

  rank_df <- rank_df[
    !is.na(rank_df$gene_symbol) & rank_df$gene_symbol != "" &
      is.finite(rank_df$logFC),
    ,
    drop = FALSE
  ]

  rank_df$P.Value_raw <- rank_df$P.Value
  rank_df$P.Value_for_rank <- rank_df$P.Value_raw
  rank_df$p_value_replaced_for_rank <- !is.finite(rank_df$P.Value_for_rank) | rank_df$P.Value_for_rank <= 0
  rank_df$P.Value_for_rank[rank_df$p_value_replaced_for_rank] <- .Machine$double.xmin

  p_replaced_df <- rank_df[rank_df$p_value_replaced_for_rank, , drop = FALSE]
  if (nrow(p_replaced_df) > 0) {
    write_csv(p_replaced_df, file.path(table_dir, paste0("step13_p_values_replaced_for_rank_", timepoint, ".csv")))
  }

  rank_df$rank_stat <- sign(rank_df$logFC) * -log10(rank_df$P.Value_for_rank)

  ## If multiple gene-level rows map to the same SYMBOL, keep the row with the largest
  ## absolute rank statistic so that fgsea receives one score per gene symbol.
  rank_df <- rank_df[order(abs(rank_df$rank_stat), decreasing = TRUE), , drop = FALSE]

  duplicated_symbols <- rank_df$gene_symbol[duplicated(rank_df$gene_symbol)]
  if (length(duplicated_symbols) > 0) {
    dup_df <- rank_df[rank_df$gene_symbol %in% duplicated_symbols, , drop = FALSE]
    write_csv(dup_df, file.path(table_dir, paste0("step13_duplicate_gene_symbols_removed_", timepoint, ".csv")))
  }

  rank_df <- rank_df[!duplicated(rank_df$gene_symbol), , drop = FALSE]
  rank_df <- rank_df[order(rank_df$rank_stat, decreasing = TRUE), , drop = FALSE]
  rank_df$rank_order <- seq_len(nrow(rank_df))

  rank_df
}

get_msigdb_hallmark_mouse <- function() {
  out <- tryCatch(
    msigdbr::msigdbr(species = "Mus musculus", collection = "H"),
    error = function(e1) {
      tryCatch(
        msigdbr::msigdbr(species = "Mus musculus", category = "H"),
        error = function(e2) {
          stop(
            "Failed to retrieve Mus musculus Hallmark gene sets from msigdbr using both collection='H' and category='H'.\n",
            "collection error: ", conditionMessage(e1), "\n",
            "category error: ", conditionMessage(e2)
          )
        }
      )
    }
  )
  out
}

enrichment_curve_df <- function(stats, pathway_genes, pathway_name, timepoint, nes, padj, pval) {
  stats <- sort(stats, decreasing = TRUE)
  pathway_genes <- unique(pathway_genes)
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
    timepoint = timepoint,
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
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

archive_current_script(archive_file)

## =========================
## 3. Packages and locked input audit
## =========================

safe_library("fgsea", bioc = TRUE)
safe_library("msigdbr")
safe_library("ggplot2")
safe_library("dplyr")

require_file(de_1w_file, "1W locked Step05 DE file")
require_file(de_4w_file, "4W locked Step05 DE file")

input_audit <- data.frame(
  input = c("de_1w_file", "de_4w_file"),
  path = c(de_1w_file, de_4w_file),
  exists = c(file.exists(de_1w_file), file.exists(de_4w_file)),
  size_bytes = c(file.info(de_1w_file)$size, file.info(de_4w_file)$size),
  stringsAsFactors = FALSE
)
write_csv(input_audit, file.path(table_dir, "step13_locked_DE_input_file_audit.csv"))

cat("Using exact locked DE files:\n")
print(input_audit)

## =========================
## 4. Rank construction
## =========================

rank_1w <- make_rank_table(de_1w_file, "1W")
rank_4w <- make_rank_table(de_4w_file, "4W")

write_csv(rank_1w, file.path(table_dir, "step13_mouse_GSEA_rank_1W.csv"))
write_csv(rank_4w, file.path(table_dir, "step13_mouse_GSEA_rank_4W.csv"))

ranks_1w <- rank_1w$rank_stat
names(ranks_1w) <- rank_1w$gene_symbol
ranks_1w <- sort(ranks_1w, decreasing = TRUE)

ranks_4w <- rank_4w$rank_stat
names(ranks_4w) <- rank_4w$gene_symbol
ranks_4w <- sort(ranks_4w, decreasing = TRUE)

## =========================
## 5. Hallmark gene sets and fgsea
## =========================

msig_raw <- get_msigdb_hallmark_mouse()

pathway_col <- pick_col(msig_raw, c("gs_name", "gs_exact_source"), what = "msigdbr pathway name")
gene_col <- pick_col(msig_raw, c("gene_symbol", "mouse_gene_symbol", "db_gene_symbol"), what = "msigdbr gene symbol")

pathways <- split(msig_raw[[gene_col]], msig_raw[[pathway_col]])
pathways <- lapply(pathways, unique)

hallmark_source <- data.frame(
  pathway = rep(names(pathways), lengths(pathways)),
  gene_symbol = unlist(pathways, use.names = FALSE),
  stringsAsFactors = FALSE
)
write_csv(hallmark_source, file.path(table_dir, "step13_mouse_Hallmark_gene_sets_used.csv"))

set.seed(fgsea_seed_1w)
fgsea_1w <- as.data.frame(fgsea::fgsea(
  pathways = pathways,
  stats = ranks_1w,
  minSize = fgsea_min_size,
  maxSize = fgsea_max_size,
  eps = fgsea_eps,
  nproc = fgsea_nproc
))

set.seed(fgsea_seed_4w)
fgsea_4w <- as.data.frame(fgsea::fgsea(
  pathways = pathways,
  stats = ranks_4w,
  minSize = fgsea_min_size,
  maxSize = fgsea_max_size,
  eps = fgsea_eps,
  nproc = fgsea_nproc
))

fgsea_1w$timepoint <- "1W"
fgsea_4w$timepoint <- "4W"

fgsea_all <- rbind(fgsea_1w, fgsea_4w)
fgsea_all <- fgsea_all[, c("timepoint", setdiff(colnames(fgsea_all), "timepoint"))]

write_csv(fgsea_1w, file.path(table_dir, "step13_fgsea_Hallmark_1W.csv"))
write_csv(fgsea_4w, file.path(table_dir, "step13_fgsea_Hallmark_4W.csv"))
write_csv(fgsea_all, file.path(table_dir, "step13_fgsea_Hallmark_all.csv"))

target_stats <- fgsea_all[fgsea_all$pathway %in% target_pathways, , drop = FALSE]
target_stats <- target_stats[order(match(target_stats$pathway, target_pathways), target_stats$timepoint), , drop = FALSE]
write_csv(target_stats, file.path(table_dir, "step13_Figure3B_target_pathway_stats.csv"))

missing_targets <- setdiff(target_pathways, unique(target_stats$pathway))
if (length(missing_targets) > 0) {
  stop("Some target pathways were not found in fgsea results: ", paste(missing_targets, collapse = ", "))
}

## =========================
## 6. Running-curve source data
## =========================

curve_list <- list()

for (tp in c("1W", "4W")) {
  ranks <- if (tp == "1W") ranks_1w else ranks_4w
  fg <- fgsea_all[fgsea_all$timepoint == tp, , drop = FALSE]

  for (pw in target_pathways) {
    row <- fg[fg$pathway == pw, , drop = FALSE]
    curve_list[[paste(tp, pw, sep = "__")]] <- enrichment_curve_df(
      stats = ranks,
      pathway_genes = pathways[[pw]],
      pathway_name = pw,
      timepoint = tp,
      nes = row$NES[1],
      padj = row$padj[1],
      pval = row$pval[1]
    )
  }
}

curve_df <- do.call(rbind, curve_list)
rownames(curve_df) <- NULL

if (is.null(curve_df) || nrow(curve_df) == 0) {
  stop("No running-curve source data could be generated for the target pathways.")
}

curve_df$timepoint <- factor(curve_df$timepoint, levels = c("1W", "4W"))
curve_df$pathway <- factor(curve_df$pathway, levels = target_pathways)

write_csv(curve_df, file.path(table_dir, "step13_Figure3B_running_curve_source_data.csv"))

## =========================
## 7. Plot Figure 3B
## =========================

label_df <- target_stats
label_df$pathway_pretty <- pretty_pathway(label_df$pathway)
label_df$timepoint <- factor(label_df$timepoint, levels = c("1W", "4W"))
label_df$panel_label <- paste0(
  label_df$pathway_pretty,
  " — ",
  label_df$timepoint,
  "\nNES = ",
  sprintf("%.3f", label_df$NES),
  "   |   FDR = ",
  format_fdr(label_df$padj)
)

panel_order <- c(
  label_df$panel_label[label_df$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" &
                         as.character(label_df$timepoint) == "1W"],
  label_df$panel_label[label_df$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" &
                         as.character(label_df$timepoint) == "4W"],
  label_df$panel_label[label_df$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" &
                         as.character(label_df$timepoint) == "1W"],
  label_df$panel_label[label_df$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" &
                         as.character(label_df$timepoint) == "4W"]
)

curve_df$pathway_pretty <- pretty_pathway(curve_df$pathway)
curve_df$timepoint_chr <- as.character(curve_df$timepoint)
curve_df$panel_key <- paste(curve_df$pathway_pretty, curve_df$timepoint_chr, sep = "__")

label_key <- data.frame(
  pathway = as.character(label_df$pathway),
  timepoint_chr = as.character(label_df$timepoint),
  pathway_pretty = label_df$pathway_pretty,
  panel_label = label_df$panel_label,
  stringsAsFactors = FALSE
)
label_key$panel_key <- paste(label_key$pathway_pretty, label_key$timepoint_chr, sep = "__")

curve_df <- merge(
  curve_df,
  label_key[, c("panel_key", "panel_label")],
  by = "panel_key",
  all.x = TRUE
)
curve_df$panel_label <- factor(curve_df$panel_label, levels = panel_order)

peak_df <- aggregate(running_ES ~ panel_label, data = curve_df, FUN = max)
colnames(peak_df)[2] <- "peak_ES"

rug_df <- curve_df[curve_df$is_hit, , drop = FALSE]
rug_df$rug_ymin <- -0.045
rug_df$rug_ymax <- 0.045

p <- ggplot(curve_df, aes(x = rank, y = running_ES)) +
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
  labs(
    x = NULL,
    y = "Running enrichment score"
  ) +
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

ggsave(png_file, p, width = 11.8, height = 8.0, dpi = 320)
ggsave(pdf_file, p, width = 11.8, height = 8.0)

## =========================
## 8. Summary and reproducibility
## =========================

summary_df <- data.frame(
  metric = c(
    "DE_file_1W",
    "DE_file_4W",
    "DE_file_selection",
    "rank_method",
    "rank_direction",
    "p_value_handling",
    "n_ranked_genes_1W",
    "n_ranked_genes_4W",
    "n_p_values_replaced_1W",
    "n_p_values_replaced_4W",
    "Hallmark_species",
    "Hallmark_collection",
    "Hallmark_pathways_tested",
    "fgsea_minSize",
    "fgsea_maxSize",
    "fgsea_eps",
    "fgsea_nproc",
    "fgsea_seed_1W",
    "fgsea_seed_4W",
    "target_pathways_plotted",
    "multiple_testing_scope"
  ),
  value = as.character(c(
    de_1w_file,
    de_4w_file,
    "exact locked Step05 filenames; no fuzzy search",
    "sign(logFC) * -log10(P.Value_for_rank)",
    "logFC > 0 indicates ACLR > Contra",
    "P.Value <= 0 or non-finite P.Value replaced with .Machine$double.xmin before rank calculation",
    length(ranks_1w),
    length(ranks_4w),
    sum(rank_1w$p_value_replaced_for_rank),
    sum(rank_4w$p_value_replaced_for_rank),
    "Mus musculus",
    "MSigDB Hallmark collection H",
    length(pathways),
    fgsea_min_size,
    fgsea_max_size,
    fgsea_eps,
    fgsea_nproc,
    fgsea_seed_1w,
    fgsea_seed_4w,
    paste(target_pathways, collapse = "; "),
    "adjusted P values from full Hallmark collection, not only the two target pathways"
  )),
  stringsAsFactors = FALSE
)

write_csv(summary_df, file.path(table_dir, "step13_Figure3B_GSEA_summary.csv"))
write_csv(
  summary_df[summary_df$metric %in% c(
    "DE_file_selection", "rank_method", "p_value_handling",
    "fgsea_minSize", "fgsea_maxSize", "fgsea_eps",
    "fgsea_nproc", "fgsea_seed_1W", "fgsea_seed_4W",
    "multiple_testing_scope"
  ), ],
  file.path(table_dir, "step13_fgsea_reproducibility_settings.csv")
)

versions <- data.frame(
  item = c("R", "fgsea", "msigdbr", "ggplot2", "dplyr", "RNGkind"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("fgsea")),
    as.character(utils::packageVersion("msigdbr")),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("dplyr")),
    paste(RNGkind(), collapse = "; ")
  ),
  stringsAsFactors = FALSE
)

write_csv(versions, file.path(table_dir, "step13_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step13_sessionInfo.txt"))

cat("\n===== STEP13 SUMMARY =====\n")
cat("\nGSEA summary:\n")
print(summary_df)

cat("\nTarget pathway statistics:\n")
print(target_stats[, c("timepoint", "pathway", "pval", "padj", "ES", "NES", "size")])

cat("\nFigure3B saved to:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")

cat("\nStep13 completed successfully.\n")
