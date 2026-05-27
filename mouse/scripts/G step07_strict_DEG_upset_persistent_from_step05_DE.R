# Step07: Construct mouse strict DEG UpSet structure and persistent gene table.
# Purpose:
# This locked-method Step07 reads the locked Step05 duplicateCorrelation limma-voom
# DE tables, defines strict DEGs as FDR < 0.05 and |logFC| > 1, calculates the 1W-only,
# 4W-only, shared-any-direction, shared-direction-consistent persistent genes, and
# shared-direction-discordant genes. Persistent genes are defined exactly as genes
# that are strict DEGs at both 1W and 4W with consistent logFC direction. The script
# generates a manuscript-ready UpSet-style Figure2A, saves all source data, summaries,
# logs, software versions, and archives this script for reproducibility.

## =========================
## 0. Define paths and locked parameters
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

de_1w_file <- file.path(base_dir, "03_DE_analysis", "step05_DE_1W_ACLR_vs_Contra_limma_voom.csv")
de_4w_file <- file.path(base_dir, "03_DE_analysis", "step05_DE_4W_ACLR_vs_Contra_limma_voom.csv")

figure_dir <- file.path(base_dir, "06_figures", "Figure2")
table_dir <- file.path(base_dir, "07_tables", "step07_strict_DEG_upset_persistent")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step07_strict_DEG_upset_persistent"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

strict_FDR <- 0.05
strict_abs_logFC <- 1

## =========================
## 1. Delete old Step07 outputs
## =========================

old_files <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  file.path(figure_dir, "Figure2A_mouse_strict_DEG_UpSet.png"),
  file.path(figure_dir, "Figure2A_mouse_strict_DEG_UpSet.pdf"),
  file.path(log_dir, paste0(step_name, "_log.txt")),
  file.path(script_dir, "step07_strict_DEG_upset_persistent.R")
)
unlink(old_files[file.exists(old_files)], force = TRUE)

## =========================
## 2. Start log and helper functions
## =========================

sink(log_file, split = TRUE)

cat("===== STEP07 STRICT DEG UPSET + PERSISTENT STRUCTURE =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Strict threshold: FDR < ", strict_FDR, " and |logFC| > ", strict_abs_logFC, "\n", sep = "")
cat("Persistent definition: strict at both 1W and 4W with consistent logFC direction.\n\n")

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
    writeLines("# Step07 archive fallback.", path, useBytes = TRUE)
    cat("Script path not detected; archive fallback saved.\n")
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

find_fdr_col <- function(df) {
  candidates <- c("FDR", "adj.P.Val", "adj_P_Val", "padj", "qvalue")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) stop("No FDR-like column found. Columns: ", paste(colnames(df), collapse = ", "))
  hit[1]
}

standardize_de <- function(df, timepoint) {
  fdr_col <- find_fdr_col(df)

  if (!"gene_id" %in% colnames(df)) stop(timepoint, ": missing gene_id column.")
  if (!"logFC" %in% colnames(df)) stop(timepoint, ": missing logFC column.")

  symbol_col <- if ("SYMBOL" %in% colnames(df)) "SYMBOL" else NA_character_
  entrez_col <- if ("ENTREZID" %in% colnames(df)) "ENTREZID" else NA_character_

  out <- data.frame(
    gene_id = as.character(df$gene_id),
    ENTREZID = if (!is.na(entrez_col)) as.character(df[[entrez_col]]) else "",
    SYMBOL = if (!is.na(symbol_col)) as.character(df[[symbol_col]]) else "",
    logFC = as.numeric(df$logFC),
    FDR = as.numeric(df[[fdr_col]]),
    timepoint = timepoint,
    stringsAsFactors = FALSE
  )

  out$ENTREZID[is.na(out$ENTREZID)] <- ""
  out$SYMBOL[is.na(out$SYMBOL)] <- ""
  out$strict_DEG <- out$FDR < strict_FDR & abs(out$logFC) > strict_abs_logFC
  out$direction <- ifelse(out$logFC > 0, "Up", "Down")
  out$fdr_source_column <- fdr_col

  out
}

## =========================
## 3. Load packages and archive script
## =========================

safe_library("ggplot2")
safe_library("cowplot")

archive_current_script(file.path(script_dir, "step07_strict_DEG_upset_persistent.R"))

## =========================
## 4. Load locked Step05 DE tables
## =========================

if (!file.exists(de_1w_file)) stop("Missing locked Step05 1W DE table: ", de_1w_file)
if (!file.exists(de_4w_file)) stop("Missing locked Step05 4W DE table: ", de_4w_file)

de1_raw <- read.csv(de_1w_file, stringsAsFactors = FALSE, check.names = FALSE)
de4_raw <- read.csv(de_4w_file, stringsAsFactors = FALSE, check.names = FALSE)

de1 <- standardize_de(de1_raw, "1W")
de4 <- standardize_de(de4_raw, "4W")

write_csv(de1, file.path(table_dir, "step07_DE_1W_standardized.csv"))
write_csv(de4, file.path(table_dir, "step07_DE_4W_standardized.csv"))

## =========================
## 5. Strict DEG and persistent structure
## =========================

s1 <- de1[de1$strict_DEG, c("gene_id", "ENTREZID", "SYMBOL", "logFC", "FDR", "direction")]
s4 <- de4[de4$strict_DEG, c("gene_id", "ENTREZID", "SYMBOL", "logFC", "FDR", "direction")]

colnames(s1) <- c("gene_id", "ENTREZID_1W", "SYMBOL_1W", "logFC_1W", "FDR_1W", "direction_1W")
colnames(s4) <- c("gene_id", "ENTREZID_4W", "SYMBOL_4W", "logFC_4W", "FDR_4W", "direction_4W")

structure <- merge(s1, s4, by = "gene_id", all = TRUE)

structure$in_1W <- !is.na(structure$direction_1W)
structure$in_4W <- !is.na(structure$direction_4W)

structure$direction_relation <- ifelse(
  structure$in_1W & structure$in_4W & structure$direction_1W == structure$direction_4W,
  "same_direction",
  ifelse(structure$in_1W & structure$in_4W & structure$direction_1W != structure$direction_4W,
         "opposite_direction", NA)
)

structure$category <- ifelse(
  structure$in_1W & structure$in_4W & structure$direction_relation == "same_direction",
  "persistent_direction_consistent",
  ifelse(
    structure$in_1W & structure$in_4W & structure$direction_relation == "opposite_direction",
    "shared_direction_discordant",
    ifelse(structure$in_1W & !structure$in_4W, "1W_only", "4W_only")
  )
)

structure$ENTREZID <- ifelse(!is.na(structure$ENTREZID_1W) & structure$ENTREZID_1W != "", structure$ENTREZID_1W, structure$ENTREZID_4W)
structure$SYMBOL <- ifelse(!is.na(structure$SYMBOL_1W) & structure$SYMBOL_1W != "", structure$SYMBOL_1W, structure$SYMBOL_4W)

structure$ENTREZID[is.na(structure$ENTREZID)] <- ""
structure$SYMBOL[is.na(structure$SYMBOL)] <- ""

write_csv(structure, file.path(table_dir, "step07_strict_DEG_structure_gene_table.csv"))

persistent <- structure[structure$category == "persistent_direction_consistent", , drop = FALSE]
persistent$direction_persistent <- persistent$direction_1W

persistent <- persistent[, c(
  "gene_id", "ENTREZID", "SYMBOL",
  "logFC_1W", "FDR_1W", "direction_1W",
  "logFC_4W", "FDR_4W", "direction_4W",
  "direction_persistent", "category"
)]

write_csv(persistent, file.path(table_dir, "step07_persistent_direction_consistent_genes.csv"))

persistent_up <- persistent[persistent$direction_persistent == "Up", , drop = FALSE]
persistent_down <- persistent[persistent$direction_persistent == "Down", , drop = FALSE]

write_csv(persistent_up, file.path(table_dir, "step07_persistent_up_genes.csv"))
write_csv(persistent_down, file.path(table_dir, "step07_persistent_down_genes.csv"))

## =========================
## 6. Summary tables
## =========================

strict_summary <- data.frame(
  metric = c(
    "tested_genes_1W",
    "tested_genes_4W",
    "strict_1W_total",
    "strict_1W_up",
    "strict_1W_down",
    "strict_4W_total",
    "strict_4W_up",
    "strict_4W_down",
    "shared_strict_any_direction",
    "shared_strict_same_direction_persistent",
    "shared_strict_opposite_direction",
    "strict_1W_only",
    "strict_4W_only",
    "persistent_up",
    "persistent_down"
  ),
  n = c(
    nrow(de1),
    nrow(de4),
    nrow(s1),
    sum(s1$direction_1W == "Up"),
    sum(s1$direction_1W == "Down"),
    nrow(s4),
    sum(s4$direction_4W == "Up"),
    sum(s4$direction_4W == "Down"),
    sum(structure$in_1W & structure$in_4W),
    nrow(persistent),
    sum(structure$category == "shared_direction_discordant"),
    sum(structure$category == "1W_only"),
    sum(structure$category == "4W_only"),
    nrow(persistent_up),
    nrow(persistent_down)
  ),
  stringsAsFactors = FALSE
)

write_csv(strict_summary, file.path(table_dir, "step07_strict_DEG_persistent_summary.csv"))

category_counts <- as.data.frame(table(structure$category), stringsAsFactors = FALSE)
colnames(category_counts) <- c("category", "n")
category_counts$category <- as.character(category_counts$category)
category_counts <- category_counts[order(match(category_counts$category, c(
  "1W_only", "4W_only", "persistent_direction_consistent", "shared_direction_discordant"
))), ]

write_csv(category_counts, file.path(table_dir, "step07_structure_category_counts.csv"))

## =========================
## 7. Figure2A UpSet-style plot
## =========================

# Use a reproducible custom UpSet-style display:
# top bars = intersection sizes; bottom matrix = 1W/4W membership.
upset_df <- data.frame(
  intersection = c("1W only", "1W ∩ 4W", "4W only"),
  n = c(
    sum(structure$category == "1W_only"),
    sum(structure$in_1W & structure$in_4W),
    sum(structure$category == "4W_only")
  ),
  in_1W = c(TRUE, TRUE, FALSE),
  in_4W = c(FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

upset_df$intersection <- factor(upset_df$intersection, levels = c("1W only", "1W ∩ 4W", "4W only"))

write_csv(upset_df, file.path(table_dir, "step07_Figure2A_upset_source_data.csv"))

bar_plot <- ggplot2::ggplot(upset_df, ggplot2::aes(x = intersection, y = n)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.35, size = 4) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(x = NULL, y = "Number of strict DEGs") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

matrix_df <- rbind(
  data.frame(intersection = upset_df$intersection, set = "Strict 1W", active = upset_df$in_1W),
  data.frame(intersection = upset_df$intersection, set = "Strict 4W", active = upset_df$in_4W)
)
matrix_df$set <- factor(matrix_df$set, levels = c("Strict 4W", "Strict 1W"))

matrix_plot <- ggplot2::ggplot(matrix_df, ggplot2::aes(x = intersection, y = set)) +
  ggplot2::geom_point(ggplot2::aes(alpha = active), size = 3) +
  ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.15), guide = "none") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

combined <- cowplot::plot_grid(bar_plot, matrix_plot, ncol = 1, rel_heights = c(3, 1), align = "v")

ggplot2::ggsave(file.path(figure_dir, "Figure2A_mouse_strict_DEG_UpSet.png"), combined, width = 5.8, height = 5.2, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "Figure2A_mouse_strict_DEG_UpSet.pdf"), combined, width = 5.8, height = 5.2)

## =========================
## 8. Software versions
## =========================

software_versions <- data.frame(
  item = c("R", "ggplot2", "cowplot"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("cowplot"))
  ),
  stringsAsFactors = FALSE
)
write_csv(software_versions, file.path(table_dir, "step07_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step07_sessionInfo.txt"))

## =========================
## 9. Console summary
## =========================

cat("\n===== STEP07 SUMMARY =====\n")

cat("\nStrict DEG and persistent summary:\n")
print(strict_summary)

cat("\nStructure category counts:\n")
print(category_counts)

cat("\nPersistent gene preview:\n")
print(head(persistent, 20))

cat("\nFigure2A saved to:\n")
cat(file.path(figure_dir, "Figure2A_mouse_strict_DEG_UpSet.png"), "\n")
cat(file.path(figure_dir, "Figure2A_mouse_strict_DEG_UpSet.pdf"), "\n")

cat("\nStep07 completed successfully.\n")

sink()
