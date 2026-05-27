# Step15: Figure 3D — Top CellAge-overlapping persistent genes by ACLR direction
# English note:
# This script generates Figure 3D using the current rebuilt mouse discovery workflow.
# It uses:
#   1) Step07 persistent direction-consistent mouse gene table;
#   2) Step11C current CellAge-overlap table produced from current Step11A CellAge cleaning
#      and current Step11B gprofiler2::gorth strict 1:1 mouse-to-human remapping.
# It does NOT read the old frozen Figure3C workspace or any archived overlap-count-specific files.
#
# Genes are ranked within ACLR-up and ACLR-down directions by:
#   combined_effect_score = -log10(FDR_1W) * abs(logFC_1W) +
#                           -log10(FDR_4W) * abs(logFC_4W)
# The plot displays the signed mean persistent logFC across 1W and 4W.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths
## =========================
base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

persistent_file <- file.path(
  base_dir, "07_tables", "step07_strict_DEG_upset_persistent",
  "step07_persistent_direction_consistent_genes.csv"
)

overlap_file <- file.path(
  base_dir, "07_tables", "step11C_intersect_current_gorth_mapping_with_CellAge_clean",
  "step11C_persistent_CellAge_overlap_genes.csv"
)

# Optional audit output from Step11D, used only if present.
step11D_audit_file <- file.path(
  base_dir, "07_tables", "step11D_audit_persistent_mapping_denominators",
  "step11D_persistent_mapping_denominator_audit_summary.csv"
)

table_dir <- file.path(base_dir, "07_tables", "step15_Figure3D_top_CellAge_overlap_by_direction_current")
figure_dir <- file.path(base_dir, "06_figures", "Figure3")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

for (d in c(table_dir, figure_dir, log_dir, script_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

step_name <- "step15_Figure3D_top_CellAge_overlap_by_direction_current"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

png_file <- file.path(figure_dir, "Figure3D_top_CellAge_overlapping_persistent_genes_by_ACLR_direction_current.png")
pdf_file <- file.path(figure_dir, "Figure3D_top_CellAge_overlapping_persistent_genes_by_ACLR_direction_current.pdf")

# Canonical manuscript copies. These intentionally overwrite previous Figure3D canonical files.
canonical_png_file <- file.path(figure_dir, "Figure3D_top_CellAge_overlapping_persistent_genes_by_ACLR_direction.png")
canonical_pdf_file <- file.path(figure_dir, "Figure3D_top_CellAge_overlapping_persistent_genes_by_ACLR_direction.pdf")

source_data_file <- file.path(table_dir, "step15_Figure3D_source_data_current.csv")
full_merged_file <- file.path(table_dir, "step15_Figure3D_full_merged_overlap_table_current.csv")
summary_file <- file.path(table_dir, "step15_Figure3D_run_summary_current.csv")
input_audit_file <- file.path(table_dir, "step15_Figure3D_input_file_audit_current.csv")
software_file <- file.path(table_dir, "step15_Figure3D_software_versions_current.csv")
session_file <- file.path(table_dir, "step15_Figure3D_sessionInfo_current.txt")
archive_file <- file.path(script_dir, "step15_Figure3D_top_CellAge_overlap_by_direction_current.R")

old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  png_file, pdf_file, canonical_png_file, canonical_pdf_file,
  log_file, software_file, session_file
)
unlink(old_outputs[file.exists(old_outputs)], recursive = TRUE, force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP15 FIGURE3D CURRENT =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Persistent input: ", persistent_file, "\n", sep = "")
cat("Current CellAge overlap input: ", overlap_file, "\n", sep = "")
cat("Frozen workspace used: FALSE\n\n")

## =========================
## 1. Helper functions and packages
## =========================
safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (cc in colnames(x)) {
    if (is.list(x[[cc]])) {
      x[[cc]] <- vapply(x[[cc]], function(v) paste(as.character(v), collapse = ";"), character(1))
    }
  }
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA
  x
}

file_audit <- function(label, file) {
  exists_now <- file.exists(file)
  data.frame(
    label = label,
    file = file,
    exists = exists_now,
    size_bytes = if (exists_now) file.info(file)$size else NA_real_,
    mtime = if (exists_now) as.character(file.info(file)$mtime) else NA_character_,
    md5 = if (exists_now) as.character(tools::md5sum(file)) else NA_character_,
    stringsAsFactors = FALSE
  )
}

find_col <- function(df, candidates, what) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop("Cannot identify ", what, " column. Checked: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

safe_library("ggplot2")
safe_library("dplyr")

## =========================
## 2. Input checks and audit
## =========================
if (!file.exists(persistent_file)) {
  stop("Missing Step07 persistent table: ", persistent_file)
}
if (!file.exists(overlap_file)) {
  stop(
    "Missing current Step11C CellAge-overlap table: ", overlap_file, "\n",
    "Please run the current Step11C script first."
  )
}

input_audit <- rbind(
  file_audit("Step07 persistent direction-consistent genes", persistent_file),
  file_audit("Step11C current CellAge overlap table", overlap_file),
  file_audit("Step11D denominator audit summary optional", step11D_audit_file)
)
write_csv(input_audit, input_audit_file)
cat("\nInput file audit:\n")
print(input_audit)

persistent <- read.csv(persistent_file, stringsAsFactors = FALSE, check.names = FALSE)
overlap <- read.csv(overlap_file, stringsAsFactors = FALSE, check.names = FALSE)

required_persistent_cols <- c(
  "gene_id", "SYMBOL", "logFC_1W", "FDR_1W", "logFC_4W", "FDR_4W", "direction_persistent"
)
miss_persistent <- setdiff(required_persistent_cols, colnames(persistent))
if (length(miss_persistent) > 0) {
  stop("Persistent table is missing required columns: ", paste(miss_persistent, collapse = ", "))
}

mouse_col <- find_col(overlap, c("input", "mouse_symbol", "SYMBOL", "symbol_mouse", "mouse_symbol_clean"), "mouse symbol in overlap table")
human_col <- find_col(overlap, c("ortholog_name", "human_symbol", "Gene symbol", "ortholog_symbol"), "human ortholog symbol in overlap table")
effect_col <- find_col(overlap, c("Senescence Effect", "CellAge_effect", "senescence_effect", "effect"), "CellAge senescence effect")

persistent$mouse_symbol <- clean_symbol(persistent$SYMBOL)
overlap$mouse_symbol <- clean_symbol(overlap[[mouse_col]])
overlap$human_symbol <- clean_symbol(overlap[[human_col]])
overlap$CellAge_effect_current <- clean_symbol(overlap[[effect_col]])

persistent_rows_total <- nrow(persistent)
persistent_rows_with_symbol <- sum(!is.na(persistent$mouse_symbol))
persistent_unique_symbols <- length(unique(na.omit(persistent$mouse_symbol)))
overlap_rows_total <- nrow(overlap)
overlap_unique_human <- length(unique(na.omit(overlap$human_symbol)))
overlap_unique_mouse <- length(unique(na.omit(overlap$mouse_symbol)))

cat("\nInput counts:\n")
cat("Step07 persistent rows: ", persistent_rows_total, "\n", sep = "")
cat("Step07 persistent rows with valid SYMBOL: ", persistent_rows_with_symbol, "\n", sep = "")
cat("Step07 persistent unique mouse symbols: ", persistent_unique_symbols, "\n", sep = "")
cat("Step11C overlap rows: ", overlap_rows_total, "\n", sep = "")
cat("Step11C unique mouse symbols: ", overlap_unique_mouse, "\n", sep = "")
cat("Step11C unique human symbols: ", overlap_unique_human, "\n\n", sep = "")

## =========================
## 3. Merge persistent DE statistics with current CellAge overlap table
## =========================
overlap_unique <- overlap %>%
  filter(!is.na(mouse_symbol) & mouse_symbol != "") %>%
  distinct(mouse_symbol, .keep_all = TRUE)

merged <- persistent %>%
  inner_join(overlap_unique, by = "mouse_symbol")

if (nrow(merged) == 0) {
  stop("No genes remained after joining Step07 persistent genes with current Step11C CellAge-overlap table.")
}

## =========================
## 4. Ranking and top-gene selection
## =========================
eps <- 1e-300
merged$FDR_1W_num <- as.numeric(merged$FDR_1W)
merged$FDR_4W_num <- as.numeric(merged$FDR_4W)
merged$logFC_1W_num <- as.numeric(merged$logFC_1W)
merged$logFC_4W_num <- as.numeric(merged$logFC_4W)

merged$FDR_1W_safe <- ifelse(is.na(merged$FDR_1W_num), NA_real_, pmax(merged$FDR_1W_num, eps))
merged$FDR_4W_safe <- ifelse(is.na(merged$FDR_4W_num), NA_real_, pmax(merged$FDR_4W_num, eps))
merged$rank_component_1W <- -log10(merged$FDR_1W_safe) * abs(merged$logFC_1W_num)
merged$rank_component_4W <- -log10(merged$FDR_4W_safe) * abs(merged$logFC_4W_num)
merged$combined_effect_score <- merged$rank_component_1W + merged$rank_component_4W
merged$mean_logFC <- (merged$logFC_1W_num + merged$logFC_4W_num) / 2

merged$plot_direction <- ifelse(
  merged$direction_persistent == "Up", "Up in ACLR",
  ifelse(merged$direction_persistent == "Down", "Down in ACLR", NA)
)
merged <- merged[!is.na(merged$plot_direction), , drop = FALSE]

if (nrow(merged) == 0) {
  stop("Merged CellAge-overlap genes do not contain Up/Down direction_persistent values.")
}

top_n_each <- 10

top_up <- merged %>%
  filter(plot_direction == "Up in ACLR") %>%
  arrange(desc(combined_effect_score), desc(mean_logFC), human_symbol, mouse_symbol) %>%
  slice_head(n = top_n_each)

top_down <- merged %>%
  filter(plot_direction == "Down in ACLR") %>%
  arrange(desc(combined_effect_score), desc(abs(mean_logFC)), human_symbol, mouse_symbol) %>%
  slice_head(n = top_n_each)

plot_df <- bind_rows(top_up, top_down)
if (nrow(plot_df) == 0) {
  stop("No genes available for Figure3D plotting.")
}

plot_df$display_gene <- ifelse(!is.na(plot_df$human_symbol) & plot_df$human_symbol != "",
                               plot_df$human_symbol, plot_df$mouse_symbol)
plot_df$plot_value <- ifelse(plot_df$plot_direction == "Up in ACLR",
                             abs(plot_df$mean_logFC), -abs(plot_df$mean_logFC))

plot_df$plot_direction <- factor(plot_df$plot_direction,
                                 levels = c("Up in ACLR", "Down in ACLR"))

up_order_top_to_bottom <- plot_df %>%
  filter(plot_direction == "Up in ACLR") %>%
  arrange(desc(plot_value), desc(combined_effect_score), display_gene) %>%
  pull(display_gene)

down_order_top_to_bottom <- plot_df %>%
  filter(plot_direction == "Down in ACLR") %>%
  arrange(plot_value, desc(combined_effect_score), display_gene) %>%
  pull(display_gene)

plot_df$facet_gene <- ifelse(
  plot_df$plot_direction == "Up in ACLR",
  paste0(plot_df$display_gene, "___UP"),
  paste0(plot_df$display_gene, "___DOWN")
)
facet_levels <- c(paste0(rev(up_order_top_to_bottom), "___UP"),
                  paste0(rev(down_order_top_to_bottom), "___DOWN"))
plot_df$facet_gene <- factor(plot_df$facet_gene, levels = facet_levels)

## =========================
## 5. Save source data
## =========================
source_df <- plot_df %>%
  mutate(display_gene = sub("___(UP|DOWN)$", "", as.character(facet_gene))) %>%
  select(
    display_gene, mouse_symbol, human_symbol, gene_id,
    plot_direction, CellAge_effect_current,
    logFC_1W, FDR_1W, logFC_4W, FDR_4W,
    mean_logFC, rank_component_1W, rank_component_4W, combined_effect_score,
    plot_value
  )

write_csv(merged, full_merged_file)
write_csv(source_df, source_data_file)

## =========================
## 6. Plot Figure 3D
## =========================
xlim_max <- max(abs(plot_df$plot_value), na.rm = TRUE) * 1.25
if (!is.finite(xlim_max) || xlim_max <= 0) xlim_max <- 1

p <- ggplot(plot_df, aes(x = plot_value, y = facet_gene, fill = plot_direction)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey30") +
  facet_wrap(~ plot_direction, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Up in ACLR" = "#D95F02", "Down in ACLR" = "#1F78B4")) +
  scale_y_discrete(labels = function(x) sub("___(UP|DOWN)$", "", x)) +
  scale_x_continuous(
    limits = c(-xlim_max, xlim_max),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    title = "Top CellAge-overlapping persistent genes by ACLR direction",
    x = "Mean persistent logFC across 1W and 4W",
    y = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 19, margin = margin(b = 10)),
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 12, color = "grey20"),
    axis.text.y = element_text(size = 12, color = "grey10"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey88"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(10, 12, 10, 12)
  )

ggsave(png_file, p, width = 11.5, height = 7.6, dpi = 320)
ggsave(pdf_file, p, width = 11.5, height = 7.6)
file.copy(png_file, canonical_png_file, overwrite = TRUE)
file.copy(pdf_file, canonical_pdf_file, overwrite = TRUE)

cat("Saved: ", png_file, "\n", sep = "")
cat("Saved: ", pdf_file, "\n", sep = "")
cat("Saved canonical copy: ", canonical_png_file, "\n", sep = "")
cat("Saved canonical copy: ", canonical_pdf_file, "\n", sep = "")

## =========================
## 7. Summary, versions, and script archive
## =========================
summary_df <- data.frame(
  metric = c(
    "persistent_input_file",
    "current_overlap_input_file",
    "frozen_workspace_used",
    "Step07_persistent_table_rows",
    "Step07_persistent_rows_with_valid_SYMBOL",
    "Step07_persistent_unique_mouse_symbols",
    "Step11C_overlap_rows",
    "Step11C_unique_mouse_symbols",
    "Step11C_unique_human_symbols",
    "n_CellAge_overlap_genes_merged_with_Step07",
    "n_up_in_ACLR_available",
    "n_down_in_ACLR_available",
    "top_n_each_direction_requested",
    "n_genes_plotted_total",
    "n_up_plotted",
    "n_down_plotted",
    "ranking_metric",
    "plot_metric",
    "facet_order",
    "within_panel_order_rule",
    "figure_png",
    "figure_pdf",
    "canonical_figure_png",
    "canonical_figure_pdf"
  ),
  value = c(
    persistent_file,
    overlap_file,
    "FALSE",
    persistent_rows_total,
    persistent_rows_with_symbol,
    persistent_unique_symbols,
    overlap_rows_total,
    overlap_unique_mouse,
    overlap_unique_human,
    nrow(merged),
    sum(merged$plot_direction == "Up in ACLR"),
    sum(merged$plot_direction == "Down in ACLR"),
    top_n_each,
    nrow(plot_df),
    sum(plot_df$plot_direction == "Up in ACLR"),
    sum(plot_df$plot_direction == "Down in ACLR"),
    "combined_effect_score = -log10(FDR_1W)*abs(logFC_1W) + -log10(FDR_4W)*abs(logFC_4W)",
    "plot_value = signed mean logFC across 1W and 4W",
    "left = Up in ACLR; right = Down in ACLR",
    "top-to-bottom ordered from stronger effect to weaker effect within each facet",
    png_file,
    pdf_file,
    canonical_png_file,
    canonical_pdf_file
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_df, summary_file)

software_versions <- data.frame(
  item = c("R", "ggplot2", "dplyr"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("dplyr"))
  ),
  stringsAsFactors = FALSE
)
write_csv(software_versions, software_file)
writeLines(capture.output(sessionInfo()), session_file)
cat("Saved: ", session_file, "\n", sep = "")

# Archive this script when run with Rscript; otherwise write a clear fallback note.
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }
  return(NA_character_)
}
script_path <- get_script_path()
if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path, warn = FALSE), archive_file, useBytes = TRUE)
} else {
  writeLines("# Script path not detected. Please manually save the executed Step15 script here.", archive_file)
}
cat("Saved script archive/fallback: ", archive_file, "\n", sep = "")

cat("\n===== STEP15 CURRENT SUMMARY =====\n")
print(summary_df)

cat("\nTop up genes plotted, top-to-bottom:\n")
print(data.frame(order_top_to_bottom = seq_along(up_order_top_to_bottom), gene = up_order_top_to_bottom), row.names = FALSE)

cat("\nTop down genes plotted, top-to-bottom:\n")
print(data.frame(order_top_to_bottom = seq_along(down_order_top_to_bottom), gene = down_order_top_to_bottom), row.names = FALSE)

cat("\nStep15 Figure3D current completed successfully.\n")
cat("Main figure:\n", png_file, "\n", pdf_file, "\n", sep = "")
cat("Canonical manuscript figure:\n", canonical_png_file, "\n", canonical_pdf_file, "\n", sep = "")

sink()

cat("\nStep15 Figure3D current completed. Please open:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
