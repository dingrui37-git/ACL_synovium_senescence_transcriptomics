# StepH2w12_current78: fastp-trimmed branch GSEA running-curve figure
# Purpose:
# Generate a Supplementary Figure-style Hallmark GSEA running-curve plot for the
# fastp-trimmed current78 pig early sensitivity branch, using the same visual style
# as the main Step20_current78 Figure4B GSEA plot.
#
# Data source:
#   StepH2w9_current78_GSEA_stability outputs
#
# Figure:
#   - 2 x 2 panels:
#       HALLMARK_INFLAMMATORY_RESPONSE at t7 and t28
#       HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION at t7 and t28
#   - Uses trimmed-branch ranked gene lists and trimmed-branch target pathway statistics
#   - Target pathway FDR values are from full Hallmark collection-level fgsea
#   - Only ONE PNG and ONE PDF are generated
#
# This script does NOT auto-archive or overwrite any manually saved R script.

options(stringsAsFactors = FALSE)

## =========================
## 0. User paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

fastp_dir <- file.path(rebuild_root, "02_pig_early_fastp_sensitivity")

# StepH2w9 output folder. If you used the original successful folder, this is correct.
# If you used the v2_symbolmap new-output folder instead, the script will try that as fallback.
stepH2w9_dir_primary <- file.path(fastp_dir, "tables", "stepH2w9_current78_GSEA_stability")
stepH2w9_dir_fallback <- file.path(fastp_dir, "tables", "stepH2w9_current78_GSEA_stability_v2_symbolmap")
stepH2w9_dir <- if (dir.exists(stepH2w9_dir_primary)) stepH2w9_dir_primary else stepH2w9_dir_fallback

out_fig_dir <- file.path(fastp_dir, "figures", "stepH2w12_current78_trimmed_GSEA_running_curves")
out_table_dir <- file.path(fastp_dir, "tables", "stepH2w12_current78_trimmed_GSEA_running_curves")
logs_dir <- file.path(out_table_dir, "logs")
objects_dir <- file.path(out_table_dir, "objects")

for (d in c(out_fig_dir, out_table_dir, logs_dir, objects_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

full_log_file <- file.path(logs_dir, "stepH2w12_current78_trimmed_GSEA_running_curves_full_log.txt")
summary_log_file <- file.path(logs_dir, "stepH2w12_current78_trimmed_GSEA_running_curves_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== StepH2w12_current78 trimmed GSEA running-curve figure =====\n")
cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("StepH2w9 directory:", stepH2w9_dir, "\n")
cat("Purpose: plot trimmed-branch target Hallmark GSEA running curves in Step20 style.\n\n")

## =========================
## 1. Packages and helpers
## =========================

required_pkgs <- c("ggplot2", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

stop_if_missing <- function(file, label) {
  if (!file.exists(file)) stop("Missing ", label, ": ", file, call. = FALSE)
}

safe_read_csv <- function(file) {
  read.csv(file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(x, file) {
  write.csv(x, file, row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")
  cat("Saved: ", file, "\n", sep = "")
}

find_existing_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

pretty_pathway <- function(x) {
  x <- as.character(x)
  x[x == "HALLMARK_INFLAMMATORY_RESPONSE"] <- "Inflammatory response"
  x[x == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT / ECM remodeling"
  x
}

pretty_contrast <- function(x) {
  x <- as.character(x)
  x[x == "ACLT_untreated_t7_vs_Control"] <- "t7"
  x[x == "ACLT_untreated_t28_vs_Control"] <- "t28"
  x
}

format_fdr <- function(x) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "e", digits = 2))
}

enrichment_curve_df <- function(stats, pathway_genes, pathway_name, contrast_label, nes, padj, pval) {
  stats <- sort(stats, decreasing = TRUE)
  hits <- names(stats) %in% pathway_genes
  N <- length(stats)
  Nh <- sum(hits)
  Nm <- N - Nh
  if (Nh == 0 || Nm == 0) {
    return(data.frame())
  }

  weights <- abs(stats)
  Phit <- ifelse(hits, weights / sum(weights[hits]), 0)
  Pmiss <- ifelse(!hits, 1 / Nm, 0)
  running <- cumsum(Phit - Pmiss)

  data.frame(
    branch = "trimmed",
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

## =========================
## 2. Input files from StepH2w9
## =========================

target_results_file <- file.path(stepH2w9_dir, "stepH2w9_current78_target_pathway_results_all_branches.csv")
pathways_rds_file <- file.path(stepH2w9_dir, "objects", "stepH2w9_current78_sus_scrofa_hallmark_pathways.rds")
rank_t7_file <- file.path(stepH2w9_dir, "stepH2w9_current78_ranked_genes_trimmed__t7.csv")
rank_t28_file <- file.path(stepH2w9_dir, "stepH2w9_current78_ranked_genes_trimmed__t28.csv")
summary_by_timepoint_file <- file.path(stepH2w9_dir, "stepH2w9_current78_GSEA_stability_summary_by_timepoint.csv")

stop_if_missing(target_results_file, "StepH2w9 target pathway results")
stop_if_missing(pathways_rds_file, "StepH2w9 Sus scrofa Hallmark pathway RDS")
stop_if_missing(rank_t7_file, "StepH2w9 trimmed t7 ranked gene list")
stop_if_missing(rank_t28_file, "StepH2w9 trimmed t28 ranked gene list")

target_results_all <- safe_read_csv(target_results_file)
pathways <- readRDS(pathways_rds_file)
rank_t7_df <- safe_read_csv(rank_t7_file)
rank_t28_df <- safe_read_csv(rank_t28_file)
stability_summary <- if (file.exists(summary_by_timepoint_file)) safe_read_csv(summary_by_timepoint_file) else data.frame()

target_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

missing_targets <- setdiff(target_pathways, names(pathways))
if (length(missing_targets) > 0) {
  stop("Target Hallmark pathways missing from pathway RDS: ", paste(missing_targets, collapse = ", "), call. = FALSE)
}

required_rank_cols <- c("gene_symbol", "rank_stat")
if (!all(required_rank_cols %in% colnames(rank_t7_df))) {
  stop("Trimmed t7 ranked list lacks required columns: ", paste(setdiff(required_rank_cols, colnames(rank_t7_df)), collapse = ", "), call. = FALSE)
}
if (!all(required_rank_cols %in% colnames(rank_t28_df))) {
  stop("Trimmed t28 ranked list lacks required columns: ", paste(setdiff(required_rank_cols, colnames(rank_t28_df)), collapse = ", "), call. = FALSE)
}

rank_to_stats <- function(rank_df) {
  stats <- as.numeric(rank_df$rank_stat)
  names(stats) <- as.character(rank_df$gene_symbol)
  keep <- !is.na(stats) & is.finite(stats) & !is.na(names(stats)) & names(stats) != ""
  stats <- stats[keep]
  stats <- stats[!duplicated(names(stats))]
  sort(stats, decreasing = TRUE)
}

stats_t7 <- rank_to_stats(rank_t7_df)
stats_t28 <- rank_to_stats(rank_t28_df)

cat("Trimmed ranked symbols:\n")
cat("t7 :", length(stats_t7), "\n")
cat("t28:", length(stats_t28), "\n\n")

target_results <- target_results_all[
  target_results_all$branch == "trimmed" &
    target_results_all$pathway %in% target_pathways &
    target_results_all$contrast %in% c("ACLT_untreated_t7_vs_Control", "ACLT_untreated_t28_vs_Control"),
  ,
  drop = FALSE
]

if (nrow(target_results) != 4) {
  stop("Expected 4 trimmed target pathway result rows, but found: ", nrow(target_results), call. = FALSE)
}

## =========================
## 3. Running curve source data
## =========================

target_running <- list()

for (contrast_label in c("ACLT_untreated_t7_vs_Control", "ACLT_untreated_t28_vs_Control")) {
  stats_vec <- if (contrast_label == "ACLT_untreated_t7_vs_Control") stats_t7 else stats_t28
  fg_this <- target_results[target_results$contrast == contrast_label, , drop = FALSE]

  for (pw in target_pathways) {
    row <- fg_this[fg_this$pathway == pw, , drop = FALSE]
    if (nrow(row) != 1) {
      stop("Target pathway not found exactly once: ", contrast_label, " / ", pw, call. = FALSE)
    }

    target_running[[paste(contrast_label, pw, sep = "__")]] <- enrichment_curve_df(
      stats = stats_vec,
      pathway_genes = pathways[[pw]],
      pathway_name = pw,
      contrast_label = contrast_label,
      nes = row$NES[1],
      padj = row$padj[1],
      pval = row$pval[1]
    )
  }
}

running_df <- dplyr::bind_rows(target_running)

if (nrow(running_df) == 0) {
  stop("Running curve source data is empty. Please check ranked gene symbols and pathway membership.", call. = FALSE)
}

## =========================
## 4. Panel labels and plot data
## =========================

label_df <- target_results
label_df$pathway_pretty <- pretty_pathway(label_df$pathway)
label_df$timepoint <- pretty_contrast(label_df$contrast)
label_df$timepoint <- factor(label_df$timepoint, levels = c("t7", "t28"))
label_df$panel_label <- paste0(
  label_df$pathway_pretty,
  " — ",
  label_df$timepoint,
  "\nNES = ",
  sprintf("%.3f", as.numeric(label_df$NES)),
  "   |   FDR = ",
  format_fdr(label_df$padj)
)

panel_order <- c(
  label_df$panel_label[label_df$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & as.character(label_df$timepoint) == "t7"],
  label_df$panel_label[label_df$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & as.character(label_df$timepoint) == "t28"],
  label_df$panel_label[label_df$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & as.character(label_df$timepoint) == "t7"],
  label_df$panel_label[label_df$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & as.character(label_df$timepoint) == "t28"]
)
panel_order <- panel_order[nchar(panel_order) > 0]

running_df$pathway_pretty <- pretty_pathway(running_df$pathway)
running_df$timepoint <- pretty_contrast(running_df$contrast)
running_df$panel_key <- paste(running_df$pathway_pretty, running_df$timepoint, sep = "__")

label_key <- data.frame(
  pathway = as.character(label_df$pathway),
  timepoint = as.character(label_df$timepoint),
  pathway_pretty = label_df$pathway_pretty,
  panel_label = label_df$panel_label,
  stringsAsFactors = FALSE
)
label_key$panel_key <- paste(label_key$pathway_pretty, label_key$timepoint, sep = "__")

running_df <- merge(
  running_df,
  label_key[, c("panel_key", "panel_label")],
  by = "panel_key",
  all.x = TRUE
)
running_df$panel_label <- factor(running_df$panel_label, levels = panel_order)

peak_df <- aggregate(running_ES ~ panel_label, data = running_df, FUN = max)
colnames(peak_df)[2] <- "peak_ES"

rug_df <- running_df[running_df$is_hit, , drop = FALSE]
rug_df$rug_ymin <- -0.045
rug_df$rug_ymax <- 0.045

## =========================
## 5. Plot in Step20 mouse-style green running-curve format
## =========================

p_combined <- ggplot2::ggplot(running_df, ggplot2::aes(x = rank, y = running_ES)) +
  ggplot2::geom_hline(
    data = peak_df,
    ggplot2::aes(yintercept = peak_ES),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = rug_df,
    ggplot2::aes(x = rank, xend = rank, y = rug_ymin, yend = rug_ymax),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.35,
    linewidth = 0.18
  ) +
  ggplot2::geom_line(color = "#00CD00", linewidth = 0.75) +
  ggplot2::facet_wrap(~ panel_label, ncol = 2, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "Running enrichment score") +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 13.5, face = "plain", lineheight = 1.15),
    panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = ggplot2::element_line(color = "#F0F0F0", linewidth = 0.25),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
  )

combined_png <- file.path(out_fig_dir, "Supplementary_fastp_trimmed_current78_Hallmark_GSEA_running_curves_mouse_style.png")
combined_pdf <- file.path(out_fig_dir, "Supplementary_fastp_trimmed_current78_Hallmark_GSEA_running_curves_mouse_style.pdf")

# Generate exactly one PNG and one PDF.
ggplot2::ggsave(combined_png, p_combined, width = 11.8, height = 8.0, dpi = 320)
ggplot2::ggsave(combined_pdf, p_combined, width = 11.8, height = 8.0)

cat("Saved plot: ", combined_png, "\n", sep = "")
cat("Saved plot: ", combined_pdf, "\n\n", sep = "")

## =========================
## 6. Save source data, stats, summary
## =========================

running_source_file <- file.path(out_table_dir, "stepH2w12_current78_trimmed_GSEA_running_curve_source_data.csv")
target_results_file_out <- file.path(out_table_dir, "stepH2w12_current78_trimmed_target_pathway_results_used_for_plot.csv")
stability_context_file <- file.path(out_table_dir, "stepH2w12_current78_GSEA_stability_summary_context.csv")

write_csv(running_df, running_source_file)
write_csv(target_results, target_results_file_out)
if (nrow(stability_summary) > 0) write_csv(stability_summary, stability_context_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "ggplot2_version",
    "dplyr_version",
    "source_step",
    "branch_plotted",
    "rank_statistic",
    "gene_identifier_for_GSEA",
    "target_pathways_note",
    "figure_role"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("dplyr")),
    stepH2w9_dir,
    "trimmed fastp sensitivity branch",
    "sign(logFC) * -log10(PValue), inherited from StepH2w9_current78",
    "pig gene symbols",
    "fgsea was run on the full Hallmark collection in StepH2w9; target pathways were extracted afterward",
    "supplementary sensitivity figure, not replacement for primary Figure4B"
  ),
  stringsAsFactors = FALSE
)

version_file <- file.path(out_table_dir, "stepH2w12_current78_versions_and_method_records.csv")
write_csv(version_df, version_file)

summary_df <- data.frame(
  metric = c(
    "stepH2w9_input_dir",
    "ranked_symbols_trimmed_t7",
    "ranked_symbols_trimmed_t28",
    "target_pathways_plotted",
    "n_target_panels",
    "trimmed_inflammatory_t7_NES",
    "trimmed_inflammatory_t7_FDR",
    "trimmed_inflammatory_t28_NES",
    "trimmed_inflammatory_t28_FDR",
    "trimmed_EMT_t7_NES",
    "trimmed_EMT_t7_FDR",
    "trimmed_EMT_t28_NES",
    "trimmed_EMT_t28_FDR",
    "plot_png",
    "plot_pdf",
    "running_source_data",
    "target_results_source_data",
    "output_fig_dir",
    "output_table_dir"
  ),
  value = c(
    stepH2w9_dir,
    length(stats_t7),
    length(stats_t28),
    paste(target_pathways, collapse = "; "),
    nlevels(running_df$panel_label),
    target_results$NES[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & target_results$contrast == "ACLT_untreated_t7_vs_Control"][1],
    target_results$padj[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & target_results$contrast == "ACLT_untreated_t7_vs_Control"][1],
    target_results$NES[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & target_results$contrast == "ACLT_untreated_t28_vs_Control"][1],
    target_results$padj[target_results$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & target_results$contrast == "ACLT_untreated_t28_vs_Control"][1],
    target_results$NES[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & target_results$contrast == "ACLT_untreated_t7_vs_Control"][1],
    target_results$padj[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & target_results$contrast == "ACLT_untreated_t7_vs_Control"][1],
    target_results$NES[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & target_results$contrast == "ACLT_untreated_t28_vs_Control"][1],
    target_results$padj[target_results$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & target_results$contrast == "ACLT_untreated_t28_vs_Control"][1],
    combined_png,
    combined_pdf,
    running_source_file,
    target_results_file_out,
    out_fig_dir,
    out_table_dir
  ),
  stringsAsFactors = FALSE
)

summary_file <- file.path(out_table_dir, "stepH2w12_current78_trimmed_GSEA_running_curves_summary.csv")
write_csv(summary_df, summary_file)

saveRDS(
  list(
    target_results = target_results,
    running_df = running_df,
    label_df = label_df,
    stability_summary = stability_summary,
    summary_df = summary_df,
    version_df = version_df
  ),
  file.path(objects_dir, "stepH2w12_current78_trimmed_GSEA_running_curves_workspace.rds")
)

writeLines(capture.output(sessionInfo()), file.path(out_table_dir, "stepH2w12_current78_sessionInfo.txt"))

cat("\n===== StepH2w12_current78 trimmed GSEA figure summary =====\n")
print(summary_df, row.names = FALSE)

cat("\nTarget pathway statistics used for plot:\n")
print(target_results[, c("branch", "contrast", "pathway", "NES", "padj", "pval")], row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

cat("\nStepH2w12_current78 trimmed GSEA running-curve figure completed successfully.\n")

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== StepH2w12_current78 trimmed GSEA running-curve figure SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(summary_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Target pathway statistics used for plot:", summary_con)
writeLines(capture.output(print(target_results[, c("branch", "contrast", "pathway", "NES", "padj", "pval")], row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This figure uses trimmed-branch target pathway statistics and ranked lists from StepH2w9_current78.", summary_con)
writeLines("It follows the same Step20 Figure4B green running-curve style.", summary_con)
writeLines("The fastp branch is a sensitivity analysis, not a replacement for the primary untrimmed pig early GSEA.", summary_con)
close(summary_con)

sink()

cat("\nStepH2w12_current78 trimmed GSEA running-curve figure completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
