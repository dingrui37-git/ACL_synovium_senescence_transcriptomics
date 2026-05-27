# StepH2w11_current78: fastp-trimmed branch signature score figure
# Purpose:
# Generate a Supplementary Figure-style plot for the fastp-trimmed current78-derived
# pig early signature scores, using the same visual style as the main Step19 Figure4A
# signature score plot.
#
# Data source:
#   StepH2w8_current78_signature_score_stability outputs
#
# Figure:
#   - Three facets: Directional score, Up-score, Down-score reoriented
#   - Groups: Control, ACLT_t7, ACLT_t28
#   - P values are taken from the trimmed-branch Wilcoxon group comparisons
#   - Only ONE PNG and ONE PDF are generated
#
# This script does NOT auto-archive or overwrite any manually saved R script.

options(stringsAsFactors = FALSE)

## =========================
## 0. Basic settings
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

fastp_dir <- file.path(rebuild_root, "02_pig_early_fastp_sensitivity")
stepH2w8_dir <- file.path(fastp_dir, "tables", "stepH2w8_current78_signature_score_stability")

out_dir <- file.path(fastp_dir, "figures", "stepH2w11_current78_trimmed_signature_score_figure")
table_out_dir <- file.path(fastp_dir, "tables", "stepH2w11_current78_trimmed_signature_score_figure")
log_dir <- file.path(table_out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "stepH2w11_current78_trimmed_signature_score_figure_full_log.txt")
summary_log_file <- file.path(log_dir, "stepH2w11_current78_trimmed_signature_score_figure_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== StepH2w11_current78 trimmed signature score figure =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Purpose: generate trimmed-branch current78 pig early signature score plot in main Figure4A style.\n\n")

## =========================
## 1. Helper functions
## =========================

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

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "p = NA",
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))
  )
}

## =========================
## 2. Input files
## =========================

scores_long_file <- file.path(stepH2w8_dir, "stepH2w8_current78_signature_scores_long.csv")
group_comparison_file <- file.path(stepH2w8_dir, "stepH2w8_current78_signature_score_group_comparison_by_branch.csv")
stability_metrics_file <- file.path(stepH2w8_dir, "stepH2w8_current78_signature_score_stability_metrics.csv")

stop_if_missing(scores_long_file, "StepH2w8 score long table")
stop_if_missing(group_comparison_file, "StepH2w8 group comparison table")
stop_if_missing(stability_metrics_file, "StepH2w8 stability metrics table")

cat("Scores long file: ", scores_long_file, "\n", sep = "")
cat("Group comparison file: ", group_comparison_file, "\n", sep = "")
cat("Stability metrics file: ", stability_metrics_file, "\n\n", sep = "")

scores_long_all <- safe_read_csv(scores_long_file)
group_stats_all <- safe_read_csv(group_comparison_file)
stability_metrics <- safe_read_csv(stability_metrics_file)

required_score_cols <- c("branch", "sample_id", "group", "score_name", "score_value")
if (!all(required_score_cols %in% colnames(scores_long_all))) {
  stop("Score long table lacks required columns: ", paste(setdiff(required_score_cols, colnames(scores_long_all)), collapse = ", "), call. = FALSE)
}
required_stat_cols <- c("branch", "comparison", "score", "wilcox_p_value")
if (!all(required_stat_cols %in% colnames(group_stats_all))) {
  stop("Group comparison table lacks required columns: ", paste(setdiff(required_stat_cols, colnames(group_stats_all)), collapse = ", "), call. = FALSE)
}

## =========================
## 3. Prepare trimmed branch plot data
## =========================

score_plot_levels <- c("directional_score", "up_score", "down_score_reoriented")
score_plot_labels <- c(
  directional_score = "Directional score",
  up_score = "Up-score",
  down_score_reoriented = "Down-score reoriented"
)

plot_long <- scores_long_all[
  scores_long_all$branch == "trimmed" & scores_long_all$score_name %in% score_plot_levels,
  ,
  drop = FALSE
]

if (nrow(plot_long) == 0) stop("No trimmed-branch score rows found for manuscript-relevant scores.", call. = FALSE)

plot_long$score_value <- as.numeric(plot_long$score_value)
plot_long$score_label <- factor(score_plot_labels[plot_long$score_name],
                                levels = unname(score_plot_labels[score_plot_levels]))
plot_long$group <- factor(plot_long$group, levels = c("Control", "ACLT_t7", "ACLT_t28"))

stats_plot <- group_stats_all[
  group_stats_all$branch == "trimmed" &
    group_stats_all$score %in% score_plot_levels &
    group_stats_all$comparison %in% c("ACLT_t7_vs_Control", "ACLT_t28_vs_Control"),
  ,
  drop = FALSE
]

if (nrow(stats_plot) != 6) {
  stop("Expected 6 trimmed-branch comparison rows for 3 scores x 2 comparisons, but found: ", nrow(stats_plot), call. = FALSE)
}

stats_plot$score_label <- factor(score_plot_labels[stats_plot$score],
                                 levels = unname(score_plot_labels[score_plot_levels]))
stats_plot$xend <- ifelse(stats_plot$comparison == "ACLT_t7_vs_Control", 2, 3)
stats_plot$xmid <- stats_plot$xend
stats_plot$group_label <- ifelse(stats_plot$comparison == "ACLT_t7_vs_Control", "ACLT_t7", "ACLT_t28")
stats_plot$p_label <- fmt_p(as.numeric(stats_plot$wilcox_p_value))

# Put p labels slightly above each score facet, matching Step19 style.
y_ranges <- aggregate(score_value ~ score_label, plot_long, function(x) diff(range(x, na.rm = TRUE)))
y_max <- aggregate(score_value ~ score_label, plot_long, max, na.rm = TRUE)
y_pos <- merge(y_max, y_ranges, by = "score_label")
colnames(y_pos) <- c("score_label", "y_max", "y_range")
y_pos$y_range[!is.finite(y_pos$y_range) | y_pos$y_range == 0] <- 1

stats_plot <- merge(stats_plot, y_pos, by = "score_label", all.x = TRUE)
stats_plot$y_label <- stats_plot$y_max + ifelse(stats_plot$group_label == "ACLT_t7", 0.10, 0.22) * stats_plot$y_range
stats_plot$y_segment <- stats_plot$y_max + ifelse(stats_plot$group_label == "ACLT_t7", 0.05, 0.17) * stats_plot$y_range

## =========================
## 4. Plot
## =========================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
suppressPackageStartupMessages(library(ggplot2))

point_fill <- c(Control = "white", ACLT_t7 = "#D95F02", ACLT_t28 = "#E69F00")
point_edge <- c(Control = "#7A7A7A", ACLT_t7 = "#7F3300", ACLT_t28 = "#8A5A00")

p_combined <- ggplot(plot_long, aes(x = group, y = score_value)) +
  geom_point(
    aes(fill = group, color = group),
    shape = 21,
    stroke = 0.55,
    size = 3.1,
    position = position_jitter(width = 0.08, height = 0, seed = 1),
    alpha = 0.95
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.44,
    linewidth = 0.35,
    color = "black"
  ) +
  geom_segment(
    data = stats_plot,
    aes(x = 1, xend = xend, y = y_segment, yend = y_segment),
    inherit.aes = FALSE,
    color = "#666666",
    linewidth = 0.35
  ) +
  geom_text(
    data = stats_plot,
    aes(x = xmid, y = y_label, label = p_label),
    inherit.aes = FALSE,
    size = 4.0
  ) +
  facet_grid(rows = vars(score_label), switch = "y", scales = "free_y") +
  scale_fill_manual(values = point_fill) +
  scale_color_manual(values = point_edge) +
  labs(
    title = "Fastp-trimmed current78-derived pig early signature scores",
    x = NULL,
    y = "Mean row-z signature score"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 19),
    panel.grid.major = element_line(color = "#E6E6E6"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = "black"),
    strip.text = element_text(size = 13),
    strip.text.y.left = element_text(angle = -90, size = 13),
    axis.text.x = element_text(size = 12, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15),
    legend.position = "none"
  )

combined_png <- file.path(out_dir, "Supplementary_fastp_trimmed_current78_signature_scores_mouse_style.png")
combined_pdf <- file.path(out_dir, "Supplementary_fastp_trimmed_current78_signature_scores_mouse_style.pdf")

# Generate exactly one PNG and one PDF.
ggsave(combined_png, p_combined, width = 7.4, height = 9.2, dpi = 320)
ggsave(combined_pdf, p_combined, width = 7.4, height = 9.2)

cat("Saved plot: ", combined_png, "\n", sep = "")
cat("Saved plot: ", combined_pdf, "\n\n", sep = "")

## =========================
## 5. Save source data and summary
## =========================

source_data_file <- file.path(table_out_dir, "stepH2w11_current78_trimmed_signature_score_plot_source_data.csv")
stats_source_file <- file.path(table_out_dir, "stepH2w11_current78_trimmed_signature_score_plot_statistics.csv")
stability_source_file <- file.path(table_out_dir, "stepH2w11_current78_signature_score_stability_metrics_used_for_context.csv")

write_csv(plot_long, source_data_file)
write_csv(stats_plot, stats_source_file)
write_csv(stability_metrics, stability_source_file)

summary_df <- data.frame(
  metric = c(
    "scores_long_input_file",
    "group_comparison_input_file",
    "branch_plotted",
    "scores_plotted",
    "n_samples_plotted",
    "n_Control",
    "n_ACLT_t7",
    "n_ACLT_t28",
    "directional_score_trimmed_t7_p",
    "directional_score_trimmed_t28_p",
    "up_score_trimmed_t7_p",
    "up_score_trimmed_t28_p",
    "down_score_reoriented_trimmed_t7_p",
    "down_score_reoriented_trimmed_t28_p",
    "plot_png",
    "plot_pdf",
    "source_data_file",
    "statistics_file",
    "output_dir"
  ),
  value = c(
    scores_long_file,
    group_comparison_file,
    "trimmed",
    paste(score_plot_levels, collapse = "; "),
    length(unique(plot_long$sample_id)),
    length(unique(plot_long$sample_id[plot_long$group == "Control"])),
    length(unique(plot_long$sample_id[plot_long$group == "ACLT_t7"])),
    length(unique(plot_long$sample_id[plot_long$group == "ACLT_t28"])),
    group_stats_all$wilcox_p_value[group_stats_all$branch == "trimmed" & group_stats_all$score == "directional_score" & group_stats_all$comparison == "ACLT_t7_vs_Control"][1],
    group_stats_all$wilcox_p_value[group_stats_all$branch == "trimmed" & group_stats_all$score == "directional_score" & group_stats_all$comparison == "ACLT_t28_vs_Control"][1],
    group_stats_all$wilcox_p_value[group_stats_all$branch == "trimmed" & group_stats_all$score == "up_score" & group_stats_all$comparison == "ACLT_t7_vs_Control"][1],
    group_stats_all$wilcox_p_value[group_stats_all$branch == "trimmed" & group_stats_all$score == "up_score" & group_stats_all$comparison == "ACLT_t28_vs_Control"][1],
    group_stats_all$wilcox_p_value[group_stats_all$branch == "trimmed" & group_stats_all$score == "down_score_reoriented" & group_stats_all$comparison == "ACLT_t7_vs_Control"][1],
    group_stats_all$wilcox_p_value[group_stats_all$branch == "trimmed" & group_stats_all$score == "down_score_reoriented" & group_stats_all$comparison == "ACLT_t28_vs_Control"][1],
    combined_png,
    combined_pdf,
    source_data_file,
    stats_source_file,
    out_dir
  ),
  stringsAsFactors = FALSE
)

summary_file <- file.path(table_out_dir, "stepH2w11_current78_trimmed_signature_score_figure_summary.csv")
write_csv(summary_df, summary_file)

cat("\n===== StepH2w11_current78 trimmed signature score figure summary =====\n")
print(summary_df, row.names = FALSE)

cat("\nTrimmed-branch statistics used for plot:\n")
print(stats_plot[, c("score", "comparison", "wilcox_p_value", "p_label")], row.names = FALSE)

cat("\nStepH2w11_current78 trimmed signature score figure completed successfully.\n")

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== StepH2w11_current78 trimmed signature score figure SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(summary_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Trimmed-branch statistics used for plot:", summary_con)
writeLines(capture.output(print(stats_plot[, c("score", "comparison", "wilcox_p_value", "p_label")], row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This figure uses trimmed-branch scores from StepH2w8_current78 and follows the main Step19 Figure4A visual style.", summary_con)
writeLines("The fastp branch is a sensitivity analysis, not a replacement for the primary untrimmed pig early analysis.", summary_con)
close(summary_con)

sink()

cat("\nStepH2w11_current78 trimmed signature score figure completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
