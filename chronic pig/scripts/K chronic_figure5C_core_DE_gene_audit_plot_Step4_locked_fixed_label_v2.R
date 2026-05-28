# Chronic Figure5C current78: gene-level lollipop/dot plot for core-gene chronic DE audit
# Purpose:
#   Generate Figure5C as a gene-level chronic DE audit plot for the fixed early-defined
#   24 core ortholog genes, using a lollipop/dot-plot design.
#
# Main display:
#   - y-axis: gene names
#   - x-axis: chronic logFC (ACLT_alone_52W vs Control_52W)
#   - color: chronic retention class
#   - shape: early signature direction (Up_in_ACLR / Down_in_ACLR)
#   - larger point / label mark: chronic strict + direction-consistent genes
#   - genes not tested in chronic DE are displayed at a dedicated "Not tested" position
#
# Interpretation:
#   - This figure does NOT redefine a chronic core set.
#   - It audits how the fixed early-defined 24 core genes behave in chronic DE.
#
# Output:
#   One PNG and one PDF only, plus source data/statistics/summary logs.
#   This script does NOT auto-archive or overwrite any manually saved R script.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_chronic_dir <- file.path(rebuild_root, "03_pig_chronic")
step5b_dir <- file.path(pig_chronic_dir, "tables", "chronic_step5B_current78_core_DE_audit_Step4_locked")

out_table_dir <- file.path(pig_chronic_dir, "tables", "chronic_figure5C_current78_core_DE_gene_audit_plot_Step4_locked_fixed_label_v2")
out_fig_dir <- file.path(pig_chronic_dir, "figures", "Figure5C_current78_core_DE_gene_audit_plot_Step4_locked_fixed_label_v2")
log_dir <- file.path(out_table_dir, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "chronic_figure5C_current78_core_DE_gene_audit_plot_Step4_locked_fixed_label_v2_full_log.txt")
summary_log_file <- file.path(log_dir, "chronic_figure5C_current78_core_DE_gene_audit_plot_Step4_locked_fixed_label_v2_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== CHRONIC FIGURE5C CURRENT78 CORE DE GENE-LEVEL AUDIT PLOT =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Purpose: generate Figure5C as a gene-level lollipop/dot plot for the chronic DE audit of the fixed early-defined 24 core genes.\n\n")

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

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

retention_label_map <- c(
  direction_consistent_strict = "Strict + direction-consistent",
  direction_consistent_FDR_only = "Direction-consistent, FDR only",
  direction_consistent_large_effect_only = "Direction-consistent, |logFC| > 1 only",
  direction_consistent_but_attenuated = "Direction-consistent but attenuated",
  opposite_or_reversed = "Opposite / reversed",
  not_tested = "Not tested in chronic DE"
)

retention_rank_map <- c(
  direction_consistent_strict = 1,
  direction_consistent_FDR_only = 2,
  direction_consistent_large_effect_only = 3,
  direction_consistent_but_attenuated = 4,
  opposite_or_reversed = 5,
  not_tested = 6
)

## =========================
## 2. Input files
## =========================

audit_table_file <- file.path(step5b_dir, "chronic_step5B_early_defined_24_core_chronic_DE_audit_table.csv")
summary_file <- file.path(step5b_dir, "chronic_step5B_early_defined_24_core_chronic_DE_audit_summary.csv")
strict_retained_file <- file.path(step5b_dir, "chronic_step5B_core_genes_chronic_strict_and_direction_consistent.csv")

stop_if_missing(audit_table_file, "Chronic Step5B audit table")
stop_if_missing(summary_file, "Chronic Step5B summary")
stop_if_missing(strict_retained_file, "Chronic Step5B strict retained gene table")

audit_df <- safe_read_csv(audit_table_file)
summary_df <- safe_read_csv(summary_file)
strict_df <- safe_read_csv(strict_retained_file)

required_cols <- c(
  "display_symbol", "pig_ensg", "early_signature_direction",
  "tested_in_chronic_DE", "chronic_logFC", "chronic_PValue", "chronic_FDR",
  "chronic_strict_DEG", "chronic_strict_and_direction_consistent",
  "chronic_retention_class"
)
missing_cols <- setdiff(required_cols, colnames(audit_df))
if (length(missing_cols) > 0) {
  stop("Audit table missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

cat("Input audit table: ", audit_table_file, "\n", sep = "")
cat("Input summary file: ", summary_file, "\n\n", sep = "")

## =========================
## 3. Prepare plot data
## =========================

plot_df <- audit_df

plot_df$display_symbol <- as.character(plot_df$display_symbol)
plot_df$display_symbol[is.na(plot_df$display_symbol) | plot_df$display_symbol == ""] <- plot_df$pig_ensg[is.na(plot_df$display_symbol) | plot_df$display_symbol == ""]
plot_df$early_signature_direction <- factor(plot_df$early_signature_direction, levels = c("Up_in_ACLR", "Down_in_ACLR"))
plot_df$tested_in_chronic_DE <- as.logical(plot_df$tested_in_chronic_DE)
plot_df$chronic_logFC <- suppressWarnings(as.numeric(plot_df$chronic_logFC))
plot_df$chronic_PValue <- suppressWarnings(as.numeric(plot_df$chronic_PValue))
plot_df$chronic_FDR <- suppressWarnings(as.numeric(plot_df$chronic_FDR))
plot_df$chronic_strict_DEG <- as.logical(plot_df$chronic_strict_DEG)
plot_df$chronic_strict_and_direction_consistent <- as.logical(plot_df$chronic_strict_and_direction_consistent)
plot_df$chronic_retention_class <- as.character(plot_df$chronic_retention_class)

# Normalize class names if needed
plot_df$chronic_retention_class[is.na(plot_df$chronic_retention_class) & !plot_df$tested_in_chronic_DE] <- "not_tested"
plot_df$retention_label <- retention_label_map[plot_df$chronic_retention_class]
plot_df$retention_rank <- retention_rank_map[plot_df$chronic_retention_class]
plot_df$retention_label[is.na(plot_df$retention_label)] <- plot_df$chronic_retention_class[is.na(plot_df$retention_label)]
plot_df$retention_rank[is.na(plot_df$retention_rank)] <- 99

tested_vals <- plot_df$chronic_logFC[plot_df$tested_in_chronic_DE & is.finite(plot_df$chronic_logFC)]
if (length(tested_vals) == 0) {
  stop("No tested chronic logFC values found in audit table.", call. = FALSE)
}

x_min <- min(tested_vals, na.rm = TRUE)
x_max <- max(tested_vals, na.rm = TRUE)
x_range <- x_max - x_min
if (!is.finite(x_range) || x_range == 0) x_range <- 1

# Dedicated plotting location for not-tested genes
not_test_x <- x_min - 0.85 * x_range
plot_df$plot_x <- ifelse(plot_df$tested_in_chronic_DE, plot_df$chronic_logFC, not_test_x)
plot_df$plot_x_label <- ifelse(plot_df$tested_in_chronic_DE, fmt_num(plot_df$chronic_logFC, 2), "Not tested")

# Text-label positions.
# - Not-tested labels are placed to the right of the not-tested point to avoid left clipping.
# - For strict + direction-consistent genes, the asterisk is placed between the point and
#   the logFC number, and the numeric label is moved farther away to avoid overlap.
plot_df$is_strict_dc <- plot_df$chronic_strict_and_direction_consistent %in% TRUE
plot_df$is_positive_tested <- plot_df$tested_in_chronic_DE & is.finite(plot_df$plot_x) & plot_df$plot_x >= 0
plot_df$is_negative_tested <- plot_df$tested_in_chronic_DE & is.finite(plot_df$plot_x) & plot_df$plot_x < 0

plot_df$star_x <- ifelse(
  plot_df$is_strict_dc & plot_df$is_positive_tested,
  plot_df$plot_x + 0.095 * x_range,
  ifelse(
    plot_df$is_strict_dc & plot_df$is_negative_tested,
    plot_df$plot_x + 0.095 * x_range,
    NA_real_
  )
)

plot_df$label_x <- ifelse(
  plot_df$tested_in_chronic_DE,
  ifelse(
    plot_df$is_positive_tested,
    plot_df$plot_x + ifelse(plot_df$is_strict_dc, 0.235 * x_range, 0.065 * x_range),
    plot_df$plot_x - 0.065 * x_range
  ),
  not_test_x + 0.12 * x_range
)
plot_df$label_hjust <- ifelse(
  plot_df$tested_in_chronic_DE,
  ifelse(plot_df$is_positive_tested, 0, 1),
  0
)

# Gene ordering:
#   1. retention class priority
#   2. Up_in_ACLR before Down_in_ACLR
#   3. within class, by chronic logFC magnitude
direction_rank <- ifelse(plot_df$early_signature_direction == "Up_in_ACLR", 1, 2)
magnitude_rank <- ifelse(is.na(plot_df$chronic_logFC), -Inf, abs(plot_df$chronic_logFC))

ord <- order(
  plot_df$retention_rank,
  direction_rank,
  -magnitude_rank,
  plot_df$display_symbol
)

plot_df <- plot_df[ord, , drop = FALSE]
plot_df$gene_label <- factor(plot_df$display_symbol, levels = rev(plot_df$display_symbol))

# Add "*" mark to strict retained genes on y-axis display if desired in a companion label column
plot_df$strict_star <- ifelse(plot_df$chronic_strict_and_direction_consistent %in% TRUE, "*", "")

# Colors
retention_colors <- c(
  "Strict + direction-consistent" = "#D95F02",
  "Direction-consistent, FDR only" = "#FDB863",
  "Direction-consistent, |logFC| > 1 only" = "#80CDC1",
  "Direction-consistent but attenuated" = "#67A9CF",
  "Opposite / reversed" = "#8073AC",
  "Not tested in chronic DE" = "#BDBDBD"
)

# Statistics summary for log
n_fixed_core <- nrow(plot_df)
n_tested <- sum(plot_df$tested_in_chronic_DE, na.rm = TRUE)
n_not_tested <- sum(!plot_df$tested_in_chronic_DE, na.rm = TRUE)
n_strict_dc <- sum(plot_df$chronic_strict_and_direction_consistent, na.rm = TRUE)
n_reversed <- sum(plot_df$chronic_retention_class == "opposite_or_reversed", na.rm = TRUE)

strict_symbols <- plot_df$display_symbol[plot_df$chronic_strict_and_direction_consistent %in% TRUE]
not_tested_symbols <- plot_df$display_symbol[!plot_df$tested_in_chronic_DE]
reversed_symbols <- plot_df$display_symbol[plot_df$chronic_retention_class == "opposite_or_reversed"]

## =========================
## 4. Plot
## =========================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
suppressPackageStartupMessages(library(ggplot2))

# Segments only for tested genes
segment_df <- plot_df[plot_df$tested_in_chronic_DE, , drop = FALSE]

# Dynamic figure height for gene labels
plot_height <- max(8.2, 0.36 * nrow(plot_df) + 2.7)

# Axis breaks: include special "Not tested" position
pretty_breaks <- pretty(c(x_min, x_max), n = 6)
x_breaks <- c(not_test_x, pretty_breaks)
x_breaks <- unique(round(x_breaks, 6))
x_labels <- ifelse(abs(x_breaks - not_test_x) < 1e-8, "Not tested", as.character(x_breaks))

# small buffer for text labels
right_buffer <- x_max + 0.82 * x_range
left_buffer <- not_test_x - 0.12 * x_range

p <- ggplot(plot_df, aes(y = gene_label)) +
  geom_vline(xintercept = 0, color = "#999999", linetype = "dashed", linewidth = 0.4) +
  geom_segment(
    data = segment_df,
    aes(x = 0, xend = plot_x, y = gene_label, yend = gene_label, color = retention_label),
    linewidth = 0.65,
    alpha = 0.9
  ) +
  geom_point(
    aes(x = plot_x, fill = retention_label, shape = early_signature_direction),
    color = "black",
    stroke = 0.45,
    size = ifelse(plot_df$chronic_strict_and_direction_consistent %in% TRUE, 3.9, 3.0)
  ) +
  geom_text(
    data = plot_df[plot_df$chronic_strict_and_direction_consistent %in% TRUE, , drop = FALSE],
    aes(x = star_x, label = "*"),
    hjust = 0.5,
    size = 5.0,
    fontface = "bold"
  ) +
  geom_text(
    aes(x = label_x, label = plot_x_label, hjust = label_hjust),
    size = 3.05,
    color = "#333333"
  ) +
  scale_fill_manual(values = retention_colors, drop = FALSE) +
  scale_color_manual(values = retention_colors, drop = FALSE) +
  scale_shape_manual(
    values = c("Up_in_ACLR" = 24, "Down_in_ACLR" = 25),
    labels = c("Up_in_ACLR" = "Early up in ACLR", "Down_in_ACLR" = "Early down in ACLR"),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    limits = c(left_buffer, right_buffer),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Chronic DE audit of the fixed early-defined 24 core genes",
    subtitle = "Asterisks mark genes that remain strict and direction-consistent in chronic DE",
    x = "Chronic logFC (ACLT_alone_52W vs Control_52W)",
    y = NULL,
    fill = "Retention class",
    color = "Retention class",
    shape = "Early signature direction"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#E6E6E6"),
    axis.text.y = element_text(size = 10.5),
    axis.text.x = element_text(size = 10.5),
    axis.title.x = element_text(size = 13),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10.5),
    legend.text = element_text(size = 9.5),
    plot.margin = margin(t = 8, r = 24, b = 8, l = 18)
  ) +
  coord_cartesian(clip = "off") +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 4)),
    color = "none",
    shape = guide_legend(override.aes = list(fill = "white", size = 4))
  )

figure_png <- file.path(out_fig_dir, "Figure5C_current78_core_DE_gene_audit_lollipop.png")
figure_pdf <- file.path(out_fig_dir, "Figure5C_current78_core_DE_gene_audit_lollipop.pdf")

ggsave(figure_png, p, width = 11.8, height = plot_height, dpi = 320)
ggsave(figure_pdf, p, width = 11.8, height = plot_height)

cat("Saved plot: ", figure_png, "\n", sep = "")
cat("Saved plot: ", figure_pdf, "\n\n", sep = "")

## =========================
## 5. Save source data and summaries
## =========================

plot_source_file <- file.path(out_table_dir, "chronic_figure5C_current78_plot_source_data.csv")
plot_stats_file <- file.path(out_table_dir, "chronic_figure5C_current78_plot_statistics.csv")
summary_out_file <- file.path(out_table_dir, "chronic_figure5C_current78_core_DE_gene_audit_plot_summary.csv")
version_file <- file.path(out_table_dir, "chronic_figure5C_current78_versions_and_method_records.csv")
support_lists_file <- file.path(out_table_dir, "chronic_figure5C_current78_gene_category_lists.csv")

write_csv(plot_df, plot_source_file)

plot_stats <- data.frame(
  metric = c(
    "fixed_early_defined_core_genes",
    "tested_in_chronic_DE",
    "not_tested_in_chronic_DE",
    "strict_and_direction_consistent",
    "opposite_or_reversed",
    "strict_retained_genes",
    "not_tested_genes",
    "reversed_genes"
  ),
  value = c(
    n_fixed_core,
    n_tested,
    n_not_tested,
    n_strict_dc,
    n_reversed,
    ifelse(length(strict_symbols) > 0, paste(strict_symbols, collapse = ", "), "none"),
    ifelse(length(not_tested_symbols) > 0, paste(not_tested_symbols, collapse = ", "), "none"),
    ifelse(length(reversed_symbols) > 0, paste(reversed_symbols, collapse = ", "), "none")
  ),
  stringsAsFactors = FALSE
)
write_csv(plot_stats, plot_stats_file)

support_lists <- data.frame(
  category = c(
    rep("strict_and_direction_consistent", length(strict_symbols)),
    rep("not_tested_in_chronic_DE", length(not_tested_symbols)),
    rep("opposite_or_reversed", length(reversed_symbols))
  ),
  gene = c(strict_symbols, not_tested_symbols, reversed_symbols),
  stringsAsFactors = FALSE
)
if (nrow(support_lists) == 0) {
  support_lists <- data.frame(category = character(0), gene = character(0), stringsAsFactors = FALSE)
}
write_csv(support_lists, support_lists_file)

summary_out <- data.frame(
  metric = c(
    "step5B_audit_input_file",
    "n_fixed_early_defined_core_genes",
    "n_tested_in_chronic_DE",
    "n_not_tested_in_chronic_DE",
    "n_strict_and_direction_consistent",
    "n_opposite_or_reversed",
    "strict_retained_genes",
    "not_tested_genes",
    "reversed_genes",
    "figure_png",
    "figure_pdf",
    "plot_source_data_file",
    "plot_statistics_file",
    "interpretation_policy",
    "output_dir"
  ),
  value = c(
    audit_table_file,
    n_fixed_core,
    n_tested,
    n_not_tested,
    n_strict_dc,
    n_reversed,
    ifelse(length(strict_symbols) > 0, paste(strict_symbols, collapse = ", "), "none"),
    ifelse(length(not_tested_symbols) > 0, paste(not_tested_symbols, collapse = ", "), "none"),
    ifelse(length(reversed_symbols) > 0, paste(reversed_symbols, collapse = ", "), "none"),
    figure_png,
    figure_pdf,
    plot_source_file,
    plot_stats_file,
    "Figure5C is a gene-level chronic DE audit plot for the fixed early-defined 24 core genes; it does not redefine a chronic core set.",
    out_table_dir
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_out, summary_out_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "ggplot2_version",
    "figure_source_step",
    "core_definition_policy",
    "plot_design",
    "figure_role"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2")),
    "Chronic Step5B current78 core DE audit, final Step4 source locked",
    "fixed early-defined 24 core ortholog genes from pig early primary-analysis workflow",
    "gene-level lollipop/dot plot with retention-class colors and early-direction shapes",
    "main Figure5C; fixed not-tested label clipping and separated strict-star/logFC labels"
  ),
  stringsAsFactors = FALSE
)
write_csv(version_df, version_file)

cat("\n===== Chronic Figure5C gene-level audit plot summary =====\n")
print(summary_out, row.names = FALSE)

cat("\nPlot statistics:\n")
print(plot_stats, row.names = FALSE)

cat("\nPreview of plot data:\n")
print(plot_df[, c("display_symbol", "pig_ensg", "early_signature_direction", "tested_in_chronic_DE", "chronic_logFC", "chronic_FDR", "chronic_retention_class", "plot_x_label", "label_x", "star_x")], row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== Chronic Figure5C current78 core DE gene-level audit plot SUMMARY TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(summary_out, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Plot statistics:", summary_con)
writeLines(capture.output(print(plot_stats, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Preview of plot data:", summary_con)
writeLines(capture.output(print(plot_df[, c("display_symbol", "pig_ensg", "early_signature_direction", "tested_in_chronic_DE", "chronic_logFC", "chronic_FDR", "chronic_retention_class", "plot_x_label", "label_x", "star_x")], row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This figure is intended to serve as the main Figure5C.", summary_con)
writeLines("It shows which fixed early-defined core genes are strict retained, attenuated, reversed, or not tested in chronic DE.", summary_con)
writeLines("It does not redefine a chronic core gene set.", summary_con)
close(summary_con)

sink()

cat("\nChronic Figure5C current78 core DE gene-level audit plot completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
