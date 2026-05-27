# Step21_current78: Figure4C validation flow barplot
# Purpose:
# Create a Figure4C-style horizontal barplot that summarizes the current78 pig validation
# flow from mouse-derived persistent CellAge-overlap genes to the final 24 core ortholog genes.
# This script follows the current formal workflow outputs and generates both a full log and
# a concise summary log for review.

options(stringsAsFactors = FALSE)

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
}
install_if_missing(c("readr", "dplyr", "ggplot2"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

find_col <- function(df, candidates, label) {
  idx <- match(candidates, colnames(df))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    stop(sprintf("Could not find column for %s. Checked: %s", label, paste(candidates, collapse = ", ")))
  }
  colnames(df)[idx[1]]
}

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early"
step18_dir <- file.path(base_dir, "tables", "step18_current78_pig_early_signature_remap")
out_dir <- file.path(base_dir, "tables", "step21_current78_Figure4C_validation_flow_barplot")
fig_dir <- file.path(base_dir, "figures", "Figure4")
log_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "step21_current78_Figure4C_validation_flow_barplot_full_log.txt")
summary_log_file <- file.path(log_dir, "step21_current78_Figure4C_validation_flow_barplot_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)
cat("===== STEP21 CURRENT78 FIGURE4C VALIDATION FLOW BARPLOT =====\n")

summary_file <- file.path(step18_dir, "step18_current78_pig_early_signature_remap_summary.csv")
gorth_raw_file <- file.path(step18_dir, "step18_current78_human_to_pig_gorth_raw.csv")
validation_file <- file.path(step18_dir, "step18_current78_pig_signature_validation_table.csv")

if (!file.exists(summary_file)) stop("Missing Step18 current78 summary file: ", summary_file)
if (!file.exists(gorth_raw_file)) stop("Missing Step18 current78 gorth raw file: ", gorth_raw_file)
if (!file.exists(validation_file)) stop("Missing Step18 current78 validation file: ", validation_file)

cat("Summary file: ", summary_file, "\n", sep = "")
cat("gorth raw file: ", gorth_raw_file, "\n", sep = "")
cat("Validation file: ", validation_file, "\n", sep = "")

summary_df_raw <- read_csv(summary_file, show_col_types = FALSE)
gorth_raw <- read_csv(gorth_raw_file, show_col_types = FALSE)
validation_df <- read_csv(validation_file, show_col_types = FALSE)

get_metric <- function(metric_name) {
  x <- summary_df_raw$value[match(metric_name, summary_df_raw$metric)]
  if (length(x) == 0 || is.na(x)) stop("Missing metric in Step18 summary: ", metric_name)
  as.numeric(x)
}

# Raw mapping stage: count unique human inputs that produced any raw gorth row
raw_input_col <- find_col(gorth_raw, c("input", "human_gene", "query", "Gene symbol"), "gorth raw input gene")
raw_mapped_unique_human <- gorth_raw %>%
  mutate(raw_input = as.character(.data[[raw_input_col]])) %>%
  filter(!is.na(raw_input), raw_input != "") %>%
  distinct(raw_input) %>%
  nrow()

n_input_human <- get_metric("current_human_CellAge_overlap_input_genes")
n_strict_1to1 <- get_metric("current_human_genes_strict_1to1_mapped_to_pig")
n_detected <- get_metric("current_pig_signature_genes_detected")
n_direction_consistent <- get_metric("direction_consistent_both_timepoints")
n_core <- get_metric("core_strict_both_t7_t28")

plot_df <- data.frame(
  step_order = 1:6,
  category = c(
    sprintf("%d human genes\n(mouse persistent ∩ CellAge)", n_input_human),
    sprintf("%d human→pig raw ortholog mapping", raw_mapped_unique_human),
    sprintf("%d strict 1:1 pig orthologs", n_strict_1to1),
    sprintf("%d detected in pig expression matrix", n_detected),
    sprintf("%d direction-consistent in pig\n(t7 and t28 direction match mouse)", n_direction_consistent),
    sprintf("%d strict core ortholog genes\n(t7 + t28 + direction-consistent)", n_core)
  ),
  n = c(n_input_human, raw_mapped_unique_human, n_strict_1to1, n_detected, n_direction_consistent, n_core),
  stringsAsFactors = FALSE
)

# top-to-bottom order like the reference figure
plot_df$category <- factor(plot_df$category, levels = rev(plot_df$category))
plot_df <- plot_df %>% arrange(desc(step_order))

# Save source data
write_csv(plot_df, file.path(out_dir, "step21_current78_Figure4C_validation_flow_source_data.csv"))

# Plot
max_n <- max(plot_df$n)
text_nudge <- max_n * 0.005
x_upper <- max_n + max(5, ceiling(max_n * 0.12))

p <- ggplot(plot_df, aes(x = n, y = category)) +
  geom_col(fill = "#5F5F5F", width = 0.75) +
  geom_text(aes(label = n), hjust = -0.05, size = 6) +
  scale_x_continuous(limits = c(0, x_upper), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Pig validation flow from mouse-derived persistent CellAge genes",
    x = "Number of genes",
    y = NULL
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#DDDDDD", linewidth = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 15),
    axis.line = element_line(color = "black"),
    plot.margin = margin(15, 20, 15, 15)
  )

png_file <- file.path(fig_dir, "Figure4C_current78_pig_validation_flow_barplot.png")
pdf_file <- file.path(fig_dir, "Figure4C_current78_pig_validation_flow_barplot.pdf")
canonical_png <- file.path(fig_dir, "Figure4C_pig_validation_flow_barplot.png")
canonical_pdf <- file.path(fig_dir, "Figure4C_pig_validation_flow_barplot.pdf")

png(png_file, width = 1800, height = 1300, res = 200)
print(p)
dev.off()

pdf(pdf_file, width = 12, height = 8.5)
print(p)
dev.off()

file.copy(png_file, canonical_png, overwrite = TRUE)
file.copy(pdf_file, canonical_pdf, overwrite = TRUE)

summary_df <- data.frame(
  metric = c(
    "current_human_CellAge_overlap_input_genes",
    "human_genes_with_raw_pig_mapping",
    "current_human_genes_strict_1to1_mapped_to_pig",
    "current_pig_signature_genes_detected",
    "direction_consistent_both_timepoints",
    "core_strict_both_t7_t28",
    "plot_png",
    "plot_pdf",
    "canonical_plot_png",
    "canonical_plot_pdf",
    "output_dir"
  ),
  value = c(
    n_input_human,
    raw_mapped_unique_human,
    n_strict_1to1,
    n_detected,
    n_direction_consistent,
    n_core,
    png_file,
    pdf_file,
    canonical_png,
    canonical_pdf,
    out_dir
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_df, file.path(out_dir, "step21_current78_Figure4C_validation_flow_barplot_summary.csv"))

cat("\n===== STEP21 CURRENT78 FIGURE4C VALIDATION FLOW SUMMARY =====\n")
print(summary_df)
cat("\nPlot source data:\n")
print(plot_df)
cat("\nStep21 current78 Figure4C validation flow barplot completed successfully.\n")

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== STEP21 CURRENT78 FIGURE4C VALIDATION FLOW SUMMARY TO SEND ME =====", summary_con)
writeLines(capture.output(print(summary_df)), summary_con)
writeLines("", summary_con)
writeLines("Plot source data:", summary_con)
writeLines(capture.output(print(plot_df)), summary_con)
close(summary_con)

sink()
cat("\nStep21 current78 Figure4C validation flow barplot completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
