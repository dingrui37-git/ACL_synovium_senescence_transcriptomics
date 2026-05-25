# Step08: Redraw Figure 2B logFC scatter with persistent-up, persistent-down, and opposite-direction classes.
# Purpose:
# This corrected Step08 does not rerun DE analysis. It uses the locked Step07 strict DEG
# structure table and redraws Figure 2B to show effect-size concordance between 1W and 4W.
# Shared strict genes are classified as Persistent_up, Persistent_down, or Opposite_direction.
# The plot includes x=0/y=0 dashed lines, a y=x reference line, summary text, selected gene
# labels, source data, logs, software versions, and an archived copy of this script.

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

step07_table <- file.path(
  base_dir, "07_tables", "step07_strict_DEG_upset_persistent",
  "step07_strict_DEG_structure_gene_table.csv"
)

figure_dir <- file.path(base_dir, "06_figures", "Figure2")
table_dir <- file.path(base_dir, "07_tables", "step08_logFC_scatter")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step08_logFC_scatter"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

old_files <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  file.path(figure_dir, "Figure2B_mouse_shared_strict_logFC_scatter.png"),
  file.path(figure_dir, "Figure2B_mouse_shared_strict_logFC_scatter.pdf"),
  log_file,
  file.path(script_dir, "step08_logFC_scatter.R")
)
unlink(old_files[file.exists(old_files)], force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP08 FIGURE2B LOGFC SCATTER REDRAW =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Input table: ", step07_table, "\n\n", sep = "")

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
  } else {
    writeLines("# Step08 archive fallback.", path, useBytes = TRUE)
  }
}

safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

safe_library("ggplot2")
safe_library("ggrepel")

archive_current_script(file.path(script_dir, "step08_logFC_scatter.R"))

if (!file.exists(step07_table)) {
  stop("Missing Step07 structure table. Please run locked Step07 first: ", step07_table)
}

structure <- read.csv(step07_table, stringsAsFactors = FALSE, check.names = FALSE)

shared <- structure[structure$in_1W & structure$in_4W, , drop = FALSE]
if (nrow(shared) == 0) stop("No shared strict genes found.")

shared$plot_class <- ifelse(
  shared$category == "shared_direction_discordant",
  "Opposite_direction",
  ifelse(shared$direction_1W == "Up" & shared$direction_4W == "Up",
         "Persistent_up",
         "Persistent_down")
)

shared$plot_class <- factor(
  shared$plot_class,
  levels = c("Persistent_up", "Persistent_down", "Opposite_direction")
)

shared$label_name <- ifelse(!is.na(shared$SYMBOL) & shared$SYMBOL != "", shared$SYMBOL, shared$gene_id)

shared$combined_score <- -log10(pmax(shared$FDR_1W, .Machine$double.xmin)) * abs(shared$logFC_1W) +
  -log10(pmax(shared$FDR_4W, .Machine$double.xmin)) * abs(shared$logFC_4W)

pearson_r <- unname(cor(shared$logFC_1W, shared$logFC_4W, method = "pearson"))
spearman_rho <- unname(cor(shared$logFC_1W, shared$logFC_4W, method = "spearman"))

n_shared <- nrow(shared)
n_persistent <- sum(shared$plot_class %in% c("Persistent_up", "Persistent_down"))
n_opposite <- sum(shared$plot_class == "Opposite_direction")
n_persistent_up <- sum(shared$plot_class == "Persistent_up")
n_persistent_down <- sum(shared$plot_class == "Persistent_down")

summary_df <- data.frame(
  metric = c(
    "shared_strict_genes",
    "persistent_same_direction",
    "opposite_direction",
    "persistent_up",
    "persistent_down",
    "pearson_r",
    "spearman_rho"
  ),
  value = c(
    n_shared,
    n_persistent,
    n_opposite,
    n_persistent_up,
    n_persistent_down,
    round(pearson_r, 4),
    round(spearman_rho, 4)
  ),
  stringsAsFactors = FALSE
)

write_csv(shared, file.path(table_dir, "step08_Figure2B_shared_strict_logFC_scatter_source_data.csv"))
write_csv(summary_df, file.path(table_dir, "step08_Figure2B_logFC_scatter_summary.csv"))

# Label selection:
# - top 10 Persistent_up
# - top 8 Persistent_down
# - all Opposite_direction genes
lab_up <- shared[shared$plot_class == "Persistent_up", , drop = FALSE]
lab_down <- shared[shared$plot_class == "Persistent_down", , drop = FALSE]
lab_opp <- shared[shared$plot_class == "Opposite_direction", , drop = FALSE]

lab_up <- lab_up[order(-lab_up$combined_score), , drop = FALSE]
lab_down <- lab_down[order(-lab_down$combined_score), , drop = FALSE]

label_df <- rbind(
  head(lab_up, 10),
  head(lab_down, 8),
  lab_opp
)
label_df <- label_df[!duplicated(label_df$gene_id), , drop = FALSE]

write_csv(label_df, file.path(table_dir, "step08_Figure2B_selected_labels.csv"))

# Position summary text.
x_min <- min(shared$logFC_1W, na.rm = TRUE)
x_max <- max(shared$logFC_1W, na.rm = TRUE)
y_min <- min(shared$logFC_4W, na.rm = TRUE)
y_max <- max(shared$logFC_4W, na.rm = TRUE)

summary_text <- paste0(
  "Shared strict genes = ", n_shared, "\n",
  "Persistent = ", n_persistent, "\n",
  "Opposite = ", n_opposite, "\n",
  "Pearson r = ", sprintf("%.2f", pearson_r), "\n",
  "Spearman rho = ", sprintf("%.2f", spearman_rho)
)

p <- ggplot2::ggplot(shared, ggplot2::aes(x = logFC_1W, y = logFC_4W)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.45, color = "grey45") +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey60") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey60") +
  ggplot2::geom_point(
    ggplot2::aes(color = plot_class),
    size = 1.55,
    alpha = 0.72
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    ggplot2::aes(label = label_name, color = plot_class),
    size = 3.0,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = x_min + 0.06 * (x_max - x_min),
    y = y_max - 0.08 * (y_max - y_min),
    label = summary_text,
    hjust = 0,
    vjust = 1,
    size = 4.0
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Persistent_up" = "#D55E00",
      "Persistent_down" = "#0072B2",
      "Opposite_direction" = "grey50"
    ),
    drop = FALSE
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    title = "Shared strict DEGs: 1W vs 4W effect-size concordance",
    x = "logFC at 1W (ACLR vs Contra)",
    y = "logFC at 4W (ACLR vs Contra)",
    color = NULL
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank()
  )

ggplot2::ggsave(file.path(figure_dir, "Figure2B_mouse_shared_strict_logFC_scatter.png"), p, width = 7.2, height = 5.8, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "Figure2B_mouse_shared_strict_logFC_scatter.pdf"), p, width = 7.2, height = 5.8)

software_versions <- data.frame(
  item = c("R", "ggplot2", "ggrepel"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("ggrepel"))
  ),
  stringsAsFactors = FALSE
)
write_csv(software_versions, file.path(table_dir, "step08_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step08_sessionInfo.txt"))

cat("\n===== STEP08 SUMMARY =====\n")
print(summary_df)

cat("\nSelected labels:\n")
print(label_df[, c("gene_id", "SYMBOL", "label_name", "logFC_1W", "logFC_4W", "plot_class", "combined_score")])

cat("\nFigure2B saved to:\n")
cat(file.path(figure_dir, "Figure2B_mouse_shared_strict_logFC_scatter.png"), "\n")
cat(file.path(figure_dir, "Figure2B_mouse_shared_strict_logFC_scatter.pdf"), "\n")

cat("\nStep08 completed successfully.\n")

sink()
