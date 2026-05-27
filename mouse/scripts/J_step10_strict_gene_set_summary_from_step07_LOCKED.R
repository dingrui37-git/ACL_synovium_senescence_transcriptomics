# Step10: Figure 2D strict DEG category summary rebuilt directly from Step07.
# English note:
# This locked revision removes the previous fallback/hard-coded count logic.
# It reads only the current Step07 strict DEG structure table and recomputes the
# Figure 2D counts from that table:
#   - Persistent = strict at both 1W and 4W with the same logFC direction
#   - 1W-only = strict only at 1W
#   - 4W-only = strict only at 4W
# Shared strict genes with opposite directions are excluded from Figure 2D and
# saved separately for audit. This script does not rerun DE analysis.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths and locked inputs
## =========================

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

step07_structure_file <- file.path(
  base_dir, "07_tables", "step07_strict_DEG_upset_persistent",
  "step07_strict_DEG_structure_gene_table.csv"
)

table_dir <- file.path(base_dir, "07_tables", "step10_strict_gene_set_summary_color_position_fixed")
figure_dir <- file.path(base_dir, "06_figures", "Figure2")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

for (d in c(table_dir, figure_dir, log_dir, script_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

step_name <- "step10_strict_gene_set_summary_color_position_fixed"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

png_file <- file.path(figure_dir, "Figure2D_strict_gene_set_summary_color_position_fixed.png")
pdf_file <- file.path(figure_dir, "Figure2D_strict_gene_set_summary_color_position_fixed.pdf")

## Canonical manuscript copies.
png_file_canonical <- file.path(figure_dir, "Figure2D_strict_gene_set_summary.png")
pdf_file_canonical <- file.path(figure_dir, "Figure2D_strict_gene_set_summary.pdf")

archive_file <- file.path(script_dir, "step10_strict_gene_set_summary_color_position_fixed.R")

## =========================
## 1. Clean only Step10 outputs
## =========================

old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  log_file,
  png_file,
  pdf_file,
  png_file_canonical,
  pdf_file_canonical,
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

cat("===== STEP10 FIGURE2D FROM CURRENT STEP07 STRUCTURE =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Base directory: ", base_dir, "\n", sep = "")
cat("Locked input: ", step07_structure_file, "\n\n", sep = "")
cat("Important: no old Step10 summary table and no hard-coded fallback counts are used.\n\n")

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
        "# Step10 archive fallback.",
        "# R could not detect the executed script path.",
        "# Please keep the externally executed script as the authoritative copy."
      ),
      path,
      useBytes = TRUE
    )
    cat("Script path not detected; archive fallback saved to: ", path, "\n", sep = "")
  }
}

safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path, row_names = FALSE) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (cc in colnames(x)) {
    if (is.list(x[[cc]])) {
      x[[cc]] <- vapply(x[[cc]], function(v) paste(as.character(v), collapse = ";"), character(1))
    }
  }
  write.csv(x, path, row.names = row_names, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

require_columns <- function(df, required_cols, label) {
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

normalize_direction <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("Up", "up", "UP", "Up in ACLR", "up_in_ACLR", "Persistent_up")] <- "Up in ACLR"
  x[x %in% c("Down", "down", "DOWN", "Down in ACLR", "down_in_ACLR", "Persistent_down")] <- "Down in ACLR"
  x[!x %in% c("Up in ACLR", "Down in ACLR")] <- NA_character_
  x
}

safe_library("ggplot2")
safe_library("grid")

archive_current_script(archive_file)

## =========================
## 3. Read current Step07 structure and rebuild counts
## =========================

if (!file.exists(step07_structure_file)) {
  stop("Missing required Step07 structure table: ", step07_structure_file)
}

structure_df <- read.csv(step07_structure_file, stringsAsFactors = FALSE, check.names = FALSE)
require_columns(
  structure_df,
  c("gene_id", "category", "in_1W", "in_4W", "direction_1W", "direction_4W"),
  "Step07 structure table"
)

cat("Step07 structure rows: ", nrow(structure_df), "\n", sep = "")
cat("Step07 categories detected:\n")
print(table(structure_df$category, useNA = "ifany"))

source_df <- structure_df

source_df$category_plot <- ifelse(
  source_df$category == "persistent_direction_consistent",
  "Persistent",
  ifelse(
    source_df$category == "1W_only",
    "1W-only",
    ifelse(source_df$category == "4W_only", "4W-only", NA_character_)
  )
)

source_df$direction_source <- ifelse(
  source_df$category %in% c("persistent_direction_consistent", "1W_only"),
  source_df$direction_1W,
  ifelse(source_df$category == "4W_only", source_df$direction_4W, NA_character_)
)
source_df$direction_plot <- normalize_direction(source_df$direction_source)

figure2d_source <- source_df[
  !is.na(source_df$category_plot) & !is.na(source_df$direction_plot),
  ,
  drop = FALSE
]

excluded_opposite <- source_df[source_df$category == "shared_direction_discordant", , drop = FALSE]
if (nrow(excluded_opposite) > 0) {
  write_csv(excluded_opposite, file.path(table_dir, "step10_shared_strict_opposite_direction_excluded_from_Figure2D.csv"))
}

if (nrow(figure2d_source) == 0) {
  stop("No Figure2D source rows could be derived from Step07 structure table.")
}

write_csv(
  figure2d_source,
  file.path(table_dir, "step10_Figure2D_source_rows_recomputed_from_step07.csv")
)

counts_long <- aggregate(
  gene_id ~ category_plot + direction_plot,
  data = figure2d_source,
  FUN = length
)
colnames(counts_long) <- c("category", "direction", "n")

all_combos <- expand.grid(
  category = c("Persistent", "1W-only", "4W-only"),
  direction = c("Down in ACLR", "Up in ACLR"),
  stringsAsFactors = FALSE
)

counts_long <- merge(all_combos, counts_long, by = c("category", "direction"), all.x = TRUE)
counts_long$n[is.na(counts_long$n)] <- 0L
counts_long$n <- as.integer(counts_long$n)

counts_long$category <- factor(counts_long$category, levels = c("Persistent", "1W-only", "4W-only"))
counts_long$direction <- factor(counts_long$direction, levels = c("Down in ACLR", "Up in ACLR"))
counts_long <- counts_long[order(counts_long$category, counts_long$direction), , drop = FALSE]

write_csv(counts_long, file.path(table_dir, "step10_counts_long_recomputed_from_step07.csv"))

down <- counts_long[counts_long$direction == "Down in ACLR", c("category", "n")]
up <- counts_long[counts_long$direction == "Up in ACLR", c("category", "n")]
wide_counts <- merge(down, up, by = "category", suffixes = c("_down", "_up"))

wide_counts$category <- as.character(wide_counts$category)
wide_counts <- wide_counts[match(c("Persistent", "1W-only", "4W-only"), wide_counts$category), , drop = FALSE]

colnames(wide_counts)[colnames(wide_counts) == "n_down"] <- "Down in ACLR"
colnames(wide_counts)[colnames(wide_counts) == "n_up"] <- "Up in ACLR"

wide_counts$x <- seq_len(nrow(wide_counts))
wide_counts$total <- wide_counts[["Down in ACLR"]] + wide_counts[["Up in ACLR"]]

write_csv(wide_counts, file.path(table_dir, "step10_counts_wide_recomputed_from_step07.csv"))

## Optional audit against the currently locked mouse numbers. This audit does not control the plot.
expected_current_counts <- data.frame(
  category = c("Persistent", "1W-only", "4W-only"),
  expected_down_in_ACLR = c(446L, 929L, 84L),
  expected_up_in_ACLR = c(970L, 1079L, 211L),
  stringsAsFactors = FALSE
)
audit_counts <- merge(
  expected_current_counts,
  wide_counts[, c("category", "Down in ACLR", "Up in ACLR", "total")],
  by = "category",
  all.x = TRUE
)
colnames(audit_counts)[colnames(audit_counts) == "Down in ACLR"] <- "observed_down_in_ACLR"
colnames(audit_counts)[colnames(audit_counts) == "Up in ACLR"] <- "observed_up_in_ACLR"
audit_counts$matches_current_locked_counts <- with(
  audit_counts,
  expected_down_in_ACLR == observed_down_in_ACLR &
    expected_up_in_ACLR == observed_up_in_ACLR
)
write_csv(audit_counts, file.path(table_dir, "step10_current_locked_count_audit_non_plot_controlling.csv"))

cat("\nRecomputed Figure2D counts from Step07:\n")
print(wide_counts)

cat("\nOpposite-direction shared strict genes excluded from Figure2D: ", nrow(excluded_opposite), "\n", sep = "")

## =========================
## 4. Manual block coordinates: down lower, up upper
## =========================

bar_half_width <- 0.36

rect_down <- data.frame(
  category = wide_counts$category,
  direction = "Down in ACLR",
  x = wide_counts$x,
  xmin = wide_counts$x - bar_half_width,
  xmax = wide_counts$x + bar_half_width,
  ymin = 0,
  ymax = wide_counts[["Down in ACLR"]],
  label_y = wide_counts[["Down in ACLR"]] / 2,
  n = wide_counts[["Down in ACLR"]],
  stringsAsFactors = FALSE
)

rect_up <- data.frame(
  category = wide_counts$category,
  direction = "Up in ACLR",
  x = wide_counts$x,
  xmin = wide_counts$x - bar_half_width,
  xmax = wide_counts$x + bar_half_width,
  ymin = wide_counts[["Down in ACLR"]],
  ymax = wide_counts$total,
  label_y = wide_counts[["Down in ACLR"]] + wide_counts[["Up in ACLR"]] / 2,
  n = wide_counts[["Up in ACLR"]],
  stringsAsFactors = FALSE
)

plot_rects <- rbind(rect_down, rect_up)
plot_rects$category <- factor(plot_rects$category, levels = c("Persistent", "1W-only", "4W-only"))
plot_rects$direction <- factor(plot_rects$direction, levels = c("Down in ACLR", "Up in ACLR"))
plot_rects$label_size <- ifelse(plot_rects$n < 120, 4.8, 5.6)

total_labels <- data.frame(
  category = wide_counts$category,
  x = wide_counts$x,
  y = wide_counts$total + max(wide_counts$total) * 0.035,
  label = paste0("Total = ", wide_counts$total),
  stringsAsFactors = FALSE
)
total_labels$category <- factor(total_labels$category, levels = c("Persistent", "1W-only", "4W-only"))

write_csv(plot_rects, file.path(table_dir, "step10_manual_rectangles_down_lower_up_upper.csv"))
write_csv(total_labels, file.path(table_dir, "step10_total_labels.csv"))

## =========================
## 5. Plot
## =========================

fill_values <- c(
  "Down in ACLR" = "#0072B2",
  "Up in ACLR" = "#D55E00"
)

p <- ggplot() +
  geom_rect(
    data = plot_rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = direction),
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    data = plot_rects[plot_rects$n > 0, , drop = FALSE],
    aes(x = x, y = label_y, label = n, size = label_size),
    color = "white",
    fontface = "plain"
  ) +
  scale_size_identity() +
  geom_text(
    data = total_labels,
    aes(x = x, y = y, label = label),
    color = "black",
    size = 5.6
  ) +
  scale_fill_manual(
    values = fill_values,
    breaks = c("Down in ACLR", "Up in ACLR"),
    labels = c("Down in ACLR", "Up in ACLR")
  ) +
  scale_x_continuous(
    breaks = wide_counts$x,
    labels = wide_counts$category,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    breaks = seq(0, ceiling(max(wide_counts$total) / 500) * 500, by = 500)
  ) +
  coord_cartesian(ylim = c(0, max(total_labels$y) * 1.05), clip = "off") +
  labs(
    title = "Strict DEG category summary across time points",
    x = NULL,
    y = "Number of genes",
    fill = NULL
  ) +
  guides(fill = guide_legend(
    reverse = FALSE,
    override.aes = list(color = "black")
  )) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.65, "cm"),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    plot.margin = margin(t = 12, r = 16, b = 10, l = 10)
  )

ggsave(png_file, p, width = 10.8, height = 7.8, dpi = 320)
ggsave(pdf_file, p, width = 10.8, height = 7.8)

ggsave(png_file_canonical, p, width = 10.8, height = 7.8, dpi = 320)
ggsave(pdf_file_canonical, p, width = 10.8, height = 7.8)

## =========================
## 6. Summary and reproducibility
## =========================

summary_df <- data.frame(
  metric = c(
    "step07_structure_file",
    "counts_source",
    "Persistent_down",
    "Persistent_up",
    "Persistent_total",
    "1W_only_down",
    "1W_only_up",
    "1W_only_total",
    "4W_only_down",
    "4W_only_up",
    "4W_only_total",
    "shared_strict_opposite_direction_excluded",
    "plot_logic"
  ),
  value = as.character(c(
    step07_structure_file,
    "recomputed directly from step07_strict_DEG_structure_gene_table.csv",
    wide_counts[wide_counts$category == "Persistent", "Down in ACLR"],
    wide_counts[wide_counts$category == "Persistent", "Up in ACLR"],
    wide_counts[wide_counts$category == "Persistent", "total"],
    wide_counts[wide_counts$category == "1W-only", "Down in ACLR"],
    wide_counts[wide_counts$category == "1W-only", "Up in ACLR"],
    wide_counts[wide_counts$category == "1W-only", "total"],
    wide_counts[wide_counts$category == "4W-only", "Down in ACLR"],
    wide_counts[wide_counts$category == "4W-only", "Up in ACLR"],
    wide_counts[wide_counts$category == "4W-only", "total"],
    nrow(excluded_opposite),
    "manual geom_rect: Down lower blue, Up upper orange"
  )),
  stringsAsFactors = FALSE
)
write_csv(summary_df, file.path(table_dir, "step10_color_position_fixed_summary.csv"))

versions <- data.frame(
  item = c("R", "ggplot2", "grid"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2")),
    "base R grid"
  ),
  stringsAsFactors = FALSE
)
write_csv(versions, file.path(table_dir, "step10_color_position_fixed_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step10_sessionInfo.txt"))

cat("\n===== STEP10 SUMMARY =====\n")
print(summary_df)

cat("\nOutput files:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat(png_file_canonical, "\n")
cat(pdf_file_canonical, "\n")

cat("\nStep10 completed successfully.\n")
