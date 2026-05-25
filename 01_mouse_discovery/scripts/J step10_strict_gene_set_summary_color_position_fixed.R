# step10_strict_gene_set_summary_color_position_fixed.R
# Purpose:
# Re-draw Figure 2D strict DEG category summary with manually fixed stacking:
#   - Down in ACLR is always the lower blue block.
#   - Up in ACLR is always the upper orange/red block.
#   - Each number is positioned at the center of its own block.
#   - Total label is positioned above each bar.
#
# This script does not rerun differential expression analysis.
# It reads existing Step10/Step07 summary tables when available and only fixes the Figure2D drawing layer.

options(stringsAsFactors = FALSE)

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

table_dir <- file.path(base_dir, "07_tables", "step10_strict_gene_set_summary_color_position_fixed")
figure_dir <- file.path(base_dir, "06_figures", "Figure2")
log_dir <- file.path(base_dir, "08_logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "step10_strict_gene_set_summary_color_position_fixed_log.txt")

png_file <- file.path(figure_dir, "Figure2D_strict_gene_set_summary_color_position_fixed.png")
pdf_file <- file.path(figure_dir, "Figure2D_strict_gene_set_summary_color_position_fixed.pdf")

## Also save a canonical copy for easy replacement if needed.
png_file_canonical <- file.path(figure_dir, "Figure2D_strict_gene_set_summary.png")
pdf_file_canonical <- file.path(figure_dir, "Figure2D_strict_gene_set_summary.pdf")

unlink(c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  log_file,
  png_file,
  pdf_file,
  png_file_canonical,
  pdf_file_canonical
)[file.exists(c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  log_file,
  png_file,
  pdf_file,
  png_file_canonical,
  pdf_file_canonical
))], recursive = TRUE, force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP10 FIGURE2D COLOR/POSITION FIX =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Base directory: ", base_dir, "\n", sep = "")

safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

normalize_category <- function(x) {
  x <- as.character(x)
  x <- gsub("_", "-", x)
  x <- trimws(x)
  x[x %in% c("Persistent", "persistent", "Persistent same direction", "persistent_direction_consistent")] <- "Persistent"
  x[x %in% c("1W-only", "1W only", "1W_only", "1Wonly")] <- "1W-only"
  x[x %in% c("4W-only", "4W only", "4W_only", "4Wonly")] <- "4W-only"
  x
}

normalize_direction <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("Down", "down", "Persistent_down", "persistent_down", "Down in ACLR", "down_in_ACLR")] <- "Down in ACLR"
  x[x %in% c("Up", "up", "Persistent_up", "persistent_up", "Up in ACLR", "up_in_ACLR")] <- "Up in ACLR"
  x
}

build_counts_from_long <- function(df) {
  colnames(df) <- make.names(colnames(df), unique = TRUE)

  category_col <- NULL
  for (cc in c("category", "Category", "gene_class", "class")) {
    cc2 <- make.names(cc)
    if (cc2 %in% colnames(df)) {
      category_col <- cc2
      break
    }
  }

  direction_col <- NULL
  for (cc in c("direction", "Direction", "regulation")) {
    cc2 <- make.names(cc)
    if (cc2 %in% colnames(df)) {
      direction_col <- cc2
      break
    }
  }

  n_col <- NULL
  for (cc in c("n", "count", "Count", "genes", "N")) {
    cc2 <- make.names(cc)
    if (cc2 %in% colnames(df)) {
      n_col <- cc2
      break
    }
  }

  if (is.null(category_col) || is.null(direction_col) || is.null(n_col)) return(NULL)

  out <- data.frame(
    category = normalize_category(df[[category_col]]),
    direction = normalize_direction(df[[direction_col]]),
    n = as.integer(df[[n_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[out$category %in% c("Persistent", "1W-only", "4W-only") &
               out$direction %in% c("Down in ACLR", "Up in ACLR"), , drop = FALSE]
  if (nrow(out) == 0) return(NULL)

  agg <- aggregate(n ~ category + direction, data = out, FUN = sum)
  agg
}

build_counts_from_wide <- function(df) {
  colnames_original <- colnames(df)
  colnames(df) <- make.names(colnames(df), unique = TRUE)

  category_col <- NULL
  for (cc in c("category", "Category")) {
    cc2 <- make.names(cc)
    if (cc2 %in% colnames(df)) {
      category_col <- cc2
      break
    }
  }
  if (is.null(category_col)) return(NULL)

  ## Find down/up columns after make.names().
  down_candidates <- make.names(c("Down in ACLR", "Down.in.ACLR", "down", "Down"))
  up_candidates <- make.names(c("Up in ACLR", "Up.in.ACLR", "up", "Up"))

  down_col <- down_candidates[down_candidates %in% colnames(df)]
  up_col <- up_candidates[up_candidates %in% colnames(df)]

  if (length(down_col) == 0 || length(up_col) == 0) return(NULL)

  out <- data.frame(
    category = normalize_category(df[[category_col]]),
    down = as.integer(df[[down_col[1]]]),
    up = as.integer(df[[up_col[1]]]),
    stringsAsFactors = FALSE
  )
  out <- out[out$category %in% c("Persistent", "1W-only", "4W-only"), , drop = FALSE]
  if (nrow(out) == 0) return(NULL)

  long <- rbind(
    data.frame(category = out$category, direction = "Down in ACLR", n = out$down, stringsAsFactors = FALSE),
    data.frame(category = out$category, direction = "Up in ACLR", n = out$up, stringsAsFactors = FALSE)
  )
  long
}

safe_library("ggplot2")
safe_library("grid")

## =========================
## 1. Read counts from existing locked summary files
## =========================

candidate_files <- c(
  file.path(base_dir, "07_tables", "step10_strict_gene_set_summary", "step10_strict_gene_set_summary_counts_wide.csv"),
  file.path(base_dir, "07_tables", "step10_strict_gene_set_summary", "step10_strict_gene_set_summary_counts_long.csv"),
  file.path(base_dir, "07_tables", "step10_strict_gene_set_summary", "step10_strict_gene_set_summary_run_summary.csv"),
  file.path(base_dir, "07_tables", "step10_strict_gene_set_summary", "step10_strict_gene_set_summary_manual_plot_rectangles.csv")
)

counts_long <- NULL
source_file <- NA_character_

for (ff in candidate_files) {
  tab <- safe_read(ff)
  if (is.null(tab)) next

  tmp <- build_counts_from_wide(tab)
  if (is.null(tmp)) tmp <- build_counts_from_long(tab)

  if (!is.null(tmp) && nrow(tmp) > 0) {
    counts_long <- tmp
    source_file <- ff
    break
  }
}

## Fallback to locked counts from current confirmed Figure2D if no table can be parsed.
## These values are the current locked mouse strict-DEG structure:
## Persistent: down=446, up=970; 1W-only: down=929, up=1079; 4W-only: down=84, up=211.
if (is.null(counts_long)) {
  counts_long <- data.frame(
    category = rep(c("Persistent", "1W-only", "4W-only"), each = 2),
    direction = rep(c("Down in ACLR", "Up in ACLR"), times = 3),
    n = c(446, 970, 929, 1079, 84, 211),
    stringsAsFactors = FALSE
  )
  source_file <- "locked_counts_fallback_from_confirmed_current_Figure2D"
}

cat("Counts source: ", source_file, "\n", sep = "")
cat("\nCounts long used:\n")
print(counts_long)

## Ensure complete combinations and order.
all_combos <- expand.grid(
  category = c("Persistent", "1W-only", "4W-only"),
  direction = c("Down in ACLR", "Up in ACLR"),
  stringsAsFactors = FALSE
)

counts_long <- merge(all_combos, counts_long, by = c("category", "direction"), all.x = TRUE)
counts_long$n[is.na(counts_long$n)] <- 0L

counts_long$category <- factor(counts_long$category, levels = c("Persistent", "1W-only", "4W-only"))
counts_long$direction <- factor(counts_long$direction, levels = c("Down in ACLR", "Up in ACLR"))
counts_long <- counts_long[order(counts_long$category, counts_long$direction), ]

write_csv(counts_long, file.path(table_dir, "step10_counts_long_used_for_color_position_fixed.csv"))

## Convert to wide.
down <- counts_long[counts_long$direction == "Down in ACLR", c("category", "n")]
up <- counts_long[counts_long$direction == "Up in ACLR", c("category", "n")]
wide_counts <- merge(down, up, by = "category", suffixes = c("_down", "_up"))
wide_counts$category <- as.character(wide_counts$category)
wide_counts <- wide_counts[match(c("Persistent", "1W-only", "4W-only"), wide_counts$category), ]

colnames(wide_counts)[colnames(wide_counts) == "n_down"] <- "Down in ACLR"
colnames(wide_counts)[colnames(wide_counts) == "n_up"] <- "Up in ACLR"

wide_counts$x <- seq_len(nrow(wide_counts))
wide_counts$total <- wide_counts[["Down in ACLR"]] + wide_counts[["Up in ACLR"]]

write_csv(wide_counts, file.path(table_dir, "step10_counts_wide_used_for_color_position_fixed.csv"))

## =========================
## 2. Manual block coordinates: down lower, up upper
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
## 3. Plot
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

## Also save canonical copy.
ggsave(png_file_canonical, p, width = 10.8, height = 7.8, dpi = 320)
ggsave(pdf_file_canonical, p, width = 10.8, height = 7.8)

summary_df <- data.frame(
  metric = c(
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
    "plot_logic"
  ),
  value = c(
    source_file,
    wide_counts[wide_counts$category == "Persistent", "Down in ACLR"],
    wide_counts[wide_counts$category == "Persistent", "Up in ACLR"],
    wide_counts[wide_counts$category == "Persistent", "total"],
    wide_counts[wide_counts$category == "1W-only", "Down in ACLR"],
    wide_counts[wide_counts$category == "1W-only", "Up in ACLR"],
    wide_counts[wide_counts$category == "1W-only", "total"],
    wide_counts[wide_counts$category == "4W-only", "Down in ACLR"],
    wide_counts[wide_counts$category == "4W-only", "Up in ACLR"],
    wide_counts[wide_counts$category == "4W-only", "total"],
    "manual geom_rect: Down lower blue, Up upper orange"
  ),
  stringsAsFactors = FALSE
)

write_csv(summary_df, file.path(table_dir, "step10_color_position_fixed_summary.csv"))

versions <- data.frame(
  item = c("R", "ggplot2"),
  version = c(as.character(getRversion()), as.character(utils::packageVersion("ggplot2"))),
  stringsAsFactors = FALSE
)
write_csv(versions, file.path(table_dir, "step10_color_position_fixed_software_versions.csv"))

cat("\n===== STEP10 COLOR/POSITION FIX SUMMARY =====\n")
print(summary_df)

cat("\nOutput files:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat(png_file_canonical, "\n")
cat(pdf_file_canonical, "\n")

cat("\nStep10 color/position fix completed successfully.\n")

sink()

cat("\nStep10 color/position fix completed. Please open:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
