# =========================================================
# English note:
# This script generates Figure 1A as a simplified three-block
# overall study workflow with manual line breaks for readability.
# Each information item is placed on its own line, matching a
# manuscript-style workflow schematic.
# Output: exactly one PNG and one PDF.
# =========================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

# -----------------------------
# 1. Output paths
# -----------------------------
out_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/06_figures/Figure1"
log_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/08_logs"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(out_dir, "Figure1A_overall_study_workflow_3block_linebreak_v2.png")
pdf_file <- file.path(out_dir, "Figure1A_overall_study_workflow_3block_linebreak_v2.pdf")
summary_file <- file.path(log_dir, "step01_Figure1A_overall_study_workflow_3block_linebreak_v2_summary_to_send.txt")

# -----------------------------
# 2. Manual text content
# -----------------------------
mouse_body <- paste(
  "GSE271903 mouse synovium, paired ACLR vs Contra at 1W and 4W",
  "Paired DE analysis at 1W and 4W",
  "Shared strict DE genes at 1W & 4W: 1,421",
  "Persistent same-directional DE genes: 1,416",
  "Persistent same-directional DE genes ∩ CellAge: 78",
  "Projected CellAge score",
  "Hallmark GSEA: inflammatory response and EMT / ECM remodeling",
  sep = "\n"
)

early_pig_body <- paste(
  "E-MTAB-6664 early pig synovium, Control vs ACLT at 1W and 4W",
  "78 human genes projected to pig ortholog space",
  "75 detected in pig; 47 direction-consistent",
  "Persistent senescence-associated signature-score validation",
  "Hallmark GSEA: inflammatory response and EMT / ECM remodeling",
  "Early-defined strict core ortholog genes: 24",
  sep = "\n"
)

chronic_pig_body <- paste(
  "GSE228848 chronic pig synovium, Control_52W vs ACLT_alone_52W",
  "Persistent senescence-associated signature-score validation",
  "Hallmark GSEA: inflammatory response and EMT / ECM remodeling",
  "Retention audit of early-defined core genes performed",
  "No chronic-specific signature redefinition",
  sep = "\n"
)

blocks <- data.frame(
  title = c(
    "Mouse discovery",
    "Persistent program validation in early pig",
    "Pig chronic extension"
  ),
  body = c(mouse_body, early_pig_body, chronic_pig_body),
  x = c(10, 10, 10),
  y = c(17.6, 10.4, 3.7),
  w = c(16.6, 16.6, 16.6),
  h = c(5.6, 5.4, 4.8),
  header_fill = c("#2F5E8F", "#B56400", "#A03C76"),
  body_fill = c("#EAF2FB", "#FAECDD", "#F7E5F0"),
  stringsAsFactors = FALSE
)

blocks$xmin <- blocks$x - blocks$w / 2
blocks$xmax <- blocks$x + blocks$w / 2
blocks$ymin <- blocks$y - blocks$h / 2
blocks$ymax <- blocks$y + blocks$h / 2

header_h <- 0.9

body_df <- data.frame(
  x = blocks$x,
  y = blocks$y - 0.2,
  label = blocks$body,
  stringsAsFactors = FALSE
)

arrows_df <- data.frame(
  x = c(10, 10),
  y = c(blocks$ymin[1] - 0.18, blocks$ymin[2] - 0.18),
  xend = c(10, 10),
  yend = c(blocks$ymax[2] + 0.18, blocks$ymax[3] + 0.18)
)

# -----------------------------
# 3. Plot
# -----------------------------
p <- ggplot() +
  geom_rect(
    data = blocks,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = blocks$body_fill,
    color = "#666666",
    linewidth = 0.6
  ) +
  geom_rect(
    data = blocks,
    aes(xmin = xmin, xmax = xmax, ymin = ymax - header_h, ymax = ymax),
    fill = blocks$header_fill,
    color = "#666666",
    linewidth = 0.6
  ) +
  geom_text(
    data = blocks,
    aes(x = x, y = ymax - header_h / 2, label = title),
    color = "white",
    fontface = "bold",
    size = 7.6,
    lineheight = 0.95
  ) +
  geom_text(
    data = body_df,
    aes(x = x, y = y, label = label),
    color = "#222222",
    size = 4.9,
    lineheight = 1.18
  ) +
  geom_segment(
    data = arrows_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#666666",
    linewidth = 0.72,
    arrow = arrow(length = unit(0.24, "cm"), type = "closed")
  ) +
  annotate(
    "text",
    x = 10,
    y = 22.9,
    label = "Overall study workflow for cross-species synovial\nsenescence analysis after ACL injury",
    fontface = "bold",
    size = 8.5,
    lineheight = 1.0
  ) +
  annotate(
    "text",
    x = 10,
    y = 21.8,
    label = "Mouse discovery of a persistent program → early pig validation → chronic pig extension",
    size = 5.1,
    color = "#444444"
  ) +
  coord_cartesian(
    xlim = c(0.7, 19.3),
    ylim = c(0.3, 23.8),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(20, 22, 20, 22)
  )

# -----------------------------
# 4. Save outputs
# -----------------------------
ggsave(
  filename = png_file,
  plot = p,
  width = 12.5,
  height = 14.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 12.5,
  height = 14.5,
  units = "in",
  bg = "white"
)

# -----------------------------
# 5. Summary log
# -----------------------------
sink(summary_file)
cat("===== STEP01 FIGURE1A 3-BLOCK LINEBREAK V2 SUMMARY =====\n")
cat("Output PNG: ", png_file, "\n", sep = "")
cat("Output PDF: ", pdf_file, "\n", sep = "")
cat("\nBlock contents:\n\n")
for (i in seq_len(nrow(blocks))) {
  cat("[", i, "] ", blocks$title[i], "\n", sep = "")
  cat(blocks$body[i], "\n\n")
}
cat("Script completed successfully.\n")
sink()

cat("Figure1A 3-block linebreak v2 completed successfully.\n")
cat("PNG: ", png_file, "\n", sep = "")
cat("PDF: ", pdf_file, "\n", sep = "")
cat("Summary log: ", summary_file, "\n", sep = "")
