# Chronic Step6 / Figure5D current78: early vs chronic summary comparison
# Purpose:
#   Generate Figure5D as a compact summary-comparison panel integrating the key
#   early-pig and chronic-pig findings under the current78 / 75-gene / 24-core framework.
#
# Figure design:
#   A summary tile plot comparing:
#     1) Signature directional score
#     2) Hallmark inflammatory response NES
#     3) Hallmark EMT / ECM remodeling NES
#     4) Core-gene status
#   across:
#     - Early 1W
#     - Early 4W
#     - Chronic 52W
#
# Interpretation:
#   - This is a summary/comparison panel, not a re-analysis panel.
#   - It reuses already completed early and chronic outputs.
#   - The chronic cohort remains an extension-validation cohort; no chronic-specific
#     signature or core set is redefined here.
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

pig_early_dir <- file.path(rebuild_root, "02_pig_early")
pig_chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

out_table_dir <- file.path(pig_chronic_dir, "tables", "chronic_step6_current78_figure5D_early_vs_chronic_summary")
out_fig_dir <- file.path(pig_chronic_dir, "figures", "Figure5D_current78_early_vs_chronic_summary")
log_dir <- file.path(out_table_dir, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

full_log_file <- file.path(log_dir, "chronic_step6_current78_figure5D_full_log.txt")
summary_log_file <- file.path(log_dir, "chronic_step6_current78_figure5D_summary_to_send_me.txt")

sink(full_log_file, split = TRUE)

cat("===== CHRONIC STEP6 CURRENT78 FIGURE5D EARLY VS CHRONIC SUMMARY =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Purpose: generate Figure5D as an early-vs-chronic summary comparison panel.\n\n")

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

get_metric_value <- function(df, metric_name) {
  hit <- df$value[df$metric == metric_name]
  if (length(hit) == 0) return(NA)
  hit[1]
}

fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_fdr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (is.na(x)) return("FDR = NA")
  if (x < 0.001) return("FDR < 0.001")
  paste0("FDR = ", formatC(x, format = "f", digits = 3))
}

pick_fill_for_effect <- function(value, fdr) {
  value <- suppressWarnings(as.numeric(value))
  fdr <- suppressWarnings(as.numeric(fdr))
  if (is.na(value)) return("#EFEFEF")
  if (is.na(fdr)) {
    if (value > 0) return("#F7E3CC")
    if (value < 0) return("#DCEAF7")
    return("#EFEFEF")
  }
  if (value > 0) {
    if (fdr < 0.05) return("#E49A4A")
    if (fdr < 0.10) return("#F2BF84")
    return("#F7E3CC")
  }
  if (value < 0) {
    if (fdr < 0.05) return("#6BAED6")
    if (fdr < 0.10) return("#A9D0EA")
    return("#DCEAF7")
  }
  "#EFEFEF"
}

pick_fill_for_core_fraction <- function(frac) {
  frac <- suppressWarnings(as.numeric(frac))
  if (is.na(frac)) return("#EFEFEF")
  if (frac >= 0.40) return("#E49A4A")
  if (frac >= 0.25) return("#F2BF84")
  if (frac >= 0.10) return("#F7E3CC")
  "#FBF3EA"
}

## =========================
## 2. Input files
## =========================

# Early
early_sig_stats_file <- file.path(pig_early_dir, "tables", "step19_current78_pig_signature_score", "step19_current78_pig_signature_score_group_comparisons.csv")
early_gsea_targets_file <- file.path(pig_early_dir, "tables", "step20_current78_pig_early_hallmark_gsea", "step20_current78_target_pathway_results_from_full_hallmark.csv")
early_step18_summary_file <- file.path(pig_early_dir, "tables", "step18_current78_pig_early_signature_remap", "step18_current78_pig_early_signature_remap_summary.csv")

# Chronic
chronic_sig_stats_file <- file.path(
  pig_chronic_dir,
  "tables",
  "chronic_step2_current78_signature_score_TMM_aligned",
  "chronic_step2_current78_TMM_aligned_signature_score_group_comparisons.csv"
)
chronic_gsea_targets_file <- file.path(pig_chronic_dir, "tables", "chronic_step4_current78_DE_GSEA", "chronic_step4_target_pathway_results_from_full_hallmark.csv")
# Use the final Step4-locked Step5B core-gene DE audit, not the older non-locked Step5B folder.
chronic_step5b_summary_file <- file.path(pig_chronic_dir, "tables", "chronic_step5B_current78_core_DE_audit_Step4_locked", "chronic_step5B_early_defined_24_core_chronic_DE_audit_summary.csv")

stop_if_missing(early_sig_stats_file, "early Step19 signature score statistics")
stop_if_missing(early_gsea_targets_file, "early Step20 target pathway results")
stop_if_missing(early_step18_summary_file, "early Step18 remap summary")
stop_if_missing(chronic_sig_stats_file, "chronic Step2 signature score statistics")
stop_if_missing(chronic_gsea_targets_file, "chronic Step4 target pathway results")
stop_if_missing(chronic_step5b_summary_file, "chronic Step5B core DE audit summary")

early_sig_stats <- safe_read_csv(early_sig_stats_file)
early_gsea <- safe_read_csv(early_gsea_targets_file)
early_step18_summary <- safe_read_csv(early_step18_summary_file)

chronic_sig_stats <- safe_read_csv(chronic_sig_stats_file)
chronic_gsea <- safe_read_csv(chronic_gsea_targets_file)
chronic_step5b_summary <- safe_read_csv(chronic_step5b_summary_file)

cat("Early signature statistics: ", early_sig_stats_file, "\n", sep = "")
cat("Early target GSEA: ", early_gsea_targets_file, "\n", sep = "")
cat("Early core summary: ", early_step18_summary_file, "\n", sep = "")
cat("Chronic signature statistics: ", chronic_sig_stats_file, "\n", sep = "")
cat("Chronic target GSEA: ", chronic_gsea_targets_file, "\n", sep = "")
cat("Chronic core audit summary: ", chronic_step5b_summary_file, "\n\n", sep = "")

## =========================
## 3. Extract summary values
## =========================

# --- Signature directional score
early_sig_t7 <- early_sig_stats[
  early_sig_stats$score == "directional_score" &
    early_sig_stats$comparison == "ACLT_t7_vs_Control",
  ,
  drop = FALSE
]
early_sig_t28 <- early_sig_stats[
  early_sig_stats$score == "directional_score" &
    early_sig_stats$comparison == "ACLT_t28_vs_Control",
  ,
  drop = FALSE
]
chronic_sig <- chronic_sig_stats[
  chronic_sig_stats$score == "directional_score" &
    chronic_sig_stats$comparison == "ACLT_alone_52W_vs_Control_52W",
  ,
  drop = FALSE
]

if (nrow(early_sig_t7) != 1 || nrow(early_sig_t28) != 1 || nrow(chronic_sig) != 1) {
  stop("Could not uniquely recover directional_score statistics for early/chronic comparisons.", call. = FALSE)
}

# --- GSEA target results
get_gsea_row <- function(df, contrast_name, pathway_name) {
  rr <- df[df$contrast == contrast_name & df$pathway == pathway_name, , drop = FALSE]
  if (nrow(rr) != 1) stop("Could not recover GSEA row: contrast=", contrast_name, ", pathway=", pathway_name, call. = FALSE)
  rr
}

early_inf_t7 <- get_gsea_row(early_gsea, "ACLT_untreated_t7_vs_Control", "HALLMARK_INFLAMMATORY_RESPONSE")
early_inf_t28 <- get_gsea_row(early_gsea, "ACLT_untreated_t28_vs_Control", "HALLMARK_INFLAMMATORY_RESPONSE")
early_emt_t7 <- get_gsea_row(early_gsea, "ACLT_untreated_t7_vs_Control", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
early_emt_t28 <- get_gsea_row(early_gsea, "ACLT_untreated_t28_vs_Control", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

chronic_inf <- get_gsea_row(chronic_gsea, "ACLT_alone_52W_vs_Control_52W", "HALLMARK_INFLAMMATORY_RESPONSE")
chronic_emt <- get_gsea_row(chronic_gsea, "ACLT_alone_52W_vs_Control_52W", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

# --- Early/core summary
early_detected <- suppressWarnings(as.numeric(get_metric_value(early_step18_summary, "current_pig_signature_genes_detected")))
early_strict_t7 <- suppressWarnings(as.numeric(get_metric_value(early_step18_summary, "strict_t7")))
early_strict_t28 <- suppressWarnings(as.numeric(get_metric_value(early_step18_summary, "strict_t28")))
early_core_24 <- suppressWarnings(as.numeric(get_metric_value(early_step18_summary, "core_strict_both_t7_t28")))

# --- Chronic core audit summary
chronic_core_tested <- suppressWarnings(as.numeric(get_metric_value(chronic_step5b_summary, "core_genes_tested_in_chronic_DE")))
chronic_core_dir <- suppressWarnings(as.numeric(get_metric_value(chronic_step5b_summary, "direction_consistent_with_early")))
chronic_core_strict_dir <- suppressWarnings(as.numeric(get_metric_value(chronic_step5b_summary, "chronic_strict_and_direction_consistent")))

if (any(is.na(c(early_detected, early_strict_t7, early_strict_t28, early_core_24,
                chronic_core_tested, chronic_core_dir, chronic_core_strict_dir)))) {
  stop("Failed to recover one or more early/chronic core-summary metrics.", call. = FALSE)
}

## =========================
## 4. Build Figure5D tile data
## =========================

tile_df <- data.frame(
  stage = rep(c("Early 1W", "Early 4W", "Chronic 52W"), each = 4),
  metric = rep(c(
    "Signature\ndirectional score",
    "Hallmark NES:\nInflammatory response",
    "Hallmark NES:\nEMT / ECM remodeling",
    "Core-gene\nstatus"
  ), times = 3),
  stringsAsFactors = FALSE
)

tile_df$stage <- factor(tile_df$stage, levels = c("Early 1W", "Early 4W", "Chronic 52W"))
tile_df$metric <- factor(tile_df$metric, levels = rev(c(
  "Signature\ndirectional score",
  "Hallmark NES:\nInflammatory response",
  "Hallmark NES:\nEMT / ECM remodeling",
  "Core-gene\nstatus"
)))

tile_df$label <- NA_character_
tile_df$fill_color <- "#EFEFEF"
tile_df$numeric_value <- NA_real_
tile_df$fdr_value <- NA_real_
tile_df$note <- NA_character_

# Fill each tile
for (i in seq_len(nrow(tile_df))) {
  st <- as.character(tile_df$stage[i])
  met <- as.character(tile_df$metric[i])

  if (met == "Signature\ndirectional score") {
    if (st == "Early 1W") {
      effect <- as.numeric(early_sig_t7$median_difference_case_minus_control[1])
      fdr <- as.numeric(early_sig_t7$BH_FDR_across_all_score_comparisons[1])
      tile_df$label[i] <- paste0("Δscore = ", fmt_num(effect), "\n", fmt_fdr(fdr))
      tile_df$fill_color[i] <- pick_fill_for_effect(effect, fdr)
      tile_df$numeric_value[i] <- effect
      tile_df$fdr_value[i] <- fdr
      tile_df$note[i] <- "Early 1W directional score vs Control"
    } else if (st == "Early 4W") {
      effect <- as.numeric(early_sig_t28$median_difference_case_minus_control[1])
      fdr <- as.numeric(early_sig_t28$BH_FDR_across_all_score_comparisons[1])
      tile_df$label[i] <- paste0("Δscore = ", fmt_num(effect), "\n", fmt_fdr(fdr))
      tile_df$fill_color[i] <- pick_fill_for_effect(effect, fdr)
      tile_df$numeric_value[i] <- effect
      tile_df$fdr_value[i] <- fdr
      tile_df$note[i] <- "Early 4W directional score vs Control"
    } else if (st == "Chronic 52W") {
      effect <- as.numeric(chronic_sig$median_difference_case_minus_control[1])
      fdr <- as.numeric(chronic_sig$BH_FDR_across_all_score_comparisons[1])
      tile_df$label[i] <- paste0("Δscore = ", fmt_num(effect), "\n", fmt_fdr(fdr))
      tile_df$fill_color[i] <- pick_fill_for_effect(effect, fdr)
      tile_df$numeric_value[i] <- effect
      tile_df$fdr_value[i] <- fdr
      tile_df$note[i] <- "Chronic 52W directional score vs Control"
    }
  }

  if (met == "Hallmark NES:\nInflammatory response") {
    if (st == "Early 1W") {
      nes <- as.numeric(early_inf_t7$NES[1]); fdr <- as.numeric(early_inf_t7$padj[1])
    } else if (st == "Early 4W") {
      nes <- as.numeric(early_inf_t28$NES[1]); fdr <- as.numeric(early_inf_t28$padj[1])
    } else {
      nes <- as.numeric(chronic_inf$NES[1]); fdr <- as.numeric(chronic_inf$padj[1])
    }
    tile_df$label[i] <- paste0("NES = ", fmt_num(nes), "\n", fmt_fdr(fdr))
    tile_df$fill_color[i] <- pick_fill_for_effect(nes, fdr)
    tile_df$numeric_value[i] <- nes
    tile_df$fdr_value[i] <- fdr
    tile_df$note[i] <- "Hallmark inflammatory response"
  }

  if (met == "Hallmark NES:\nEMT / ECM remodeling") {
    if (st == "Early 1W") {
      nes <- as.numeric(early_emt_t7$NES[1]); fdr <- as.numeric(early_emt_t7$padj[1])
    } else if (st == "Early 4W") {
      nes <- as.numeric(early_emt_t28$NES[1]); fdr <- as.numeric(early_emt_t28$padj[1])
    } else {
      nes <- as.numeric(chronic_emt$NES[1]); fdr <- as.numeric(chronic_emt$padj[1])
    }
    tile_df$label[i] <- paste0("NES = ", fmt_num(nes), "\n", fmt_fdr(fdr))
    tile_df$fill_color[i] <- pick_fill_for_effect(nes, fdr)
    tile_df$numeric_value[i] <- nes
    tile_df$fdr_value[i] <- fdr
    tile_df$note[i] <- "Hallmark EMT / ECM remodeling"
  }

  if (met == "Core-gene\nstatus") {
    if (st == "Early 1W") {
      frac <- early_strict_t7 / early_detected
      tile_df$label[i] <- paste0("Strict at 1W\n", early_strict_t7, " / ", early_detected)
      tile_df$fill_color[i] <- pick_fill_for_core_fraction(frac)
      tile_df$numeric_value[i] <- frac
      tile_df$note[i] <- "Among the fixed 75 pig signature genes, genes strict at early 1W"
    } else if (st == "Early 4W") {
      frac <- early_strict_t28 / early_detected
      tile_df$label[i] <- paste0("Strict at 4W\n", early_strict_t28, " / ", early_detected)
      tile_df$fill_color[i] <- pick_fill_for_core_fraction(frac)
      tile_df$numeric_value[i] <- frac
      tile_df$note[i] <- "Among the fixed 75 pig signature genes, genes strict at early 4W"
    } else {
      frac <- chronic_core_strict_dir / early_core_24
      tile_df$label[i] <- paste0(
        "Strict retained\n", chronic_core_strict_dir, " / ", early_core_24,
        "\n(", chronic_core_dir, "/", chronic_core_tested, " dir-cons.)"
      )
      tile_df$fill_color[i] <- pick_fill_for_core_fraction(frac)
      tile_df$numeric_value[i] <- frac
      tile_df$note[i] <- "Among the fixed early-defined 24 core genes, strict + direction-consistent in chronic DE"
    }
  }
}

## =========================
## 5. Plot
## =========================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
suppressPackageStartupMessages(library(ggplot2))

p <- ggplot(tile_df, aes(x = stage, y = metric)) +
  geom_tile(aes(fill = fill_color), color = "black", linewidth = 0.35, width = 0.95, height = 0.95, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(aes(label = label), size = 4.1, lineheight = 0.95) +
  labs(
    title = "Figure 5D. Early-to-chronic comparison of the current78-derived pig synovial program",
    subtitle = "Fixed 75-gene signature and fixed early-defined 24 core genes projected from early to chronic pig synovium",
    x = NULL,
    y = NULL,
    caption = paste(
      "Signature row: median directional-score difference (case - control) with FDR.",
      "Hallmark rows: NES with FDR from full Hallmark fgsea.",
      "Core-gene row: early strict counts among 75 signature genes, and chronic strict retained counts among the fixed early-defined 24 core genes.",
      sep = "\n"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11.5),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11.5),
    plot.caption = element_text(size = 9.5, hjust = 0),
    legend.position = "none"
  )

figure_png <- file.path(out_fig_dir, "Figure5D_current78_early_vs_chronic_summary_tileplot.png")
figure_pdf <- file.path(out_fig_dir, "Figure5D_current78_early_vs_chronic_summary_tileplot.pdf")

ggsave(figure_png, p, width = 11.4, height = 6.4, dpi = 320)
ggsave(figure_pdf, p, width = 11.4, height = 6.4)

cat("Saved plot: ", figure_png, "\n", sep = "")
cat("Saved plot: ", figure_pdf, "\n\n", sep = "")

## =========================
## 6. Save source data and summaries
## =========================

source_data_file <- file.path(out_table_dir, "chronic_step6_current78_figure5D_plot_source_data.csv")
summary_file <- file.path(out_table_dir, "chronic_step6_current78_figure5D_summary.csv")
version_file <- file.path(out_table_dir, "chronic_step6_current78_versions_and_method_records.csv")
comparison_values_file <- file.path(out_table_dir, "chronic_step6_current78_early_vs_chronic_extracted_values.csv")

write_csv(tile_df, source_data_file)

comparison_values <- data.frame(
  category = c(
    "early_signature_directional_score_t7",
    "early_signature_directional_score_t28",
    "chronic_signature_directional_score",
    "early_inflammatory_NES_t7",
    "early_inflammatory_NES_t28",
    "chronic_inflammatory_NES",
    "early_EMT_NES_t7",
    "early_EMT_NES_t28",
    "chronic_EMT_NES",
    "early_signature_detected_genes",
    "early_strict_t7",
    "early_strict_t28",
    "early_core_genes",
    "chronic_core_genes_tested",
    "chronic_core_direction_consistent",
    "chronic_core_strict_direction_consistent"
  ),
  value = c(
    early_sig_t7$median_difference_case_minus_control[1],
    early_sig_t28$median_difference_case_minus_control[1],
    chronic_sig$median_difference_case_minus_control[1],
    early_inf_t7$NES[1],
    early_inf_t28$NES[1],
    chronic_inf$NES[1],
    early_emt_t7$NES[1],
    early_emt_t28$NES[1],
    chronic_emt$NES[1],
    early_detected,
    early_strict_t7,
    early_strict_t28,
    early_core_24,
    chronic_core_tested,
    chronic_core_dir,
    chronic_core_strict_dir
  ),
  stringsAsFactors = FALSE
)
write_csv(comparison_values, comparison_values_file)

summary_df <- data.frame(
  metric = c(
    "early_signature_stats_file",
    "early_gsea_targets_file",
    "early_step18_summary_file",
    "chronic_signature_stats_file",
    "chronic_gsea_targets_file",
    "chronic_step5B_summary_file",
    "figure_png",
    "figure_pdf",
    "plot_source_data_file",
    "extracted_values_file",
    "interpretation_policy",
    "output_dir"
  ),
  value = c(
    early_sig_stats_file,
    early_gsea_targets_file,
    early_step18_summary_file,
    chronic_sig_stats_file,
    chronic_gsea_targets_file,
    chronic_step5b_summary_file,
    figure_png,
    figure_pdf,
    source_data_file,
    comparison_values_file,
    "Figure5D is a summary-comparison panel integrating completed early and chronic outputs under the fixed current78 / 75-gene / 24-core framework. It does not redefine a chronic signature or chronic core.",
    out_table_dir
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_df, summary_file)

version_df <- data.frame(
  item = c(
    "R_version",
    "ggplot2_version",
    "early_signature_source_step",
    "early_GSEA_source_step",
    "chronic_signature_source_step",
    "chronic_GSEA_source_step",
    "core_comparison_policy"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2")),
    "Step19 current78 pig early signature score",
    "Step20 current78 pig early full-Hallmark GSEA",
    "Chronic Step2 current78 TMM-aligned signature score",
    "Chronic Step4 current78 edgeR QLF DE + full-Hallmark GSEA",
    "Compare early strict counts among the fixed 75-gene signature with chronic strict retained counts among the fixed early-defined 24 core genes; no chronic core redefinition."
  ),
  stringsAsFactors = FALSE
)
write_csv(version_df, version_file)

cat("\n===== Chronic Step6 Figure5D summary =====\n")
print(summary_df, row.names = FALSE)

cat("\nExtracted comparison values:\n")
print(comparison_values, row.names = FALSE)

cat("\nTile plot source data:\n")
print(tile_df, row.names = FALSE)

cat("\nVersion and method records:\n")
print(version_df, row.names = FALSE)

summary_con <- file(summary_log_file, open = "wt")
writeLines("===== Chronic Step6 current78 Figure5D early vs chronic summary TO SEND ME =====", summary_con)
writeLines(paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), summary_con)
writeLines("", summary_con)
writeLines("Main run summary:", summary_con)
writeLines(capture.output(print(summary_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Extracted comparison values:", summary_con)
writeLines(capture.output(print(comparison_values, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Tile plot source data:", summary_con)
writeLines(capture.output(print(tile_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Version and method records:", summary_con)
writeLines(capture.output(print(version_df, row.names = FALSE)), summary_con)
writeLines("", summary_con)
writeLines("Interpretation checkpoint:", summary_con)
writeLines("This figure summarizes completed early and chronic analyses under the fixed current78 / 75-gene / 24-core framework.", summary_con)
writeLines("The chronic cohort remains an extension-validation dataset.", summary_con)
writeLines("No chronic-specific signature or chronic core set is redefined in this figure.", summary_con)
close(summary_con)

sink()

cat("\nChronic Step6 / Figure5D current78 early-vs-chronic summary completed.\n")
cat("Summary to send me: ", summary_log_file, "\n", sep = "")
