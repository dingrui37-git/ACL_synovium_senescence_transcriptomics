options(stringsAsFactors = FALSE)
required <- c("ggplot2", "ragg", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

out_root <- "D:/workspace/ACLsenescence2_reviewer2_figure_reorganization"
out_dir <- file.path(out_root, "supplementary_figure5_directional_null")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(p) utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
fmt_p <- function(x) ifelse(x < 1e-4, "<0.0001", sprintf("%.4f", x))
fmt_fdr <- function(x) ifelse(x < 1e-4, "<0.0001", sprintf("%.3g", x))
export_ggplot <- function(p, base, width_mm, height_mm) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(base, ".svg"), width = w, height = h); print(p); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h); print(p); grDevices::dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(p); grDevices::dev.off()
  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = 300); print(p); grDevices::dev.off()
}

in_dir <- "D:/workspace/ACLsenescence2_reviewer1_directional_score_null"
summ <- read_csv(file.path(in_dir, "directional_score_matched_null_summary.csv"))
perm <- read_csv(file.path(in_dir, "directional_score_matched_null_all_permutations.csv"))
summ$display <- c("Early pig | 1W", "Early pig | 4W", "Chronic pig | 52W")
perm$display <- ifelse(
  perm$dataset == "pig_early" & perm$comparison == "ACLT_t7_vs_Control", "Early pig | 1W",
  ifelse(perm$dataset == "pig_early" & perm$comparison == "ACLT_t28_vs_Control", "Early pig | 4W", "Chronic pig | 52W")
)
lev <- c("Early pig | 1W", "Early pig | 4W", "Chronic pig | 52W")
summ$display <- factor(summ$display, levels = lev)
perm$display <- factor(perm$display, levels = lev)
ann <- summ
ann$label <- paste0("Observed difference = ", sprintf("%.2f", ann$observed_median_difference),
                    "\nP = ", fmt_p(ann$empirical_p_one_sided_activation),
                    "; BH-FDR = ", fmt_fdr(ann$BH_FDR_one_sided_across_three_comparisons))
ann$label_x <- 1.04

p <- ggplot2::ggplot(perm, ggplot2::aes(x = null_median_difference)) +
  ggplot2::geom_density(fill = "#BFD7EA", colour = "#2166AC", linewidth = 0.55, adjust = 1.05) +
  ggplot2::geom_vline(data = summ, ggplot2::aes(xintercept = null_q025), colour = "#7F8C8D", linetype = "dashed", linewidth = 0.35) +
  ggplot2::geom_vline(data = summ, ggplot2::aes(xintercept = null_q975), colour = "#7F8C8D", linetype = "dashed", linewidth = 0.35) +
  ggplot2::geom_vline(data = summ, ggplot2::aes(xintercept = 0), colour = "#404040", linetype = "dotted", linewidth = 0.35) +
  ggplot2::geom_vline(data = summ, ggplot2::aes(xintercept = observed_median_difference), colour = "#B2182B", linewidth = 0.85) +
  ggplot2::geom_label(data = ann, ggplot2::aes(x = label_x, y = Inf, label = label), hjust = 1, vjust = 1.15, lineheight = 0.95, size = 2.35, colour = "#222222", fill = "white", alpha = 0.82, linewidth = 0.15, label.padding = grid::unit(0.12, "lines")) +
  ggplot2::facet_wrap(~ display, nrow = 1, scales = "free_y") +
  ggplot2::scale_x_continuous(name = "Median directional-score difference (case - control)", breaks = seq(-0.5, 1.0, by = 0.5), expand = c(0.01, 0.01)) +
  ggplot2::coord_cartesian(xlim = c(-0.50, 1.08)) +
  ggplot2::labs(y = "Matched-null density", title = "Signature-specific directional activation in pig") +
  ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
                 strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
                 strip.text = ggplot2::element_text(face = "bold", size = 8),
                 axis.title = ggplot2::element_text(size = 8.5),
                 axis.text = ggplot2::element_text(size = 7.5),
                 panel.spacing = grid::unit(0.75, "lines"))

base <- file.path(out_dir, "SupplementaryFigure5_DirectionalScoreMatchedNull")
export_ggplot(p, base, width_mm = 190, height_mm = 82)
make_subset_plot <- function(label_subset) {
  ss <- summ[summ$display %in% label_subset, , drop = FALSE]
  pp <- perm[perm$display %in% label_subset, , drop = FALSE]
  aa <- ann[ann$display %in% label_subset, , drop = FALSE]
  ss$display <- factor(as.character(ss$display), levels = label_subset)
  pp$display <- factor(as.character(pp$display), levels = label_subset)
  aa$display <- factor(as.character(aa$display), levels = label_subset)
  title_text <- if (length(label_subset) == 1L) "Directional activation | Pig 52W" else "Signature-specific directional activation in pig"
  ggplot2::ggplot(pp, ggplot2::aes(x = null_median_difference)) +
    ggplot2::geom_density(fill = "#BFD7EA", colour = "#2166AC", linewidth = 0.55, adjust = 1.05) +
    ggplot2::geom_vline(data = ss, ggplot2::aes(xintercept = null_q025), colour = "#7F8C8D", linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_vline(data = ss, ggplot2::aes(xintercept = null_q975), colour = "#7F8C8D", linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_vline(data = ss, ggplot2::aes(xintercept = 0), colour = "#404040", linetype = "dotted", linewidth = 0.35) +
    ggplot2::geom_vline(data = ss, ggplot2::aes(xintercept = observed_median_difference), colour = "#B2182B", linewidth = 0.85) +
    ggplot2::geom_label(data = aa, ggplot2::aes(x = label_x, y = Inf, label = label), hjust = 1, vjust = 1.15, lineheight = 0.95, size = 2.35, colour = "#222222", fill = "white", alpha = 0.82, linewidth = 0.15, label.padding = grid::unit(0.12, "lines")) +
    ggplot2::facet_wrap(~ display, nrow = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(name = "Median directional-score difference (case - control)", breaks = seq(-0.5, 1.0, by = 0.5), expand = c(0.01, 0.01)) +
    ggplot2::coord_cartesian(xlim = c(-0.50, 1.08)) +
    ggplot2::labs(y = "Matched-null density", title = title_text) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
                   strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
                   strip.text = ggplot2::element_text(face = "bold", size = 8),
                   axis.title = ggplot2::element_text(size = 8.5), axis.text = ggplot2::element_text(size = 7.5),
                   panel.spacing = grid::unit(0.75, "lines"))
}
export_ggplot(make_subset_plot(c("Early pig | 1W", "Early pig | 4W")),
              file.path(out_dir, "SupplementaryFigure5_DirectionalScoreMatchedNull_EarlyPig"),
              width_mm = 145, height_mm = 82)
export_ggplot(make_subset_plot("Chronic pig | 52W"),
              file.path(out_dir, "SupplementaryFigure5_DirectionalScoreMatchedNull_ChronicPig"),
              width_mm = 105, height_mm = 82)
utils::write.csv(summ, file.path(out_dir, "SupplementaryFigure5_plot_summary.csv"), row.names = FALSE)
utils::write.csv(perm, file.path(out_dir, "SupplementaryFigure5_plot_source_data.csv"), row.names = FALSE)
writeLines(c("Supplementary Figure 5: matched-gene-set directional-score null control.",
             "Shaded density: 10,000 matched random-gene-set effects.",
             "Dashed lines: null 95% interval; dotted line: zero; red line: observed signature effect.",
             "BH-FDR is across the three prespecified pig comparisons."),
           file.path(out_dir, "directional_score_null_plot_contract.txt"))
