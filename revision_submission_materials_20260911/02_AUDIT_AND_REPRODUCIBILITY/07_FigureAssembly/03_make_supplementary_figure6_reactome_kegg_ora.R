options(stringsAsFactors = FALSE)
required <- c("ggplot2", "ragg", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

out_root <- "D:/workspace/ACLsenescence2_reviewer2_figure_reorganization"
out_dir <- file.path(out_root, "supplementary_figure6_reactome_kegg_ora")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(p) utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
export_ggplot <- function(p, base, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(base, ".svg"), width = w, height = h); print(p); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h); print(p); grDevices::dev.off()
  ragg::agg_tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); print(p); grDevices::dev.off()
  ragg::agg_png(paste0(base, ".png"), width = w, height = h, units = "in", res = 300); print(p); grDevices::dev.off()
}

all_results <- read_csv("D:/workspace/ACLsenescence2_reviewer3_persistent_ORA/persistent_ORA_all_Reactome_KEGG_results.csv")
sig <- all_results[all_results$gene_list %in% c("persistent_up", "persistent_down") &
                    all_results$primary_FDR_status == "FDR_supported", , drop = FALSE]
sig$primary_FDR <- as.numeric(sig$primary_FDR)
sig$enrichment_ratio <- as.numeric(sig$enrichment_ratio)
sig$overlap_size <- as.numeric(sig$overlap_size)
sig$direction_label <- ifelse(sig$gene_list == "persistent_up", "Persistent up", "Persistent down")
sig <- sig[order(sig$gene_list, sig$database, sig$primary_FDR), , drop = FALSE]
grp <- interaction(sig$gene_list, sig$database, drop = TRUE)
top <- do.call(rbind, lapply(split(sig, grp), function(z) utils::head(z, 8L)))
top$pathway_label <- as.character(top$pathway_name)
too_long <- nchar(top$pathway_label) > 58
top$pathway_label[too_long] <- paste0(substr(top$pathway_label[too_long], 1, 55), "...")
top$pathway_label <- factor(top$pathway_label, levels = rev(unique(top$pathway_label)))
top$neglog10_FDR <- -log10(top$primary_FDR)
top$database <- factor(top$database, levels = c("Reactome", "KEGG"))
top$direction_label <- factor(top$direction_label, levels = c("Persistent up", "Persistent down"))
top$panel <- factor(paste(top$direction_label, top$database, sep = " | "),
                    levels = c("Persistent up | Reactome", "Persistent up | KEGG",
                               "Persistent down | Reactome", "Persistent down | KEGG"))

p <- ggplot2::ggplot(top, ggplot2::aes(x = enrichment_ratio, y = pathway_label,
                                      size = overlap_size, colour = neglog10_FDR)) +
  ggplot2::geom_point(alpha = 0.9) +
  ggplot2::facet_wrap(~ panel, ncol = 2, scales = "free") +
  ggplot2::scale_colour_viridis_c(option = "C", name = expression(-log[10] * "(FDR)")) +
  ggplot2::scale_size_continuous(range = c(2.1, 6.6), name = "Overlapping genes") +
  ggplot2::labs(x = "Enrichment ratio (observed / expected)", y = NULL,
                title = "Broader pathway context for direction-consistent persistent genes") +
  ggplot2::theme_classic(base_size = 8.5) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
                 strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
                 strip.text = ggplot2::element_text(face = "bold", size = 8),
                 axis.text.y = ggplot2::element_text(size = 6.6),
                 axis.text.x = ggplot2::element_text(size = 7),
                 axis.title.x = ggplot2::element_text(size = 8.3),
                 panel.spacing = grid::unit(0.65, "lines"),
                 legend.title = ggplot2::element_text(size = 7.8),
                 legend.text = ggplot2::element_text(size = 7.2))
base <- file.path(out_dir, "SupplementaryFigure6_Reactome_KEGG_ORA")
export_ggplot(p, base, width_mm = 220, height_mm = 210)
utils::write.csv(top, file.path(out_dir, "SupplementaryFigure6_plot_source_top8_per_facet.csv"), row.names = FALSE)
utils::write.csv(sig, file.path(out_dir, "SupplementaryFigure6_all_FDR_supported_results.csv"), row.names = FALSE)
utils::write.csv(all_results, file.path(out_dir, "SupplementaryFigure6_all_tested_results.csv"), row.names = FALSE)
writeLines(c("Plot scope: up to eight pathways with the lowest primary FDR among FDR-supported pathways per database and direction; no significant stratum would be displayed without points.",
             "Supplementary Figure 6: Reactome/KEGG broader pathway ORA for persistent genes.",
             "The complete results for all tested Reactome/KEGG pathways, including nominal P values and BH-adjusted FDR values, are supplied as source data in SupplementaryFigure6_all_tested_results.csv.",
             "Primary ORA FDR: within each persistent-gene direction, BH jointly across 937 Reactome and 350 KEGG pathways.",
             "Persistent-up input: 970 genes; persistent-down input: 446 genes; the same database-annotated subset of the common 14,184-gene eligibility universe was used for both directions."),
           file.path(out_dir, "SupplementaryFigure6_plot_contract.txt"))
