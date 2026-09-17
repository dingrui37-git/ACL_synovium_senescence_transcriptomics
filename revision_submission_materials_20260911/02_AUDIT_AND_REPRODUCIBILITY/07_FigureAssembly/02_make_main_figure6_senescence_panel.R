options(stringsAsFactors = FALSE)
required <- c("ComplexHeatmap", "circlize", "ragg", "svglite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

out_root <- "D:/workspace/ACLsenescence2_reviewer2_figure_reorganization"
out_dir <- file.path(out_root, "figure6_senescence_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
read_csv <- function(p) utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)

export_heatmap <- function(ht, base_path, width_mm, height_mm, draw_args = list()) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  draw_one <- function() do.call(ComplexHeatmap::draw, c(list(object = ht), draw_args))
  svglite::svglite(paste0(base_path, ".svg"), width = w, height = h); draw_one(); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(base_path, ".pdf"), width = w, height = h); draw_one(); grDevices::dev.off()
  ragg::agg_tiff(paste0(base_path, ".tiff"), width = w, height = h, units = "in", res = 600, compression = "lzw"); draw_one(); grDevices::dev.off()
  ragg::agg_png(paste0(base_path, ".png"), width = w, height = h, units = "in", res = 300); draw_one(); grDevices::dev.off()
}

panel_dir <- "D:/workspace/ACLsenescence2_reviewer2_curated_senescence_panel"
panel8 <- read_csv(file.path(panel_dir, "curated_8set_panel_all_results.csv"))
hm_mouse <- read_csv("E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step13_Figure3B_mouse_Hallmark_GSEA/step13_fgsea_Hallmark_all.csv")
hm_early <- read_csv("E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/step20_current78_fgsea_full_hallmark_combined.csv")
hm_chronic <- read_csv("E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step4_current78_DE_GSEA/chronic_step4_fgsea_full_hallmark_ACLT_alone_52W_vs_Control_52W.csv")
if (!all(table(hm_mouse$timepoint) == 50L) ||
    !all(table(hm_early$contrast) == 50L) || nrow(hm_chronic) != 50L) {
  stop("The contextual readouts must come from complete 50-set Hallmark analyses.", call. = FALSE)
}

make_hallmark <- function(x, dataset, map_col = NULL, map_time = NULL, fixed_comparison = NULL) {
  y <- x[x$pathway %in% c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"), , drop = FALSE]
  if (!is.null(map_time)) y$comparison <- unname(map_time[as.character(y$timepoint)])
  else if (!is.null(map_col)) y$comparison <- unname(map_col[as.character(y$contrast)])
  else y$comparison <- fixed_comparison
  data.frame(dataset = dataset, comparison = as.character(y$comparison),
             pathway = as.character(y$pathway), NES = as.numeric(y$NES),
             pval = as.numeric(y$pval), FDR = as.numeric(y$padj), size = as.numeric(y$size),
             FDR_family = "BH across complete 50-set Hallmark collection within comparison",
             FDR_source = "Existing full Hallmark GSEA result",
             source = "Existing Hallmark GSEA", stringsAsFactors = FALSE)
}

hm <- rbind(
  make_hallmark(hm_mouse, "mouse_discovery", map_time = c("1W" = "mouse_1W", "4W" = "mouse_4W")),
  make_hallmark(hm_early, "pig_early", map_col = c(
    "ACLT_untreated_t7_vs_Control" = "pig_early_1W",
    "ACLT_untreated_t28_vs_Control" = "pig_early_4W"
  )),
  make_hallmark(hm_chronic, "pig_chronic", fixed_comparison = "pig_chronic_52W")
)

path8 <- data.frame(dataset = panel8$dataset, comparison = panel8$comparison,
                    pathway = panel8$pathway, NES = as.numeric(panel8$NES),
                    pval = as.numeric(panel8$pval), FDR = as.numeric(panel8$primary_FDR),
                    FDR_family = "BH across five direct and three mechanistic prespecified entries within comparison",
                    FDR_source = "Curated eight-set panel GSEA result",
                    size = as.numeric(panel8$fgsea_size),
                    source = "Prespecified senescence panel GSEA", stringsAsFactors = FALSE)
pathway_order <- c(
  "REACTOME_CELLULAR_SENESCENCE",
  "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
  "FRIDMAN_SENESCENCE_UP", "FRIDMAN_SENESCENCE_DN", "SenMayo",
  "HALLMARK_P53_PATHWAY", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)
path10 <- rbind(path8, hm)
path10$pathway <- as.character(path10$pathway)
path10 <- path10[path10$pathway %in% pathway_order, , drop = FALSE]
path10$pathway <- factor(path10$pathway, levels = pathway_order)
path10 <- path10[order(path10$comparison, path10$pathway), , drop = FALSE]
path10$pathway <- as.character(path10$pathway)
friendly <- c("Cellular Senescence", "Reactome SASP", "Fridman UP", "Fridman DN",
              "SenMayo", "p53", "E2F", "G2M",
              "Hallmark Inflammatory Response", "Hallmark EMT")
path10$display_name <- friendly[match(path10$pathway, pathway_order)]
path10$tier <- ifelse(path10$pathway %in% pathway_order[1:5], "Direct senescence/SASP",
                      ifelse(path10$pathway %in% pathway_order[6:8], "Supporting mechanism",
                             "Existing Hallmark context"))
if (any(!is.finite(path10$FDR))) stop("Missing FDR in Figure 6 pathway input.", call. = FALSE)

## Reconstruct mouse score from the same 75 mouse symbols used for the pig projection.
sig <- read_csv("E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv")
mouse_score_one <- function(tp) {
  voom <- read_csv(sprintf("E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/03_DE_analysis/step05_voom_logCPM_%s.csv", tp))
  gid <- as.character(voom[[1]])
  mat <- as.matrix(voom[, -1, drop = FALSE]); storage.mode(mat) <- "numeric"
  de <- read_csv(sprintf("E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/03_DE_analysis/step05_DE_%s_ACLR_vs_Contra_limma_voom.csv", tp))
  mapper <- stats::setNames(as.character(de$SYMBOL), as.character(de$gene_id))
  symbol <- unname(mapper[gid])
  keep <- !is.na(symbol) & symbol %in% sig$mouse_symbol
  map <- data.frame(gene_id = gid[keep], mouse_symbol = symbol[keep], stringsAsFactors = FALSE)
  map <- map[!duplicated(map$mouse_symbol), , drop = FALSE]
  map <- merge(map, sig[, c("mouse_symbol", "signature_direction")], by = "mouse_symbol", sort = FALSE)
  if (nrow(map) != 75L) stop("Mouse ", tp, " signature mapping has ", nrow(map), " genes; expected 75.")
  z <- t(scale(t(mat[match(map$gene_id, gid), , drop = FALSE]))); z[!is.finite(z)] <- NA_real_
  signs <- ifelse(map$signature_direction == "Up_in_ACLR", 1, -1)
  score <- colMeans(z * signs, na.rm = TRUE)
  control <- score[grepl("L$", names(score))]; case <- score[grepl("R$", names(score))]
  p <- stats::wilcox.test(case, control, paired = TRUE, exact = FALSE)$p.value
  data.frame(comparison = paste0("mouse_", tp),
             directional_score_delta_median = stats::median(case) - stats::median(control),
             median_control = stats::median(control), median_case = stats::median(case), pval = p,
              test = "paired Wilcoxon signed-rank", n_control = length(control), n_case = length(case),
              score_FDR_primary_score_analysis = NA_real_,
              score_FDR_family = "Descriptive mouse derivation reference; no inferential FDR",
              score_role = "Descriptive derivation-species reference score; 75 mouse genes corresponding to measurable one-to-one pig orthologues",
              stringsAsFactors = FALSE)
}

mouse_score <- rbind(mouse_score_one("1W"), mouse_score_one("4W"))
pig_early <- read_csv("E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_score_group_comparisons.csv")
pig_chronic <- read_csv("E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_score_group_comparisons.csv")
pig_early <- pig_early[pig_early$score == "directional_score", , drop = FALSE]
pig_chronic <- pig_chronic[pig_chronic$score == "directional_score", , drop = FALSE]
pig_score <- rbind(
  data.frame(comparison = c("pig_early_1W", "pig_early_4W"),
             directional_score_delta_median = pig_early$median_difference_case_minus_control,
             median_control = pig_early$median_control, median_case = pig_early$median_case,
             pval = pig_early$wilcox_p_value, test = "Wilcoxon rank-sum",
             n_control = pig_early$n_control, n_case = pig_early$n_case,
             score_FDR_primary_score_analysis = pig_early$BH_FDR_across_all_score_comparisons,
             score_FDR_family = "Original early pig score analysis: BH across 12 score-comparison rows (six score definitions by two time contrasts)",
             score_role = "Cross-species validation score; projected 75 pig orthologues"),
  data.frame(comparison = "pig_chronic_52W",
             directional_score_delta_median = pig_chronic$median_difference_case_minus_control,
             median_control = pig_chronic$median_control, median_case = pig_chronic$median_case,
             pval = pig_chronic$wilcox_p_value, test = "Wilcoxon rank-sum",
             n_control = pig_chronic$n_control, n_case = pig_chronic$n_case,
             score_FDR_primary_score_analysis = pig_chronic$BH_FDR_across_all_score_comparisons,
             score_FDR_family = "Original chronic pig score analysis: BH across six score definitions",
             score_role = "Cross-species validation score; projected 75 pig orthologues")
)
comparison_order <- c("mouse_1W", "mouse_4W", "pig_early_1W", "pig_early_4W", "pig_chronic_52W")
score_stats <- rbind(mouse_score, pig_score)
score_stats <- score_stats[match(comparison_order, score_stats$comparison), , drop = FALSE]
pig_score_idx <- grepl("^pig_", score_stats$comparison)
## Mouse scores reuse genes and directions derived from the mouse discovery data;
## they are descriptive reference values rather than independent validation tests.
## For continuity with the original Figure 4A/5A score analyses, pig score cells
## retain the dataset-specific primary-score BH-FDR already reported in those
## analyses rather than introducing a new three-contrast FDR family in Figure 6.
score_stats$score_inference_role <- ifelse(pig_score_idx,
                                           "Inferential pig validation",
                                           "Descriptive mouse derivation reference")
score_stats$display_comparison <- c("Mouse 1W", "Mouse 4W", "Pig 1W", "Pig 4W", "Pig 52W")

score_mat <- matrix(score_stats$directional_score_delta_median, nrow = 1,
                    dimnames = list("Directional score\n(mouse reference; pig validation)", score_stats$display_comparison))
score_fdr <- matrix(score_stats$score_FDR_primary_score_analysis, nrow = 1,
                    dimnames = dimnames(score_mat))
row_names <- friendly
path_mat <- matrix(NA_real_, nrow = length(row_names), ncol = length(comparison_order),
                   dimnames = list(row_names, score_stats$display_comparison))
path_fdr <- path_mat
path_family <- matrix(NA_character_, nrow = length(row_names), ncol = length(comparison_order),
                      dimnames = dimnames(path_mat))
for (i in seq_len(nrow(path10))) {
  rr <- match(path10$display_name[i], row_names)
  cc <- match(path10$comparison[i], comparison_order)
  path_mat[rr, cc] <- path10$NES[i]
  path_fdr[rr, cc] <- path10$FDR[i]
  path_family[rr, cc] <- path10$FDR_family[i]
}

utils::write.csv(score_stats, file.path(out_dir, "Figure6_directional_score_statistics.csv"), row.names = FALSE)
utils::write.csv(path10, file.path(out_dir, "Figure6_10set_GSEA_results_with_separate_FDR_families.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene_set = rownames(score_mat), score_mat, check.names = FALSE),
                 file.path(out_dir, "Figure6_directional_score_matrix.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene_set = rownames(path_mat), path_mat, check.names = FALSE),
                 file.path(out_dir, "Figure6_NES_matrix_numeric.csv"), row.names = FALSE)
disp_path <- matrix(paste0(sprintf("%.2f", path_mat), ifelse(path_fdr < 0.05, "*", "")),
                    nrow = nrow(path_mat), dimnames = dimnames(path_mat))
score_stars <- ifelse(is.finite(score_fdr) & score_fdr < 0.05, "*", "")
disp_score <- matrix(paste0(sprintf("%.2f", score_mat), score_stars),
                     nrow = 1, dimnames = dimnames(score_mat))
utils::write.csv(data.frame(gene_set = rownames(path_mat), disp_path, check.names = FALSE),
                 file.path(out_dir, "Figure6_NES_matrix_display_with_FDR_stars.csv"), row.names = FALSE)
utils::write.csv(data.frame(gene_set = rownames(score_mat), disp_score, check.names = FALSE),
                 file.path(out_dir, "Figure6_directional_score_matrix_display_with_FDR_stars.csv"), row.names = FALSE)

score_lim <- max(1, ceiling(max(abs(score_mat), na.rm = TRUE) * 10) / 10)
path_lim <- max(3.2, ceiling(max(abs(path_mat), na.rm = TRUE) * 10) / 10)
score_col <- circlize::colorRamp2(c(-score_lim, 0, score_lim), c("#2166AC", "#F7F7F7", "#B2182B"))
path_col <- circlize::colorRamp2(c(-path_lim, 0, path_lim), c("#2166AC", "#F7F7F7", "#B2182B"))
score_cell <- function(j, i, x, y, w, h, fill) {
  grid::grid.text(sprintf("%.2f", score_mat[i, j]), x = x, y = y,
                  gp = grid::gpar(fontsize = 9, col = ifelse(abs(score_mat[i, j]) > score_lim * .62, "white", "#222222")))
  if (is.finite(score_fdr[i, j]) && score_fdr[i, j] < .05)
    grid::grid.text("*", x = x + .38 * w, y = y + .30 * h,
                    gp = grid::gpar(fontsize = 10, fontface = "bold", col = ifelse(score_mat[i, j] > score_lim * .62, "white", "#222222")))
}
path_cell <- function(j, i, x, y, w, h, fill) {
  grid::grid.text(sprintf("%.2f", path_mat[i, j]), x = x, y = y,
                  gp = grid::gpar(fontsize = 7.8, col = ifelse(abs(path_mat[i, j]) > path_lim * .62, "white", "#222222")))
  if (is.finite(path_fdr[i, j]) && path_fdr[i, j] < .05)
    grid::grid.text("*", x = x + .38 * w, y = y + .30 * h,
                    gp = grid::gpar(fontsize = 8.7, fontface = "bold", col = ifelse(path_mat[i, j] > path_lim * .62, "white", "#222222")))
}

ht_score <- ComplexHeatmap::Heatmap(score_mat, name = "Score difference", col = score_col,
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 9, fontface = "bold"), cell_fun = score_cell,
  rect_gp = grid::gpar(col = "white", lwd = 1.1), heatmap_legend_param = list(title = "Score difference"))
ht_path <- ComplexHeatmap::Heatmap(path_mat, name = "NES", col = path_col,
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 8.5, fontface = "bold"),
  row_names_gp = grid::gpar(fontsize = 7.7), row_names_side = "left",
  row_split = factor(c(rep("Direct senescence/SASP", 5), rep("Supporting mechanism", 3),
                       rep("Existing Hallmark context", 2)),
                     levels = c("Direct senescence/SASP", "Supporting mechanism",
                                "Existing Hallmark context")),
  row_title_rot = 90, row_title_gp = grid::gpar(fontsize = 8.2, fontface = "bold"),
  cell_fun = path_cell, rect_gp = grid::gpar(col = "white", lwd = 1.1),
  heatmap_legend_param = list(title = "NES"))
ht_combined <- get("%v%", envir = asNamespace("ComplexHeatmap"))(ht_score, ht_path)
export_heatmap(ht_combined,
  file.path(out_dir, "Figure6_senescence_related_panel_NES_heatmap"),
  width_mm = 225, height_mm = 188,
  draw_args = list(heatmap_legend_side = "right", annotation_legend_side = "right",
                   merge_legends = FALSE, padding = grid::unit(c(6, 8, 6, 8), "mm")))

source_long <- rbind(
  data.frame(comparison = score_stats$comparison, comparison_label = score_stats$display_comparison,
             row_label = "Directional score (mouse reference; pig validation)", value_type = "Directional score median difference",
             value = score_stats$directional_score_delta_median,
             FDR = score_stats$score_FDR_primary_score_analysis,
             FDR_family = score_stats$score_FDR_family,
             inference_role = score_stats$score_inference_role,
             score_role = score_stats$score_role,
             source = "75-gene mouse reference score and pig projected validation score"),
  data.frame(comparison = rep(comparison_order, times = nrow(path_mat)),
             comparison_label = rep(colnames(path_mat), times = nrow(path_mat)),
             row_label = rep(rownames(path_mat), each = length(comparison_order)),
             value_type = "GSEA NES", value = as.vector(t(path_mat)), FDR = as.vector(t(path_fdr)),
             FDR_family = as.vector(t(path_family)),
             inference_role = NA_character_,
             score_role = NA_character_,
             source = "Eight-set senescence panel or complete 50-set Hallmark GSEA")
)
utils::write.csv(source_long, file.path(out_dir, "Figure6_source_data_long.csv"), row.names = FALSE)
writeLines(c(
  "Rows: directional score (mouse reference; pig validation), five direct senescence/SASP gene sets, three supporting mechanistic gene sets and two existing Hallmark context readouts.",
  "Columns: Mouse 1W; Mouse 4W; Pig 1W; Pig 4W; Pig 52W.",
  "GSEA FDR: the five direct and three mechanistic prespecified entries use BH within comparison across eight; inflammation and EMT retain FDR from the complete 50-set Hallmark analysis.",
  "Directional-score row: mouse is a descriptive derivation-species reference score based on the 75 mouse genes corresponding to the measurable one-to-one pig orthologues; pig is the projected validation score. Mouse scores use paired Wilcoxon tests but are not treated as independent validation tests because gene membership and directions were derived from mouse data.",
  "Pig directional-score cells retain the original dataset-specific primary score-analysis BH-FDR: early pig used all 12 score-comparison rows (six score definitions by two time contrasts), and chronic pig used all six score definitions. Mouse score cells have no inferential FDR or asterisk.",
  "Asterisk: FDR < 0.05 within the applicable inferential family. Color is centered at zero.",
  "The former Figure 5D information-integration panel was intentionally omitted."
), file.path(out_dir, "Figure6_contract_and_QA.txt"))
utils::write.csv(data.frame(metric = c("directional_score_rows", "GSEA_rows", "comparisons",
                                "n_GSEA_FDR_families", "directional_score_FDR_significant",
                                "mouse_score_inferential", "pig_score_FDR_family_definition",
                                "eight_set_FDR_significant", "complete_50_Hallmark_FDR_significant",
                                "source_data_rows"),
                      value = c(1, 10, 5, 2, sum(score_fdr < .05, na.rm = TRUE),
                                "No", "Early pig: 12 rows; chronic pig: 6 rows",
                                sum(path_fdr[seq_len(8), ] < .05, na.rm = TRUE),
                                sum(path_fdr[9:10, ] < .05, na.rm = TRUE), nrow(source_long))),
          file.path(out_dir, "Figure6_QA_summary.csv"), row.names = FALSE)
