# Step14: Figure 3C four-step bar plot based on current remapped mouse-to-human ortholog mapping and CellAge intersection
# Purpose:
# This step redraws Figure 3C using the current reproducible pipeline outputs and explicitly separates:
#   1) Step07 strict persistent mouse gene entries
#   2) Step07 persistent genes with valid mouse symbols for ortholog mapping
#   3) Current Step11B strict 1:1 human orthologs
#   4) Current Step11C CellAge-overlapping genes
#
# Rationale:
# Step07 can define persistent genes at the gene-entry level. Some entries may lack a usable mouse symbol
# and therefore cannot be submitted to gprofiler2::gorth for symbol-level mouse-to-human ortholog mapping.
# This script therefore keeps the official Step07 entry count as the first denominator and reports the
# mappable-symbol count as a separate second bar instead of labelling it as "persistent mouse genes".
#
# Current inputs:
#   Step07 persistent mouse genes:
#     07_tables/step07_strict_DEG_upset_persistent/step07_persistent_direction_consistent_genes.csv
#   Step11A cleaned CellAge file:
#     07_tables/step11A_CellAge_raw_to_clean/step11A_CellAge_clean_271903.csv
#   Step11B current remapped strict 1:1 mouse-human ortholog mapping:
#     07_tables/step11B_current_gorth_strict_1to1_remap/step11B_current_mouse_human_ortholog_mapping_strict_1to1.csv
#   Step11C current CellAge overlap:
#     07_tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_genes.csv
#
# This script does not use the old frozen Figure3C workspace.

options(stringsAsFactors = FALSE)

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

step07_dir <- file.path(base_dir, "07_tables", "step07_strict_DEG_upset_persistent")
step11A_dir <- file.path(base_dir, "07_tables", "step11A_CellAge_raw_to_clean")
step11B_dir <- file.path(base_dir, "07_tables", "step11B_current_gorth_strict_1to1_remap")
step11C_dir <- file.path(base_dir, "07_tables", "step11C_intersect_current_gorth_mapping_with_CellAge_clean")
step11D_dir <- file.path(base_dir, "07_tables", "step11D_audit_persistent_mapping_denominators")

table_dir <- file.path(base_dir, "07_tables", "step14_Figure3C_current_mapping_framework_barplot")
figure_dir <- file.path(base_dir, "06_figures", "Figure3")
log_dir <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

for (d in c(table_dir, figure_dir, log_dir, script_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

step_name <- "step14_Figure3C_current_mapping_framework_barplot"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

png_file <- file.path(figure_dir, "Figure3C_persistent_mouse_to_human_CellAge_framework_barplot_current_4step.png")
pdf_file <- file.path(figure_dir, "Figure3C_persistent_mouse_to_human_CellAge_framework_barplot_current_4step.pdf")

# Canonical manuscript copies. These will overwrite previous Figure3C canonical files.
png_file_canonical <- file.path(figure_dir, "Figure3C_persistent_mouse_to_human_CellAge_framework.png")
pdf_file_canonical <- file.path(figure_dir, "Figure3C_persistent_mouse_to_human_CellAge_framework.pdf")

# Clean only current Step14 outputs, not upstream files.
old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  log_file,
  file.path(figure_dir, "Figure3C_persistent_mouse_to_human_CellAge_framework_barplot_current.png"),
  file.path(figure_dir, "Figure3C_persistent_mouse_to_human_CellAge_framework_barplot_current.pdf"),
  png_file,
  pdf_file,
  png_file_canonical,
  pdf_file_canonical
)
unlink(old_outputs[file.exists(old_outputs)], recursive = TRUE, force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP14 FIGURE3C CURRENT MAPPING BARPLOT: FOUR-STEP DENOMINATOR-AWARE VERSION =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Base directory: ", base_dir, "\n\n", sep = "")

safe_library <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (cc in colnames(x)) {
    if (is.list(x[[cc]])) {
      x[[cc]] <- vapply(x[[cc]], function(v) paste(as.character(v), collapse = ";"), character(1))
    }
  }
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

find_existing_file <- function(candidates, label) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Missing ", label, ". Checked:\n", paste(candidates, collapse = "\n"))
  }
  normalizePath(hit[1], winslash = "/", mustWork = TRUE)
}

# MD5 is used only as a reproducibility fingerprint for the exact input file contents.
get_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  as.character(out)
}

file_audit <- function(path, label) {
  info <- file.info(path)
  data.frame(
    label = label,
    file = normalizePath(path, winslash = "/", mustWork = TRUE),
    size_bytes = as.numeric(info$size),
    mtime = as.character(info$mtime),
    ctime = as.character(info$ctime),
    md5 = get_md5(path),
    stringsAsFactors = FALSE
  )
}

find_col <- function(df, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    if (required) {
      stop("Could not find ", label, " column. Checked: ", paste(candidates, collapse = ", "),
           "\nAvailable columns: ", paste(colnames(df), collapse = ", "))
    }
    return(NA_character_)
  }
  hit[1]
}

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None", "none")] <- NA_character_
  x
}

safe_library("ggplot2")

# Current input files.
persistent_file <- find_existing_file(
  c(
    file.path(step07_dir, "step07_persistent_direction_consistent_genes.csv"),
    file.path(base_dir, "07_tables", "step07_persistent_direction_consistent_genes.csv")
  ),
  "Step07 persistent direction-consistent gene table"
)

persistent_summary_file <- find_existing_file(
  c(
    file.path(step07_dir, "step07_strict_DEG_persistent_summary.csv"),
    file.path(step07_dir, "step07_strict_DEG_summary.csv"),
    file.path(base_dir, "07_tables", "step07_strict_DEG_persistent_summary.csv")
  ),
  "Step07 persistent summary table"
)

cellage_clean_file <- find_existing_file(
  c(
    file.path(step11A_dir, "step11A_CellAge_clean_271903.csv"),
    file.path(step11A_dir, "step11A_CellAge_clean_symbol_list.csv")
  ),
  "Step11A cleaned CellAge table"
)

mapping_file <- find_existing_file(
  c(file.path(step11B_dir, "step11B_current_mouse_human_ortholog_mapping_strict_1to1.csv")),
  "Step11B current strict 1:1 mouse-human ortholog mapping table"
)

overlap_file <- find_existing_file(
  c(file.path(step11C_dir, "step11C_persistent_CellAge_overlap_genes.csv")),
  "Step11C current persistent-CellAge overlap table"
)

step11D_summary_file <- file.path(step11D_dir, "step11D_persistent_mapping_denominator_audit_summary.csv")
step11D_available <- file.exists(step11D_summary_file)

cat("Input files used:\n")
cat("Persistent genes: ", persistent_file, "\n", sep = "")
cat("Persistent summary: ", persistent_summary_file, "\n", sep = "")
cat("CellAge clean: ", cellage_clean_file, "\n", sep = "")
cat("Current strict 1:1 mapping: ", mapping_file, "\n", sep = "")
cat("Current CellAge overlap: ", overlap_file, "\n", sep = "")
cat("Step11D denominator audit available: ", step11D_available, "\n\n", sep = "")

persistent_df <- read.csv(persistent_file, check.names = FALSE, stringsAsFactors = FALSE)
persistent_summary_df <- read.csv(persistent_summary_file, check.names = FALSE, stringsAsFactors = FALSE)
cellage_clean <- read.csv(cellage_clean_file, check.names = FALSE, stringsAsFactors = FALSE)
mapping_df <- read.csv(mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
overlap_df <- read.csv(overlap_file, check.names = FALSE, stringsAsFactors = FALSE)
step11D_summary <- if (step11D_available) read.csv(step11D_summary_file, check.names = FALSE, stringsAsFactors = FALSE) else NULL

get_metric_value <- function(df, metric_name, numeric = TRUE, fallback = NA) {
  if (is.null(df) || !all(c("metric", "value") %in% colnames(df))) return(fallback)
  hit <- df$value[df$metric == metric_name]
  if (length(hit) == 0) return(fallback)
  hit <- hit[1]
  if (numeric) {
    suppressWarnings(as.numeric(hit))
  } else {
    as.character(hit)
  }
}

# Step07 official entry count: use table rows as the primary definition.
# If Step11D has already audited it, cross-check rather than silently relying on a label.
persistent_entry_total <- nrow(persistent_df)
step11D_persistent_rows <- get_metric_value(step11D_summary, "Step07_persistent_table_rows", numeric = TRUE, fallback = NA_real_)
if (!is.na(step11D_persistent_rows) && step11D_persistent_rows != persistent_entry_total) {
  warning("Step11D persistent row count differs from Step07 persistent table rows. ",
          "Step11D = ", step11D_persistent_rows, "; Step07 table rows = ", persistent_entry_total,
          ". Plot will use Step07 table rows.")
}

persistent_mouse_col <- find_col(
  persistent_df,
  c("input", "mouse_symbol", "SYMBOL", "symbol", "gene_symbol", "Gene symbol", "gene", "label_final"),
  "persistent mouse gene symbol"
)

mapping_mouse_col <- find_col(
  mapping_df,
  c("input", "mouse_symbol", "query", "Gene symbol", "gene_symbol"),
  "mapping mouse input"
)

mapping_human_col <- find_col(
  mapping_df,
  c("ortholog_name", "human_symbol", "name", "ortholog", "Gene symbol"),
  "mapping human ortholog"
)

cellage_symbol_col <- find_col(
  cellage_clean,
  c("Gene symbol", "Gene.Symbol", "gene_symbol", "CellAge_human_symbol", "symbol", "SYMBOL"),
  "CellAge human gene symbol"
)

overlap_human_col <- find_col(
  overlap_df,
  c("ortholog_name", "Gene symbol", "human_symbol", "CellAge_human_symbol", "gene_symbol"),
  "overlap human gene symbol"
)

persistent_symbols_all <- clean_symbol(persistent_df[[persistent_mouse_col]])
persistent_mappable_symbols <- unique(na.omit(persistent_symbols_all))
persistent_mappable_total <- length(persistent_mappable_symbols)
persistent_missing_symbol_rows <- sum(is.na(persistent_symbols_all))

mapping_mouse <- unique(na.omit(clean_symbol(mapping_df[[mapping_mouse_col]])))
mapping_human <- unique(na.omit(clean_symbol(mapping_df[[mapping_human_col]])))
cellage_symbols <- unique(na.omit(clean_symbol(cellage_clean[[cellage_symbol_col]])))
overlap_symbols <- unique(na.omit(clean_symbol(overlap_df[[overlap_human_col]])))

strict_ortholog_total <- length(mapping_human)
cellage_unique_total <- length(cellage_symbols)
overlap_total <- length(overlap_symbols)

# Row-level audit: join Step07 persistent entries to the current mapping by mouse symbol.
map_pair <- unique(data.frame(
  mouse_symbol_clean = clean_symbol(mapping_df[[mapping_mouse_col]]),
  human_symbol_clean = clean_symbol(mapping_df[[mapping_human_col]]),
  stringsAsFactors = FALSE
))
map_pair <- map_pair[!is.na(map_pair$mouse_symbol_clean) & !is.na(map_pair$human_symbol_clean), , drop = FALSE]

persistent_row_audit <- data.frame(
  step07_row_id = seq_len(nrow(persistent_df)),
  mouse_symbol_raw = as.character(persistent_df[[persistent_mouse_col]]),
  mouse_symbol_clean = persistent_symbols_all,
  stringsAsFactors = FALSE
)
persistent_row_audit <- merge(
  persistent_row_audit,
  map_pair,
  by = "mouse_symbol_clean",
  all.x = TRUE,
  sort = FALSE
)
persistent_row_audit$has_valid_mouse_symbol <- !is.na(persistent_row_audit$mouse_symbol_clean)
persistent_row_audit$has_strict_1to1_human_ortholog <- !is.na(persistent_row_audit$human_symbol_clean)
persistent_row_audit$is_CellAge_overlap <- persistent_row_audit$human_symbol_clean %in% overlap_symbols

row_level_strict_ortholog_total <- sum(persistent_row_audit$has_strict_1to1_human_ortholog, na.rm = TRUE)
row_level_overlap_total <- sum(persistent_row_audit$is_CellAge_overlap, na.rm = TRUE)
unique_mouse_symbols_with_overlap <- length(unique(na.omit(persistent_row_audit$mouse_symbol_clean[persistent_row_audit$is_CellAge_overlap])))

# Save rows without usable symbols: this explains why 1416 entries become 1409 mappable symbols.
missing_symbol_rows_df <- persistent_df[is.na(persistent_symbols_all), , drop = FALSE]
if (nrow(missing_symbol_rows_df) > 0) {
  missing_symbol_rows_df$step07_row_id <- which(is.na(persistent_symbols_all))
}
write_csv(missing_symbol_rows_df, file.path(table_dir, "step14_Step07_persistent_rows_without_mappable_symbol.csv"))
write_csv(persistent_row_audit, file.path(table_dir, "step14_Step07_persistent_rows_with_mapping_and_CellAge_status.csv"))

# Audit recomputed overlap from Step11A and Step11B.
recomputed_overlap_symbols <- intersect(mapping_human, cellage_symbols)
recomputed_overlap_n <- length(recomputed_overlap_symbols)
step11C_matches_recomputed <- identical(sort(overlap_symbols), sort(recomputed_overlap_symbols))
if (!step11C_matches_recomputed) {
  warning("Step11C overlap symbols differ from recomputed Step11A ∩ Step11B symbols. ",
          "Step11C = ", overlap_total, "; recomputed = ", recomputed_overlap_n,
          ". Plot will use Step11C overlap count.")
}

if (persistent_mappable_total != length(mapping_mouse)) {
  warning("Mappable persistent mouse symbols differ from Step11B mapping mouse inputs. ",
          "Mappable Step07 symbols = ", persistent_mappable_total,
          "; Step11B mapped mouse inputs = ", length(mapping_mouse),
          ". This may occur if some mappable symbols failed strict 1:1 mapping. This is expected, but should be recorded.")
}

# Four-step Figure3C data.
plot_df <- data.frame(
  category = c(
    "Strict persistent\nmouse gene entries",
    "Mappable mouse\nsymbols",
    "Strict 1:1 human\northologs",
    "CellAge\noverlap"
  ),
  n = c(
    persistent_entry_total,
    persistent_mappable_total,
    strict_ortholog_total,
    overlap_total
  ),
  denominator = persistent_entry_total,
  interpretation = c(
    "Step07 persistent direction-consistent gene entries",
    "Step07 persistent entries with non-missing mouse symbols available for ortholog mapping",
    "Unique strict one-to-one human ortholog symbols from current Step11B remapping",
    "Unique human ortholog symbols overlapping current Step11A cleaned CellAge symbols"
  ),
  stringsAsFactors = FALSE
)

plot_df$percent_of_persistent_entries <- plot_df$n / persistent_entry_total * 100
plot_df$label <- paste0(
  format(plot_df$n, big.mark = ","),
  "\n(",
  ifelse(abs(plot_df$percent_of_persistent_entries - 100) < 1e-8,
         "100",
         sprintf("%.1f", plot_df$percent_of_persistent_entries)),
  "%)"
)
plot_df$category <- factor(plot_df$category, levels = plot_df$category)

write_csv(plot_df, file.path(table_dir, "step14_Figure3C_current_4step_barplot_source_data.csv"))

input_audit <- rbind(
  file_audit(persistent_file, "Step07 persistent direction-consistent genes"),
  file_audit(persistent_summary_file, "Step07 persistent summary table"),
  file_audit(cellage_clean_file, "Step11A cleaned CellAge table"),
  file_audit(mapping_file, "Step11B current strict 1:1 mapping"),
  file_audit(overlap_file, "Step11C current CellAge overlap"),
  if (step11D_available) file_audit(step11D_summary_file, "Step11D denominator audit summary") else NULL
)
write_csv(input_audit, file.path(table_dir, "step14_Figure3C_current_4step_input_file_audit.csv"))

summary_df <- data.frame(
  metric = c(
    "persistent_input_file",
    "persistent_summary_input_file",
    "CellAge_clean_input_file",
    "current_mapping_input_file",
    "current_overlap_input_file",
    "Step11D_audit_input_file",
    "Step07_persistent_gene_entries",
    "Step07_persistent_rows_with_nonmissing_mouse_symbol",
    "Step07_persistent_rows_missing_mouse_symbol",
    "Step07_persistent_unique_mappable_mouse_symbols",
    "current_strict_1to1_mapping_rows",
    "current_strict_1to1_unique_mouse_inputs_with_strict_1to1_ortholog",
    "current_strict_1to1_unique_human_orthologs",
    "row_level_persistent_rows_with_strict_1to1_ortholog",
    "CellAge_clean_unique_symbols",
    "Step11C_unique_human_CellAge_overlap_symbols",
    "row_level_persistent_rows_with_CellAge_overlap",
    "unique_mouse_symbols_with_CellAge_overlap",
    "recomputed_overlap_from_Step11A_and_Step11B",
    "Step11C_matches_recomputed_overlap_symbols",
    "mappable_symbol_percent_of_persistent_entries",
    "strict_1to1_human_ortholog_percent_of_persistent_entries",
    "CellAge_overlap_percent_of_persistent_entries",
    "strict_1to1_human_ortholog_percent_of_mappable_symbols",
    "CellAge_overlap_percent_of_mappable_symbols",
    "workspace_used",
    "method_note",
    "plot_file_png",
    "plot_file_pdf",
    "canonical_plot_file_png",
    "canonical_plot_file_pdf"
  ),
  value = c(
    persistent_file,
    persistent_summary_file,
    cellage_clean_file,
    mapping_file,
    overlap_file,
    if (step11D_available) step11D_summary_file else "not available",
    persistent_entry_total,
    nrow(persistent_df) - persistent_missing_symbol_rows,
    persistent_missing_symbol_rows,
    persistent_mappable_total,
    nrow(mapping_df),
    length(mapping_mouse),
    strict_ortholog_total,
    row_level_strict_ortholog_total,
    cellage_unique_total,
    overlap_total,
    row_level_overlap_total,
    unique_mouse_symbols_with_overlap,
    recomputed_overlap_n,
    step11C_matches_recomputed,
    sprintf("%.1f", persistent_mappable_total / persistent_entry_total * 100),
    sprintf("%.1f", strict_ortholog_total / persistent_entry_total * 100),
    sprintf("%.1f", overlap_total / persistent_entry_total * 100),
    sprintf("%.1f", strict_ortholog_total / persistent_mappable_total * 100),
    sprintf("%.1f", overlap_total / persistent_mappable_total * 100),
    "FALSE",
    "Figure3C generated from current Step07/Step11A/Step11B/Step11C outputs; Step07 persistent entries and mappable symbols are shown separately; no frozen workspace was used.",
    png_file,
    pdf_file,
    png_file_canonical,
    pdf_file_canonical
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_df, file.path(table_dir, "step14_Figure3C_current_4step_barplot_run_summary.csv"))

p <- ggplot(plot_df, aes(x = category, y = n)) +
  geom_col(width = 0.68, fill = "#333333", color = "#333333") +
  geom_text(
    aes(label = label),
    vjust = -0.55,
    size = 5.4,
    lineheight = 1.05,
    color = "black"
  ) +
  scale_y_continuous(
    limits = c(0, max(plot_df$n) * 1.22),
    breaks = seq(0, ceiling(max(plot_df$n) / 500) * 500, by = 500),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Mapping of persistent mouse genes into the human CellAge framework",
    subtitle = "Persistent entries lacking mappable mouse symbols are shown before ortholog filtering",
    x = NULL,
    y = "Number of genes"
  ) +
  theme_bw(base_size = 17) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 21, margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 12.5, color = "grey35", margin = margin(b = 12)),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 13.2, color = "grey25", lineheight = 0.95),
    axis.text.y = element_text(size = 14.5, color = "grey25"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.45),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(t = 10, r = 18, b = 12, l = 12)
  )

ggsave(png_file, p, width = 12.2, height = 7.1, dpi = 320)
ggsave(pdf_file, p, width = 12.2, height = 7.1)
ggsave(png_file_canonical, p, width = 12.2, height = 7.1, dpi = 320)
ggsave(pdf_file_canonical, p, width = 12.2, height = 7.1)

versions <- data.frame(
  item = c("R", "ggplot2"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("ggplot2"))
  ),
  stringsAsFactors = FALSE
)
write_csv(versions, file.path(table_dir, "step14_Figure3C_current_4step_barplot_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step14_Figure3C_current_4step_barplot_sessionInfo.txt"))

# Archive this script if running by Rscript; otherwise save a marker.
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
  }
  NA_character_
}
script_path <- get_script_path()
archive_path <- file.path(script_dir, "step14_Figure3C_current_mapping_framework_barplot_4step.R")
if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path), archive_path)
} else {
  writeLines("# Script path not detected. Please manually save the executed Step14 four-step script here.", archive_path)
}
cat("Saved/checked script archive: ", archive_path, "\n", sep = "")

cat("\n===== STEP14 FIGURE3C CURRENT FOUR-STEP BARPLOT SUMMARY =====\n")
print(summary_df)

cat("\nInput file audit:\n")
print(input_audit)

cat("\nPlot source data:\n")
print(plot_df)

cat("\nRows without mappable mouse symbol from Step07 persistent table:\n")
print(missing_symbol_rows_df)

cat("\nFigure3C four-step barplot saved to:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
cat("\nCanonical copies saved to:\n")
cat(png_file_canonical, "\n")
cat(pdf_file_canonical, "\n")

cat("\nStep14 Figure3C current four-step mapping barplot completed successfully.\n")

sink()

cat("\nStep14 Figure3C current four-step mapping barplot completed. Please open:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
