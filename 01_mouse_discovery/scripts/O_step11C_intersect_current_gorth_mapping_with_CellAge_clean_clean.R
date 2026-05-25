# Step11C: Intersect current remapped strict 1:1 mouse-human ortholog mapping with cleaned CellAge table.
# English note:
# This script reconstructs the persistent ∩ CellAge human gene set for Figure 3C/3D
# using two reproducible source-data tables:
#   1) Step11A cleaned CellAge human gene table
#   2) Step11B current gorth-derived strict 1:1 mouse-to-human ortholog mapping table
#
# Locked method:
# - Do not rerun online gprofiler2::gorth in this step.
# - Use the current Step11B remapped strict 1:1 ortholog table as input.
# - Do not use any precomputed CellAge-overlap table as input.
# - Recompute the overlap by joining current remapped human ortholog symbols with cleaned CellAge symbols.
#
options(stringsAsFactors = FALSE)

base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery"

cellage_clean_file <- file.path(
  base_dir, "07_tables", "step11A_CellAge_raw_to_clean",
  "step11A_CellAge_clean_271903.csv"
)

current_mapping_file <- file.path(
  base_dir, "07_tables", "step11B_current_gorth_strict_1to1_remap",
  "step11B_current_mouse_human_ortholog_mapping_strict_1to1.csv"
)

table_dir <- file.path(base_dir, "07_tables", "step11C_intersect_current_gorth_mapping_with_CellAge_clean")
log_dir <- file.path(base_dir, "08_logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step11C_intersect_current_gorth_mapping_with_CellAge_clean"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  log_file
)
unlink(old_outputs[file.exists(old_outputs)], recursive = TRUE, force = TRUE)

sink(log_file, split = TRUE)

cat("===== STEP11C INTERSECT CURRENT GORTH MAPPING WITH CELLAGE CLEAN =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("CellAge clean input:\n", cellage_clean_file, "\n", sep = "")
cat("Current remapped mapping input:\n", current_mapping_file, "\n", sep = "")
cat("Overlap will be recomputed directly from current Step11B mapping and Step11A cleaned CellAge symbols.\n\n")

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

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "-", ".")] <- NA_character_
  x
}

if (!file.exists(cellage_clean_file)) {
  stop("Missing Step11A CellAge clean table: ", cellage_clean_file)
}
if (!file.exists(current_mapping_file)) {
  stop("Missing Step11B current remapped mapping table: ", current_mapping_file)
}

cellage_clean <- read.csv(cellage_clean_file, stringsAsFactors = FALSE, check.names = FALSE)
current_mapping <- read.csv(current_mapping_file, stringsAsFactors = FALSE, check.names = FALSE)

required_cellage_cols <- c("Gene symbol", "Senescence Effect", "CellAge_symbol")
missing_cellage <- setdiff(required_cellage_cols, colnames(cellage_clean))
if (length(missing_cellage) > 0) {
  stop("CellAge clean table is missing required columns: ", paste(missing_cellage, collapse = ", "))
}

required_mapping_cols <- c("input", "input_ensg", "ortholog_name", "ortholog_ensg", "n_human_orthologs", "n_mouse_inputs")
missing_mapping <- setdiff(required_mapping_cols, colnames(current_mapping))
if (length(missing_mapping) > 0) {
  stop("Current remapped mapping table is missing required columns: ", paste(missing_mapping, collapse = ", "))
}

cellage_clean$CellAge_symbol <- clean_symbol(cellage_clean$CellAge_symbol)
cellage_clean <- cellage_clean[!is.na(cellage_clean$CellAge_symbol) & cellage_clean$CellAge_symbol != "", , drop = FALSE]

## Enforce one CellAge row per symbol.
cellage_duplicate_symbols <- cellage_clean$CellAge_symbol[duplicated(cellage_clean$CellAge_symbol)]
if (length(cellage_duplicate_symbols) > 0) {
  dup <- cellage_clean[cellage_clean$CellAge_symbol %in% cellage_duplicate_symbols, , drop = FALSE]
  write_csv(dup, file.path(table_dir, "step11C_CellAge_clean_duplicate_symbols_detected.csv"))
  stop("CellAge clean table still contains duplicated CellAge_symbol values.")
}

current_mapping$input <- clean_symbol(current_mapping$input)
current_mapping$ortholog_name <- clean_symbol(current_mapping$ortholog_name)
current_mapping$n_human_orthologs <- as.numeric(current_mapping$n_human_orthologs)
current_mapping$n_mouse_inputs <- as.numeric(current_mapping$n_mouse_inputs)

mapping_bad <- current_mapping[
  is.na(current_mapping$input) | is.na(current_mapping$ortholog_name) |
    current_mapping$n_human_orthologs != 1 | current_mapping$n_mouse_inputs != 1,
  ,
  drop = FALSE
]

if (nrow(mapping_bad) > 0) {
  write_csv(mapping_bad, file.path(table_dir, "step11C_current_mapping_failed_validation_rows.csv"))
  stop("Current remapped mapping failed strict 1:1 validation.")
}

## Recompute overlap: human ortholog symbols x CellAge symbols.
overlap <- merge(
  current_mapping,
  cellage_clean,
  by.x = "ortholog_name",
  by.y = "CellAge_symbol",
  all = FALSE,
  sort = FALSE
)

overlap <- overlap[order(overlap$ortholog_name, overlap$input), , drop = FALSE]

## Save source inputs and recomputed overlap.
write_csv(cellage_clean, file.path(table_dir, "step11C_input_CellAge_clean_271903.csv"))
write_csv(current_mapping, file.path(table_dir, "step11C_input_current_mouse_human_ortholog_mapping_strict_1to1.csv"))
write_csv(overlap, file.path(table_dir, "step11C_persistent_CellAge_overlap_genes.csv"))

## Compact table for manuscript source data.
compact_cols <- intersect(
  c(
    "input", "input_ensg", "ortholog_name", "ortholog_ensg",
    "Entrez ID", "Gene symbol", "Gene name",
    "Cancer Cell", "Type of senescence", "Senescence Effect", "Reference",
    "n_human_orthologs", "n_mouse_inputs"
  ),
  colnames(overlap)
)
write_csv(overlap[, compact_cols, drop = FALSE], file.path(table_dir, "step11C_persistent_CellAge_overlap_genes_compact.csv"))

effect_counts <- as.data.frame(table(overlap[["Senescence Effect"]]), stringsAsFactors = FALSE)
colnames(effect_counts) <- c("Senescence_Effect", "n")
write_csv(effect_counts, file.path(table_dir, "step11C_overlap_Senescence_Effect_counts.csv"))

summary_df <- data.frame(
  metric = c(
    "CellAge_clean_input_file",
    "current_mapping_input_file",
    "CellAge_clean_unique_symbols",
    "current_strict_1to1_mapping_rows",
    "current_mapping_unique_mouse_inputs",
    "current_mapping_unique_human_orthologs",
    "recomputed_persistent_CellAge_overlap_genes",
    "main_overlap_output",
    "compact_overlap_output",
    "method_note"
  ),
  value = c(
    cellage_clean_file,
    current_mapping_file,
    length(unique(cellage_clean$CellAge_symbol)),
    nrow(current_mapping),
    length(unique(current_mapping$input)),
    length(unique(current_mapping$ortholog_name)),
    nrow(overlap),
    file.path(table_dir, "step11C_persistent_CellAge_overlap_genes.csv"),
    file.path(table_dir, "step11C_persistent_CellAge_overlap_genes_compact.csv"),
    "recomputed by intersecting Step11B current remapped strict 1:1 human ortholog symbols with Step11A cleaned CellAge symbols; no hard-coded expected overlap count was used"
  ),
  stringsAsFactors = FALSE
)

write_csv(summary_df, file.path(table_dir, "step11C_intersection_summary.csv"))

software_versions <- data.frame(
  item = c("R"),
  version = c(as.character(getRversion())),
  stringsAsFactors = FALSE
)
write_csv(software_versions, file.path(table_dir, "step11C_software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step11C_sessionInfo.txt"))

cat("\n===== STEP11C SUMMARY =====\n")
print(summary_df)

cat("\nSenescence Effect counts among recomputed overlap genes:\n")
print(effect_counts)

cat("\nFirst 20 recomputed overlap genes:\n")
print(head(overlap[, compact_cols, drop = FALSE], 20))

cat("\nStep11C completed successfully.\n")
cat("Main overlap output:\n", file.path(table_dir, "step11C_persistent_CellAge_overlap_genes.csv"), "\n", sep = "")

sink()

cat("\nStep11C completed.\n")
cat("Main overlap output:\n", file.path(table_dir, "step11C_persistent_CellAge_overlap_genes.csv"), "\n", sep = "")
