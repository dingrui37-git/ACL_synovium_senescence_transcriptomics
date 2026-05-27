# Step18_current78: Remap current mouse-derived 78 human CellAge-overlap genes to pig and validate in pig early QLF DE
# Purpose:
# Starting from the current 78 human CellAge-overlapping persistent genes from mouse discovery,
# this script performs a de novo human-to-pig gprofiler2::gorth mapping, keeps strict 1:1
# human-pig orthologs, checks whether pig orthologs are detected in the pig early gene count
# matrix, and evaluates direction consistency and strict validation in pig early t7/t28 QLF
# differential expression results. All outputs are saved under the pig early project folder
# for manuscript-ready reproducibility.

## =========================
## 0. Packages
## =========================

pkg_needed <- c("gprofiler2", "dplyr", "readr", "tools")
for (pkg in pkg_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
suppressPackageStartupMessages({
  library(gprofiler2)
  library(dplyr)
  library(readr)
  library(tools)
})

## =========================
## 1. Paths
## =========================

project_dir <- "E:/R/ACLsenescence2/rebuild_submission"
mouse_dir   <- file.path(project_dir, "02_mouse_discovery")
pig_dir     <- file.path(project_dir, "02_pig_early")

expected_current78_human_genes <- 78L

pick_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

mouse_current78_candidates <- c(
  file.path(
    mouse_dir,
    "07_tables",
    "step11C_intersect_current_gorth_mapping_with_CellAge_clean",
    "step11C_persistent_CellAge_overlap_genes.csv"
  ),
  file.path(
    mouse_dir,
    "07_tables",
    "step11C_intersect_current_gorth_mapping_with_CellAge_clean",
    "step11C_persistent_CellAge_overlap_genes_compact.csv"
  ),
  file.path(
    mouse_dir,
    "tables",
    "step11C_intersect_current_gorth_mapping_with_CellAge_clean",
    "step11C_persistent_CellAge_overlap_genes.csv"
  ),
  file.path(
    mouse_dir,
    "tables",
    "step11C_intersect_current_gorth_mapping_with_CellAge_clean",
    "step11C_persistent_CellAge_overlap_genes_compact.csv"
  )
)

mouse_current78_file <- pick_first_existing(mouse_current78_candidates)

mouse_persistent_file <- file.path(
  mouse_dir,
  "07_tables",
  "step07_strict_DEG_upset_persistent",
  "step07_persistent_direction_consistent_genes.csv"
)

pig_count_file <- file.path(pig_dir, "tables", "step16_pig_early_gene_count_matrix.csv")
pig_de_t7_file <- file.path(pig_dir, "tables", "step17_pig_early_DE_t7_vs_CON_QLF.csv")
pig_de_t28_file <- file.path(pig_dir, "tables", "step17_pig_early_DE_t28_vs_CON_QLF.csv")

fallback_pig_dir <- "E:/R/ACLsenescence2 LD/rebuild_submission/02_pig_early"
if (!file.exists(pig_count_file)) {
  pig_count_file <- file.path(fallback_pig_dir, "tables", "step16_pig_early_gene_count_matrix.csv")
}
if (!file.exists(pig_de_t7_file)) {
  pig_de_t7_file <- file.path(fallback_pig_dir, "tables", "step17_pig_early_DE_t7_vs_CON_QLF.csv")
}
if (!file.exists(pig_de_t28_file)) {
  pig_de_t28_file <- file.path(fallback_pig_dir, "tables", "step17_pig_early_DE_t28_vs_CON_QLF.csv")
}

out_dir    <- file.path(pig_dir, "tables", "step18_current78_pig_early_signature_remap")
log_dir    <- file.path(pig_dir, "logs")
script_dir <- file.path(pig_dir, "scripts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "step18_current78_pig_early_signature_remap_log.txt")

old_outputs <- list.files(out_dir, pattern = "^step18_current78_", full.names = TRUE)
if (length(old_outputs) > 0) file.remove(old_outputs)

sink(log_file, split = TRUE)
cat("===== STEP18_CURRENT78 PIG EARLY SIGNATURE REMAP AND VALIDATION =====\n")
cat("Run time:", format(Sys.time()), "\n")
cat("Project directory:", project_dir, "\n")
cat("Mouse current78 input:", mouse_current78_file, "\n")
cat("Expected current human CellAge-overlap genes:", expected_current78_human_genes, "\n")
cat("Mouse current78 candidate paths checked:\n")
cat(paste0("  - ", mouse_current78_candidates, collapse = "\n"), "\n", sep = "")
cat("Pig count matrix:", pig_count_file, "\n")
cat("Pig DE t7:", pig_de_t7_file, "\n")
cat("Pig DE t28:", pig_de_t28_file, "\n\n")

## =========================
## 2. Helper functions
## =========================

stop_if_missing <- function(path, label) {
  if (length(path) != 1 || is.na(path) || !file.exists(path)) {
    stop(label, " does not exist. Checked path: ", path, call. = FALSE)
  }
}

read_csv_keep <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null", "None")] <- NA
  x
}

find_col <- function(df, candidates, label) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop("Could not find ", label, ". Checked: ", paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  hit[1]
}

file_audit <- function(path, label) {
  info <- file.info(path)
  data.frame(
    label = label,
    file = normalizePath(path, winslash = "/", mustWork = FALSE),
    size_bytes = as.numeric(info$size),
    mtime = as.character(info$mtime),
    ctime = as.character(info$ctime),
    md5 = if (file.exists(path)) tools::md5sum(path)[[1]] else NA_character_,
    stringsAsFactors = FALSE
  )
}

safe_pkg_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
}

safe_gprofiler_version <- function(org) {
  out <- tryCatch(
    capture.output(print(gprofiler2::get_version_info(organism = org))),
    error = function(e) paste("get_version_info failed:", conditionMessage(e))
  )
  paste(out, collapse = " | ")
}

infer_direction_from_persistent <- function(df) {
  lower_names <- tolower(colnames(df))
  explicit_candidates <- colnames(df)[lower_names %in% c(
    "signature_direction", "persistent_direction", "mouse_direction", "direction", "direction_consistent"
  )]
  if (length(explicit_candidates) > 0) {
    raw <- as.character(df[[explicit_candidates[1]]])
    raw_low <- tolower(raw)
    dir <- ifelse(grepl("up|positive|higher", raw_low), "Up_in_ACLR",
                  ifelse(grepl("down|negative|lower", raw_low), "Down_in_ACLR", NA_character_))
    if (sum(!is.na(dir)) > 0) return(dir)
  }

  logfc_1w_candidates <- c("logFC_1W", "logFC.1W", "logFC_1w", "logFC_t1", "logFC_1")
  logfc_4w_candidates <- c("logFC_4W", "logFC.4W", "logFC_4w", "logFC_t4", "logFC_2")
  c1 <- logfc_1w_candidates[logfc_1w_candidates %in% colnames(df)]
  c4 <- logfc_4w_candidates[logfc_4w_candidates %in% colnames(df)]
  if (length(c1) > 0) {
    lf1 <- suppressWarnings(as.numeric(df[[c1[1]]]))
    dir <- ifelse(lf1 > 0, "Up_in_ACLR", ifelse(lf1 < 0, "Down_in_ACLR", NA_character_))
    if (length(c4) > 0) {
      lf4 <- suppressWarnings(as.numeric(df[[c4[1]]]))
      dir[!is.na(lf1) & !is.na(lf4) & sign(lf1) != sign(lf4)] <- NA_character_
    }
    return(dir)
  }
  stop("Could not infer mouse persistent direction. Need direction or logFC columns.", call. = FALSE)
}

get_de_core <- function(de_df, label) {
  gene_col <- find_col(de_df, c("gene_id", "pig_ensg", "GeneID", "genes", "gene", "ENSEMBL", "ensembl_gene_id"), paste0(label, " pig gene id column"))
  logfc_col <- find_col(de_df, c("logFC", "LogFC", "log2FC", "log2FoldChange"), paste0(label, " logFC column"))
  fdr_col <- find_col(de_df, c("FDR", "adj.P.Val", "adj.P.Val.", "padj", "adj_p_value", "adj.P.Value"), paste0(label, " FDR column"))
  p_candidates <- c("PValue", "P.Value", "P value", "P.Value.", "pvalue", "P")
  p_col <- p_candidates[p_candidates %in% colnames(de_df)]
  out <- data.frame(
    pig_ensg = as.character(de_df[[gene_col]]),
    logFC = suppressWarnings(as.numeric(de_df[[logfc_col]])),
    FDR = suppressWarnings(as.numeric(de_df[[fdr_col]])),
    stringsAsFactors = FALSE
  )
  if (length(p_col) > 0) out$PValue <- suppressWarnings(as.numeric(de_df[[p_col[1]]]))
  out
}

direction_match <- function(signature_direction, pig_logFC) {
  ifelse(
    is.na(signature_direction) | is.na(pig_logFC),
    NA,
    ifelse(signature_direction == "Up_in_ACLR", pig_logFC > 0,
           ifelse(signature_direction == "Down_in_ACLR", pig_logFC < 0, NA))
  )
}

## =========================
## 3. Check inputs and audit
## =========================

stop_if_missing(mouse_current78_file, "Current 78 human CellAge-overlap file")
stop_if_missing(mouse_persistent_file, "Mouse persistent file")
stop_if_missing(pig_count_file, "Pig early count matrix")
stop_if_missing(pig_de_t7_file, "Pig early t7 QLF DE table")
stop_if_missing(pig_de_t28_file, "Pig early t28 QLF DE table")

input_audit <- bind_rows(
  file_audit(mouse_current78_file, "Current mouse-derived 78 human CellAge-overlap genes"),
  file_audit(mouse_persistent_file, "Step07 mouse persistent direction-consistent genes"),
  file_audit(pig_count_file, "Pig early Step16 gene count matrix"),
  file_audit(pig_de_t7_file, "Pig early Step17 t7 vs CON QLF DE"),
  file_audit(pig_de_t28_file, "Pig early Step17 t28 vs CON QLF DE")
)
write.csv(input_audit, file.path(out_dir, "step18_current78_input_file_audit.csv"), row.names = FALSE)
cat("Input file audit:\n")
print(input_audit)

## =========================
## 4. Read current 78 human genes and mouse directions
## =========================

current78_df <- read_csv_keep(mouse_current78_file)
persistent_df <- read_csv_keep(mouse_persistent_file)

human_col <- find_col(
  current78_df,
  c("ortholog_name", "human_gene", "human_symbol", "Gene symbol", "Gene.symbol", "CellAge_human_symbol", "SYMBOL"),
  "human gene symbol column in current78 file"
)
mouse_col_current <- find_col(
  current78_df,
  c("input", "mouse_gene", "mouse_symbol", "mouse_input", "Mouse symbol", "SYMBOL_mouse"),
  "mouse source symbol column in current78 file"
)

current78_core <- current78_df %>%
  mutate(human_gene = clean_symbol(.data[[human_col]]), mouse_symbol = clean_symbol(.data[[mouse_col_current]])) %>%
  filter(!is.na(human_gene)) %>%
  distinct(human_gene, .keep_all = TRUE)

if (nrow(current78_core) != expected_current78_human_genes) {
  write.csv(
    current78_core,
    file.path(out_dir, "step18_current78_input_gene_count_failed_audit.csv"),
    row.names = FALSE
  )
  stop(
    "Current CellAge-overlap input gene count mismatch. Expected ",
    expected_current78_human_genes,
    " unique human genes, but found ",
    nrow(current78_core),
    ". This usually means an outdated Step11C file was selected or the input table columns were parsed incorrectly. ",
    "Audit saved to: ",
    file.path(out_dir, "step18_current78_input_gene_count_failed_audit.csv"),
    call. = FALSE
  )
}

mouse_col_persistent <- find_col(
  persistent_df,
  c("input", "mouse_gene", "mouse_symbol", "SYMBOL", "symbol", "gene_symbol", "label_final"),
  "mouse symbol column in Step07 persistent file"
)
persistent_df$mouse_symbol <- clean_symbol(persistent_df[[mouse_col_persistent]])
persistent_df$signature_direction <- infer_direction_from_persistent(persistent_df)

persistent_dir <- persistent_df %>%
  filter(!is.na(mouse_symbol), !is.na(signature_direction)) %>%
  select(mouse_symbol, signature_direction) %>%
  distinct(mouse_symbol, .keep_all = TRUE)

current78_annot <- current78_core %>% left_join(persistent_dir, by = "mouse_symbol")

if (any(is.na(current78_annot$signature_direction))) {
  missing_dir <- current78_annot %>% filter(is.na(signature_direction))
  write.csv(missing_dir, file.path(out_dir, "step18_current78_genes_missing_mouse_direction.csv"), row.names = FALSE)
  stop("Some current78 genes do not have mouse persistent direction. See step18_current78_genes_missing_mouse_direction.csv", call. = FALSE)
}

write.csv(current78_annot, file.path(out_dir, "step18_current78_human_genes_with_mouse_direction.csv"), row.names = FALSE)
cat("\nCurrent78 human genes with mouse direction:\n")
cat("n_current78_unique_human_genes =", nrow(current78_core), "\n")
cat("n_current78_with_direction =", nrow(current78_annot), "\n")
print(table(current78_annot$signature_direction))

## =========================
## 5. Human-to-pig gorth remapping
## =========================

human_query <- sort(unique(current78_annot$human_gene))
cat("\nRunning gprofiler2::gorth human -> pig...\n")
gorth_raw <- gprofiler2::gorth(
  query = human_query,
  source_organism = "hsapiens",
  target_organism = "sscrofa",
  mthreshold = Inf,
  filter_na = FALSE
)
write.csv(gorth_raw, file.path(out_dir, "step18_current78_human_to_pig_gorth_raw.csv"), row.names = FALSE)
if (nrow(gorth_raw) == 0) stop("gorth returned zero rows.", call. = FALSE)

ortholog_col <- find_col(gorth_raw, c("ortholog_name", "pig_symbol", "target"), "pig ortholog symbol column from gorth")
ortholog_ensg_col <- find_col(gorth_raw, c("ortholog_ensg", "pig_ensg", "target_ensg"), "pig ortholog Ensembl column from gorth")

gorth_counted <- gorth_raw %>%
  mutate(human_gene = clean_symbol(input), pig_symbol = clean_symbol(.data[[ortholog_col]]), pig_ensg = clean_symbol(.data[[ortholog_ensg_col]])) %>%
  filter(!is.na(human_gene), !is.na(pig_ensg)) %>%
  group_by(human_gene) %>% mutate(n_pig_orthologs = n_distinct(pig_ensg)) %>% ungroup() %>%
  group_by(pig_ensg) %>% mutate(n_human_inputs = n_distinct(human_gene)) %>% ungroup()

gorth_strict <- gorth_counted %>%
  filter(n_pig_orthologs == 1, n_human_inputs == 1) %>%
  distinct(human_gene, pig_ensg, .keep_all = TRUE)

write.csv(gorth_counted, file.path(out_dir, "step18_current78_human_to_pig_gorth_counted_pairs.csv"), row.names = FALSE)
write.csv(gorth_strict, file.path(out_dir, "step18_current78_human_to_pig_strict_1to1_pairs.csv"), row.names = FALSE)

## =========================
## 6. Check pig expression matrix presence
## =========================

count_df <- read_csv_keep(pig_count_file)
count_gene_col <- find_col(count_df, c("gene_id", "pig_ensg", "GeneID", "gene", "ENSEMBL"), "pig count matrix gene id column")
expr_gene_ids <- unique(as.character(count_df[[count_gene_col]]))

mapping_with_presence <- current78_annot %>%
  left_join(gorth_strict %>% select(human_gene, pig_ensg, pig_symbol, n_pig_orthologs, n_human_inputs), by = "human_gene") %>%
  mutate(
    mapped_to_pig_strict_1to1 = !is.na(pig_ensg),
    present_in_step16_counts = mapped_to_pig_strict_1to1 & pig_ensg %in% expr_gene_ids,
    included_in_current78_pig_signature = present_in_step16_counts & !is.na(signature_direction)
  )

write.csv(mapping_with_presence, file.path(out_dir, "step18_current78_human_to_pig_strict_1to1_with_presence.csv"), row.names = FALSE)

pig_signature <- mapping_with_presence %>%
  filter(included_in_current78_pig_signature) %>%
  distinct(pig_ensg, .keep_all = TRUE) %>%
  arrange(signature_direction, human_gene)

write.csv(pig_signature, file.path(out_dir, "step18_current78_pig_signature_gene_table.csv"), row.names = FALSE)
write.csv(pig_signature %>% select(human_gene, mouse_symbol, pig_ensg, pig_symbol, signature_direction), file.path(out_dir, "step18_current78_pig_signature_gene_list_all.csv"), row.names = FALSE)
write.csv(pig_signature %>% filter(signature_direction == "Up_in_ACLR"), file.path(out_dir, "step18_current78_pig_signature_gene_list_up.csv"), row.names = FALSE)
write.csv(pig_signature %>% filter(signature_direction == "Down_in_ACLR"), file.path(out_dir, "step18_current78_pig_signature_gene_list_down.csv"), row.names = FALSE)

## =========================
## 7. Pig early gene-level validation
## =========================

de_t7 <- get_de_core(read_csv_keep(pig_de_t7_file), "t7")
de_t28 <- get_de_core(read_csv_keep(pig_de_t28_file), "t28")

validation <- pig_signature %>%
  select(human_gene, mouse_symbol, pig_ensg, pig_symbol, signature_direction) %>%
  left_join(de_t7 %>% rename(logFC_t7 = logFC, FDR_t7 = FDR), by = "pig_ensg") %>%
  left_join(de_t28 %>% rename(logFC_t28 = logFC, FDR_t28 = FDR), by = "pig_ensg") %>%
  mutate(
    direction_consistent_t7 = direction_match(signature_direction, logFC_t7),
    direction_consistent_t28 = direction_match(signature_direction, logFC_t28),
    direction_consistent_both_timepoints = direction_consistent_t7 %in% TRUE & direction_consistent_t28 %in% TRUE,
    direction_consistent_any_timepoint = direction_consistent_t7 %in% TRUE | direction_consistent_t28 %in% TRUE,
    strict_t7 = direction_consistent_t7 %in% TRUE & !is.na(FDR_t7) & !is.na(logFC_t7) & FDR_t7 < 0.05 & abs(logFC_t7) > 1,
    strict_t28 = direction_consistent_t28 %in% TRUE & !is.na(FDR_t28) & !is.na(logFC_t28) & FDR_t28 < 0.05 & abs(logFC_t28) > 1,
    core_strict_both_t7_t28 = strict_t7 & strict_t28
  ) %>%
  arrange(desc(core_strict_both_t7_t28), desc(direction_consistent_both_timepoints), signature_direction, human_gene)

write.csv(validation, file.path(out_dir, "step18_current78_pig_signature_validation_table.csv"), row.names = FALSE)
core_genes <- validation %>% filter(core_strict_both_t7_t28)
write.csv(core_genes, file.path(out_dir, "step18_current78_pig_core_strict_both_t7_t28_genes.csv"), row.names = FALSE)

## =========================
## 8. Versions and summary
## =========================

version_df <- data.frame(
  item = c("R.version.string", "platform", "gprofiler2_package_version", "gprofiler2_base_url", "gprofiler2_hsapiens_version_info", "gprofiler2_sscrofa_version_info", "gorth_source_organism", "gorth_target_organism", "gorth_mthreshold", "gorth_filter_na", "strict_1to1_rule"),
  value = c(
    R.version.string,
    R.version$platform,
    safe_pkg_version("gprofiler2"),
    tryCatch(gprofiler2::get_base_url(), error = function(e) paste("get_base_url failed:", conditionMessage(e))),
    safe_gprofiler_version("hsapiens"),
    safe_gprofiler_version("sscrofa"),
    "hsapiens",
    "sscrofa",
    "Inf",
    "FALSE",
    "n_pig_orthologs == 1 and n_human_inputs == 1"
  ),
  stringsAsFactors = FALSE
)
write.csv(version_df, file.path(out_dir, "step18_current78_mapping_and_software_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "step18_current78_sessionInfo.txt"))

summary_df <- data.frame(
  metric = c(
    "current_human_CellAge_overlap_input_genes",
    "current_human_genes_with_mouse_direction",
    "gorth_raw_rows",
    "current_human_genes_strict_1to1_mapped_to_pig",
    "current_pig_orthologs_present_in_step16_counts",
    "current_pig_signature_genes_detected",
    "current_pig_signature_up",
    "current_pig_signature_down",
    "direction_consistent_both_timepoints",
    "direction_consistent_any_timepoint",
    "strict_t7",
    "strict_t28",
    "core_strict_both_t7_t28",
    "output_dir"
  ),
  value = c(
    length(unique(current78_annot$human_gene)),
    sum(!is.na(current78_annot$signature_direction)),
    nrow(gorth_raw),
    sum(mapping_with_presence$mapped_to_pig_strict_1to1),
    sum(mapping_with_presence$present_in_step16_counts),
    nrow(pig_signature),
    sum(pig_signature$signature_direction == "Up_in_ACLR"),
    sum(pig_signature$signature_direction == "Down_in_ACLR"),
    sum(validation$direction_consistent_both_timepoints, na.rm = TRUE),
    sum(validation$direction_consistent_any_timepoint, na.rm = TRUE),
    sum(validation$strict_t7, na.rm = TRUE),
    sum(validation$strict_t28, na.rm = TRUE),
    sum(validation$core_strict_both_t7_t28, na.rm = TRUE),
    normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "step18_current78_pig_early_signature_remap_summary.csv"), row.names = FALSE)

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
  NA_character_
}
script_path <- get_script_path()
archive_path <- file.path(script_dir, "step18_current78_pig_early_signature_remap.R")
if (!is.na(script_path) && file.exists(script_path)) {
  writeLines(readLines(script_path, warn = FALSE), archive_path)
} else {
  writeLines("# Script path not detected. Please manually save the executed script here for reproducibility.", archive_path)
}

cat("\n===== STEP18_CURRENT78 SUMMARY =====\n")
print(summary_df)
cat("\nSignature direction counts:\n")
print(table(pig_signature$signature_direction))
cat("\nCore genes preview:\n")
print(head(core_genes[, c("human_gene", "mouse_symbol", "pig_ensg", "pig_symbol", "signature_direction", "logFC_t7", "FDR_t7", "logFC_t28", "FDR_t28"), drop = FALSE], 20))
cat("\nVersion records:\n")
print(version_df[version_df$item %in% c("gprofiler2_package_version", "gprofiler2_base_url", "gorth_source_organism", "gorth_target_organism", "strict_1to1_rule"), ])
cat("\nStep18 current78 pig early signature remap completed successfully.\n")
cat("Main summary: ", file.path(out_dir, "step18_current78_pig_early_signature_remap_summary.csv"), "\n", sep = "")

sink()
cat("\nStep18 current78 pig early signature remap completed.\n")
cat("Main summary: ", file.path(out_dir, "step18_current78_pig_early_signature_remap_summary.csv"), "\n", sep = "")
