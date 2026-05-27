# StepSINGSENS_01: current78-derived early pig singscore sensitivity analysis.
# This script computes rank-based singscore sensitivity scores using the full genome-wide
# TMM-normalized logCPM matrix regenerated from the raw gene-level count matrix, without
# filterByExpr before rankGenes. It uses the locked current78-derived pig early signature
# expectation: 75 detected pig signature genes, including 65 up-signature genes and 10
# down-signature genes. Outputs include score tables, Wilcoxon group comparisons,
# correlations against the main z-score based signature scores, audit tables, methods text,
# and a summary log for reproducibility.

options(stringsAsFactors = FALSE)

method_version <- "2026-05-10_FULL_LOGCPM_v2_robust_signature_column_detection_no_filtered_rank_background"
step_name <- "StepSINGSENS_01_pig_early_current78_singscore_sensitivity"

input_count_file <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv"
signature_table_file <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv"
main_scores_file <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv"
main_signature_logcpm_file <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_logCPM_matrix.csv"

output_root <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/singscore_sensitivity"

dirs <- file.path(output_root, c("tables", "source_data", "figures", "logs", "objects", "scripts"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
log_file <- file.path(output_root, "logs", paste0(step_name, "_summary_log.txt"))

run_status <- "FAILED"

clean_previous_outputs <- function(root_dir) {
  if (!dir.exists(root_dir)) return(integer(0))
  targets <- list.files(root_dir, pattern = "^StepSINGSENS_01", recursive = TRUE, full.names = TRUE)
  targets <- targets[!grepl("summary_log\\.txt$", targets)]
  if (length(targets) > 0) {
    ok <- file.remove(targets)
    return(sum(ok))
  }
  0L
}

sink(log_file, split = TRUE)
on.exit({
  cat("\n============================================================\n")
  cat(step_name, "finished with status:", run_status, "\n")
  cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Summary log saved to:\n", log_file, "\n", sep = "")
  cat("============================================================\n")
  sink()
}, add = TRUE)

cat("============================================================\n")
cat(step_name, "\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Method version:", method_version, "\n")
cat("============================================================\n\n")

cat("Cleaning previous StepSINGSENS_01 outputs from output_root subfolders...\n")
removed_n <- clean_previous_outputs(output_root)
cat("Number of previous StepSINGSENS_01 files removed:", removed_n, "\n\n")

required_pkgs <- c("edgeR", "limma", "singscore", "ggplot2", "dplyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
       ". Please install them before rerunning this script.")
}

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(singscore)
  library(ggplot2)
  library(dplyr)
})

cat("Required packages loaded.\n")
cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
cat("limma version:", as.character(packageVersion("limma")), "\n")
cat("singscore version:", as.character(packageVersion("singscore")), "\n")
cat("ggplot2 version:", as.character(packageVersion("ggplot2")), "\n")
cat("dplyr version:", as.character(packageVersion("dplyr")), "\n\n")

for (f in c(input_count_file, signature_table_file, main_scores_file, main_signature_logcpm_file)) {
  if (!file.exists(f)) stop("Required input file does not exist: ", f)
}

sample_meta <- data.frame(
  sample_id = c(paste0("CON", 1:6), paste0("INJS", 1:6), paste0("INJL", 1:6)),
  group = factor(c(rep("Control", 6), rep("ACLT_1W", 6), rep("ACLT_4W", 6)),
                 levels = c("Control", "ACLT_1W", "ACLT_4W")),
  time_point = c(rep("t0", 6), rep("t7", 6), rep("t28", 6)),
  sample_number = rep(1:6, 3),
  stringsAsFactors = FALSE
)

read_csv_keep_names <- function(file) {
  read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
}

make_unique_rownames <- function(x) {
  make.unique(as.character(x), sep = "__dup")
}

# Infer the best gene identifier column by overlap with a reference gene universe.
infer_gene_column <- function(df, reference_genes = NULL, preferred_patterns = c("^gene_id$", "pig.*gene", "pig.*ensembl", "ensembl", "gene")) {
  if (ncol(df) == 0) return(NA_character_)
  cn <- colnames(df)
  non_numeric_cols <- cn[!vapply(df, is.numeric, logical(1))]
  candidate_cols <- unique(c(
    cn[Reduce(`|`, lapply(preferred_patterns, function(p) grepl(p, cn, ignore.case = TRUE)))],
    non_numeric_cols,
    cn[1]
  ))
  candidate_cols <- candidate_cols[candidate_cols %in% cn]
  if (length(candidate_cols) == 0) return(cn[1])
  if (!is.null(reference_genes)) {
    overlap_n <- vapply(candidate_cols, function(cc) {
      vals <- unique(as.character(df[[cc]]))
      sum(vals %in% reference_genes, na.rm = TRUE)
    }, integer(1))
    if (max(overlap_n, na.rm = TRUE) > 0) return(candidate_cols[which.max(overlap_n)])
  }
  candidate_cols[1]
}

standardize_direction <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[^a-z0-9_+\\-]", "_", y)
  out <- rep(NA_character_, length(y))
  out[grepl("down|dn|negative|repress|decreas|minus|\\-1|^-|low", y)] <- "down"
  out[grepl("up|positive|induc|increas|plus|\\+1|^\\+|high", y)] <- "up"
  # If a string contains both, use stricter exact terms where possible.
  both <- which(grepl("up", y) & grepl("down", y))
  if (length(both) > 0) out[both] <- NA_character_
  out
}

infer_direction_column <- function(df) {
  cn <- colnames(df)
  preferred <- cn[grepl("direction|status|regulation|mouse|persistent|signature|class", cn, ignore.case = TRUE)]
  candidate_cols <- unique(c(preferred, cn[!vapply(df, is.numeric, logical(1))]))
  candidate_cols <- candidate_cols[candidate_cols %in% cn]
  if (length(candidate_cols) == 0) return(NA_character_)
  score <- vapply(candidate_cols, function(cc) {
    z <- standardize_direction(df[[cc]])
    sum(!is.na(z)) + 100L * as.integer(any(z == "up", na.rm = TRUE) & any(z == "down", na.rm = TRUE))
  }, integer(1))
  if (max(score, na.rm = TRUE) <= 0) return(NA_character_)
  candidate_cols[which.max(score)]
}

read_matrix_with_gene_ids <- function(file, reference_genes = NULL) {
  df <- read_csv_keep_names(file)
  gene_col <- infer_gene_column(df, reference_genes = reference_genes)
  if (is.na(gene_col) || !(gene_col %in% colnames(df))) {
    stop("Could not infer gene identifier column for matrix file: ", file)
  }
  gene_ids <- as.character(df[[gene_col]])
  data_cols <- setdiff(colnames(df), gene_col)
  mat_df <- df[, data_cols, drop = FALSE]
  numeric_ok <- vapply(mat_df, function(z) all(suppressWarnings(!is.na(as.numeric(z))) | is.na(z)), logical(1))
  mat_df <- mat_df[, numeric_ok, drop = FALSE]
  mat <- as.matrix(data.frame(lapply(mat_df, as.numeric), check.names = FALSE))
  colnames(mat) <- colnames(mat_df)
  rownames(mat) <- make_unique_rownames(gene_ids)
  list(matrix = mat, gene_col = gene_col, raw_df = df)
}

find_score_column <- function(df, patterns) {
  cn <- colnames(df)
  for (p in patterns) {
    hit <- cn[grepl(p, cn, ignore.case = TRUE)]
    hit <- hit[vapply(df[hit], is.numeric, logical(1))]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

extract_total_score <- function(score_df, expected_samples) {
  score_df <- as.data.frame(score_df, stringsAsFactors = FALSE)
  score_col <- find_score_column(score_df, c("^TotalScore$", "total.*score", "score"))
  if (is.na(score_col)) {
    numeric_cols <- colnames(score_df)[vapply(score_df, is.numeric, logical(1))]
    if (length(numeric_cols) == 0) stop("Could not find a numeric score column in singscore output.")
    score_col <- numeric_cols[1]
  }
  scores <- as.numeric(score_df[[score_col]])
  sample_ids <- rownames(score_df)
  if (is.null(sample_ids) || any(is.na(sample_ids)) || any(sample_ids == "") || length(sample_ids) != length(scores)) {
    sample_ids <- expected_samples
  }
  if (length(sample_ids) != length(scores)) {
    stop("singscore output sample number does not match expected sample number.")
  }
  data.frame(sample_id = sample_ids, score = scores, stringsAsFactors = FALSE)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

safe_wilcox <- function(score_df, score_col, case_group) {
  control <- score_df[[score_col]][score_df$group == "Control"]
  case <- score_df[[score_col]][score_df$group == case_group]
  out <- data.frame(
    score = score_col,
    comparison = paste0(case_group, " vs Control"),
    n_control = length(control),
    n_case = length(case),
    median_control = ifelse(length(control) > 0, median(control, na.rm = TRUE), NA_real_),
    median_case = ifelse(length(case) > 0, median(case, na.rm = TRUE), NA_real_),
    median_case_minus_control = ifelse(length(control) > 0 && length(case) > 0,
                                       median(case, na.rm = TRUE) - median(control, na.rm = TRUE), NA_real_),
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )
  if (length(control) > 0 && length(case) > 0) {
    out$p_value <- suppressWarnings(stats::wilcox.test(case, control, exact = FALSE)$p.value)
  }
  out
}

cat("Input files:\n")
cat("Count matrix:", input_count_file, "\n")
cat("Signature table:", signature_table_file, "\n")
cat("Main z-score scores:", main_scores_file, "\n")
cat("Main signature logCPM matrix:", main_signature_logcpm_file, "\n\n")

# -----------------------------------------------------------------------------
# 1) Regenerate full genome-wide TMM logCPM matrix without filterByExpr.
# -----------------------------------------------------------------------------
count_df <- read_csv_keep_names(input_count_file)
cat("Raw count table dimension:", nrow(count_df), "rows x", ncol(count_df), "columns\n")
cat("First 10 columns:\n")
print(head(colnames(count_df), 10))

count_gene_col <- infer_gene_column(count_df, reference_genes = NULL, preferred_patterns = c("^gene_id$", "gene"))
if (is.na(count_gene_col)) stop("Could not infer gene_id column in count matrix.")
count_gene_ids <- as.character(count_df[[count_gene_col]])
count_sample_cols <- setdiff(colnames(count_df), count_gene_col)
count_sample_cols <- intersect(sample_meta$sample_id, count_sample_cols)
if (length(count_sample_cols) != 18) {
  stop("Expected 18 core sample columns in count matrix, but found: ", length(count_sample_cols),
       ". Found columns: ", paste(count_sample_cols, collapse = ", "))
}

counts <- as.matrix(count_df[, count_sample_cols, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- make_unique_rownames(count_gene_ids)
# Current count matrix should not contain duplicate gene IDs, but keep a safety check.
if (any(grepl("__dup", rownames(counts), fixed = TRUE))) {
  cat("WARNING: Duplicate gene IDs were detected and made unique for matrix operations.\n")
}
counts <- counts[, sample_meta$sample_id, drop = FALSE]
cat("Validated count matrix dimension:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Detect and remove duplicated artificial rows only if duplicated rownames exist after make.unique? no collapse.
# edgeR accepts unique row names.
dge <- edgeR::DGEList(counts = counts, group = sample_meta$group)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
full_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

cat("Full genome-wide TMM logCPM matrix generated.\n")
cat("filterByExpr applied before rankGenes: FALSE\n")
cat("full_logCPM dimension:", nrow(full_logcpm), "genes x", ncol(full_logcpm), "samples\n\n")

write.csv(data.frame(gene_id = rownames(full_logcpm), full_logcpm, check.names = FALSE),
          file.path(output_root, "source_data", "StepSINGSENS_01_full_genomewide_TMM_logCPM_matrix_no_filterByExpr.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 2) Read current78 signature table and infer pig gene IDs/directions robustly.
# -----------------------------------------------------------------------------
sig_df <- read_csv_keep_names(signature_table_file)
cat("Signature table dimension:", nrow(sig_df), "rows x", ncol(sig_df), "columns\n")
cat("Signature table column names:\n")
print(colnames(sig_df))

sig_gene_col <- infer_gene_column(sig_df, reference_genes = rownames(full_logcpm))
sig_direction_col <- infer_direction_column(sig_df)
cat("Inferred signature gene ID column:", sig_gene_col, "\n")
cat("Inferred signature direction column:", sig_direction_col, "\n")

if (is.na(sig_gene_col) || !(sig_gene_col %in% colnames(sig_df))) {
  stop("Could not infer pig gene ID column from signature table. Please inspect StepSINGSENS_01_signature_table_column_audit.csv.")
}
if (is.na(sig_direction_col) || !(sig_direction_col %in% colnames(sig_df))) {
  stop("Could not infer up/down direction column from signature table. Please inspect StepSINGSENS_01_signature_table_column_audit.csv.")
}

column_audit <- data.frame(
  column = colnames(sig_df),
  is_numeric = vapply(sig_df, is.numeric, logical(1)),
  overlap_with_full_logcpm_gene_ids = vapply(colnames(sig_df), function(cc) sum(unique(as.character(sig_df[[cc]])) %in% rownames(full_logcpm), na.rm = TRUE), integer(1)),
  standardized_direction_nonNA = vapply(colnames(sig_df), function(cc) sum(!is.na(standardize_direction(sig_df[[cc]]))), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(column_audit,
          file.path(output_root, "tables", "StepSINGSENS_01_signature_table_column_audit.csv"),
          row.names = FALSE)

sig_work <- data.frame(
  row_index = seq_len(nrow(sig_df)),
  pig_gene_id = as.character(sig_df[[sig_gene_col]]),
  signature_direction = standardize_direction(sig_df[[sig_direction_col]]),
  raw_direction_value = as.character(sig_df[[sig_direction_col]]),
  stringsAsFactors = FALSE
)
sig_work <- sig_work[!is.na(sig_work$pig_gene_id) & sig_work$pig_gene_id != "", , drop = FALSE]
sig_work$detected_in_full_logcpm <- sig_work$pig_gene_id %in% rownames(full_logcpm)

# Deduplicate by pig_gene_id while preserving direction. If duplicated with same direction, keep one.
dup_summary <- sig_work %>%
  group_by(pig_gene_id) %>%
  summarise(
    n_rows = n(),
    directions = paste(sort(unique(na.omit(signature_direction))), collapse = ";"),
    detected_in_full_logcpm = any(detected_in_full_logcpm),
    .groups = "drop"
  )
write.csv(dup_summary,
          file.path(output_root, "tables", "StepSINGSENS_01_signature_gene_duplicate_audit.csv"),
          row.names = FALSE)

conflicting <- dup_summary[grepl(";", dup_summary$directions), , drop = FALSE]
if (nrow(conflicting) > 0) {
  write.csv(conflicting,
            file.path(output_root, "tables", "StepSINGSENS_01_signature_gene_direction_conflicts.csv"),
            row.names = FALSE)
  stop("Conflicting up/down directions were detected for duplicated pig_gene_id values. See StepSINGSENS_01_signature_gene_direction_conflicts.csv.")
}

sig_unique <- sig_work %>%
  filter(!is.na(signature_direction), signature_direction %in% c("up", "down")) %>%
  arrange(row_index) %>%
  distinct(pig_gene_id, .keep_all = TRUE)

signature_detection_audit <- sig_unique %>%
  mutate(
    detected_in_count_matrix = pig_gene_id %in% rownames(counts),
    detected_in_full_logcpm = pig_gene_id %in% rownames(full_logcpm)
  ) %>%
  arrange(signature_direction, pig_gene_id)
write.csv(signature_detection_audit,
          file.path(output_root, "tables", "StepSINGSENS_01_current78_signature_detection_audit.csv"),
          row.names = FALSE)

up_genes <- signature_detection_audit$pig_gene_id[signature_detection_audit$signature_direction == "up" & signature_detection_audit$detected_in_full_logcpm]
down_genes <- signature_detection_audit$pig_gene_id[signature_detection_audit$signature_direction == "down" & signature_detection_audit$detected_in_full_logcpm]
all_detected_genes <- c(up_genes, down_genes)

observed_total <- length(all_detected_genes)
observed_up <- length(up_genes)
observed_down <- length(down_genes)

cat("\nCurrent78 signature detection from full genome-wide TMM logCPM matrix:\n")
cat("Detected total:", observed_total, "\n")
cat("Detected up:", observed_up, "\n")
cat("Detected down:", observed_down, "\n\n")

if (!(observed_total == 75 && observed_up == 65 && observed_down == 10)) {
  final_bad <- data.frame(
    metric = c("run_status", "reason", "observed_total", "observed_up", "observed_down",
               "expected_total", "expected_up", "expected_down", "signature_gene_col", "signature_direction_col"),
    value = c("FAILED", "Detected current78 signature counts do not match locked expectation.",
              observed_total, observed_up, observed_down, 75, 65, 10, sig_gene_col, sig_direction_col),
    stringsAsFactors = FALSE
  )
  write.csv(final_bad,
            file.path(output_root, "tables", "StepSINGSENS_01_failed_signature_count_summary.csv"),
            row.names = FALSE)
  stop("The detected current78 signature counts do not match the locked expectation of 75 total, 65 up, 10 down. Observed: total=",
       observed_total, ", up=", observed_up, ", down=", observed_down,
       ". Audit tables were written under: ", file.path(output_root, "tables"))
}

write.csv(data.frame(gene_id = up_genes, signature_direction = "up"),
          file.path(output_root, "source_data", "StepSINGSENS_01_current78_up_genes_used.csv"),
          row.names = FALSE)
write.csv(data.frame(gene_id = down_genes, signature_direction = "down"),
          file.path(output_root, "source_data", "StepSINGSENS_01_current78_down_genes_used.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 3) Audit main z-score signature matrices without failing when structure differs.
# -----------------------------------------------------------------------------
main_sig_matrix_info <- tryCatch({
  tmp <- read_matrix_with_gene_ids(main_signature_logcpm_file, reference_genes = rownames(full_logcpm))
  main_sig_genes <- rownames(tmp$matrix)
  if (length(main_sig_genes) == 0) main_sig_genes <- character(0)
  audit <- data.frame(
    gene_id = main_sig_genes,
    in_current78_signature_table = main_sig_genes %in% signature_detection_audit$pig_gene_id,
    in_full_logcpm = main_sig_genes %in% rownames(full_logcpm),
    stringsAsFactors = FALSE
  )
  write.csv(audit,
            file.path(output_root, "tables", "StepSINGSENS_01_main_signature_logCPM_gene_audit.csv"),
            row.names = FALSE)
  list(gene_col = tmp$gene_col, n_genes = length(main_sig_genes), status = "parsed")
}, error = function(e) {
  writeLines(conditionMessage(e),
             con = file.path(output_root, "logs", "StepSINGSENS_01_main_signature_logCPM_parse_warning.txt"))
  list(gene_col = NA_character_, n_genes = NA_integer_, status = paste0("not_parsed: ", conditionMessage(e)))
})
cat("Main signature logCPM audit status:", main_sig_matrix_info$status, "\n")
cat("Main signature logCPM inferred gene column:", main_sig_matrix_info$gene_col, "\n")
cat("Main signature logCPM genes parsed:", main_sig_matrix_info$n_genes, "\n\n")

# -----------------------------------------------------------------------------
# 4) singscore rankGenes and simpleScore.
# -----------------------------------------------------------------------------
cat("Running singscore::rankGenes on full genome-wide TMM logCPM matrix...\n")
ranked <- singscore::rankGenes(full_logcpm)
saveRDS(ranked, file.path(output_root, "objects", "StepSINGSENS_01_singscore_ranked_genes_full_logCPM.rds"))
cat("rankGenes completed.\n\n")

score_directional_raw <- singscore::simpleScore(ranked, upSet = up_genes, downSet = down_genes,
                                                centerScore = TRUE, knownDirection = TRUE)
score_up_raw <- singscore::simpleScore(ranked, upSet = up_genes,
                                       centerScore = TRUE, knownDirection = TRUE)
score_down_as_up_raw <- singscore::simpleScore(ranked, upSet = down_genes,
                                               centerScore = TRUE, knownDirection = TRUE)

directional_df <- extract_total_score(score_directional_raw, expected_samples = colnames(full_logcpm))
up_df <- extract_total_score(score_up_raw, expected_samples = colnames(full_logcpm))
down_raw_df <- extract_total_score(score_down_as_up_raw, expected_samples = colnames(full_logcpm))

singscore_by_sample <- sample_meta %>%
  left_join(directional_df %>% rename(directional_singscore = score), by = "sample_id") %>%
  left_join(up_df %>% rename(up_singscore = score), by = "sample_id") %>%
  left_join(down_raw_df %>% rename(down_raw_singscore = score), by = "sample_id") %>%
  mutate(down_reoriented_singscore = -down_raw_singscore)

write.csv(singscore_by_sample,
          file.path(output_root, "tables", "StepSINGSENS_01_current78_singscore_scores_by_sample.csv"),
          row.names = FALSE)
cat("singscore scores by sample:\n")
print(singscore_by_sample)
cat("\n")

# -----------------------------------------------------------------------------
# 5) Compare with main z-score based signature scores.
# -----------------------------------------------------------------------------
main_scores <- read_csv_keep_names(main_scores_file)
cat("Main z-score score table columns:\n")
print(colnames(main_scores))

main_sample_col <- colnames(main_scores)[tolower(colnames(main_scores)) %in% c("sample_id", "sample", "sampleid")]
if (length(main_sample_col) == 0) {
  # choose column with maximal overlap with sample IDs
  overlap <- vapply(colnames(main_scores), function(cc) sum(as.character(main_scores[[cc]]) %in% sample_meta$sample_id), integer(1))
  if (max(overlap) > 0) main_sample_col <- colnames(main_scores)[which.max(overlap)] else stop("Could not infer sample_id column in main z-score scores file.")
} else {
  main_sample_col <- main_sample_col[1]
}

main_scores$sample_id_for_join <- as.character(main_scores[[main_sample_col]])
main_directional_col <- find_score_column(main_scores, c("^directional_score$", "directional.*score"))
main_up_col <- find_score_column(main_scores, c("^up_score$", "up.*score"))
main_down_reoriented_col <- find_score_column(main_scores, c("down.*reorient.*score", "down.*reoriented", "down_reoriented_score"))

cat("Inferred main score sample column:", main_sample_col, "\n")
cat("Inferred main directional score column:", main_directional_col, "\n")
cat("Inferred main up score column:", main_up_col, "\n")
cat("Inferred main down reoriented score column:", main_down_reoriented_col, "\n\n")

merged_scores <- singscore_by_sample %>%
  left_join(main_scores, by = c("sample_id" = "sample_id_for_join"), suffix = c("", "_main"))
write.csv(merged_scores,
          file.path(output_root, "tables", "StepSINGSENS_01_current78_singscore_merged_with_main_zscore_scores.csv"),
          row.names = FALSE)

cor_rows <- list()
if (!is.na(main_directional_col)) {
  cor_rows[[length(cor_rows) + 1]] <- data.frame(
    comparison = "main_zscore_directional_score_vs_directional_singscore",
    main_score_column = main_directional_col,
    singscore_column = "directional_singscore",
    spearman_rho = safe_cor(as.numeric(merged_scores[[main_directional_col]]), merged_scores$directional_singscore),
    n_samples = sum(is.finite(as.numeric(merged_scores[[main_directional_col]])) & is.finite(merged_scores$directional_singscore)),
    stringsAsFactors = FALSE
  )
}
if (!is.na(main_up_col)) {
  cor_rows[[length(cor_rows) + 1]] <- data.frame(
    comparison = "main_zscore_up_score_vs_up_singscore",
    main_score_column = main_up_col,
    singscore_column = "up_singscore",
    spearman_rho = safe_cor(as.numeric(merged_scores[[main_up_col]]), merged_scores$up_singscore),
    n_samples = sum(is.finite(as.numeric(merged_scores[[main_up_col]])) & is.finite(merged_scores$up_singscore)),
    stringsAsFactors = FALSE
  )
}
if (!is.na(main_down_reoriented_col)) {
  cor_rows[[length(cor_rows) + 1]] <- data.frame(
    comparison = "main_zscore_down_reoriented_score_vs_down_reoriented_singscore",
    main_score_column = main_down_reoriented_col,
    singscore_column = "down_reoriented_singscore",
    spearman_rho = safe_cor(as.numeric(merged_scores[[main_down_reoriented_col]]), merged_scores$down_reoriented_singscore),
    n_samples = sum(is.finite(as.numeric(merged_scores[[main_down_reoriented_col]])) & is.finite(merged_scores$down_reoriented_singscore)),
    stringsAsFactors = FALSE
  )
}
correlation_summary <- if (length(cor_rows) > 0) bind_rows(cor_rows) else data.frame(
  comparison = character(0), main_score_column = character(0), singscore_column = character(0),
  spearman_rho = numeric(0), n_samples = integer(0)
)
write.csv(correlation_summary,
          file.path(output_root, "tables", "StepSINGSENS_01_current78_singscore_vs_main_zscore_correlation.csv"),
          row.names = FALSE)
cat("Correlation summary against main z-score based scores:\n")
print(correlation_summary)
cat("\n")

# -----------------------------------------------------------------------------
# 6) Wilcoxon group comparisons for rank-based singscore scores.
# -----------------------------------------------------------------------------
score_cols <- c("directional_singscore", "up_singscore", "down_reoriented_singscore")
wilcox_results <- bind_rows(lapply(score_cols, function(sc) {
  bind_rows(
    safe_wilcox(singscore_by_sample, sc, "ACLT_1W"),
    safe_wilcox(singscore_by_sample, sc, "ACLT_4W")
  )
})) %>%
  mutate(p_adj_BH_within_singscore_sensitivity = p.adjust(p_value, method = "BH"))
write.csv(wilcox_results,
          file.path(output_root, "tables", "StepSINGSENS_01_current78_singscore_group_comparison_wilcox.csv"),
          row.names = FALSE)
cat("Wilcoxon group comparison results:\n")
print(wilcox_results)
cat("\n")

# -----------------------------------------------------------------------------
# 7) Simple figures for audit/submission support.
# -----------------------------------------------------------------------------
plot_df <- singscore_by_sample
p_dir <- ggplot(plot_df, aes(x = group, y = directional_singscore)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, size = 2.5) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Directional singscore", title = "Current78-derived early pig signature")
ggsave(file.path(output_root, "figures", "StepSINGSENS_01_directional_singscore_by_group.pdf"),
       p_dir, width = 4.6, height = 3.8, useDingbats = FALSE)
ggsave(file.path(output_root, "figures", "StepSINGSENS_01_directional_singscore_by_group.png"),
       p_dir, width = 4.6, height = 3.8, dpi = 450)

if (!is.na(main_directional_col)) {
  scatter_df <- merged_scores
  p_cor <- ggplot(scatter_df, aes(x = .data[[main_directional_col]], y = directional_singscore, label = sample_id)) +
    geom_point(size = 2.5) +
    theme_classic(base_size = 11) +
    labs(x = "Main z-score directional score", y = "Directional singscore",
         title = "z-score score vs rank-based singscore")
  ggsave(file.path(output_root, "figures", "StepSINGSENS_01_directional_zscore_vs_singscore_scatter.pdf"),
         p_cor, width = 4.6, height = 3.8, useDingbats = FALSE)
  ggsave(file.path(output_root, "figures", "StepSINGSENS_01_directional_zscore_vs_singscore_scatter.png"),
         p_cor, width = 4.6, height = 3.8, dpi = 450)
}

# -----------------------------------------------------------------------------
# 8) Method parameters and methods text.
# -----------------------------------------------------------------------------
method_parameters <- data.frame(
  parameter = c(
    "method_version",
    "input_count_file",
    "signature_table_file",
    "main_scores_file",
    "main_signature_logcpm_file",
    "output_root",
    "signature_gene_id_column_inferred",
    "signature_direction_column_inferred",
    "rankGenes_input",
    "filterByExpr_applied_before_rankGenes",
    "normalization_method",
    "logCPM_prior_count",
    "singscore_rank_function",
    "singscore_score_function",
    "centerScore",
    "knownDirection",
    "directional_singscore_definition",
    "down_reoriented_singscore_definition",
    "expected_signature_total",
    "expected_signature_up",
    "expected_signature_down"
  ),
  value = c(
    method_version,
    input_count_file,
    signature_table_file,
    main_scores_file,
    main_signature_logcpm_file,
    output_root,
    sig_gene_col,
    sig_direction_col,
    "full genome-wide TMM-normalized logCPM matrix regenerated from raw count matrix",
    "FALSE",
    "edgeR::calcNormFactors(method='TMM')",
    "1",
    "singscore::rankGenes",
    "singscore::simpleScore",
    "TRUE",
    "TRUE",
    "simpleScore(upSet=65 up genes, downSet=10 down genes, centerScore=TRUE, knownDirection=TRUE)",
    "-simpleScore(upSet=10 down genes, centerScore=TRUE, knownDirection=TRUE)",
    "75",
    "65",
    "10"
  ),
  stringsAsFactors = FALSE
)
write.csv(method_parameters,
          file.path(output_root, "tables", "StepSINGSENS_01_method_parameters.csv"),
          row.names = FALSE)

methods_text <- c(
  "Singscore sensitivity analysis was performed for the current78-derived early pig signature using a full genome-wide TMM-normalized logCPM matrix regenerated from the raw gene-level count matrix.",
  "No filterByExpr filtering was applied before singscore::rankGenes, so that the rank background represented the full measured gene-level expression space rather than the PCA/MDS-filtered expression matrix.",
  "The locked detected pig signature comprised 75 genes, including 65 up-signature genes and 10 down-signature genes.",
  "Gene ranks were calculated using singscore::rankGenes, and scores were calculated using singscore::simpleScore with centerScore = TRUE and knownDirection = TRUE.",
  "The directional singscore was calculated using the up-signature genes as upSet and the down-signature genes as downSet.",
  "The down-signature score was additionally calculated by scoring down-signature genes as an upSet and multiplying the resulting score by -1 to obtain a biologically reoriented down-signature score.",
  "Two-sided Wilcoxon rank-sum tests were used for targeted comparisons of ACLT_1W or ACLT_4W against Control. Spearman correlations were used to compare rank-based singscore outputs with the main z-score-based signature scores."
)
writeLines(methods_text,
           con = file.path(output_root, "tables", "StepSINGSENS_01_methods_text_current78_singscore_sensitivity.txt"))

# -----------------------------------------------------------------------------
# 9) Final summary.
# -----------------------------------------------------------------------------
get_cor_value <- function(pattern) {
  if (nrow(correlation_summary) == 0) return(NA_real_)
  hit <- grepl(pattern, correlation_summary$comparison)
  if (!any(hit)) return(NA_real_)
  correlation_summary$spearman_rho[which(hit)[1]]
}
get_wilcox_p <- function(score, comp) {
  hit <- wilcox_results$score == score & wilcox_results$comparison == comp
  if (!any(hit)) return(NA_real_)
  wilcox_results$p_value[which(hit)[1]]
}

final_summary <- data.frame(
  metric = c(
    "run_status",
    "method_version",
    "rankGenes_input",
    "filterByExpr_applied_before_rankGenes",
    "signature_detected_total",
    "signature_detected_up",
    "signature_detected_down",
    "samples",
    "groups",
    "directional_zscore_vs_singscore_spearman_rho",
    "up_zscore_vs_singscore_spearman_rho",
    "down_reoriented_zscore_vs_singscore_spearman_rho",
    "directional_singscore_ACLT_1W_vs_Control_p",
    "directional_singscore_ACLT_4W_vs_Control_p",
    "main_signature_logCPM_audit_status",
    "main_signature_logCPM_genes_parsed"
  ),
  value = c(
    "SUCCESS",
    method_version,
    "full genome-wide TMM-normalized logCPM matrix regenerated from raw count matrix",
    "FALSE",
    observed_total,
    observed_up,
    observed_down,
    nrow(sample_meta),
    paste(paste(names(table(sample_meta$group)), as.integer(table(sample_meta$group))), collapse = "; "),
    get_cor_value("directional"),
    get_cor_value("up_score"),
    get_cor_value("down_reoriented"),
    get_wilcox_p("directional_singscore", "ACLT_1W vs Control"),
    get_wilcox_p("directional_singscore", "ACLT_4W vs Control"),
    main_sig_matrix_info$status,
    main_sig_matrix_info$n_genes
  ),
  stringsAsFactors = FALSE
)
write.csv(final_summary,
          file.path(output_root, "tables", "StepSINGSENS_01_final_summary_for_review.csv"),
          row.names = FALSE)

cat("Final summary for review:\n")
print(final_summary)
cat("\n")

cat("Key output files:\n")
cat("1) Scores by sample:\n")
cat(file.path(output_root, "tables", "StepSINGSENS_01_current78_singscore_scores_by_sample.csv"), "\n")
cat("2) Correlation with main z-score scores:\n")
cat(file.path(output_root, "tables", "StepSINGSENS_01_current78_singscore_vs_main_zscore_correlation.csv"), "\n")
cat("3) Wilcoxon group comparisons:\n")
cat(file.path(output_root, "tables", "StepSINGSENS_01_current78_singscore_group_comparison_wilcox.csv"), "\n")
cat("4) Signature detection audit:\n")
cat(file.path(output_root, "tables", "StepSINGSENS_01_current78_signature_detection_audit.csv"), "\n")
cat("5) Full genome-wide TMM logCPM matrix used for rankGenes:\n")
cat(file.path(output_root, "source_data", "StepSINGSENS_01_full_genomewide_TMM_logCPM_matrix_no_filterByExpr.csv"), "\n")
cat("6) Final summary:\n")
cat(file.path(output_root, "tables", "StepSINGSENS_01_final_summary_for_review.csv"), "\n\n")

cat("Session information:\n")
print(sessionInfo())

run_status <- "SUCCESS"
cat("\n============================================================\n")
cat(step_name, "completed successfully.\n")
cat("Finished at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
