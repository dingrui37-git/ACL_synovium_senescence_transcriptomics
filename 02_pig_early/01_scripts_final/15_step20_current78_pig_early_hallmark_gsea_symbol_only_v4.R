# ============================================================
# IMMEDIATE CONSOLE CHECKPOINT
# If you see this line, the correct Step20 v4 script is being sourced.
# ============================================================
cat("\n>>> STEP20 SYMBOL-ONLY V4 SCRIPT STARTED SUCCESSFULLY <<<\n")
cat("Script run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")

# ============================================================
# Step20_current78: Pig early Hallmark GSEA using full Sus scrofa Hallmark collection
# Purpose:
#   Starting from the Step17 pig early edgeR QLF DE tables, this step builds
#   preranked gene lists for ACLT untreated t7 vs Control and ACLT untreated
#   t28 vs Control, runs fgsea against the full MSigDB Hallmark collection for
#   Sus scrofa, and then extracts the two prespecified target pathways:
#     - HALLMARK_INFLAMMATORY_RESPONSE
#     - HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
#
# Method note:
#   This current version runs fgsea on the full Hallmark collection rather than
#   only on the two target pathways, so adjusted P values are retained from the
#   full Hallmark-level analysis.
#
# Important:
#   This is the official replacement version for Step20; it keeps the original analysis logic and explicitly records fgsea seed and nproc settings for reproducibility.
#   It does NOT auto-archive itself and does NOT write an empty script placeholder,
#   to avoid overwriting manually saved R scripts.
#   Corrected identifier policy: formal fgsea ranked vectors use pig gene symbols;
#   unresolved missing symbols are excluded rather than replaced with gene_id fallback.
#   V4 additionally fixes zero-row audit table handling.
# ============================================================

## =========================
## 0. User paths
## =========================
base_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_pig_early"

tables_dir <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures", "Figure4B_current78")
out_dir <- file.path(tables_dir, "step20_current78_pig_early_hallmark_gsea")
logs_dir <- file.path(out_dir, "logs")
objects_dir <- file.path(out_dir, "objects")

for (d in c(out_dir, logs_dir, objects_dir, figures_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## Close any leftover normal sink from a previous interrupted Step20 run.
## Do not touch message sinks here. Restart R if message output looks redirected.
while (sink.number() > 0) {
  sink(NULL)
}

full_log_file <- file.path(logs_dir, "step20_current78_pig_early_hallmark_gsea_full_log.txt")
summary_log_file <- file.path(logs_dir, "step20_current78_pig_early_hallmark_gsea_summary_to_send_me.txt")

## Fix accidental line break typo guard if copied manually
if (!exists("summary_log_file")) {
  summary_log_file <- file.path(logs_dir, "step20_current78_pig_early_hallmark_gsea_summary_to_send_me.txt")
}

sink(full_log_file, split = TRUE)
cat("===== STEP20_CURRENT78 PIG EARLY HALLMARK GSEA =====\n")
cat("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Base directory:", base_dir, "\n")
cat("Method: full Sus scrofa Hallmark collection fgsea; target pathways extracted after full run.\n")

## =========================
## 1. Packages
## =========================
required_pkgs <- c("fgsea", "msigdbr", "ggplot2", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}

## =========================
## 2. Inputs
## =========================
t7_file <- file.path(tables_dir, "step17_pig_early_DE_t7_vs_CON_QLF.csv")
t28_file <- file.path(tables_dir, "step17_pig_early_DE_t28_vs_CON_QLF.csv")

if (!file.exists(t7_file)) stop("Missing t7 QLF DE table: ", t7_file)
if (!file.exists(t28_file)) stop("Missing t28 QLF DE table: ", t28_file)

cat("\nInput DE tables:\n")
cat("t7 :", t7_file, "\n")
cat("t28:", t28_file, "\n")

de_t7 <- read.csv(t7_file, stringsAsFactors = FALSE, check.names = FALSE)
de_t28 <- read.csv(t28_file, stringsAsFactors = FALSE, check.names = FALSE)

## =========================
## 3. Helper functions
## =========================
find_col <- function(df, candidates, required = TRUE, label = "column") {
  hits <- candidates[candidates %in% colnames(df)]
  if (length(hits) > 0) return(hits[1])
  if (required) {
    stop("Could not find ", label, ". Checked: ", paste(candidates, collapse = ", "))
  }
  NA_character_
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

safe_p <- function(p) {
  p <- as.numeric(p)
  p[is.na(p)] <- NA_real_
  p[p <= 0] <- .Machine$double.xmin
  p
}

build_gene_id_symbol_map_from_available_de_tables <- function(de_list) {
  map_parts <- list()
  for (nm in names(de_list)) {
    df <- de_list[[nm]]
    gene_id_col <- find_col(df, c("gene_id", "pig_ensg", "GeneID", "Geneid", "feature_id"), TRUE, paste0(nm, " gene id column"))
    symbol_col <- find_col(df, c("gene_name", "gene_symbol", "SYMBOL", "symbol", "pig_symbol", "external_gene_name"), FALSE, paste0(nm, " gene symbol column"))
    if (is.na(symbol_col)) next
    tmp <- data.frame(
      source_table = nm,
      gene_id = as.character(df[[gene_id_col]]),
      gene_symbol = clean_symbol(df[[symbol_col]]),
      stringsAsFactors = FALSE
    )
    tmp <- tmp[!is.na(tmp$gene_id) & nzchar(tmp$gene_id) &
                 !is.na(tmp$gene_symbol) & nzchar(tmp$gene_symbol), , drop = FALSE]
    if (nrow(tmp) > 0) map_parts[[nm]] <- tmp
  }

  if (length(map_parts) == 0) {
    empty_audit <- data.frame(
      source_table = character(),
      gene_id = character(),
      gene_symbol = character(),
      stringsAsFactors = FALSE
    )
    return(list(map = setNames(character(0), character(0)), audit = empty_audit))
  }

  audit <- unique(do.call(rbind, map_parts))
  # Keep only unambiguous gene_id -> one symbol mappings.
  n_symbols <- stats::aggregate(gene_symbol ~ gene_id, data = audit, FUN = function(x) length(unique(x)))
  one_to_one_ids <- n_symbols$gene_id[n_symbols$gene_symbol == 1]
  audit_unique <- audit[audit$gene_id %in% one_to_one_ids, c("gene_id", "gene_symbol"), drop = FALSE]
  audit_unique <- unique(audit_unique)
  map <- setNames(audit_unique$gene_symbol, audit_unique$gene_id)

  list(map = map, audit = audit)
}

rank_one_table <- function(df, contrast_label, gene_id_to_symbol_map = NULL) {
  gene_id_col <- find_col(df, c("gene_id", "pig_ensg", "GeneID", "Geneid", "feature_id"), TRUE, "gene id column")
  symbol_col <- find_col(df, c("gene_name", "gene_symbol", "SYMBOL", "symbol", "pig_symbol", "external_gene_name"), FALSE, "gene symbol column")
  logfc_col <- find_col(df, c("logFC", "log2FC", "log2FoldChange"), TRUE, "logFC column")
  p_col <- find_col(df, c("PValue", "P.Value", "pvalue", "P.value", "P"), TRUE, "P value column")
  fdr_col <- find_col(df, c("FDR", "adj.P.Val", "padj", "qvalue"), FALSE, "FDR column")

  gene_id <- as.character(df[[gene_id_col]])
  gene_symbol_from_de <- if (!is.na(symbol_col)) clean_symbol(df[[symbol_col]]) else rep(NA_character_, nrow(df))
  gene_symbol <- gene_symbol_from_de

  identifier_source <- rep(NA_character_, nrow(df))
  direct_symbol <- !is.na(gene_symbol) & nzchar(gene_symbol)
  identifier_source[direct_symbol] <- "pig_gene_symbol_from_DE_table"

  missing_symbol <- is.na(gene_symbol) | !nzchar(gene_symbol)
  n_symbol_filled_by_mapping <- 0L

  if (!is.null(gene_id_to_symbol_map) && length(gene_id_to_symbol_map) > 0 && any(missing_symbol)) {
    missing_idx <- which(missing_symbol)
    mapped_symbol <- rep(NA_character_, length(missing_idx))
    lookup_ids <- gene_id[missing_idx]
    can_lookup <- !is.na(lookup_ids) & nzchar(lookup_ids)
    if (any(can_lookup)) {
      mapped_symbol[can_lookup] <- unname(gene_id_to_symbol_map[lookup_ids[can_lookup]])
    }
    mapped_symbol <- clean_symbol(mapped_symbol)
    valid_mapped <- !is.na(mapped_symbol) & nzchar(mapped_symbol)

    if (any(valid_mapped)) {
      fill_idx <- missing_idx[valid_mapped]
      gene_symbol[fill_idx] <- mapped_symbol[valid_mapped]
      identifier_source[fill_idx] <- "gene_id_to_symbol_mapping_from_available_annotation"
      n_symbol_filled_by_mapping <- length(fill_idx)
    }
  }

  logfc <- suppressWarnings(as.numeric(df[[logfc_col]]))
  pval_original <- suppressWarnings(as.numeric(df[[p_col]]))
  pval <- safe_p(df[[p_col]])
  fdr <- if (!is.na(fdr_col)) suppressWarnings(as.numeric(df[[fdr_col]])) else rep(NA_real_, nrow(df))

  rank_stat <- sign(logfc) * (-log10(pval))

  rank_df_all <- data.frame(
    contrast = contrast_label,
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    gene_symbol_original = gene_symbol_from_de,
    identifier_source = identifier_source,
    logFC = logfc,
    PValue_original = pval_original,
    PValue = pval,
    FDR = fdr,
    rank_stat = rank_stat,
    stringsAsFactors = FALSE
  )

  missing_or_invalid_symbol <- is.na(rank_df_all$gene_symbol) | !nzchar(rank_df_all$gene_symbol)
  invalid_numeric <- is.na(rank_df_all$logFC) | is.na(rank_df_all$PValue) |
    rank_df_all$PValue <= 0 | is.na(rank_df_all$rank_stat)

  keep <- !missing_or_invalid_symbol & !invalid_numeric
  removed_df <- rank_df_all[!keep, , drop = FALSE]
  ## Use character(nrow(.)) so zero-row audit tables do not trigger
  ## `$<-.data.frame` replacement-length errors.
  removed_df$removal_reason <- character(nrow(removed_df))
  if (nrow(removed_df) > 0) {
    removed_missing_symbol <- missing_or_invalid_symbol[!keep]
    removed_invalid_numeric <- invalid_numeric[!keep]
    removed_df$removal_reason[removed_missing_symbol] <- "no_valid_pig_gene_symbol_after_available_annotation_mapping"
    removed_df$removal_reason[removed_invalid_numeric & is.na(removed_df$removal_reason)] <- "missing_logFC_or_PValue_or_rank_stat"
  }

  rank_df <- rank_df_all[keep, , drop = FALSE]
  if (nrow(rank_df) == 0) {
    stop("No valid ranked genes remained for ", contrast_label,
         " after enforcing pig gene-symbol identifier space.")
  }

  # Keep one row per pig gene symbol using the largest absolute rank statistic.
  ord <- order(-abs(rank_df$rank_stat), rank_df$gene_symbol)
  rank_df <- rank_df[ord, , drop = FALSE]

  is_dup <- duplicated(rank_df$gene_symbol)
  duplicated_symbol_rows_removed <- sum(is_dup)
  duplicated_removed_df <- rank_df[is_dup, , drop = FALSE]
  if (nrow(duplicated_removed_df) > 0) {
    duplicated_removed_df$removal_reason <- "duplicate_pig_gene_symbol_lower_absolute_rank_stat"
  }

  rank_df <- rank_df[!is_dup, , drop = FALSE]
  rank_df <- rank_df[order(rank_df$rank_stat, decreasing = TRUE), , drop = FALSE]
  rownames(rank_df) <- NULL

  stats_vec <- rank_df$rank_stat
  names(stats_vec) <- rank_df$gene_symbol
  stats_vec <- sort(stats_vec, decreasing = TRUE)

  list(
    rank_df = rank_df,
    stats_vec = stats_vec,
    removed_invalid_rows = nrow(removed_df),
    removed_no_valid_symbol = sum(missing_or_invalid_symbol),
    removed_invalid_numeric = sum(!missing_or_invalid_symbol & invalid_numeric),
    n_symbol_filled_by_mapping = n_symbol_filled_by_mapping,
    duplicated_symbol_rows_removed = duplicated_symbol_rows_removed,
    removed_df = removed_df,
    duplicated_removed_df = duplicated_removed_df,
    n_ranked_symbols = length(stats_vec)
  )
}

clean_fgsea_out <- function(fgsea_result, contrast_label) {
  out <- as.data.frame(fgsea_result, stringsAsFactors = FALSE)
  if (nrow(out) == 0) {
    return(data.frame())
  }
  if ("leadingEdge" %in% colnames(out)) {
    out$leadingEdge_size <- vapply(out$leadingEdge, length, integer(1))
    out$leadingEdge_genes <- vapply(out$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    out$leadingEdge <- NULL
  }
  out$contrast <- contrast_label
  out <- out[, c("contrast", setdiff(colnames(out), "contrast")), drop = FALSE]
  out
}

file_audit <- function(path, label) {
  info <- file.info(path)
  data.frame(
    label = label,
    file = path,
    exists = file.exists(path),
    size_bytes = if (file.exists(path)) info$size else NA_real_,
    mtime = if (file.exists(path)) as.character(info$mtime) else NA_character_,
    stringsAsFactors = FALSE
  )
}

## =========================
## 4. Ranked gene lists
## =========================
## Identifier policy for formal Hallmark GSEA:
##   Use pig gene symbols as the ranked-vector identifiers to match msigdbr
##   Sus scrofa Hallmark gene set membership. Missing symbols are first mapped
##   from available DE-table annotation when an unambiguous gene_id -> symbol
##   mapping exists. Records still lacking a valid pig gene symbol are excluded
##   from the formal fgsea ranked vector instead of using gene_id as fallback.
gene_id_symbol_mapping <- build_gene_id_symbol_map_from_available_de_tables(
  list(t7 = de_t7, t28 = de_t28)
)
cat("\nGene ID to symbol mapping from available DE-table annotation:\n")
cat("n mapping rows:", nrow(gene_id_symbol_mapping$audit), "\n")
cat("n unambiguous gene_id symbols:", length(gene_id_symbol_mapping$map), "\n")

cat("\nBuilding ranked gene list for t7...\n")
rank_t7 <- rank_one_table(
  de_t7,
  "ACLT_untreated_t7_vs_Control",
  gene_id_to_symbol_map = gene_id_symbol_mapping$map
)
cat("t7 ranked gene list built successfully.\n")

cat("\nBuilding ranked gene list for t28...\n")
rank_t28 <- rank_one_table(
  de_t28,
  "ACLT_untreated_t28_vs_Control",
  gene_id_to_symbol_map = gene_id_symbol_mapping$map
)
cat("t28 ranked gene list built successfully.\n")

cat("\nRanked symbols:\n")
cat("t7 :", rank_t7$n_ranked_symbols, "\n")
cat("t28:", rank_t28$n_ranked_symbols, "\n")
cat("t7 symbol-filled-by-mapping:", rank_t7$n_symbol_filled_by_mapping, "\n")
cat("t28 symbol-filled-by-mapping:", rank_t28$n_symbol_filled_by_mapping, "\n")
cat("t7 removed-no-valid-symbol-or-invalid:", rank_t7$removed_invalid_rows, "\n")
cat("t28 removed-no-valid-symbol-or-invalid:", rank_t28$removed_invalid_rows, "\n")

## =========================
## 5. MSigDB Hallmark collection for Sus scrofa
## =========================
# msigdbr changed its argument names across versions. Try the current and legacy forms.
msig_df <- tryCatch(
  msigdbr::msigdbr(species = "Sus scrofa", collection = "H"),
  error = function(e) {
    msigdbr::msigdbr(species = "Sus scrofa", category = "H")
  }
)

if (nrow(msig_df) == 0) stop("msigdbr returned no Sus scrofa Hallmark gene sets.")

## Record MSigDB database version when available in the msigdbr output.
## Some msigdbr versions include a db_version column; older versions may not.
msigdb_database_version <- "not_available_in_msigdbr_output"
msigdb_database_version_source <- "db_version column absent"
if ("db_version" %in% colnames(msig_df)) {
  dbv <- unique(stats::na.omit(as.character(msig_df$db_version)))
  dbv <- dbv[nzchar(dbv)]
  if (length(dbv) > 0) {
    msigdb_database_version <- paste(dbv, collapse = "; ")
    msigdb_database_version_source <- "msig_df$db_version"
  } else {
    msigdb_database_version <- "db_version column present but empty"
    msigdb_database_version_source <- "msig_df$db_version"
  }
}
cat("MSigDB database version:", msigdb_database_version, "\n")
cat("MSigDB database version source:", msigdb_database_version_source, "\n")

pathway_col <- find_col(msig_df, c("gs_name", "gs_exact_source"), TRUE, "MSigDB pathway name column")
gene_col <- find_col(msig_df, c("gene_symbol", "db_gene_symbol"), TRUE, "MSigDB gene symbol column")

hallmark_members <- data.frame(
  pathway = as.character(msig_df[[pathway_col]]),
  gene_symbol = as.character(msig_df[[gene_col]]),
  stringsAsFactors = FALSE
)
hallmark_members <- hallmark_members[!is.na(hallmark_members$pathway) & !is.na(hallmark_members$gene_symbol), , drop = FALSE]
hallmark_members <- unique(hallmark_members)

pathways <- split(hallmark_members$gene_symbol, hallmark_members$pathway)
pathways <- lapply(pathways, unique)
pathways <- pathways[order(names(pathways))]

fgsea_min_size <- 10
fgsea_max_size <- 500
fgsea_eps <- 0
fgsea_nproc <- 1
fgsea_seed_t7 <- 1
fgsea_seed_t28 <- 1

target_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

missing_targets <- setdiff(target_pathways, names(pathways))
if (length(missing_targets) > 0) {
  stop("Target Hallmark pathways missing from Sus scrofa msigdbr result: ", paste(missing_targets, collapse = ", "))
}

cat("\nHallmark gene sets loaded:", length(pathways), "\n")
cat("fgsea parameters: minSize=", fgsea_min_size, ", maxSize=", fgsea_max_size, ", eps=", fgsea_eps, ", nproc=", fgsea_nproc, "\n", sep = "")
cat("fgsea seeds: t7=", fgsea_seed_t7, ", t28=", fgsea_seed_t28, "\n", sep = "")

## =========================
## 6. Run fgsea on full Hallmark collection
## =========================
cat("\nfgsea reproducibility settings:\n")
cat("t7 seed:", fgsea_seed_t7, "\n")
cat("t28 seed:", fgsea_seed_t28, "\n")
cat("nproc:", fgsea_nproc, "\n")

set.seed(fgsea_seed_t7)
fg_t7 <- fgsea::fgsea(
  pathways = pathways,
  stats = rank_t7$stats_vec,
  minSize = fgsea_min_size,
  maxSize = fgsea_max_size,
  eps = fgsea_eps,
  nproc = fgsea_nproc
)

set.seed(fgsea_seed_t28)
fg_t28 <- fgsea::fgsea(
  pathways = pathways,
  stats = rank_t28$stats_vec,
  minSize = fgsea_min_size,
  maxSize = fgsea_max_size,
  eps = fgsea_eps,
  nproc = fgsea_nproc
)

fg_t7_out <- clean_fgsea_out(fg_t7, "ACLT_untreated_t7_vs_Control")
fg_t28_out <- clean_fgsea_out(fg_t28, "ACLT_untreated_t28_vs_Control")
fg_all <- rbind(fg_t7_out, fg_t28_out)
fg_all <- fg_all[order(fg_all$contrast, fg_all$padj, -abs(fg_all$NES)), , drop = FALSE]

target_results <- fg_all[fg_all$pathway %in% target_pathways, , drop = FALSE]
target_results <- target_results[order(target_results$pathway, target_results$contrast), , drop = FALSE]

## =========================
## 7. Coverage and mouse-style running-curve source data for targets
## =========================
coverage_df <- data.frame(
  pathway = names(pathways),
  n_genes_in_pathway = vapply(pathways, length, integer(1)),
  n_genes_present_in_t7_rank = vapply(pathways, function(gs) sum(gs %in% names(rank_t7$stats_vec)), integer(1)),
  n_genes_present_in_t28_rank = vapply(pathways, function(gs) sum(gs %in% names(rank_t28$stats_vec)), integer(1)),
  stringsAsFactors = FALSE
)
coverage_df <- coverage_df[order(coverage_df$pathway), , drop = FALSE]

enrichment_curve_df <- function(stats, pathway_genes, pathway_name, contrast_label, nes, padj, pval) {
  stats <- sort(stats, decreasing = TRUE)
  hits <- names(stats) %in% pathway_genes
  N <- length(stats)
  Nh <- sum(hits)
  Nm <- N - Nh
  if (Nh == 0 || Nm == 0) return(data.frame())
  weights <- abs(stats)
  Phit <- ifelse(hits, weights / sum(weights[hits]), 0)
  Pmiss <- ifelse(!hits, 1 / Nm, 0)
  running <- cumsum(Phit - Pmiss)
  data.frame(
    contrast = contrast_label,
    pathway = pathway_name,
    rank = seq_along(stats),
    running_ES = running,
    is_hit = hits,
    gene_symbol = names(stats),
    NES = nes,
    padj = padj,
    pval = pval,
    stringsAsFactors = FALSE
  )
}

target_running <- list()
for (contrast_label in c("ACLT_untreated_t7_vs_Control", "ACLT_untreated_t28_vs_Control")) {
  stats_vec <- if (contrast_label == "ACLT_untreated_t7_vs_Control") rank_t7$stats_vec else rank_t28$stats_vec
  fg_this <- target_results[target_results$contrast == contrast_label, , drop = FALSE]
  for (pw in target_pathways) {
    row <- fg_this[fg_this$pathway == pw, , drop = FALSE]
    if (nrow(row) != 1) stop("Target pathway not found exactly once in target_results: ", contrast_label, " / ", pw)
    target_running[[paste(contrast_label, pw, sep = "__")]] <- enrichment_curve_df(
      stats = stats_vec,
      pathway_genes = pathways[[pw]],
      pathway_name = pw,
      contrast_label = contrast_label,
      nes = row$NES[1],
      padj = row$padj[1],
      pval = row$pval[1]
    )
  }
}
running_df <- dplyr::bind_rows(target_running)

## =========================
## 8. Mouse Figure3B-style 2 x 2 target enrichment curve plot
## =========================
pretty_pathway <- function(x) {
  x <- as.character(x)
  x[x == "HALLMARK_INFLAMMATORY_RESPONSE"] <- "Inflammatory response"
  x[x == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT / ECM remodeling"
  x
}
pretty_contrast <- function(x) {
  x <- as.character(x)
  x[x == "ACLT_untreated_t7_vs_Control"] <- "t7"
  x[x == "ACLT_untreated_t28_vs_Control"] <- "t28"
  x
}
format_fdr <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 2))
}

label_df <- target_results
label_df$pathway_pretty <- pretty_pathway(label_df$pathway)
label_df$timepoint <- pretty_contrast(label_df$contrast)
label_df$timepoint <- factor(label_df$timepoint, levels = c("t7", "t28"))
label_df$panel_label <- paste0(
  label_df$pathway_pretty,
  " — ",
  label_df$timepoint,
  "\nNES = ",
  sprintf("%.3f", label_df$NES),
  "   |   FDR = ",
  format_fdr(label_df$padj)
)

panel_order <- c(
  label_df$panel_label[label_df$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & as.character(label_df$timepoint) == "t7"],
  label_df$panel_label[label_df$pathway == "HALLMARK_INFLAMMATORY_RESPONSE" & as.character(label_df$timepoint) == "t28"],
  label_df$panel_label[label_df$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & as.character(label_df$timepoint) == "t7"],
  label_df$panel_label[label_df$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" & as.character(label_df$timepoint) == "t28"]
)

## Guard against empty labels if factor/coercion changes.
panel_order <- panel_order[nchar(panel_order) > 0]

running_df$pathway_pretty <- pretty_pathway(running_df$pathway)
running_df$timepoint <- pretty_contrast(running_df$contrast)
running_df$panel_key <- paste(running_df$pathway_pretty, running_df$timepoint, sep = "__")

label_key <- data.frame(
  pathway = as.character(label_df$pathway),
  timepoint = as.character(label_df$timepoint),
  pathway_pretty = label_df$pathway_pretty,
  panel_label = label_df$panel_label,
  stringsAsFactors = FALSE
)
label_key$panel_key <- paste(label_key$pathway_pretty, label_key$timepoint, sep = "__")

running_df <- merge(
  running_df,
  label_key[, c("panel_key", "panel_label")],
  by = "panel_key",
  all.x = TRUE
)
running_df$panel_label <- factor(running_df$panel_label, levels = panel_order)

peak_df <- aggregate(running_ES ~ panel_label, data = running_df, FUN = max)
colnames(peak_df)[2] <- "peak_ES"
rug_df <- running_df[running_df$is_hit, , drop = FALSE]
rug_df$rug_ymin <- -0.045
rug_df$rug_ymax <- 0.045

plot_files <- data.frame(plot = character(), file = character(), stringsAsFactors = FALSE)
p_combined <- ggplot2::ggplot(running_df, ggplot2::aes(x = rank, y = running_ES)) +
  ggplot2::geom_hline(
    data = peak_df,
    ggplot2::aes(yintercept = peak_ES),
    inherit.aes = FALSE,
    color = "red",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "red",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  ggplot2::geom_segment(
    data = rug_df,
    ggplot2::aes(x = rank, xend = rank, y = rug_ymin, yend = rug_ymax),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.35,
    linewidth = 0.18
  ) +
  ggplot2::geom_line(color = "#00CD00", linewidth = 0.75) +
  ggplot2::facet_wrap(~ panel_label, ncol = 2, scales = "free_y") +
  ggplot2::labs(x = NULL, y = "Running enrichment score") +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 13.5, face = "plain", lineheight = 1.15),
    panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.4),
    panel.grid.minor = ggplot2::element_line(color = "#F0F0F0", linewidth = 0.25),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
  )

combined_png <- file.path(figures_dir, "Figure4B_current78_pig_early_Hallmark_GSEA_running_curves_mouse_style.png")
combined_pdf <- file.path(figures_dir, "Figure4B_current78_pig_early_Hallmark_GSEA_running_curves_mouse_style.pdf")
ggplot2::ggsave(combined_png, p_combined, width = 11.8, height = 8.0, dpi = 320)
ggplot2::ggsave(combined_pdf, p_combined, width = 11.8, height = 8.0)
plot_files <- rbind(plot_files, data.frame(plot = "Figure4B_current78_pig_early_Hallmark_GSEA_running_curves_mouse_style", file = combined_png, stringsAsFactors = FALSE))
plot_files <- rbind(plot_files, data.frame(plot = "Figure4B_current78_pig_early_Hallmark_GSEA_running_curves_mouse_style", file = combined_pdf, stringsAsFactors = FALSE))

## Also save separate one-panel plots with the same green-curve style for flexible figure layout.
for (contrast_label in unique(running_df$contrast)) {
  for (pw in target_pathways) {
    dfp <- running_df[running_df$contrast == contrast_label & running_df$pathway == pw, , drop = FALSE]
    if (nrow(dfp) == 0) next
    peak_one <- aggregate(running_ES ~ panel_label, data = dfp, FUN = max)
    colnames(peak_one)[2] <- "peak_ES"
    rug_one <- dfp[dfp$is_hit, , drop = FALSE]
    rug_one$rug_ymin <- -0.045
    rug_one$rug_ymax <- 0.045
    p_one <- ggplot2::ggplot(dfp, ggplot2::aes(x = rank, y = running_ES)) +
      ggplot2::geom_hline(data = peak_one, ggplot2::aes(yintercept = peak_ES), inherit.aes = FALSE, color = "red", linetype = "dashed", linewidth = 0.45) +
      ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.35) +
      ggplot2::geom_segment(data = rug_one, ggplot2::aes(x = rank, xend = rank, y = rug_ymin, yend = rug_ymax), inherit.aes = FALSE, color = "black", alpha = 0.35, linewidth = 0.18) +
      ggplot2::geom_line(color = "#00CD00", linewidth = 0.75) +
      ggplot2::facet_wrap(~ panel_label, ncol = 1, scales = "free_y") +
      ggplot2::labs(x = NULL, y = "Running enrichment score") +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 13.5, face = "plain", lineheight = 1.15),
        panel.grid.major = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.4),
        panel.grid.minor = ggplot2::element_line(color = "#F0F0F0", linewidth = 0.25),
        axis.title.y = ggplot2::element_text(size = 14),
        axis.text = ggplot2::element_text(size = 10),
        axis.title.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
      )
    safe_name <- gsub("[^A-Za-z0-9]+", "_", paste0(contrast_label, "_", pw, "_mouse_style"))
    png_file <- file.path(figures_dir, paste0("Figure4B_current78_", safe_name, ".png"))
    pdf_file <- file.path(figures_dir, paste0("Figure4B_current78_", safe_name, ".pdf"))
    ggplot2::ggsave(png_file, p_one, width = 5.9, height = 4.0, dpi = 320)
    ggplot2::ggsave(pdf_file, p_one, width = 5.9, height = 4.0)
    plot_files <- rbind(plot_files, data.frame(plot = safe_name, file = png_file, stringsAsFactors = FALSE))
    plot_files <- rbind(plot_files, data.frame(plot = safe_name, file = pdf_file, stringsAsFactors = FALSE))
  }
}

## =========================
## 9. Save outputs
## =========================
write.csv(rank_t7$rank_df, file.path(out_dir, "step20_current78_ranked_genes_t7_vs_Control.csv"), row.names = FALSE)
write.csv(rank_t28$rank_df, file.path(out_dir, "step20_current78_ranked_genes_t28_vs_Control.csv"), row.names = FALSE)
write.csv(rank_t7$removed_df, file.path(out_dir, "step20_current78_ranked_genes_t7_removed_no_valid_symbol_or_invalid_values.csv"), row.names = FALSE)
write.csv(rank_t28$removed_df, file.path(out_dir, "step20_current78_ranked_genes_t28_removed_no_valid_symbol_or_invalid_values.csv"), row.names = FALSE)
write.csv(rank_t7$duplicated_removed_df, file.path(out_dir, "step20_current78_ranked_genes_t7_removed_duplicate_symbols.csv"), row.names = FALSE)
write.csv(rank_t28$duplicated_removed_df, file.path(out_dir, "step20_current78_ranked_genes_t28_removed_duplicate_symbols.csv"), row.names = FALSE)
write.csv(gene_id_symbol_mapping$audit, file.path(out_dir, "step20_current78_available_gene_id_to_symbol_mapping_audit.csv"), row.names = FALSE)
write.csv(hallmark_members, file.path(out_dir, "step20_current78_sus_scrofa_hallmark_membership.csv"), row.names = FALSE)
write.csv(coverage_df, file.path(out_dir, "step20_current78_hallmark_gene_set_coverage.csv"), row.names = FALSE)
write.csv(fg_t7_out, file.path(out_dir, "step20_current78_fgsea_full_hallmark_t7_vs_Control.csv"), row.names = FALSE)
write.csv(fg_t28_out, file.path(out_dir, "step20_current78_fgsea_full_hallmark_t28_vs_Control.csv"), row.names = FALSE)
write.csv(fg_all, file.path(out_dir, "step20_current78_fgsea_full_hallmark_combined.csv"), row.names = FALSE)
write.csv(target_results, file.path(out_dir, "step20_current78_target_pathway_results_from_full_hallmark.csv"), row.names = FALSE)
write.csv(running_df, file.path(out_dir, "step20_current78_target_pathway_running_curve_source_data.csv"), row.names = FALSE)
write.csv(plot_files, file.path(out_dir, "step20_current78_target_pathway_plot_files.csv"), row.names = FALSE)

version_df <- data.frame(
  item = c(
    "R_version",
    "fgsea_package_version",
    "msigdbr_package_version",
    "ggplot2_package_version",
    "dplyr_package_version",
    "msigdbr_species",
    "msigdbr_collection",
    "MSigDB_database_version",
    "MSigDB_database_version_source",
    "msigdbr_output_columns",
    "fgsea_minSize",
    "fgsea_maxSize",
    "fgsea_eps",
    "fgsea_nproc",
    "fgsea_seed_t7",
    "fgsea_seed_t28",
    "rank_statistic",
    "ranked_vector_identifier_policy",
    "missing_symbol_handling",
    "target_pathways_note"
  ),
  value = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("fgsea")),
    as.character(utils::packageVersion("msigdbr")),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("dplyr")),
    "Sus scrofa",
    "MSigDB Hallmark collection H",
    msigdb_database_version,
    msigdb_database_version_source,
    paste(colnames(msig_df), collapse = "; "),
    fgsea_min_size,
    fgsea_max_size,
    fgsea_eps,
    fgsea_nproc,
    fgsea_seed_t7,
    fgsea_seed_t28,
    "sign(logFC) * -log10(PValue)",
    "Pig gene symbols are used as formal fgsea ranked-vector identifiers to match msigdbr Sus scrofa Hallmark membership.",
    "Records lacking pig gene symbols are mapped from available DE-table gene_id-to-symbol annotation when unambiguous; unresolved records are excluded and audited, not replaced by gene_id fallback.",
    "fgsea was run on the full Hallmark collection; target pathways were extracted afterward."
  ),
  stringsAsFactors = FALSE
)
write.csv(version_df, file.path(out_dir, "step20_current78_software_versions_and_parameters.csv"), row.names = FALSE)

fgsea_reproducibility_df <- data.frame(
  contrast = c("ACLT_untreated_t7_vs_Control", "ACLT_untreated_t28_vs_Control"),
  seed = c(fgsea_seed_t7, fgsea_seed_t28),
  nproc = c(fgsea_nproc, fgsea_nproc),
  minSize = c(fgsea_min_size, fgsea_min_size),
  maxSize = c(fgsea_max_size, fgsea_max_size),
  eps = c(fgsea_eps, fgsea_eps),
  fgsea_implementation = "fgsea default multilevel implementation; no fixed nperm specified",
  stringsAsFactors = FALSE
)
write.csv(fgsea_reproducibility_df, file.path(out_dir, "step20_current78_fgsea_reproducibility_settings.csv"), row.names = FALSE)

input_audit <- rbind(
  file_audit(t7_file, "Step17 t7 vs Control QLF DE table"),
  file_audit(t28_file, "Step17 t28 vs Control QLF DE table")
)
write.csv(input_audit, file.path(out_dir, "step20_current78_input_file_audit.csv"), row.names = FALSE)

summary_df <- data.frame(
  metric = c(
    "t7_DE_input_file",
    "t28_DE_input_file",
    "n_ranked_symbols_t7",
    "n_ranked_symbols_t28",
    "n_symbol_filled_by_mapping_t7",
    "n_symbol_filled_by_mapping_t28",
    "n_removed_no_valid_symbol_or_invalid_t7",
    "n_removed_no_valid_symbol_or_invalid_t28",
    "n_duplicate_symbols_removed_t7",
    "n_duplicate_symbols_removed_t28",
    "identifier_policy",
    "n_full_hallmark_gene_sets_loaded",
    "MSigDB_database_version",
    "MSigDB_database_version_source",
    "n_full_hallmark_gene_sets_tested_t7",
    "n_full_hallmark_gene_sets_tested_t28",
    "fgsea_minSize",
    "fgsea_maxSize",
    "fgsea_eps",
    "fgsea_nproc",
    "fgsea_seed_t7",
    "fgsea_seed_t28",
    "rank_statistic",
    "target_pathways",
    "full_collection_FDR_used",
    "output_dir"
  ),
  value = c(
    t7_file,
    t28_file,
    rank_t7$n_ranked_symbols,
    rank_t28$n_ranked_symbols,
    rank_t7$n_symbol_filled_by_mapping,
    rank_t28$n_symbol_filled_by_mapping,
    rank_t7$removed_invalid_rows,
    rank_t28$removed_invalid_rows,
    rank_t7$duplicated_symbol_rows_removed,
    rank_t28$duplicated_symbol_rows_removed,
    "pig gene symbol only; missing symbols mapped from available annotation when unambiguous; unresolved records excluded, no gene_id fallback in formal ranked vector",
    length(pathways),
    msigdb_database_version,
    msigdb_database_version_source,
    nrow(fg_t7_out),
    nrow(fg_t28_out),
    fgsea_min_size,
    fgsea_max_size,
    fgsea_eps,
    fgsea_nproc,
    fgsea_seed_t7,
    fgsea_seed_t28,
    "sign(logFC) * -log10(PValue)",
    paste(target_pathways, collapse = "; "),
    TRUE,
    out_dir
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "step20_current78_pig_early_hallmark_gsea_summary.csv"), row.names = FALSE)

save(
  de_t7, de_t28,
  rank_t7, rank_t28,
  hallmark_members, pathways,
  fg_t7, fg_t28, fg_t7_out, fg_t28_out, fg_all,
  target_results, coverage_df, running_df,
  version_df, fgsea_reproducibility_df, summary_df,
  file = file.path(objects_dir, "step20_current78_pig_early_hallmark_gsea_workspace.RData")
)

## =========================
## 10. Summary-to-send log
## =========================
summary_lines <- c(
  "===== STEP20_CURRENT78 SUMMARY TO SEND =====",
  paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Main run summary:",
  paste(capture.output(print(summary_df)), collapse = "\n"),
  "",
  "Target pathway results from full Hallmark fgsea:",
  paste(capture.output(print(target_results)), collapse = "\n"),
  "",
  "Hallmark coverage for target pathways:",
  paste(capture.output(print(coverage_df[coverage_df$pathway %in% target_pathways, , drop = FALSE])), collapse = "\n"),
  "",
  "Version and parameter records:",
  paste(capture.output(print(version_df)), collapse = "\n"),
  "",
  "fgsea reproducibility settings:",
  paste(capture.output(print(fgsea_reproducibility_df)), collapse = "\n"),
  "",
  "Primary output files:",
  paste("full_combined_fgsea:", file.path(out_dir, "step20_current78_fgsea_full_hallmark_combined.csv")),
  paste("target_pathway_results:", file.path(out_dir, "step20_current78_target_pathway_results_from_full_hallmark.csv")),
  paste("ranked_t7:", file.path(out_dir, "step20_current78_ranked_genes_t7_vs_Control.csv")),
  paste("ranked_t28:", file.path(out_dir, "step20_current78_ranked_genes_t28_vs_Control.csv")),
  paste("removed_t7:", file.path(out_dir, "step20_current78_ranked_genes_t7_removed_no_valid_symbol_or_invalid_values.csv")),
  paste("removed_t28:", file.path(out_dir, "step20_current78_ranked_genes_t28_removed_no_valid_symbol_or_invalid_values.csv")),
  paste("running_curve_source_data:", file.path(out_dir, "step20_current78_target_pathway_running_curve_source_data.csv")),
  paste("summary_csv:", file.path(out_dir, "step20_current78_pig_early_hallmark_gsea_summary.csv")),
  paste("fgsea_reproducibility_settings:", file.path(out_dir, "step20_current78_fgsea_reproducibility_settings.csv")),
  "",
  "Plots:",
  paste(capture.output(print(plot_files)), collapse = "\n"),
  "",
  "Interpretation checkpoint:",
  "Please verify that fgsea was run on the full Hallmark collection and that target pathway FDR values are from the full collection-level analysis.",
  "Current parameters should be minSize = 10, maxSize = 500, eps = 0, nproc = 1, seed_t7 = 1, seed_t28 = 1."
)
writeLines(summary_lines, summary_log_file, useBytes = TRUE)

cat("\n===== STEP20_CURRENT78 SUMMARY =====\n")
print(summary_df)
cat("\nTarget pathway results from full Hallmark fgsea:\n")
print(target_results)
cat("\nStep20 current78 pig early Hallmark GSEA completed successfully.\n")
cat("Summary-to-send log: ", summary_log_file, "\n", sep = "")

sink()

cat("\nStep20 current78 pig early Hallmark GSEA completed.\n")
cat("Summary-to-send log: ", summary_log_file, "\n", sep = "")
