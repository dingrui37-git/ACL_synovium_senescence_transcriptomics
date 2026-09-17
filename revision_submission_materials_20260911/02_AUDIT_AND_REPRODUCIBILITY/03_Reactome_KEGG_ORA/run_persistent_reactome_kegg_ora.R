options(stringsAsFactors = FALSE, timeout = 180)

## Reviewer-requested complementary pathway analysis.
## Formal input: mouse genes that are strict DEGs at both 1W and 4W and retain
## the same direction. The CellAge-filtered 78-gene subset is deliberately not
## used for ORA because that would introduce circular pathway interpretation.

out_dir <- "D:/workspace/ACLsenescence2_reviewer3_persistent_ORA"
tmp_dir <- file.path(out_dir, "tmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = tmp_dir, TEMP = tmp_dir, TMP = tmp_dir)

input_dir <- "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step07_strict_DEG_upset_persistent"
persistent_file <- file.path(input_dir, "step07_persistent_direction_consistent_genes.csv")
de_1w_file <- file.path(input_dir, "step07_DE_1W_standardized.csv")
de_4w_file <- file.path(input_dir, "step07_DE_4W_standardized.csv")

required_packages <- c("msigdbr", "KEGGREST")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

min_gene_set_size <- 10L
max_gene_set_size <- 500L
strict_fdr <- 0.05
strict_abs_logfc <- 1
ora_alternative <- "over-representation; one-sided hypergeometric upper tail"
primary_multiplicity <- paste(
  "BH within direction across all tested mouse-native Reactome and KEGG pathways pooled;",
  "up- and down-regulated persistent lists are separate prespecified families"
)

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  message("Saved: ", path)
}

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "none", "None")] <- NA_character_
  sub("\\.0$", "", x)
}

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL", "none", "None")] <- NA_character_
  x
}

file_audit <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
  info <- file.info(path)
  data.frame(
    label = label,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    md5 = unname(tools::md5sum(path)),
    stringsAsFactors = FALSE
  )
}

required_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(!is.finite(x))) return(NA_real_)
  stats::median(x[is.finite(x)])
}

## 1. Freeze and validate the analysis input.
input_audit <- do.call(rbind, list(
  file_audit(persistent_file, "Step07 direction-consistent persistent genes"),
  file_audit(de_1w_file, "Step07 standardized 1W DE universe"),
  file_audit(de_4w_file, "Step07 standardized 4W DE universe")
))
write_csv(input_audit, file.path(out_dir, "persistent_ORA_input_file_audit.csv"))

persistent <- utils::read.csv(persistent_file, stringsAsFactors = FALSE, check.names = FALSE)
de_1w <- utils::read.csv(de_1w_file, stringsAsFactors = FALSE, check.names = FALSE)
de_4w <- utils::read.csv(de_4w_file, stringsAsFactors = FALSE, check.names = FALSE)

required_columns(
  persistent,
  c("ENTREZID", "SYMBOL", "logFC_1W", "FDR_1W", "direction_1W",
    "logFC_4W", "FDR_4W", "direction_4W", "direction_persistent", "category"),
  "Persistent table"
)
required_columns(de_1w, c("ENTREZID", "SYMBOL", "logFC", "FDR"), "1W DE table")
required_columns(de_4w, c("ENTREZID", "SYMBOL", "logFC", "FDR"), "4W DE table")

for (nm in c("logFC_1W", "FDR_1W", "logFC_4W", "FDR_4W")) {
  persistent[[nm]] <- suppressWarnings(as.numeric(persistent[[nm]]))
}
persistent$ENTREZID <- clean_id(persistent$ENTREZID)
persistent$SYMBOL <- clean_symbol(persistent$SYMBOL)
persistent$direction_persistent <- trimws(as.character(persistent$direction_persistent))

de_1w$ENTREZID <- clean_id(de_1w$ENTREZID)
de_4w$ENTREZID <- clean_id(de_4w$ENTREZID)
de_1w$logFC <- suppressWarnings(as.numeric(de_1w$logFC))
de_4w$logFC <- suppressWarnings(as.numeric(de_4w$logFC))
de_1w$FDR <- suppressWarnings(as.numeric(de_1w$FDR))
de_4w$FDR <- suppressWarnings(as.numeric(de_4w$FDR))

valid_de_1w <- !is.na(de_1w$ENTREZID) & is.finite(de_1w$logFC) & is.finite(de_1w$FDR)
valid_de_4w <- !is.na(de_4w$ENTREZID) & is.finite(de_4w$logFC) & is.finite(de_4w$FDR)
background_entrez <- sort(intersect(unique(de_1w$ENTREZID[valid_de_1w]),
                                    unique(de_4w$ENTREZID[valid_de_4w])))
background_symbol_rows <- rbind(
  data.frame(gene_id = de_1w$ENTREZID, gene_symbol = clean_symbol(de_1w$SYMBOL)),
  data.frame(gene_id = de_4w$ENTREZID, gene_symbol = clean_symbol(de_4w$SYMBOL))
)
background_symbol_rows <- unique(background_symbol_rows[
  !is.na(background_symbol_rows$gene_id) & !is.na(background_symbol_rows$gene_symbol), , drop = FALSE
])
background_symbol_map <- vapply(
  split(background_symbol_rows$gene_symbol, background_symbol_rows$gene_id),
  function(x) paste(sort(unique(x)), collapse = ";"),
  character(1)
)

if (anyDuplicated(persistent$ENTREZID)) {
  stop("Persistent input contains duplicated Entrez IDs.", call. = FALSE)
}
if (nrow(persistent) != 1416L) {
  stop("Expected 1,416 persistent genes, found ", nrow(persistent), ".", call. = FALSE)
}
if (!all(persistent$category == "persistent_direction_consistent")) {
  stop("Persistent category audit failed.", call. = FALSE)
}
if (!all(persistent$FDR_1W < strict_fdr & persistent$FDR_4W < strict_fdr &
         abs(persistent$logFC_1W) > strict_abs_logfc &
         abs(persistent$logFC_4W) > strict_abs_logfc)) {
  stop("Persistent strict-DE threshold audit failed.", call. = FALSE)
}
direction_expected <- ifelse(persistent$logFC_1W > 0 & persistent$logFC_4W > 0, "Up",
                             ifelse(persistent$logFC_1W < 0 & persistent$logFC_4W < 0, "Down", NA_character_))
if (any(is.na(direction_expected)) || !all(direction_expected == persistent$direction_persistent)) {
  stop("Persistent direction-consistency audit failed.", call. = FALSE)
}
if (!all(persistent$ENTREZID %in% background_entrez)) {
  stop("One or more persistent genes are absent from the two-timepoint DE universe.", call. = FALSE)
}

persistent$mean_abs_logFC <- rowMeans(abs(cbind(persistent$logFC_1W, persistent$logFC_4W)))
persistent$delta_logFC_4W_minus_1W <- persistent$logFC_4W - persistent$logFC_1W
persistent$abs_effect_stronger_at_4W <- abs(persistent$logFC_4W) > abs(persistent$logFC_1W)

gene_lists <- list(
  persistent_up = persistent$ENTREZID[persistent$direction_persistent == "Up"],
  persistent_down = persistent$ENTREZID[persistent$direction_persistent == "Down"],
  persistent_all = persistent$ENTREZID
)
gene_list_roles <- c(
  persistent_up = "primary direction-stratified ORA",
  persistent_down = "primary direction-stratified ORA",
  persistent_all = "secondary combined-list context"
)

write_csv(persistent, file.path(out_dir, "persistent_ORA_input_genes_1416.csv"))
write_csv(persistent[persistent$direction_persistent == "Up", , drop = FALSE],
          file.path(out_dir, "persistent_ORA_input_genes_up_970.csv"))
write_csv(persistent[persistent$direction_persistent == "Down", , drop = FALSE],
          file.path(out_dir, "persistent_ORA_input_genes_down_446.csv"))

## 2. Load species-native Reactome and KEGG pathway definitions.
## The common eligibility universe is defined before database annotation is
## applied; each database then uses its annotated subset as the same testing
## universe for persistent-up and persistent-down lists.
reactome_raw <- msigdbr::msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M2",
  subcollection = "CP:REACTOME"
)
if (nrow(reactome_raw) == 0L) stop("No mouse-native Reactome gene sets returned.", call. = FALSE)

reactome_membership <- unique(data.frame(
  database = "Reactome",
  database_release = as.character(reactome_raw$db_version),
  pathway_id = as.character(reactome_raw$gs_id),
  pathway = as.character(reactome_raw$gs_name),
  pathway_name = as.character(reactome_raw$gs_description),
  gene_id = clean_id(reactome_raw$ncbi_gene),
  gene_symbol = clean_symbol(reactome_raw$gene_symbol),
  stringsAsFactors = FALSE
))
reactome_membership <- reactome_membership[!is.na(reactome_membership$gene_id) &
                                             !is.na(reactome_membership$pathway), , drop = FALSE]

kegg_release_text <- as.character(KEGGREST::keggInfo("kegg"))
kegg_release_line <- strsplit(kegg_release_text, "\n", fixed = TRUE)[[1]]
kegg_release_line <- kegg_release_line[grepl("^[[:space:]]*pathway", kegg_release_line)]
kegg_release_date <- if (length(kegg_release_line) > 0L) {
  sub(".*([0-9]{4}/[0-9]{2}/[0-9]{2}).*", "\\1", trimws(kegg_release_line[1]))
} else {
  "release_not_reported"
}
kegg_release <- gsub("/", "-", kegg_release_date, fixed = TRUE)
kegg_pathway_names <- KEGGREST::keggList("pathway", "mmu")
kegg_links <- KEGGREST::keggLink("pathway", "mmu")

kegg_membership <- data.frame(
  database = "KEGG",
  database_release = kegg_release,
  pathway_id = sub("^path:", "", unname(kegg_links)),
  pathway = sub("^path:", "", unname(kegg_links)),
  gene_id = sub("^mmu:", "", names(kegg_links)),
  stringsAsFactors = FALSE
)
kegg_name_map <- data.frame(
  pathway_id = names(kegg_pathway_names),
  pathway_name = sub(" - Mus musculus \\(house mouse\\)$", "", unname(kegg_pathway_names)),
  stringsAsFactors = FALSE
)
kegg_membership <- merge(kegg_membership, kegg_name_map, by = "pathway_id", all.x = TRUE, sort = FALSE)
kegg_membership$gene_symbol <- unname(background_symbol_map[kegg_membership$gene_id])
kegg_membership <- unique(kegg_membership[, c(
  "database", "database_release", "pathway_id", "pathway", "pathway_name", "gene_id", "gene_symbol"
)])

membership_all <- rbind(reactome_membership, kegg_membership)
membership_all <- membership_all[order(membership_all$database, membership_all$pathway, membership_all$gene_id), ]
write_csv(membership_all, file.path(out_dir, "persistent_ORA_pathway_membership_snapshot.csv"))
writeLines(kegg_release_text, file.path(out_dir, "persistent_ORA_KEGG_release_info.txt"))

## 3. One-sided hypergeometric ORA.
run_database_ora <- function(database_name, membership) {
  database_genes <- unique(membership$gene_id)
  analysis_universe <- sort(intersect(background_entrez, database_genes))
  pathways_all <- split(membership$gene_id, membership$pathway)
  pathways_all <- lapply(pathways_all, function(x) sort(unique(intersect(x, analysis_universe))))
  pathway_sizes <- lengths(pathways_all)
  tested_names <- names(pathways_all)[pathway_sizes >= min_gene_set_size &
                                       pathway_sizes <= max_gene_set_size]
  pathways <- pathways_all[tested_names]

  metadata <- unique(membership[, c("pathway", "pathway_id", "pathway_name", "database_release")])
  metadata <- metadata[match(tested_names, metadata$pathway), , drop = FALSE]

  universe_table <- data.frame(
    database = database_name,
    gene_id = analysis_universe,
    gene_symbol = unname(background_symbol_map[analysis_universe]),
    persistent_gene = analysis_universe %in% persistent$ENTREZID,
    persistent_direction = persistent$direction_persistent[match(analysis_universe, persistent$ENTREZID)],
    stringsAsFactors = FALSE
  )
  write_csv(universe_table, file.path(out_dir, paste0("persistent_ORA_universe_", database_name, ".csv")))

  output <- list()
  output_i <- 0L
  for (list_name in names(gene_lists)) {
    query <- sort(unique(intersect(gene_lists[[list_name]], analysis_universe)))
    N <- length(analysis_universe)
    n <- length(query)

    rows <- lapply(seq_along(pathways), function(i) {
      members <- pathways[[i]]
      overlap <- intersect(query, members)
      K <- length(members)
      k <- length(overlap)
      p <- stats::phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
      expected <- n * K / N
      enrichment_ratio <- if (expected > 0) k / expected else NA_real_
      odds_ratio <- ((k + 0.5) * (N - K - n + k + 0.5)) /
        ((n - k + 0.5) * (K - k + 0.5))

      pg <- persistent[match(overlap, persistent$ENTREZID), , drop = FALSE]
      pg <- pg[order(-pg$mean_abs_logFC, pg$SYMBOL), , drop = FALSE]
      data.frame(
        database = database_name,
        database_release = metadata$database_release[i],
        gene_list = list_name,
        analysis_role = unname(gene_list_roles[list_name]),
        pathway_id = metadata$pathway_id[i],
        pathway = metadata$pathway[i],
        pathway_name = metadata$pathway_name[i],
        universe_size = N,
        query_size_in_universe = n,
        pathway_size_in_universe = K,
        overlap_size = k,
        expected_overlap = expected,
        enrichment_ratio = enrichment_ratio,
        odds_ratio_haldane = odds_ratio,
        pval = p,
        median_logFC_1W_overlap = safe_median(pg$logFC_1W),
        median_logFC_4W_overlap = safe_median(pg$logFC_4W),
        median_abs_logFC_1W_overlap = safe_median(abs(pg$logFC_1W)),
        median_abs_logFC_4W_overlap = safe_median(abs(pg$logFC_4W)),
        median_delta_logFC_4W_minus_1W = safe_median(pg$delta_logFC_4W_minus_1W),
        proportion_overlap_stronger_abs_effect_at_4W = if (nrow(pg) > 0L) mean(pg$abs_effect_stronger_at_4W) else NA_real_,
        overlap_gene_ids = paste(pg$ENTREZID, collapse = ";"),
        overlap_gene_symbols = paste(pg$SYMBOL, collapse = ";"),
        stringsAsFactors = FALSE
      )
    })
    result <- do.call(rbind, rows)
    result$FDR_within_database <- stats::p.adjust(result$pval, method = "BH")
    output_i <- output_i + 1L
    output[[output_i]] <- result
  }

  result_all <- do.call(rbind, output)
  audit <- data.frame(
    database = database_name,
    database_release = paste(unique(membership$database_release), collapse = "; "),
    background_genes_tested_at_both_timepoints = length(background_entrez),
    database_annotated_background_genes = length(analysis_universe),
    gene_sets_loaded = length(pathways_all),
    gene_sets_tested_size_10_to_500 = length(pathways),
    persistent_all_in_database_universe = length(intersect(gene_lists$persistent_all, analysis_universe)),
    persistent_up_in_database_universe = length(intersect(gene_lists$persistent_up, analysis_universe)),
    persistent_down_in_database_universe = length(intersect(gene_lists$persistent_down, analysis_universe)),
    min_gene_set_size = min_gene_set_size,
    max_gene_set_size = max_gene_set_size,
    stringsAsFactors = FALSE
  )
  list(result = result_all, audit = audit)
}

reactome_run <- run_database_ora("Reactome", reactome_membership)
kegg_run <- run_database_ora("KEGG", kegg_membership)

ora <- rbind(reactome_run$result, kegg_run$result)
run_audit <- rbind(reactome_run$audit, kegg_run$audit)

## Primary family: pool all tested Reactome and KEGG pathways within each
## prespecified direction. The combined 1,416-gene list remains secondary.
ora$FDR_joint_Reactome_KEGG <- NA_real_
for (list_name in unique(ora$gene_list)) {
  idx <- which(ora$gene_list == list_name)
  ora$FDR_joint_Reactome_KEGG[idx] <- stats::p.adjust(ora$pval[idx], method = "BH")
}
ora$primary_test <- ora$gene_list %in% c("persistent_up", "persistent_down")
ora$primary_FDR <- ifelse(ora$primary_test, ora$FDR_joint_Reactome_KEGG, NA_real_)
ora$primary_FDR_status <- ifelse(
  !ora$primary_test, "secondary_combined_list",
  ifelse(ora$primary_FDR < 0.05, "FDR_supported", ifelse(ora$pval < 0.05, "nominal_only", "not_significant"))
)

ora <- ora[order(match(ora$gene_list, c("persistent_up", "persistent_down", "persistent_all")),
                 ora$FDR_joint_Reactome_KEGG, ora$pval, ora$database, ora$pathway), ]
rownames(ora) <- NULL

significant_primary <- ora[ora$primary_test & !is.na(ora$primary_FDR) & ora$primary_FDR < 0.05, , drop = FALSE]
significant_within_database <- ora[!is.na(ora$FDR_within_database) & ora$FDR_within_database < 0.05, , drop = FALSE]

top_primary <- do.call(rbind, lapply(split(significant_primary, list(significant_primary$gene_list,
                                                                     significant_primary$database), drop = TRUE),
                                     function(x) head(x[order(x$primary_FDR, -x$enrichment_ratio, x$pathway), ], 20L)))
if (is.null(top_primary)) top_primary <- significant_primary[0, , drop = FALSE]

summary_by_family <- do.call(rbind, lapply(split(ora, list(ora$gene_list, ora$database), drop = TRUE), function(x) {
  data.frame(
    gene_list = x$gene_list[1],
    analysis_role = x$analysis_role[1],
    database = x$database[1],
    n_tested_pathways = nrow(x),
    n_nominal_p_lt_0.05 = sum(x$pval < 0.05, na.rm = TRUE),
    n_within_database_FDR_lt_0.05 = sum(x$FDR_within_database < 0.05, na.rm = TRUE),
    n_joint_Reactome_KEGG_FDR_lt_0.05 = sum(x$FDR_joint_Reactome_KEGG < 0.05, na.rm = TRUE),
    min_joint_FDR = min(x$FDR_joint_Reactome_KEGG, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
summary_by_family <- summary_by_family[order(match(summary_by_family$gene_list,
                                                   c("persistent_up", "persistent_down", "persistent_all")),
                                             summary_by_family$database), ]

input_summary <- data.frame(
  metric = c(
    "persistent_genes_total", "persistent_genes_up", "persistent_genes_down",
    "unique_persistent_entrez", "unique_persistent_symbols",
    "background_entrez_tested_at_both_1W_and_4W",
    "persistent_up_stronger_abs_effect_at_4W", "persistent_down_stronger_abs_effect_at_4W",
    "persistent_up_median_abs_logFC_1W", "persistent_up_median_abs_logFC_4W",
    "persistent_down_median_abs_logFC_1W", "persistent_down_median_abs_logFC_4W",
    "persistent_definition", "primary_ORA_gene_lists", "secondary_ORA_gene_list",
    "primary_FDR_scope", "within_database_FDR_role", "ORA_test",
    "gene_set_size_filter", "Reactome_source", "KEGG_source",
    "CellAge_78_used_as_ORA_input"
  ),
  value = c(
    nrow(persistent), sum(persistent$direction_persistent == "Up"),
    sum(persistent$direction_persistent == "Down"),
    length(unique(persistent$ENTREZID)), length(unique(stats::na.omit(persistent$SYMBOL))),
    length(background_entrez),
    paste0(sum(persistent$direction_persistent == "Up" & persistent$abs_effect_stronger_at_4W),
           "/", sum(persistent$direction_persistent == "Up")),
    paste0(sum(persistent$direction_persistent == "Down" & persistent$abs_effect_stronger_at_4W),
           "/", sum(persistent$direction_persistent == "Down")),
    safe_median(abs(persistent$logFC_1W[persistent$direction_persistent == "Up"])),
    safe_median(abs(persistent$logFC_4W[persistent$direction_persistent == "Up"])),
    safe_median(abs(persistent$logFC_1W[persistent$direction_persistent == "Down"])),
    safe_median(abs(persistent$logFC_4W[persistent$direction_persistent == "Down"])),
    "FDR < 0.05 and |logFC| > 1 at both 1W and 4W, with the same logFC direction",
    "persistent_up; persistent_down", "persistent_all",
    primary_multiplicity,
    "secondary multiplicity-scope result",
    ora_alternative,
    paste0(min_gene_set_size, "-", max_gene_set_size,
           " genes in the database-annotated subset of the common eligibility universe"),
    paste0("MSigDB mouse-native M2:CP:REACTOME ", paste(unique(reactome_membership$database_release), collapse = ";")),
    paste0("live KEGG REST mouse pathway snapshot; ", kegg_release),
    "No; excluded to avoid circular enrichment after CellAge-based selection"
  ),
  stringsAsFactors = FALSE
)

write_csv(ora, file.path(out_dir, "persistent_ORA_all_Reactome_KEGG_results.csv"))
write_csv(significant_primary, file.path(out_dir, "persistent_ORA_primary_FDR_significant_results.csv"))
write_csv(significant_within_database, file.path(out_dir, "persistent_ORA_within_database_FDR_significant_results.csv"))
write_csv(top_primary, file.path(out_dir, "persistent_ORA_top20_primary_results_per_direction_database.csv"))
write_csv(summary_by_family, file.path(out_dir, "persistent_ORA_result_counts_summary.csv"))
write_csv(run_audit, file.path(out_dir, "persistent_ORA_database_and_universe_audit.csv"))
write_csv(input_summary, file.path(out_dir, "persistent_ORA_analysis_specification.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "persistent_ORA_sessionInfo.txt"))

cat("\n===== Persistent Reactome/KEGG ORA completed =====\n")
cat("Persistent genes:", nrow(persistent), "(up", sum(persistent$direction_persistent == "Up"),
    "; down", sum(persistent$direction_persistent == "Down"), ")\n")
print(run_audit, row.names = FALSE)
cat("\nResult counts:\n")
print(summary_by_family, row.names = FALSE)
cat("\nPrimary joint-FDR-significant pathways:", nrow(significant_primary), "\n")
