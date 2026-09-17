# Build a supplementary table contrasting persistent, 1-week-only and 4-week-only
# mouse strict-DE programs with representative genes and Reactome/KEGG ORA.
#
# Locked mouse DEG definition (inherited from Step07): FDR < 0.05 and |logFC| > 1.
# Persistent = strict DE at both time points with the same logFC direction.
# 1W-only / 4W-only = strict DE at only the indicated time point.
# ORA is direction-stratified and uses the common two-timepoint tested universe,
# with database-specific annotation subsets and 10-500 gene pathway size limits.
# Reactome and KEGG pathways are pooled for BH correction within each of the six
# temporal-class-by-direction lists (1,287 tested pathways per list).

options(stringsAsFactors = FALSE, scipen = 999)

out_dir <- "D:/workspace/ACLsenescence2_temporal_program_summary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

structure_file <- paste0(
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/",
  "step07_strict_DEG_upset_persistent/step07_strict_DEG_structure_gene_table.csv"
)
de1_file <- paste0(
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/",
  "step07_strict_DEG_upset_persistent/step07_DE_1W_standardized.csv"
)
de4_file <- paste0(
  "E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/",
  "step07_strict_DEG_upset_persistent/step07_DE_4W_standardized.csv"
)
ora_dir <- "D:/workspace/ACLsenescence2_reviewer3_persistent_ORA"
membership_file <- file.path(ora_dir, "persistent_ORA_pathway_membership_snapshot.csv")
universe_R_file <- file.path(ora_dir, "persistent_ORA_universe_Reactome.csv")
universe_K_file <- file.path(ora_dir, "persistent_ORA_universe_KEGG.csv")

required <- c(structure_file, de1_file, de4_file, membership_file,
              universe_R_file, universe_K_file)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = " | "))

read_utf8_csv <- function(path) {
  # The archived pathway snapshot contains a few non-UTF8 description bytes;
  # base R's native reader preserves the complete CSV, whereas forcing UTF-8
  # can truncate the file at the first invalid byte.
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA"))
}

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

structure <- read_utf8_csv(structure_file)
de1 <- read_utf8_csv(de1_file)
de4 <- read_utf8_csv(de4_file)
membership <- read_utf8_csv(membership_file)
universe_R <- read_utf8_csv(universe_R_file)
universe_K <- read_utf8_csv(universe_K_file)

for (nm in c("ENTREZID", "ENTREZID_1W", "ENTREZID_4W")) {
  if (nm %in% names(structure)) structure[[nm]] <- clean_id(structure[[nm]])
}
structure$SYMBOL <- clean_symbol(structure$SYMBOL)
for (nm in c("logFC_1W", "FDR_1W", "logFC_4W", "FDR_4W")) {
  structure[[nm]] <- safe_num(structure[[nm]])
}
structure$category <- as.character(structure$category)
if (!"direction_persistent" %in% names(structure)) {
  # The Step07 structure table stores the shared direction as direction_1W;
  # for the consistent persistent category this is also the persistent direction.
  structure$direction_persistent <- ifelse(
    structure$category == "persistent_direction_consistent",
    structure$direction_1W,
    NA_character_
  )
}

target_classes <- c("persistent", "1-week-only", "4-week-only")
class_source <- c(
  persistent = "persistent_direction_consistent",
  `1-week-only` = "1W_only",
  `4-week-only` = "4W_only"
)

# Derive the same two-timepoint eligibility universe used by the primary ORA.
de1$ENTREZID <- clean_id(de1$ENTREZID)
de4$ENTREZID <- clean_id(de4$ENTREZID)
common_eligibility <- sort(intersect(na.omit(unique(de1$ENTREZID)),
                                     na.omit(unique(de4$ENTREZID))))
if (length(common_eligibility) != 14184L) {
  warning("Common eligibility universe is ", length(common_eligibility),
          " rather than the audited 14,184 genes.")
}

# Direction-stratified query lists and representative genes.
gene_lists <- list()
gene_rows <- list()

for (tc in target_classes) {
  src <- class_source[[tc]]
  x <- structure[structure$category == src, , drop = FALSE]
  if (tc == "persistent") {
    x$direction <- ifelse(x$direction_persistent == "Up", "Up", "Down")
    x$effect_timepoint <- "both"
    x$effect_logFC <- rowMeans(cbind(x$logFC_1W, x$logFC_4W), na.rm = TRUE)
    x$effect_FDR <- pmax(x$FDR_1W, x$FDR_4W)
    x$rank_score <- rowMeans(cbind(
      -log10(pmax(x$FDR_1W, .Machine$double.xmin)) * abs(x$logFC_1W),
      -log10(pmax(x$FDR_4W, .Machine$double.xmin)) * abs(x$logFC_4W)
    ), na.rm = TRUE)
  } else if (tc == "1-week-only") {
    x$direction <- ifelse(x$direction_1W == "Up", "Up", "Down")
    x$effect_timepoint <- "1W"
    x$effect_logFC <- x$logFC_1W
    x$effect_FDR <- x$FDR_1W
    x$rank_score <- -log10(pmax(x$FDR_1W, .Machine$double.xmin)) * abs(x$logFC_1W)
  } else {
    x$direction <- ifelse(x$direction_4W == "Up", "Up", "Down")
    x$effect_timepoint <- "4W"
    x$effect_logFC <- x$logFC_4W
    x$effect_FDR <- x$FDR_4W
    x$rank_score <- -log10(pmax(x$FDR_4W, .Machine$double.xmin)) * abs(x$logFC_4W)
  }
  for (dr in c("Up", "Down")) {
    key <- paste(tc, dr, sep = "__")
    y <- x[x$direction == dr & !is.na(x$ENTREZID), , drop = FALSE]
    y <- y[!duplicated(y$ENTREZID), , drop = FALSE]
    gene_lists[[key]] <- sort(unique(y$ENTREZID))
    # Up to ten representatives per direction, ranked using the original DE
    # significance/effect-size convention (mean across both time points for persistent).
    y <- y[order(-y$rank_score, y$effect_FDR, -abs(y$effect_logFC), y$SYMBOL), , drop = FALSE]
    nrep <- min(10L, nrow(y))
    if (nrep > 0L) {
      yy <- y[seq_len(nrep), , drop = FALSE]
      gene_rows[[key]] <- data.frame(
        temporal_class = tc,
        direction = dr,
        gene_rank = seq_len(nrep),
        gene_id = yy$ENTREZID,
        gene_symbol = yy$SYMBOL,
        effect_timepoint = yy$effect_timepoint,
        logFC_1W = yy$logFC_1W,
        FDR_1W = yy$FDR_1W,
        logFC_4W = yy$logFC_4W,
        FDR_4W = yy$FDR_4W,
        ranking_score = yy$rank_score,
        stringsAsFactors = FALSE
      )
    }
  }
}

representative_genes <- do.call(rbind, gene_rows)
rownames(representative_genes) <- NULL
write_csv(representative_genes,
          file.path(out_dir, "Supplementary_Table_temporal_programs_representative_genes.csv"))

# Prepare database-specific tested pathways using the exact existing membership
# snapshot and the same 10-500-gene filter as the persistent ORA.
membership$gene_id <- clean_id(membership$gene_id)
membership$pathway <- as.character(membership$pathway)
membership$pathway_name <- as.character(membership$pathway_name)
membership$database <- as.character(membership$database)

db_objects <- list(
  Reactome = list(membership = membership[membership$database == "Reactome", , drop = FALSE],
                  universe = clean_id(universe_R$gene_id)),
  KEGG = list(membership = membership[membership$database == "KEGG", , drop = FALSE],
              universe = clean_id(universe_K$gene_id))
)

tested_pathways <- list()
for (db in names(db_objects)) {
  mm <- db_objects[[db]]$membership
  u <- sort(unique(na.omit(db_objects[[db]]$universe)))
  split_members <- split(mm$gene_id, mm$pathway)
  split_members <- lapply(split_members, function(z) sort(unique(intersect(na.omit(z), u))))
  sizes <- lengths(split_members)
  keep <- names(split_members)[sizes >= 10L & sizes <= 500L]
  pp <- split_members[keep]
  meta <- unique(mm[, c("pathway", "pathway_id", "pathway_name", "database_release"), drop = FALSE])
  meta <- meta[match(names(pp), meta$pathway), , drop = FALSE]
  meta$pathway_size_in_universe <- lengths(pp)
  tested_pathways[[db]] <- list(pathways = pp, metadata = meta, universe = u)
}

# One-sided hypergeometric ORA for each of six temporal-class-by-direction lists.
ora_rows <- list()
row_i <- 0L
for (key in names(gene_lists)) {
  bits <- strsplit(key, "__", fixed = TRUE)[[1]]
  tc <- bits[1]
  dr <- bits[2]
  for (db in names(tested_pathways)) {
    obj <- tested_pathways[[db]]
    query <- sort(unique(intersect(gene_lists[[key]], obj$universe)))
    N <- length(obj$universe)
    n <- length(query)
    for (i in seq_along(obj$pathways)) {
      members <- obj$pathways[[i]]
      K <- length(members)
      overlap <- intersect(query, members)
      k <- length(overlap)
      p <- if (n == 0L || k == 0L) 1 else stats::phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
      row_i <- row_i + 1L
      ov_sym <- membership$gene_symbol[match(overlap, membership$gene_id)]
      ov_sym <- ov_sym[!is.na(ov_sym) & ov_sym != ""]
      ora_rows[[row_i]] <- data.frame(
        temporal_class = tc,
        direction = dr,
        database = db,
        database_release = obj$metadata$database_release[i],
        pathway_id = obj$metadata$pathway_id[i],
        pathway = obj$metadata$pathway[i],
        pathway_name = obj$metadata$pathway_name[i],
        universe_size = N,
        query_size_in_universe = n,
        pathway_size_in_universe = K,
        overlap_size = k,
        expected_overlap = n * K / N,
        enrichment_ratio = if (n * K == 0) NA_real_ else (k / n) / (K / N),
        p_value = p,
        overlap_gene_ids = paste(overlap, collapse = ";"),
        overlap_gene_symbols = paste(ov_sym, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
}
ora <- do.call(rbind, ora_rows)
rownames(ora) <- NULL
ora$FDR_BH_joint_Reactome_KEGG <- NA_real_
for (key in unique(paste(ora$temporal_class, ora$direction, sep = "__"))) {
  idx <- which(paste(ora$temporal_class, ora$direction, sep = "__") == key)
  ora$FDR_BH_joint_Reactome_KEGG[idx] <- stats::p.adjust(ora$p_value[idx], method = "BH")
}
ora$FDR_status <- ifelse(ora$FDR_BH_joint_Reactome_KEGG < 0.05,
                         "FDR_supported",
                         ifelse(ora$p_value < 0.05, "nominal_only", "not_significant"))
ora <- ora[order(match(ora$temporal_class, target_classes),
                 match(ora$direction, c("Up", "Down")),
                 match(ora$database, c("Reactome", "KEGG")),
                 ora$FDR_BH_joint_Reactome_KEGG, ora$p_value, ora$pathway_name), , drop = FALSE]
write_csv(ora, file.path(out_dir, "Supplementary_Table_temporal_programs_all_ORA_results.csv"))

# Display up to five FDR-supported pathways per database and direction. If a
# stratum has no FDR-supported pathway, display its top nominal pathways and
# label them explicitly rather than implying FDR support.
top_rows <- list()
ti <- 0L
for (tc in target_classes) for (dr in c("Up", "Down")) for (db in c("Reactome", "KEGG")) {
  z <- ora[ora$temporal_class == tc & ora$direction == dr & ora$database == db, , drop = FALSE]
  supported <- z[z$FDR_status == "FDR_supported", , drop = FALSE]
  if (nrow(supported) > 0L) {
    z <- supported[order(supported$FDR_BH_joint_Reactome_KEGG,
                        -supported$enrichment_ratio, supported$pathway_name), , drop = FALSE]
    basis <- "FDR-supported"
  } else {
    z <- z[order(z$p_value, -z$enrichment_ratio, z$pathway_name), , drop = FALSE]
    basis <- "top nominal (no FDR-supported pathway in this stratum)"
  }
  z <- head(z, 5L)
  if (nrow(z) > 0L) {
    z$display_rank <- seq_len(nrow(z))
    z$display_basis <- basis
    ti <- ti + 1L
    top_rows[[ti]] <- z
  }
}
top_pathways <- do.call(rbind, top_rows)
rownames(top_pathways) <- NULL
top_pathways <- top_pathways[, c(
  "temporal_class", "direction", "database", "display_rank", "display_basis",
  "pathway_id", "pathway_name", "pathway_size_in_universe", "query_size_in_universe",
  "overlap_size", "expected_overlap", "enrichment_ratio", "p_value",
  "FDR_BH_joint_Reactome_KEGG", "overlap_gene_symbols"
)]
write_csv(top_pathways,
          file.path(out_dir, "Supplementary_Table_temporal_programs_top_pathways.csv"))

# Compact class-level summary for the supplement and for manuscript drafting.
summary_rows <- list()
for (tc in target_classes) {
  src <- class_source[[tc]]
  x <- structure[structure$category == src, , drop = FALSE]
  direction_vector <- if (tc == "persistent") x$direction_persistent else if (tc == "1-week-only") x$direction_1W else x$direction_4W
  n_up <- sum(direction_vector == "Up", na.rm = TRUE)
  n_down <- sum(direction_vector == "Down", na.rm = TRUE)
  pathway_text <- function(dr, db) {
    z <- top_pathways[top_pathways$temporal_class == tc & top_pathways$direction == dr & top_pathways$database == db, , drop = FALSE]
    if (nrow(z) == 0L) return("None")
    paste(sprintf("%s [FDR=%s]", z$pathway_name, format.pval(z$FDR_BH_joint_Reactome_KEGG, digits = 3, eps = 1e-300)), collapse = "; ")
  }
  gene_text <- function(dr) {
    z <- representative_genes[representative_genes$temporal_class == tc & representative_genes$direction == dr, , drop = FALSE]
    if (nrow(z) == 0L) return("None")
    paste(z$gene_symbol, collapse = ", ")
  }
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    temporal_class = tc,
    strict_gene_count = nrow(x),
    up_count = n_up,
    down_count = n_down,
    representative_up_genes = gene_text("Up"),
    representative_down_genes = gene_text("Down"),
    top_Reactome_up = pathway_text("Up", "Reactome"),
    top_Reactome_down = pathway_text("Down", "Reactome"),
    top_KEGG_up = pathway_text("Up", "KEGG"),
    top_KEGG_down = pathway_text("Down", "KEGG"),
    stringsAsFactors = FALSE
  )
}
summary_table <- do.call(rbind, summary_rows)
write_csv(summary_table,
          file.path(out_dir, "Supplementary_Table_temporal_programs_summary.csv"))

# Counts and audit information.
audit <- data.frame(
  metric = c(
    "strict_DE_threshold", "persistent_definition", "one_week_only_definition",
    "four_week_only_definition", "common_eligibility_universe_genes",
    "Reactome_annotated_universe_genes", "KEGG_annotated_universe_genes",
    "Reactome_tested_pathways", "KEGG_tested_pathways", "pooled_tested_pathways_per_list",
    "BH_family", "representative_genes_per_class_direction"
  ),
  value = c(
    "FDR < 0.05 and |logFC| > 1",
    "strict DE at 1W and 4W with consistent direction",
    "strict DE at 1W and not strict DE at 4W",
    "strict DE at 4W and not strict DE at 1W",
    length(common_eligibility), length(tested_pathways$Reactome$universe),
    length(tested_pathways$KEGG$universe), length(tested_pathways$Reactome$pathways),
    length(tested_pathways$KEGG$pathways), length(tested_pathways$Reactome$pathways) + length(tested_pathways$KEGG$pathways),
    "BH within each temporal class × direction across pooled Reactome + KEGG pathways",
    "up to 10 per direction"
  ), stringsAsFactors = FALSE
)
write_csv(audit, file.path(out_dir, "Supplementary_Table_temporal_programs_QA.csv"))

session_lines <- tryCatch(capture.output(sessionInfo()), error = function(e) {
  c(paste0("sessionInfo() unavailable: ", conditionMessage(e)),
    paste0("R.version: ", R.version.string))
})
writeLines(session_lines, file.path(out_dir, "Supplementary_Table_temporal_programs_sessionInfo.txt"), useBytes = TRUE)
cat("Wrote temporal program supplement to", out_dir, "\n")
cat("Classes:", paste(target_classes, collapse = ", "), "\n")
cat("Gene counts:\n")
print(summary_table[, c("temporal_class", "strict_gene_count", "up_count", "down_count")])
cat("ORA rows:", nrow(ora), "\n")
cat("FDR-supported rows:", sum(ora$FDR_status == "FDR_supported"), "\n")
