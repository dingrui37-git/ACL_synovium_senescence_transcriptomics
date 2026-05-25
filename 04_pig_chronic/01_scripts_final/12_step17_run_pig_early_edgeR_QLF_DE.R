# ============================================================
# Step 17: pig-early differential expression (edgeR QLF)
# Project root: E:/R/ACLsenescence2
#
# Purpose:
#   - use Step 16 gene-level counts
#   - run pig early DE for:
#       1) ACLT_untreated_t7  vs CON_t0
#       2) ACLT_untreated_t28 vs CON_t0
#   - method: edgeR glmQLFit + glmQLFTest
#   - save full tables / summaries / workspace into rebuild_submission
# ============================================================

setwd("E:/R/ACLsenescence2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rebuild_root <- file.path(project_root, "rebuild_submission")

pig_dir <- file.path(rebuild_root, "02_pig_early")
scripts_dir <- file.path(pig_dir, "scripts")
objects_dir <- file.path(pig_dir, "objects")
tables_dir <- file.path(pig_dir, "tables")
logs_dir <- file.path(pig_dir, "logs")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is not installed. Please install/load edgeR before running Step 17.")
}

normalize_slash <- function(x) {
  gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
}

safe_parse_gtf_gene_map <- function(gtf_path) {
  lines <- tryCatch(
    readLines(gzfile(gtf_path, open = "rt"), warn = FALSE),
    error = function(e1) {
      tryCatch(readLines(gtf_path, warn = FALSE), error = function(e2) character(0))
    }
  )

  if (length(lines) == 0) {
    return(data.frame(gene_id = character(0), gene_name = character(0), stringsAsFactors = FALSE))
  }

  lines <- lines[!startsWith(lines, "#")]
  gene_lines <- lines[grepl("\tgene\t", lines, fixed = TRUE)]

  if (length(gene_lines) == 0) {
    return(data.frame(gene_id = character(0), gene_name = character(0), stringsAsFactors = FALSE))
  }

  fields <- strsplit(gene_lines, "\t", fixed = TRUE)
  attr_col <- vapply(fields, function(x) if (length(x) >= 9) x[9] else "", character(1))

  gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attr_col, perl = TRUE)
  gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", attr_col, perl = TRUE)

  gene_id[gene_id == attr_col] <- NA_character_
  gene_name[gene_name == attr_col] <- NA_character_
  gene_name[is.na(gene_name) | gene_name == ""] <- gene_id[is.na(gene_name) | gene_name == ""]

  out <- data.frame(
    gene_id = gene_id,
    gene_name = gene_name,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene_id) & out$gene_id != "", , drop = FALSE]
  out <- out[!duplicated(out$gene_id), , drop = FALSE]
  out
}

count_file <- file.path(tables_dir, "step16_pig_early_gene_count_matrix.csv")
manifest_file <- file.path(tables_dir, "step09_pig_early_core_sample_fastq_manifest.csv")
canon_file <- file.path(tables_dir, "step14_pig_early_canonical_reference_manifest.csv")

needed <- c(count_file, manifest_file, canon_file)
if (!all(file.exists(needed))) {
  stop("Required Step 09/14/16 files are missing. Please complete Step 16 first.")
}

count_df <- read.csv(count_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8", check.names = FALSE)
manifest_df <- read.csv(manifest_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
canon_df <- read.csv(canon_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

gtf_rel <- canon_df$relative_path[canon_df$item == "annotation_gtf"]
if (length(gtf_rel) != 1 || is.na(gtf_rel)) {
  stop("Canonical GTF was not found in Step 14 manifest.")
}
gtf_abs <- normalize_slash(file.path(project_root, gtf_rel))
if (!file.exists(gtf_abs)) {
  stop("Canonical GTF does not exist on disk: ", gtf_abs)
}

gene_map_df <- safe_parse_gtf_gene_map(gtf_abs)

if (!("gene_id" %in% colnames(count_df))) {
  stop("Step 16 count matrix must contain a gene_id column.")
}

sample_ids <- manifest_df$sample_id
if (!all(sample_ids %in% colnames(count_df))) {
  missing_samples <- sample_ids[!sample_ids %in% colnames(count_df)]
  stop("Some manifest samples are missing in Step 16 count matrix: ", paste(missing_samples, collapse = ", "))
}

count_mat <- as.matrix(count_df[, sample_ids, drop = FALSE])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- count_df$gene_id

meta_df <- manifest_df[, c("sample_id", "core_group")]
meta_df$sample_id <- as.character(meta_df$sample_id)
meta_df$core_group <- as.character(meta_df$core_group)

output_t7 <- file.path(tables_dir, "step17_pig_early_DE_t7_vs_CON_QLF.csv")
output_t28 <- file.path(tables_dir, "step17_pig_early_DE_t28_vs_CON_QLF.csv")
output_summary <- file.path(tables_dir, "step17_pig_early_DE_run_summary.csv")
output_filter <- file.path(tables_dir, "step17_pig_early_DE_filter_stats.csv")
output_gene_map <- file.path(tables_dir, "step17_pig_early_gene_id_name_map.csv")
output_workspace <- file.path(objects_dir, "step17_pig_early_edgeR_QLF_DE_workspace.RData")
output_log <- file.path(logs_dir, "step17_pig_early_edgeR_QLF_DE_log.txt")

step17_targets <- c(
  output_t7, output_t28, output_summary, output_filter,
  output_gene_map, output_workspace, output_log
)

cleanup_outputs <- function() {
  existing <- step17_targets[file.exists(step17_targets)]
  if (length(existing) > 0) unlink(existing, force = TRUE)
}

cleanup_outputs()

run_one_contrast <- function(target_group, ref_group = "CON_t0") {
  keep_samples <- meta_df$core_group %in% c(ref_group, target_group)
  sub_meta <- meta_df[keep_samples, , drop = FALSE]
  sub_counts <- count_mat[, sub_meta$sample_id, drop = FALSE]

  group <- factor(sub_meta$core_group, levels = c(ref_group, target_group))
  design <- model.matrix(~ group)

  dge <- edgeR::DGEList(counts = sub_counts, samples = sub_meta)
  keep_genes <- edgeR::filterByExpr(dge, design = design)

  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)

  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  tab <- data.frame(
    gene_id = rownames(tab),
    as.data.frame(tab, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )

  if (nrow(gene_map_df) > 0) {
    m_idx <- match(tab$gene_id, gene_map_df$gene_id)
    tab$gene_name <- gene_map_df$gene_name[m_idx]
  } else {
    tab$gene_name <- NA_character_
  }

  tab$contrast <- paste0(target_group, "_vs_", ref_group)
  tab$target_group <- target_group
  tab$reference_group <- ref_group
  tab$significant_FDR05 <- tab$FDR < 0.05
  tab$strict_sig <- tab$FDR < 0.05 & abs(tab$logFC) > 1
  tab$direction <- ifelse(
    tab$logFC > 0, paste0("Up_in_", target_group),
    ifelse(tab$logFC < 0, paste0("Down_in_", target_group), "No_change")
  )

  summary_row <- data.frame(
    contrast = paste0(target_group, "_vs_", ref_group),
    target_group = target_group,
    reference_group = ref_group,
    n_genes_before_filter = nrow(sub_counts),
    n_genes_after_filter = nrow(dge$counts),
    n_samples = ncol(sub_counts),
    n_target = sum(group == target_group),
    n_reference = sum(group == ref_group),
    n_FDR05 = sum(tab$FDR < 0.05, na.rm = TRUE),
    n_strict = sum(tab$strict_sig, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    dge = dge,
    design = design,
    fit = fit,
    qlf = qlf,
    table = tab,
    summary = summary_row
  )
}

start_time <- Sys.time()

res_t7 <- tryCatch(
  run_one_contrast("ACLT_untreated_t7", "CON_t0"),
  error = function(e) {
    cleanup_outputs()
    stop("Step 17 failed for t7 vs control. Error: ", conditionMessage(e))
  }
)

res_t28 <- tryCatch(
  run_one_contrast("ACLT_untreated_t28", "CON_t0"),
  error = function(e) {
    cleanup_outputs()
    stop("Step 17 failed for t28 vs control. Error: ", conditionMessage(e))
  }
)

end_time <- Sys.time()
elapsed_min <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)

run_summary_df <- rbind(res_t7$summary, res_t28$summary)
run_summary_df$elapsed_minutes_total <- elapsed_min

filter_stats_df <- data.frame(
  metric = c(
    "project_root",
    "count_matrix_file",
    "gtf_relative",
    "n_input_genes_total",
    "n_samples_total",
    "t7_genes_after_filter",
    "t28_genes_after_filter",
    "elapsed_minutes_total"
  ),
  value = c(
    project_root,
    count_file,
    gtf_rel,
    nrow(count_mat),
    ncol(count_mat),
    nrow(res_t7$dge$counts),
    nrow(res_t28$dge$counts),
    elapsed_min
  ),
  stringsAsFactors = FALSE
)

write.csv(res_t7$table, output_t7, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(res_t28$table, output_t28, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(run_summary_df, output_summary, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(filter_stats_df, output_filter, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(gene_map_df, output_gene_map, row.names = FALSE, fileEncoding = "UTF-8")

save(
  count_mat,
  meta_df,
  gene_map_df,
  res_t7,
  res_t28,
  run_summary_df,
  filter_stats_df,
  file = output_workspace
)

log_lines <- c(
  "Step 17 completed successfully.",
  paste("Project root:", project_root),
  paste("Input count matrix:", count_file),
  paste("GTF:", gtf_abs),
  paste("Total genes in Step 16 count matrix:", nrow(count_mat)),
  paste("Total samples:", ncol(count_mat)),
  paste("t7 genes after filter:", nrow(res_t7$dge$counts)),
  paste("t28 genes after filter:", nrow(res_t28$dge$counts)),
  paste("Elapsed minutes total:", elapsed_min)
)

writeLines(log_lines, output_log, useBytes = TRUE)

writeLines(
  c(
    "# Step 17 run script archived automatically.",
    "setwd(\"E:/R/ACLsenescence2\")"
  ),
  con = file.path(scripts_dir, "step17_run_pig_early_edgeR_QLF_DE.R"),
  useBytes = TRUE
)

writeLines(
  c(
    "# Step 17 check script archived automatically.",
    "setwd(\"E:/R/ACLsenescence2\")"
  ),
  con = file.path(scripts_dir, "step17_check_pig_early_edgeR_QLF_DE.R"),
  useBytes = TRUE
)

message("Step 17 finished successfully.")
