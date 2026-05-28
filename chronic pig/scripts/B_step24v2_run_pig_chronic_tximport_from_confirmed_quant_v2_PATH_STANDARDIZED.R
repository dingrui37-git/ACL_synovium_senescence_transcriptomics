# ============================================================
# Step 24 v2: pig-chronic tx2gene + tximport from confirmed quant folder
# Project root: E:/R/ACLsenescence2
#
# Fix:
#   - quant_dir is standardized to rebuild_submission/raw data/GSE228848_synovium_quant
#     to match Step23/CHRONQC/CHRONLOCK audit scripts.
#   - chronic Salmon files in this folder are named like:
#       *_quant.sf.txt.gz
#     rather than plain quant.sf / quant.sf.gz
#   - file matching now supports:
#       quant.sf
#       quant.sf.gz
#       quant.sf.txt
#       quant.sf.txt.gz
# ============================================================

setwd("E:/R/ACLsenescence2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
rebuild_root <- file.path(project_root, "rebuild_submission")
chronic_dir <- file.path(rebuild_root, "03_pig_chronic")

scripts_dir <- file.path(chronic_dir, "scripts")
objects_dir <- file.path(chronic_dir, "objects")
tables_dir <- file.path(chronic_dir, "tables")
logs_dir <- file.path(chronic_dir, "logs")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("tximport", quietly = TRUE)) {
  stop("tximport is not installed. Please install/load tximport before running Step 24 v2.")
}

normalize_slash <- function(x) {
  gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
}

extract_attr <- function(attr_vec, key) {
  pat <- paste0('.*', key, ' "([^"]+)".*')
  out <- sub(pat, "\\1", attr_vec, perl = TRUE)
  out[out == attr_vec] <- NA_character_
  out
}

quant_dir <- file.path(rebuild_root, "raw data", "GSE228848_synovium_quant")
gtf_file <- file.path(project_root, "reference", "Sus_scrofa_Ensembl115", "Sus_scrofa.Sscrofa11.1.115.gtf.gz")

if (!dir.exists(quant_dir)) {
  stop("Confirmed chronic quant folder does not exist: ", quant_dir)
}
if (!file.exists(gtf_file)) {
  stop("Local pig Ensembl 115 GTF not found: ", gtf_file)
}

all_files <- list.files(
  quant_dir,
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

all_files <- normalize_slash(all_files)
all_files <- all_files[file.exists(all_files)]

is_salmon_quant_name <- function(x) {
  grepl("quant\\.sf(\\.txt)?(\\.gz)?$", basename(x), ignore.case = TRUE)
}

quant_files <- all_files[vapply(all_files, is_salmon_quant_name, logical(1))]

if (length(quant_files) == 0) {
  stop(
    "No Salmon quant files were found in: ", quant_dir,
    ". Expected names like quant.sf, quant.sf.gz, or *_quant.sf.txt.gz"
  )
}

extract_quant_sample_id <- function(path) {
  b <- basename(path)
  b <- sub("\\.gz$", "", b, ignore.case = TRUE)
  b <- sub("\\.txt$", "", b, ignore.case = TRUE)
  b <- sub("\\.sf$", "", b, ignore.case = TRUE)
  b <- sub("_quant$", "", b, ignore.case = TRUE)
  b
}

sample_ids <- vapply(quant_files, extract_quant_sample_id, character(1))

quant_inventory_df <- data.frame(
  sample_id = sample_ids,
  absolute_path = quant_files,
  relative_path = sub(paste0("^", normalize_slash(project_root), "/?"), "", quant_files),
  filename = basename(quant_files),
  parent_dir = basename(dirname(quant_files)),
  size_bytes = file.info(quant_files)$size,
  stringsAsFactors = FALSE
)
quant_inventory_df$size_mb <- round(quant_inventory_df$size_bytes / 1024^2, 3)
quant_inventory_df <- quant_inventory_df[order(quant_inventory_df$sample_id, quant_inventory_df$relative_path), , drop = FALSE]

if (anyDuplicated(quant_inventory_df$sample_id)) {
  dup_ids <- unique(quant_inventory_df$sample_id[duplicated(quant_inventory_df$sample_id)])
  stop("Duplicated sample IDs detected after filename parsing: ", paste(dup_ids, collapse = ", "))
}

# ---------- build tx2gene from local GTF ----------
gtf_df <- tryCatch(
  read.delim(
    gzfile(gtf_file),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "#",
    stringsAsFactors = FALSE,
    fill = TRUE
  ),
  error = function(e1) {
    read.delim(
      gtf_file,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "#",
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  }
)

if (ncol(gtf_df) < 9) {
  stop("GTF parsing failed: fewer than 9 columns detected.")
}

colnames(gtf_df)[1:9] <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

tx_rows <- gtf_df[gtf_df$feature %in% c("transcript", "exon"), c("feature", "attribute"), drop = FALSE]
tx_rows$transcript_id <- extract_attr(tx_rows$attribute, "transcript_id")
tx_rows$gene_id <- extract_attr(tx_rows$attribute, "gene_id")

tx2gene_df <- unique(tx_rows[, c("transcript_id", "gene_id"), drop = FALSE])
tx2gene_df <- tx2gene_df[!is.na(tx2gene_df$transcript_id) & tx2gene_df$transcript_id != "", , drop = FALSE]
tx2gene_df <- tx2gene_df[!is.na(tx2gene_df$gene_id) & tx2gene_df$gene_id != "", , drop = FALSE]
tx2gene_df <- tx2gene_df[order(tx2gene_df$transcript_id), , drop = FALSE]

if (nrow(tx2gene_df) == 0) {
  stop("tx2gene construction failed: no transcript-gene pairs extracted from GTF.")
}

gene_rows <- gtf_df[gtf_df$feature == "gene", "attribute", drop = FALSE]
gene_annotation_df <- data.frame(
  gene_id = extract_attr(gene_rows$attribute, "gene_id"),
  gene_name = extract_attr(gene_rows$attribute, "gene_name"),
  gene_biotype = extract_attr(gene_rows$attribute, "gene_biotype"),
  stringsAsFactors = FALSE
)
gene_annotation_df$gene_name[is.na(gene_annotation_df$gene_name) | gene_annotation_df$gene_name == ""] <- gene_annotation_df$gene_id[is.na(gene_annotation_df$gene_name) | gene_annotation_df$gene_name == ""]
gene_annotation_df <- unique(gene_annotation_df)
gene_annotation_df <- gene_annotation_df[!is.na(gene_annotation_df$gene_id) & gene_annotation_df$gene_id != "", , drop = FALSE]
gene_annotation_df <- gene_annotation_df[!duplicated(gene_annotation_df$gene_id), , drop = FALSE]
gene_annotation_df <- gene_annotation_df[order(gene_annotation_df$gene_id), , drop = FALSE]

# ---------- tximport ----------
files_named <- quant_inventory_df$absolute_path
names(files_named) <- quant_inventory_df$sample_id

txi <- tximport::tximport(
  files = files_named,
  type = "salmon",
  tx2gene = tx2gene_df,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "no"
)

counts_mat <- txi$counts
abundance_mat <- txi$abundance
length_mat <- txi$length

counts_df <- data.frame(
  gene_id = rownames(counts_mat),
  as.data.frame(counts_mat, stringsAsFactors = FALSE, check.names = FALSE),
  stringsAsFactors = FALSE
)

abundance_df <- data.frame(
  gene_id = rownames(abundance_mat),
  as.data.frame(abundance_mat, stringsAsFactors = FALSE, check.names = FALSE),
  stringsAsFactors = FALSE
)

length_df <- data.frame(
  gene_id = rownames(length_mat),
  as.data.frame(length_mat, stringsAsFactors = FALSE, check.names = FALSE),
  stringsAsFactors = FALSE
)

run_summary_df <- data.frame(
  metric = c(
    "project_root",
    "quant_dir",
    "gtf_file",
    "tximport_type",
    "ignoreTxVersion",
    "countsFromAbundance",
    "n_quant_files",
    "n_unique_samples",
    "n_tx2gene_rows",
    "n_gene_annotation_rows",
    "counts_n_genes",
    "counts_n_samples",
    "abundance_n_genes",
    "abundance_n_samples",
    "length_n_genes",
    "length_n_samples"
  ),
  value = c(
    project_root,
    quant_dir,
    gtf_file,
    "salmon",
    TRUE,
    "no",
    nrow(quant_inventory_df),
    length(unique(quant_inventory_df$sample_id)),
    nrow(tx2gene_df),
    nrow(gene_annotation_df),
    nrow(counts_mat),
    ncol(counts_mat),
    nrow(abundance_mat),
    ncol(abundance_mat),
    nrow(length_mat),
    ncol(length_mat)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  quant_inventory_df,
  file.path(tables_dir, "step24_pig_chronic_quant_inventory.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  tx2gene_df,
  file.path(tables_dir, "step24_pig_chronic_tx2gene.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  gene_annotation_df,
  file.path(tables_dir, "step24_pig_chronic_gene_annotation_from_gtf.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  counts_df,
  file.path(tables_dir, "step24_pig_chronic_gene_level_estimated_counts_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  abundance_df,
  file.path(tables_dir, "step24_pig_chronic_gene_level_abundance_tpm_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  length_df,
  file.path(tables_dir, "step24_pig_chronic_gene_level_length_matrix.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  run_summary_df,
  file.path(tables_dir, "step24_pig_chronic_tximport_run_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

save(
  quant_inventory_df, tx2gene_df, gene_annotation_df,
  txi, counts_mat, abundance_mat, length_mat,
  counts_df, abundance_df, length_df, run_summary_df,
  file = file.path(objects_dir, "step24_pig_chronic_tximport_workspace.RData")
)

writeLines(
  c(
    "Step 24 v2 completed successfully.",
    paste("Project root:", project_root),
    paste("Quant directory:", quant_dir),
    paste("GTF file:", gtf_file),
    paste("Quant files:", nrow(quant_inventory_df)),
    paste("Unique samples:", length(unique(quant_inventory_df$sample_id))),
    paste("tx2gene rows:", nrow(tx2gene_df)),
    paste("Counts matrix dim:", paste(nrow(counts_mat), "x", ncol(counts_mat))),
    paste("Abundance matrix dim:", paste(nrow(abundance_mat), "x", ncol(abundance_mat))),
    paste("Length matrix dim:", paste(nrow(length_mat), "x", ncol(length_mat)))
  ),
  file.path(logs_dir, "step24_pig_chronic_tximport_log.txt"),
  useBytes = TRUE
)

writeLines(
  c(
    "# Step 24 v2 run script archived automatically.",
    "setwd(\"E:/R/ACLsenescence2\")"
  ),
  con = file.path(scripts_dir, "step24_run_pig_chronic_tximport_from_confirmed_quant_v2.R"),
  useBytes = TRUE
)

writeLines(
  c(
    "# Step 24 v2 check script archived automatically.",
    "setwd(\"E:/R/ACLsenescence2\")"
  ),
  con = file.path(scripts_dir, "step24_check_pig_chronic_tximport_from_confirmed_quant_v2.R"),
  useBytes = TRUE
)

message("Step 24 v2 finished successfully.")
