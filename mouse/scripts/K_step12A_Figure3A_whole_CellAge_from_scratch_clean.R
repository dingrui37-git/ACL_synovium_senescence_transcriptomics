# Step12A: Figure 3A whole CellAge projected score from scratch.
# Purpose:
# Generate Figure 3A using a fully reproducible from-scratch workflow:
# CellAge raw file -> gprofiler2::gorth human-to-mouse strict 1:1 ortholog mapping
# -> projected mouse CellAge gene set -> 48-sample TMM logCPM -> row-z -> score.
#
# Locked method:
# 1. Figure3A is a broad CellAge database-level projection, not the persistent ∩ CellAge overlap score.
# 2. CellAge genes are read from the raw CellAge file under 02_mouse_discovery/00_raw_data/cellAge.
# 3. Human-to-mouse ortholog mapping is performed from scratch using gprofiler2::gorth.
# 4. Strict 1:1 orthologs are retained: one human input has one mouse ortholog, and that mouse ortholog maps back to one human input within the gorth output.
# 5. Step04 cleaned count matrices for 1W and 4W are combined into one 48-sample discovery count matrix.
# 6. Genes with any missing count are removed before edgeR::DGEList; no imputation is performed.
# 7. TMM normalization is performed using edgeR::calcNormFactors().
# 8. logCPM is obtained using edgeR::cpm(log = TRUE, prior.count = 1).
# 9. Row-wise z-score is calculated across all 48 mouse discovery samples.
# 10. Sample-level CellAge scores are mean row-z scores for CellAge all / induces / inhibits.
# 11. Within each time point, ACLR vs Contra is tested using paired Wilcoxon test.

options(stringsAsFactors = FALSE)

## =========================
## 0. Paths
## =========================

project_root <- "E:/R/ACLsenescence2"
base_dir <- file.path(project_root, "rebuild_submission", "02_mouse_discovery")
data_raw_dir <- file.path(base_dir, "00_raw_data", "cellAge")

expr_1w_file <- file.path(base_dir, "02_expression_matrix", "step04_expr_1W.csv")
expr_4w_file <- file.path(base_dir, "02_expression_matrix", "step04_expr_4W.csv")
anno_1w_file <- file.path(base_dir, "01_metadata", "step04_anno_1W.csv")
anno_4w_file <- file.path(base_dir, "01_metadata", "step04_anno_4W.csv")

figure_dir <- file.path(base_dir, "06_figures", "Figure3")
table_dir  <- file.path(base_dir, "07_tables", "step12A_Figure3A_whole_CellAge_from_scratch")
log_dir    <- file.path(base_dir, "08_logs")
script_dir <- file.path(base_dir, "09_scripts")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

step_name <- "step12A_Figure3A_whole_CellAge_from_scratch"
log_file <- file.path(log_dir, paste0(step_name, "_log.txt"))

png_file <- file.path(figure_dir, "Figure3A_whole_CellAge_projected_scores_from_scratch.png")
pdf_file <- file.path(figure_dir, "Figure3A_whole_CellAge_projected_scores_from_scratch.pdf")

## =========================
## 1. Remove old Step12A from-scratch outputs only
## =========================

old_outputs <- c(
  list.files(table_dir, full.names = TRUE, recursive = TRUE),
  png_file,
  pdf_file,
  log_file,
  file.path(script_dir, "step12A_Figure3A_whole_CellAge_from_scratch.R")
)
unlink(old_outputs[file.exists(old_outputs)], recursive = TRUE, force = TRUE)

## =========================
## 2. Helpers
## =========================

sink(log_file, split = TRUE)

cat("===== STEP12A FIGURE3A WHOLE CELLAGE FROM SCRATCH =====\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Raw CellAge directory: ", data_raw_dir, "\n", sep = "")
cat("Expected raw CellAge file: csv/tsv/txt/xls/xlsx with CellAge/cellage in filename, e.g. cellage3.tsv.\n")
cat("Expression input 1W: ", expr_1w_file, "\n", sep = "")
cat("Expression input 4W: ", expr_4w_file, "\n", sep = "")

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    return(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE))
  }
  NA_character_
}

archive_current_script <- function(path) {
  src <- get_script_path()
  if (!is.na(src) && file.exists(src)) {
    writeLines(readLines(src, warn = FALSE, encoding = "UTF-8"), path, useBytes = TRUE)
  } else {
    writeLines("# Step12A from-scratch archive fallback.", path, useBytes = TRUE)
  }
}

safe_library <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  cat("Saved: ", path, "\n", sep = "")
}

pick_col <- function(df, candidates, required = TRUE, what = NULL) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Cannot find required column", if (!is.null(what)) paste0(" for ", what) else "",
         ". Checked: ", paste(candidates, collapse = ", "))
  }
  NA_character_
}

clean_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "NULL", "null", "na")] <- ""
  x
}

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", format(round(p, 3), nsmall = 3))
}

paired_wilcox <- function(df) {
  wide <- reshape(
    df[, c("mouse_id", "treatment", "score")],
    idvar = "mouse_id",
    timevar = "treatment",
    direction = "wide"
  )
  wide <- wide[complete.cases(wide[, c("score.Contra", "score.ACLR")]), , drop = FALSE]
  delta <- wide$score.ACLR - wide$score.Contra

  wt <- if (nrow(wide) >= 3) {
    suppressWarnings(wilcox.test(wide$score.ACLR, wide$score.Contra, paired = TRUE))
  } else {
    NULL
  }

  data.frame(
    n_pairs = nrow(wide),
    mean_contra = mean(wide$score.Contra, na.rm = TRUE),
    mean_aclr = mean(wide$score.ACLR, na.rm = TRUE),
    mean_delta_aclr_minus_contra = mean(delta, na.rm = TRUE),
    median_delta_aclr_minus_contra = median(delta, na.rm = TRUE),
    wilcox_p = if (is.null(wt)) NA_real_ else wt$p.value,
    p_label = fmt_p(if (is.null(wt)) NA_real_ else wt$p.value),
    stringsAsFactors = FALSE
  )
}

read_table_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext %in% c("tsv", "txt")) {
    return(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (ext %in% c("xls", "xlsx")) {
    safe_library("readxl")
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop("Unsupported file extension: ", path)
}

infer_cellage_version <- function(path) {
  bn <- basename(path)
  hit <- regmatches(bn, regexpr("(?i)cellage[_ -]*[0-9]+", bn, perl = TRUE))
  if (length(hit) == 1 && nchar(hit) > 0) return(hit)
  tools::file_path_sans_ext(bn)
}

file_audit_row <- function(path, label) {
  info <- file.info(path)
  data.frame(
    item = c(
      paste0(label, "_path"),
      paste0(label, "_basename"),
      paste0(label, "_inferred_version_from_filename"),
      paste0(label, "_size_bytes"),
      paste0(label, "_modified_time"),
      paste0(label, "_md5")
    ),
    value = c(
      normalizePath(path, winslash = "/", mustWork = FALSE),
      basename(path),
      infer_cellage_version(path),
      as.character(info$size),
      format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
      unname(tools::md5sum(path))
    ),
    stringsAsFactors = FALSE
  )
}

get_gprofiler_version_info_safe <- function() {
  out <- tryCatch(
    {
      if ("get_version_info" %in% getNamespaceExports("gprofiler2")) {
        tmp <- gprofiler2::get_version_info(organism = "hsapiens")
        paste(capture.output(print(tmp)), collapse = " | ")
      } else {
        "gprofiler2::get_version_info not available in this installed package version"
      }
    },
    error = function(e) paste0("gprofiler version info unavailable: ", conditionMessage(e))
  )
  out
}

read_count_matrix <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if ("gene_id" %in% colnames(df)) {
    row_ids <- as.character(df$gene_id)
    df$gene_id <- NULL
  } else {
    row_ids <- as.character(df[[1]])
    df[[1]] <- NULL
  }

  if ("gene_symbol" %in% colnames(df)) df$gene_symbol <- NULL
  if ("SYMBOL" %in% colnames(df)) df$SYMBOL <- NULL
  if ("symbol" %in% colnames(df)) df$symbol <- NULL

  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- row_ids
  mat
}

archive_current_script(file.path(script_dir, "step12A_Figure3A_whole_CellAge_from_scratch.R"))

safe_library("edgeR", bioc = TRUE)
safe_library("ggplot2")
safe_library("gprofiler2")
safe_library("clusterProfiler", bioc = TRUE)
safe_library("org.Mm.eg.db", bioc = TRUE)

gprofiler_package_version <- as.character(utils::packageVersion("gprofiler2"))
gprofiler_version_info <- get_gprofiler_version_info_safe()

cat("\nMapping tool version audit:\n")
cat("gprofiler2 package version: ", gprofiler_package_version, "\n", sep = "")
cat("gprofiler2 version info: ", gprofiler_version_info, "\n", sep = "")

## =========================
## 3. Load raw CellAge file
## =========================

cellage_candidates <- list.files(
  data_raw_dir,
  pattern = "(CellAge|cellage).*[.](csv|tsv|txt|xls|xlsx)$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(cellage_candidates) == 0) {
  stop("No CellAge file found under: ", data_raw_dir,
       "\nExpected a csv/tsv/txt/xls/xlsx file with CellAge in the filename.")
}

## Prefer files likely containing gene records, not old score outputs.
cellage_candidates <- cellage_candidates[!grepl("Figure3A|score|coverage|ortholog|mapping", basename(cellage_candidates), ignore.case = TRUE)]
if (length(cellage_candidates) == 0) {
  stop("Only derived CellAge files were found, but no raw CellAge file was found under: ", data_raw_dir)
}

cellage_file <- cellage_candidates[1]
cat("\nUsing raw CellAge file:\n", cellage_file, "\n", sep = "")

cellage_file_audit <- file_audit_row(cellage_file, "CellAge_raw_file")
cat("CellAge inferred version from filename: ", infer_cellage_version(cellage_file), "\n", sep = "")
cat("CellAge file modified time: ", cellage_file_audit$value[cellage_file_audit$item == "CellAge_raw_file_modified_time"], "\n", sep = "")
cat("CellAge file md5: ", cellage_file_audit$value[cellage_file_audit$item == "CellAge_raw_file_md5"], "\n", sep = "")
write_csv(cellage_file_audit, file.path(table_dir, "step12A_CellAge_raw_file_audit.csv"))

cellage_raw <- read_table_auto(cellage_file)
cat("CellAge raw table dimensions: ", paste(dim(cellage_raw), collapse = " x "), "\n", sep = "")

human_symbol_col <- pick_col(
  cellage_raw,
  c("Gene symbol", "Gene Symbol", "gene_symbol", "Symbol", "symbol", "Gene", "gene"),
  what = "CellAge human gene symbol"
)

effect_col <- pick_col(
  cellage_raw,
  c("Senescence Effect", "CellAge_effect", "senescence_effect", "Effect", "effect"),
  required = FALSE,
  what = "CellAge senescence effect"
)

cellage_clean <- cellage_raw
cellage_clean$CellAge_human_symbol <- clean_symbol(cellage_clean[[human_symbol_col]])

if (!is.na(effect_col)) {
  cellage_clean$CellAge_effect <- clean_symbol(cellage_clean[[effect_col]])
} else {
  cellage_clean$CellAge_effect <- "Unclear"
}

cellage_clean <- cellage_clean[cellage_clean$CellAge_human_symbol != "", , drop = FALSE]
cellage_clean$CellAge_effect <- ifelse(
  tolower(cellage_clean$CellAge_effect) == "induces",
  "Induces",
  ifelse(tolower(cellage_clean$CellAge_effect) == "inhibits", "Inhibits", "Unclear")
)

## Keep unique human gene + effect rows. This preserves genes annotated to different effects if present.
cellage_genes <- unique(cellage_clean[, c("CellAge_human_symbol", "CellAge_effect")])
colnames(cellage_genes) <- c("human_symbol", "CellAge_effect")

write_csv(cellage_clean, file.path(table_dir, "step12A_CellAge_raw_cleaned_used.csv"))
write_csv(cellage_genes, file.path(table_dir, "step12A_CellAge_human_genes_by_effect_used.csv"))

cat("\nCellAge human gene/effect summary:\n")
print(table(cellage_genes$CellAge_effect, useNA = "ifany"))
cat("Unique human symbols: ", length(unique(cellage_genes$human_symbol)), "\n", sep = "")

## =========================
## 4. From-scratch human-to-mouse gorth mapping and strict 1:1 filtering
## =========================

query_symbols <- sort(unique(cellage_genes$human_symbol))

cat("\nRunning gprofiler2::gorth human-to-mouse mapping from scratch...\n")
cat("Number of human query symbols: ", length(query_symbols), "\n", sep = "")

gorth_raw <- gprofiler2::gorth(
  query = query_symbols,
  source_organism = "hsapiens",
  target_organism = "mmusculus",
  mthreshold = Inf,
  filter_na = FALSE
)

write_csv(gorth_raw, file.path(table_dir, "step12A_gorth_human_to_mouse_raw.csv"))

cat("gorth raw output dimensions: ", paste(dim(gorth_raw), collapse = " x "), "\n", sep = "")
cat("gorth raw output columns: ", paste(colnames(gorth_raw), collapse = ", "), "\n", sep = "")

if (!all(c("input", "ortholog_name") %in% colnames(gorth_raw))) {
  stop("gorth output lacks required columns input and ortholog_name.")
}

gorth_valid <- gorth_raw[
  !is.na(gorth_raw$input) & gorth_raw$input != "" &
    !is.na(gorth_raw$ortholog_name) & gorth_raw$ortholog_name != "",
  ,
  drop = FALSE
]

## Count unique ortholog relationships within the returned table.
n_mouse_orthologs <- aggregate(
  ortholog_name ~ input,
  data = unique(gorth_valid[, c("input", "ortholog_name")]),
  FUN = function(x) length(unique(x))
)
colnames(n_mouse_orthologs) <- c("input", "n_mouse_orthologs")

n_human_inputs <- aggregate(
  input ~ ortholog_name,
  data = unique(gorth_valid[, c("input", "ortholog_name")]),
  FUN = function(x) length(unique(x))
)
colnames(n_human_inputs) <- c("ortholog_name", "n_human_inputs")

gorth_counted <- merge(gorth_valid, n_mouse_orthologs, by = "input", all.x = TRUE)
gorth_counted <- merge(gorth_counted, n_human_inputs, by = "ortholog_name", all.x = TRUE)

gorth_strict <- gorth_counted[
  gorth_counted$n_mouse_orthologs == 1 & gorth_counted$n_human_inputs == 1,
  ,
  drop = FALSE
]
gorth_strict <- gorth_strict[!duplicated(paste(gorth_strict$input, gorth_strict$ortholog_name, sep = "|||")), , drop = FALSE]

write_csv(gorth_counted, file.path(table_dir, "step12A_gorth_human_to_mouse_all_counted_pairs.csv"))
write_csv(gorth_strict, file.path(table_dir, "step12A_gorth_human_to_mouse_strict_1to1_pairs.csv"))

## Join CellAge effects back to strict mapping.
map_std <- merge(
  cellage_genes,
  gorth_strict,
  by.x = "human_symbol",
  by.y = "input",
  all.x = FALSE,
  all.y = FALSE
)

map_std$mouse_symbol <- clean_symbol(map_std$ortholog_name)
map_std$score_type <- ifelse(
  map_std$CellAge_effect == "Induces",
  "CellAge_induces",
  ifelse(map_std$CellAge_effect == "Inhibits", "CellAge_inhibits", "CellAge_all_only")
)

write_csv(map_std, file.path(table_dir, "step12A_CellAge_human_to_mouse_strict_1to1_map_from_scratch.csv"))

## =========================
## 5. Load expression and metadata, build 48-sample TMM logCPM
## =========================

required_files <- c(expr_1w_file, expr_4w_file, anno_1w_file, anno_4w_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

expr_1w <- read_count_matrix(expr_1w_file)
expr_4w <- read_count_matrix(expr_4w_file)

anno_1w <- read.csv(anno_1w_file, stringsAsFactors = FALSE)
anno_4w <- read.csv(anno_4w_file, stringsAsFactors = FALSE)

anno_1w$time_point <- "1W"
anno_4w$time_point <- "4W"
anno_1w$title <- anno_1w$sample_id
anno_4w$title <- anno_4w$sample_id

if (!("sex" %in% colnames(anno_1w)) && "Sex" %in% colnames(anno_1w)) anno_1w$sex <- anno_1w$Sex
if (!("sex" %in% colnames(anno_4w)) && "Sex" %in% colnames(anno_4w)) anno_4w$sex <- anno_4w$Sex

common_genes <- intersect(rownames(expr_1w), rownames(expr_4w))
expr_1w_common <- expr_1w[common_genes, , drop = FALSE]
expr_4w_common <- expr_4w[common_genes, , drop = FALSE]

expr_all <- cbind(expr_1w_common, expr_4w_common)

meta_all <- rbind(anno_1w, anno_4w)
meta_all <- meta_all[match(colnames(expr_all), meta_all$title), , drop = FALSE]

if (!all(colnames(expr_all) == meta_all$title)) {
  mismatch <- data.frame(expr_col = colnames(expr_all), meta_title = meta_all$title, stringsAsFactors = FALSE)
  write_csv(mismatch, file.path(table_dir, "step12A_expression_metadata_alignment_mismatch.csv"))
  stop("Expression matrix and metadata are not aligned.")
}

na_gene_idx <- apply(expr_all, 1, function(x) any(is.na(x)))
na_gene_ids <- rownames(expr_all)[na_gene_idx]
if (length(na_gene_ids) > 0) {
  write_csv(
    data.frame(gene_id = na_gene_ids, reason = "removed_before_DGEList_due_to_any_NA_count"),
    file.path(table_dir, "step12A_genes_removed_before_DGEList_due_to_NA_counts.csv")
  )
  expr_all <- expr_all[!na_gene_idx, , drop = FALSE]
}

if (any(is.na(expr_all))) stop("NA counts remain after removing rows with missing values.")

dge_all <- edgeR::DGEList(counts = expr_all)
dge_all <- edgeR::calcNormFactors(dge_all)
logcpm_all <- edgeR::cpm(dge_all, log = TRUE, prior.count = 1)

## =========================
## 6. Build genome-wide mouse symbol -> gene_id map
## =========================

## The Step04 expression matrices may only contain st_gene_id row identifiers and no SYMBOL column.
## Therefore, build a genome-wide gene_id -> mouse SYMBOL map by searching current project
## tables and, when needed, mapping Entrez IDs to SYMBOL using org.Mm.eg.db + clusterProfiler::bitr.

read_expression_annotation_direct <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  gene_id <- if ("gene_id" %in% colnames(df)) as.character(df$gene_id) else as.character(df[[1]])

  symbol_col <- NA_character_
  for (cc in c("gene_symbol", "SYMBOL", "symbol", "Symbol", "Gene symbol", "Gene Symbol")) {
    if (cc %in% colnames(df)) {
      symbol_col <- cc
      break
    }
  }

  if (is.na(symbol_col)) {
    return(data.frame(gene_id = gene_id, mouse_symbol = NA_character_, source = basename(path), stringsAsFactors = FALSE))
  }

  data.frame(
    gene_id = gene_id,
    mouse_symbol = clean_symbol(df[[symbol_col]]),
    source = basename(path),
    stringsAsFactors = FALSE
  )
}

guess_delim_reader <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  if (ext %in% c("tsv", "txt")) return(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  if (ext %in% c("xls", "xlsx")) {
    safe_library("readxl")
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE, check.names = FALSE))
  }
  NULL
}

build_symbol_map_from_candidate <- function(path, expr_gene_ids) {
  out_empty <- data.frame(gene_id = character(0), mouse_symbol = character(0), source = character(0), stringsAsFactors = FALSE)

  tab <- tryCatch(guess_delim_reader(path), error = function(e) NULL)
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0 || ncol(tab) < 2) return(out_empty)

  ## Remove zero-length or unnamed columns defensively.
  if (is.null(colnames(tab)) || any(is.na(colnames(tab)))) return(out_empty)
  colnames(tab) <- make.unique(as.character(colnames(tab)))

  ## Identify the column matching current expression gene IDs.
  match_counts <- vapply(seq_along(tab), function(j) {
    x <- tab[[j]]
    if (length(x) != nrow(tab)) return(0L)
    sum(as.character(x) %in% expr_gene_ids, na.rm = TRUE)
  }, integer(1))
  names(match_counts) <- colnames(tab)

  if (length(match_counts) == 0 || all(is.na(match_counts)) || max(match_counts, na.rm = TRUE) < 100) {
    return(out_empty)
  }

  id_col <- names(match_counts)[which.max(match_counts)]
  if (length(id_col) != 1 || is.na(id_col) || !(id_col %in% colnames(tab))) return(out_empty)
  if (length(tab[[id_col]]) != nrow(tab)) return(out_empty)

  ## Direct SYMBOL column if available.
  symbol_col <- NA_character_
  for (cc in c("gene_symbol", "SYMBOL", "symbol", "Symbol", "Gene symbol", "Gene Symbol", "mgi_symbol", "MGI_symbol")) {
    if (cc %in% colnames(tab) && length(tab[[cc]]) == nrow(tab)) {
      symbol_col <- cc
      break
    }
  }

  if (!is.na(symbol_col)) {
    gene_ids_tmp <- as.character(tab[[id_col]])
    symbols_tmp <- clean_symbol(tab[[symbol_col]])
    if (length(gene_ids_tmp) == length(symbols_tmp) && length(gene_ids_tmp) == nrow(tab)) {
      tmp <- data.frame(
        gene_id = gene_ids_tmp,
        mouse_symbol = symbols_tmp,
        source = rep(path, length(gene_ids_tmp)),
        stringsAsFactors = FALSE
      )
      tmp <- tmp[tmp$gene_id %in% expr_gene_ids & tmp$mouse_symbol != "", , drop = FALSE]
      if (nrow(tmp) > 0) return(unique(tmp))
    }
  }

  ## Otherwise try an Entrez ID column and map via org.Mm.eg.db.
  entrez_candidates <- c("Entrez ID", "EntrezID", "ENTREZID", "entrez_id", "entrez", "GeneID")
  entrez_col <- NA_character_
  for (cc in entrez_candidates) {
    if (cc %in% colnames(tab) && cc != id_col && length(tab[[cc]]) == nrow(tab)) {
      vals <- clean_symbol(tab[[cc]])
      if (length(vals) == nrow(tab) && sum(vals != "") > 100) {
        entrez_col <- cc
        break
      }
    }
  }

  if (!is.na(entrez_col)) {
    gene_ids_tmp <- as.character(tab[[id_col]])
    entrez_tmp <- clean_symbol(tab[[entrez_col]])
    if (length(gene_ids_tmp) != length(entrez_tmp) || length(gene_ids_tmp) != nrow(tab)) return(out_empty)

    tmp0 <- data.frame(
      gene_id = gene_ids_tmp,
      ENTREZID = entrez_tmp,
      stringsAsFactors = FALSE
    )
    tmp0 <- tmp0[tmp0$gene_id %in% expr_gene_ids & tmp0$ENTREZID != "", , drop = FALSE]
    if (nrow(tmp0) == 0) return(out_empty)

    bitr_map <- tryCatch(
      clusterProfiler::bitr(
        unique(tmp0$ENTREZID),
        fromType = "ENTREZID",
        toType = "SYMBOL",
        OrgDb = org.Mm.eg.db::org.Mm.eg.db
      ),
      error = function(e) NULL
    )
    if (is.null(bitr_map) || nrow(bitr_map) == 0) return(out_empty)

    tmp <- merge(tmp0, bitr_map, by = "ENTREZID", all.x = FALSE)
    tmp <- data.frame(
      gene_id = tmp$gene_id,
      mouse_symbol = clean_symbol(tmp$SYMBOL),
      source = rep(path, nrow(tmp)),
      stringsAsFactors = FALSE
    )
    tmp <- tmp[tmp$mouse_symbol != "", , drop = FALSE]
    if (nrow(tmp) > 0) return(unique(tmp))
  }

  out_empty
}

expr_gene_ids <- rownames(expr_all)

## First try direct annotations inside Step04 matrices.
anno_expr_1w <- read_expression_annotation_direct(expr_1w_file)
anno_expr_4w <- read_expression_annotation_direct(expr_4w_file)
expr_symbol_map <- unique(rbind(anno_expr_1w, anno_expr_4w))
expr_symbol_map <- expr_symbol_map[
  !is.na(expr_symbol_map$mouse_symbol) &
    expr_symbol_map$mouse_symbol != "" &
    expr_symbol_map$gene_id %in% expr_gene_ids,
  ,
  drop = FALSE
]

## If Step04 matrices lack symbols, search project/raw files for genome-wide annotation.
annotation_sources_used <- data.frame(source = character(0), n_mapped_rows = integer(0), stringsAsFactors = FALSE)

if (nrow(expr_symbol_map) == 0) {
  candidate_files <- unique(c(
    list.files(file.path(base_dir, "02_expression_matrix"), pattern = "[.](csv|tsv|txt|xls|xlsx)$", recursive = TRUE, full.names = TRUE),
    list.files(file.path(base_dir, "07_tables"), pattern = "[.](csv|tsv|txt|xls|xlsx)$", recursive = TRUE, full.names = TRUE),
    list.files(file.path(base_dir, "00_raw_data"), pattern = "[.](csv|tsv|txt|xls|xlsx)$", recursive = TRUE, full.names = TRUE),
    list.files(data_raw_dir, pattern = "[.](csv|tsv|txt|xls|xlsx)$", recursive = TRUE, full.names = TRUE)
  ))

  ## Avoid derived score outputs but allow raw expression/annotation files.
  candidate_files <- candidate_files[!grepl("Figure3A|score|coverage|paired_stats|projected|CellAge_human_to_mouse", basename(candidate_files), ignore.case = TRUE)]

  recovered_maps <- list()
  for (ff in candidate_files) {
    tmp <- build_symbol_map_from_candidate(ff, expr_gene_ids)
    if (nrow(tmp) > 0) {
      recovered_maps[[ff]] <- tmp
      annotation_sources_used <- rbind(
        annotation_sources_used,
        data.frame(source = ff, n_mapped_rows = nrow(tmp), stringsAsFactors = FALSE)
      )
    }
  }

  if (length(recovered_maps) > 0) {
    expr_symbol_map <- unique(do.call(rbind, recovered_maps))
  }
}

expr_symbol_map <- expr_symbol_map[
  !is.na(expr_symbol_map$mouse_symbol) &
    expr_symbol_map$mouse_symbol != "" &
    expr_symbol_map$gene_id %in% rownames(expr_all),
  ,
  drop = FALSE
]

write_csv(annotation_sources_used, file.path(table_dir, "step12A_gene_symbol_annotation_sources_used.csv"))
write_csv(expr_symbol_map, file.path(table_dir, "step12A_genomewide_gene_id_to_mouse_symbol_map.csv"))

if (nrow(expr_symbol_map) == 0) {
  stop(
    "Cannot build genome-wide mouse_symbol -> gene_id map. ",
    "No Step04 symbol column was found, and no project/raw annotation file with st_gene_id plus SYMBOL/Entrez mapping could be recovered. ",
    "Please provide the original GSE271903 count matrix or a gene annotation table containing st_gene_id and Entrez/SYMBOL."
  )
}

symbol_counts <- aggregate(gene_id ~ mouse_symbol, data = expr_symbol_map, FUN = function(x) length(unique(x)))
ambiguous_symbols <- symbol_counts$mouse_symbol[symbol_counts$gene_id > 1]

if (length(ambiguous_symbols) > 0) {
  ambiguous_df <- expr_symbol_map[expr_symbol_map$mouse_symbol %in% ambiguous_symbols, , drop = FALSE]
  write_csv(ambiguous_df, file.path(table_dir, "step12A_ambiguous_mouse_symbols_multiple_gene_ids_excluded.csv"))
  expr_symbol_map <- expr_symbol_map[!(expr_symbol_map$mouse_symbol %in% ambiguous_symbols), , drop = FALSE]
}

map_gene <- merge(map_std, expr_symbol_map, by = "mouse_symbol", all.x = TRUE)
unmapped <- map_gene[is.na(map_gene$gene_id) | map_gene$gene_id == "", , drop = FALSE]
if (nrow(unmapped) > 0) {
  write_csv(unmapped, file.path(table_dir, "step12A_projected_CellAge_mouse_symbols_not_detected_or_unmapped.csv"))
}

map_gene_detected <- map_gene[!is.na(map_gene$gene_id) & map_gene$gene_id != "", , drop = FALSE]
map_gene_detected <- map_gene_detected[!duplicated(paste(map_gene_detected$human_symbol, map_gene_detected$gene_id, map_gene_detected$CellAge_effect, sep = "|||")), , drop = FALSE]

all_genes <- unique(map_gene_detected$gene_id)
induces_genes <- unique(map_gene_detected$gene_id[map_gene_detected$score_type == "CellAge_induces"])
inhibits_genes <- unique(map_gene_detected$gene_id[map_gene_detected$score_type == "CellAge_inhibits"])
all_only_genes <- unique(map_gene_detected$gene_id[map_gene_detected$score_type == "CellAge_all_only"])

if (length(all_genes) == 0) stop("No detected whole-CellAge projected mouse genes found in expression matrix.")
if (length(induces_genes) == 0) stop("No detected CellAge_induces genes found in expression matrix.")
if (length(inhibits_genes) == 0) stop("No detected CellAge_inhibits genes found in expression matrix.")

logcpm_cellage <- logcpm_all[all_genes, , drop = FALSE]
z_cellage <- row_zscore(logcpm_cellage)

## =========================
## 7. Compute sample-level scores
## =========================

score_module <- function(z_mat, genes) {
  present <- intersect(genes, rownames(z_mat))
  if (length(present) == 0) return(rep(NA_real_, ncol(z_mat)))
  colMeans(z_mat[present, , drop = FALSE], na.rm = TRUE)
}

scores_wide <- data.frame(
  title = meta_all$title,
  mouse_id = meta_all$mouse_id,
  treatment = meta_all$treatment,
  time_point = meta_all$time_point,
  Sex = meta_all$sex,
  CellAge_all = score_module(z_cellage, all_genes),
  CellAge_induces = score_module(z_cellage, induces_genes),
  CellAge_inhibits = score_module(z_cellage, inhibits_genes),
  stringsAsFactors = FALSE
)

score_long <- reshape(
  scores_wide,
  varying = c("CellAge_all", "CellAge_induces", "CellAge_inhibits"),
  v.names = "score",
  timevar = "score_type",
  times = c("CellAge_all", "CellAge_induces", "CellAge_inhibits"),
  direction = "long"
)
rownames(score_long) <- NULL

score_long$score_label <- c(
  CellAge_all = "CellAge all",
  CellAge_induces = "CellAge induces",
  CellAge_inhibits = "CellAge inhibits"
)[score_long$score_type]

score_long$time_point <- factor(score_long$time_point, levels = c("1W", "4W"))
score_long$treatment <- factor(score_long$treatment, levels = c("Contra", "ACLR"))
score_long$score_label <- factor(score_long$score_label, levels = c("CellAge all", "CellAge induces", "CellAge inhibits"))

## =========================
## 8. Paired statistics
## =========================

stats_list <- list()
for (score_label_i in levels(score_long$score_label)) {
  for (time_i in levels(score_long$time_point)) {
    sub <- score_long[score_long$score_label == score_label_i & score_long$time_point == time_i, , drop = FALSE]
    st <- paired_wilcox(sub)
    st$score_label <- score_label_i
    st$time_point <- time_i
    st$score_type <- unique(as.character(sub$score_type))
    st$y_max <- max(sub$score, na.rm = TRUE)
    st$y_label <- st$y_max + 0.08
    stats_list[[paste(score_label_i, time_i, sep = "_")]] <- st
  }
}
stats <- do.call(rbind, stats_list)
rownames(stats) <- NULL
stats <- stats[, c(
  "score_label", "time_point", "score_type", "n_pairs",
  "mean_contra", "mean_aclr", "mean_delta_aclr_minus_contra",
  "median_delta_aclr_minus_contra", "wilcox_p", "p_label", "y_max", "y_label"
)]

## =========================
## 9. Save source data
## =========================

coverage <- data.frame(
  gene_set = c("CellAge_all", "CellAge_induces", "CellAge_inhibits"),
  n_mouse_projected_genes = c(
    length(unique(map_std$mouse_symbol)),
    length(unique(map_std$mouse_symbol[map_std$score_type == "CellAge_induces"])),
    length(unique(map_std$mouse_symbol[map_std$score_type == "CellAge_inhibits"]))
  ),
  n_genes_present_in_expr = c(length(all_genes), length(induces_genes), length(inhibits_genes)),
  stringsAsFactors = FALSE
)

detected <- data.frame(
  score_type = c("CellAge_all", "CellAge_induces", "CellAge_inhibits"),
  genes_in_set = coverage$n_mouse_projected_genes,
  genes_detected_in_expr_all = coverage$n_genes_present_in_expr,
  expression_scale = "edgeR_TMM_logCPM_48_samples",
  row_z_scope = "all_48_mouse_discovery_samples",
  stringsAsFactors = FALSE
)

write_csv(map_gene_detected, file.path(table_dir, "step12A_whole_CellAge_projected_mouse_mapping_detected_from_scratch.csv"))
write_csv(coverage, file.path(table_dir, "step12A_CellAge_score_gene_coverage.csv"))
write_csv(detected, file.path(table_dir, "step12A_detected_gene_summary.csv"))
write_csv(scores_wide, file.path(table_dir, "step12A_Figure3A_CellAge_scores_by_sample.csv"))
write_csv(score_long, file.path(table_dir, "step12A_Figure3A_CellAge_scores_long.csv"))
write_csv(stats, file.path(table_dir, "step12A_Figure3A_CellAge_paired_stats.csv"))

logcpm_export <- as.data.frame(logcpm_cellage)
logcpm_export$gene_id <- rownames(logcpm_export)
write_csv(logcpm_export, file.path(table_dir, "step12A_TMM_logCPM_whole_CellAge_projected_genes_48samples.csv"))

z_export <- as.data.frame(z_cellage)
z_export$gene_id <- rownames(z_export)
write_csv(z_export, file.path(table_dir, "step12A_row_z_whole_CellAge_projected_genes_48samples.csv"))

## =========================
## 10. Plot Figure3A
## =========================

p <- ggplot() +
  geom_line(
    data = score_long,
    aes(x = treatment, y = score, group = interaction(score_label, time_point, mouse_id)),
    color = "#B8B8B8",
    linewidth = 0.65
  ) +
  geom_point(
    data = score_long[score_long$treatment == "Contra", , drop = FALSE],
    aes(x = treatment, y = score),
    shape = 21,
    fill = "white",
    color = "#7A7A7A",
    stroke = 0.5,
    size = 3.1
  ) +
  geom_point(
    data = score_long[score_long$treatment == "ACLR", , drop = FALSE],
    aes(x = treatment, y = score),
    shape = 21,
    fill = "#D95F02",
    color = "#7F3300",
    stroke = 0.5,
    size = 3.1
  ) +
  geom_text(
    data = stats,
    aes(x = 1.5, y = y_label, label = p_label),
    size = 4.7
  ) +
  facet_grid(rows = vars(score_label), cols = vars(time_point), switch = "y") +
  labs(
    title = "Projected CellAge scores in paired synovium samples",
    x = NULL,
    y = "Mean row-z CellAge score"
  ) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    panel.grid.major = element_line(color = "#E6E6E6"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0", color = "black"),
    strip.text = element_text(size = 13),
    strip.text.y.left = element_text(angle = -90, size = 13),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )

ggsave(png_file, p, width = 12.5, height = 10.5, dpi = 320)
ggsave(pdf_file, p, width = 12.5, height = 10.5)

## =========================
## 11. Save versions and summary
## =========================

summary_df <- data.frame(
  metric = c(
    "score_scope",
    "CellAge_raw_file",
    "ortholog_mapping",
    "expression_input",
    "normalization",
    "log_expression",
    "row_z_scope",
    "CellAge_all_projected_mouse_symbols",
    "CellAge_all_detected_gene_ids",
    "CellAge_induces_detected_gene_ids",
    "CellAge_inhibits_detected_gene_ids",
    "expr_all_genes_after_NA_row_removal",
    "expr_all_samples",
    "genes_removed_before_DGEList_due_to_NA_counts",
    "ambiguous_mouse_symbols_excluded",
    "score_rows_long"
  ),
  value = c(
    "whole CellAge projected mouse gene set, not the persistent intersect CellAge gene subset",
    cellage_file,
    "gprofiler2::gorth hsapiens to mmusculus strict 1:1 from scratch",
    "Step04 cleaned count matrices, 1W + 4W combined",
    "edgeR::calcNormFactors default TMM",
    "edgeR::cpm(log = TRUE, prior.count = 1)",
    "all 48 mouse discovery samples",
    coverage$n_mouse_projected_genes[coverage$gene_set == "CellAge_all"],
    coverage$n_genes_present_in_expr[coverage$gene_set == "CellAge_all"],
    coverage$n_genes_present_in_expr[coverage$gene_set == "CellAge_induces"],
    coverage$n_genes_present_in_expr[coverage$gene_set == "CellAge_inhibits"],
    nrow(expr_all),
    ncol(expr_all),
    length(na_gene_ids),
    length(ambiguous_symbols),
    nrow(score_long)
  ),
  stringsAsFactors = FALSE
)

write_csv(summary_df, file.path(table_dir, "step12A_summary.csv"))

versions <- data.frame(
  item = c("R", "edgeR", "ggplot2", "gprofiler2"),
  version = c(
    as.character(getRversion()),
    as.character(utils::packageVersion("edgeR")),
    as.character(utils::packageVersion("ggplot2")),
    as.character(utils::packageVersion("gprofiler2"))
  ),
  stringsAsFactors = FALSE
)
write_csv(versions, file.path(table_dir, "step12A_software_versions.csv"))

version_audit <- rbind(
  cellage_file_audit,
  data.frame(
    item = c(
      "gprofiler2_package_version",
      "gprofiler2_version_info",
      "gorth_source_organism",
      "gorth_target_organism",
      "gorth_mthreshold",
      "gorth_filter_na",
      "edgeR_package_version",
      "clusterProfiler_package_version",
      "org.Mm.eg.db_package_version"
    ),
    value = c(
      gprofiler_package_version,
      gprofiler_version_info,
      "hsapiens",
      "mmusculus",
      "Inf",
      "FALSE",
      as.character(utils::packageVersion("edgeR")),
      as.character(utils::packageVersion("clusterProfiler")),
      as.character(utils::packageVersion("org.Mm.eg.db"))
    ),
    stringsAsFactors = FALSE
  )
)
write_csv(version_audit, file.path(table_dir, "step12A_CellAge_and_mapping_version_audit.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "step12A_sessionInfo.txt"))

cat("\n===== STEP12A WHOLE CELLAGE FROM SCRATCH SUMMARY =====\n")
cat("\nFigure output files:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")

cat("\nSummary:\n")
print(summary_df)

cat("\nCoverage:\n")
print(coverage)

cat("\nPaired statistics:\n")
print(stats)

cat("\nScore table dimensions:\n")
print(dim(score_long))

cat("\nCellAge and mapping version audit:\n")
print(version_audit)

cat("\nStep12A whole CellAge from scratch completed successfully.\n")

sink()

cat("\nStep12A whole CellAge from scratch completed. Please open:\n")
cat(png_file, "\n")
cat(pdf_file, "\n")
