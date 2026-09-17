## Compatibility audit of deposited GSE228848 Salmon transcript IDs
## against the Ensembl release 115 tx2gene table used for tximport.

options(stringsAsFactors = FALSE)

quant_dir <- "E:/R/ACLsenescence2/rebuild_submission/raw data/GSE228848_synovium_quant"
tx2gene_file <- "E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_tx2gene.csv"
out_dir <- "D:/workspace/ACLsenescence2_chronic_salmon_methods_revision"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

quant_files <- sort(list.files(
  quant_dir,
  pattern = "quant[.]sf.*[.]gz$",
  full.names = TRUE,
  ignore.case = TRUE
))
stopifnot(length(quant_files) == 96L)

## This mirrors the transcript-version handling used by tximport with
## ignoreTxVersion = TRUE for Ensembl identifiers such as ENSSSCT... .4.
strip_version <- function(x) sub("[.][0-9]+$", "", x)

read_quant <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  tab <- read.delim(
    con,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    colClasses = c("character", "NULL", "NULL", "NULL", "numeric"),
    check.names = FALSE
  )
  list(ids = tab[[1L]], num_reads = tab[[2L]])
}

quant_by_file <- lapply(quant_files, read_quant)
ids_by_file <- lapply(quant_by_file, `[[`, "ids")
num_reads_by_file <- lapply(quant_by_file, `[[`, "num_reads")
raw_unique_by_file <- lapply(ids_by_file, unique)
normalized_unique_by_file <- lapply(raw_unique_by_file, strip_version)

deposited_unique_ids <- sort(unique(unlist(normalized_unique_by_file, use.names = FALSE)))
tx2gene <- read.csv(tx2gene_file, check.names = FALSE)
tx2gene_ids <- sort(unique(strip_version(tx2gene$transcript_id)))

mapped <- deposited_unique_ids %in% tx2gene_ids
per_file <- data.frame(
  file = basename(quant_files),
  n_rows = vapply(ids_by_file, length, integer(1)),
  n_unique_ids_after_version_normalization = vapply(normalized_unique_by_file, length, integer(1)),
  n_mapped_to_release115_tx2gene = vapply(
    normalized_unique_by_file,
    function(x) sum(x %in% tx2gene_ids),
    integer(1)
  ),
  mapping_percentage = vapply(
    normalized_unique_by_file,
    function(x) 100 * mean(x %in% tx2gene_ids),
    numeric(1)
  ),
  total_NumReads = vapply(num_reads_by_file, sum, numeric(1), na.rm = TRUE),
  mapped_NumReads = vapply(
    seq_along(num_reads_by_file),
    function(i) {
      m <- strip_version(ids_by_file[[i]]) %in% tx2gene_ids
      sum(num_reads_by_file[[i]][m], na.rm = TRUE)
    },
    numeric(1)
  )
)
per_file$unmapped_NumReads <- per_file$total_NumReads - per_file$mapped_NumReads
per_file$NumReads_mapping_percentage <- 100 * per_file$mapped_NumReads / per_file$total_NumReads

summary <- data.frame(
  metric = c(
    "n_quant_files",
    "unique_deposited_transcript_ids_after_version_normalization",
    "unique_ids_mapped_to_ensembl115_tx2gene",
    "unique_ids_unmapped",
    "mapping_percentage",
    "unmapped_percentage",
    "NumReads_mapping_percentage_mean",
    "NumReads_mapping_percentage_min",
    "NumReads_mapping_percentage_max",
    "NumReads_mapping_percentage_weighted_all_files",
    "per_file_mapping_percentage_min",
    "per_file_mapping_percentage_max"
  ),
  value = c(
    length(quant_files),
    length(deposited_unique_ids),
    sum(mapped),
    sum(!mapped),
    100 * mean(mapped),
    100 * mean(!mapped),
    mean(per_file$NumReads_mapping_percentage),
    min(per_file$NumReads_mapping_percentage),
    max(per_file$NumReads_mapping_percentage),
    100 * sum(per_file$mapped_NumReads) / sum(per_file$total_NumReads),
    min(per_file$mapping_percentage),
    max(per_file$mapping_percentage)
  )
)

write.csv(
  summary,
  file.path(out_dir, "transcript_mapping_compatibility_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  per_file,
  file.path(out_dir, "transcript_mapping_compatibility_per_file.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  data.frame(
    transcript_id = deposited_unique_ids,
    mapped_to_ensembl115_tx2gene = mapped
  ),
  file.path(out_dir, "transcript_mapping_compatibility_ids.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  c(
    "GSE228848 transcript-ID compatibility audit",
    paste0("Quantification files: ", length(quant_files)),
    paste0("Unique deposited transcript IDs after version normalization: ", length(deposited_unique_ids)),
    paste0("Mapped to Ensembl release 115 tx2gene: ", sum(mapped)),
    paste0("Unmapped: ", sum(!mapped)),
    sprintf("Compatibility: %.5f%%", 100 * mean(mapped)),
    "Version suffixes were removed before matching, consistent with tximport ignoreTxVersion = TRUE."
  ),
  file.path(out_dir, "transcript_mapping_compatibility_audit.txt"),
  useBytes = TRUE
)
