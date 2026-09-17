options(stringsAsFactors = FALSE)

## Reviewer 1 negative-control analysis for the mouse-defined directional score.
## Outputs are isolated under D:/workspace and do not overwrite frozen project results.

n_perm <- 10000L
random_seed_by_dataset <- c(pig_early = 20260901L, pig_chronic = 20260902L)
root <- "E:/R/ACLsenescence2/rebuild_submission"
out_dir <- "D:/workspace/ACLsenescence2_reviewer1_directional_score_null"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
orthology_mapping_path <- file.path(out_dir, "orthology_eligibility_mapping.csv")
orthology_raw_path <- file.path(out_dir, "orthology_gorth_raw.csv")

if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is required.", call. = FALSE)
}
if (!requireNamespace("gprofiler2", quietly = TRUE)) {
  stop("gprofiler2 is required to construct the orthologue-eligible control universe.", call. = FALSE)
}

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

clean_id <- function(x) sub("\\.[0-9]+$", "", trimws(as.character(x)))

row_z <- function(mat) {
  mu <- rowMeans(mat)
  s <- apply(mat, 1L, stats::sd)
  z <- sweep(mat, 1L, mu, "-")
  z <- sweep(z, 1L, s, "/")
  z[!is.finite(z)] <- NA_real_
  z
}

make_bins <- function(x, n_bins) {
  br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n_bins + 1L),
                               na.rm = TRUE, names = FALSE, type = 8))
  if (length(br) < 2L) return(rep("1", length(x)))
  as.character(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE))
}

prepare_data <- function(count_file, manifest_file, group_map, expected_groups) {
  counts <- read_csv(count_file)
  manifest <- read_csv(manifest_file)
  gene_candidates <- c("gene_id", "Geneid", "gene", "GeneID", "ensembl_gene_id", "pig_gene_id")
  gene_col <- gene_candidates[gene_candidates %in% names(counts)][1]
  if (is.na(gene_col)) stop("No gene ID column in ", count_file, call. = FALSE)

  sample_candidates <- c("sample_id", "sample", "SampleID", "geo_accession", "GSM", "run", "Run")
  sample_col <- sample_candidates[sample_candidates %in% names(manifest)][1]
  group_candidates <- c("core_group", "group", "treatment", "condition", "title", "sample_title")
  group_col <- group_candidates[group_candidates %in% names(manifest)][1]
  if (is.na(sample_col) || is.na(group_col)) stop("Manifest columns unresolved: ", manifest_file, call. = FALSE)

  sample_id <- trimws(as.character(manifest[[sample_col]]))
  group_raw <- trimws(as.character(manifest[[group_col]]))
  group <- unname(group_map[group_raw])
  keep <- !is.na(group) & sample_id %in% names(counts)
  sample_id <- sample_id[keep]
  group <- group[keep]
  ord <- order(match(group, expected_groups), sample_id)
  sample_id <- sample_id[ord]
  group <- group[ord]
  if (!all(expected_groups %in% group)) stop("Expected groups missing from ", manifest_file, call. = FALSE)

  ids <- clean_id(counts[[gene_col]])
  mat <- as.matrix(data.frame(lapply(counts[, sample_id, drop = FALSE], function(v) as.numeric(as.character(v))),
                              check.names = FALSE))
  rownames(mat) <- ids
  mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]
  if (anyDuplicated(rownames(mat))) mat <- rowsum(mat, rownames(mat), reorder = FALSE)
  if (anyNA(mat) || any(mat < 0)) stop("Invalid counts in ", count_file, call. = FALSE)

  dge <- edgeR::DGEList(counts = mat, group = group)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
  list(logcpm = logcpm, group = group, sample_id = sample_id)
}

build_orthology_eligibility <- function(gene_ids, mapping_path, raw_path, batch_size = 1000L) {
  gene_ids <- sort(unique(clean_id(gene_ids)))
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]
  if (file.exists(mapping_path)) {
    cached <- read_csv(mapping_path)
    required <- c("input_ensg", "ortholog_ensg", "n_human_orthologs", "one_to_one_eligible")
    if (all(required %in% names(cached)) && all(gene_ids %in% cached$input_ensg)) {
      return(cached[match(gene_ids, cached$input_ensg), , drop = FALSE])
    }
  }

  batches <- split(gene_ids, ceiling(seq_along(gene_ids) / batch_size))
  raw_list <- vector("list", length(batches))
  for (bi in seq_along(batches)) {
    query <- batches[[bi]]
    message(sprintf("g:Profiler orthology batch %d/%d (%d pig genes)", bi, length(batches), length(query)))
    res <- NULL
    last_error <- ""
    for (attempt in seq_len(3L)) {
      res <- tryCatch(
        gprofiler2::gorth(query = query, source_organism = "sscrofa",
                          target_organism = "hsapiens", mthreshold = Inf,
                          filter_na = TRUE),
        error = function(e) {
          last_error <<- conditionMessage(e)
          NULL
        }
      )
      if (!is.null(res)) break
    }
    if (is.null(res)) {
      stop("g:Profiler orthology query failed for batch ", bi, ": ", last_error, call. = FALSE)
    }
    if (nrow(res) > 0L) raw_list[[bi]] <- res
  }

  raw <- do.call(rbind, raw_list[!vapply(raw_list, is.null, logical(1))])
  if (is.null(raw) || nrow(raw) == 0L) {
    stop("g:Profiler returned no pig-to-human orthology mappings.", call. = FALSE)
  }
  if (!all(c("input_ensg", "ortholog_ensg") %in% names(raw))) {
    stop("Unexpected g:Profiler gorth columns; input_ensg/ortholog_ensg are required.", call. = FALSE)
  }
  raw$input_ensg <- clean_id(raw$input_ensg)
  raw$ortholog_ensg <- clean_id(raw$ortholog_ensg)
  raw <- raw[!is.na(raw$input_ensg) & nzchar(raw$input_ensg) &
               !is.na(raw$ortholog_ensg) & nzchar(raw$ortholog_ensg), , drop = FALSE]
  raw <- raw[!duplicated(raw[c("input_ensg", "ortholog_ensg")]), , drop = FALSE]
  utils::write.csv(raw, raw_path, row.names = FALSE)

  targets <- split(raw$ortholog_ensg, raw$input_ensg)
  n_targets <- vapply(gene_ids, function(g) {
    if (!g %in% names(targets)) 0L else length(unique(targets[[g]]))
  }, integer(1))
  ortholog_ensg <- vapply(gene_ids, function(g) {
    if (n_targets[match(g, gene_ids)] != 1L) NA_character_ else unique(targets[[g]])[1L]
  }, character(1))
  mapping <- data.frame(
    input_ensg = gene_ids,
    ortholog_ensg = ortholog_ensg,
    n_human_orthologs = n_targets,
    one_to_one_eligible = n_targets == 1L,
    mapping_status = ifelse(n_targets == 0L, "unmapped",
                            ifelse(n_targets == 1L, "one_to_one", "one_to_many")),
    stringsAsFactors = FALSE
  )
  utils::write.csv(mapping, mapping_path, row.names = FALSE)
  mapping
}

matched_null <- function(dataset_name, dat, signature, comparisons, orthology, random_seed) {
  sig <- signature[signature$pig_ensg %in% rownames(dat$logcpm), , drop = FALSE]
  if (nrow(sig) != 75L || sum(sig$direction == "Up_in_ACLR") != 65L ||
      sum(sig$direction == "Down_in_ACLR") != 10L) {
    stop(dataset_name, ": signature is not 75/65/10 after matching.", call. = FALSE)
  }

  orthology <- orthology[match(rownames(dat$logcpm), orthology$input_ensg), , drop = FALSE]
  if (any(is.na(orthology$input_ensg))) stop(dataset_name, ": orthology map does not cover all measured genes.", call. = FALSE)
  eligible_ortholog_genes <- orthology$input_ensg[orthology$one_to_one_eligible]
  if (!all(sig$pig_ensg %in% eligible_ortholog_genes)) {
    stop(dataset_name, ": one or more signature genes are not one-to-one orthologue-eligible.", call. = FALSE)
  }

  gene_mean <- rowMeans(dat$logcpm)
  gene_sd <- apply(dat$logcpm, 1L, stats::sd)
  eligible <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
  universe <- data.frame(gene = rownames(dat$logcpm), mean = gene_mean, sd = gene_sd,
                         stringsAsFactors = FALSE)
  universe <- universe[eligible & universe$gene %in% eligible_ortholog_genes &
                         !(universe$gene %in% sig$pig_ensg), , drop = FALSE]
  sig_feat <- data.frame(gene = sig$pig_ensg, direction = sig$direction,
                         mean = gene_mean[sig$pig_ensg], sd = gene_sd[sig$pig_ensg],
                         stringsAsFactors = FALSE)

  chosen_bins <- NA_integer_
  for (nb in c(5L, 4L, 3L, 2L)) {
    all_mean <- c(universe$mean, sig_feat$mean)
    all_sd <- c(universe$sd, sig_feat$sd)
    mb <- make_bins(all_mean, nb)
    sb <- make_bins(all_sd, nb)
    universe$base_stratum <- paste(mb[seq_len(nrow(universe))], sb[seq_len(nrow(universe))], sep = "_")
    sig_idx <- nrow(universe) + seq_len(nrow(sig_feat))
    sig_feat$base_stratum <- paste(mb[sig_idx], sb[sig_idx], sep = "_")
    need <- table(sig_feat$base_stratum)
    have <- table(universe$base_stratum)
    feasible <- all(names(need) %in% names(have)) && all(have[names(need)] >= need)
    if (feasible) { chosen_bins <- nb; break }
  }
  if (is.na(chosen_bins)) stop(dataset_name, ": unable to construct matched strata.", call. = FALSE)

  z_all <- row_z(dat$logcpm)
  sign_obs <- ifelse(sig_feat$direction == "Up_in_ACLR", 1, -1)
  observed_score <- colMeans(z_all[sig_feat$gene, , drop = FALSE] * sign_obs)

  effect <- function(score, case_group, control_group) {
    stats::median(score[dat$group == case_group]) - stats::median(score[dat$group == control_group])
  }
  observed_effect <- vapply(seq_len(nrow(comparisons)), function(i) {
    effect(observed_score, comparisons$case[i], comparisons$control[i])
  }, numeric(1))

  set.seed(random_seed)
  null_effect <- matrix(NA_real_, nrow = n_perm, ncol = nrow(comparisons))
  colnames(null_effect) <- comparisons$comparison
  strata <- unique(sig_feat$base_stratum)
  for (b in seq_len(n_perm)) {
    genes <- character(0)
    signs <- numeric(0)
    for (st in strata) {
      ssub <- sig_feat[sig_feat$base_stratum == st, , drop = FALSE]
      cand <- universe$gene[universe$base_stratum == st]
      ## Draw the required total number jointly so all genes remain unique;
      ## then assign the prespecified up/down counts observed in this stratum.
      draw <- sample(cand, nrow(ssub), replace = FALSE)
      n_up <- sum(ssub$direction == "Up_in_ACLR")
      n_down <- sum(ssub$direction == "Down_in_ACLR")
      draw_sign <- sample(c(rep(1, n_up), rep(-1, n_down)), length(draw), replace = FALSE)
      genes <- c(genes, draw)
      signs <- c(signs, draw_sign)
    }
    null_score <- colMeans(z_all[genes, , drop = FALSE] * signs)
    null_effect[b, ] <- vapply(seq_len(nrow(comparisons)), function(i) {
      effect(null_score, comparisons$case[i], comparisons$control[i])
    }, numeric(1))
  }

  result <- data.frame(
    dataset = dataset_name,
    comparison = comparisons$comparison,
    n_control = vapply(comparisons$control, function(g) sum(dat$group == g), integer(1)),
    n_case = vapply(comparisons$case, function(g) sum(dat$group == g), integer(1)),
    observed_median_difference = observed_effect,
    null_mean = colMeans(null_effect),
    null_sd = apply(null_effect, 2L, stats::sd),
    null_q025 = apply(null_effect, 2L, stats::quantile, probs = 0.025),
    null_q50 = apply(null_effect, 2L, stats::quantile, probs = 0.5),
    null_q975 = apply(null_effect, 2L, stats::quantile, probs = 0.975),
    empirical_p_one_sided_activation = vapply(seq_along(observed_effect), function(i) {
      (1 + sum(null_effect[, i] >= observed_effect[i])) / (n_perm + 1)
    }, numeric(1)),
    empirical_p_two_sided = vapply(seq_along(observed_effect), function(i) {
      (1 + sum(abs(null_effect[, i]) >= abs(observed_effect[i]))) / (n_perm + 1)
    }, numeric(1)),
    observed_percentile_in_null = vapply(seq_along(observed_effect), function(i) {
      mean(null_effect[, i] < observed_effect[i])
    }, numeric(1)),
    n_permutations = n_perm,
    n_expression_variability_bins = chosen_bins,
    random_seed = random_seed,
    control_universe_definition = "Measured pig genes satisfying the same pig-side one-to-one pig-to-human orthologue eligibility criterion used for projected-signature construction, excluding the 75 signature genes",
    matching_definition = "Joint mean-expression x across-sample SD strata with prespecified up/down counts within each stratum",
    variability_definition = "Gene-level SD across all samples within the corresponding cohort",
    sampling_within_permutation = "Without replacement within each mean-expression x SD stratum",
    sampling_across_permutations = "Independent permutations; genes may recur across permutations",
    one_sided_test_definition = "Prespecified injury-concordant activation direction (Delta directional score > 0)",
    stringsAsFactors = FALSE
  )

  null_long <- do.call(rbind, lapply(seq_len(ncol(null_effect)), function(i) {
    data.frame(dataset = dataset_name, comparison = colnames(null_effect)[i],
               permutation = seq_len(n_perm), null_median_difference = null_effect[, i])
  }))
  diagnostics <- rbind(
    data.frame(dataset = dataset_name, set = "observed_signature", n = nrow(sig_feat),
               mean_expression = mean(sig_feat$mean), median_expression = stats::median(sig_feat$mean),
               mean_gene_sd = mean(sig_feat$sd), median_gene_sd = stats::median(sig_feat$sd)),
    data.frame(dataset = dataset_name, set = "eligible_control_universe", n = nrow(universe),
               mean_expression = mean(universe$mean), median_expression = stats::median(universe$mean),
               mean_gene_sd = mean(universe$sd), median_gene_sd = stats::median(universe$sd))
  )
  strata_levels <- unique(sig_feat$base_stratum)
  sig_tab <- table(factor(sig_feat$base_stratum, levels = strata_levels))
  up_tab <- table(factor(sig_feat$base_stratum[sig_feat$direction == "Up_in_ACLR"], levels = strata_levels))
  down_tab <- table(factor(sig_feat$base_stratum[sig_feat$direction == "Down_in_ACLR"], levels = strata_levels))
  universe_tab <- table(factor(universe$base_stratum, levels = strata_levels))
  list(result = result, null = null_long, diagnostics = diagnostics,
       strata = data.frame(dataset = dataset_name, base_stratum = strata_levels,
                           signature_n = as.integer(sig_tab), signature_up_n = as.integer(up_tab),
                           signature_down_n = as.integer(down_tab), universe_n = as.integer(universe_tab),
                           expression_variability_bins = chosen_bins,
                           stringsAsFactors = FALSE))
}

signature_path <- file.path(root, "02_pig_early", "tables", "step18_current78_pig_early_signature_remap",
                            "step18_current78_pig_signature_gene_table.csv")
sig_raw <- read_csv(signature_path)
signature <- data.frame(pig_ensg = clean_id(sig_raw$pig_ensg),
                        direction = trimws(sig_raw$signature_direction), stringsAsFactors = FALSE)
signature <- signature[!duplicated(signature$pig_ensg) &
                         signature$direction %in% c("Up_in_ACLR", "Down_in_ACLR"), , drop = FALSE]

early <- prepare_data(
  file.path(root, "02_pig_early", "tables", "step16_pig_early_gene_count_matrix.csv"),
  file.path(root, "02_pig_early", "tables", "step09_pig_early_core_sample_fastq_manifest.csv"),
  c(CON_t0 = "Control", ACLT_untreated_t7 = "ACLT_t7", ACLT_untreated_t28 = "ACLT_t28"),
  c("Control", "ACLT_t7", "ACLT_t28")
)
early_cmp <- data.frame(comparison = c("ACLT_t7_vs_Control", "ACLT_t28_vs_Control"),
                        control = "Control", case = c("ACLT_t7", "ACLT_t28"))

chronic <- prepare_data(
  file.path(root, "03_pig_chronic", "tables", "step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv"),
  file.path(root, "03_pig_chronic", "tables", "step25v3_pig_chronic_main_comparison_manifest.csv"),
  c(Control_52W = "Control_52W", ACLT_alone_52W = "ACLT_alone_52W"),
  c("Control_52W", "ACLT_alone_52W")
)
chronic_cmp <- data.frame(comparison = "ACLT_alone_52W_vs_Control_52W",
                          control = "Control_52W", case = "ACLT_alone_52W")

all_measured_genes <- unique(c(rownames(early$logcpm), rownames(chronic$logcpm)))
orthology <- build_orthology_eligibility(all_measured_genes, orthology_mapping_path, orthology_raw_path)

early_out <- matched_null("pig_early", early, signature, early_cmp, orthology,
                          random_seed_by_dataset[["pig_early"]])
chronic_out <- matched_null("pig_chronic", chronic, signature, chronic_cmp, orthology,
                            random_seed_by_dataset[["pig_chronic"]])
summary <- rbind(early_out$result, chronic_out$result)
summary$BH_FDR_one_sided_across_three_comparisons <- stats::p.adjust(summary$empirical_p_one_sided_activation, method = "BH")
summary$BH_FDR_two_sided_across_three_comparisons <- stats::p.adjust(summary$empirical_p_two_sided, method = "BH")

utils::write.csv(summary, file.path(out_dir, "directional_score_matched_null_summary.csv"), row.names = FALSE)
utils::write.csv(rbind(early_out$null, chronic_out$null),
                 file.path(out_dir, "directional_score_matched_null_all_permutations.csv"), row.names = FALSE)
utils::write.csv(rbind(early_out$diagnostics, chronic_out$diagnostics),
                 file.path(out_dir, "directional_score_matched_null_diagnostics.csv"), row.names = FALSE)
utils::write.csv(rbind(early_out$strata, chronic_out$strata),
                 file.path(out_dir, "directional_score_matched_null_strata_audit.csv"), row.names = FALSE)
session_lines <- c(
  R = R.version.string,
  paste0("edgeR=", as.character(utils::packageVersion("edgeR"))),
  paste0("gprofiler2=", as.character(utils::packageVersion("gprofiler2"))),
  paste0("n_perm=", n_perm),
  paste0("random_seed_pig_early=", random_seed_by_dataset[["pig_early"]]),
  paste0("random_seed_pig_chronic=", random_seed_by_dataset[["pig_chronic"]]),
  paste0("orthology_source=sscrofa_to_hsapiens_gProfiler_gorth"),
  paste0("gprofiler_query_date=2026-09-03"),
  paste0("orthology_mapping_path=", orthology_mapping_path),
  paste0("sd_scope=all_samples_within_corresponding_cohort"),
  paste0("sampling_within_permutation=without_replacement"),
  paste0("sampling_across_permutations=independent_permutations")
)
writeLines(session_lines, file.path(out_dir, "directional_score_matched_null_sessionInfo.txt"))
print(summary)
