options(stringsAsFactors = FALSE)

out_dir <- "D:/workspace/ACLsenescence2_reviewer3_persistent_ORA"
all_file <- file.path(out_dir, "persistent_ORA_all_Reactome_KEGG_results.csv")
sig_file <- file.path(out_dir, "persistent_ORA_primary_FDR_significant_results.csv")
audit_file <- file.path(out_dir, "persistent_ORA_database_and_universe_audit.csv")
input_file <- file.path(out_dir, "persistent_ORA_input_genes_1416.csv")
membership_file <- file.path(out_dir, "persistent_ORA_pathway_membership_snapshot.csv")

required <- c(all_file, sig_file, audit_file, input_file, membership_file)
if (!all(file.exists(required))) {
  stop("Missing required ORA output: ", paste(required[!file.exists(required)], collapse = "; "), call. = FALSE)
}

ora <- utils::read.csv(all_file, stringsAsFactors = FALSE, check.names = FALSE)
sig <- utils::read.csv(sig_file, stringsAsFactors = FALSE, check.names = FALSE)
audit <- utils::read.csv(audit_file, stringsAsFactors = FALSE, check.names = FALSE)
input <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
membership <- utils::read.csv(membership_file, stringsAsFactors = FALSE, check.names = FALSE)

checks <- list()
add_check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    status = if (isTRUE(pass)) "PASS" else "FAIL",
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

add_check("input_rows", nrow(input) == 1416L, nrow(input))
add_check("input_unique_entrez", length(unique(input$ENTREZID)) == 1416L,
          length(unique(input$ENTREZID)))
add_check("input_direction_counts",
          sum(input$direction_persistent == "Up") == 970L &&
            sum(input$direction_persistent == "Down") == 446L,
          paste0("up=", sum(input$direction_persistent == "Up"),
                 "; down=", sum(input$direction_persistent == "Down")))
add_check("input_strict_thresholds",
          all(input$FDR_1W < 0.05 & input$FDR_4W < 0.05 &
                abs(input$logFC_1W) > 1 & abs(input$logFC_4W) > 1),
          "FDR<0.05 and |logFC|>1 at both timepoints")

expected_rows <- 3L * sum(audit$gene_sets_tested_size_10_to_500)
add_check("all_result_rows", nrow(ora) == expected_rows,
          paste0("observed=", nrow(ora), "; expected=", expected_rows))
add_check("result_unique_key",
          !anyDuplicated(ora[, c("database", "gene_list", "pathway")]),
          paste0("duplicate_keys=", sum(duplicated(ora[, c("database", "gene_list", "pathway")]))))
add_check("membership_unique_key",
          !anyDuplicated(membership[, c("database", "pathway", "gene_id")]),
          paste0("duplicate_rows=", sum(duplicated(membership[, c("database", "pathway", "gene_id")]))))

p_recalc <- stats::phyper(
  ora$overlap_size - 1L,
  ora$pathway_size_in_universe,
  ora$universe_size - ora$pathway_size_in_universe,
  ora$query_size_in_universe,
  lower.tail = FALSE
)
max_p_diff <- max(abs(p_recalc - ora$pval), na.rm = TRUE)
add_check("hypergeometric_p_recalculation", max_p_diff < 1e-12,
          format(max_p_diff, scientific = TRUE))

within_recalc <- rep(NA_real_, nrow(ora))
for (key in unique(paste(ora$database, ora$gene_list, sep = "||"))) {
  idx <- which(paste(ora$database, ora$gene_list, sep = "||") == key)
  within_recalc[idx] <- stats::p.adjust(ora$pval[idx], method = "BH")
}
max_within_diff <- max(abs(within_recalc - ora$FDR_within_database), na.rm = TRUE)
add_check("within_database_BH_recalculation", max_within_diff < 1e-12,
          format(max_within_diff, scientific = TRUE))

joint_recalc <- rep(NA_real_, nrow(ora))
for (gene_list in unique(ora$gene_list)) {
  idx <- which(ora$gene_list == gene_list)
  joint_recalc[idx] <- stats::p.adjust(ora$pval[idx], method = "BH")
}
max_joint_diff <- max(abs(joint_recalc - ora$FDR_joint_Reactome_KEGG), na.rm = TRUE)
add_check("joint_Reactome_KEGG_BH_recalculation", max_joint_diff < 1e-12,
          format(max_joint_diff, scientific = TRUE))

count_ids <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x) | x == "", 0L, lengths(strsplit(x, ";", fixed = TRUE)))
}
overlap_count_diff <- max(abs(count_ids(ora$overlap_gene_ids) - ora$overlap_size))
add_check("overlap_gene_count", overlap_count_diff == 0L, overlap_count_diff)
add_check("contingency_bounds",
          all(ora$overlap_size <= ora$query_size_in_universe &
                ora$overlap_size <= ora$pathway_size_in_universe &
                ora$query_size_in_universe <= ora$universe_size &
                ora$pathway_size_in_universe <= ora$universe_size),
          "all k<=n,K and n,K<=N")
add_check("probability_ranges",
          all(is.finite(ora$pval) & ora$pval >= 0 & ora$pval <= 1 &
                is.finite(ora$FDR_within_database) & ora$FDR_within_database >= 0 &
                ora$FDR_within_database <= 1 &
                is.finite(ora$FDR_joint_Reactome_KEGG) & ora$FDR_joint_Reactome_KEGG >= 0 &
                ora$FDR_joint_Reactome_KEGG <= 1),
          "all P and FDR values in [0,1]")
add_check("primary_significant_count", nrow(sig) == 101L,
          paste0("observed=", nrow(sig), "; up=", sum(sig$gene_list == "persistent_up"),
                 "; down=", sum(sig$gene_list == "persistent_down")))

sig_key <- paste(sig$database, sig$gene_list, sig$pathway, sep = "||")
expected_sig <- ora[ora$primary_test & ora$primary_FDR < 0.05, , drop = FALSE]
expected_sig_key <- paste(expected_sig$database, expected_sig$gene_list, expected_sig$pathway, sep = "||")
add_check("significant_file_exact_subset",
          setequal(sig_key, expected_sig_key) && nrow(sig) == nrow(expected_sig),
          paste0("saved=", nrow(sig), "; recomputed=", nrow(expected_sig)))
add_check("database_releases_recorded",
          any(grepl("2026.1.Mm", audit$database_release, fixed = TRUE)) &&
            any(grepl("2026-08-31", audit$database_release, fixed = TRUE)),
          paste(audit$database, audit$database_release, collapse = "; "))
add_check("pathway_names_complete", all(!is.na(ora$pathway_name) & nzchar(ora$pathway_name)),
          paste0("missing=", sum(is.na(ora$pathway_name) | !nzchar(ora$pathway_name))))

qa <- do.call(rbind, checks)
utils::write.csv(qa, file.path(out_dir, "persistent_ORA_QA_checks.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c(
  "Persistent Reactome/KEGG ORA QA report",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Checks passed: ", sum(qa$status == "PASS"), "/", nrow(qa)),
  paste0("Overall status: ", if (all(qa$status == "PASS")) "PASS" else "FAIL"),
  "",
  capture.output(print(qa, row.names = FALSE))
), file.path(out_dir, "persistent_ORA_QA_report.txt"))

print(qa, row.names = FALSE)
if (!all(qa$status == "PASS")) stop("ORA QA failed.", call. = FALSE)
cat("ALL_QA_CHECKS_PASS\n")
