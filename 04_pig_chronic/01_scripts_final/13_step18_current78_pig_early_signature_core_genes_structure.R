===== CHRONIC STEP2 CURRENT78 TMM-ALIGNED SIGNATURE SCORE =====
Run time: 2026-05-16 10:48:37 CEST
Method version: 2026-05-16_chronic_step2_current78_TMM_aligned_with_pig_early_Step19
Interpretation: chronic pig extension validation using fixed pig early current78 75-gene signature.
Normalization: edgeR TMM-normalized logCPM; no filterByExpr before scoring.
Primary score: directional_score = (sum z_up + sum -z_down)/(n_up+n_down).

警告信息:
程序包‘ggplot2’是用R版本4.5.3 来建造的 
Required packages loaded.
edgeR version: 4.8.2 
ggplot2 version: 4.0.2 

Signature file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
Chronic main manifest: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
Chronic main count matrix: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv
Step1 detectability summary: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step1_current78_signature_core_detectability_audit/chronic_step1_current78_signature_core_detectability_summary.csv exists=TRUE

Fixed signature genes: 75
Up genes: 65
Down genes: 10

Chronic main count matrix dimensions used: 22438 genes x 24 samples
Sample group counts:

   Control_52W ACLT_alone_52W 
            12             12 

Normalization summary:
                               metric                                  value
1            input_count_matrix_genes                                  22438
2          input_count_matrix_samples                                     24
3 filterByExpr_applied_before_scoring                                  FALSE
4                  edgeR_DGEList_used                                   TRUE
5                   TMM_normalization edgeR::calcNormFactors(method = 'TMM')
6                     logCPM_function                 edgeR::cpm(log = TRUE)
7                         prior_count                                      1
8                normalized_lib_sizes                                   TRUE
9            estimated_counts_rounded                                  FALSE

TMM library summary:
                                            sample_id          group raw_estimated_count_sum TMM_norm_factor effective_library_size detected_genes_count_gt0
GSM7140486_CON1_synovium     GSM7140486_CON1_synovium    Control_52W                 7108677       1.1173152                7942633                    16731
GSM7140491_CON2_synovium     GSM7140491_CON2_synovium    Control_52W                 8959039       1.0410712                9326997                    16672
GSM7140495_CON3_synovium     GSM7140495_CON3_synovium    Control_52W                 9719015       1.1126848               10814201                    16993
GSM7140500_CON4_synovium     GSM7140500_CON4_synovium    Control_52W                 8391254       1.0728326                9002412                    16672
GSM7140502_CON5_synovium     GSM7140502_CON5_synovium    Control_52W                12170231       0.9537912               11607859                    16538
GSM7140507_CON6_synovium     GSM7140507_CON6_synovium    Control_52W                16132991       1.1147778               17984701                    17390
GSM7140511_CON7_synovium     GSM7140511_CON7_synovium    Control_52W                11920148       1.1091978               13221802                    17215
GSM7140516_CON8_synovium     GSM7140516_CON8_synovium    Control_52W                 9996314       1.1221556               11217420                    16810
GSM7140518_CON9_synovium     GSM7140518_CON9_synovium    Control_52W                12285504       1.1005048               13520257                    17058
GSM7140523_CON10_synovium   GSM7140523_CON10_synovium    Control_52W                13966921       1.0709385               14957713                    16999
GSM7140527_CON11_synovium   GSM7140527_CON11_synovium    Control_52W                16354930       0.8895768               14548967                    16694
GSM7140532_CON12_synovium   GSM7140532_CON12_synovium    Control_52W                11603076       1.0880467               12624688                    16604
GSM7140487_ACLT1_synovium   GSM7140487_ACLT1_synovium ACLT_alone_52W                10334100       0.8887394                9184322                    16660
GSM7140488_ACLT2_synovium   GSM7140488_ACLT2_synovium ACLT_alone_52W                10521457       1.0006256               10528039                    16880
GSM7140489_ACLT3_synovium   GSM7140489_ACLT3_synovium ACLT_alone_52W                18669780       0.9213546               17201488                    16495
GSM7140490_ACLT4_synovium   GSM7140490_ACLT4_synovium ACLT_alone_52W                 9298837       0.9904546                9210076                    16367
GSM7140492_ACLT5_synovium   GSM7140492_ACLT5_synovium ACLT_alone_52W                 9066422       0.9958939                9029194                    16613
GSM7140493_ACLT6_synovium   GSM7140493_ACLT6_synovium ACLT_alone_52W                10831988       0.9033712                9785306                    16291
GSM7140494_ACLT7_synovium   GSM7140494_ACLT7_synovium ACLT_alone_52W                10627443       0.9275873                9857881                    16710
GSM7140496_ACLT8_synovium   GSM7140496_ACLT8_synovium ACLT_alone_52W                12200497       1.0253431               12509695                    16912
GSM7140497_ACLT9_synovium   GSM7140497_ACLT9_synovium ACLT_alone_52W                 7706256       1.0095390                7779765                    16517
GSM7140498_ACLT10_synovium GSM7140498_ACLT10_synovium ACLT_alone_52W                 9007486       1.0418801                9384721                    16511
GSM7140499_ACLT11_synovium GSM7140499_ACLT11_synovium ACLT_alone_52W                 8387249       0.7973198                6687320                    15741
GSM7140501_ACLT12_synovium GSM7140501_ACLT12_synovium ACLT_alone_52W                10569945       0.8191519                8658390                    15555

Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5A_current78_chronic_signature_score_TMM_aligned/Figure5A_current78_chronic_signature_scores_TMM_aligned.png
Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5A_current78_chronic_signature_score_TMM_aligned/Figure5A_current78_chronic_signature_scores_TMM_aligned.pdf

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_gene_metadata_used.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_all_gene_TMM_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_TMM_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_zscore_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_scores_by_sample.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_scores_long.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_score_group_comparisons.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_sample_metadata_used.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_plot_files.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_normalization_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_library_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_input_file_audit.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_score_run_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_versions_and_method_records.csv

===== Chronic Step2 current78 TMM-aligned signature score summary =====
                                                 metric
                                         method_version
                                   signature_input_file
                             chronic_main_manifest_file
                               chronic_main_counts_file
                                 n_signature_genes_used
                                   n_up_signature_genes
                                 n_down_signature_genes
                                              n_samples
                                          n_Control_52W
                                       n_ACLT_alone_52W
                                            score_scale
                            filterByExpr_before_scoring
                                      TMM_normalization
                               estimated_counts_rounded
                   primary_directional_score_definition
                audit_directional_score_sum_means_saved
                                             group_test
                   directional_score_median_Control_52W
                directional_score_median_ACLT_alone_52W
 directional_score_median_difference_case_minus_control
                       directional_score_wilcox_p_value
  directional_score_BH_FDR_across_all_score_comparisons
                                       main_score_table
                                 group_comparison_table
                                             figure_png
                                             figure_pdf
                                             output_dir
                                                                                                                                                                                      value
                                                                                                                       2026-05-16_chronic_step2_current78_TMM_aligned_with_pig_early_Step19
                                        E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
                                                                             E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
                                                             E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv
                                                                                                                                                                                         75
                                                                                                                                                                                         65
                                                                                                                                                                                         10
                                                                                                                                                                                         24
                                                                                                                                                                                         12
                                                                                                                                                                                         12
                                                                                                    edgeR TMM-normalized logCPM; row-wise z-score across 24 chronic main-comparison samples
                                                                                                                                                                                      FALSE
                                                                               edgeR::calcNormFactors(method = 'TMM'); edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
                                                                                                                                                                                      FALSE
                                                                                                                                                     (sum up z + sum -down z)/(n_up+n_down)
                                                                                                                       TRUE; directional_score_sum_means = up_score + down_score_reoriented
                                                                                                          two-sided Wilcoxon rank-sum test for ACLT_alone_52W vs Control_52W; exact = FALSE
                                                                                                                                                                         -0.139433548042642
                                                                                                                                                                          0.162449627277128
                                                                                                                                                                           0.30188317531977
                                                                                                                                                                       0.000900935596360016
                                                                                                                                                                        0.00401900971132003
        E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_scores_by_sample.csv
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_score_group_comparisons.csv
           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5A_current78_chronic_signature_score_TMM_aligned/Figure5A_current78_chronic_signature_scores_TMM_aligned.png
           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5A_current78_chronic_signature_score_TMM_aligned/Figure5A_current78_chronic_signature_scores_TMM_aligned.pdf
                                                                           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned

Signature direction counts:

Down_in_ACLR   Up_in_ACLR 
          10           65 

Group comparison statistics:
                    comparison                       score n_control n_case median_control median_case mean_control  mean_case median_difference_case_minus_control
 ACLT_alone_52W_vs_Control_52W           directional_score        12     12    -0.13943355   0.1624496   -0.1677294  0.1677294                            0.3018832
 ACLT_alone_52W_vs_Control_52W directional_score_sum_means        12     12    -0.30068253   0.3215883   -0.2964636  0.2964636                            0.6222708
 ACLT_alone_52W_vs_Control_52W                    up_score        12     12    -0.12959273   0.1460881   -0.1748195  0.1748195                            0.2756808
 ACLT_alone_52W_vs_Control_52W       down_score_reoriented        12     12    -0.05067284   0.1158261   -0.1216441  0.1216441                            0.1664990
 ACLT_alone_52W_vs_Control_52W              down_score_raw        12     12     0.05067284  -0.1158261    0.1216441 -0.1216441                           -0.1664990
 ACLT_alone_52W_vs_Control_52W      total_score_unoriented        12     12    -0.11677319   0.1151519   -0.1352910  0.1352910                            0.2319250
 mean_difference_case_minus_control wilcox_p_value BH_FDR_across_all_score_comparisons
                          0.3354589   0.0009009356                         0.004019010
                          0.5929272   0.0020095049                         0.004019010
                          0.3496390   0.0016520395                         0.004019010
                          0.2432882   0.1748533069                         0.174853307
                         -0.2432882   0.1748533069                         0.174853307
                          0.2705820   0.0051079054                         0.007661858

Version and method records:
                     item                                                                                                                        value
                R_version                                                                                                                        4.5.2
            edgeR_version                                                                                                                        4.8.2
          ggplot2_version                                                                                                                        4.0.2
     signature_definition                                             current78-derived pig early 75-gene signature; 65 Up_in_ACLR and 10 Down_in_ACLR
              score_input                                                       GSE228848 chronic main-comparison tximport gene-level estimated counts
            normalization edgeR::DGEList; edgeR::calcNormFactors(method = 'TMM'); edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
 filtering_before_scoring                                                   No filterByExpr before scoring, to avoid losing predefined signature genes
      row_standardization                                                              row-wise z-score across 24 chronic main-comparison samples only
    primary_score_formula                                                                   directional_score = (sum z_up + sum -z_down)/(n_up+n_down)
    chronic_analysis_role                                                            extension validation only; no chronic signature/core redefinition

Chronic Step2 current78 TMM-aligned signature score completed successfully.

Chronic Step2 current78 TMM-aligned signature score completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/logs/chronic_step2_current78_TMM_aligned_signature_score_summary_to_send_me.txt
===== CHRONIC STEP3 CURRENT78 TMM-ALIGNED SINGSCORE SENSITIVITY =====
Run time: 2026-05-16 13:47:34 CEST
Interpretation: singscore sensitivity analysis using fixed pig early current78 75-gene signature.
Input transformation: edgeR TMM-normalized logCPM; no filterByExpr before ranking/scoring.

Signature file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
Chronic main manifest: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
Chronic main count matrix: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv
Updated TMM-aligned Step2 z-score scores: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_scores_by_sample.csv

Normalization summary:
                                    metric                                  value
1                 input_count_matrix_genes                                  22438
2               input_count_matrix_samples                                     24
3    filterByExpr_applied_before_singscore                                  FALSE
4                       edgeR_DGEList_used                                   TRUE
5                        TMM_normalization edgeR::calcNormFactors(method = 'TMM')
6                          logCPM_function                 edgeR::cpm(log = TRUE)
7                              prior_count                                      1
8                     normalized_lib_sizes                                   TRUE
9                 estimated_counts_rounded                                  FALSE
10 n_signature_genes_available_for_scoring                                     75
11                    n_up_signature_genes                                     65
12                  n_down_signature_genes                                     10

TMM library summary:
                                            sample_id          group raw_estimated_count_sum TMM_norm_factor effective_library_size detected_genes_count_gt0
GSM7140486_CON1_synovium     GSM7140486_CON1_synovium    Control_52W                 7108677       1.1173152                7942633                    16731
GSM7140491_CON2_synovium     GSM7140491_CON2_synovium    Control_52W                 8959039       1.0410712                9326997                    16672
GSM7140495_CON3_synovium     GSM7140495_CON3_synovium    Control_52W                 9719015       1.1126848               10814201                    16993
GSM7140500_CON4_synovium     GSM7140500_CON4_synovium    Control_52W                 8391254       1.0728326                9002412                    16672
GSM7140502_CON5_synovium     GSM7140502_CON5_synovium    Control_52W                12170231       0.9537912               11607859                    16538
GSM7140507_CON6_synovium     GSM7140507_CON6_synovium    Control_52W                16132991       1.1147778               17984701                    17390
GSM7140511_CON7_synovium     GSM7140511_CON7_synovium    Control_52W                11920148       1.1091978               13221802                    17215
GSM7140516_CON8_synovium     GSM7140516_CON8_synovium    Control_52W                 9996314       1.1221556               11217420                    16810
GSM7140518_CON9_synovium     GSM7140518_CON9_synovium    Control_52W                12285504       1.1005048               13520257                    17058
GSM7140523_CON10_synovium   GSM7140523_CON10_synovium    Control_52W                13966921       1.0709385               14957713                    16999
GSM7140527_CON11_synovium   GSM7140527_CON11_synovium    Control_52W                16354930       0.8895768               14548967                    16694
GSM7140532_CON12_synovium   GSM7140532_CON12_synovium    Control_52W                11603076       1.0880467               12624688                    16604
GSM7140487_ACLT1_synovium   GSM7140487_ACLT1_synovium ACLT_alone_52W                10334100       0.8887394                9184322                    16660
GSM7140488_ACLT2_synovium   GSM7140488_ACLT2_synovium ACLT_alone_52W                10521457       1.0006256               10528039                    16880
GSM7140489_ACLT3_synovium   GSM7140489_ACLT3_synovium ACLT_alone_52W                18669780       0.9213546               17201488                    16495
GSM7140490_ACLT4_synovium   GSM7140490_ACLT4_synovium ACLT_alone_52W                 9298837       0.9904546                9210076                    16367
GSM7140492_ACLT5_synovium   GSM7140492_ACLT5_synovium ACLT_alone_52W                 9066422       0.9958939                9029194                    16613
GSM7140493_ACLT6_synovium   GSM7140493_ACLT6_synovium ACLT_alone_52W                10831988       0.9033712                9785306                    16291
GSM7140494_ACLT7_synovium   GSM7140494_ACLT7_synovium ACLT_alone_52W                10627443       0.9275873                9857881                    16710
GSM7140496_ACLT8_synovium   GSM7140496_ACLT8_synovium ACLT_alone_52W                12200497       1.0253431               12509695                    16912
GSM7140497_ACLT9_synovium   GSM7140497_ACLT9_synovium ACLT_alone_52W                 7706256       1.0095390                7779765                    16517
GSM7140498_ACLT10_synovium GSM7140498_ACLT10_synovium ACLT_alone_52W                 9007486       1.0418801                9384721                    16511
GSM7140499_ACLT11_synovium GSM7140499_ACLT11_synovium ACLT_alone_52W                 8387249       0.7973198                6687320                    15741
GSM7140501_ACLT12_synovium GSM7140501_ACLT12_synovium ACLT_alone_52W                10569945       0.8191519                8658390                    15555

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_signature_gene_metadata_used.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_TMM_aligned_genomewide_logCPM_matrix_for_singscore.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_TMM_aligned_normalization_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_TMM_aligned_library_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_scores_by_sample.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_scores_long.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_group_comparisons.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_vs_zscore_correlation.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_TMM_aligned_sample_metadata_used.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_TMM_aligned_plot_files.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_TMM_aligned_input_file_audit.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_run_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_TMM_aligned_versions_and_method_records.csv

===== Chronic Step3 current78 TMM-aligned singscore summary =====
                                                        metric
                                                method_version
                                          signature_input_file
                                    chronic_main_manifest_file
                                      chronic_main_counts_file
                         updated_TMM_aligned_step2_scores_file
                                        n_signature_genes_used
                                          n_up_signature_genes
                                        n_down_signature_genes
                                     n_genomewide_genes_ranked
                                                     n_samples
                                                 n_Control_52W
                                              n_ACLT_alone_52W
                                               singscore_scale
                                 filterByExpr_before_singscore
                                             TMM_normalization
                                      estimated_counts_rounded
                              directional_singscore_definition
                                  down_singscore_reorientation
                                                    group_test
                      directional_singscore_median_Control_52W
                   directional_singscore_median_ACLT_alone_52W
    directional_singscore_median_difference_case_minus_control
                          directional_singscore_wilcox_p_value
 directional_singscore_BH_FDR_across_all_singscore_comparisons
          TMM_aligned_directional_zscore_vs_singscore_spearman
                                          main_singscore_table
                                        group_comparison_table
                            zscore_singscore_correlation_table
                                                    figure_png
                                                    figure_pdf
                                                    output_dir
                                                                                                                                                                                            value
                                                                                                                          2026-05-16_chronic_step3_current78_TMM_aligned_singscore_sensitivity_v1
                                              E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
                                                                                   E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
                                                                   E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv
              E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_scores_by_sample.csv
                                                                                                                                                                                               75
                                                                                                                                                                                               65
                                                                                                                                                                                               10
                                                                                                                                                                                            22438
                                                                                                                                                                                               24
                                                                                                                                                                                               12
                                                                                                                                                                                               12
                                                                                                                       genome-wide sample-wise ranks from TMM-normalized logCPM; centerScore=TRUE
                                                                                                                                                                                            FALSE
                                                                                     edgeR::calcNormFactors(method = 'TMM'); edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
                                                                                                                                                                                            FALSE
                                                                                                   simpleScore(rankData, upSet=65 genes, downSet=10 genes, centerScore=TRUE, knownDirection=TRUE)
                                                                                                                            down genes scored as upSet to obtain raw score, then multiplied by -1
                                                                                                                two-sided Wilcoxon rank-sum test for ACLT_alone_52W vs Control_52W; exact = FALSE
                                                                                                                                                                                 0.10422263662957
                                                                                                                                                                                0.144578150066651
                                                                                                                                                                               0.0403555134370819
                                                                                                                                                                              0.00725995498819057
                                                                                                                                                                               0.0145199099763811
                                                                                                                                                                                0.732173913043478
        E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_scores_by_sample.csv
       E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_group_comparisons.csv
   E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/chronic_step3_current78_TMM_aligned_singscore_vs_zscore_correlation.csv
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Supplementary_chronic_current78_singscore_TMM_aligned/Supplementary_chronic_current78_TMM_aligned_singscore_scores_mouse_style.png
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Supplementary_chronic_current78_singscore_TMM_aligned/Supplementary_chronic_current78_TMM_aligned_singscore_scores_mouse_style.pdf
                                                                           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned

Normalization summary:
                                  metric                                  value
                input_count_matrix_genes                                  22438
              input_count_matrix_samples                                     24
   filterByExpr_applied_before_singscore                                  FALSE
                      edgeR_DGEList_used                                   TRUE
                       TMM_normalization edgeR::calcNormFactors(method = 'TMM')
                         logCPM_function                 edgeR::cpm(log = TRUE)
                             prior_count                                      1
                    normalized_lib_sizes                                   TRUE
                estimated_counts_rounded                                  FALSE
 n_signature_genes_available_for_scoring                                     75
                    n_up_signature_genes                                     65
                  n_down_signature_genes                                     10

Signature direction counts:

Down_in_ACLR   Up_in_ACLR 
          10           65 

Singscore group comparison statistics:
                    comparison                     score n_control n_case median_control  median_case mean_control    mean_case median_difference_case_minus_control
 ACLT_alone_52W_vs_Control_52W     directional_singscore        12     12     0.10422264  0.144578150   0.10184728  0.146780468                           0.04035551
 ACLT_alone_52W_vs_Control_52W              up_singscore        12     12     0.12671042  0.149780470   0.12614392  0.145571413                           0.02307005
 ACLT_alone_52W_vs_Control_52W down_singscore_reoriented        12     12    -0.02915998  0.003054218  -0.02429664  0.001209054                           0.03221420
 ACLT_alone_52W_vs_Control_52W        down_singscore_raw        12     12     0.02915998 -0.003054218   0.02429664 -0.001209054                          -0.03221420
 mean_difference_case_minus_control wilcox_p_value BH_FDR_across_all_singscore_comparisons
                         0.04493318    0.007259955                             0.014519910
                         0.01942749    0.001353941                             0.005415766
                         0.02550569    0.099877403                             0.099877403
                        -0.02550569    0.099877403                             0.099877403

TMM-aligned z-score vs singscore correlation:
                                      comparison n_samples  spearman   pearson
     TMM_aligned_directional_zscore_vs_singscore        24 0.7321739 0.7014277
              TMM_aligned_up_zscore_vs_singscore        24 0.9408696 0.9402832
 TMM_aligned_down_reoriented_zscore_vs_singscore        24 0.8017391 0.8513580

Version and method records:
                     item                                                                                  value
                R_version                                                                                  4.5.2
            edgeR_version                                                                                  4.8.2
        singscore_version                                                                                 1.30.0
          ggplot2_version                                                                                  4.0.2
     signature_definition       current78-derived pig early 75-gene signature; 65 Up_in_ACLR and 10 Down_in_ACLR
              score_input                 GSE228848 chronic main-comparison tximport gene-level estimated counts
               rank_input        genome-wide TMM-normalized logCPM matrix for 24 chronic main-comparison samples
 filtering_before_ranking No filterByExpr before singscore ranking/scoring, to retain predefined signature genes
              centerScore                                                                                   TRUE
           knownDirection                                                                                   TRUE
    chronic_analysis_role                      sensitivity analysis only; no chronic signature/core redefinition

Chronic Step3 current78 TMM-aligned singscore sensitivity completed successfully.

Chronic Step3 current78 TMM-aligned singscore sensitivity completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step3_current78_singscore_sensitivity_TMM_aligned/logs/chronic_step3_current78_TMM_aligned_singscore_sensitivity_summary_to_send_me.txt
===== CHRONIC STEP5 CURRENT78 TMM-ALIGNED CORE HEATMAP =====
Run time: 2026-05-16 14:55:51 CEST
Method version: 2026-05-16_chronic_step5_current78_core_heatmap_TMM_aligned_v1
Interpretation: chronic visualization of the fixed early-defined 24 core genes.
This step does NOT redefine a chronic core set.

Required packages loaded.
edgeR version: 4.8.2 
pheatmap version: 1.0.13 
RColorBrewer version: 1.1.3 

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_core_gene_file_candidates.csv
Core gene file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_gene_table.csv
Manifest file: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
Count matrix file: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv

Fixed early-defined core genes loaded: 24

null device 
          1 
null device 
          1 
Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned/Figure5C_current78_main24_core_heatmap_TMM_aligned.png
Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned/Figure5C_current78_main24_core_heatmap_TMM_aligned.pdf

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_fixed_early_defined_24_core_gene_metadata.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_all_gene_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_core_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_core_heatmap_zscore_matrix_unclipped.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_core_heatmap_zscore_matrix_clipped_for_plot.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_sample_metadata_used.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_core_gene_order_used_in_heatmap.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_normalization_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_library_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_plot_files.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_input_file_audit.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_core_heatmap_run_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/chronic_step5_TMM_aligned_versions_and_method_records.csv

===== Chronic Step5 current78 TMM-aligned core heatmap summary =====
                                        metric
                                method_version
                          core_gene_input_file
                                 manifest_file
                                    count_file
              n_fixed_early_defined_core_genes
 n_core_genes_detected_in_chronic_count_matrix
                                     n_samples
                                 n_Control_52W
                              n_ACLT_alone_52W
                           heatmap_input_scale
                   filterByExpr_before_heatmap
                      estimated_counts_rounded
                     heatmap_value_source_data
                         heatmap_value_plotted
                      zscore_clipping_for_plot
                            row_label_priority
                                row_order_rule
                             column_order_rule
                                  cluster_rows
                                  cluster_cols
                                    figure_png
                                    figure_pdf
                                    output_dir
                                                                                                                                                                   value
                                                                                                          2026-05-16_chronic_step5_current78_core_heatmap_TMM_aligned_v1
                     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_gene_table.csv
                                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv
                                                                                                                                                                      24
                                                                                                                                                                      24
                                                                                                                                                                      24
                                                                                                                                                                      12
                                                                                                                                                                      12
                                                                          edgeR TMM-normalized logCPM based on chronic main-comparison gene-level estimated count matrix
                                                                                                                                                                   FALSE
                                                                                                                                                                   FALSE
                                                                                                                    unclipped row-wise z-score across 24 chronic samples
                                                                                                                row-wise z-score clipped to [-2, 2] for heatmap plotting
                                                                                                           TRUE; source data include both unclipped and clipped matrices
                                                                                                                                      human_gene > pig_symbol > pig_ensg
                                                              Up_in_ACLR first, then Down_in_ACLR; within each group ordered by mean(case-control) TMM-logCPM difference
                                                                                                                  Control_52W samples first, then ACLT_alone_52W samples
                                                                                                                                                                   FALSE
                                                                                                                                                                   FALSE
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned/Figure5C_current78_main24_core_heatmap_TMM_aligned.png
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned/Figure5C_current78_main24_core_heatmap_TMM_aligned.pdf
                                                           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned

Normalization summary:
                      metric                                  value
    input_count_matrix_genes                                  22438
  input_count_matrix_samples                                     24
          edgeR_DGEList_used                                   TRUE
           TMM_normalization edgeR::calcNormFactors(method = 'TMM')
             logCPM_function                 edgeR::cpm(log = TRUE)
                 prior_count                                      1
        normalized_lib_sizes                                   TRUE
 filterByExpr_before_heatmap                                  FALSE
    estimated_counts_rounded                                  FALSE
          n_fixed_core_genes                                     24
 n_fixed_core_genes_detected                                     24

Core gene order preview:
 row_label           pig_ensg pig_symbol human_gene signature_direction mean_diff_case_minus_control_TMM_logCPM
      SDC1 ENSSSCG00000008601       SDC1       SDC1          Up_in_ACLR                              2.29904586
     SPARC ENSSSCG00000017082      SPARC      SPARC          Up_in_ACLR                              1.40856902
     LOXL2 ENSSSCG00000033608      LOXL2      LOXL2          Up_in_ACLR                              1.20370946
  SERPINE1 ENSSSCG00000025698   SERPINE1   SERPINE1          Up_in_ACLR                              1.02104059
     HTRA1 ENSSSCG00000010703      HTRA1      HTRA1          Up_in_ACLR                              0.71409080
     HIF1A ENSSSCG00000005096      HIF1A      HIF1A          Up_in_ACLR                              0.47509865
   RACGAP1 ENSSSCG00000000217    RACGAP1    RACGAP1          Up_in_ACLR                              0.46739836
     ESPL1 ENSSSCG00000000265      ESPL1      ESPL1          Up_in_ACLR                              0.45202430
     FOXM1 ENSSSCG00000000739      FOXM1      FOXM1          Up_in_ACLR                              0.39562084
     RUNX1 ENSSSCG00000035537      RUNX1      RUNX1          Up_in_ACLR                              0.36739886
      CDK1 ENSSSCG00000010214       CDK1       CDK1          Up_in_ACLR                              0.34695774
     KIFC1 ENSSSCG00000001510      KIFC1      KIFC1          Up_in_ACLR                              0.25240081
       TTK ENSSSCG00000004466        TTK        TTK          Up_in_ACLR                              0.25227359
     CKAP2 ENSSSCG00000009378      CKAP2      CKAP2          Up_in_ACLR                              0.24778850
     CCNB1 ENSSSCG00000029326      CCNB1      CCNB1          Up_in_ACLR                              0.22572418
      NEK2 ENSSSCG00000015604       NEK2       NEK2          Up_in_ACLR                              0.20096781
     KIF2C ENSSSCG00000003931      KIF2C      KIF2C          Up_in_ACLR                              0.15844107
      BUB1 ENSSSCG00000030469       BUB1       BUB1          Up_in_ACLR                              0.12895294
    KIF20A ENSSSCG00000014326     KIF20A     KIF20A          Up_in_ACLR                              0.07577240
    ENTPD7 ENSSSCG00000010540     ENTPD7     ENTPD7          Up_in_ACLR                             -0.01845322
   TNFSF15 ENSSSCG00000059646    TNFSF15    TNFSF15          Up_in_ACLR                             -0.23093169
      FASN ENSSSCG00000029944       FASN       FASN        Down_in_ACLR                             -0.29044648
       CBS ENSSSCG00000040779        CBS        CBS        Down_in_ACLR                             -1.07709796
     PPARG ENSSSCG00000011579      PPARG      PPARG        Down_in_ACLR                             -1.61707980

Version and method records:
                     item                                                                                            value
                R_version                                                                                            4.5.2
            edgeR_version                                                                                            4.8.2
         pheatmap_version                                                                                           1.0.13
     RColorBrewer_version                                                                                            1.1.3
    fixed_core_definition              fixed early-defined 24 core ortholog genes from pig early primary-analysis workflow
            heatmap_scale                                                                      edgeR TMM-normalized logCPM
 filtering_before_heatmap                       No filterByExpr before heatmap extraction, to retain predefined core genes
            heatmap_value row-wise z-score across 24 chronic main-comparison samples; clipped to [-2, 2] for plotting only
   display_label_priority                                                               human_gene > pig_symbol > pig_ensg
            analysis_role                               chronic extension visualization only; no chronic core redefinition

Chronic Step5 current78 TMM-aligned core heatmap completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned/logs/chronic_step5_current78_TMM_aligned_core_heatmap_summary_to_send_me.txt
===== CHRONIC STEP5 CURRENT78 TMM-ALIGNED CORE HEATMAP, PIG EARLY FIGURE4D STYLE =====
Run time: 2026-05-16 15:03:26 CEST
Method version: 2026-05-16_chronic_step5_TMM_aligned_core_heatmap_pig_early_style_v1
Interpretation: chronic visualization of fixed early-defined 24 core genes.
This step does NOT redefine a chronic core set.

Required packages loaded.
edgeR version: 4.8.2 
pheatmap version: 1.0.13 
RColorBrewer version: 1.1.3 

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_core_gene_file_candidates.csv
Core gene file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_gene_table.csv
Manifest file: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
Count matrix file: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv

Fixed early-defined core genes loaded: 24

Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style.pdf
Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style.png

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_fixed_early_defined_24_core_gene_metadata.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_all_gene_TMM_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_core_TMM_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_core_heatmap_zscore_matrix_unclipped.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_core_heatmap_zscore_matrix_clipped_for_plot.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_sample_metadata_used.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_normalization_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_TMM_normalization_factors.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_plot_files.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_input_file_audit.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_core_heatmap_run_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/chronic_step5_pig_early_style_versions_and_method_records.csv

===== Chronic Step5 current78 TMM-aligned pig-early-style core heatmap summary =====
                                        metric
                                method_version
                          core_gene_input_file
                                 manifest_file
                                    count_file
              n_fixed_early_defined_core_genes
 n_core_genes_detected_in_chronic_count_matrix
                                     n_samples
                                 n_Control_52W
                              n_ACLT_alone_52W
                           normalization_scope
                          normalization_method
                   filterByExpr_before_heatmap
                      estimated_counts_rounded
                           heatmap_color_scale
                     heatmap_value_source_data
                         heatmap_value_plotted
                          zscore_clip_for_plot
                        display_label_priority
                                row_clustering
                             column_clustering
                     row_dendrogram_treeheight
                             column_order_rule
                                      pdf_file
                                      png_file
                                    output_dir
                                          note
                                                                                                                                                                                                   value
                                                                                                                                    2026-05-16_chronic_step5_TMM_aligned_core_heatmap_pig_early_style_v1
                                                     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_gene_table.csv
                                                                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
                                                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv
                                                                                                                                                                                                      24
                                                                                                                                                                                                      24
                                                                                                                                                                                                      24
                                                                                                                                                                                                      12
                                                                                                                                                                                                      12
                                       Full chronic main-comparison gene-level estimated count matrix across the 24 plotted samples; 24 core genes are extracted after TMM-normalized logCPM calculation
                                                     edgeR::DGEList(full count matrix) -> edgeR::calcNormFactors(method = 'TMM') -> edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
                                                                                                                                                                                                   FALSE
                                                                                                                                                                                                   FALSE
                                                                                                                                             rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red
                                                                                                                                                    unclipped row-wise z-score across 24 chronic samples
                                                                                                                                               row-wise z-score clipped to [-2, 2] for pheatmap plotting
                                                                                                                                                                                                 [-2, 2]
                                                                                                                       pig_symbol > human_gene > pig_ensg, matching pig early Figure4D visual convention
                                                                                                                                                              TRUE; euclidean distance, complete linkage
                                                                                                                                                                                                   FALSE
                                                                                                                                                                                                      55
                                                                                                                                                  Control_52W samples first, then ACLT_alone_52W samples
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style.pdf
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style.png
                                                                           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style
         Style-aligned chronic Step5: TMM normalization is estimated from the full chronic main-comparison count matrix before extracting fixed 24 core genes; pig early Figure4D heatmap style is used.

Normalization summary:
                      metric                                  value
    input_count_matrix_genes                                  22438
  input_count_matrix_samples                                     24
          edgeR_DGEList_used                                   TRUE
           TMM_normalization edgeR::calcNormFactors(method = 'TMM')
             logCPM_function                 edgeR::cpm(log = TRUE)
                 prior_count                                      1
        normalized_lib_sizes                                   TRUE
 filterByExpr_before_heatmap                                  FALSE
    estimated_counts_rounded                                  FALSE
          n_fixed_core_genes                                     24
 n_fixed_core_genes_detected                                     24

Core genes preview:
 display_label human_gene           pig_ensg pig_symbol        MouseDirection signature_direction
           CBS        CBS ENSSSCG00000040779        CBS mouse_persistent_down        Down_in_ACLR
          FASN       FASN ENSSSCG00000029944       FASN mouse_persistent_down        Down_in_ACLR
         PPARG      PPARG ENSSSCG00000011579      PPARG mouse_persistent_down        Down_in_ACLR
          BUB1       BUB1 ENSSSCG00000030469       BUB1   mouse_persistent_up          Up_in_ACLR
         CCNB1      CCNB1 ENSSSCG00000029326      CCNB1   mouse_persistent_up          Up_in_ACLR
          CDK1       CDK1 ENSSSCG00000010214       CDK1   mouse_persistent_up          Up_in_ACLR
         CKAP2      CKAP2 ENSSSCG00000009378      CKAP2   mouse_persistent_up          Up_in_ACLR
        ENTPD7     ENTPD7 ENSSSCG00000010540     ENTPD7   mouse_persistent_up          Up_in_ACLR
         ESPL1      ESPL1 ENSSSCG00000000265      ESPL1   mouse_persistent_up          Up_in_ACLR
         FOXM1      FOXM1 ENSSSCG00000000739      FOXM1   mouse_persistent_up          Up_in_ACLR
         HIF1A      HIF1A ENSSSCG00000005096      HIF1A   mouse_persistent_up          Up_in_ACLR
         HTRA1      HTRA1 ENSSSCG00000010703      HTRA1   mouse_persistent_up          Up_in_ACLR
        KIF20A     KIF20A ENSSSCG00000014326     KIF20A   mouse_persistent_up          Up_in_ACLR
         KIF2C      KIF2C ENSSSCG00000003931      KIF2C   mouse_persistent_up          Up_in_ACLR
         KIFC1      KIFC1 ENSSSCG00000001510      KIFC1   mouse_persistent_up          Up_in_ACLR
         LOXL2      LOXL2 ENSSSCG00000033608      LOXL2   mouse_persistent_up          Up_in_ACLR
          NEK2       NEK2 ENSSSCG00000015604       NEK2   mouse_persistent_up          Up_in_ACLR
       RACGAP1    RACGAP1 ENSSSCG00000000217    RACGAP1   mouse_persistent_up          Up_in_ACLR
         RUNX1      RUNX1 ENSSSCG00000035537      RUNX1   mouse_persistent_up          Up_in_ACLR
          SDC1       SDC1 ENSSSCG00000008601       SDC1   mouse_persistent_up          Up_in_ACLR
      SERPINE1   SERPINE1 ENSSSCG00000025698   SERPINE1   mouse_persistent_up          Up_in_ACLR
         SPARC      SPARC ENSSSCG00000017082      SPARC   mouse_persistent_up          Up_in_ACLR
       TNFSF15    TNFSF15 ENSSSCG00000059646    TNFSF15   mouse_persistent_up          Up_in_ACLR
           TTK        TTK ENSSSCG00000004466        TTK   mouse_persistent_up          Up_in_ACLR

Sample annotation:
                  sample_id      group_raw          group
   GSM7140486_CON1_synovium    Control_52W    Control_52W
   GSM7140491_CON2_synovium    Control_52W    Control_52W
   GSM7140495_CON3_synovium    Control_52W    Control_52W
   GSM7140500_CON4_synovium    Control_52W    Control_52W
   GSM7140502_CON5_synovium    Control_52W    Control_52W
   GSM7140507_CON6_synovium    Control_52W    Control_52W
   GSM7140511_CON7_synovium    Control_52W    Control_52W
   GSM7140516_CON8_synovium    Control_52W    Control_52W
   GSM7140518_CON9_synovium    Control_52W    Control_52W
  GSM7140523_CON10_synovium    Control_52W    Control_52W
  GSM7140527_CON11_synovium    Control_52W    Control_52W
  GSM7140532_CON12_synovium    Control_52W    Control_52W
  GSM7140487_ACLT1_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140488_ACLT2_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140489_ACLT3_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140490_ACLT4_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140492_ACLT5_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140493_ACLT6_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140494_ACLT7_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140496_ACLT8_synovium ACLT_alone_52W ACLT_alone_52W
  GSM7140497_ACLT9_synovium ACLT_alone_52W ACLT_alone_52W
 GSM7140498_ACLT10_synovium ACLT_alone_52W ACLT_alone_52W
 GSM7140499_ACLT11_synovium ACLT_alone_52W ACLT_alone_52W
 GSM7140501_ACLT12_synovium ACLT_alone_52W ACLT_alone_52W

Version and method records:
                     item                                                                                            value
                R_version                                                                                            4.5.2
            edgeR_version                                                                                            4.8.2
         pheatmap_version                                                                                           1.0.13
     RColorBrewer_version                                                                                            1.1.3
    fixed_core_definition              fixed early-defined 24 core ortholog genes from pig early primary-analysis workflow
            heatmap_scale                                                                      edgeR TMM-normalized logCPM
 filtering_before_heatmap                       No filterByExpr before heatmap extraction, to retain predefined core genes
            heatmap_value row-wise z-score across 24 chronic main-comparison samples; clipped to [-2, 2] for plotting only
   display_label_priority                                                               pig_symbol > human_gene > pig_ensg
             visual_style               pig early Figure4D Step22 style: rev(RdYlBu), row dendrogram, no column clustering
            analysis_role                               chronic extension visualization only; no chronic core redefinition

Chronic Step5 current78 TMM-aligned pig-early-style core heatmap completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style/logs/chronic_step5_TMM_aligned_pig_early_style_core_heatmap_summary_to_send_me.txt
===== CHRONIC STEP5 CURRENT78 TMM-ALIGNED CORE HEATMAP, PIG EARLY FIGURE4D STYLE =====
Run time: 2026-05-16 15:07:17 CEST
Method version: 2026-05-16_chronic_step5_TMM_aligned_core_heatmap_pig_early_style_compact_v2
Interpretation: chronic visualization of fixed early-defined 24 core genes.
This step does NOT redefine a chronic core set.

Required packages loaded.
edgeR version: 4.8.2 
pheatmap version: 1.0.13 
RColorBrewer version: 1.1.3 

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_core_gene_file_candidates.csv
Core gene file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_gene_table.csv
Manifest file: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
Count matrix file: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv

Fixed early-defined core genes loaded: 24

Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact.pdf
Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact.png

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_fixed_early_defined_24_core_gene_metadata.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_all_gene_TMM_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_core_TMM_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_core_heatmap_zscore_matrix_unclipped.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_core_heatmap_zscore_matrix_clipped_for_plot.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_sample_metadata_used.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_normalization_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_TMM_normalization_factors.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_plot_files.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_input_file_audit.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_core_heatmap_run_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/chronic_step5_pig_early_style_versions_and_method_records.csv

===== Chronic Step5 current78 TMM-aligned pig-early-style core heatmap summary =====
                                        metric
                                method_version
                          core_gene_input_file
                                 manifest_file
                                    count_file
              n_fixed_early_defined_core_genes
 n_core_genes_detected_in_chronic_count_matrix
                                     n_samples
                                 n_Control_52W
                              n_ACLT_alone_52W
                           normalization_scope
                          normalization_method
                   filterByExpr_before_heatmap
                      estimated_counts_rounded
                           heatmap_color_scale
                     heatmap_value_source_data
                         heatmap_value_plotted
                          zscore_clip_for_plot
                        display_label_priority
                                row_clustering
                             column_clustering
                     row_dendrogram_treeheight
                             column_order_rule
                                      pdf_file
                                      png_file
                                    output_dir
                                          note
                                                                                                                                                                                                                   value
                                                                                                                                            2026-05-16_chronic_step5_TMM_aligned_core_heatmap_pig_early_style_compact_v2
                                                                     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_gene_table.csv
                                                                                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_manifest.csv
                                                                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step25v3_pig_chronic_main_comparison_gene_level_counts_matrix.csv
                                                                                                                                                                                                                      24
                                                                                                                                                                                                                      24
                                                                                                                                                                                                                      24
                                                                                                                                                                                                                      12
                                                                                                                                                                                                                      12
                                                       Full chronic main-comparison gene-level estimated count matrix across the 24 plotted samples; 24 core genes are extracted after TMM-normalized logCPM calculation
                                                                     edgeR::DGEList(full count matrix) -> edgeR::calcNormFactors(method = 'TMM') -> edgeR::cpm(log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
                                                                                                                                                                                                                   FALSE
                                                                                                                                                                                                                   FALSE
                                                                                                                                                             rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red
                                                                                                                                                                    unclipped row-wise z-score across 24 chronic samples
                                                                                                                                                               row-wise z-score clipped to [-2, 2] for pheatmap plotting
                                                                                                                                                                                                                 [-2, 2]
                                                                                                                                       pig_symbol > human_gene > pig_ensg, matching pig early Figure4D visual convention
                                                                                                                                                                              TRUE; euclidean distance, complete linkage
                                                                                                                                                                                                                   FALSE
                                                                                                                                                                                                                      55
                                                                                                                                                                  Control_52W samples first, then ACLT_alone_52W samples
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact.pdf
 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact/Figure5C_current78_main24_core_heatmap_TMM_aligned_pig_early_style_compact.png
                                                                                   E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact
                         Style-aligned chronic Step5: TMM normalization is estimated from the full chronic main-comparison count matrix before extracting fixed 24 core genes; pig early Figure4D heatmap style is used.

Normalization summary:
                      metric                                  value
    input_count_matrix_genes                                  22438
  input_count_matrix_samples                                     24
          edgeR_DGEList_used                                   TRUE
           TMM_normalization edgeR::calcNormFactors(method = 'TMM')
             logCPM_function                 edgeR::cpm(log = TRUE)
                 prior_count                                      1
        normalized_lib_sizes                                   TRUE
 filterByExpr_before_heatmap                                  FALSE
    estimated_counts_rounded                                  FALSE
          n_fixed_core_genes                                     24
 n_fixed_core_genes_detected                                     24

Core genes preview:
 display_label human_gene           pig_ensg pig_symbol        MouseDirection signature_direction
           CBS        CBS ENSSSCG00000040779        CBS mouse_persistent_down        Down_in_ACLR
          FASN       FASN ENSSSCG00000029944       FASN mouse_persistent_down        Down_in_ACLR
         PPARG      PPARG ENSSSCG00000011579      PPARG mouse_persistent_down        Down_in_ACLR
          BUB1       BUB1 ENSSSCG00000030469       BUB1   mouse_persistent_up          Up_in_ACLR
         CCNB1      CCNB1 ENSSSCG00000029326      CCNB1   mouse_persistent_up          Up_in_ACLR
          CDK1       CDK1 ENSSSCG00000010214       CDK1   mouse_persistent_up          Up_in_ACLR
         CKAP2      CKAP2 ENSSSCG00000009378      CKAP2   mouse_persistent_up          Up_in_ACLR
        ENTPD7     ENTPD7 ENSSSCG00000010540     ENTPD7   mouse_persistent_up          Up_in_ACLR
         ESPL1      ESPL1 ENSSSCG00000000265      ESPL1   mouse_persistent_up          Up_in_ACLR
         FOXM1      FOXM1 ENSSSCG00000000739      FOXM1   mouse_persistent_up          Up_in_ACLR
         HIF1A      HIF1A ENSSSCG00000005096      HIF1A   mouse_persistent_up          Up_in_ACLR
         HTRA1      HTRA1 ENSSSCG00000010703      HTRA1   mouse_persistent_up          Up_in_ACLR
        KIF20A     KIF20A ENSSSCG00000014326     KIF20A   mouse_persistent_up          Up_in_ACLR
         KIF2C      KIF2C ENSSSCG00000003931      KIF2C   mouse_persistent_up          Up_in_ACLR
         KIFC1      KIFC1 ENSSSCG00000001510      KIFC1   mouse_persistent_up          Up_in_ACLR
         LOXL2      LOXL2 ENSSSCG00000033608      LOXL2   mouse_persistent_up          Up_in_ACLR
          NEK2       NEK2 ENSSSCG00000015604       NEK2   mouse_persistent_up          Up_in_ACLR
       RACGAP1    RACGAP1 ENSSSCG00000000217    RACGAP1   mouse_persistent_up          Up_in_ACLR
         RUNX1      RUNX1 ENSSSCG00000035537      RUNX1   mouse_persistent_up          Up_in_ACLR
          SDC1       SDC1 ENSSSCG00000008601       SDC1   mouse_persistent_up          Up_in_ACLR
      SERPINE1   SERPINE1 ENSSSCG00000025698   SERPINE1   mouse_persistent_up          Up_in_ACLR
         SPARC      SPARC ENSSSCG00000017082      SPARC   mouse_persistent_up          Up_in_ACLR
       TNFSF15    TNFSF15 ENSSSCG00000059646    TNFSF15   mouse_persistent_up          Up_in_ACLR
           TTK        TTK ENSSSCG00000004466        TTK   mouse_persistent_up          Up_in_ACLR

Sample annotation:
                  sample_id      group_raw          group plot_label
   GSM7140486_CON1_synovium    Control_52W    Control_52W       CON1
   GSM7140491_CON2_synovium    Control_52W    Control_52W       CON2
   GSM7140495_CON3_synovium    Control_52W    Control_52W       CON3
   GSM7140500_CON4_synovium    Control_52W    Control_52W       CON4
   GSM7140502_CON5_synovium    Control_52W    Control_52W       CON5
   GSM7140507_CON6_synovium    Control_52W    Control_52W       CON6
   GSM7140511_CON7_synovium    Control_52W    Control_52W       CON7
   GSM7140516_CON8_synovium    Control_52W    Control_52W       CON8
   GSM7140518_CON9_synovium    Control_52W    Control_52W       CON9
  GSM7140523_CON10_synovium    Control_52W    Control_52W      CON10
  GSM7140527_CON11_synovium    Control_52W    Control_52W      CON11
  GSM7140532_CON12_synovium    Control_52W    Control_52W      CON12
  GSM7140487_ACLT1_synovium ACLT_alone_52W ACLT_alone_52W      ACLT1
  GSM7140488_ACLT2_synovium ACLT_alone_52W ACLT_alone_52W      ACLT2
  GSM7140489_ACLT3_synovium ACLT_alone_52W ACLT_alone_52W      ACLT3
  GSM7140490_ACLT4_synovium ACLT_alone_52W ACLT_alone_52W      ACLT4
  GSM7140492_ACLT5_synovium ACLT_alone_52W ACLT_alone_52W      ACLT5
  GSM7140493_ACLT6_synovium ACLT_alone_52W ACLT_alone_52W      ACLT6
  GSM7140494_ACLT7_synovium ACLT_alone_52W ACLT_alone_52W      ACLT7
  GSM7140496_ACLT8_synovium ACLT_alone_52W ACLT_alone_52W      ACLT8
  GSM7140497_ACLT9_synovium ACLT_alone_52W ACLT_alone_52W      ACLT9
 GSM7140498_ACLT10_synovium ACLT_alone_52W ACLT_alone_52W     ACLT10
 GSM7140499_ACLT11_synovium ACLT_alone_52W ACLT_alone_52W     ACLT11
 GSM7140501_ACLT12_synovium ACLT_alone_52W ACLT_alone_52W     ACLT12

Version and method records:
                     item                                                                                                               value
                R_version                                                                                                               4.5.2
            edgeR_version                                                                                                               4.8.2
         pheatmap_version                                                                                                              1.0.13
     RColorBrewer_version                                                                                                               1.1.3
    fixed_core_definition                                 fixed early-defined 24 core ortholog genes from pig early primary-analysis workflow
            heatmap_scale                                                                                         edgeR TMM-normalized logCPM
 filtering_before_heatmap                                          No filterByExpr before heatmap extraction, to retain predefined core genes
            heatmap_value                    row-wise z-score across 24 chronic main-comparison samples; clipped to [-2, 2] for plotting only
   display_label_priority                                                                                  pig_symbol > human_gene > pig_ensg
             visual_style                  pig early Figure4D Step22 style: rev(RdYlBu), row dendrogram, no column clustering; compact layout
      compact_plot_labels Column display labels shortened to CON1-CON12 and ACLT1-ACLT12 for plotting; raw sample IDs retained in source data
            analysis_role                                                  chronic extension visualization only; no chronic core redefinition

Chronic Step5 current78 TMM-aligned pig-early-style core heatmap completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5_current78_core_heatmap_TMM_aligned_pig_early_style_compact/logs/chronic_step5_TMM_aligned_pig_early_style_core_heatmap_summary_to_send_me.txt
===== CHRONIC STEP6 CURRENT78 FIGURE5D EARLY VS CHRONIC SUMMARY =====
Run time: 2026-05-16 15:46:35 CEST
Purpose: generate Figure5D as an early-vs-chronic summary comparison panel.

Early signature statistics: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_score_group_comparisons.csv
Early target GSEA: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/step20_current78_target_pathway_results_from_full_hallmark.csv
Early core summary: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_early_signature_remap_summary.csv
Chronic signature statistics: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_score_group_comparisons.csv
Chronic target GSEA: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step4_current78_DE_GSEA/chronic_step4_target_pathway_results_from_full_hallmark.csv
Chronic core audit summary: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5B_current78_core_DE_audit/chronic_step5B_early_defined_24_core_chronic_DE_audit_summary.csv

Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5D_current78_early_vs_chronic_summary/Figure5D_current78_early_vs_chronic_summary_tileplot.png
Saved plot: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5D_current78_early_vs_chronic_summary/Figure5D_current78_early_vs_chronic_summary_tileplot.pdf

Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary/chronic_step6_current78_figure5D_plot_source_data.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary/chronic_step6_current78_early_vs_chronic_extracted_values.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary/chronic_step6_current78_figure5D_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary/chronic_step6_current78_versions_and_method_records.csv

===== Chronic Step6 Figure5D summary =====
                       metric
   early_signature_stats_file
      early_gsea_targets_file
    early_step18_summary_file
 chronic_signature_stats_file
    chronic_gsea_targets_file
  chronic_step5B_summary_file
                   figure_png
                   figure_pdf
        plot_source_data_file
        extracted_values_file
        interpretation_policy
                   output_dir
                                                                                                                                                                                                     value
                                                E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_score_group_comparisons.csv
                                        E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/step20_current78_target_pathway_results_from_full_hallmark.csv
                                              E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_early_signature_remap_summary.csv
                E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score_TMM_aligned/chronic_step2_current78_TMM_aligned_signature_score_group_comparisons.csv
                                                  E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step4_current78_DE_GSEA/chronic_step4_target_pathway_results_from_full_hallmark.csv
                                     E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step5B_current78_core_DE_audit/chronic_step5B_early_defined_24_core_chronic_DE_audit_summary.csv
                                        E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5D_current78_early_vs_chronic_summary/Figure5D_current78_early_vs_chronic_summary_tileplot.png
                                        E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/figures/Figure5D_current78_early_vs_chronic_summary/Figure5D_current78_early_vs_chronic_summary_tileplot.pdf
                              E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary/chronic_step6_current78_figure5D_plot_source_data.csv
                      E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary/chronic_step6_current78_early_vs_chronic_extracted_values.csv
 Figure5D is a summary-comparison panel integrating completed early and chronic outputs under the fixed current78 / 75-gene / 24-core framework. It does not redefine a chronic signature or chronic core.
                                                                                    E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary

Extracted comparison values:
                                 category      value
     early_signature_directional_score_t7  0.9703711
    early_signature_directional_score_t28  0.6914711
      chronic_signature_directional_score  0.3018832
                early_inflammatory_NES_t7  1.3623498
               early_inflammatory_NES_t28  1.0118053
                 chronic_inflammatory_NES -1.0259994
                         early_EMT_NES_t7  1.6329280
                        early_EMT_NES_t28  1.5134305
                          chronic_EMT_NES  2.1476422
           early_signature_detected_genes 75.0000000
                          early_strict_t7 34.0000000
                         early_strict_t28 25.0000000
                         early_core_genes 24.0000000
                chronic_core_genes_tested 20.0000000
        chronic_core_direction_consistent 19.0000000
 chronic_core_strict_direction_consistent  4.0000000

Tile plot source data:
       stage                               metric                                      label fill_color numeric_value    fdr_value
    Early 1W         Signature\ndirectional score                 Δscore = 0.97\nFDR = 0.005    #E49A4A     0.9703711 5.074868e-03
    Early 1W Hallmark NES:\nInflammatory response                    NES = 1.36\nFDR = 0.030    #E49A4A     1.3623498 3.045636e-02
    Early 1W  Hallmark NES:\nEMT / ECM remodeling                    NES = 1.63\nFDR = 0.001    #E49A4A     1.6329280 1.123435e-03
    Early 1W                    Core-gene\nstatus                      Strict at 1W\n34 / 75    #E49A4A     0.4533333           NA
    Early 4W         Signature\ndirectional score                 Δscore = 0.69\nFDR = 0.005    #E49A4A     0.6914711 5.074868e-03
    Early 4W Hallmark NES:\nInflammatory response                    NES = 1.01\nFDR = 0.528    #F7E3CC     1.0118053 5.278643e-01
    Early 4W  Hallmark NES:\nEMT / ECM remodeling                    NES = 1.51\nFDR = 0.003    #E49A4A     1.5134305 3.483966e-03
    Early 4W                    Core-gene\nstatus                      Strict at 4W\n25 / 75    #F2BF84     0.3333333           NA
 Chronic 52W         Signature\ndirectional score                 Δscore = 0.30\nFDR = 0.004    #E49A4A     0.3018832 4.019010e-03
 Chronic 52W Hallmark NES:\nInflammatory response                   NES = -1.03\nFDR = 0.415    #DCEAF7    -1.0259994 4.154154e-01
 Chronic 52W  Hallmark NES:\nEMT / ECM remodeling                    NES = 2.15\nFDR < 0.001    #E49A4A     2.1476422 3.364377e-09
 Chronic 52W                    Core-gene\nstatus Strict retained\n4 / 24\n(19/20 dir-cons.)    #F7E3CC     0.1666667           NA
                                                                                     note
                                                    Early 1W directional score vs Control
                                                           Hallmark inflammatory response
                                                            Hallmark EMT / ECM remodeling
                         Among the fixed 75 pig signature genes, genes strict at early 1W
                                                    Early 4W directional score vs Control
                                                           Hallmark inflammatory response
                                                            Hallmark EMT / ECM remodeling
                         Among the fixed 75 pig signature genes, genes strict at early 4W
                                                 Chronic 52W directional score vs Control
                                                           Hallmark inflammatory response
                                                            Hallmark EMT / ECM remodeling
 Among the fixed early-defined 24 core genes, strict + direction-consistent in chronic DE

Version and method records:
                          item
                     R_version
               ggplot2_version
   early_signature_source_step
        early_GSEA_source_step
 chronic_signature_source_step
      chronic_GSEA_source_step
        core_comparison_policy
                                                                                                                                                                        value
                                                                                                                                                                        4.5.2
                                                                                                                                                                        4.0.2
                                                                                                                                   Step19 current78 pig early signature score
                                                                                                                                Step20 current78 pig early full-Hallmark GSEA
                                                                                                                          Chronic Step2 current78 TMM-aligned signature score
                                                                                                                    Chronic Step4 current78 edgeR QLF DE + full-Hallmark GSEA
 Compare early strict counts among the fixed 75-gene signature with chronic strict retained counts among the fixed early-defined 24 core genes; no chronic core redefinition.

Chronic Step6 / Figure5D current78 early-vs-chronic summary completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step6_current78_figure5D_early_vs_chronic_summary/logs/chronic_step6_current78_figure5D_summary_to_send_me.txt
===== STEP20_CURRENT78 PIG EARLY HALLMARK GSEA =====
Run time: 2026-05-16 20:21:30 
Base directory: E:/R/ACLsenescence2/rebuild_submission/02_pig_early 
Method: full Sus scrofa Hallmark collection fgsea; target pathways extracted after full run.

Input DE tables:
t7 : E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t7_vs_CON_QLF.csv 
t28: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t28_vs_CON_QLF.csv 

Gene ID to symbol mapping from available DE-table annotation:
n mapping rows: 28235 
n unambiguous gene_id symbols: 14304 
===== STEP20_CURRENT78 PIG EARLY HALLMARK GSEA =====
Run time: 2026-05-16 20:28:52 
Base directory: E:/R/ACLsenescence2/rebuild_submission/02_pig_early 
Method: full Sus scrofa Hallmark collection fgsea; target pathways extracted after full run.

Input DE tables:
t7 : E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t7_vs_CON_QLF.csv 
t28: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t28_vs_CON_QLF.csv 

Gene ID to symbol mapping from available DE-table annotation:
n mapping rows: 28235 
n unambiguous gene_id symbols: 14304 

Building ranked gene list for t7...

===== Step H2w6: trimmed branch edgeR QLF DE =====
[1] TRUE
2026-05-16 20:38:56 | Removed previous Step H2w6 output files: 13 
2026-05-16 20:38:56 | Project root: E:/R/ACLsenescence2 
2026-05-16 20:38:56 | Branch dir: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity 
2026-05-16 20:38:56 | Counts file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_gene_level_counts_matrix.csv 
2026-05-16 20:38:56 | Sample info file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_sample_info_for_DE.csv 
2026-05-16 20:38:56 | edgeR version: 4.8.2 
2026-05-16 20:38:57 | Gene-id column detected: gene_id 
2026-05-16 20:38:57 | Input genes: 35682 
2026-05-16 20:38:57 | Input samples: 18 
2026-05-16 20:38:57 | Sample groups: 6; 6; 6 
2026-05-16 20:38:57 | Starting contrast: ACLT_t7_vs_Control_t0 
2026-05-16 20:39:05 | Finished contrast: ACLT_t7_vs_Control_t0 | tested genes: 14458 | FDR<0.05: 6088 | strict: 2038 
2026-05-16 20:39:05 | Starting contrast: ACLT_t28_vs_Control_t0 
2026-05-16 20:39:12 | Finished contrast: ACLT_t28_vs_Control_t0 | tested genes: 14217 | FDR<0.05: 6005 | strict: 1897 

===== Step H2w6 contrast summary =====
             comparison n_control_samples n_treatment_samples n_tested_genes_after_filterByExpr n_FDR_lt_0_05 n_strict_FDR_0_05_abslogFC_1 n_strict_up n_strict_down
  ACLT_t7_vs_Control_t0                 6                   6                             14458          6088                         2038         908          1130
 ACLT_t28_vs_Control_t0                 6                   6                             14217          6005                         1897         672          1225

===== Step H2w6 overall summary =====
                           metric                                                                                                                              value
                     project_root                                                                                                                E:/R/ACLsenescence2
                       branch_dir                                                              E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity
                      counts_file      E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_gene_level_counts_matrix.csv
                 sample_info_file            E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_sample_info_for_DE.csv
                h2w5_summary_file E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w5_v2_trimmed_featureCounts_overall_summary.csv
                    edgeR_version                                                                                                                              4.8.2
                        n_samples                                                                                                                                 18
                    n_input_genes                                                                                                                              35682
                      n_contrasts                                                                                                                                  2
                n_t7_tested_genes                                                                                                                              14458
                 n_t7_FDR_lt_0_05                                                                                                                               6088
  n_t7_strict_FDR_0_05_abslogFC_1                                                                                                                               2038
               n_t28_tested_genes                                                                                                                              14217
                n_t28_FDR_lt_0_05                                                                                                                               6005
 n_t28_strict_FDR_0_05_abslogFC_1                                                                                                                               1897
                    batch_ok_h2w6                                                                                                                               TRUE
              sample_library_file  E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_count_library_summary_18samples.csv
                 run_summary_file         E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_edgeR_QLF_DE_run_summary.csv
               method_params_file      E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_edgeR_QLF_method_parameters.csv
                     summary_file        E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_edgeR_QLF_overall_summary.csv
                       config_rds                E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/objects/stepH2w6_trimmed_edgeR_QLF_config.rds
                         log_file                   E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/logs/stepH2w6_trimmed_edgeR_QLF_DE_log.txt
                           status                                                                                                   Step H2w6 completed successfully

输出文件：
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_edgeR_QLF_DE_run_summary.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_edgeR_QLF_overall_summary.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_edgeR_QLF_method_parameters.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_count_library_summary_18samples.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/logs/stepH2w6_trimmed_edgeR_QLF_DE_log.txt 

结果：H2w6 成功。trimmed branch 已完成 edgeR QLF DE。下一步可以比较 untrimmed 与 trimmed 分支的 DE/score/GSEA 稳定性。

===== End of Step H2w6 =====
[1] TRUE
2026-05-16 20:40:12 | ===== Step H2w7 started ===== 
2026-05-16 20:40:12 | project_root: E:/R/ACLsenescence2 
2026-05-16 20:40:12 | untrimmed_tables_dir: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables 
2026-05-16 20:40:12 | trimmed_tables_dir: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables 
2026-05-16 20:40:17 | Wrote candidate file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_file_selection_candidates.csv 
2026-05-16 20:40:17 | Wrote selected file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_file_selection_selected_files.csv 
2026-05-16 20:40:17 | Comparing contrast: ACLT_t7_vs_Control_t0 
2026-05-16 20:40:17 |   untrimmed_file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t7_vs_CON_QLF.csv 
2026-05-16 20:40:17 |   trimmed_file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t7_vs_Control_t0_QLF.csv 
2026-05-16 20:40:19 |   wrote gene comparison: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_stability_ACLT_t7_vs_Control_t0_gene_level_comparison.csv 
2026-05-16 20:40:19 |   wrote strict membership: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_stability_ACLT_t7_vs_Control_t0_strict_membership.csv 
2026-05-16 20:40:19 | Comparing contrast: ACLT_t28_vs_Control_t0 
2026-05-16 20:40:19 |   untrimmed_file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t28_vs_CON_QLF.csv 
2026-05-16 20:40:19 |   trimmed_file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t28_vs_Control_t0_QLF.csv 
2026-05-16 20:40:21 |   wrote gene comparison: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_stability_ACLT_t28_vs_Control_t0_gene_level_comparison.csv 
2026-05-16 20:40:21 |   wrote strict membership: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_stability_ACLT_t28_vs_Control_t0_strict_membership.csv 

===== Step H2w7 DE stability summary =====
               contrast                                                                                    untrimmed_file
  ACLT_t7_vs_Control_t0  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t7_vs_CON_QLF.csv
 ACLT_t28_vs_Control_t0 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t28_vs_CON_QLF.csv
                                                                                                                    trimmed_file n_untrimmed_genes n_trimmed_genes n_common_genes
  E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t7_vs_Control_t0_QLF.csv             14233           14458          14228
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t28_vs_Control_t0_QLF.csv             14002           14217          13993
 pearson_logFC_common spearman_logFC_common pearson_minuslog10P_common n_strict_untrimmed n_strict_trimmed n_strict_intersection n_strict_union strict_jaccard
             0.998779              0.998895                   0.996638               1990             2038                  1935           2093       0.924510
             0.999072              0.998902                   0.996977               1838             1897                  1805           1930       0.935233
 pct_untrimmed_strict_recovered_in_trimmed pct_trimmed_strict_also_in_untrimmed direction_agreement_pct_common direction_agreement_pct_strict_union median_abs_logFC_diff_common
                                    97.236                               94.946                         98.925                                  100                      0.01238
                                    98.205                               95.150                         98.835                                  100                      0.01209
 p95_abs_logFC_diff_common top100_FDR_overlap top500_FDR_overlap
                  0.047724                 95                476
                  0.044982                 94                477

===== Step H2w7 interpretation flags =====
               contrast                stability_flag spearman_logFC_common strict_jaccard direction_agreement_pct_strict_union
  ACLT_t7_vs_Control_t0 stable_with_minor_differences              0.998895       0.924510                                  100
 ACLT_t28_vs_Control_t0 stable_with_minor_differences              0.998902       0.935233                                  100

===== Step H2w7 overall summary =====
                                   metric                                                                                                                                            value
                             project_root                                                                                                                              E:/R/ACLsenescence2
                               branch_dir                                                                            E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity
                     untrimmed_tables_dir                                                                                       E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables
                       trimmed_tables_dir                                                                     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables
                     n_contrasts_compared                                                                                                                                                2
                         min_common_genes                                                                                                                                            13993
                min_spearman_logFC_common                                                                                                                                         0.998895
                       min_strict_jaccard                                                                                                                                          0.92451
 min_direction_agreement_pct_strict_union                                                                                                                                              100
                            batch_ok_h2w7                                                                                                                                             TRUE
                           candidate_file                           E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_file_selection_candidates.csv
                            selected_file                       E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_file_selection_selected_files.csv
                             summary_file              E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_trimmed_vs_untrimmed_DE_stability_summary.csv
                      interpretation_file E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_trimmed_vs_untrimmed_DE_stability_interpretation_flags.csv
                                 log_file                    E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/logs/stepH2w7_trimmed_vs_untrimmed_DE_stability_log.txt
                                   status                                                                                                                 Step H2w7 completed successfully

输出文件：
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_trimmed_vs_untrimmed_DE_stability_summary.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_trimmed_vs_untrimmed_DE_stability_interpretation_flags.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_file_selection_candidates.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_DE_file_selection_selected_files.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w7_trimmed_vs_untrimmed_DE_stability_overall_summary.csv 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/objects/stepH2w7_trimmed_vs_untrimmed_DE_stability_config.rds 
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/logs/stepH2w7_trimmed_vs_untrimmed_DE_stability_log.txt 

结果：H2w7 成功。已完成 trimmed 与 untrimmed 分支的 DE 稳定性比较。下一步可以进入 signature score / GSEA 稳定性比较。
2026-05-16 20:40:21 | ===== Step H2w7 finished ===== 

===== End of Step H2w7 =====
===== StepH2w9_current78 GSEA stability =====
Run time: 2026-05-16 20:44:30 CEST
Method: full Sus scrofa Hallmark fgsea; minSize=10, maxSize=500, eps=0; target pathways extracted afterward.

Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_trimmed_DE_candidates_t7.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_trimmed_DE_candidates_t28.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_selected_DE_files.csv
Selected DE files:
     branch                      contrast timepoint                                                                                                                         de_file
1 untrimmed  ACLT_untreated_t7_vs_Control        t7                                E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t7_vs_CON_QLF.csv
2 untrimmed ACLT_untreated_t28_vs_Control       t28                               E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t28_vs_CON_QLF.csv
3   trimmed  ACLT_untreated_t7_vs_Control        t7  E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t7_vs_Control_t0_QLF.csv
4   trimmed ACLT_untreated_t28_vs_Control       t28 E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t28_vs_Control_t0_QLF.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_gene_id_to_symbol_map_from_untrimmed_DE.csv
Hallmark gene sets loaded: 50
Unique Hallmark symbols: 4204
Gene ID -> symbol map rows from untrimmed DE: 14304
fgsea parameters: minSize=10, maxSize=500, eps=0, nproc=1

Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_ranked_gene_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_ranked_genes_untrimmed__t7.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_ranked_genes_untrimmed__t28.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_ranked_genes_trimmed__t7.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_ranked_genes_trimmed__t28.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_pathway_source_and_parameters.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_hallmark_coverage_by_branch.csv
Running fgsea for untrimmed__t7: ranked symbols=14229, pathways passing size filter=50
  |                                                                                                                                                                                                    |                                                                                                                                                                                            |   0%  |                                                                                                                                                                                                    |=======                                                                                                                                                                                     |   4%  |                                                                                                                                                                                                    |==============                                                                                                                                                                              |   8%  |                                                                                                                                                                                                    |======================                                                                                                                                                                      |  12%  |                                                                                                                                                                                                    |=============================                                                                                                                                                               |  15%  |                                                                                                                                                                                                    |====================================                                                                                                                                                        |  19%  |                                                                                                                                                                                                    |===========================================                                                                                                                                                 |  23%  |                                                                                                                                                                                                    |===================================================                                                                                                                                         |  27%  |                                                                                                                                                                                                    |==========================================================                                                                                                                                  |  31%  |                                                                                                                                                                                                    |=================================================================                                                                                                                           |  35%  |                                                                                                                                                                                                    |========================================================================                                                                                                                    |  38%  |                                                                                                                                                                                                    |================================================================================                                                                                                            |  42%  |                                                                                                                                                                                                    |=======================================================================================                                                                                                     |  46%  |                                                                                                                                                                                                    |==============================================================================================                                                                                              |  50%  |                                                                                                                                                                                                    |=====================================================================================================                                                                                       |  54%  |                                                                                                                                                                                                    |============================================================================================================                                                                                |  58%  |                                                                                                                                                                                                    |====================================================================================================================                                                                        |  62%  |                                                                                                                                                                                                    |===========================================================================================================================                                                                 |  65%  |                                                                                                                                                                                                    |==================================================================================================================================                                                          |  69%  |                                                                                                                                                                                                    |=========================================================================================================================================                                                   |  73%  |                                                                                                                                                                                                    |=================================================================================================================================================                                           |  77%  |                                                                                                                                                                                                    |========================================================================================================================================================                                    |  81%  |                                                                                                                                                                                                    |===============================================================================================================================================================                             |  85%  |                                                                                                                                                                                                    |======================================================================================================================================================================                      |  88%  |                                                                                                                                                                                                    |==============================================================================================================================================================================              |  92%  |                                                                                                                                                                                                    |=====================================================================================================================================================================================       |  96%  |                                                                                                                                                                                                    |============================================================================================================================================================================================| 100%

Running fgsea for untrimmed__t28: ranked symbols=13998, pathways passing size filter=50
  |                                                                                                                                                                                                    |                                                                                                                                                                                            |   0%  |                                                                                                                                                                                                    |=========                                                                                                                                                                                   |   5%  |                                                                                                                                                                                                    |===================                                                                                                                                                                         |  10%  |                                                                                                                                                                                                    |============================                                                                                                                                                                |  15%  |                                                                                                                                                                                                    |======================================                                                                                                                                                      |  20%  |                                                                                                                                                                                                    |===============================================                                                                                                                                             |  25%  |                                                                                                                                                                                                    |========================================================                                                                                                                                    |  30%  |                                                                                                                                                                                                    |==================================================================                                                                                                                          |  35%  |                                                                                                                                                                                                    |===========================================================================                                                                                                                 |  40%  |                                                                                                                                                                                                    |=====================================================================================                                                                                                       |  45%  |                                                                                                                                                                                                    |==============================================================================================                                                                                              |  50%  |                                                                                                                                                                                                    |=======================================================================================================                                                                                     |  55%  |                                                                                                                                                                                                    |=================================================================================================================                                                                           |  60%  |                                                                                                                                                                                                    |==========================================================================================================================                                                                  |  65%  |                                                                                                                                                                                                    |====================================================================================================================================                                                        |  70%  |                                                                                                                                                                                                    |=============================================================================================================================================                                               |  75%  |                                                                                                                                                                                                    |======================================================================================================================================================                                      |  80%  |                                                                                                                                                                                                    |================================================================================================================================================================                            |  85%  |                                                                                                                                                                                                    |=========================================================================================================================================================================                   |  90%  |                                                                                                                                                                                                    |===================================================================================================================================================================================         |  95%  |                                                                                                                                                                                                    |============================================================================================================================================================================================| 100%

Running fgsea for trimmed__t7: ranked symbols=14235, pathways passing size filter=50
  |                                                                                                                                                                                                    |                                                                                                                                                                                            |   0%  |                                                                                                                                                                                                    |=======                                                                                                                                                                                     |   4%  |                                                                                                                                                                                                    |==============                                                                                                                                                                              |   7%  |                                                                                                                                                                                                    |=====================                                                                                                                                                                       |  11%  |                                                                                                                                                                                                    |============================                                                                                                                                                                |  15%  |                                                                                                                                                                                                    |===================================                                                                                                                                                         |  19%  |                                                                                                                                                                                                    |==========================================                                                                                                                                                  |  22%  |                                                                                                                                                                                                    |=================================================                                                                                                                                           |  26%  |                                                                                                                                                                                                    |========================================================                                                                                                                                    |  30%  |                                                                                                                                                                                                    |===============================================================                                                                                                                             |  33%  |                                                                                                                                                                                                    |======================================================================                                                                                                                      |  37%  |                                                                                                                                                                                                    |=============================================================================                                                                                                               |  41%  |                                                                                                                                                                                                    |====================================================================================                                                                                                        |  44%  |                                                                                                                                                                                                    |===========================================================================================                                                                                                 |  48%  |                                                                                                                                                                                                    |=================================================================================================                                                                                           |  52%  |                                                                                                                                                                                                    |========================================================================================================                                                                                    |  56%  |                                                                                                                                                                                                    |===============================================================================================================                                                                             |  59%  |                                                                                                                                                                                                    |======================================================================================================================                                                                      |  63%  |                                                                                                                                                                                                    |=============================================================================================================================                                                               |  67%  |                                                                                                                                                                                                    |====================================================================================================================================                                                        |  70%  |                                                                                                                                                                                                    |===========================================================================================================================================                                                 |  74%  |                                                                                                                                                                                                    |==================================================================================================================================================                                          |  78%  |                                                                                                                                                                                                    |=========================================================================================================================================================                                   |  81%  |                                                                                                                                                                                                    |================================================================================================================================================================                            |  85%  |                                                                                                                                                                                                    |=======================================================================================================================================================================                     |  89%  |                                                                                                                                                                                                    |==============================================================================================================================================================================              |  93%  |                                                                                                                                                                                                    |=====================================================================================================================================================================================       |  96%  |                                                                                                                                                                                                    |============================================================================================================================================================================================| 100%

Running fgsea for trimmed__t28: ranked symbols=14056, pathways passing size filter=50
  |                                                                                                                                                                                                    |                                                                                                                                                                                            |   0%  |                                                                                                                                                                                                    |=========                                                                                                                                                                                   |   5%  |                                                                                                                                                                                                    |===================                                                                                                                                                                         |  10%  |                                                                                                                                                                                                    |============================                                                                                                                                                                |  15%  |                                                                                                                                                                                                    |======================================                                                                                                                                                      |  20%  |                                                                                                                                                                                                    |===============================================                                                                                                                                             |  25%  |                                                                                                                                                                                                    |========================================================                                                                                                                                    |  30%  |                                                                                                                                                                                                    |==================================================================                                                                                                                          |  35%  |                                                                                                                                                                                                    |===========================================================================                                                                                                                 |  40%  |                                                                                                                                                                                                    |=====================================================================================                                                                                                       |  45%  |                                                                                                                                                                                                    |==============================================================================================                                                                                              |  50%  |                                                                                                                                                                                                    |=======================================================================================================                                                                                     |  55%  |                                                                                                                                                                                                    |=================================================================================================================                                                                           |  60%  |                                                                                                                                                                                                    |==========================================================================================================================                                                                  |  65%  |                                                                                                                                                                                                    |====================================================================================================================================                                                        |  70%  |                                                                                                                                                                                                    |=============================================================================================================================================                                               |  75%  |                                                                                                                                                                                                    |======================================================================================================================================================                                      |  80%  |                                                                                                                                                                                                    |================================================================================================================================================================                            |  85%  |                                                                                                                                                                                                    |=========================================================================================================================================================================                   |  90%  |                                                                                                                                                                                                    |===================================================================================================================================================================================         |  95%  |                                                                                                                                                                                                    |============================================================================================================================================================================================| 100%

Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_fgsea_full_hallmark_all_branches.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_target_pathway_results_all_branches.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_GSEA_trimmed_vs_untrimmed_all_common_pathways.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_GSEA_target_pathways_stability.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_GSEA_stability_summary_by_timepoint.csv
null device 
          1 
null device 
          1 
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_versions_and_method_records.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_GSEA_stability_overall_summary.csv

===== StepH2w9_current78 overall summary =====
                                   metric
                     untrimmed_t7_DE_file
                    untrimmed_t28_DE_file
                       trimmed_t7_DE_file
                      trimmed_t28_DE_file
         n_full_hallmark_gene_sets_loaded
            n_fgsea_pathways_untrimmed_t7
           n_fgsea_pathways_untrimmed_t28
              n_fgsea_pathways_trimmed_t7
             n_fgsea_pathways_trimmed_t28
                            fgsea_minSize
                            fgsea_maxSize
                                fgsea_eps
                           rank_statistic
                 gene_identifier_for_GSEA
                 full_collection_FDR_used
                          target_pathways
     min_spearman_NES_all_common_pathways
        max_abs_delta_NES_target_pathways
 target_direction_agreement_all_contrasts
  target_pathways_all_found_all_contrasts
                   overall_stability_flag
                         selected_de_file
                          symbol_map_file
                      ranked_summary_file
                           fgsea_all_file
                          comparison_file
                   target_comparison_file
                summary_by_timepoint_file
                                 plot_png
                                 plot_pdf
                               output_dir
                                                                                                                                                                               value
                                                                                    E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t7_vs_CON_QLF.csv
                                                                                   E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t28_vs_CON_QLF.csv
                                                      E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t7_vs_Control_t0_QLF.csv
                                                     E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w6_trimmed_DE_ACLT_t28_vs_Control_t0_QLF.csv
                                                                                                                                                                                  50
                                                                                                                                                                                  50
                                                                                                                                                                                  50
                                                                                                                                                                                  50
                                                                                                                                                                                  50
                                                                                                                                                                                  10
                                                                                                                                                                                 500
                                                                                                                                                                                   0
                                                                                                                                                        sign(logFC) * -log10(PValue)
                                                                                                                                                                    pig gene symbols
                                                                                                                                                                                TRUE
                                                                                                          HALLMARK_INFLAMMATORY_RESPONSE; HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
                                                                                                                                                                             0.98982
                                                                                                                                                                            0.016429
                                                                                                                                                                                TRUE
                                                                                                                                                                                TRUE
                                                                                                                                                                       highly_stable
                             E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_selected_DE_files.csv
       E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_gene_id_to_symbol_map_from_untrimmed_DE.csv
                           E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_ranked_gene_summary.csv
              E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_fgsea_full_hallmark_all_branches.csv
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_GSEA_trimmed_vs_untrimmed_all_common_pathways.csv
                E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_GSEA_target_pathways_stability.csv
           E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/stepH2w9_current78_GSEA_stability_summary_by_timepoint.csv
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/figures/stepH2w9_current78_GSEA_NES_trimmed_vs_untrimmed_scatter.png
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/figures/stepH2w9_current78_GSEA_NES_trimmed_vs_untrimmed_scatter.pdf
                                                                      E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability

===== StepH2w9_current78 summary by timepoint =====
 timepoint                      contrast n_common_pathways spearman_NES_all_common_pathways pearson_NES_all_common_pathways max_abs_delta_NES_all_common_pathways
       t28 ACLT_untreated_t28_vs_Control                50                        0.9898199                       0.9818575                              2.023222
        t7  ACLT_untreated_t7_vs_Control                50                        0.9978872                       0.9702097                              1.816289
 mean_abs_delta_NES_all_common_pathways n_target_pathways max_abs_delta_NES_target_pathways target_direction_agreement target_pathways_all_found stability_flag
                             0.06169354                 2                        0.01642932                       TRUE                      TRUE  highly_stable
                             0.09026789                 2                        0.01197492                       TRUE                      TRUE  highly_stable

===== StepH2w9_current78 target pathway comparison =====
 timepoint                                    pathway NES_untrimmed NES_trimmed   delta_NES padj_untrimmed padj_trimmed direction_same
       t28 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION      1.513431   1.5298598  0.01642932    0.003483966  0.001588196           TRUE
       t28             HALLMARK_INFLAMMATORY_RESPONSE      1.011805   0.9980778 -0.01372747    0.527864313  0.597014925           TRUE
        t7 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION      1.632928   1.6209531 -0.01197492    0.001123435  0.001670631           TRUE
        t7             HALLMARK_INFLAMMATORY_RESPONSE      1.362350   1.3741244  0.01177454    0.030456364  0.040314643           TRUE

===== Version and method records =====
                     item                                                                                   value
                R_version                                                                                   4.5.2
    fgsea_package_version                                                                                  1.36.2
  msigdbr_package_version                                                                                  26.1.0
  ggplot2_package_version                                                                                   4.0.2
    dplyr_package_version                                                                                   1.1.4
          msigdbr_species                                                                              Sus scrofa
       msigdbr_collection                                                            MSigDB Hallmark collection H
            fgsea_minSize                                                                                      10
            fgsea_maxSize                                                                                     500
                fgsea_eps                                                                                       0
              fgsea_nproc                                                                                       1
           rank_statistic                                                            sign(logFC) * -log10(PValue)
 gene_identifier_for_GSEA                                     pig gene symbols, matching current Step20_current78
     target_pathways_note fgsea was run on the full Hallmark collection; target pathways were extracted afterward
              branch_role        fastp trimming sensitivity analysis, not replacement for main untrimmed analysis

StepH2w9_current78 completed.

StepH2w9_current78 GSEA stability completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/tables/stepH2w9_current78_GSEA_stability/logs/stepH2w9_current78_GSEA_stability_summary_to_send_me.txt

===== DIAGNOSE STEP20 SOURCE ISSUE =====
Run time: 2026-05-16 20:50:23 
Target script:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea_symbol_only_corrected_v2.R

1) File existence and metadata
file.exists: FALSE 
Target script does not exist.

2) Search for possible Step20 symbol-only scripts under 02_pig_early
                                                                                                                             file  size               mtime
  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/scripts/step20_current78_pig_early_hallmark_gsea_symbol_only_corrected_v2.R 35643 2026-05-16 20:28:42
     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/scripts/step20_current78_pig_early_hallmark_gsea_symbol_only_corrected.R 34186 2026-05-16 20:20:26
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/scripts/old Step20_current78_pig_early_hallmark_gsea_add_MSigDB_db_version.R 27910 2026-05-14 14:51:28

3) Check current output log timestamps
full_log exists: TRUE 
                                                                                                                                                                size               mtime
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/logs/step20_current78_pig_early_hallmark_gsea_full_log.txt 72567 2026-05-16 20:28:54
summary_log exists: TRUE 
                                                                                                                                                                          size               mtime
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/logs/step20_current78_pig_early_hallmark_gsea_summary_to_send_me.txt 18206 2026-05-14 14:52:46

4) Attempt to source target script with echo=TRUE
If the script is complete, you should now see Step20 progress messages.
If nothing appears or it returns immediately, the target script is likely not the full Step20 script.


===== DIAGNOSE STEP20 SOURCE ISSUE =====
Run time: 2026-05-16 20:55:50 
Target script:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea_symbol_only_corrected_v2.R

1) File existence and metadata
file.exists: FALSE 
Target script does not exist.

2) Search for possible Step20 symbol-only scripts under 02_pig_early
                                                                                                                             file  size               mtime
  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/scripts/step20_current78_pig_early_hallmark_gsea_symbol_only_corrected_v2.R 35643 2026-05-16 20:28:42
     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/scripts/step20_current78_pig_early_hallmark_gsea_symbol_only_corrected.R 34186 2026-05-16 20:20:26
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/scripts/old Step20_current78_pig_early_hallmark_gsea_add_MSigDB_db_version.R 27910 2026-05-14 14:51:28

3) Check current output log timestamps
full_log exists: TRUE 
                                                                                                                                                                size               mtime
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/logs/step20_current78_pig_early_hallmark_gsea_full_log.txt 74661 2026-05-16 20:28:54
summary_log exists: TRUE 
                                                                                                                                                                          size               mtime
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step20_current78_pig_early_hallmark_gsea/logs/step20_current78_pig_early_hallmark_gsea_summary_to_send_me.txt 18206 2026-05-14 14:52:46

4) Attempt to source target script with echo=TRUE
If the script is complete, you should now see Step20 progress messages.
If nothing appears or it returns immediately, the target script is likely not the full Step20 script.


>>> STEP20 SYMBOL-ONLY V3 SCRIPT STARTED SUCCESSFULLY <<<
Script run time: 2026-05-16 20:58:22
===== STEP20_CURRENT78 PIG EARLY HALLMARK GSEA =====
Run time: 2026-05-16 20:58:22 
Base directory: E:/R/ACLsenescence2/rebuild_submission/02_pig_early 
Method: full Sus scrofa Hallmark collection fgsea; target pathways extracted after full run.

Input DE tables:
t7 : E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t7_vs_CON_QLF.csv 
t28: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step17_pig_early_DE_t28_vs_CON_QLF.csv 

Gene ID to symbol mapping from available DE-table annotation:
n mapping rows: 28235 
n unambiguous gene_id symbols: 14304 

Building ranked gene list for t7...
