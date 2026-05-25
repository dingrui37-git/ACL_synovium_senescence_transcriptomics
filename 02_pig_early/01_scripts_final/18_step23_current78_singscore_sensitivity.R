============================================================
StepSINGSENS_01_pig_early_current78_singscore_sensitivity 
Started at: 2026-05-10 15:03:36 
Method version: 2026-05-10_FULL_LOGCPM_v2_robust_signature_column_detection_no_filtered_rank_background 
============================================================

Cleaning previous StepSINGSENS_01 outputs from output_root subfolders...
Number of previous StepSINGSENS_01 files removed: 2 

Required packages loaded.
edgeR version: 4.8.2 
limma version: 3.66.0 
singscore version: 1.30.0 
ggplot2 version: 4.0.2 
dplyr version: 1.1.4 

Input files:
Count matrix: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv 
Signature table: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv 
Main z-score scores: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv 
Main signature logCPM matrix: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_logCPM_matrix.csv 

Raw count table dimension: 35682 rows x 19 columns
First 10 columns:
 [1] "gene_id" "CON1"    "CON2"    "CON3"    "CON4"    "CON5"    "CON6"    "INJS1"   "INJS2"   "INJS3"  
Validated count matrix dimension: 35682 genes x 18 samples
Full genome-wide TMM logCPM matrix generated.
filterByExpr applied before rankGenes: FALSE
full_logCPM dimension: 35682 genes x 18 samples

Signature table dimension: 75 rows x 26 columns
Signature table column names:
 [1] "ortholog_name"                       "input"                               "input_number"                        "input_ensg"                         
 [5] "ensg_number"                         "ortholog_ensg"                       "description"                         "n_human_orthologs"                  
 [9] "n_mouse_inputs"                      "Entrez ID"                           "Gene symbol"                         "Gene name"                          
[13] "Cancer Cell"                         "Type of senescence"                  "Senescence Effect"                   "Reference"                          
[17] "human_gene"                          "mouse_symbol"                        "signature_direction"                 "pig_ensg"                           
[21] "pig_symbol"                          "n_pig_orthologs"                     "n_human_inputs"                      "mapped_to_pig_strict_1to1"          
[25] "present_in_step16_counts"            "included_in_current78_pig_signature"
Inferred signature gene ID column: pig_ensg 
Inferred signature direction column: signature_direction 

Current78 signature detection from full genome-wide TMM logCPM matrix:
Detected total: 75 
Detected up: 65 
Detected down: 10 

Main signature logCPM audit status: parsed 
Main signature logCPM inferred gene column:  
Main signature logCPM genes parsed: 0 

Running singscore::rankGenes on full genome-wide TMM logCPM matrix...
rankGenes completed.

singscore scores by sample:
   sample_id   group time_point sample_number directional_singscore up_singscore down_raw_singscore down_reoriented_singscore
1       CON1 Control         t0             1          -0.007644539    0.2294434         0.23708791               -0.23708791
2       CON2 Control         t0             2           0.052693846    0.2432635         0.19056963               -0.19056963
3       CON3 Control         t0             3           0.071125683    0.2421534         0.17102770               -0.17102770
4       CON4 Control         t0             4           0.017753908    0.2436734         0.22591949               -0.22591949
5       CON5 Control         t0             5           0.042727680    0.2666624         0.22393474               -0.22393474
6       CON6 Control         t0             6           0.060906308    0.2389479         0.17804160               -0.17804160
7      INJS1 ACLT_1W         t7             1           0.138206753    0.3069189         0.16871216               -0.16871216
8      INJS2 ACLT_1W         t7             2           0.114855517    0.3022910         0.18743552               -0.18743552
9      INJS3 ACLT_1W         t7             3           0.087697421    0.3079037         0.22020632               -0.22020632
10     INJS4 ACLT_1W         t7             4           0.113023066    0.3084817         0.19545862               -0.19545862
11     INJS5 ACLT_1W         t7             5           0.161524069    0.3052110         0.14368693               -0.14368693
12     INJS6 ACLT_1W         t7             6           0.102742722    0.2962222         0.19347948               -0.19347948
13     INJL1 ACLT_4W        t28             1           0.130921784    0.2676144         0.13669264               -0.13669264
14     INJL2 ACLT_4W        t28             2           0.102156856    0.2851912         0.18303431               -0.18303431
15     INJL3 ACLT_4W        t28             3           0.211321139    0.2787690         0.06744786               -0.06744786
16     INJL4 ACLT_4W        t28             4           0.209302387    0.2841286         0.07482619               -0.07482619
17     INJL5 ACLT_4W        t28             5           0.096210240    0.2804920         0.18428179               -0.18428179
18     INJL6 ACLT_4W        t28             6           0.083986387    0.2819736         0.19798722               -0.19798722

Main z-score score table columns:
[1] "sample_id"                   "group"                       "group_raw"                   "total_score_unoriented"      "up_score"                   
[6] "down_score_raw"              "down_score_reoriented"       "directional_score"           "directional_score_sum_means"
Inferred main score sample column: sample_id 
Inferred main directional score column: directional_score 
Inferred main up score column: up_score 
Inferred main down reoriented score column: down_score_reoriented 

Correlation summary against main z-score based scores:
                                                      comparison     main_score_column          singscore_column spearman_rho n_samples
1         main_zscore_directional_score_vs_directional_singscore     directional_score     directional_singscore    0.6532508        18
2                           main_zscore_up_score_vs_up_singscore              up_score              up_singscore    0.8988648        18
3 main_zscore_down_reoriented_score_vs_down_reoriented_singscore down_score_reoriented down_reoriented_singscore    0.7337461        18

Wilcoxon group comparison results:
                      score         comparison n_control n_case median_control median_case median_case_minus_control     p_value p_adj_BH_within_singscore_sensitivity
1     directional_singscore ACLT_1W vs Control         6      6     0.04771076   0.1139393                0.06622853 0.005074868                           0.007612302
2     directional_singscore ACLT_4W vs Control         6      6     0.04771076   0.1165393                0.06882856 0.005074868                           0.007612302
3              up_singscore ACLT_1W vs Control         6      6     0.24270843   0.3060650                0.06335652 0.005074868                           0.007612302
4              up_singscore ACLT_4W vs Control         6      6     0.24270843   0.2812328                0.03852439 0.005074868                           0.007612302
5 down_reoriented_singscore ACLT_1W vs Control         6      6    -0.20725219  -0.1904575                0.01679468 0.297953062                           0.297953062
6 down_reoriented_singscore ACLT_4W vs Control         6      6    -0.20725219  -0.1598635                0.04738871 0.092695803                           0.111234963

Final summary for review:
                                             metric                                                                                   value
1                                        run_status                                                                                 SUCCESS
2                                    method_version 2026-05-10_FULL_LOGCPM_v2_robust_signature_column_detection_no_filtered_rank_background
3                                   rankGenes_input         full genome-wide TMM-normalized logCPM matrix regenerated from raw count matrix
4             filterByExpr_applied_before_rankGenes                                                                                   FALSE
5                          signature_detected_total                                                                                      75
6                             signature_detected_up                                                                                      65
7                           signature_detected_down                                                                                      10
8                                           samples                                                                                      18
9                                            groups                                                         Control 6; ACLT_1W 6; ACLT_4W 6
10     directional_zscore_vs_singscore_spearman_rho                                                                       0.653250773993808
11              up_zscore_vs_singscore_spearman_rho                                                                       0.898864809081527
12 down_reoriented_zscore_vs_singscore_spearman_rho                                                                        0.73374613003096
13       directional_singscore_ACLT_1W_vs_Control_p                                                                     0.00507486809794025
14       directional_singscore_ACLT_4W_vs_Control_p                                                                     0.00507486809794025
15               main_signature_logCPM_audit_status                                                                                  parsed
16               main_signature_logCPM_genes_parsed                                                                                       0

Key output files:
1) Scores by sample:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/singscore_sensitivity/tables/StepSINGSENS_01_current78_singscore_scores_by_sample.csv 
2) Correlation with main z-score scores:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/singscore_sensitivity/tables/StepSINGSENS_01_current78_singscore_vs_main_zscore_correlation.csv 
3) Wilcoxon group comparisons:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/singscore_sensitivity/tables/StepSINGSENS_01_current78_singscore_group_comparison_wilcox.csv 
4) Signature detection audit:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/singscore_sensitivity/tables/StepSINGSENS_01_current78_signature_detection_audit.csv 
5) Full genome-wide TMM logCPM matrix used for rankGenes:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/singscore_sensitivity/source_data/StepSINGSENS_01_full_genomewide_TMM_logCPM_matrix_no_filterByExpr.csv 
6) Final summary:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/singscore_sensitivity/tables/StepSINGSENS_01_final_summary_for_review.csv 

Session information:
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.1.4      ggplot2_4.0.2    singscore_1.30.0 edgeR_4.8.2      limma_3.66.0    

loaded via a namespace (and not attached):
 [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 gtable_0.3.6                Biobase_2.70.0              lattice_0.22-7              vctrs_0.7.2                
 [7] tools_4.5.2                 generics_0.1.4              stats4_4.5.2                tibble_3.3.1                AnnotationDbi_1.72.0        RSQLite_2.4.5              
[13] blob_1.3.0                  pkgconfig_2.0.3             Matrix_1.7-4                RColorBrewer_1.1-3          S7_0.2.1                    S4Vectors_0.48.0           
[19] graph_1.88.1                lifecycle_1.0.5             stringr_1.6.0               compiler_4.5.2              farver_2.1.2                Biostrings_2.78.0          
[25] statmod_1.5.1               Seqinfo_1.0.0               pillar_1.11.1               crayon_1.5.3                tidyr_1.3.2                 cachem_1.1.0               
[31] DelayedArray_0.36.0         abind_1.4-8                 tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7               reshape2_1.4.5             
[37] purrr_1.2.1                 labeling_0.4.3              fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                   SparseArray_1.10.8         
[43] magrittr_2.0.4              S4Arrays_1.10.1             XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0             withr_3.0.2                
[49] scales_1.4.0                bit64_4.6.0-1               XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0           bit_4.6.0                  
[55] png_0.1-8                   memoise_2.0.1               GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                 Rcpp_1.1.1                 
[61] xtable_1.8-4                glue_1.8.0                  DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0             plyr_1.8.9                 
[67] R6_2.6.1                    MatrixGenerics_1.22.0      

============================================================
StepSINGSENS_01_pig_early_current78_singscore_sensitivity completed successfully.
Finished at: 2026-05-10 15:03:46 
============================================================
============================================================
StepFQC_01_pig_early_featureCounts_assignment_level_QC 
Started at: 2026-05-10 15:58:19 
Method version: 2026-05-10_featureCounts_assignment_QC_and_strandSpecific0_justification 
============================================================

Required packages loaded.
dplyr version: 1.1.4 
ggplot2 version: 4.0.2 
Input directory:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables

Input files found:
                                      label                                                                                                                path exists
                               bam_manifest                   E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_bam_manifest_used.csv   TRUE
                   featureCounts_annotation            E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_annotation.csv   TRUE
                   featureCounts_full_table            E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_full_table.csv   TRUE
                  featureCounts_run_summary           E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_run_summary.csv   TRUE
                   featureCounts_stat_table            E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_stat_table.csv   TRUE
            featureCounts_summary_by_sample     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_summary_by_sample.csv   TRUE
                          gene_count_matrix                   E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv   TRUE
              strandness_failures_stepBfix3              E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_failures.csv   TRUE
  strandness_featureCounts_detail_stepBfix3  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_featureCounts_detail.csv   TRUE
 strandness_featureCounts_summary_stepBfix3 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_featureCounts_summary.csv   TRUE
       strandness_overall_summary_stepBfix3       E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_overall_summary.csv   TRUE
      strandness_recommended_mode_stepBfix3      E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_recommended_mode.csv   TRUE

featureCounts stat table dimension: 14 rows x 19 columns
 [1] "category"  "CON1.bam"  "CON2.bam"  "CON3.bam"  "CON4.bam"  "CON5.bam"  "CON6.bam"  "INJS1.bam" "INJS2.bam" "INJS3.bam"

Gene count matrix dimension: 35682 rows x 19 columns
 [1] "gene_id" "CON1"    "CON2"    "CON3"    "CON4"    "CON5"    "CON6"    "INJS1"   "INJS2"   "INJS3"  

Expected 18 core samples:
 sample_id   group time_point sample_number
      CON1 Control         t0             1
      CON2 Control         t0             2
      CON3 Control         t0             3
      CON4 Control         t0             4
      CON5 Control         t0             5
      CON6 Control         t0             6
     INJS1 ACLT_1W         t7             1
     INJS2 ACLT_1W         t7             2
     INJS3 ACLT_1W         t7             3
     INJS4 ACLT_1W         t7             4
     INJS5 ACLT_1W         t7             5
     INJS6 ACLT_1W         t7             6
     INJL1 ACLT_4W        t28             1
     INJL2 ACLT_4W        t28             2
     INJL3 ACLT_4W        t28             3
     INJL4 ACLT_4W        t28             4
     INJL5 ACLT_4W        t28             5
     INJL6 ACLT_4W        t28             6

Group counts:

Control ACLT_1W ACLT_4W 
      6       6       6 

Parsed featureCounts stat long table preview:
 sample_id source_column                        status   count               status_standard
      CON1      CON1.bam                      Assigned 9802241                      Assigned
      CON1      CON1.bam          Unassigned_Ambiguity  298563          Unassigned_Ambiguity
      CON1      CON1.bam            Unassigned_Chimera       0            Unassigned_Chimera
      CON1      CON1.bam          Unassigned_Duplicate       0          Unassigned_Duplicate
      CON1      CON1.bam     Unassigned_FragmentLength       0     Unassigned_FragmentLength
      CON1      CON1.bam     Unassigned_MappingQuality       0     Unassigned_MappingQuality
      CON1      CON1.bam       Unassigned_MultiMapping       0       Unassigned_MultiMapping
      CON1      CON1.bam         Unassigned_NoFeatures 2264531         Unassigned_NoFeatures
      CON1      CON1.bam           Unassigned_NonSplit       0           Unassigned_NonSplit
      CON1      CON1.bam Unassigned_Overlapping_Length       0 Unassigned_Overlapping_Length
      CON1      CON1.bam          Unassigned_Read_Type       0          Unassigned_Read_Type
      CON1      CON1.bam          Unassigned_Secondary       0          Unassigned_Secondary
      CON1      CON1.bam          Unassigned_Singleton  232330          Unassigned_Singleton
      CON1      CON1.bam           Unassigned_Unmapped 2726316           Unassigned_Unmapped
      CON2      CON2.bam                      Assigned 7182805                      Assigned
      CON2      CON2.bam          Unassigned_Ambiguity  247018          Unassigned_Ambiguity
      CON2      CON2.bam            Unassigned_Chimera       0            Unassigned_Chimera
      CON2      CON2.bam          Unassigned_Duplicate       0          Unassigned_Duplicate
      CON2      CON2.bam     Unassigned_FragmentLength       0     Unassigned_FragmentLength
      CON2      CON2.bam     Unassigned_MappingQuality       0     Unassigned_MappingQuality

FeatureCounts assignment-level QC by sample preview:
 sample_id   group time_point sample_number total_fragments_considered assigned_fragments unassigned_fragments assignment_rate_pct parse_status Assigned Unassigned_Ambiguity
      CON1 Control         t0             1                   15323981            9802241              5521740              63.967         PASS  9802241               298563
      CON2 Control         t0             2                   11663160            7182805              4480355              61.585         PASS  7182805               247018
      CON3 Control         t0             3                   18546356           10664100              7882256              57.500         PASS 10664100               390603
      CON4 Control         t0             4                   19476871           13210214              6266657              67.825         PASS 13210214               436398
      CON5 Control         t0             5                   25325474           17038498              8286976              67.278         PASS 17038498               550013
      CON6 Control         t0             6                   13948648           10394313              3554335              74.518         PASS 10394313               315342
     INJS1 ACLT_1W         t7             1                   18282532           14050912              4231620              76.854         PASS 14050912               388232
     INJS2 ACLT_1W         t7             2                   17594194           13326409              4267785              75.743         PASS 13326409               383712
     INJS3 ACLT_1W         t7             3                   15479722           12059998              3419724              77.908         PASS 12059998               353362
     INJS4 ACLT_1W         t7             4                   11395236            8857217              2538019              77.727         PASS  8857217               262628
     INJS5 ACLT_1W         t7             5                   15945383           12116184              3829199              75.986         PASS 12116184               339148
     INJS6 ACLT_1W         t7             6                    9714098            7482593              2231505              77.028         PASS  7482593               206175
     INJL1 ACLT_4W        t28             1                   14071395            9737507              4333888              69.201         PASS  9737507               290820
     INJL2 ACLT_4W        t28             2                   12287052            9393742              2893310              76.452         PASS  9393742               243881
     INJL3 ACLT_4W        t28             3                   15410230           12121875              3288355              78.661         PASS 12121875               302662
     INJL4 ACLT_4W        t28             4                   17903438           13604350              4299088              75.987         PASS 13604350               329585
     INJL5 ACLT_4W        t28             5                   15491146           11646492              3844654              75.182         PASS 11646492               321235
     INJL6 ACLT_4W        t28             6                   12223978            9446634              2777344              77.280         PASS  9446634               253789
 Unassigned_Chimera Unassigned_Duplicate Unassigned_FragmentLength Unassigned_MappingQuality Unassigned_MultiMapping Unassigned_NoFeatures Unassigned_NonSplit
                  0                    0                         0                         0                       0               2264531                   0
                  0                    0                         0                         0                       0               1643653                   0
                  0                    0                         0                         0                       0               2657541                   0
                  0                    0                         0                         0                       0               2650994                   0
                  0                    0                         0                         0                       0               3172325                   0
                  0                    0                         0                         0                       0               1661647                   0
                  0                    0                         0                         0                       0               1686895                   0
                  0                    0                         0                         0                       0               1559967                   0
                  0                    0                         0                         0                       0               1752476                   0
                  0                    0                         0                         0                       0               1172299                   0
                  0                    0                         0                         0                       0               1567243                   0
                  0                    0                         0                         0                       0               1031086                   0
                  0                    0                         0                         0                       0               1633941                   0
                  0                    0                         0                         0                       0               1204251                   0
                  0                    0                         0                         0                       0               1552420                   0
                  0                    0                         0                         0                       0               1564309                   0
                  0                    0                         0                         0                       0               1582710                   0
                  0                    0                         0                         0                       0               1477961                   0
 Unassigned_Overlapping_Length Unassigned_Read_Type Unassigned_Secondary Unassigned_Singleton Unassigned_Unmapped
                             0                    0                    0               232330             2726316
                             0                    0                    0               204973             2384711
                             0                    0                    0               280317             4553795
                             0                    0                    0               266355             2912910
                             0                    0                    0               349781             4214857
                             0                    0                    0               274230             1303116
                             0                    0                    0               318222             1838271
                             0                    0                    0               243640             2080466
                             0                    0                    0               276866             1037020
                             0                    0                    0               250036              853056
                             0                    0                    0               285361             1637447
                             0                    0                    0               149413              844831
                             0                    0                    0               216604             2192523
                             0                    0                    0               250261             1194917
                             0                    0                    0               263729             1169544
                             0                    0                    0               224276             2180918
                             0                    0                    0               229120             1711589
                             0                    0                    0               246224              799370

Count matrix consistency audit preview:
 sample_id   group time_point sample_number assigned_from_count_matrix_colsum assigned_fragments total_fragments_considered assignment_rate_pct
      CON1 Control         t0             1                           9802241            9802241                   15323981              63.967
      CON2 Control         t0             2                           7182805            7182805                   11663160              61.585
      CON3 Control         t0             3                          10664100           10664100                   18546356              57.500
      CON4 Control         t0             4                          13210214           13210214                   19476871              67.825
      CON5 Control         t0             5                          17038498           17038498                   25325474              67.278
      CON6 Control         t0             6                          10394313           10394313                   13948648              74.518
     INJS1 ACLT_1W         t7             1                          14050912           14050912                   18282532              76.854
     INJS2 ACLT_1W         t7             2                          13326409           13326409                   17594194              75.743
     INJS3 ACLT_1W         t7             3                          12059998           12059998                   15479722              77.908
     INJS4 ACLT_1W         t7             4                           8857217            8857217                   11395236              77.727
     INJS5 ACLT_1W         t7             5                          12116184           12116184                   15945383              75.986
     INJS6 ACLT_1W         t7             6                           7482593            7482593                    9714098              77.028
     INJL1 ACLT_4W        t28             1                           9737507            9737507                   14071395              69.201
     INJL2 ACLT_4W        t28             2                           9393742            9393742                   12287052              76.452
     INJL3 ACLT_4W        t28             3                          12121875           12121875                   15410230              78.661
     INJL4 ACLT_4W        t28             4                          13604350           13604350                   17903438              75.987
     INJL5 ACLT_4W        t28             5                          11646492           11646492                   15491146              75.182
     INJL6 ACLT_4W        t28             6                           9446634            9446634                   12223978              77.280
 assigned_difference_count_matrix_minus_featureCounts_stat matches_assigned_fragments_exactly
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
============================================================
StepFQC_01_pig_early_featureCounts_assignment_level_QC 
Started at: 2026-05-10 16:08:31 
Method version: 2026-05-10_featureCounts_assignment_QC_and_strandSpecific0_justification_v2_vectorized_sample_id_fix 
============================================================

Required packages loaded.
dplyr version: 1.1.4 
ggplot2 version: 4.0.2 
Input directory:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables

Input files found:
                                      label                                                                                                                path exists
                               bam_manifest                   E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_bam_manifest_used.csv   TRUE
                   featureCounts_annotation            E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_annotation.csv   TRUE
                   featureCounts_full_table            E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_full_table.csv   TRUE
                  featureCounts_run_summary           E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_run_summary.csv   TRUE
                   featureCounts_stat_table            E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_stat_table.csv   TRUE
            featureCounts_summary_by_sample     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_featurecounts_summary_by_sample.csv   TRUE
                          gene_count_matrix                   E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv   TRUE
              strandness_failures_stepBfix3              E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_failures.csv   TRUE
  strandness_featureCounts_detail_stepBfix3  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_featureCounts_detail.csv   TRUE
 strandness_featureCounts_summary_stepBfix3 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_featureCounts_summary.csv   TRUE
       strandness_overall_summary_stepBfix3       E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_overall_summary.csv   TRUE
      strandness_recommended_mode_stepBfix3      E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepBfix3_pig_early_strandness_recommended_mode.csv   TRUE

featureCounts stat table dimension: 14 rows x 19 columns
 [1] "category"  "CON1.bam"  "CON2.bam"  "CON3.bam"  "CON4.bam"  "CON5.bam"  "CON6.bam"  "INJS1.bam" "INJS2.bam" "INJS3.bam"

Gene count matrix dimension: 35682 rows x 19 columns
 [1] "gene_id" "CON1"    "CON2"    "CON3"    "CON4"    "CON5"    "CON6"    "INJS1"   "INJS2"   "INJS3"  

Expected 18 core samples:
 sample_id   group time_point sample_number
      CON1 Control         t0             1
      CON2 Control         t0             2
      CON3 Control         t0             3
      CON4 Control         t0             4
      CON5 Control         t0             5
      CON6 Control         t0             6
     INJS1 ACLT_1W         t7             1
     INJS2 ACLT_1W         t7             2
     INJS3 ACLT_1W         t7             3
     INJS4 ACLT_1W         t7             4
     INJS5 ACLT_1W         t7             5
     INJS6 ACLT_1W         t7             6
     INJL1 ACLT_4W        t28             1
     INJL2 ACLT_4W        t28             2
     INJL3 ACLT_4W        t28             3
     INJL4 ACLT_4W        t28             4
     INJL5 ACLT_4W        t28             5
     INJL6 ACLT_4W        t28             6

Group counts:

Control ACLT_1W ACLT_4W 
      6       6       6 

Parsed featureCounts stat long table preview:
 sample_id source_column                        status   count               status_standard
      CON1      CON1.bam                      Assigned 9802241                      Assigned
      CON1      CON1.bam          Unassigned_Ambiguity  298563          Unassigned_Ambiguity
      CON1      CON1.bam            Unassigned_Chimera       0            Unassigned_Chimera
      CON1      CON1.bam          Unassigned_Duplicate       0          Unassigned_Duplicate
      CON1      CON1.bam     Unassigned_FragmentLength       0     Unassigned_FragmentLength
      CON1      CON1.bam     Unassigned_MappingQuality       0     Unassigned_MappingQuality
      CON1      CON1.bam       Unassigned_MultiMapping       0       Unassigned_MultiMapping
      CON1      CON1.bam         Unassigned_NoFeatures 2264531         Unassigned_NoFeatures
      CON1      CON1.bam           Unassigned_NonSplit       0           Unassigned_NonSplit
      CON1      CON1.bam Unassigned_Overlapping_Length       0 Unassigned_Overlapping_Length
      CON1      CON1.bam          Unassigned_Read_Type       0          Unassigned_Read_Type
      CON1      CON1.bam          Unassigned_Secondary       0          Unassigned_Secondary
      CON1      CON1.bam          Unassigned_Singleton  232330          Unassigned_Singleton
      CON1      CON1.bam           Unassigned_Unmapped 2726316           Unassigned_Unmapped
      CON2      CON2.bam                      Assigned 7182805                      Assigned
      CON2      CON2.bam          Unassigned_Ambiguity  247018          Unassigned_Ambiguity
      CON2      CON2.bam            Unassigned_Chimera       0            Unassigned_Chimera
      CON2      CON2.bam          Unassigned_Duplicate       0          Unassigned_Duplicate
      CON2      CON2.bam     Unassigned_FragmentLength       0     Unassigned_FragmentLength
      CON2      CON2.bam     Unassigned_MappingQuality       0     Unassigned_MappingQuality

FeatureCounts assignment-level QC by sample preview:
 sample_id   group time_point sample_number total_fragments_considered assigned_fragments unassigned_fragments assignment_rate_pct parse_status Assigned Unassigned_Ambiguity
      CON1 Control         t0             1                   15323981            9802241              5521740              63.967         PASS  9802241               298563
      CON2 Control         t0             2                   11663160            7182805              4480355              61.585         PASS  7182805               247018
      CON3 Control         t0             3                   18546356           10664100              7882256              57.500         PASS 10664100               390603
      CON4 Control         t0             4                   19476871           13210214              6266657              67.825         PASS 13210214               436398
      CON5 Control         t0             5                   25325474           17038498              8286976              67.278         PASS 17038498               550013
      CON6 Control         t0             6                   13948648           10394313              3554335              74.518         PASS 10394313               315342
     INJS1 ACLT_1W         t7             1                   18282532           14050912              4231620              76.854         PASS 14050912               388232
     INJS2 ACLT_1W         t7             2                   17594194           13326409              4267785              75.743         PASS 13326409               383712
     INJS3 ACLT_1W         t7             3                   15479722           12059998              3419724              77.908         PASS 12059998               353362
     INJS4 ACLT_1W         t7             4                   11395236            8857217              2538019              77.727         PASS  8857217               262628
     INJS5 ACLT_1W         t7             5                   15945383           12116184              3829199              75.986         PASS 12116184               339148
     INJS6 ACLT_1W         t7             6                    9714098            7482593              2231505              77.028         PASS  7482593               206175
     INJL1 ACLT_4W        t28             1                   14071395            9737507              4333888              69.201         PASS  9737507               290820
     INJL2 ACLT_4W        t28             2                   12287052            9393742              2893310              76.452         PASS  9393742               243881
     INJL3 ACLT_4W        t28             3                   15410230           12121875              3288355              78.661         PASS 12121875               302662
     INJL4 ACLT_4W        t28             4                   17903438           13604350              4299088              75.987         PASS 13604350               329585
     INJL5 ACLT_4W        t28             5                   15491146           11646492              3844654              75.182         PASS 11646492               321235
     INJL6 ACLT_4W        t28             6                   12223978            9446634              2777344              77.280         PASS  9446634               253789
 Unassigned_Chimera Unassigned_Duplicate Unassigned_FragmentLength Unassigned_MappingQuality Unassigned_MultiMapping Unassigned_NoFeatures Unassigned_NonSplit
                  0                    0                         0                         0                       0               2264531                   0
                  0                    0                         0                         0                       0               1643653                   0
                  0                    0                         0                         0                       0               2657541                   0
                  0                    0                         0                         0                       0               2650994                   0
                  0                    0                         0                         0                       0               3172325                   0
                  0                    0                         0                         0                       0               1661647                   0
                  0                    0                         0                         0                       0               1686895                   0
                  0                    0                         0                         0                       0               1559967                   0
                  0                    0                         0                         0                       0               1752476                   0
                  0                    0                         0                         0                       0               1172299                   0
                  0                    0                         0                         0                       0               1567243                   0
                  0                    0                         0                         0                       0               1031086                   0
                  0                    0                         0                         0                       0               1633941                   0
                  0                    0                         0                         0                       0               1204251                   0
                  0                    0                         0                         0                       0               1552420                   0
                  0                    0                         0                         0                       0               1564309                   0
                  0                    0                         0                         0                       0               1582710                   0
                  0                    0                         0                         0                       0               1477961                   0
 Unassigned_Overlapping_Length Unassigned_Read_Type Unassigned_Secondary Unassigned_Singleton Unassigned_Unmapped
                             0                    0                    0               232330             2726316
                             0                    0                    0               204973             2384711
                             0                    0                    0               280317             4553795
                             0                    0                    0               266355             2912910
                             0                    0                    0               349781             4214857
                             0                    0                    0               274230             1303116
                             0                    0                    0               318222             1838271
                             0                    0                    0               243640             2080466
                             0                    0                    0               276866             1037020
                             0                    0                    0               250036              853056
                             0                    0                    0               285361             1637447
                             0                    0                    0               149413              844831
                             0                    0                    0               216604             2192523
                             0                    0                    0               250261             1194917
                             0                    0                    0               263729             1169544
                             0                    0                    0               224276             2180918
                             0                    0                    0               229120             1711589
                             0                    0                    0               246224              799370

Count matrix consistency audit preview:
 sample_id   group time_point sample_number assigned_from_count_matrix_colsum assigned_fragments total_fragments_considered assignment_rate_pct
      CON1 Control         t0             1                           9802241            9802241                   15323981              63.967
      CON2 Control         t0             2                           7182805            7182805                   11663160              61.585
      CON3 Control         t0             3                          10664100           10664100                   18546356              57.500
      CON4 Control         t0             4                          13210214           13210214                   19476871              67.825
      CON5 Control         t0             5                          17038498           17038498                   25325474              67.278
      CON6 Control         t0             6                          10394313           10394313                   13948648              74.518
     INJS1 ACLT_1W         t7             1                          14050912           14050912                   18282532              76.854
     INJS2 ACLT_1W         t7             2                          13326409           13326409                   17594194              75.743
     INJS3 ACLT_1W         t7             3                          12059998           12059998                   15479722              77.908
     INJS4 ACLT_1W         t7             4                           8857217            8857217                   11395236              77.727
     INJS5 ACLT_1W         t7             5                          12116184           12116184                   15945383              75.986
     INJS6 ACLT_1W         t7             6                           7482593            7482593                    9714098              77.028
     INJL1 ACLT_4W        t28             1                           9737507            9737507                   14071395              69.201
     INJL2 ACLT_4W        t28             2                           9393742            9393742                   12287052              76.452
     INJL3 ACLT_4W        t28             3                          12121875           12121875                   15410230              78.661
     INJL4 ACLT_4W        t28             4                          13604350           13604350                   17903438              75.987
     INJL5 ACLT_4W        t28             5                          11646492           11646492                   15491146              75.182
     INJL6 ACLT_4W        t28             6                           9446634            9446634                   12223978              77.280
 assigned_difference_count_matrix_minus_featureCounts_stat matches_assigned_fragments_exactly
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE
                                                         0                               TRUE

Strand-specific mode audit:
Recommended mode raw extraction: NA 
stepBfix3 strandness failures rows: 54 

Strandness featureCounts summary table:
 strandSpecific n_samples_total n_success n_failed mean_assignment_rate_pct median_assignment_rate_pct min_assignment_rate_pct max_assignment_rate_pct
              0              18         0       18                       NA                         NA                      NA                      NA
              1              18         0       18                       NA                         NA                      NA                      NA
              2              18         0       18                       NA                         NA                      NA                      NA

Strandness overall summary table:
                            metric                                                                                                      value
                      project_root                                                                                        E:/R/ACLsenescence2
                 bam_manifest_file E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/stepB_pig_early_bam_manifest_for_strandness.csv
                          gtf_file                      E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz
                n_expected_samples                                                                                                         18
                        n_bam_rows                                                                                                         18
                      n_total_runs                                                                                                         54
                 n_successful_runs                                                                                                          0
                     n_failed_runs                                                                                                         54
 best_mode_by_mean_assignment_rate                                                                                                       <NA>
                            status                                                                                    Step B fix v3 completed

Strandness recommended mode table:
 strandSpecific n_samples_total n_success n_failed mean_assignment_rate_pct median_assignment_rate_pct min_assignment_rate_pct max_assignment_rate_pct recommended_strandSpecific
             NA              18         0       18                       NA                         NA                      NA                      NA                         NA
                                          recommendation_rule
 No successful featureCounts run; inspect failure table first

Group-level featureCounts assignment QC summary:
# A tibble: 3 × 12
  group   n_samples mean_total_fragments_con…¹ median_total_fragmen…² mean_assigned_fragme…³ median_assigned_frag…⁴ mean_unassigned_frag…⁵ median_unassigned_fr…⁶ mean_assignment_rate…⁷
  <fct>       <int>                      <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>
1 Control         6                  17380748.              16935168.              11382028.              10529206.               5998720.               5894198.                   65.4
2 ACLT_1W         6                  14735194.              15712552.              11315552.              12088091                3419642                3624462.                   76.9
3 ACLT_4W         6                  14564540.              14740812.              10991767.              10692000.               3572773.               3566504.                   75.5
# ℹ abbreviated names: ¹​mean_total_fragments_considered, ²​median_total_fragments_considered, ³​mean_assigned_fragments, ⁴​median_assigned_fragments, ⁵​mean_unassigned_fragments,
#   ⁶​median_unassigned_fragments, ⁷​mean_assignment_rate_pct
# ℹ 3 more variables: median_assignment_rate_pct <dbl>, min_assignment_rate_pct <dbl>, max_assignment_rate_pct <dbl>

Final summary for review:
                                                metric                                                                                                value
                                            run_status                                                                                              SUCCESS
                                        method_version 2026-05-10_featureCounts_assignment_QC_and_strandSpecific0_justification_v2_vectorized_sample_id_fix
                                             input_dir                                           E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables
                                           output_root                      E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC
                                 expected_core_samples                                                                                                   18
                     featureCounts_stat_samples_parsed                                                                                                   18
                    missing_featureCounts_stat_samples                                                                                                 none
                      extra_featureCounts_stat_samples                                                                                                 none
           all_core_featureCounts_stat_samples_present                                                                                                 TRUE
                              featureCounts_parse_PASS                                                                                                 TRUE
                       mean_total_fragments_considered                                                                                         15560160.778
                               mean_assigned_fragments                                                                                         11229782.444
                             mean_unassigned_fragments                                                                                          4330378.333
                              mean_assignment_rate_pct                                                                                               72.593
                               min_assignment_rate_pct                                                                                                 57.5
                               max_assignment_rate_pct                                                                                               78.661
                          count_matrix_samples_present                                                                                                   18
                          count_matrix_missing_samples                                                                                                 none
 count_matrix_assigned_counts_match_featureCounts_stat                                                                                                 TRUE
             strandness_stepBfix3_recommended_mode_raw                                                                                           not_parsed
                    strandness_stepBfix3_failures_rows                                                                                                   54

Key output files:
1) FeatureCounts assignment-level QC by sample:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_featureCounts_assignment_level_QC_by_sample.csv 
2) FeatureCounts status-wide table:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_featureCounts_assignment_status_wide.csv 
3) Count matrix vs featureCounts Assigned audit:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_count_matrix_vs_featureCounts_assigned_audit.csv 
4) Group-level featureCounts assignment QC summary:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_group_featureCounts_assignment_QC_summary.csv 
5) Strand-specific mode comparison combined table:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_stepBfix3_strandSpecific_mode_comparison_combined.csv 
6) Strand-specific recommended mode table:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_stepBfix3_strandSpecific_recommended_mode.csv 
7) Methods text:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_methods_text_featureCounts_assignment_QC.txt 
8) Final summary table:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/tables/StepFQC_01_final_summary_for_review.csv 
9) Assignment-rate figure:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/featurecounts_QC/figures/StepFQC_01_featureCounts_assignment_rate_by_sample.pdf 

Session information:
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.1.4      ggplot2_4.0.2    singscore_1.30.0 edgeR_4.8.2      limma_3.66.0    

loaded via a namespace (and not attached):
 [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 gtable_0.3.6                Biobase_2.70.0              lattice_0.22-7              vctrs_0.7.2                
 [7] tools_4.5.2                 generics_0.1.4              stats4_4.5.2                tibble_3.3.1                AnnotationDbi_1.72.0        RSQLite_2.4.5              
[13] blob_1.3.0                  pkgconfig_2.0.3             Matrix_1.7-4                RColorBrewer_1.1-3          S7_0.2.1                    S4Vectors_0.48.0           
[19] graph_1.88.1                lifecycle_1.0.5             stringr_1.6.0               compiler_4.5.2              farver_2.1.2                Biostrings_2.78.0          
[25] statmod_1.5.1               Seqinfo_1.0.0               pillar_1.11.1               crayon_1.5.3                tidyr_1.3.2                 cachem_1.1.0               
[31] DelayedArray_0.36.0         abind_1.4-8                 tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7               reshape2_1.4.5             
[37] purrr_1.2.1                 labeling_0.4.3              fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                   SparseArray_1.10.8         
[43] magrittr_2.0.4              S4Arrays_1.10.1             XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0             withr_3.0.2                
[49] scales_1.4.0                bit64_4.6.0-1               XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0           bit_4.6.0                  
[55] png_0.1-8                   memoise_2.0.1               GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                 Rcpp_1.1.1                 
[61] xtable_1.8-4                glue_1.8.0                  DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0             plyr_1.8.9                 
[67] R6_2.6.1                    MatrixGenerics_1.22.0      
============================================================
StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison
Started at: 2026-05-10 16:34:26 
Method version: 2026-05-10_clean_strandSpecific_0_1_2_featureCounts_comparison 
============================================================

Cleaning previous StepSTRAND_01 outputs from output_root subfolders...
Number of previous StepSTRAND_01 files removed: 0 

Required packages loaded.
Rsubread version: 2.24.0 
dplyr version: 1.1.4 
ggplot2 version: 4.0.2 

Script archive note: --file argument was not detected. If running interactively, manually save this script in:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/scripts

Input BAM directory:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align

GTF file candidate:
E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz

Output root:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison

Expected 18 core samples:
 sample_id   group time_point sample_number
      CON1 Control         t0             1
      CON2 Control         t0             2
      CON3 Control         t0             3
      CON4 Control         t0             4
      CON5 Control         t0             5
      CON6 Control         t0             6
     INJS1 ACLT_1W         t7             1
     INJS2 ACLT_1W         t7             2
     INJS3 ACLT_1W         t7             3
     INJS4 ACLT_1W         t7             4
     INJS5 ACLT_1W         t7             5
     INJS6 ACLT_1W         t7             6
     INJL1 ACLT_4W        t28             1
     INJL2 ACLT_4W        t28             2
     INJL3 ACLT_4W        t28             3
     INJL4 ACLT_4W        t28             4
     INJL5 ACLT_4W        t28             5
     INJL6 ACLT_4W        t28             6

Group counts:

Control ACLT_1W ACLT_4W 
      6       6       6 


ERROR:
GTF .gz file does not exist: E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz

============================================================
StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison finished with status: FAILED 
Finished at: 2026-05-10 16:34:27 
Summary log saved to:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/logs/StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison_summary_log.txt
============================================================

Caught error:
GTF .gz file does not exist: E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz 

============================================================
StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison finished with status: FAILED 
Finished at: 2026-05-10 16:34:27 
Summary log saved to:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/logs/StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison_summary_log.txt
============================================================
Cleaning previous StepSTRAND_01 outputs from output_root subfolders...
Number of previous StepSTRAND_01 files removed: 1 

============================================================
StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison
Started at: 2026-05-10 16:38:16 
Method version: 2026-05-10_clean_strandSpecific_0_1_2_featureCounts_comparison_v2_safe_logging 
============================================================

Required packages loaded.
Rsubread version: 2.24.0 
dplyr version: 1.1.4 
ggplot2 version: 4.0.2 

Script archive note: --file argument was not detected. If running interactively, manually save this script in:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/scripts

Input BAM directory:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align

GTF file candidate:
E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz

Output root:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison

Expected 18 core samples:
 sample_id   group time_point sample_number
      CON1 Control         t0             1
      CON2 Control         t0             2
      CON3 Control         t0             3
      CON4 Control         t0             4
      CON5 Control         t0             5
      CON6 Control         t0             6
     INJS1 ACLT_1W         t7             1
     INJS2 ACLT_1W         t7             2
     INJS3 ACLT_1W         t7             3
     INJS4 ACLT_1W         t7             4
     INJS5 ACLT_1W         t7             5
     INJS6 ACLT_1W         t7             6
     INJL1 ACLT_4W        t28             1
     INJL2 ACLT_4W        t28             2
     INJL3 ACLT_4W        t28             3
     INJL4 ACLT_4W        t28             4
     INJL5 ACLT_4W        t28             5
     INJL6 ACLT_4W        t28             6

Group counts:

Control ACLT_1W ACLT_4W 
      6       6       6 


ERROR:
GTF .gz file does not exist: E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz

============================================================
StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison finished with status: FAILED 
Finished at: 2026-05-10 16:38:17 
Summary log saved to:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/logs/StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison_summary_log.txt
============================================================
[1] FALSE
============================================================
StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison 
Started at: 2026-05-10 16:55:33 
Method version: 2026-05-10_clean_strandSpecific_0_1_2_featureCounts_comparison_v4_fixed_gtf_path 
============================================================

Required packages loaded.
Rsubread version: 2.24.0 
dplyr version: 1.1.4 
ggplot2 version: 4.0.2 
tidyr version: 1.3.2 

Script archive note: --file argument was not detected. If running interactively, manually save this script in:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/scripts 

Input BAM directory:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align

Locked count matrix file:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv

Output root:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison

Expected 18 core samples:
 sample_id   group time_point sample_number
      CON1 Control         t0             1
      CON2 Control         t0             2
      CON3 Control         t0             3
      CON4 Control         t0             4
      CON5 Control         t0             5
      CON6 Control         t0             6
     INJS1 ACLT_1W         t7             1
     INJS2 ACLT_1W         t7             2
     INJS3 ACLT_1W         t7             3
     INJS4 ACLT_1W         t7             4
     INJS5 ACLT_1W         t7             5
     INJS6 ACLT_1W         t7             6
     INJL1 ACLT_4W        t28             1
     INJL2 ACLT_4W        t28             2
     INJL3 ACLT_4W        t28             3
     INJL4 ACLT_4W        t28             4
     INJL5 ACLT_4W        t28             5
     INJL6 ACLT_4W        t28             6

Group counts:

Control ACLT_1W ACLT_4W 
      6       6       6 

BAM file check:
 sample_id                                                                                 bam_path bam_exists
      CON1  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/CON1.bam       TRUE
      CON2  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/CON2.bam       TRUE
      CON3  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/CON3.bam       TRUE
      CON4  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/CON4.bam       TRUE
      CON5  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/CON5.bam       TRUE
      CON6  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/CON6.bam       TRUE
     INJS1 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJS1.bam       TRUE
     INJS2 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJS2.bam       TRUE
     INJS3 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJS3.bam       TRUE
     INJS4 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJS4.bam       TRUE
     INJS5 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJS5.bam       TRUE
     INJS6 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJS6.bam       TRUE
     INJL1 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJL1.bam       TRUE
     INJL2 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJL2.bam       TRUE
     INJL3 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJL3.bam       TRUE
     INJL4 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJL4.bam       TRUE
     INJL5 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJL5.bam       TRUE
     INJL6 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/bam/step15J_rsubread_align/INJL6.bam       TRUE

Fixed GTF file supplied by user:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf

Selected fixed GTF candidate:
                                                                                                       candidate_path                       basename priority_score
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf Sus_scrofa.Sscrofa11.1.115.gtf            999
              selection_reason
 User-confirmed fixed GTF path

Working GTF used by featureCounts:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/source_data/gtf_working_copy/Sus_scrofa.Sscrofa11.1.115.gtf

Locked count matrix dimensions: 35682 rows x 19 columns
Locked count matrix sample column sums:
 sample_id locked_assigned_count
      CON1               9802241
      CON2               7182805
      CON3              10664100
      CON4              13210214
      CON5              17038498
      CON6              10394313
     INJS1              14050912
     INJS2              13326409
     INJS3              12059998
     INJS4               8857217
     INJS5              12116184
     INJS6               7482593
     INJL1               9737507
     INJL2               9393742
     INJL3              12121875
     INJL4              13604350
     INJL5              11646492
     INJL6               9446634


------------------------------------------------------------
Running Rsubread::featureCounts with strandSpecific = 0 
------------------------------------------------------------

        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       Rsubread 2.24.0

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 18 BAM files                                     ||
||                                                                            ||
||                           CON1.bam                                         ||
||                           CON2.bam                                         ||
||                           CON3.bam                                         ||
||                           CON4.bam                                         ||
||                           CON5.bam                                         ||
||                           CON6.bam                                         ||
||                           INJS1.bam                                        ||
||                           INJS2.bam                                        ||
||                           INJS3.bam                                        ||
||                           INJS4.bam                                        ||
||                           INJS5.bam                                        ||
||                           INJS6.bam                                        ||
||                           INJL1.bam                                        ||
||                           INJL2.bam                                        ||
||                           INJL3.bam                                        ||
||                           INJL4.bam                                        ||
||                           INJL5.bam                                        ||
||                           INJL6.bam                                        ||
||                                                                            ||
||              Paired-end : yes                                              ||
||        Count read pairs : yes                                              ||
||              Annotation : Sus_scrofa.Sscrofa11.1.115.gtf (GTF)             ||
||      Dir for temp files : .                                                ||
||                 Threads : 4                                                ||
||                   Level : meta-feature level                               ||
||      Multimapping reads : counted                                          ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file Sus_scrofa.Sscrofa11.1.115.gtf ...                    ||
||    Features : 567727                                                       ||
||    Meta-features : 35682                                                   ||
||    Chromosomes/contigs : 302                                               ||
||                                                                            ||
|| Process BAM file CON1.bam...                                               ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15323981                                             ||
||    Successfully assigned alignments : 9802241 (64.0%)                      ||
||    Running time : 0.32 minutes                                             ||
||                                                                            ||
|| Process BAM file CON2.bam...                                               ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 11663160                                             ||
||    Successfully assigned alignments : 7182805 (61.6%)                      ||
||    Running time : 0.25 minutes                                             ||
||                                                                            ||
|| Process BAM file CON3.bam...                                               ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 18546356                                             ||
||    Successfully assigned alignments : 10664100 (57.5%)                     ||
||    Running time : 0.37 minutes                                             ||
||                                                                            ||
|| Process BAM file CON4.bam...                                               ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 19476871                                             ||
||    Successfully assigned alignments : 13210214 (67.8%)                     ||
||    Running time : 0.41 minutes                                             ||
||                                                                            ||
|| Process BAM file CON5.bam...                                               ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 25325474                                             ||
||    Successfully assigned alignments : 17038498 (67.3%)                     ||
||    Running time : 0.55 minutes                                             ||
||                                                                            ||
|| Process BAM file CON6.bam...                                               ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 13948648                                             ||
||    Successfully assigned alignments : 10394313 (74.5%)                     ||
||    Running time : 0.33 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS1.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 18282532                                             ||
||    Successfully assigned alignments : 14050912 (76.9%)                     ||
||    Running time : 0.74 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS2.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 17594194                                             ||
||    Successfully assigned alignments : 13326409 (75.7%)                     ||
||    Running time : 0.57 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS3.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15479722                                             ||
||    Successfully assigned alignments : 12059998 (77.9%)                     ||
||    Running time : 0.51 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS4.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 11395236                                             ||
||    Successfully assigned alignments : 8857217 (77.7%)                      ||
||    Running time : 0.42 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS5.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15945383                                             ||
||    Successfully assigned alignments : 12116184 (76.0%)                     ||
||    Running time : 0.60 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS6.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 9714098                                              ||
||    Successfully assigned alignments : 7482593 (77.0%)                      ||
||    Running time : 0.36 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL1.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 14071395                                             ||
||    Successfully assigned alignments : 9737507 (69.2%)                      ||
||    Running time : 0.51 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL2.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 12287052                                             ||
||    Successfully assigned alignments : 9393742 (76.5%)                      ||
||    Running time : 0.47 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL3.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15410230                                             ||
||    Successfully assigned alignments : 12121875 (78.7%)                     ||
||    Running time : 0.56 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL4.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 17903438                                             ||
||    Successfully assigned alignments : 13604350 (76.0%)                     ||
||    Running time : 0.64 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL5.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15491146                                             ||
||    Successfully assigned alignments : 11646492 (75.2%)                     ||
||    Running time : 0.57 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL6.bam...                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 12223978                                             ||
||    Successfully assigned alignments : 9446634 (77.3%)                      ||
||    Running time : 0.52 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
\\============================================================================//

Completed strandSpecific = 0 
 sample_id   group assigned_fragments total_fragments_considered assignment_rate_pct parse_status
      CON1 Control            9802241                   15323981            63.96667         PASS
      CON2 Control            7182805                   11663160            61.58541         PASS
      CON3 Control           10664100                   18546356            57.49971         PASS
      CON4 Control           13210214                   19476871            67.82513         PASS
      CON5 Control           17038498                   25325474            67.27810         PASS
      CON6 Control           10394313                   13948648            74.51843         PASS
     INJS1 ACLT_1W           14050912                   18282532            76.85430         PASS
     INJS2 ACLT_1W           13326409                   17594194            75.74322         PASS
     INJS3 ACLT_1W           12059998                   15479722            77.90836         PASS
     INJS4 ACLT_1W            8857217                   11395236            77.72737         PASS
     INJS5 ACLT_1W           12116184                   15945383            75.98553         PASS
     INJS6 ACLT_1W            7482593                    9714098            77.02818         PASS
     INJL1 ACLT_4W            9737507                   14071395            69.20072         PASS
     INJL2 ACLT_4W            9393742                   12287052            76.45237         PASS
     INJL3 ACLT_4W           12121875                   15410230            78.66122         PASS
     INJL4 ACLT_4W           13604350                   17903438            75.98736         PASS
     INJL5 ACLT_4W           11646492                   15491146            75.18160         PASS
     INJL6 ACLT_4W            9446634                   12223978            77.27954         PASS

------------------------------------------------------------
Running Rsubread::featureCounts with strandSpecific = 1 
------------------------------------------------------------

        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       Rsubread 2.24.0

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 18 BAM files                                     ||
||                                                                            ||
||                           CON1.bam                                         ||
||                           CON2.bam                                         ||
||                           CON3.bam                                         ||
||                           CON4.bam                                         ||
||                           CON5.bam                                         ||
||                           CON6.bam                                         ||
||                           INJS1.bam                                        ||
||                           INJS2.bam                                        ||
||                           INJS3.bam                                        ||
||                           INJS4.bam                                        ||
||                           INJS5.bam                                        ||
||                           INJS6.bam                                        ||
||                           INJL1.bam                                        ||
||                           INJL2.bam                                        ||
||                           INJL3.bam                                        ||
||                           INJL4.bam                                        ||
||                           INJL5.bam                                        ||
||                           INJL6.bam                                        ||
||                                                                            ||
||              Paired-end : yes                                              ||
||        Count read pairs : yes                                              ||
||              Annotation : Sus_scrofa.Sscrofa11.1.115.gtf (GTF)             ||
||      Dir for temp files : .                                                ||
||                 Threads : 4                                                ||
||                   Level : meta-feature level                               ||
||      Multimapping reads : counted                                          ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file Sus_scrofa.Sscrofa11.1.115.gtf ...                    ||
||    Features : 567727                                                       ||
||    Meta-features : 35682                                                   ||
||    Chromosomes/contigs : 302                                               ||
||                                                                            ||
|| Process BAM file CON1.bam...                                               ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15323981                                             ||
||    Successfully assigned alignments : 5158304 (33.7%)                      ||
||    Running time : 0.55 minutes                                             ||
||                                                                            ||
|| Process BAM file CON2.bam...                                               ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 11663160                                             ||
||    Successfully assigned alignments : 3791753 (32.5%)                      ||
||    Running time : 0.45 minutes                                             ||
||                                                                            ||
|| Process BAM file CON3.bam...                                               ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 18546356                                             ||
||    Successfully assigned alignments : 5651411 (30.5%)                      ||
||    Running time : 0.64 minutes                                             ||
||                                                                            ||
|| Process BAM file CON4.bam...                                               ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 19476871                                             ||
||    Successfully assigned alignments : 7008182 (36.0%)                      ||
||    Running time : 0.70 minutes                                             ||
||                                                                            ||
|| Process BAM file CON5.bam...                                               ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 25325474                                             ||
||    Successfully assigned alignments : 9027626 (35.6%)                      ||
||    Running time : 0.88 minutes                                             ||
||                                                                            ||
|| Process BAM file CON6.bam...                                               ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 13948648                                             ||
||    Successfully assigned alignments : 5508401 (39.5%)                      ||
||    Running time : 0.34 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS1.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 18282532                                             ||
||    Successfully assigned alignments : 7432599 (40.7%)                      ||
||    Running time : 0.62 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS2.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 17594194                                             ||
||    Successfully assigned alignments : 7050950 (40.1%)                      ||
||    Running time : 0.59 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS3.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15479722                                             ||
||    Successfully assigned alignments : 6382251 (41.2%)                      ||
||    Running time : 0.52 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS4.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 11395236                                             ||
||    Successfully assigned alignments : 4695553 (41.2%)                      ||
||    Running time : 0.35 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS5.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15945383                                             ||
||    Successfully assigned alignments : 6408503 (40.2%)                      ||
||    Running time : 0.51 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS6.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 9714098                                              ||
||    Successfully assigned alignments : 3947788 (40.6%)                      ||
||    Running time : 0.31 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL1.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 14071395                                             ||
||    Successfully assigned alignments : 5140837 (36.5%)                      ||
||    Running time : 0.42 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL2.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 12287052                                             ||
||    Successfully assigned alignments : 4951766 (40.3%)                      ||
||    Running time : 0.40 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL3.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15410230                                             ||
||    Successfully assigned alignments : 6360706 (41.3%)                      ||
||    Running time : 0.48 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL4.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 17903438                                             ||
||    Successfully assigned alignments : 7138640 (39.9%)                      ||
||    Running time : 0.58 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL5.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15491146                                             ||
||    Successfully assigned alignments : 6149246 (39.7%)                      ||
||    Running time : 0.53 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL6.bam...                                              ||
||    Strand specific : stranded                                              ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 12223978                                             ||
||    Successfully assigned alignments : 4975048 (40.7%)                      ||
||    Running time : 0.38 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
\\============================================================================//

Completed strandSpecific = 1 
 sample_id   group assigned_fragments total_fragments_considered assignment_rate_pct parse_status
      CON1 Control            5158304                   15323981            33.66164         PASS
      CON2 Control            3791753                   11663160            32.51051         PASS
      CON3 Control            5651411                   18546356            30.47181         PASS
      CON4 Control            7008182                   19476871            35.98207         PASS
      CON5 Control            9027626                   25325474            35.64642         PASS
      CON6 Control            5508401                   13948648            39.49057         PASS
     INJS1 ACLT_1W            7432599                   18282532            40.65410         PASS
     INJS2 ACLT_1W            7050950                   17594194            40.07544         PASS
     INJS3 ACLT_1W            6382251                   15479722            41.22975         PASS
     INJS4 ACLT_1W            4695553                   11395236            41.20628         PASS
     INJS5 ACLT_1W            6408503                   15945383            40.19034         PASS
     INJS6 ACLT_1W            3947788                    9714098            40.63978         PASS
     INJL1 ACLT_4W            5140837                   14071395            36.53395         PASS
     INJL2 ACLT_4W            4951766                   12287052            40.30068         PASS
     INJL3 ACLT_4W            6360706                   15410230            41.27587         PASS
     INJL4 ACLT_4W            7138640                   17903438            39.87301         PASS
     INJL5 ACLT_4W            6149246                   15491146            39.69523         PASS
     INJL6 ACLT_4W            4975048                   12223978            40.69909         PASS

------------------------------------------------------------
Running Rsubread::featureCounts with strandSpecific = 2 
------------------------------------------------------------

        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       Rsubread 2.24.0

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 18 BAM files                                     ||
||                                                                            ||
||                           CON1.bam                                         ||
||                           CON2.bam                                         ||
||                           CON3.bam                                         ||
||                           CON4.bam                                         ||
||                           CON5.bam                                         ||
||                           CON6.bam                                         ||
||                           INJS1.bam                                        ||
||                           INJS2.bam                                        ||
||                           INJS3.bam                                        ||
||                           INJS4.bam                                        ||
||                           INJS5.bam                                        ||
||                           INJS6.bam                                        ||
||                           INJL1.bam                                        ||
||                           INJL2.bam                                        ||
||                           INJL3.bam                                        ||
||                           INJL4.bam                                        ||
||                           INJL5.bam                                        ||
||                           INJL6.bam                                        ||
||                                                                            ||
||              Paired-end : yes                                              ||
||        Count read pairs : yes                                              ||
||              Annotation : Sus_scrofa.Sscrofa11.1.115.gtf (GTF)             ||
||      Dir for temp files : .                                                ||
||                 Threads : 4                                                ||
||                   Level : meta-feature level                               ||
||      Multimapping reads : counted                                          ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file Sus_scrofa.Sscrofa11.1.115.gtf ...                    ||
||    Features : 567727                                                       ||
||    Meta-features : 35682                                                   ||
||    Chromosomes/contigs : 302                                               ||
||                                                                            ||
|| Process BAM file CON1.bam...                                               ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15323981                                             ||
||    Successfully assigned alignments : 5162430 (33.7%)                      ||
||    Running time : 0.46 minutes                                             ||
||                                                                            ||
|| Process BAM file CON2.bam...                                               ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 11663160                                             ||
||    Successfully assigned alignments : 3802380 (32.6%)                      ||
||    Running time : 0.44 minutes                                             ||
||                                                                            ||
|| Process BAM file CON3.bam...                                               ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 18546356                                             ||
||    Successfully assigned alignments : 5647070 (30.4%)                      ||
||    Running time : 0.64 minutes                                             ||
||                                                                            ||
|| Process BAM file CON4.bam...                                               ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 19476871                                             ||
||    Successfully assigned alignments : 7027846 (36.1%)                      ||
||    Running time : 0.68 minutes                                             ||
||                                                                            ||
|| Process BAM file CON5.bam...                                               ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 25325474                                             ||
||    Successfully assigned alignments : 9018916 (35.6%)                      ||
||    Running time : 0.96 minutes                                             ||
||                                                                            ||
|| Process BAM file CON6.bam...                                               ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 13948648                                             ||
||    Successfully assigned alignments : 5493857 (39.4%)                      ||
||    Running time : 0.45 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS1.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 18282532                                             ||
||    Successfully assigned alignments : 7356098 (40.2%)                      ||
||    Running time : 0.62 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS2.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 17594194                                             ||
||    Successfully assigned alignments : 7011135 (39.8%)                      ||
||    Running time : 0.55 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS3.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15479722                                             ||
||    Successfully assigned alignments : 6350816 (41.0%)                      ||
||    Running time : 0.56 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS4.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 11395236                                             ||
||    Successfully assigned alignments : 4669958 (41.0%)                      ||
||    Running time : 0.52 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS5.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15945383                                             ||
||    Successfully assigned alignments : 6363541 (39.9%)                      ||
||    Running time : 0.53 minutes                                             ||
||                                                                            ||
|| Process BAM file INJS6.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 9714098                                              ||
||    Successfully assigned alignments : 3936646 (40.5%)                      ||
||    Running time : 0.36 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL1.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 14071395                                             ||
||    Successfully assigned alignments : 5123729 (36.4%)                      ||
||    Running time : 0.48 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL2.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 12287052                                             ||
||    Successfully assigned alignments : 4916502 (40.0%)                      ||
||    Running time : 0.45 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL3.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15410230                                             ||
||    Successfully assigned alignments : 6342431 (41.2%)                      ||
||    Running time : 0.54 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL4.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 17903438                                             ||
||    Successfully assigned alignments : 7118316 (39.8%)                      ||
||    Running time : 0.61 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL5.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 15491146                                             ||
||    Successfully assigned alignments : 6130614 (39.6%)                      ||
||    Running time : 0.77 minutes                                             ||
||                                                                            ||
|| Process BAM file INJL6.bam...                                              ||
||    Strand specific : reversely stranded                                    ||
||    Paired-end reads are included.                                          ||
||    Total alignments : 12223978                                             ||
||    Successfully assigned alignments : 4975535 (40.7%)                      ||
||    Running time : 0.46 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
\\============================================================================//

Completed strandSpecific = 2 
 sample_id   group assigned_fragments total_fragments_considered assignment_rate_pct parse_status
      CON1 Control            5162430                   15323981            33.68857         PASS
      CON2 Control            3802380                   11663160            32.60163         PASS
      CON3 Control            5647070                   18546356            30.44841         PASS
      CON4 Control            7027846                   19476871            36.08303         PASS
      CON5 Control            9018916                   25325474            35.61203         PASS
      CON6 Control            5493857                   13948648            39.38630         PASS
     INJS1 ACLT_1W            7356098                   18282532            40.23566         PASS
     INJS2 ACLT_1W            7011135                   17594194            39.84914         PASS
     INJS3 ACLT_1W            6350816                   15479722            41.02668         PASS
     INJS4 ACLT_1W            4669958                   11395236            40.98167         PASS
     INJS5 ACLT_1W            6363541                   15945383            39.90836         PASS
     INJS6 ACLT_1W            3936646                    9714098            40.52508         PASS
     INJL1 ACLT_4W            5123729                   14071395            36.41237         PASS
     INJL2 ACLT_4W            4916502                   12287052            40.01368         PASS
     INJL3 ACLT_4W            6342431                   15410230            41.15728         PASS
     INJL4 ACLT_4W            7118316                   17903438            39.75949         PASS
     INJL5 ACLT_4W            6130614                   15491146            39.57495         PASS
     INJL6 ACLT_4W            4975535                   12223978            40.70308         PASS

Strand-specific mode summary:
# A tibble: 3 × 10
  strandSpecific n_samples_total n_parse_pass mean_assigned_fragments median_assigned_frag…¹ mean_total_fragments…² mean_assignment_rate…³ median_assignment_ra…⁴ min_assignment_rate_…⁵
           <dbl>           <int>        <int>                   <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>
1              0              18           18               11229782.              11155296               15560161.                   72.6                   75.9                   57.5
2              1              18           18                5932198                5900328.              15560161.                   38.3                   40.0                   30.5
3              2              18           18                5913768.               5888842               15560161.                   38.2                   39.8                   30.4
# ℹ abbreviated names: ¹​median_assigned_fragments, ²​mean_total_fragments_considered, ³​mean_assignment_rate_pct, ⁴​median_assignment_rate_pct, ⁵​min_assignment_rate_pct
# ℹ 1 more variable: max_assignment_rate_pct <dbl>

Recommended strandSpecific mode:
 recommended_strandSpecific                                                        recommendation_rule
                          0 Mode with the highest mean assignment_rate_pct across the 18 core samples.
                                                                                                         method_note
 All modes were run using the same BAM files, GTF annotation and featureCounts parameters except for strandSpecific.

strandSpecific = 0 rerun vs locked count matrix audit:
 sample_id   group time_point sample_number assigned_from_fc_counts_colsum locked_count_matrix_colsum difference_mode0_rerun_minus_locked mode0_rerun_matches_locked_count_matrix
      CON1 Control         t0             1                        9802241                    9802241                                   0                                    TRUE
      CON2 Control         t0             2                        7182805                    7182805                                   0                                    TRUE
      CON3 Control         t0             3                       10664100                   10664100                                   0                                    TRUE
      CON4 Control         t0             4                       13210214                   13210214                                   0                                    TRUE
      CON5 Control         t0             5                       17038498                   17038498                                   0                                    TRUE
      CON6 Control         t0             6                       10394313                   10394313                                   0                                    TRUE
     INJS1 ACLT_1W         t7             1                       14050912                   14050912                                   0                                    TRUE
     INJS2 ACLT_1W         t7             2                       13326409                   13326409                                   0                                    TRUE
     INJS3 ACLT_1W         t7             3                       12059998                   12059998                                   0                                    TRUE
     INJS4 ACLT_1W         t7             4                        8857217                    8857217                                   0                                    TRUE
     INJS5 ACLT_1W         t7             5                       12116184                   12116184                                   0                                    TRUE
     INJS6 ACLT_1W         t7             6                        7482593                    7482593                                   0                                    TRUE
     INJL1 ACLT_4W        t28             1                        9737507                    9737507                                   0                                    TRUE
     INJL2 ACLT_4W        t28             2                        9393742                    9393742                                   0                                    TRUE
     INJL3 ACLT_4W        t28             3                       12121875                   12121875                                   0                                    TRUE
     INJL4 ACLT_4W        t28             4                       13604350                   13604350                                   0                                    TRUE
     INJL5 ACLT_4W        t28             5                       11646492                   11646492                                   0                                    TRUE
     INJL6 ACLT_4W        t28             6                        9446634                    9446634                                   0                                    TRUE

Final summary for review:
                                                          metric
                                                      run_status
                                                  method_version
                                               bam_files_present
                                           gtf_original_selected
                                                gtf_working_used
                        strandSpecific0_mean_assignment_rate_pct
                        strandSpecific1_mean_assignment_rate_pct
                        strandSpecific2_mean_assignment_rate_pct
                                      recommended_strandSpecific
                                             interpretation_flag
 strandSpecific0_rerun_assigned_counts_match_locked_count_matrix
                                                     output_root
                                                                                                                                                value
                                                                                                                                              SUCCESS
                                                                     2026-05-10_clean_strandSpecific_0_1_2_featureCounts_comparison_v4_fixed_gtf_path
                                                                                                                                                   18
                                 E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf
 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/source_data/gtf_working_copy/Sus_scrofa.Sscrofa11.1.115.gtf
                                                                                                                                               72.594
                                                                                                                                               38.341
                                                                                                                                                38.22
                                                                                                                                                    0
                                                                                                    strandSpecific_0_has_highest_mean_assignment_rate
                                                                                                                                                 TRUE
                                                             E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison

Key output files:
1) Assignment rate by sample:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/tables/StepSTRAND_01_strandSpecific_assignment_rate_by_sample.csv
2) Mode-level summary:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/tables/StepSTRAND_01_strandSpecific_mode_summary.csv
3) Recommended mode:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/tables/StepSTRAND_01_recommended_strandSpecific_mode.csv
4) strandSpecific = 0 rerun vs locked count matrix audit:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/tables/StepSTRAND_01_strandSpecific0_rerun_vs_locked_count_matrix_audit.csv
5) GTF candidates detected:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/tables/StepSTRAND_01_gtf_candidates_detected.csv
6) Assignment-rate figure:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/figures/StepSTRAND_01_strandSpecific_assignment_rate_comparison.pdf
7) Final summary:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/tables/StepSTRAND_01_final_summary_for_review.csv

Session information:
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_1.3.2      Rsubread_2.24.0  dplyr_1.1.4      ggplot2_4.0.2    singscore_1.30.0 edgeR_4.8.2      limma_3.66.0    

loaded via a namespace (and not attached):
 [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 gtable_0.3.6                Biobase_2.70.0              lattice_0.22-7              vctrs_0.7.2                
 [7] tools_4.5.2                 generics_0.1.4              parallel_4.5.2              stats4_4.5.2                tibble_3.3.1                AnnotationDbi_1.72.0       
[13] RSQLite_2.4.5               blob_1.3.0                  pkgconfig_2.0.3             Matrix_1.7-4                RColorBrewer_1.1-3          S7_0.2.1                   
[19] S4Vectors_0.48.0            graph_1.88.1                lifecycle_1.0.5             stringr_1.6.0               compiler_4.5.2              farver_2.1.2               
[25] Biostrings_2.78.0           statmod_1.5.1               Seqinfo_1.0.0               pillar_1.11.1               crayon_1.5.3                cachem_1.1.0               
[31] DelayedArray_0.36.0         abind_1.4-8                 tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7               reshape2_1.4.5             
[37] purrr_1.2.1                 labeling_0.4.3              fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                   SparseArray_1.10.8         
[43] magrittr_2.0.4              S4Arrays_1.10.1             XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0             withr_3.0.2                
[49] scales_1.4.0                bit64_4.6.0-1               XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0           bit_4.6.0                  
[55] png_0.1-8                   memoise_2.0.1               GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                 Rcpp_1.1.1                 
[61] xtable_1.8-4                glue_1.8.0                  DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0             plyr_1.8.9                 
[67] R6_2.6.1                    MatrixGenerics_1.22.0      

============================================================
StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison finished with status: SUCCESS 
Finished at: 2026-05-10 17:25:10 
Summary log saved to:
E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplement/strandSpecific_comparison/logs/StepSTRAND_01_pig_early_strandSpecific_featureCounts_comparison_summary_log.txt
============================================================
============================================================
StepCHRONQC_01_chronic_pig_Salmon_tximport_QC
Started at: 2026-05-10 19:16:44 
Method version: 2026-05-10_chronic_pig_Salmon_tximport_quantification_audit_v1 
============================================================

Script archive note: --file argument was not detected. If running interactively, manually save this script in:
E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/scripts

Required packages loaded.
dplyr version: 1.1.4 
ggplot2 version: 4.0.2 

Input paths:
                      label                                                                                                                                                  path
1                 quant_dir                                                                              E:/R/ACLsenescence2/rebuild_submission/raw data/GSE228848_synovium_quant
2                tables_dir                                                                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables
3             metadata_file E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score/chronic_step2_current78_sample_metadata_used.csv
4 tximport_run_summary_file                                              E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_tximport_run_summary.csv
5                  gtf_file                                  E:/R/ACLsenescence2/rebuild_submission/02_pig_early_fastp_sensitivity/reference_cache/Sus_scrofa.Sscrofa11.1.115.gtf
  exists
1   TRUE
2   TRUE
3   TRUE
4   TRUE
5   TRUE

Quant-like files found: 96 
Tissue counts from filename:

cartilage  synovium 
       48        48 
Group counts from filename:

    ACLT_alone_52W        Control_52W Reconstruction_52W         Repair_52W 
                24                 24                 24                 24 

Synovium-candidate quant files: 48 

    ACLT_alone_52W        Control_52W Reconstruction_52W         Repair_52W 
                12                 12                 12                 12 
Auditing quant file 1 of 48 : GSM7140486_CON1_synovium_quant.sf.txt.gz 
Auditing quant file 2 of 48 : GSM7140487_ACLT1_synovium_quant.sf.txt.gz 
Auditing quant file 3 of 48 : GSM7140488_ACLT2_synovium_quant.sf.txt.gz 
Auditing quant file 4 of 48 : GSM7140489_ACLT3_synovium_quant.sf.txt.gz 
Auditing quant file 5 of 48 : GSM7140490_ACLT4_synovium_quant.sf.txt.gz 
Auditing quant file 6 of 48 : GSM7140491_CON2_synovium_quant.sf.txt.gz 
Auditing quant file 7 of 48 : GSM7140492_ACLT5_synovium_quant.sf.txt.gz 
Auditing quant file 8 of 48 : GSM7140493_ACLT6_synovium_quant.sf.txt.gz 
Auditing quant file 9 of 48 : GSM7140494_ACLT7_synovium_quant.sf.txt.gz 
Auditing quant file 10 of 48 : GSM7140495_CON3_synovium_quant.sf.txt.gz 
Auditing quant file 11 of 48 : GSM7140496_ACLT8_synovium_quant.sf.txt.gz 
Auditing quant file 12 of 48 : GSM7140497_ACLT9_synovium_quant.sf.txt.gz 
Auditing quant file 13 of 48 : GSM7140498_ACLT10_synovium_quant.sf.txt.gz 
Auditing quant file 14 of 48 : GSM7140499_ACLT11_synovium_quant.sf.txt.gz 
Auditing quant file 15 of 48 : GSM7140500_CON4_synovium_quant.sf.txt.gz 
Auditing quant file 16 of 48 : GSM7140501_ACLT12_synovium_quant.sf.txt.gz 
Auditing quant file 17 of 48 : GSM7140502_CON5_synovium_quant.sf.txt.gz 
Auditing quant file 18 of 48 : GSM7140503_RECON1_synovium_quant.sf.txt.gz 
Auditing quant file 19 of 48 : GSM7140504_RECON2_synovium_quant.sf.txt.gz 
Auditing quant file 20 of 48 : GSM7140505_RECON3_synovium_quant.sf.txt.gz 
Auditing quant file 21 of 48 : GSM7140506_RECON4_synovium_quant.sf.txt.gz 
Auditing quant file 22 of 48 : GSM7140507_CON6_synovium_quant.sf.txt.gz 
Auditing quant file 23 of 48 : GSM7140508_RECON5_synovium_quant.sf.txt.gz 
Auditing quant file 24 of 48 : GSM7140509_RECON6_synovium_quant.sf.txt.gz 
Auditing quant file 25 of 48 : GSM7140510_RECON7_synovium_quant.sf.txt.gz 
Auditing quant file 26 of 48 : GSM7140511_CON7_synovium_quant.sf.txt.gz 
Auditing quant file 27 of 48 : GSM7140512_RECON8_synovium_quant.sf.txt.gz 
Auditing quant file 28 of 48 : GSM7140513_RECON9_synovium_quant.sf.txt.gz 
Auditing quant file 29 of 48 : GSM7140514_RECON10_synovium_quant.sf.txt.gz 
Auditing quant file 30 of 48 : GSM7140515_RECON11_synovium_quant.sf.txt.gz 
Auditing quant file 31 of 48 : GSM7140516_CON8_synovium_quant.sf.txt.gz 
Auditing quant file 32 of 48 : GSM7140517_RECON12_synovium_quant.sf.txt.gz 
Auditing quant file 33 of 48 : GSM7140518_CON9_synovium_quant.sf.txt.gz 
Auditing quant file 34 of 48 : GSM7140519_REPAIR1_synovium_quant.sf.txt.gz 
Auditing quant file 35 of 48 : GSM7140520_REPAIR2_synovium_quant.sf.txt.gz 
Auditing quant file 36 of 48 : GSM7140521_REPAIR3_synovium_quant.sf.txt.gz 
Auditing quant file 37 of 48 : GSM7140522_REPAIR4_synovium_quant.sf.txt.gz 
Auditing quant file 38 of 48 : GSM7140523_CON10_synovium_quant.sf.txt.gz 
Auditing quant file 39 of 48 : GSM7140524_REPAIR5_synovium_quant.sf.txt.gz 
Auditing quant file 40 of 48 : GSM7140525_REPAIR6_synovium_quant.sf.txt.gz 
Auditing quant file 41 of 48 : GSM7140526_REPAIR7_synovium_quant.sf.txt.gz 
Auditing quant file 42 of 48 : GSM7140527_CON11_synovium_quant.sf.txt.gz 
Auditing quant file 43 of 48 : GSM7140528_REPAIR8_synovium_quant.sf.txt.gz 
Auditing quant file 44 of 48 : GSM7140529_REPAIR9_synovium_quant.sf.txt.gz 
Auditing quant file 45 of 48 : GSM7140530_REPAIR10_synovium_quant.sf.txt.gz 
Auditing quant file 46 of 48 : GSM7140531_REPAIR11_synovium_quant.sf.txt.gz 
Auditing quant file 47 of 48 : GSM7140532_CON12_synovium_quant.sf.txt.gz 
Auditing quant file 48 of 48 : GSM7140533_REPAIR12_synovium_quant.sf.txt.gz 

Quant.sf integrity summary:
                              metric    value
1              quant_files_all_found       96
2     synovium_candidate_quant_files       48
3              quant_files_read_PASS       48
4  quant_files_standard_columns_PASS       48
5                  min_n_transcripts    46295
6                  max_n_transcripts    46295
7                 mean_n_transcripts    46295
8                        min_tpm_sum  1000000
9                        max_tpm_sum  1000000
10                      mean_tpm_sum  1000000
11                  min_numreads_sum  7483049
12                  max_numreads_sum 19992966
13                 mean_numreads_sum 13321649

Detected tximport/gene-level files:
                      label                                                                                                                   path exists
1               counts_file E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_estimated_counts_matrix.csv   TRUE
2            abundance_file    E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_abundance_tpm_matrix.csv   TRUE
3     effective_length_file           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_length_matrix.csv   TRUE
4      gene_annotation_file           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_annotation_from_gtf.csv   TRUE
5              tx2gene_file                            E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_tx2gene.csv   TRUE
6 tximport_run_summary_file               E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_tximport_run_summary.csv   TRUE

Gene-level matrix dimensions:
            matrix                                                                                                                   path n_genes_rows n_columns_total n_sample_columns
1           counts E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_estimated_counts_matrix.csv        22438              97               96
2    abundance_TPM    E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_abundance_tpm_matrix.csv        22438              97               96
3 effective_length           E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_length_matrix.csv        22438              97               96
  first_column
1      gene_id
2      gene_id
3      gene_id

Matrix sample-column consistency:
                  comparison identical_sample_columns counts_n other_n
1        counts_vs_abundance                     TRUE       96      96
2 counts_vs_effective_length                     TRUE       96      96

Gene-level matrix sample QC preview:
                     sample_id geo_accession label_from_sample_id total_estimated_counts detected_genes_count_gt0 detected_genes_count_ge1 tissue_from_sample_id
1    GSM7140438_CON1_cartilage    GSM7140438       CON1_cartilage               14591099                    16445                    16443             cartilage
2   GSM7140439_ACLT1_cartilage    GSM7140439      ACLT1_cartilage               15096743                    17006                    17000             cartilage
3   GSM7140440_ACLT2_cartilage    GSM7140440      ACLT2_cartilage               16784162                    16406                    16403             cartilage
4   GSM7140441_ACLT3_cartilage    GSM7140441      ACLT3_cartilage               21257502                    16094                    16092             cartilage
5   GSM7140442_ACLT4_cartilage    GSM7140442      ACLT4_cartilage               11437664                    16523                    16517             cartilage
6    GSM7140443_CON2_cartilage    GSM7140443       CON2_cartilage               13261603                    15686                    15685             cartilage
7   GSM7140444_ACLT5_cartilage    GSM7140444      ACLT5_cartilage               12477842                    16035                    16033             cartilage
8   GSM7140445_ACLT6_cartilage    GSM7140445      ACLT6_cartilage               14792199                    16188                    16186             cartilage
9   GSM7140446_ACLT7_cartilage    GSM7140446      ACLT7_cartilage               15318184                    16155                    16150             cartilage
10   GSM7140447_CON3_cartilage    GSM7140447       CON3_cartilage               16951157                    16604                    16600             cartilage
11  GSM7140448_ACLT8_cartilage    GSM7140448      ACLT8_cartilage               15483878                    16303                    16297             cartilage
12  GSM7140449_ACLT9_cartilage    GSM7140449      ACLT9_cartilage               16007927                    17038                    17035             cartilage
13 GSM7140450_ACLT10_cartilage    GSM7140450     ACLT10_cartilage               15004740                    16171                    16166             cartilage
14 GSM7140451_ACLT11_cartilage    GSM7140451     ACLT11_cartilage               16617065                    16502                    16499             cartilage
15   GSM7140452_CON4_cartilage    GSM7140452       CON4_cartilage               16618734                    16438                    16437             cartilage
16 GSM7140453_ACLT12_cartilage    GSM7140453     ACLT12_cartilage               14097209                    16048                    16044             cartilage
17   GSM7140454_CON5_cartilage    GSM7140454       CON5_cartilage               16031334                    16661                    16654             cartilage
18 GSM7140455_RECON1_cartilage    GSM7140455     RECON1_cartilage               13751988                    16436                    16432             cartilage
19 GSM7140456_RECON2_cartilage    GSM7140456     RECON2_cartilage                8845046                    16255                    16251             cartilage
20 GSM7140457_RECON3_cartilage    GSM7140457     RECON3_cartilage               19699395                    16461                    16456             cartilage
   inferred_group_from_sample_id tpm_sum_gene_level genes_with_tpm_gt0 effective_length_na_n effective_length_nonpositive_n effective_length_median
1                    Control_52W           993345.1              16445                     0                              0                1999.090
2                 ACLT_alone_52W           991537.3              17006                     0                              0                2045.647
3                 ACLT_alone_52W           992778.7              16406                     0                              0                2017.151
4                 ACLT_alone_52W           992511.1              16094                     0                              0                2032.960
5                 ACLT_alone_52W           991249.4              16523                     0                              0                1974.590
6                    Control_52W           993701.5              15686                     0                              0                1988.202
7                 ACLT_alone_52W           992246.1              16035                     0                              0                2061.716
8                 ACLT_alone_52W           993117.0              16188                     0                              0                1907.822
9                 ACLT_alone_52W           993973.5              16155                     0                              0                1919.220
10                   Control_52W           992636.2              16604                     0                              0                2040.076
11                ACLT_alone_52W           993768.0              16303                     0                              0                2013.006
12                ACLT_alone_52W           991883.1              17038                     0                              0                2092.524
13                ACLT_alone_52W           993311.1              16171                     0                              0                2048.118
14                ACLT_alone_52W           992865.4              16502                     0                              0                1970.960
15                   Control_52W           992464.1              16438                     0                              0                2060.209
16                ACLT_alone_52W           992709.7              16048                     0                              0                1805.490
17                   Control_52W           990791.4              16661                     0                              0                1992.369
18            Reconstruction_52W           992738.1              16436                     0                              0                1972.800
19            Reconstruction_52W           991763.2              16255                     0                              0                1931.476
20            Reconstruction_52W           991486.9              16461                     0                              0                2002.846

Group counts inferred from matrix sample IDs:
  inferred_group_from_sample_id n_samples
1                ACLT_alone_52W        24
2                   Control_52W        24
3            Reconstruction_52W        24
4                    Repair_52W        24

Matrix samples vs quant files audit preview:
   matrix_geo_accession            matrix_sample_id quant_file basename label_from_filename tissue_from_filename inferred_group_from_filename has_matching_synovium_quant_file
1            GSM7140438   GSM7140438_CON1_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
2            GSM7140439  GSM7140439_ACLT1_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
3            GSM7140440  GSM7140440_ACLT2_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
4            GSM7140441  GSM7140441_ACLT3_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
5            GSM7140442  GSM7140442_ACLT4_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
6            GSM7140443   GSM7140443_CON2_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
7            GSM7140444  GSM7140444_ACLT5_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
8            GSM7140445  GSM7140445_ACLT6_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
9            GSM7140446  GSM7140446_ACLT7_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
10           GSM7140447   GSM7140447_CON3_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
11           GSM7140448  GSM7140448_ACLT8_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
12           GSM7140449  GSM7140449_ACLT9_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
13           GSM7140450 GSM7140450_ACLT10_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
14           GSM7140451 GSM7140451_ACLT11_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
15           GSM7140452   GSM7140452_CON4_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
16           GSM7140453 GSM7140453_ACLT12_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
17           GSM7140454   GSM7140454_CON5_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
18           GSM7140455 GSM7140455_RECON1_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
19           GSM7140456 GSM7140456_RECON2_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE
20           GSM7140457 GSM7140457_RECON3_cartilage       <NA>     <NA>                <NA>                 <NA>                         <NA>                            FALSE

Metadata file loaded. Columns:
[1] "sample_id" "group_raw" "group"    
                                                                                                                                          metadata_file exists n_rows n_cols
1 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score/chronic_step2_current78_sample_metadata_used.csv   TRUE     24      3
  inferred_sample_col inferred_group_col
1           sample_id              group

Metadata vs matrix sample audit summary:
Metadata rows: 24 
Matched by exact matrix column: 24 
Matched by GEO accession: 24 

tx2gene transcript matching audit:
                                                                                 tx2gene_file tx2gene_exists
1 E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_tx2gene.csv           TRUE
                                                                                                    quant_file_used quant_n_transcripts tx2gene_n_transcripts matched_with_version
1 E:/R/ACLsenescence2/rebuild_submission/raw data/GSE228848_synovium_quant/GSM7140486_CON1_synovium_quant.sf.txt.gz               46295                 60440                    0
  matched_after_stripping_version match_rate_with_version_pct match_rate_after_stripping_version_pct
1                           45845                           0                                 99.028

tximport run summary preview:
                   metric                                                                                 value
1            project_root                                                                   E:/R/ACLsenescence2
2               quant_dir                                 E:/R/ACLsenescence2/data_raw/GSE228848_synovium_quant
3                gtf_file E:/R/ACLsenescence2/reference/Sus_scrofa_Ensembl115/Sus_scrofa.Sscrofa11.1.115.gtf.gz
4           tximport_type                                                                                salmon
5         ignoreTxVersion                                                                                  TRUE
6     countsFromAbundance                                                                                    no
7           n_quant_files                                                                                    96
8        n_unique_samples                                                                                    96
9          n_tx2gene_rows                                                                                 60440
10 n_gene_annotation_rows                                                                                 35682
11         counts_n_genes                                                                                 22438
12       counts_n_samples                                                                                    96
13      abundance_n_genes                                                                                 22438
14    abundance_n_samples                                                                                    96
15         length_n_genes                                                                                 22438
16       length_n_samples                                                                                    96

Final summary for review:
                                                metric                                                                               value
1                                           run_status                                                                             SUCCESS
2                                       method_version                      2026-05-10_chronic_pig_Salmon_tximport_quantification_audit_v1
3                                quant_files_all_found                                                                                  96
4                       synovium_candidate_quant_files                                                                                  48
5                                quant_files_read_PASS                                                                                  48
6                    quant_files_standard_columns_PASS                                                                                  48
7                                  counts_matrix_genes                                                                               22438
8                                counts_matrix_samples                                                                                  96
9                            abundance_matrix_detected                                                                                TRUE
10                    effective_length_matrix_detected                                                                                TRUE
11        counts_vs_abundance_sample_columns_identical                                                                                TRUE
12 counts_vs_effective_length_sample_columns_identical                                                                                TRUE
13                         mean_total_estimated_counts                                                                        13777276.977
14                          min_total_estimated_counts                                                                         7108677.116
15                          max_total_estimated_counts                                                                        21575315.985
16                       mean_detected_genes_count_gt0                                                                           16476.521
17                        min_detected_genes_count_gt0                                                                               15555
18                        max_detected_genes_count_gt0                                                                               17390
19                             mean_gene_level_TPM_sum                                                                          986368.266
20                              min_gene_level_TPM_sum                                                                          973187.294
21                              max_gene_level_TPM_sum                                                                          994784.807
22                                  metadata_file_rows                                                                                  24
23      metadata_vs_matrix_matched_by_exact_sample_col                                                                                  24
24         metadata_vs_matrix_matched_by_geo_accession                                                                                  24
25      tx2gene_match_rate_after_stripping_version_pct                                                                              99.028
26                         tximport_run_summary_exists                                                                                TRUE
27                                         output_root E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC

Key output files:
1) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_quant_file_manifest_all.csv
2) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_quant_file_manifest_synovium_candidate.csv
3) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_quant_sf_integrity_QC_by_file.csv
4) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_quant_sf_integrity_summary.csv
5) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_gene_level_matrix_dimensions.csv
6) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_matrix_sample_column_consistency.csv
7) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_gene_level_matrix_sample_QC.csv
8) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_matrix_samples_vs_quant_files_audit.csv
9) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_metadata_vs_matrix_samples_audit.csv
10) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_tx2gene_transcript_matching_audit.csv
11) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_methods_text_Salmon_tximport_QC.txt
12) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/Salmon_tximport_QC/tables/StepCHRONQC_01_final_summary_for_review.csv

Session information:
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_1.3.2      Rsubread_2.24.0  dplyr_1.1.4      ggplot2_4.0.2    singscore_1.30.0 edgeR_4.8.2      limma_3.66.0    

loaded via a namespace (and not attached):
 [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 gtable_0.3.6                Biobase_2.70.0              lattice_0.22-7              vctrs_0.7.2                
 [7] tools_4.5.2                 generics_0.1.4              parallel_4.5.2              stats4_4.5.2                tibble_3.3.1                AnnotationDbi_1.72.0       
[13] RSQLite_2.4.5               blob_1.3.0                  pkgconfig_2.0.3             Matrix_1.7-4                RColorBrewer_1.1-3          S7_0.2.1                   
[19] S4Vectors_0.48.0            graph_1.88.1                lifecycle_1.0.5             stringr_1.6.0               compiler_4.5.2              farver_2.1.2               
[25] Biostrings_2.78.0           statmod_1.5.1               Seqinfo_1.0.0               pillar_1.11.1               crayon_1.5.3                cachem_1.1.0               
[31] DelayedArray_0.36.0         abind_1.4-8                 tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7               reshape2_1.4.5             
[37] purrr_1.2.1                 labeling_0.4.3              fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                   SparseArray_1.10.8         
[43] magrittr_2.0.4              S4Arrays_1.10.1             XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0             withr_3.0.2                
[49] scales_1.4.0                bit64_4.6.0-1               XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0           bit_4.6.0                  
[55] png_0.1-8                   memoise_2.0.1               GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                 Rcpp_1.1.1                 
[61] xtable_1.8-4                glue_1.8.0                  DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0             plyr_1.8.9                 
[67] R6_2.6.1                    MatrixGenerics_1.22.0      
============================================================
StepCHRONLOCK_01_lock_synovium_only_48_matrices
Started at: 2026-05-10 19:23:33 
Method version: 2026-05-10_lock_chronic_pig_synovium_only_48_matrices_from_step24_tximport_outputs 
============================================================

Required packages loaded.
dplyr version: 1.1.4 
ggplot2 version: 4.0.2 

Input paths:
                  label                                                                                                                                                  path exists
1           counts_file                                E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_estimated_counts_matrix.csv   TRUE
2        abundance_file                                   E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_abundance_tpm_matrix.csv   TRUE
3 effective_length_file                                          E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_gene_level_length_matrix.csv   TRUE
4 tximport_summary_file                                              E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_pig_chronic_tximport_run_summary.csv   TRUE
5      metadata_24_file E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/chronic_step2_current78_signature_score/chronic_step2_current78_sample_metadata_used.csv   TRUE
6             quant_dir                                                                              E:/R/ACLsenescence2/rebuild_submission/raw data/GSE228848_synovium_quant   TRUE

All-sample matrix dimensions:
            matrix n_gene_rows n_columns_total n_sample_columns first_column
1 estimated_counts       22438              97               96      gene_id
2    abundance_TPM       22438              97               96      gene_id
3 effective_length       22438              97               96      gene_id

All-sample metadata group/tissue counts inferred from matrix columns:
           
            ACLT_alone_52W Control_52W Reconstruction_52W Repair_52W
  cartilage             12          12                 12         12
  synovium              12          12                 12         12

Synovium-only sample metadata:
                      sample_id geo_accession sample_label   tissue              group sample_number
1      GSM7140486_CON1_synovium    GSM7140486         CON1 synovium        Control_52W             1
2      GSM7140491_CON2_synovium    GSM7140491         CON2 synovium        Control_52W             2
3      GSM7140495_CON3_synovium    GSM7140495         CON3 synovium        Control_52W             3
4      GSM7140500_CON4_synovium    GSM7140500         CON4 synovium        Control_52W             4
5      GSM7140502_CON5_synovium    GSM7140502         CON5 synovium        Control_52W             5
6      GSM7140507_CON6_synovium    GSM7140507         CON6 synovium        Control_52W             6
7      GSM7140511_CON7_synovium    GSM7140511         CON7 synovium        Control_52W             7
8      GSM7140516_CON8_synovium    GSM7140516         CON8 synovium        Control_52W             8
9      GSM7140518_CON9_synovium    GSM7140518         CON9 synovium        Control_52W             9
10    GSM7140523_CON10_synovium    GSM7140523        CON10 synovium        Control_52W            10
11    GSM7140527_CON11_synovium    GSM7140527        CON11 synovium        Control_52W            11
12    GSM7140532_CON12_synovium    GSM7140532        CON12 synovium        Control_52W            12
13    GSM7140487_ACLT1_synovium    GSM7140487        ACLT1 synovium     ACLT_alone_52W             1
14    GSM7140488_ACLT2_synovium    GSM7140488        ACLT2 synovium     ACLT_alone_52W             2
15    GSM7140489_ACLT3_synovium    GSM7140489        ACLT3 synovium     ACLT_alone_52W             3
16    GSM7140490_ACLT4_synovium    GSM7140490        ACLT4 synovium     ACLT_alone_52W             4
17    GSM7140492_ACLT5_synovium    GSM7140492        ACLT5 synovium     ACLT_alone_52W             5
18    GSM7140493_ACLT6_synovium    GSM7140493        ACLT6 synovium     ACLT_alone_52W             6
19    GSM7140494_ACLT7_synovium    GSM7140494        ACLT7 synovium     ACLT_alone_52W             7
20    GSM7140496_ACLT8_synovium    GSM7140496        ACLT8 synovium     ACLT_alone_52W             8
21    GSM7140497_ACLT9_synovium    GSM7140497        ACLT9 synovium     ACLT_alone_52W             9
22   GSM7140498_ACLT10_synovium    GSM7140498       ACLT10 synovium     ACLT_alone_52W            10
23   GSM7140499_ACLT11_synovium    GSM7140499       ACLT11 synovium     ACLT_alone_52W            11
24   GSM7140501_ACLT12_synovium    GSM7140501       ACLT12 synovium     ACLT_alone_52W            12
25   GSM7140503_RECON1_synovium    GSM7140503       RECON1 synovium Reconstruction_52W             1
26   GSM7140504_RECON2_synovium    GSM7140504       RECON2 synovium Reconstruction_52W             2
27   GSM7140505_RECON3_synovium    GSM7140505       RECON3 synovium Reconstruction_52W             3
28   GSM7140506_RECON4_synovium    GSM7140506       RECON4 synovium Reconstruction_52W             4
29   GSM7140508_RECON5_synovium    GSM7140508       RECON5 synovium Reconstruction_52W             5
30   GSM7140509_RECON6_synovium    GSM7140509       RECON6 synovium Reconstruction_52W             6
31   GSM7140510_RECON7_synovium    GSM7140510       RECON7 synovium Reconstruction_52W             7
32   GSM7140512_RECON8_synovium    GSM7140512       RECON8 synovium Reconstruction_52W             8
33   GSM7140513_RECON9_synovium    GSM7140513       RECON9 synovium Reconstruction_52W             9
34  GSM7140514_RECON10_synovium    GSM7140514      RECON10 synovium Reconstruction_52W            10
35  GSM7140515_RECON11_synovium    GSM7140515      RECON11 synovium Reconstruction_52W            11
36  GSM7140517_RECON12_synovium    GSM7140517      RECON12 synovium Reconstruction_52W            12
37  GSM7140519_REPAIR1_synovium    GSM7140519      REPAIR1 synovium         Repair_52W             1
38  GSM7140520_REPAIR2_synovium    GSM7140520      REPAIR2 synovium         Repair_52W             2
39  GSM7140521_REPAIR3_synovium    GSM7140521      REPAIR3 synovium         Repair_52W             3
40  GSM7140522_REPAIR4_synovium    GSM7140522      REPAIR4 synovium         Repair_52W             4
41  GSM7140524_REPAIR5_synovium    GSM7140524      REPAIR5 synovium         Repair_52W             5
42  GSM7140525_REPAIR6_synovium    GSM7140525      REPAIR6 synovium         Repair_52W             6
43  GSM7140526_REPAIR7_synovium    GSM7140526      REPAIR7 synovium         Repair_52W             7
44  GSM7140528_REPAIR8_synovium    GSM7140528      REPAIR8 synovium         Repair_52W             8
45  GSM7140529_REPAIR9_synovium    GSM7140529      REPAIR9 synovium         Repair_52W             9
46 GSM7140530_REPAIR10_synovium    GSM7140530     REPAIR10 synovium         Repair_52W            10
47 GSM7140531_REPAIR11_synovium    GSM7140531     REPAIR11 synovium         Repair_52W            11
48 GSM7140533_REPAIR12_synovium    GSM7140533     REPAIR12 synovium         Repair_52W            12

Synovium-only group counts:

    ACLT_alone_52W        Control_52W Reconstruction_52W         Repair_52W 
                12                 12                 12                 12 

Quant-file audit for locked 48 synovium matrix:
All quant-like files found: 96 
Synovium quant-like files found: 48 
Locked synovium matrix samples with matching synovium quant files: 48 

Metadata 24 audit:
Metadata rows: 24 
Rows in all 96 matrix: 24 
Rows in locked synovium 48 matrix: 24 
Rows in main 24 matrix: 24 

Locked synovium 48 sample QC preview:
                   sample_id geo_accession sample_label   tissue          group sample_number total_estimated_counts detected_genes_count_gt0 detected_genes_count_ge1
1   GSM7140486_CON1_synovium    GSM7140486         CON1 synovium    Control_52W             1                7108677                    16731                    16728
2   GSM7140491_CON2_synovium    GSM7140491         CON2 synovium    Control_52W             2                8959039                    16672                    16664
3   GSM7140495_CON3_synovium    GSM7140495         CON3 synovium    Control_52W             3                9719015                    16993                    16984
4   GSM7140500_CON4_synovium    GSM7140500         CON4 synovium    Control_52W             4                8391254                    16672                    16671
5   GSM7140502_CON5_synovium    GSM7140502         CON5 synovium    Control_52W             5               12170231                    16538                    16532
6   GSM7140507_CON6_synovium    GSM7140507         CON6 synovium    Control_52W             6               16132991                    17390                    17383
7   GSM7140511_CON7_synovium    GSM7140511         CON7 synovium    Control_52W             7               11920148                    17215                    17210
8   GSM7140516_CON8_synovium    GSM7140516         CON8 synovium    Control_52W             8                9996314                    16810                    16803
9   GSM7140518_CON9_synovium    GSM7140518         CON9 synovium    Control_52W             9               12285504                    17058                    17054
10 GSM7140523_CON10_synovium    GSM7140523        CON10 synovium    Control_52W            10               13966921                    16999                    16997
11 GSM7140527_CON11_synovium    GSM7140527        CON11 synovium    Control_52W            11               16354930                    16694                    16688
12 GSM7140532_CON12_synovium    GSM7140532        CON12 synovium    Control_52W            12               11603076                    16604                    16597
13 GSM7140487_ACLT1_synovium    GSM7140487        ACLT1 synovium ACLT_alone_52W             1               10334100                    16660                    16658
14 GSM7140488_ACLT2_synovium    GSM7140488        ACLT2 synovium ACLT_alone_52W             2               10521457                    16880                    16877
15 GSM7140489_ACLT3_synovium    GSM7140489        ACLT3 synovium ACLT_alone_52W             3               18669780                    16495                    16493
16 GSM7140490_ACLT4_synovium    GSM7140490        ACLT4 synovium ACLT_alone_52W             4                9298837                    16367                    16362
17 GSM7140492_ACLT5_synovium    GSM7140492        ACLT5 synovium ACLT_alone_52W             5                9066422                    16613                    16605
18 GSM7140493_ACLT6_synovium    GSM7140493        ACLT6 synovium ACLT_alone_52W             6               10831988                    16291                    16289
19 GSM7140494_ACLT7_synovium    GSM7140494        ACLT7 synovium ACLT_alone_52W             7               10627443                    16710                    16701
20 GSM7140496_ACLT8_synovium    GSM7140496        ACLT8 synovium ACLT_alone_52W             8               12200497                    16912                    16905
   tpm_sum_gene_level genes_with_tpm_gt0 effective_length_na_n effective_length_nonpositive_n effective_length_median
1            981863.0              16731                     0                              0                2068.981
2            979800.5              16672                     0                              0                2070.901
3            978064.3              16993                     0                              0                2001.266
4            981488.1              16672                     0                              0                2120.246
5            973601.7              16538                     0                              0                2014.147
6            978012.3              17390                     0                              0                2095.801
7            979716.8              17215                     0                              0                2075.162
8            976774.5              16810                     0                              0                2150.801
9            975982.7              17058                     0                              0                2064.551
10           973436.5              16999                     0                              0                2066.693
11           979477.6              16694                     0                              0                2043.190
12           978679.6              16604                     0                              0                2056.304
13           980475.3              16660                     0                              0                2037.955
14           979505.1              16880                     0                              0                2062.312
15           979341.4              16495                     0                              0                2065.769
16           976777.6              16367                     0                              0                2087.923
17           977392.2              16613                     0                              0                2090.186
18           975265.8              16291                     0                              0                1959.066
19           986966.1              16710                     0                              0                2012.780
20           978661.9              16912                     0                              0                2080.269

Locked synovium 48 group QC summary:
# A tibble: 4 × 15
  group       n_samples mean_total_estimated…¹ median_total_estimat…² min_total_estimated_…³ max_total_estimated_…⁴ mean_detected_genes_…⁵ median_detected_gene…⁶ min_detected_genes_c…⁷
  <chr>           <int>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <int>
1 ACLT_alone…        12              10601788.              10427778.               7706256.              18669780.                 16438.                 16514                   15555
2 Control_52W        12              11550675.              11761612.               7108677.              16354930.                 16865.                 16770.                  16538
3 Reconstruc…        12              13830394.              14124377.               8900906.              17336700.                 16498.                 16566.                  16005
4 Repair_52W         12              14067516.              13727374.              11246102.              17197856.                 16693.                 16764                   16198
# ℹ abbreviated names: ¹​mean_total_estimated_counts, ²​median_total_estimated_counts, ³​min_total_estimated_counts, ⁴​max_total_estimated_counts, ⁵​mean_detected_genes_count_gt0,
#   ⁶​median_detected_genes_count_gt0, ⁷​min_detected_genes_count_gt0
# ℹ 6 more variables: max_detected_genes_count_gt0 <int>, mean_gene_level_TPM_sum <dbl>, min_gene_level_TPM_sum <dbl>, max_gene_level_TPM_sum <dbl>, total_effective_length_na_n <int>,
#   total_effective_length_nonpositive_n <int>

Main comparison 24 group QC summary:
# A tibble: 2 × 15
  group       n_samples mean_total_estimated…¹ median_total_estimat…² min_total_estimated_…³ max_total_estimated_…⁴ mean_detected_genes_…⁵ median_detected_gene…⁶ min_detected_genes_c…⁷
  <chr>           <int>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <dbl>                  <int>
1 ACLT_alone…        12              10601788.              10427778.               7706256.              18669780.                 16438.                 16514                   15555
2 Control_52W        12              11550675.              11761612.               7108677.              16354930.                 16865.                 16770.                  16538
# ℹ abbreviated names: ¹​mean_total_estimated_counts, ²​median_total_estimated_counts, ³​min_total_estimated_counts, ⁴​max_total_estimated_counts, ⁵​mean_detected_genes_count_gt0,
#   ⁶​median_detected_genes_count_gt0, ⁷​min_detected_genes_count_gt0
# ℹ 6 more variables: max_detected_genes_count_gt0 <int>, mean_gene_level_TPM_sum <dbl>, min_gene_level_TPM_sum <dbl>, max_gene_level_TPM_sum <dbl>, total_effective_length_na_n <int>,
#   total_effective_length_nonpositive_n <int>

Final summary for review:
                                             metric                                                                                                 value
1                                        run_status                                                                                               SUCCESS
2                                    method_version                    2026-05-10_lock_chronic_pig_synovium_only_48_matrices_from_step24_tximport_outputs
3                           all_step24_matrix_genes                                                                                                 22438
4                         all_step24_matrix_samples                                                                                                    96
5                      locked_synovium_matrix_genes                                                                                                 22438
6                    locked_synovium_matrix_samples                                                                                                    48
7                      locked_synovium_group_counts                               ACLT_alone_52W 12; Control_52W 12; Reconstruction_52W 12; Repair_52W 12
8                      main_comparison_matrix_genes                                                                                                 22438
9                    main_comparison_matrix_samples                                                                                                    24
10                     main_comparison_group_counts                                                                     ACLT_alone_52W 12; Control_52W 12
11                       synovium_quant_files_found                                                                                                    48
12 locked_synovium_samples_with_matching_quant_file                                                                                                    48
13                                 metadata_24_rows                                                                                                    24
14         metadata_24_in_locked_synovium_48_matrix                                                                                                    24
15                    metadata_24_in_main_24_matrix                                                                                                    24
16  counts_abundance_length_gene_id_order_identical                                                                                                  TRUE
17 counts_abundance_length_sample_columns_identical                                                                                                  TRUE
18          mean_synovium_48_total_estimated_counts                                                                                          12512593.483
19           min_synovium_48_total_estimated_counts                                                                                           7108677.116
20           max_synovium_48_total_estimated_counts                                                                                          18669779.773
21        mean_synovium_48_detected_genes_count_gt0                                                                                             16623.167
22         min_synovium_48_detected_genes_count_gt0                                                                                                 15555
23         max_synovium_48_detected_genes_count_gt0                                                                                                 17390
24                                      output_root E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock

Key output files:
1) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_synovium_48_gene_level_estimated_counts_matrix.csv
2) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_synovium_48_gene_level_abundance_tpm_matrix.csv
3) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_synovium_48_gene_level_effective_length_matrix.csv
4) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_synovium_48_sample_metadata_locked.csv
5) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_gene_level_estimated_counts_matrix.csv
6) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_gene_level_abundance_tpm_matrix.csv
7) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_gene_level_effective_length_matrix.csv
8) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_sample_metadata_locked.csv
9) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_locked_matrix_dimensions.csv
10) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_final_summary_for_review.csv

Session information:
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8    LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_1.3.2      Rsubread_2.24.0  dplyr_1.1.4      ggplot2_4.0.2    singscore_1.30.0 edgeR_4.8.2      limma_3.66.0    

loaded via a namespace (and not attached):
 [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 gtable_0.3.6                Biobase_2.70.0              lattice_0.22-7              vctrs_0.7.2                
 [7] tools_4.5.2                 generics_0.1.4              parallel_4.5.2              stats4_4.5.2                tibble_3.3.1                AnnotationDbi_1.72.0       
[13] RSQLite_2.4.5               blob_1.3.0                  pkgconfig_2.0.3             Matrix_1.7-4                RColorBrewer_1.1-3          S7_0.2.1                   
[19] S4Vectors_0.48.0            graph_1.88.1                lifecycle_1.0.5             stringr_1.6.0               compiler_4.5.2              farver_2.1.2               
[25] Biostrings_2.78.0           statmod_1.5.1               Seqinfo_1.0.0               pillar_1.11.1               crayon_1.5.3                cachem_1.1.0               
[31] DelayedArray_0.36.0         abind_1.4-8                 tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7               reshape2_1.4.5             
[37] purrr_1.2.1                 labeling_0.4.3              fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                   SparseArray_1.10.8         
[43] magrittr_2.0.4              S4Arrays_1.10.1             XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0             withr_3.0.2                
[49] scales_1.4.0                bit64_4.6.0-1               XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0           bit_4.6.0                  
[55] png_0.1-8                   memoise_2.0.1               GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                 Rcpp_1.1.1                 
[61] xtable_1.8-4                glue_1.8.0                  DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0             plyr_1.8.9                 
[67] R6_2.6.1                    MatrixGenerics_1.22.0      
============================================================
StepCHRONMDS_01_chronic_pig_main24_PCA_MDS
Started at: 2026-05-10 21:03:02 
Method version: 2026-05-10_chronic_pig_main24_PCA_MDS_unified_with_mouse_and_early_pig 
============================================================

Script archive note: --file argument was not detected. If running interactively, manually save this script in:
E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/scripts 

Required packages loaded.
edgeR version: 4.8.2 
limma version: 3.66.0 
ggplot2 version: 4.0.2 
ggrepel available for label plots: TRUE 

Input files:
Counts matrix: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_gene_level_estimated_counts_matrix.csv 
Metadata: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_sample_metadata_locked.csv 
Output root: E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA 

Counts table dimension: 22438 rows x 25 columns
First 10 count table columns:
 [1] "gene_id"                  "GSM7140486_CON1_synovium" "GSM7140491_CON2_synovium"
 [4] "GSM7140495_CON3_synovium" "GSM7140500_CON4_synovium" "GSM7140502_CON5_synovium"
 [7] "GSM7140507_CON6_synovium" "GSM7140511_CON7_synovium" "GSM7140516_CON8_synovium"
[10] "GSM7140518_CON9_synovium"

Metadata dimension: 24 rows x 6 columns
Metadata columns:
[1] "sample_id"     "geo_accession" "sample_label"  "tissue"        "group"         "sample_number"
Inferred gene ID column: gene_id 
Inferred metadata sample column: sample_id 
Inferred metadata group column: group 

Metadata sample count: 24 
Count matrix sample columns: 24 
Missing metadata samples in count matrix: none 
Extra count matrix samples not in metadata: none 

Validated analysis count matrix dimension: 22438 genes x 24 samples
Group counts:

   Control_52W ACLT_alone_52W 
            12             12 

Filtering and normalization summary:
                          metric
1                 method_version
2             input_count_matrix
3                 input_metadata
4      genes_before_filterByExpr
5       genes_after_filterByExpr
6  genes_removed_by_filterByExpr
7                        samples
8                         groups
9           normalization_method
10            logCPM_prior_count
11                    PCA_method
12                    MDS_method
13                     MDS_input
14                 MDS_top_genes
15            MDS_gene_selection
16   Salmon_tximport_counts_note
                                                                                                                                                                          value
1                                                                                                        2026-05-10_chronic_pig_main24_PCA_MDS_unified_with_mouse_and_early_pig
2  E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_gene_level_estimated_counts_matrix.csv
3              E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_sample_metadata_locked.csv
4                                                                                                                                                                         22438
5                                                                                                                                                                         14403
6                                                                                                                                                                          8035
7                                                                                                                                                                            24
8                                                                                                                                             Control_52W 12; ACLT_alone_52W 12
9                                                                                                                                                                     edgeR_TMM
10                                                                                                                                                                            1
11                                                                                                                    prcomp(t(filtered_TMM_logCPM), center=TRUE, scale.=FALSE)
12                                                                                                        limma::plotMDS(filtered_TMM_logCPM, top=500, gene.selection='common')
13                                                                                                                                        filtered TMM-normalized logCPM matrix
14                                                                                                                                                                          500
15                                                                                                                                                                       common
16                                                                                      Estimated counts from Salmon/tximport were used as numeric values and were not rounded.

Library summary after TMM:
                                            sample_id          group raw_estimated_count_sum
GSM7140486_CON1_synovium     GSM7140486_CON1_synovium    Control_52W                 7108677
GSM7140491_CON2_synovium     GSM7140491_CON2_synovium    Control_52W                 8959039
GSM7140495_CON3_synovium     GSM7140495_CON3_synovium    Control_52W                 9719015
GSM7140500_CON4_synovium     GSM7140500_CON4_synovium    Control_52W                 8391254
GSM7140502_CON5_synovium     GSM7140502_CON5_synovium    Control_52W                12170231
GSM7140507_CON6_synovium     GSM7140507_CON6_synovium    Control_52W                16132991
GSM7140511_CON7_synovium     GSM7140511_CON7_synovium    Control_52W                11920148
GSM7140516_CON8_synovium     GSM7140516_CON8_synovium    Control_52W                 9996314
GSM7140518_CON9_synovium     GSM7140518_CON9_synovium    Control_52W                12285504
GSM7140523_CON10_synovium   GSM7140523_CON10_synovium    Control_52W                13966921
GSM7140527_CON11_synovium   GSM7140527_CON11_synovium    Control_52W                16354930
GSM7140532_CON12_synovium   GSM7140532_CON12_synovium    Control_52W                11603076
GSM7140487_ACLT1_synovium   GSM7140487_ACLT1_synovium ACLT_alone_52W                10334100
GSM7140488_ACLT2_synovium   GSM7140488_ACLT2_synovium ACLT_alone_52W                10521457
GSM7140489_ACLT3_synovium   GSM7140489_ACLT3_synovium ACLT_alone_52W                18669780
GSM7140490_ACLT4_synovium   GSM7140490_ACLT4_synovium ACLT_alone_52W                 9298837
GSM7140492_ACLT5_synovium   GSM7140492_ACLT5_synovium ACLT_alone_52W                 9066422
GSM7140493_ACLT6_synovium   GSM7140493_ACLT6_synovium ACLT_alone_52W                10831988
GSM7140494_ACLT7_synovium   GSM7140494_ACLT7_synovium ACLT_alone_52W                10627443
GSM7140496_ACLT8_synovium   GSM7140496_ACLT8_synovium ACLT_alone_52W                12200497
GSM7140497_ACLT9_synovium   GSM7140497_ACLT9_synovium ACLT_alone_52W                 7706256
GSM7140498_ACLT10_synovium GSM7140498_ACLT10_synovium ACLT_alone_52W                 9007486
GSM7140499_ACLT11_synovium GSM7140499_ACLT11_synovium ACLT_alone_52W                 8387249
GSM7140501_ACLT12_synovium GSM7140501_ACLT12_synovium ACLT_alone_52W                10569945
                           filtered_estimated_count_sum norm_factor_TMM effective_library_size
GSM7140486_CON1_synovium                        7095768       1.1101154                7877121
GSM7140491_CON2_synovium                        8947879       1.0398987                9304888
GSM7140495_CON3_synovium                        9703857       1.0996272               10670625
GSM7140500_CON4_synovium                        8379766       1.0719755                8982903
GSM7140502_CON5_synovium                       12159517       0.9604690               11678839
GSM7140507_CON6_synovium                       16100387       1.1194556               18023669
GSM7140511_CON7_synovium                       11902507       1.1107140               13220282
GSM7140516_CON8_synovium                        9983675       1.1262162               11243776
GSM7140518_CON9_synovium                       12268106       1.0968281               13456003
GSM7140523_CON10_synovium                      13947970       1.0776925               15031622
GSM7140527_CON11_synovium                      16340716       0.8937062               14603799
GSM7140532_CON12_synovium                      11587087       1.0909146               12640522
GSM7140487_ACLT1_synovium                      10309179       0.8850194                9123823
GSM7140488_ACLT2_synovium                      10495846       0.9955605               10449249
GSM7140489_ACLT3_synovium                      18653370       0.9271038               17293610
GSM7140490_ACLT4_synovium                       9289874       0.9943548                9237430
GSM7140492_ACLT5_synovium                       9055576       0.9985963                9042865
GSM7140493_ACLT6_synovium                      10822558       0.9057636                9802679
GSM7140494_ACLT7_synovium                      10613490       0.9235909                9802523
GSM7140496_ACLT8_synovium                      12169089       1.0304370               12539480
GSM7140497_ACLT9_synovium                       7697083       1.0103289                7776585
GSM7140498_ACLT10_synovium                      8996313       1.0381669                9339674
GSM7140499_ACLT11_synovium                      8380296       0.7937897                6652193
GSM7140501_ACLT12_synovium                     10509803       0.8144152                8559343
                           detected_genes_count_gt0_raw detected_genes_count_gt0_filtered
GSM7140486_CON1_synovium                          16731                             14385
GSM7140491_CON2_synovium                          16672                             14382
GSM7140495_CON3_synovium                          16993                             14385
GSM7140500_CON4_synovium                          16672                             14388
GSM7140502_CON5_synovium                          16538                             14384
GSM7140507_CON6_synovium                          17390                             14398
GSM7140511_CON7_synovium                          17215                             14394
GSM7140516_CON8_synovium                          16810                             14387
GSM7140518_CON9_synovium                          17058                             14389
GSM7140523_CON10_synovium                         16999                             14396
GSM7140527_CON11_synovium                         16694                             14379
GSM7140532_CON12_synovium                         16604                             14368
GSM7140487_ACLT1_synovium                         16660                             14376
GSM7140488_ACLT2_synovium                         16880                             14391
GSM7140489_ACLT3_synovium                         16495                             14379
GSM7140490_ACLT4_synovium                         16367                             14363
GSM7140492_ACLT5_synovium                         16613                             14386
GSM7140493_ACLT6_synovium                         16291                             14361
GSM7140494_ACLT7_synovium                         16710                             14381
GSM7140496_ACLT8_synovium                         16912                             14393
GSM7140497_ACLT9_synovium                         16517                             14384
GSM7140498_ACLT10_synovium                        16511                             14373
GSM7140499_ACLT11_synovium                        15741                             14293
GSM7140501_ACLT12_synovium                        15555                             14248

PCA variance explained by first five PCs:
   PC variance_fraction variance_percent
1 PC1        0.26311098        26.311098
2 PC2        0.12713243        12.713243
3 PC3        0.09648993         9.648993
4 PC4        0.07997444         7.997444
5 PC5        0.04872964         4.872964

MDS object diagnostic:
MDS function used: limma::plotMDS
MDS input object: filtered TMM-normalized logCPM matrix
length(mds_obj$x): 24 
length(mds_obj$y): 24 
length(colnames(logcpm)): 24 
MDS coordinate table preview:
                    sample_id       MDS1        MDS2          group
1    GSM7140486_CON1_synovium -0.7323405  0.63818484    Control_52W
2    GSM7140491_CON2_synovium -0.8205781 -0.43598238    Control_52W
3    GSM7140495_CON3_synovium -0.6862623  0.58739419    Control_52W
4    GSM7140500_CON4_synovium -0.7545569  0.18701754    Control_52W
5    GSM7140502_CON5_synovium -0.3690633  0.23804792    Control_52W
6    GSM7140507_CON6_synovium -1.3714194 -0.31041161    Control_52W
7    GSM7140511_CON7_synovium -1.0466060 -0.23171280    Control_52W
8    GSM7140516_CON8_synovium -1.1225941 -0.11587839    Control_52W
9    GSM7140518_CON9_synovium -0.9816477  0.16076970    Control_52W
10  GSM7140523_CON10_synovium -1.4262350 -0.05077962    Control_52W
11  GSM7140527_CON11_synovium  1.1462793 -0.10076299    Control_52W
12  GSM7140532_CON12_synovium -1.3715130  0.61880493    Control_52W
13  GSM7140487_ACLT1_synovium  1.5263342  1.33739032 ACLT_alone_52W
14  GSM7140488_ACLT2_synovium  0.8825123  1.22727402 ACLT_alone_52W
15  GSM7140489_ACLT3_synovium  1.0117488 -1.03664886 ACLT_alone_52W
16  GSM7140490_ACLT4_synovium  0.7051015 -0.91733837 ACLT_alone_52W
17  GSM7140492_ACLT5_synovium  0.7612967  0.04826654 ACLT_alone_52W
18  GSM7140493_ACLT6_synovium  0.5191446 -1.11764446 ACLT_alone_52W
19  GSM7140494_ACLT7_synovium  0.5207194 -0.45811418 ACLT_alone_52W
20  GSM7140496_ACLT8_synovium  0.7931422  1.32682022 ACLT_alone_52W
21  GSM7140497_ACLT9_synovium  0.6704738  0.67457053 ACLT_alone_52W
22 GSM7140498_ACLT10_synovium  0.2420734 -0.62816862 ACLT_alone_52W
23 GSM7140499_ACLT11_synovium  1.4144023 -0.21934301 ACLT_alone_52W
24 GSM7140501_ACLT12_synovium  0.4895878 -1.42175548 ACLT_alone_52W

Final summary for review:
                      metric
1                 run_status
2             method_version
3          input_counts_file
4        input_metadata_file
5                output_root
6  genes_before_filterByExpr
7   genes_after_filterByExpr
8                    samples
9                     groups
10      normalization_method
11           PCA_PC1_percent
12           PCA_PC2_percent
13              MDS_function
14                 MDS_input
15                  MDS1_min
16                  MDS1_max
17                  MDS2_min
18                  MDS2_max
19  estimated_counts_rounded
20       group_Control_52W_n
21    group_ACLT_alone_52W_n
                                                                                                                                                                          value
1                                                                                                                                                                       SUCCESS
2                                                                                                        2026-05-10_chronic_pig_main24_PCA_MDS_unified_with_mouse_and_early_pig
3  E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_gene_level_estimated_counts_matrix.csv
4              E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/tables/step24_chronic_synovium_only_matrix_lock/tables/StepCHRONLOCK_01_main_24_sample_metadata_locked.csv
5                                                                                                      E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA
6                                                                                                                                                                         22438
7                                                                                                                                                                         14403
8                                                                                                                                                                            24
9                                                                                                                                             Control_52W 12; ACLT_alone_52W 12
10                                                                                                                                                                    edgeR_TMM
11                                                                                                                                                                    26.311098
12                                                                                                                                                                    12.713243
13                                                                                                                                                               limma::plotMDS
14                                                                                                                                        filtered TMM-normalized logCPM matrix
15                                                                                                                                                                    -1.426235
16                                                                                                                                                                     1.526334
17                                                                                                                                                                    -1.421755
18                                                                                                                                                                     1.337390
19                                                                                                                                                                        FALSE
20                                                                                                                                                                           12
21                                                                                                                                                                           12

Key output files:
1) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/source_data/StepCHRONMDS_01_main24_filtered_TMM_logCPM_matrix.csv
2) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/source_data/StepCHRONMDS_01_main24_PCA_coordinates.csv
3) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/source_data/StepCHRONMDS_01_main24_MDS_coordinates.csv
4) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/tables/StepCHRONMDS_01_main24_PCA_variance_explained.csv
5) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/tables/StepCHRONMDS_01_main24_library_summary_after_TMM.csv
6) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_PCA_unlabeled.pdf
7) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_PCA_labeled.pdf
8) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_MDS_unlabeled.pdf
9) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_MDS_labeled.pdf
10) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/tables/StepCHRONMDS_01_final_summary_for_review.csv

Session information:
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_1.3.2      Rsubread_2.24.0  dplyr_1.1.4      ggplot2_4.0.2    singscore_1.30.0
[6] edgeR_4.8.2      limma_3.66.0    

loaded via a namespace (and not attached):
 [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 gtable_0.3.6               
 [4] ggrepel_0.9.6               Biobase_2.70.0              lattice_0.22-7             
 [7] vctrs_0.7.2                 tools_4.5.2                 generics_0.1.4             
[10] parallel_4.5.2              stats4_4.5.2                tibble_3.3.1               
[13] AnnotationDbi_1.72.0        RSQLite_2.4.5               blob_1.3.0                 
[16] pkgconfig_2.0.3             Matrix_1.7-4                RColorBrewer_1.1-3         
[19] S7_0.2.1                    S4Vectors_0.48.0            graph_1.88.1               
[22] lifecycle_1.0.5             stringr_1.6.0               compiler_4.5.2             
[25] farver_2.1.2                Biostrings_2.78.0           statmod_1.5.1              
[28] Seqinfo_1.0.0               pillar_1.11.1               crayon_1.5.3               
[31] cachem_1.1.0                DelayedArray_0.36.0         abind_1.4-8                
[34] tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7              
[37] reshape2_1.4.5              purrr_1.2.1                 labeling_0.4.3             
[40] fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                  
[43] SparseArray_1.10.8          magrittr_2.0.4              S4Arrays_1.10.1            
[46] XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0            
[49] withr_3.0.2                 scales_1.4.0                bit64_4.6.0-1              
[52] XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0          
[55] bit_4.6.0                   png_0.1-8                   memoise_2.0.1              
[58] GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                
[61] Rcpp_1.1.1                  xtable_1.8-4                glue_1.8.0                 
[64] DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0            
[67] plyr_1.8.9                  R6_2.6.1                    MatrixGenerics_1.22.0      
Figure1A 3-block linebreak v2 completed successfully.
PNG: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/06_figures/Figure1/Figure1A_overall_study_workflow_3block_linebreak_v2.png
PDF: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/06_figures/Figure1/Figure1A_overall_study_workflow_3block_linebreak_v2.pdf
Summary log: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/08_logs/step01_Figure1A_overall_study_workflow_3block_linebreak_v2_summary_to_send.txt
_1.2.1            locfit_1.5-9.12             stringi_1.8.7              
[37] reshape2_1.4.5              purrr_1.2.1                 labeling_0.4.3             
[40] fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                  
[43] SparseArray_1.10.8          magrittr_2.0.4              S4Arrays_1.10.1            
[46] XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0            
[49] withr_3.0.2                 scales_1.4.0                bit64_4.6.0-1              
[52] XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0          
[55] bit_4.6.0                   png_0.1-8                   memoise_2.0.1              
[58] GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                
[61] Rcpp_1.1.1                  xtable_1.8-4                glue_1.8.0                 
[64] DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0            
[67] plyr_1.8.9                  R6_2.6.1                    MatrixGenerics_1.22.0      
Figure1A 3-block linebreak v2 completed successfully.
PNG: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/06_figures/Figure1/Figure1A_overall_study_workflow_3block_linebreak_v2.png
PDF: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/06_figures/Figure1/Figure1A_overall_study_workflow_3block_linebreak_v2.pdf
Summary log: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/08_logs/step01_Figure1A_overall_study_workflow_3block_linebreak_v2_summary_to_send.txt
lement/MDS_PCA/tables/StepCHRONMDS_01_main24_PCA_variance_explained.csv
5) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/tables/StepCHRONMDS_01_main24_library_summary_after_TMM.csv
6) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_PCA_unlabeled.pdf
7) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_PCA_labeled.pdf
8) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_MDS_unlabeled.pdf
9) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/figures/StepCHRONMDS_01_main24_MDS_labeled.pdf
10) E:/R/ACLsenescence2/rebuild_submission/03_pig_chronic/supplement/MDS_PCA/tables/StepCHRONMDS_01_final_summary_for_review.csv

Session information:
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tidyr_1.3.2      Rsubread_2.24.0  dplyr_1.1.4      ggplot2_4.0.2    singscore_1.30.0
[6] edgeR_4.8.2      limma_3.66.0    

loaded via a namespace (and not attached):
 [1] KEGGREST_1.50.0             SummarizedExperiment_1.40.0 gtable_0.3.6               
 [4] ggrepel_0.9.6               Biobase_2.70.0              lattice_0.22-7             
 [7] vctrs_0.7.2                 tools_4.5.2                 generics_0.1.4             
[10] parallel_4.5.2              stats4_4.5.2                tibble_3.3.1               
[13] AnnotationDbi_1.72.0        RSQLite_2.4.5               blob_1.3.0                 
[16] pkgconfig_2.0.3             Matrix_1.7-4                RColorBrewer_1.1-3         
[19] S7_0.2.1                    S4Vectors_0.48.0            graph_1.88.1               
[22] lifecycle_1.0.5             stringr_1.6.0               compiler_4.5.2             
[25] farver_2.1.2                Biostrings_2.78.0           statmod_1.5.1              
[28] Seqinfo_1.0.0               pillar_1.11.1               crayon_1.5.3               
[31] cachem_1.1.0                DelayedArray_0.36.0         abind_1.4-8                
[34] tidyselect_1.2.1            locfit_1.5-9.12             stringi_1.8.7              
[37] reshape2_1.4.5              purrr_1.2.1                 labeling_0.4.3             
[40] fastmap_1.2.0               grid_4.5.2                  cli_3.6.5                  
[43] SparseArray_1.10.8          magrittr_2.0.4              S4Arrays_1.10.1            
[46] XML_3.99-0.20               utf8_1.2.6                  GSEABase_1.72.0            
[49] withr_3.0.2                 scales_1.4.0                bit64_4.6.0-1              
[52] XVector_0.50.0              httr_1.4.7                  matrixStats_1.5.0          
[55] bit_4.6.0                   png_0.1-8                   memoise_2.0.1              
[58] GenomicRanges_1.62.1        IRanges_2.44.0              rlang_1.1.7                
[61] Rcpp_1.1.1                  xtable_1.8-4                glue_1.8.0                 
[64] DBI_1.2.3                   BiocGenerics_0.56.0         annotate_1.88.0            
[67] plyr_1.8.9                  R6_2.6.1                    MatrixGenerics_1.22.0      
Figure1A 3-block linebreak v2 completed successfully.
PNG: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/06_figures/Figure1/Figure1A_overall_study_workflow_3block_linebreak_v2.png
PDF: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/06_figures/Figure1/Figure1A_overall_study_workflow_3block_linebreak_v2.pdf
Summary log: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/08_logs/step01_Figure1A_overall_study_workflow_3block_linebreak_v2_summary_to_send.txt
