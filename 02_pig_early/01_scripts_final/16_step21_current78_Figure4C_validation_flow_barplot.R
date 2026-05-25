===== STEP21 CURRENT78 FIGURE4C CORE HEATMAP =====
Validation file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_validation_table.csv
Summary file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_early_signature_remap_summary.csv
Count file: E:/R/ACLsenescence2 LD/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv
Manifest file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step09_pig_early_core_sample_fastq_manifest.csv
Step19 score annotation file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv
indexing step18_current78_pig_signature_validation_table.csv [=====================================================================================] 1.38GB/s, eta:  0s                                                                                                                                                                       indexing step18_current78_pig_early_signature_remap_summary.csv [=================================================================================] 63.61MB/s, eta:  0s                                                                                                                                                                       Core genes identified: 24
indexing step16_pig_early_gene_count_matrix.csv [==================================================================================================] 2.15GB/s, eta:  0s                                                                                                                                                                       No gene symbol column found in count matrix; this is acceptable because Step21 uses pig_ensg as the matching key.
indexing step19_current78_pig_signature_scores_by_sample.csv [===================================================================================] 198.77MB/s, eta:  0s                                                                                                                                                                       Sample groups assigned from Step19 score table.
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_progress_78_to_24_summary.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_gene_table.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_heatmap_zscore_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_pig_core_ortholog_sample_annotation.csv
null device 
          1 
null device 
          1 
[1] TRUE
[1] TRUE
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/step21_current78_Figure4C_core_heatmap_summary.csv

===== STEP21 CURRENT78 FIGURE4C SUMMARY =====
                                      metric
1  current_human_CellAge_overlap_input_genes
2       current_pig_signature_genes_detected
3       direction_consistent_both_timepoints
4                    core_strict_both_t7_t28
5           core_genes_found_in_count_matrix
6                                   sample_n
7                                  control_n
8                                  aclt_t7_n
9                                 aclt_t28_n
10                               heatmap_png
11                               heatmap_pdf
12                     canonical_heatmap_png
13                     canonical_heatmap_pdf
14                                output_dir
                                                                                                                        value
1                                                                                                                          78
2                                                                                                                          75
3                                                                                                                          47
4                                                                                                                          24
5                                                                                                                          24
6                                                                                                                          18
7                                                                                                                           6
8                                                                                                                           6
9                                                                                                                           6
10 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_current78_pig_early_core_ortholog_heatmap.png
11 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_current78_pig_early_core_ortholog_heatmap.pdf
12           E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_pig_early_core_ortholog_heatmap.png
13           E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_pig_early_core_ortholog_heatmap.pdf
14                          E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap

Core genes preview:
# A tibble: 20 × 9
   human_gene mouse_symbol pig_ensg           pig_symbol signature_direction logFC_t7      FDR_t7 logFC_t28     FDR_t28
   <chr>      <chr>        <chr>              <chr>      <chr>                  <dbl>       <dbl>     <dbl>       <dbl>
 1 CBS        Cbs          ENSSSCG00000040779 CBS        Down_in_ACLR           -1.90 0.000653        -3.11 0.0000993  
 2 FASN       Fasn         ENSSSCG00000029944 FASN       Down_in_ACLR           -2.11 0.000494        -3.93 0.000000106
 3 PPARG      Pparg        ENSSSCG00000011579 PPARG      Down_in_ACLR           -1.33 0.00102         -2.33 0.0000336  
 4 BUB1       Bub1         ENSSSCG00000030469 BUB1       Up_in_ACLR              3.88 0.000000141      2.46 0.0000135  
 5 CCNB1      Ccnb1        ENSSSCG00000029326 CCNB1      Up_in_ACLR              3.81 0.000000125      2.07 0.0000103  
 6 CDK1       Cdk1         ENSSSCG00000010214 CDK1       Up_in_ACLR              3.18 0.00000328       2.06 0.0000756  
 7 CKAP2      Ckap2        ENSSSCG00000009378 CKAP2      Up_in_ACLR              3.34 0.000000324      2.00 0.0000158  
 8 ENTPD7     Entpd7       ENSSSCG00000010540 ENTPD7     Up_in_ACLR              1.70 0.0000369        1.50 0.0000135  
 9 ESPL1      Espl1        ENSSSCG00000000265 ESPL1      Up_in_ACLR              3.29 0.00000156       1.61 0.00163    
10 FOXM1      Foxm1        ENSSSCG00000000739 FOXM1      Up_in_ACLR              3.53 0.000000226      1.87 0.0000401  
11 HIF1A      Hif1a        ENSSSCG00000005096 HIF1A      Up_in_ACLR              1.06 0.000295         1.06 0.00126    
12 HTRA1      Htra1        ENSSSCG00000010703 HTRA1      Up_in_ACLR              1.38 0.0000203        1.27 0.0000243  
13 KIF20A     Kif20a       ENSSSCG00000014326 KIF20A     Up_in_ACLR              3.87 0.000000122      2.22 0.00000243 
14 KIF2C      Kif2c        ENSSSCG00000003931 KIF2C      Up_in_ACLR              3.31 0.000000122      1.73 0.00000496 
15 KIFC1      Kifc1        ENSSSCG00000001510 KIFC1      Up_in_ACLR              2.66 0.00000738       1.22 0.00130    
16 LOXL2      Loxl2        ENSSSCG00000033608 LOXL2      Up_in_ACL   human_gene mouse_symbol pig_ensg           pig_symbol signature_direction logFC_t7      FDR_t7 logFC_t28     FDR_t28
   <chr>      <chr>        <chr>              <chr>      <chr>                  <dbl>       <dbl>     <dbl>       <dbl>
 1 CBS        Cbs          ENSSSCG00000040779 CBS        Down_in_ACLR           -1.90 0.000653        -3.11 0.0000993  
 2 FASN       Fasn         ENSSSCG00000029944 FASN       Down_in_ACLR           -2.11 0.000494        -3.93 0.000000106
 3 PPARG      Pparg        ENSSSCG00000011579 PPARG    
Sample annotation:
   sample_id    group
1       CON1  Control
2       CON2  Control
3       CON3  Control
4       CON4  Control
5       CON5  Control
6       CON6  Control
7      INJS1  ACLT_t7
8      INJS2  ACLT_t7
9      INJS3  ACLT_t7
10     INJS4  ACLT_t7
11     INJS5  ACLT_t7
12     INJS6  ACLT_t7
13     INJL1 ACLT_t28
14     INJL2 ACLT_t28
15     INJL3 ACLT_t28
16     INJL4 ACLT_t28
17     INJL5 ACLT_t28
18     INJL6 ACLT_t28

Step21 current78 Figure4C heatmap completed successfully.
000324      2.00 0.0000158  
 8 ENTPD7     Entpd7       ENSSSCG00000010540 ENTPD7     Up_in_ACLR              1.70 0.0000369        1.50 0.0000135  
 9 ESPL1      Espl1        ENSSSCG00000000265 ESPL1      Up_in_ACLR              3.29 0.00000156       1.61 0.00163    
10 FOXM1      Foxm1        ENSSSCG00000000739 FOXM1      Up_in_ACLR              3.53 0.000000226      1.87 0.0000401  
11 HIF1A      Hif1a        ENSSSCG00000005096 HIF1A      Up_in_ACLR              1.06 0.000295         1.06 0.00126    
12 HTRA1      Htra1        ENSSSCG00000010703 HTRA1      Up_in_ACLR              1.38 0.0000203        1.27 0.0000243  
13 KIF20A     Kif20a       ENSSSCG00000014326 KIF20A     Up_in_ACLR              3.87 0.000000122      2.22 0.00000243 
14 KIF2C      Kif2c        human_gene mouse_symbol pig_ensg           pig_symbol signature_direction logFC_t7      FDR_t7 logFC_t28     FDR_t28
   <chr>      <chr>        <chr>              <chr>      <chr>                  <dbl>       <dbl>     <dbl>       <dbl>
 1 CBS        Cbs          ENSSSCG00000040779 CBS        Down_in_ACLR           -1.90 0.000653        -3.11 0.0000993  
 2 FASN       Fasn         ENSSSCG00000029944 FASN       Down_in_ACLR           -2.11 0.000494        -3.93 0.000000106
 3 PPARG      Pparg        ENSSSCG00000011579 PPARG      Down_in_ACLR           -1.33 0.00102         -2.33 0.0000336  
 4 BUB1       Bub1         ENSSSCG00000030469 BUB1       Up_in_ACLR              3.88 0.000000141      2.46 0.0000135  
 5 CCNB1      Ccnb1        ENSSSCG00000029326 CCNB1      Up_in_ACLR              3.81 0.000000125
Sample annotation:
   sample_id    group
1       CON1  Control
2       CON2  Control
3       CON3  Control
4       CON4  Control
5       CON5  Control
6       CON6  Control
7      INJS1  ACLT_t7
8      INJS2  ACLT_t7
9      INJS3  ACLT_t7
10     INJS4  ACLT_t7
11     INJS5  ACLT_t7
12     INJS6  ACLT_t7
13     INJL1 ACLT_t28
14     INJL2 ACLT_t28
15     INJL3 ACLT_t28
16     INJL4 ACLT_t28
17     INJL5 ACLT_t28
18     INJL6 ACLT_t28

Step21 current78 Figure4C heatmap completed successfully.

Step21 current78 Figure4C heatmap completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_core_heatmap/logs/step21_current78_Figure4C_core_heatmap_summary_to_send_me.txt
===== STEP21 CURRENT78 FIGURE4C VALIDATION FLOW BARPLOT =====
Summary file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_early_signature_remap_summary.csv
gorth raw file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_human_to_pig_gorth_raw.csv
Validation file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_validation_table.csv
indexing step18_current78_pig_early_signature_remap_summary.csv [=================================================================================] 35.15MB/s, eta:  0s                                                                                                                                                                       indexing step18_current78_human_to_pig_gorth_raw.csv [===========================================================================================] 918.22MB/s, eta:  0s                                                                                                                                                                       indexing step18_current78_pig_signature_validation_table.csv [=====================================================================================] 1.52GB/s, eta:  0s                                                                                                                                                                       Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_validation_flow_barplot/step21_current78_Figure4C_validation_flow_source_data.csv
null device 
          1 
null device 
          1 
[1] TRUE
[1] TRUE
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_validation_flow_barplot/step21_current78_Figure4C_validation_flow_barplot_summary.csv

===== STEP21 CURRENT78 FIGURE4C VALIDATION FLOW SUMMARY =====
                                          metric
1      current_human_CellAge_overlap_input_genes
2               human_genes_with_raw_pig_mapping
3  current_human_genes_strict_1to1_mapped_to_pig
4           current_pig_signature_genes_detected
5           direction_consistent_both_timepoints
6                        core_strict_both_t7_t28
7                                       plot_png
8                                       plot_pdf
9                             canonical_plot_png
10                            canonical_plot_pdf
11                                    output_dir
                                                                                                                    value
1                                                                                                                      78
2                                                                                                                      78
3                                                                                                                      75
4                                                                                                                      75
5                                                                                                                      47
6                                                                                                                      24
7  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_current78_pig_validation_flow_barplot.png
8  E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_current78_pig_validation_flow_barplot.pdf
9            E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_pig_validation_flow_barplot.png
10           E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4C_pig_validation_flow_barplot.pdf
11           E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_validation_flow_barplot

Plot source data:
  step_order                                                           category  n
1          6   24 strict core ortholog genes\n(t7 + t28 + direction-consistent) 24
2          5 47 direction-consistent in pig\n(t7 and t28 direction match mouse) 47
3          4                               75 detected in pig expression matrix 75
4          3                                        75 strict 1:1 pig orthologs 75
5          2                                  78 human→pig raw ortholog mapping 78
6          1                       78 human genes\n(mouse persistent ∩ CellAge) 78

Step21 current78 Figure4C validation flow barplot completed successfully.

Step21 current78 Figure4C validation flow barplot completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step21_current78_Figure4C_validation_flow_barplot/logs/step21_current78_Figure4C_validation_flow_barplot_summary_to_send_me.txt
===== STEP22 CURRENT78 FIGURE4D CORE HEATMAP REFERENCE STYLE =====
Validation file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_validation_table.csv
Step18 summary file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_early_signature_remap_summary.csv
Step19 score file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv
Count file: E:/R/ACLsenescence2 LD/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv
Step11C overlap file: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_79_genes.csv
indexing step18_current78_pig_signature_validation_table.csv [===================================================================================] 897.17MB/s, eta:  0s                                                                                                                                                                       indexing step18_current78_pig_early_signature_remap_summary.csv [=================================================================================] 58.08MB/s, eta:  0s                                                                                                                                                                       indexing step19_current78_pig_signature_scores_by_sample.csv [===================================================================================] 279.22MB/s, eta:  0s                                                                                                                                                                       indexing step11C_persistent_CellAge_overlap_79_genes.csv [=========================================================================================] 1.25GB/s, eta:  0s                                                                                                                                                                       indexing step16_pig_early_gene_count_matrix.csv [==------------------------------------------------------------------------------------------------] 2.15GB/s, eta:  0sindexing step16_pig_early_gene_count_matrix.csv [=================================================================================================] 51.03MB/s, eta:  0s                                                                                                                                                                       null device 
          1 
null device 
          1 
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_reference_style/step22_current78_Figure4D_core_ortholog_gene_table.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_reference_style/step22_current78_Figure4D_core_heatmap_zscore_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_reference_style/step22_current78_Figure4D_core_sample_annotation.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_reference_style/step22_current78_Figure4D_core_heatmap_reference_style_summary.csv

===== STEP22 CURRENT78 FIGURE4D SUMMARY =====
                metric                                                                                                                            value
1 current78_core_genes                                                                                                                               24
2            control_n                                                                                                                                6
3            aclt_1w_n                                                                                                                                6
4            aclt_4w_n                                                                                                                                6
5             png_file E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4D_current78_core_ortholog_heatmap_reference_style.png
6             pdf_file E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4D_current78_core_ortholog_heatmap_reference_style.pdf
7           output_dir                E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_reference_style
8                 note                                   This script intentionally writes only one PNG and one PDF, with no duplicate canonical copies.

Core genes preview:
# A tibble: 20 × 9
   display_label human_gene pig_ensg           CellAgeEffect MouseDirection        logFC_t7      FDR_t7 logFC_t28     FDR_t28
   <chr>         <chr>      <chr>              <chr>         <chr>                    <dbl>       <dbl>     <dbl>       <dbl>
 1 FASN          FASN       ENSSSCG00000029944 Induces       mouse_persistent_down    -2.11 0.000494        -3.93 0.000000106
 2 PPARG         PPARG      ENSSSCG00000011579 Induces       mouse_persistent_down    -1.33 0.00102         -2.33 0.0000336  
 3 CBS           CBS        ENSSSCG00000040779 Inhibits      mouse_persistent_down    -1.90 0.000653        -3.11 0.0000993  
 4 ENTPD7        ENTPD7     ENSSSCG00000010540 Induces       mouse_persistent_up       1.70 0.0000369        1.50 0.0000135  
 5 HTRA1         HTRA1      ENSSSCG00000010703 Induces       mouse_persistent_up       1.38 0.0000203        1.27 0.0000243  
 6 RUNX1         RUNX1      ENSSSCG00000035537 Induces       mouse_persistent_up       1.44 0.000000141      1.28 0.0000478  
 7 SERPINE1      SERPINE1   ENSSSCG00000025698 Induces       mouse_persistent_up       1.26 0.00725          1.38 0.00573    
 8 SPARC         SPARC      ENSSSCG00000017082 Induces       mouse_persistent_up       1.21 0.00721          1.58 0.00323    
 9 TNFSF15       TNFSF15    ENSSSCG00000059646 Induces       mouse_persistent_up       1.84 0.00302          1.44 0.00933    
10 BUB1          BUB1       ENSSSCG00000030469 Inhibits      mouse_persistent_up       3.88 0.000000141      2.46 0.0000135  
11 CCNB1         CCNB1      ENSSSCG00000029326 Inhibits      mouse_persistent_up       3.81 0.000000125      2.07 0.0000103  
12 CDK1          CDK1       ENSSSCG00000010214 Inhibits      mouse_persistent_up       3.18 0.00000328       2.06 0.0000756  
13 CKAP2         CKAP2      ENSSSCG00000009378 Inhibits         display_label human_gene pig_ensg           CellAgeEffect MouseDirection        logFC_t7      FDR_t7 logFC_t28     FDR_t28
   <chr>         <chr>      <chr>              <chr>         <chr>                    <dbl>       <dbl>     <dbl>       <dbl>
 1 FASN          FASN       ENSSSCG00000029944 Induces       mouse_persistent_down    -2.11 0.000494        -3.93 0.000000106
 2 PPARG         PPARG      ENSSSCG00000011579 Induces       mouse_persistent_down    -1.33 0.00102         -2.33 0.0000336  
 3 CBS           CBS        ENSSSCG00000040779 Inhibits      mouse_persistent_down    -1.90 0.000653        -3.11 0.0000993  
 4 ENTPD7        ENTPD7     ENSSSCG00000010540 Induces       mouse_persistent_up       1.70 0.0000369        1.50 0.0000135  
 5 HTRA1         HTRA1      ENSSSCG00000010703 Induces       mouse_persistent_up       1.38 0.0000203        1.27 0.0000243  
 6 RUNX1         RUNX1      ENSSSCG00000035537 Induces       mouse
Sample annotation:
# A tibble: 18 × 2
   sample_id group            
   <chr>     <fct>            
 1 CON1      Control          
 2 CON2      Control          
 3 CON3      Control          
 4 CON4      Control          
 5 CON5      Control          
 6 CON6      Control          
 7 INJS1     ACLT-untreated-1W
 8 INJS2     ACLT-untreated-1W
 9 INJS3     ACLT-untreated-1W
10 INJS4     ACLT-untreated-1W
11 INJS5     ACLT-untreated-1W
12 INJS6     ACLT-untreated-1W
13 INJL1     ACLT-untreated-4W
14 INJL2     ACLT-untreated-4W
15 INJL3     ACLT-untreated-4W
16 INJL4     ACLT-untreated-4W
17 INJL5     ACLT-untreated-4W
18 INJL6     ACLT-untreated-4W

Step22 current78 Figure4D core heatmap completed successfully.

Step22 current78 Figure4D core heatmap completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_reference_style/logs/step22_current78_Figure4D_core_heatmap_reference_style_summary_to_send_me.txt
===== STEP22 CURRENT78 FIGURE4D CORE HEATMAP STEP09 STYLE =====
Run time: 2026-05-06 21:35:44 CEST
Style: rev(RColorBrewer::RdYlBu), z-score clipped to [-2, 2], expanded canvas.

Validation file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_validation_table.csv
Step18 summary file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_early_signature_remap_summary.csv
Step19 score file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv
Count file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv
Step11C overlap file: E:/R/ACLsenescence2/rebuild_submission/02_mouse_discovery/07_tables/step11C_intersect_current_gorth_mapping_with_CellAge_clean/step11C_persistent_CellAge_overlap_79_genes.csv

indexing step18_current78_pig_signature_validation_table.csv [=====================================================================================] 1.16GB/s, eta:  0s                                                                                                                                                                       indexing step18_current78_pig_early_signature_remap_summary.csv [=================================================================================] 63.61MB/s, eta:  0s                                                                                                                                                                       indexing step19_current78_pig_signature_scores_by_sample.csv [===================================================================================] 279.22MB/s, eta:  0s                                                                                                                                                                       indexing step11C_persistent_CellAge_overlap_79_genes.csv [=========================================================================================] 1.61GB/s, eta:  0s                                                                                                                                                                       indexing step16_pig_early_gene_count_matrix.csv [=-------------------------------------------------------------------------------------------------] 1.50GB/s, eta:  0sindexing step16_pig_early_gene_count_matrix.csv [===========================================================================================------] 11.79MB/s, eta:  0sindexing step16_pig_early_gene_count_matrix.csv [=================================================================================================] 12.38MB/s, eta:  0s                                                                                                                                                                       Core genes identified: 24
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_step09_style/step22_current78_Figure4D_core_ortholog_gene_table.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_step09_style/step22_current78_Figure4D_core_heatmap_zscore_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_step09_style/step22_current78_Figure4D_core_sample_annotation.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_step09_style/step22_current78_Figure4D_core_heatmap_step09_style_summary.csv

===== STEP22 CURRENT78 FIGURE4D SUMMARY =====
                 metric                                                                                                                         value
1  current78_core_genes                                                                                                                            24
2             control_n                                                                                                                             6
3             aclt_1w_n                                                                                                                             6
4             aclt_4w_n                                                                                                                             6
5   heatmap_color_scale                                                                   rev(RColorBrewer::RdYlBu): deep blue -> lemon yellow -> red
6           zscore_clip                                                                                                                       [-2, 2]
7              pdf_file E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4D_current78_core_ortholog_heatmap_step09_style.pdf
8              png_file E:/R/ACLsenescence2/rebuild_submission/02_pig_early/figures/Figure4/Figure4D_current78_core_ortholog_heatmap_step09_style.png
9            output_dir                E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_step09_style
10                 note                                                This script writes exactly one PDF and one PNG. No duplicate canonical copies.

Core genes preview:
# A tibble: 20 × 9
   display_label human_gene pig_ensg           CellAgeEffect MouseDirection        logFC_t7      FDR_t7 logFC_t28     FDR_t28
   <chr>         <chr>      <chr>              <chr>         <chr>                    <dbl>       <dbl>     <dbl>       <dbl>
 1 FASN          FASN       ENSSSCG00000029944 Induces       mouse_persistent_down    -2.11 0.000494        -3.93 0.000000106
 2 PPARG         PPARG      ENSSSCG00000011579 Induces       mouse_persistent_down    -1.33 0.00102         -2.33 0.0000336  
 3 CBS           CBS        ENSSSCG00000040779 Inhibits      mouse_persistent_down    -1.90 0.000653        -3.11 0.0000993  
 4 ENTPD7        ENTPD7     ENSSSCG00000010540 Induces       mouse_persistent_up       1.70 0.0000369        1.50 0.0000135  
 5 HTRA1         HTRA1      ENSSSCG00000010703 Induces       mouse_persistent_up       1.38 0.0000203        1.27 0.0000243  
 6 RUNX1         RUNX1      ENSSSCG00000035537 Induces       mouse_persistent_up       1.44 0.000000141      1.28 0.0000478  
 7 SERPINE1      SERPINE1   ENSSSCG00000025698 Induces       mouse_persistent_up       1.26 0.00725          1.38 0.00573    
 8 SPARC         SPARC      ENSSSCG00000017082 Induces       mouse_persistent_up       1.21 0.00721          1.58 0.00323    
 9 TNFSF15       TNFSF15    ENSSSCG00000059646 Induces       mouse_persistent_up       1.84 0.00302          1.44 0.00933    
10 BUB1          BUB1       ENSSSCG00000030469 Inhibits      mouse_persistent_up       3.88 0.000000141      2.46 0.0000135  
11 CCNB1         CCNB1      ENSSSCG00000029326 Inhibits      mouse_persistent_up       3.81 0.000000125      2.07 0.0000103  
12 CDK1          CDK1       ENSSSCG00000010214 Inhibits      mouse_persistent_up       3.18 0.00000328       2.06 0.0000756  
13 CKAP2         CKAP2      ENSSSCG00000009378 Inhibits         display_label human_gene pig_ensg           CellAgeEffect MouseDirection        logFC_t7      FDR_t7 logFC_t28     FDR_t28
   <chr>         <chr>      <chr>              <chr>         <chr>                    <dbl>       <dbl>     <dbl>       <dbl>
 1 FASN          FASN       ENSSSCG00000029944 Induces       mouse_persistent_down    -2.11 0.000494        -3.93 0.000000106
 2 PPARG         PPARG      ENSSSCG00000011579 Induces       mouse_persistent_down    -1.33 0.00102         -2.33 0.0000336  
 3 CBS           CBS        ENSSSCG00000040779 Inhibits      mouse_persistent_down    -1.90 0.000653        -3.11 0.0000993  
 4 ENTPD7        ENTPD7     ENSSSCG00000010540 Induces       mouse_persistent_up       1.70 0.0000369        1.50 0.0000135  
 5 HTRA1         HTRA1      ENSSSCG00000010703 Induces       mouse_persistent_up       1.38 0.0000203        1.27 0.0000243  
 6 RUNX1         RUNX1      ENSSSCG00000035537 Induces       mouse
Sample annotation:
# A tibble: 18 × 2
   sample_id group            
   <chr>     <fct>            
 1 CON1      Control          
 2 CON2      Control          
 3 CON3      Control          
 4 CON4      Control          
 5 CON5      Control          
 6 CON6      Control          
 7 INJS1     ACLT-untreated-1W
 8 INJS2     ACLT-untreated-1W
 9 INJS3     ACLT-untreated-1W
10 INJS4     ACLT-untreated-1W
11 INJS5     ACLT-untreated-1W
12 INJS6     ACLT-untreated-1W
13 INJL1     ACLT-untreated-4W
14 INJL2     ACLT-untreated-4W
15 INJL3     ACLT-untreated-4W
16 INJL4     ACLT-untreated-4W
17 INJL5     ACLT-untreated-4W
18 INJL6     ACLT-untreated-4W

Step22 current78 Figure4D core heatmap completed successfully.

Step22 current78 Figure4D core heatmap completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step22_current78_Figure4D_core_heatmap_step09_style/logs/step22_current78_Figure4D_core_heatmap_step09_style_summary_to_send_me.txt
===== STEP23 CURRENT78 SINGSCORE SENSITIVITY ANALYSIS =====
Run time: 2026-05-06 21:55:54 CEST

Signature file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
Count matrix file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv
Step19 score annotation file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv

indexing step18_current78_pig_signature_gene_table.csv [=========================================================================================] 763.53MB/s, eta:  0s                                                                                                                                                                       indexing step16_pig_early_gene_count_matrix.csv [=============================================================-------------------------------------] 2.15GB/s, eta:  0sindexing step16_pig_early_gene_count_matrix.csv [================================================================================================] 318.73MB/s, eta:  0s                                                                                                                                                                       indexing step19_current78_pig_signature_scores_by_sample.csv [===================================================================================] 133.26MB/s, eta:  0s                                                                                                                                                                       Signature genes: 75
Up genes: 65
Down genes: 10

Genome-wide logCPM matrix dim:
[1] 35682    18
===== STEP23 CURRENT78 SINGSCORE SENSITIVITY ANALYSIS =====
Run time: 2026-05-06 21:58:03 CEST

Signature file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
Count matrix file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv
Step19 score annotation file: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv

indexing step18_current78_pig_signature_gene_table.csv [===========================================================================================] 1.92GB/s, eta:  0s                                                                                                                                                                       indexing step16_pig_early_gene_count_matrix.csv [================================================================================================--] 2.15GB/s, eta:  0sindexing step16_pig_early_gene_count_matrix.csv [==================================================================================================] 2.15GB/s, eta:  0s                                                                                                                                                                       indexing step19_current78_pig_signature_scores_by_sample.csv [===================================================================================] 234.55MB/s, eta:  0s                                                                                                                                                                       Signature genes: 75
Up genes: 65
Down genes: 10

Genome-wide logCPM matrix dim:
[1] 35682    18
null device 
          1 
null device 
          1 
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_signature_gene_metadata.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_scores_by_sample.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_scores_long.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_group_comparisons.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_signature_gene_logCPM_matrix.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_versions_and_parameters.csv
Saved: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_sensitivity_summary.csv

===== STEP23_CURRENT78 SINGSCORE SUMMARY =====
                      metric
1       signature_input_file
2          count_matrix_file
3     step19_annotation_file
4          n_signature_genes
5       n_up_signature_genes
6     n_down_signature_genes
7  n_genomewide_genes_ranked
8                  n_samples
9                  n_Control
10                 n_ACLT_t7
11                n_ACLT_t28
12          singscore_method
13    simpleScore_parameters
14     main_scores_by_sample
15    group_comparison_table
16  supplementary_figure_png
17  supplementary_figure_pdf
18                output_dir
                                                                                                                                                                                value
1                                 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step18_current78_pig_early_signature_remap/step18_current78_pig_signature_gene_table.csv
2                                                                                   E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step16_pig_early_gene_count_matrix.csv
3                                 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/tables/step19_current78_pig_signature_score/step19_current78_pig_signature_scores_by_sample.csv
4                                                                                                                                                                                  75
5                                                                                                                                                                                  65
6                                                                                                                                                                                  10
7                                                                                                                                                                               35682
8                                                                                                                                                                                  18
9                                                                                                                                                                                   6
10                                                                                                                                                                                  6
11                                                                                                                                                                                  6
12                                                                                                        rank-based single-sample scoring using genome-wide within-sample gene ranks
13                                                                                                                rankGenes(); simpleScore(centerScore = TRUE, knownDirection = TRUE)
14                     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_scores_by_sample.csv
15                    E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/step23_current78_singscore_group_comparisons.csv
16 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/figures/Supplementary_Figure_S9_current78_singscore_sensitivity.png
17 E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/figures/Supplementary_Figure_S9_current78_singscore_sensitivity.pdf
18                                                                     E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity

Singscore group comparison statistics:
           comparison                     score n_control n_case median_control median_case mean_control  mean_case median_difference_case_minus_control
1  ACLT_t7_vs_Control     directional_singscore         6      6     0.04771076   0.1139393   0.03959381  0.1196749                           0.06622853
2 ACLT_t28_vs_Control     directional_singscore         6      6     0.04771076   0.1165393   0.03959381  0.1389831                           0.06882856
3  ACLT_t7_vs_Control              up_singscore         6      6     0.24270843   0.3060650   0.24402399  0.3045048                           0.06335652
4 ACLT_t28_vs_Control              up_singscore         6      6     0.24270843   0.2812328   0.24402399  0.2796948                           0.03852439
5  ACLT_t7_vs_Control down_singscore_reoriented         6      6    -0.20725219  -0.1904575  -0.20443018 -0.1848298                           0.01679468
6 ACLT_t28_vs_Control down_singscore_reoriented         6      6    -0.20725219  -0.1598635  -0.20443018 -0.1407117                           0.04738871
7  ACLT_t7_vs_Control        down_singscore_raw         6      6     0.20725219   0.1904575   0.20443018  0.1848298                          -0.01679468
8 ACLT_t28_vs_Control        down_singscore_raw         6      6     0.20725219   0.1598635   0.20443018  0.1407117                          -0.04738871
  mean_difference_case_minus_control wilcox_p_value BH_FDR_across_all_singscore_comparisons
1                         0.08008111    0.005074868                              0.01014974
2                         0.09938932    0.005074868                              0.01014974
3                         0.06048077    0.005074868                              0.01014974
4                         0.03567081    0.005074868                              0.01014974
5                         0.01960034    0.297953062                              0.29795306
6                         0.06371851    0.092695803                              0.12359440
7                        -0.01960034    0.297953062                              0.29795306
8                        -0.06371851    0.092695803                              0.12359440

Version and parameter records:
                       item                                                               value
1                 R_version                                                               4.5.2
2 singscore_package_version                                                              1.30.0
3     edgeR_package_version                                                               4.8.2
4   ggplot2_package_version                                                               4.0.2
5     dplyr_package_version                                                               1.1.4
6  ranked_expression_matrix   genome-wide TMM-normalized logCPM matrix from Step16 count matrix
7           rankGenes_input                                 singscore::rankGenes(logCPM_matrix)
8    simpleScore_parameters              simpleScore(centerScore = TRUE, knownDirection = TRUE)
9                group_test two-sided Wilcoxon rank-sum test for ACLT_t7 or ACLT_t28 vs Control

Step23 current78 singscore sensitivity analysis completed successfully.

Step23 current78 singscore sensitivity completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/logs/step23_current78_singscore_sensitivity_summary_to_send_me.txt
8                              0.01014974
4                         0.03567081    0.005074868                              0.01014974
5                         0.01960034    0.297953062                              0.29795306
6                         0.06371851    0.092695803                              0.12359440
7                        -0.01960034    0.297953062                              0.29795306
8                        -0.06371851    0.092695803                              0.12359440

Version and parameter records:
                       item                                                               value
1                 R_version                                                               4.5.2
2 singscore_package_version                                                              1.30.0
3     edgeR_package_version                                                               4.8.2
4   ggplot2_package_version                                                               4.0.2
5     dplyr_package_version                                                               1.1.4
6  ranked_expression_matrix   genome-wide TMM-normalized logCPM matrix from Step16 count matrix
7           rankGenes_input                                 singscore::rankGenes(logCPM_matrix)
8    simpleScore_parameters              simpleScore(centerScore = TRUE, knownDirection = TRUE)
9                group_test two-sided Wilcoxon rank-sum test for ACLT_t7 or ACLT_t28 vs Control

Step23 current78 singscore sensitivity analysis completed successfully.

Step23 current78 singscore sensitivity completed.
Summary to send me: E:/R/ACLsenescence2/rebuild_submission/02_pig_early/supplementary files/step23_current78_singscore_sensitivity/logs/step23_current78_singscore_sensitivity_summary_to_send_me.txt
