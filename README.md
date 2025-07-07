
### Datasets and code complementary to research study:

Title: *Variation in metabolic scaling reveals aerobic constraints with temperature in fishes."*

Authors: Krista Kraskura\*, Christopher L. Jerde, Erika J. Eliason

Contact Info: Krista Kraskura, kkraskura\@towson.edu

------------------------------------------------------------------------

### About the repository:

#### /Data

*Files:*

-   /Ecology groupings.xlsx
    -   *the details for species ecology and morphology group
        assignemnt.*
-   /Kraskura_species_ecologies_mar2022.csv
    -   *the finalized dataset with species assigned ecology and
        morphology group assignments.*
-   /Fish_AMR_temp_dataset_mar2022.csv
    -   *the finalized dataset with Active and Resting Metabolic Rate
        data (AMR, RMR) to obtain aerobic scopes*
-   /Fish_RMR_temp_dataset_mar2022.csv
    -   *the finalized dataset with Resting Metabolic Rate data (RMR)*
-   /Summary_scaling_b_mmr_rmr_jul192022.csv
    -   *the dataset for Figure 1 in the main text. collated dataset
        with literature-reported scaling relationships of RMR and MMR in
        fishes.*

#### /R

***Execution:** files "get_data_temp.R",* *"get_data_phylo_matrix.R",
"mixed_model_outputs.R" contain functions that are used in other
scripts. The "colors_themes.R" is called in in each file to obtain the
same dataset, colors schemes, and libraries. Each R script has a
description on the top rows.*

*Files:*

-   /get_data_temp.R

-   /get_data_phylo_matrix.R

-   /mixed_model_outputs.R

-   /colors_themes.R

<!-- -->

-   /phylo_mixed_model.R
-   /non_scaling_figures.R
-   /species_ecologies_level_model_plots.R
-   /phylo_model_update_wEcol.R

<!-- -->

-   /data_tables_summaries_ms_export.Rmd

-   /figure1_scaling_exp_sum_mmr_rmr.R

#### /Data_exports

Exported figures from various R scripts.

#### /Figures

Exported figures from various R scripts.

------------------------------------------------------------------------

**Session Info:**
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 15.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils    
[5] datasets  methods   base     

other attached packages:
 [1] flextable_0.9.8      systemfonts_1.2.3   
 [3] reshape2_1.4.4       cowplot_1.1.3       
 [5] ggformat2_0.1.0      ggimage_0.3.3       
 [7] TDbook_0.0.6         viridisLite_0.4.2   
 [9] pryr_0.1.6           ggtree_3.10.1       
[11] rotl_3.1.0           ape_5.8             
[13] evolvability_2.0.0   car_3.1-2           
[15] carData_3.0-5        emmeans_1.10.1      
[17] lme4_1.1-35.3        Matrix_1.6-3        
[19] weathermetrics_1.2.2 patchwork_1.2.0     
[21] colorBlindness_0.1.9 lubridate_1.9.3     
[23] forcats_1.0.0        stringr_1.5.1       
[25] dplyr_1.1.4          purrr_1.0.4         
[27] readr_2.1.5          tidyr_1.3.1         
[29] tibble_3.2.1         tidyverse_2.0.0     
[31] ggh4x_0.3.1          here_1.0.1          
[33] ggpubr_0.6.0         ggplot2_3.5.2       

loaded via a namespace (and not attached):
 [1] sandwich_3.1-0         
 [2] rlang_1.1.6            
 [3] magrittr_2.0.3         
 [4] multcomp_1.4-25        
 [5] compiler_4.3.2         
 [6] vctrs_0.6.5            
 [7] pkgconfig_2.0.3        
 [8] crayon_1.5.3           
 [9] fastmap_1.1.1          
[10] ellipsis_0.3.2         
[11] backports_1.4.1        
[12] magick_2.8.3           
[13] labeling_0.4.3         
[14] rmarkdown_2.26         
[15] tzdb_0.5.0             
[16] nloptr_2.0.3           
[17] ragg_1.3.1             
[18] xfun_0.43              
[19] cachem_1.0.8           
[20] aplot_0.2.2            
[21] jsonlite_2.0.0         
[22] progress_1.2.3         
[23] uuid_1.2-0             
[24] broom_1.0.5            
[25] parallel_4.3.2         
[26] prettyunits_1.2.0      
[27] R6_2.6.1               
[28] stringi_1.8.7          
[29] RColorBrewer_1.1-3     
[30] boot_1.3-30            
[31] estimability_1.5.1     
[32] knitr_1.45             
[33] Rcpp_1.0.14            
[34] zoo_1.8-12             
[35] rentrez_1.2.3          
[36] splines_4.3.2          
[37] timechange_0.3.0       
[38] tidyselect_1.2.1       
[39] rstudioapi_0.17.1      
[40] dichromat_2.0-0.1      
[41] abind_1.4-5            
[42] codetools_0.2-20       
[43] lattice_0.22-6         
[44] plyr_1.8.9             
[45] treeio_1.26.0          
[46] withr_3.0.2            
[47] askpass_1.2.1          
[48] evaluate_0.23          
[49] coda_0.19-4.1          
[50] gridGraphics_0.5-1     
[51] survival_3.6-4         
[52] zip_2.3.1              
[53] xml2_1.3.6             
[54] pillar_1.10.2          
[55] rsconnect_1.2.2        
[56] ggfun_0.1.4            
[57] generics_0.1.3         
[58] rprojroot_2.0.4        
[59] hms_1.1.3              
[60] scales_1.4.0           
[61] tidytree_0.4.6         
[62] minqa_1.2.6            
[63] xtable_1.8-4           
[64] rncl_0.8.7             
[65] glue_1.8.0             
[66] gdtools_0.4.2          
[67] lazyeval_0.2.2         
[68] tools_4.3.2            
[69] data.table_1.17.0      
[70] ggsignif_0.6.4         
[71] fs_1.6.6               
[72] mvtnorm_1.3-3          
[73] XML_3.99-0.16.1        
[74] grid_4.3.2             
[75] colorspace_2.1-1       
[76] nlme_3.1-164           
[77] cli_3.6.5              
[78] textshaping_0.3.7      
[79] officer_0.6.9          
[80] fontBitstreamVera_0.1.1
[81] gtable_0.3.6           
[82] rstatix_0.7.2          
[83] yulab.utils_0.1.4      
[84] fontquiver_0.2.1       
[85] digest_0.6.35          
[86] ggplotify_0.1.2        
[87] TH.data_1.1-2          
[88] farver_2.1.2           
[89] htmltools_0.5.8.1      
[90] memoise_2.0.1          
[91] lifecycle_1.0.4        
[92] httr_1.4.7             
[93] openssl_2.3.2          
[94] fontLiberation_0.1.0   
[95] MASS_7.3-60  
