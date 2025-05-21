# KK_etal_synthesis_fish_temp_scaling

Last fully executed (date): February 21, 2024

Files run: 

- non_scaling_figures.R
- phylo_mixed_models.R

________________________________________
Session info:

R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.1.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggimage_0.3.3        TDbook_0.0.6         pryr_0.1.6           ggtree_3.10.0       
 [5] rotl_3.1.0           ape_5.7-1            evolvability_2.0.0   car_3.1-2           
 [9] carData_3.0-5        emmeans_1.8.9        lme4_1.1-35.1        Matrix_1.6-3        
[13] lubridate_1.9.3      stringr_1.5.1        dplyr_1.1.4          purrr_1.0.2         
[17] readr_2.1.4          tidyr_1.3.1          tibble_3.2.1         tidyverse_2.0.0     
[21] ggformat2_0.1.0      weathermetrics_1.2.2 here_1.0.1           forcats_1.0.0       
[25] cowplot_1.1.1        ggpubr_0.6.0         ggplot2_3.4.4       

loaded via a namespace (and not attached):
 [1] sandwich_3.0-2     rlang_1.1.3        magrittr_2.0.3     multcomp_1.4-25    compiler_4.3.2    
 [6] systemfonts_1.0.5  vctrs_0.6.5        pkgconfig_2.0.3    crayon_1.5.2       fastmap_1.1.1     
[11] magick_2.8.1       backports_1.4.1    labeling_0.4.3     utf8_1.2.4         rmarkdown_2.25    
[16] tzdb_0.4.0         nloptr_2.0.3       ragg_1.2.6         xfun_0.41          cachem_1.0.8      
[21] aplot_0.2.2        jsonlite_1.8.8     progress_1.2.3     broom_1.0.5        parallel_4.3.2    
[26] prettyunits_1.2.0  R6_2.5.1           stringi_1.8.3      boot_1.3-28.1      estimability_1.4.1
[31] Rcpp_1.0.12        knitr_1.45         zoo_1.8-12         rentrez_1.2.3      splines_4.3.2     
[36] timechange_0.2.0   tidyselect_1.2.0   rstudioapi_0.15.0  abind_1.4-5        yaml_2.3.7        
[41] codetools_0.2-19   curl_5.2.0         lattice_0.22-5     treeio_1.26.0      withr_3.0.0       
[46] coda_0.19-4        evaluate_0.23      gridGraphics_0.5-1 survival_3.5-7     pillar_1.9.0      
[51] ggfun_0.1.3        generics_0.1.3     rprojroot_2.0.4    hms_1.1.3          munsell_0.5.0     
[56] scales_1.3.0       tidytree_0.4.5     minqa_1.2.6        xtable_1.8-4       rncl_0.8.7        
[61] glue_1.7.0         lazyeval_0.2.2     tools_4.3.2        ggsignif_0.6.4     fs_1.6.3          
[66] mvtnorm_1.2-4      XML_3.99-0.15      grid_4.3.2         colorspace_2.1-0   patchwork_1.1.3   
[71] nlme_3.1-164       cli_3.6.2          textshaping_0.3.7  fansi_1.0.6        gtable_0.3.4      
[76] rstatix_0.7.2      yulab.utils_0.1.0  digest_0.6.33      pbkrtest_0.5.2     ggplotify_0.1.2   
[81] TH.data_1.1-2      farver_2.1.1       memoise_2.0.1      htmltools_0.5.7    lifecycle_1.0.4   
[86] httr_1.4.7         MASS_7.3-60 