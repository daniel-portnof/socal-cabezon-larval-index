# socal-cabezon-larval-index
Using sdmTMB to derive a larval abundance index from CalCOFI data to be compared with existing federal stock metrics for cabezon in Southern California.

R version 4.4.3 (2025-02-28 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] maps_3.4.3              rnaturalearthdata_1.0.0 rnaturalearth_1.1.0     sf_1.0-23              
 [5] sdmTMB_1.0.0            here_1.0.2              r4ss_1.44.0             lubridate_1.9.4        
 [9] forcats_1.0.0           stringr_1.6.0           dplyr_1.1.4             purrr_1.0.4            
[13] readr_2.1.6             tidyr_1.3.1             tibble_3.2.1            ggplot2_4.0.1          
[17] tidyverse_2.0.0        

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       TMB_1.9.19         xfun_0.51          lattice_0.22-6     tzdb_0.5.0        
 [6] Rdpack_2.6.6       vctrs_0.6.5        tools_4.4.3        generics_0.1.3     proxy_0.4-27      
[11] pkgconfig_2.0.3    Matrix_1.7-2       KernSmooth_2.23-26 RColorBrewer_1.1-3 S7_0.2.1          
[16] assertthat_0.2.1   lifecycle_1.0.4    fmesher_0.6.1      compiler_4.4.3     farver_2.1.2      
[21] textshaping_1.0.0  htmltools_0.5.8.1  class_7.3-23       yaml_2.3.10        crayon_1.5.3      
[26] pillar_1.10.1      rsconnect_1.6.1    classInt_0.4-11    reformulas_0.4.4   nlme_3.1-167      
[31] tidyselect_1.2.1   digest_0.6.37      stringi_1.8.7      labeling_0.4.3     splines_4.4.3     
[36] rprojroot_2.1.1    fastmap_1.2.0      grid_4.4.3         cli_3.6.4          magrittr_2.0.3    
[41] e1071_1.7-16       corpcor_1.6.10     withr_3.0.2        scales_1.4.0       timechange_0.3.0  
[46] rmarkdown_2.29     hms_1.1.3          kableExtra_1.4.0   coda_0.19-4.1      evaluate_1.0.3    
[51] knitr_1.50         rbibutils_2.4.1    viridisLite_0.4.2  mgcv_1.9-1         rlang_1.1.5       
[56] Rcpp_1.0.14        glue_1.8.0         DBI_1.2.3          xml2_1.3.8         svglite_2.2.2     
[61] rstudioapi_0.17.1  R6_2.6.1           systemfonts_1.3.2  units_1.0-0       
