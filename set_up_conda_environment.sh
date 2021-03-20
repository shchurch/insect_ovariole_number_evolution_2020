### This code was written by SHC in Jan 2021
### This gives the commands used to build the R environment from which
###    all analyses were performed

# build an R environment with base R
conda create -n r-4.0.3 -c conda-forge r-base=4.0.3

# activate the environment
conda activate r-4.0.3

# install R packages
conda install -c conda-forge r-devtools
R
install.packages("ape")
install.packages("dplyr")
install.packages("plyr")
install.packages("nlme")
install.packages("phytools")
install.packages("ggplot2")
install.packages("BAMMtools")
install.packages("coda")
install.packages("OUwie")
# on cluster, outside of R
#conda install -c conda-forge r-rmpfr
#conda install -c anaconda gmp
install.packages("corHMM")
install.packages("grid")
install.packages("gridExtra")
install.packages("gtable")
install.packages("geiger")
install.packages("RColorBrewer")
install.packages("rmarkdown")
install.packages("bookdown")
install.packages("kableExtra")
install.packages("ggplotify")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")

devtools::install_github("mwpennell/arbutus")

sessionInfo()
#	R version 4.0.3 (2020-10-10)
#	Platform: x86_64-apple-darwin13.4.0 (64-bit)
#	Running under: macOS Catalina 10.15.7
#	
#	Matrix products: default
#	BLAS/LAPACK: /opt/anaconda3/envs/r-4.0.3/lib/libopenblasp-r0.3.12.dylib
#	
#	locale:
#	[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#	
#	attached base packages:
#	[1] grid      stats     graphics  grDevices utils     datasets  methods  
#	[8] base     
#	
#	other attached packages:
#	 [1] arbutus_0.1        ggtree_2.4.1       ggplotify_0.0.5    kableExtra_1.3.1  
#	 [5] bookdown_0.21      rmarkdown_2.6      gtable_0.3.0       gridExtra_2.3     
#	 [9] corHMM_2.6         GenSA_1.1.7        OUwie_2.5          RColorBrewer_1.1-2
#	[13] geiger_2.0.7       nloptr_1.2.2.2     corpcor_1.6.9      coda_0.19-4       
#	[17] BAMMtools_2.1.7    ggplot2_3.3.3      phytools_0.7-70    maps_3.3.0        
#	[21] nlme_3.1-151       plyr_1.8.6         dplyr_1.0.3        ape_5.4-1         
#	
#	loaded via a namespace (and not attached):
#	 [1] subplex_1.6             bitops_1.0-6            webshot_0.5.2          
#	 [4] httr_1.4.2              numDeriv_2016.8-1.1     tools_4.0.3            
#	 [7] R6_2.5.0                KernSmooth_2.23-18      lazyeval_0.2.2         
#	[10] colorspace_2.0-0        nnet_7.3-14             withr_2.4.0            
#	[13] tidyselect_1.1.0        mnormt_2.0.2            phangorn_2.5.5         
#	[16] compiler_4.0.3          rvest_0.3.6             expm_0.999-6           
#	[19] xml2_1.3.2              caTools_1.18.1          scales_1.1.1           
#	[22] mvtnorm_1.1-1           quadprog_1.5-8          stringr_1.4.0          
#	[25] digest_0.6.27           pkgconfig_2.0.3         htmltools_0.5.1        
#	[28] parallelly_1.23.0       plotrix_3.7-8           lhs_1.1.1              
#	[31] rlang_0.4.10            rstudioapi_0.13         phylolm_2.6.2          
#	[34] gridGraphics_0.5-1      generics_0.1.0          combinat_0.0-8         
#	[37] jsonlite_1.7.2          gtools_3.8.2            RCurl_1.98-1.2         
#	[40] magrittr_2.0.1          patchwork_1.1.1         interp_1.0-33          
#	[43] Matrix_1.3-2            Rcpp_1.0.6              munsell_0.5.0          
#	[46] viridis_0.5.1           lifecycle_0.2.0         scatterplot3d_0.3-41   
#	[49] stringi_1.5.3           clusterGeneration_1.3.7 MASS_7.3-53            
#	[52] gplots_3.1.1            parallel_4.0.3          listenv_0.8.0          
#	[55] crayon_1.3.4            deldir_0.2-9            lattice_0.20-41        
#	[58] tmvnsim_1.0-2           knitr_1.30              pillar_1.4.7           
#	[61] igraph_1.2.6            future.apply_1.7.0      codetools_0.2-18       
#	[64] paleotree_3.3.25        fastmatch_1.1-0         glue_1.4.2             
#	[67] evaluate_0.14           BiocManager_1.30.10     deSolve_1.28           
#	[70] treeio_1.14.3           png_0.1-7               vctrs_0.3.6            
#	[73] tidyr_1.1.2             purrr_0.3.4             future_1.21.0          
#	[76] xfun_0.20               Rmpfr_0.8-2             tidytree_0.3.3         
#	[79] viridisLite_0.3.0       tibble_3.0.5            aplot_0.0.6            
#	[82] rvcheck_0.1.8           gmp_0.6-2               globals_0.14.0         
#	[85] ellipsis_0.3.1
