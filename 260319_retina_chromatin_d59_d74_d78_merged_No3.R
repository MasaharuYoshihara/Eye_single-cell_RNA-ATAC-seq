# GSE183684

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Signac)
library(GenomicRanges)

ATAC_merged <- readRDS("260301_retina_chromatin_merged_No2.rds")

set.seed(1)

ATAC_merged <- RunTFIDF(ATAC_merged)
ATAC_merged <- FindTopFeatures(ATAC_merged, min.cutoff = 'q0')
ATAC_merged <- RunSVD(ATAC_merged)

DepthCor(ATAC_merged)

ATAC_merged <- RunUMAP(object = ATAC_merged, reduction = 'lsi', dims = 2:30)
ATAC_merged <- FindNeighbors(object = ATAC_merged, reduction = 'lsi', dims = 2:30)
ATAC_merged <- FindClusters(object = ATAC_merged, verbose = FALSE, algorithm = 3)
DimPlot(object = ATAC_merged, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(ATAC_merged)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
ATAC_merged[['RNA']] <- CreateAssayObject(counts = gene.activities)
ATAC_merged <- NormalizeData(
  object = ATAC_merged,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC_merged$nCount_RNA)
)

DefaultAssay(ATAC_merged) <- 'RNA'

saveRDS(ATAC_merged, "260301_retina_chromatin_ATAC_merged_No3.rds")

sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8   
# [3] LC_MONETARY=Japanese_Japan.utf8 LC_NUMERIC=C                   
# [5] LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#   [1] future_1.69.0        GenomicRanges_1.58.0 GenomeInfoDb_1.42.3 
# [4] IRanges_2.40.1       S4Vectors_0.44.0     BiocGenerics_0.52.0 
# [7] Signac_1.16.0        ggplot2_4.0.0        patchwork_1.3.2     
# [10] Seurat_5.3.0         SeuratObject_5.0.2   sp_2.2-0            
# [13] dplyr_1.1.4         
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.18.0       jsonlite_1.9.0         
# [4] magrittr_2.0.3          spatstat.utils_3.1-2    farver_2.1.2           
# [7] zlibbioc_1.52.0         vctrs_0.6.5             ROCR_1.0-12            
# [10] spatstat.explore_3.3-4  Rsamtools_2.22.0        RcppRoll_0.3.1         
# [13] htmltools_0.5.8.1       sctransform_0.4.1       parallelly_1.45.1      
# [16] KernSmooth_2.23-24      htmlwidgets_1.6.4       ica_1.0-3              
# [19] plyr_1.8.9              plotly_4.12.0           zoo_1.8-13             
# [22] igraph_2.1.4            mime_0.12               lifecycle_1.0.5        
# [25] pkgconfig_2.0.3         Matrix_1.7-3            R6_2.6.1               
# [28] fastmap_1.2.0           GenomeInfoDbData_1.2.13 fitdistrplus_1.2-6     
# [31] shiny_1.13.0            digest_0.6.37           tensor_1.5.1           
# [34] RSpectra_0.16-2         irlba_2.3.5.1           labeling_0.4.3         
# [37] progressr_0.18.0        spatstat.sparse_3.1-0   httr_1.4.8             
# [40] polyclip_1.10-7         abind_1.4-8             compiler_4.4.1         
# [43] withr_3.0.2             S7_0.2.0                BiocParallel_1.40.2    
# [46] fastDummies_1.7.5       MASS_7.3-60.2           tools_4.4.1            
# [49] lmtest_0.9-40           otel_0.2.0              httpuv_1.6.15          
# [52] future.apply_1.20.2     goftest_1.2-3           glue_1.8.0             
# [55] nlme_3.1-164            promises_1.5.0          grid_4.4.1             
# [58] Rtsne_0.17              cluster_2.1.6           reshape2_1.4.4         
# [61] generics_0.1.4          gtable_0.3.6            spatstat.data_3.1-9    
# [64] tidyr_1.3.1             data.table_1.17.0       XVector_0.46.0         
# [67] spatstat.geom_3.3-5     RcppAnnoy_0.0.22        ggrepel_0.9.6          
# [70] RANN_2.6.2              pillar_1.11.1           stringr_1.6.0          
# [73] spam_2.11-1             RcppHNSW_0.6.0          later_1.4.1            
# [76] splines_4.4.1           lattice_0.22-6          survival_3.6-4         
# [79] deldir_2.0-4            tidyselect_1.2.1        Biostrings_2.74.1      
# [82] miniUI_0.1.2            pbapply_1.7-4           gridExtra_2.3          
# [85] scattermore_1.2         matrixStats_1.5.0       stringi_1.8.4          
# [88] UCSC.utils_1.2.0        lazyeval_0.2.2          codetools_0.2-20       
# [91] tibble_3.2.1            cli_3.6.4               uwot_0.2.3             
# [94] xtable_1.8-8            reticulate_1.41.0       dichromat_2.0-0.1      
# [97] Rcpp_1.0.14             globals_0.19.0          spatstat.random_3.3-2  
# [100] png_0.1-8               spatstat.univar_3.1-1   parallel_4.4.1         
# [103] dotCall64_1.2           bitops_1.0-9            listenv_0.10.0         
# [106] viridisLite_0.4.3       scales_1.4.0            ggridges_0.5.7         
# [109] purrr_1.0.4             crayon_1.5.3            rlang_1.1.5            
# [112] cowplot_1.2.0           fastmatch_1.1-8