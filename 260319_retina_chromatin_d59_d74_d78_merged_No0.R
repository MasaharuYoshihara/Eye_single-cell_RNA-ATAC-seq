# GSE183684

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Signac)
library(EnsDb.Hsapiens.v86)

# RNA

D59 <- ReadMtx(
  mtx = "GSM5567526_d59_matrix.mtx.gz",
  features = "GSM5567526_d59_features.tsv.gz",
  cells = "GSM5567526_d59_barcodes.tsv.gz",
  feature.column = 2
)
RNA_d59 <- CreateSeuratObject(counts = D59, min.cells = 3, min.features = 200)

D74 <- ReadMtx(
  mtx = "GSM5567527_d74_matrix.mtx.gz",
  features = "GSM5567527_d74_features.tsv.gz",
  cells = "GSM5567527_d74_barcodes.tsv.gz",
  feature.column = 2
)
RNA_d74 <- CreateSeuratObject(counts = D74, min.cells = 3, min.features = 200)

D78 <- ReadMtx(
  mtx = "GSM5567528_d78_matrix.mtx.gz",
  features = "GSM5567528_d78_features.tsv.gz",
  cells = "GSM5567528_d78_barcodes.tsv.gz",
  feature.column = 2
)
RNA_d78 <- CreateSeuratObject(counts = D78, min.cells = 3, min.features = 200)

RNA_d59$dataset <- 'day 59'
RNA_d74$dataset <- 'day 74'
RNA_d78$dataset <- 'day 78'

RNA_merged <- merge(
  x = RNA_d59,
  y = list(RNA_d74, RNA_d78),
  add.cell.ids = c("d59", "d74", "d78")
)

# ATAC
ATAC_d59 <- readRDS("241108_retina_chromatin_d59_No1.rds")
ATAC_d74 <- readRDS("241108_retina_chromatin_d74_No1.rds")
ATAC_d78 <- readRDS("241108_retina_chromatin_d78_No1.rds")

ATAC_d59@assays$peaks@fragments <- list()
ATAC_d74@assays$peaks@fragments <- list()
ATAC_d78@assays$peaks@fragments <- list()

frag_d59 <- CreateFragmentObject(
  path = "C:/Users/ana_emb/Documents/Watabe/GSM5567518_d59_fragments.tsv.gz",
  cells = colnames(ATAC_d59)
)
frag_d74 <- CreateFragmentObject(
  path = "C:/Users/ana_emb/Documents/Watabe/GSM5567519_d74_fragments.tsv.gz",
  cells = colnames(ATAC_d74)
)
frag_d78 <- CreateFragmentObject(
  path = "C:/Users/ana_emb/Documents/Watabe/GSM5567520_d78_fragments.tsv.gz",
  cells = colnames(ATAC_d78)
)

ATAC_d59@assays$peaks@fragments <- list(frag_d59)
ATAC_d74@assays$peaks@fragments <- list(frag_d74)
ATAC_d78@assays$peaks@fragments <- list(frag_d78)

ATAC_merged <- merge(
  x = ATAC_d59,
  y = list(ATAC_d74, ATAC_d78),
  add.cell.ids = c("d59", "d74", "d78")
)

## saveRDS(ATAC_merged, "260228_retina_chromatin_ATAC_merged.rds")
## saveRDS(RNA_merged, "260228_retina_chromatin_RNA_merged.rds")
# RNA_merged <- readRDS("260228_retina_chromatin_RNA_merged.rds")
# ATAC_merged <- readRDS("260228_retina_chromatin_ATAC_merged.rds")

sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Japanese_Japan.utf8  LC_CTYPE=Japanese_Japan.utf8    LC_MONETARY=Japanese_Japan.utf8
# [4] LC_NUMERIC=C                    LC_TIME=Japanese_Japan.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] EnsDb.Hsapiens.v86_2.99.0 ensembldb_2.30.0          AnnotationFilter_1.30.0  
# [4] GenomicFeatures_1.58.0    AnnotationDbi_1.68.0      Biobase_2.66.0           
# [7] GenomicRanges_1.58.0      GenomeInfoDb_1.42.3       IRanges_2.40.1           
# [10] S4Vectors_0.44.0          BiocGenerics_0.52.0       Signac_1.16.0            
# [13] ggplot2_4.0.0             patchwork_1.3.2           Seurat_5.3.0             
# [16] SeuratObject_5.0.2        sp_2.2-0                  dplyr_1.1.4              
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.4.1               later_1.4.1                
# [4] BiocIO_1.16.0               bitops_1.0-9                tibble_3.2.1               
# [7] polyclip_1.10-7             XML_3.99-0.22               fastDummies_1.7.5          
# [10] lifecycle_1.0.5             globals_0.19.0              lattice_0.22-6             
# [13] MASS_7.3-60.2               magrittr_2.0.3              plotly_4.12.0              
# [16] yaml_2.3.12                 httpuv_1.6.15               otel_0.2.0                 
# [19] sctransform_0.4.1           spam_2.11-1                 spatstat.sparse_3.1-0      
# [22] reticulate_1.41.0           cowplot_1.2.0               pbapply_1.7-4              
# [25] DBI_1.3.0                   RColorBrewer_1.1-3          abind_1.4-8                
# [28] zlibbioc_1.52.0             Rtsne_0.17                  purrr_1.0.4                
# [31] RCurl_1.98-1.17             GenomeInfoDbData_1.2.13     ggrepel_0.9.6              
# [34] irlba_2.3.5.1               listenv_0.10.0              spatstat.utils_3.1-2       
# [37] goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.3-2      
# [40] fitdistrplus_1.2-6          parallelly_1.45.1           codetools_0.2-20           
# [43] DelayedArray_0.32.0         RcppRoll_0.3.1              tidyselect_1.2.1           
# [46] UCSC.utils_1.2.0            farver_2.1.2                matrixStats_1.5.0          
# [49] spatstat.explore_3.3-4      GenomicAlignments_1.42.0    jsonlite_1.9.0             
# [52] progressr_0.18.0            ggridges_0.5.7              survival_3.6-4             
# [55] tools_4.4.1                 ica_1.0-3                   Rcpp_1.0.14                
# [58] glue_1.8.0                  gridExtra_2.3               SparseArray_1.6.2          
# [61] MatrixGenerics_1.18.1       withr_3.0.2                 fastmap_1.2.0              
# [64] digest_0.6.37               R6_2.6.1                    mime_0.12                  
# [67] scattermore_1.2             tensor_1.5.1                dichromat_2.0-0.1          
# [70] spatstat.data_3.1-9         RSQLite_2.4.6               tidyr_1.3.1                
# [73] generics_0.1.4              data.table_1.17.0           rtracklayer_1.66.0         
# [76] httr_1.4.8                  htmlwidgets_1.6.4           S4Arrays_1.6.0             
# [79] uwot_0.2.3                  pkgconfig_2.0.3             gtable_0.3.6               
# [82] blob_1.3.0                  lmtest_0.9-40               S7_0.2.0                   
# [85] XVector_0.46.0              htmltools_0.5.8.1           dotCall64_1.2              
# [88] ProtGenerics_1.38.0         scales_1.4.0                png_0.1-8                  
# [91] spatstat.univar_3.1-1       rstudioapi_0.18.0           reshape2_1.4.4             
# [94] rjson_0.2.23                nlme_3.1-164                curl_7.0.0                 
# [97] zoo_1.8-13                  cachem_1.1.0                stringr_1.6.0              
# [100] KernSmooth_2.23-24          parallel_4.4.1              miniUI_0.1.2               
# [103] restfulr_0.0.16             pillar_1.11.1               grid_4.4.1                 
# [106] vctrs_0.6.5                 RANN_2.6.2                  promises_1.5.0             
# [109] xtable_1.8-8                cluster_2.1.6               cli_3.6.4                  
# [112] compiler_4.4.1              Rsamtools_2.22.0            rlang_1.1.5                
# [115] crayon_1.5.3                future.apply_1.20.2         plyr_1.8.9                 
# [118] stringi_1.8.4               viridisLite_0.4.3           deldir_2.0-4               
# [121] BiocParallel_1.40.2         Biostrings_2.74.1           lazyeval_0.2.2             
# [124] spatstat.geom_3.3-5         Matrix_1.7-3                RcppHNSW_0.6.0             
# [127] bit64_4.6.0-1               future_1.69.0               KEGGREST_1.46.0            
# [130] shiny_1.13.0                SummarizedExperiment_1.36.0 ROCR_1.0-12                
# [133] igraph_2.1.4                memoise_2.0.1               fastmatch_1.1-8            
# [136] bit_4.6.0