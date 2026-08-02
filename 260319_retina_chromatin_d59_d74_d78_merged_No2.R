# GSE183684

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

ATAC_merged <- readRDS("260228_retina_chromatin_ATAC_merged.rds")

str(ATAC_merged)
head(ATAC_merged)

ATAC_merged[['peaks']]
# ChromatinAssay data with 210972 features for 27440 cells
# Variable features: 0 
# Genome: 
#   Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 3 

grange <- granges(ATAC_merged)

peaks.keep <- seqnames(granges(ATAC_merged)) %in% standardChromosomes(granges(ATAC_merged))
ATAC_merged <- ATAC_merged[as.vector(peaks.keep), ]

# compute nucleosome signal score per cell
ATAC_merged <- NucleosomeSignal(object = ATAC_merged)

# compute TSS enrichment score per cell
ATAC_merged <- TSSEnrichment(object = ATAC_merged)

# not run
# peak_ranges should be a set of genomic ranges spanning the set of peaks to be quantified per cell
peak_matrix <- FeatureMatrix(
  fragments = Fragments(ATAC_merged),
  features = grange
)

# saveRDS(ATAC_merged, "260301_retina_chromatin_merged_No2_sub1.rds")
# saveRDS(peak_matrix, "260301_retina_chromatin_merged_No2_sub2.rds")
ATAC_merged <- readRDS("260301_retina_chromatin_merged_No2_sub1.rds")

# # not run
# total_fragments <- CountFragments('GSM5567518_d59_fragments.tsv.gz')
# rownames(total_fragments) <- total_fragments$CB
# ATAC_merged$fragments <- total_fragments[colnames(ATAC_merged), "frequency_count"]

Fragments(ATAC_merged)
# [[1]]
# A Fragment object for 8071 cells
# 
# [[2]]
# A Fragment object for 7095 cells
# 
# [[3]]
# A Fragment object for 12274 cells

frag_list <- Fragments(ATAC_merged)

sample_ids <- c("d59", "d74", "d78")

frag_counts_list <- lapply(seq_along(frag_list), function(i) {
  
  counts <- CountFragments(frag_list[[i]]@path)
  
  # Add correct biological prefix
  counts$CB <- paste0(sample_ids[i], "_", counts$CB)
  
  return(counts)
})

total_fragments <- do.call(rbind, frag_counts_list)

rownames(total_fragments) <- total_fragments$CB

ATAC_merged$fragments <- 
  total_fragments[colnames(ATAC_merged), "frequency_count"]

sum(is.na(ATAC_merged$fragments))
# [1] 0
summary(ATAC_merged$fragments)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1001    6335   12298   15726   21206  311659 

ATAC_merged <- FRiP(
  object = ATAC_merged,
  assay = "peaks",
  total.fragments = "fragments"
)

# saveRDS(ATAC_merged, "260301_retina_chromatin_merged_No2_sub3.rds")

# add blacklist ratio
ATAC_merged$blacklist_ratio <- FractionCountsInRegion(
  object = ATAC_merged, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

DensityScatter(ATAC_merged, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

ATAC_merged$nucleosome_group <- ifelse(ATAC_merged$nucleosome_signal > 1.5, 'NS > 1.5', 'NS < 1.5')
FragmentHistogram(object = ATAC_merged, group.by = 'nucleosome_group')
# Warning messages:
#   1: Removed 575 rows containing non-finite outside the scale range
# (`stat_bin()`). 
# 2: Removed 4 rows containing missing values or values outside the scale range
# (`geom_bar()`). 

VlnPlot(
  object = ATAC_merged,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
  pt.size = 0.1,
  ncol = 5
)

ATAC_merged_QC <- subset(
  x = ATAC_merged,
  subset = nCount_peaks > 2000 &
    nCount_peaks < 30000 &
    FRiP > 0.2 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 1.5 &
    TSS.enrichment > 2.5
)
ATAC_merged_QC
# An object of class Seurat 
# 210972 features across 22231 samples within 1 assay 
# Active assay: peaks (210972 features, 0 variable features)
# 2 layers present: counts, data

# saveRDS(object = ATAC_merged_QC, file = "260301_retina_chromatin_merged_No2.rds")

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
#   [1] patchwork_1.3.2      ggplot2_4.0.0        GenomicRanges_1.58.0
# [4] GenomeInfoDb_1.42.3  IRanges_2.40.1       S4Vectors_0.44.0    
# [7] BiocGenerics_0.52.0  Seurat_5.3.0         SeuratObject_5.0.2  
# [10] sp_2.2-0             Signac_1.16.0       
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.18.0       jsonlite_1.9.0         
# [4] magrittr_2.0.3          ggbeeswarm_0.7.3        spatstat.utils_3.1-2   
# [7] farver_2.1.2            zlibbioc_1.52.0         vctrs_0.6.5            
# [10] ROCR_1.0-12             spatstat.explore_3.3-4  Rsamtools_2.22.0       
# [13] RcppRoll_0.3.1          htmltools_0.5.8.1       sctransform_0.4.1      
# [16] parallelly_1.45.1       KernSmooth_2.23-24      htmlwidgets_1.6.4      
# [19] ica_1.0-3               plyr_1.8.9              plotly_4.12.0          
# [22] zoo_1.8-13              igraph_2.1.4            mime_0.12              
# [25] lifecycle_1.0.5         pkgconfig_2.0.3         Matrix_1.7-3           
# [28] R6_2.6.1                fastmap_1.2.0           GenomeInfoDbData_1.2.13
# [31] fitdistrplus_1.2-6      future_1.69.0           shiny_1.13.0           
# [34] digest_0.6.37           tensor_1.5.1            RSpectra_0.16-2        
# [37] irlba_2.3.5.1           labeling_0.4.3          progressr_0.18.0       
# [40] spatstat.sparse_3.1-0   httr_1.4.8              polyclip_1.10-7        
# [43] abind_1.4-8             compiler_4.4.1          withr_3.0.2            
# [46] S7_0.2.0                BiocParallel_1.40.2     fastDummies_1.7.5      
# [49] MASS_7.3-60.2           tools_4.4.1             vipor_0.4.7            
# [52] lmtest_0.9-40           otel_0.2.0              beeswarm_0.4.0         
# [55] httpuv_1.6.15           future.apply_1.20.2     goftest_1.2-3          
# [58] glue_1.8.0              nlme_3.1-164            promises_1.5.0         
# [61] grid_4.4.1              Rtsne_0.17              cluster_2.1.6          
# [64] reshape2_1.4.4          generics_0.1.4          gtable_0.3.6           
# [67] spatstat.data_3.1-9     tidyr_1.3.1             data.table_1.17.0      
# [70] XVector_0.46.0          spatstat.geom_3.3-5     RcppAnnoy_0.0.22       
# [73] ggrepel_0.9.6           RANN_2.6.2              pillar_1.11.1          
# [76] stringr_1.6.0           spam_2.11-1             RcppHNSW_0.6.0         
# [79] later_1.4.1             splines_4.4.1           dplyr_1.1.4            
# [82] lattice_0.22-6          survival_3.6-4          deldir_2.0-4           
# [85] tidyselect_1.2.1        Biostrings_2.74.1       miniUI_0.1.2           
# [88] pbapply_1.7-4           gridExtra_2.3           scattermore_1.2        
# [91] matrixStats_1.5.0       stringi_1.8.4           UCSC.utils_1.2.0       
# [94] lazyeval_0.2.2          codetools_0.2-20        tibble_3.2.1           
# [97] cli_3.6.4               uwot_0.2.3              xtable_1.8-8           
# [100] reticulate_1.41.0       dichromat_2.0-0.1       Rcpp_1.0.14            
# [103] globals_0.19.0          spatstat.random_3.3-2   png_0.1-8              
# [106] ggrastr_1.0.2           spatstat.univar_3.1-1   parallel_4.4.1         
# [109] dotCall64_1.2           bitops_1.0-9            listenv_0.10.0         
# [112] viridisLite_0.4.3       scales_1.4.0            ggridges_0.5.7         
# [115] purrr_1.0.4             crayon_1.5.3            rlang_1.1.5            
# [118] cowplot_1.2.0           fastmatch_1.1-8