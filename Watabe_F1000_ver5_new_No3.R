# Watabe et al. F1000Research ver5
# This file is Watabe_F1000_ver5_new_No3.R for scATACseq analysis
# Use Watabe_F1000_ver5_new_No1.R for scRNAseq processing
# Use Watabe_F1000_ver5_new_No2.R for scATACseq processing

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Signac)
library(GenomicRanges)
library(scales)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(RSQLite)
library(ggseqlogo)
set.seed(123)

ATAC_merged <- readRDS("260715_ATAC_annotated.rds") # From Watabe_F1000_ver5_new_No2.R
# str(ATAC_merged)
# DimPlot(ATAC_merged,group.by = 'predicted.id')

RNA_merged <- readRDS("Thomas_et_al_scRNAseq_query.rds")

celltypes <- levels(RNA_merged)
color_vector <- hue_pal()(length(celltypes))
rna_colors <- setNames(color_vector, celltypes)
show_col(rna_colors)
ATAC_merged$celltype <- factor(ATAC_merged$predicted.id, levels = names(rna_colors))

pfm_db <- JASPAR2024()
conn <- dbConnect(RSQLite::SQLite(), db(pfm_db))
tfbs <- TFBSTools::getMatrixSet(conn, list(species = "Homo sapiens", collection = "CORE"))
ATAC_merged <- AddMotifs(
  object = ATAC_merged,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = tfbs,
  assay = "peaks"
)
DefaultAssay(ATAC_merged) <- 'peaks'
Idents(object = ATAC_merged) <- "predicted.id"
ATAC_merged <- SortIdents(ATAC_merged)
saveRDS(ATAC_merged, "260716_ATAC_addmotif.rds")

cell_types <- c("Amacrine cells", "Astrocytes",
                "Bipolar cells", "Cone precursors",
                "Cones", "Horizontal cells",
                "Microglia", "Muller glia",
                "Retinal ganglion cells", "Retinal progenitors",
                "Rods", "Transitory")

da_peaks_RPC <- FindMarkers(
      object = ATAC_merged,
      ident.1 = "Retinal progenitors",
      assay = "peaks",
      only.pos = TRUE,
      test.use = 'LR',
      min.pct = 0.05,
      latent.vars = 'nCount_peaks'
)
saveRDS(da_peaks_RPC, "260716_da_peaks_RPC.rds")

ATAC_merged <-readRDS("260716_ATAC_addmotif.rds")

DimPlot(ATAC_merged, group.by = 'predicted.id')
# save image as DimPlot_ATAC_annotated_260716.tiff W800 H600

da_peaks_RPC_Transitory <- FindMarkers(
  object = ATAC_merged,
  ident.1 = c("Retinal progenitors", "Transitory"),
  assay = "peaks",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
saveRDS(da_peaks_RPC_Transitory, "260716_da_peaks_RPC_Transitory.rds")

# sessionInfo()
# R version 4.5.2 (2025-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# LAPACK version 3.12.1
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.utf8 
# [2] LC_CTYPE=English_United Kingdom.utf8   
# [3] LC_MONETARY=English_United Kingdom.utf8
# [4] LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils    
# [6] datasets  methods   base     
# 
# other attached packages:
#   [1] ggseqlogo_0.2.2                  
# [2] RSQLite_3.53.2                   
# [3] BSgenome.Hsapiens.UCSC.hg38_1.4.5
# [4] BSgenome_1.78.0                  
# [5] rtracklayer_1.70.1               
# [6] BiocIO_1.20.0                    
# [7] Biostrings_2.78.0                
# [8] XVector_0.50.0                   
# [9] GenomeInfoDb_1.46.2              
# [10] TFBSTools_1.48.0                 
# [11] JASPAR2024_0.99.7                
# [12] BiocFileCache_3.0.0              
# [13] dbplyr_2.6.0                     
# [14] scales_1.4.0                     
# [15] GenomicRanges_1.62.1             
# [16] Seqinfo_1.0.0                    
# [17] IRanges_2.44.0                   
# [18] S4Vectors_0.48.1                 
# [19] BiocGenerics_0.56.0              
# [20] generics_0.1.4                   
# [21] Signac_1.17.1                    
# [22] ggplot2_4.0.2                    
# [23] patchwork_1.3.2                  
# [24] Seurat_5.4.0                     
# [25] SeuratObject_5.3.0               
# [26] sp_2.2-1                         
# [27] dplyr_1.2.0                      
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.23            splines_4.5.2              
# [3] later_1.4.8                 filelock_1.0.3             
# [5] bitops_1.0-9                tibble_3.3.1               
# [7] polyclip_1.10-7             XML_3.99-0.23              
# [9] DirichletMultinomial_1.52.0 fastDummies_1.7.6          
# [11] httr2_1.2.2                 lifecycle_1.0.5            
# [13] pwalign_1.6.0               globals_0.19.1             
# [15] lattice_0.22-9              MASS_7.3-65                
# [17] magrittr_2.0.4              plotly_4.12.0              
# [19] yaml_2.3.12                 httpuv_1.6.16              
# [21] otel_0.2.0                  sctransform_0.4.3          
# [23] spam_2.11-3                 spatstat.sparse_3.1-0      
# [25] reticulate_1.45.0           cowplot_1.2.0              
# [27] pbapply_1.7-4               DBI_1.3.0                  
# [29] RColorBrewer_1.1-3          abind_1.4-8                
# [31] Rtsne_0.17                  purrr_1.2.1                
# [33] RCurl_1.98-1.19             rappdirs_0.3.4             
# [35] ggrepel_0.9.7               irlba_2.3.7                
# [37] listenv_0.10.1              spatstat.utils_3.2-1       
# [39] seqLogo_1.76.0              goftest_1.2-3              
# [41] RSpectra_0.16-2             spatstat.random_3.4-4      
# [43] fitdistrplus_1.2-6          parallelly_1.46.1          
# [45] codetools_0.2-20            DelayedArray_0.36.1        
# [47] RcppRoll_0.3.2              tidyselect_1.2.1           
# [49] UCSC.utils_1.6.1            farver_2.1.2               
# [51] matrixStats_1.5.0           spatstat.explore_3.7-0     
# [53] GenomicAlignments_1.46.0    jsonlite_2.0.0             
# [55] progressr_0.19.0            ggridges_0.5.7             
# [57] survival_3.8-6              tools_4.5.2                
# [59] TFMPvalue_1.0.0             ica_1.0-3                  
# [61] Rcpp_1.1.1                  glue_1.8.0                 
# [63] gridExtra_2.3               SparseArray_1.10.10        
# [65] MatrixGenerics_1.22.0       withr_3.0.3                
# [67] fastmap_1.2.0               caTools_1.18.3             
# [69] digest_0.6.39               R6_2.6.1                   
# [71] mime_0.13                   scattermore_1.2            
# [73] gtools_3.9.5                tensor_1.5.1               
# [75] dichromat_2.0-0.1           spatstat.data_3.1-9        
# [77] cigarillo_1.0.0             tidyr_1.3.2                
# [79] data.table_1.18.2.1         httr_1.4.8                 
# [81] htmlwidgets_1.6.4           S4Arrays_1.10.1            
# [83] uwot_0.2.4                  pkgconfig_2.0.3            
# [85] gtable_0.3.6                blob_1.3.0                 
# [87] lmtest_0.9-40               S7_0.2.1                   
# [89] htmltools_0.5.9             dotCall64_1.2              
# [91] Biobase_2.70.0              png_0.1-8                  
# [93] spatstat.univar_3.1-6       rstudioapi_0.19.0          
# [95] reshape2_1.4.5              rjson_0.2.23               
# [97] nlme_3.1-168                curl_7.1.0                 
# [99] cachem_1.1.0                zoo_1.8-15                 
# [101] stringr_1.6.0               KernSmooth_2.23-26         
# [103] parallel_4.5.2              miniUI_0.1.2               
# [105] restfulr_0.0.17             pillar_1.11.1              
# [107] grid_4.5.2                  vctrs_0.7.1                
# [109] RANN_2.6.2                  promises_1.5.0             
# [111] xtable_1.8-8                cluster_2.1.8.2            
# [113] cli_3.6.6                   compiler_4.5.2             
# [115] Rsamtools_2.26.0            rlang_1.2.0                
# [117] crayon_1.5.3                future.apply_1.20.2        
# [119] labeling_0.4.3              plyr_1.8.9                 
# [121] stringi_1.8.7               viridisLite_0.4.3          
# [123] deldir_2.0-4                BiocParallel_1.44.0        
# [125] lazyeval_0.2.2              spatstat.geom_3.7-0        
# [127] Matrix_1.7-4                RcppHNSW_0.6.0             
# [129] sparseMatrixStats_1.22.0    bit64_4.8.2                
# [131] future_1.70.0               shiny_1.13.0               
# [133] SummarizedExperiment_1.40.0 ROCR_1.0-12                
# [135] igraph_2.2.2                memoise_2.0.1              
# [137] fastmatch_1.1-8             bit_4.6.0

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Signac)
library(GenomicRanges)
library(scales)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(RSQLite)
library(ggseqlogo)
set.seed(123)

ATAC_merged <-readRDS("260716_ATAC_addmotif.rds")
da_peaks_RPC <- readRDS("260716_da_peaks_RPC.rds")
top.da.peak_RPC <- rownames(da_peaks_RPC[da_peaks_RPC$p_val < 0.005 & da_peaks_RPC$pct.1 > 0.2, ])
enriched.motifs <- FindMotifs(
  object = ATAC_merged,
  features = top.da.peak_RPC,
  assay = "peaks"
)
head(enriched.motifs)
# motif observed background percent.observed
# MA2328.1 MA2328.1     1662      10069         59.29361
# MA1959.2 MA1959.2     1626       9880         58.00928
# MA1513.2 MA1513.2     1778      11058         63.43204
# MA0746.3 MA0746.3     1726      11234         61.57688
# MA0742.2 MA0742.2     1977      13270         70.53157
# MA0516.3 MA0516.3     2035      13910         72.60078
# percent.background fold.enrichment pvalue motif.name
# MA2328.1            25.1725        2.355492      0      ZBED4
# MA1959.2            24.7000        2.348554      0       KLF7
# MA1513.2            27.6450        2.294521      0      KLF15
# MA0746.3            28.0850        2.192518      0        SP3
# MA0742.2            33.1750        2.126046      0      KLF12
# MA0516.3            34.7750        2.087729      0        SP2
# p.adjust
# MA2328.1        0
# MA1959.2        0
# MA1513.2        0
# MA0746.3        0
# MA0742.2        0
# MA0516.3        0
motifs_of_interest <- c("MA0700.3", "MA0069.1", "MA0718.2", "MA0726.2")
# MA0700.3 = LHX2 (former MA0700.2)
# MA0069.1 = PAX6
# MA0718.2 = RAX (former A0718.1)
# MA0726.2 = VSX2 (former MA0726.1)
# ID changed from JASPAR2020
enrichment_moi <- enriched.motifs[motifs_of_interest, , drop = FALSE]
# motif observed background percent.observed
# MA0700.3 MA0700.3      975       5887        34.784160
# MA0069.1 MA0069.1      200       1535         7.135212
# MA0718.2 MA0718.2      975       5887        34.784160
# MA0726.2 MA0726.2      352       1771        12.557974
# percent.background fold.enrichment        pvalue
# MA0700.3            14.7175        2.363456 9.787746e-169
# MA0069.1             3.8375        1.859339  9.437479e-18
# MA0718.2            14.7175        2.363456 9.787746e-169
# MA0726.2             4.4275        2.836358  2.026034e-74
# motif.name      p.adjust
# MA0700.3       LHX2 9.272601e-168
# MA0069.1       PAX6  1.003690e-17
# MA0718.2        RAX 9.272601e-168
# MA0726.2       VSX2  5.834979e-74
write.csv(enrichment_moi, "260719_enrichment_moi_RPC.csv")
DefaultAssay(ATAC_merged) <- "peaks"
MotifPlot(
  object = ATAC_merged,
  motifs = head(rownames(enriched.motifs))
)
MotifPlot(
  object = ATAC_merged,
  motifs = c("LHX2", "PAX6", "RAX", "VSX2")
)
# save image as RPC_enriched_motif_260719.tiff W800 H600

ATAC_merged <-readRDS("260716_ATAC_addmotif.rds")
da_peaks_RPC_Transitory <- readRDS("260716_da_peaks_RPC_Transitory.rds")
top.da.peak_RPC_Transitory <- rownames(da_peaks_RPC_Transitory[da_peaks_RPC_Transitory$p_val < 0.005 & da_peaks_RPC_Transitory$pct.1 > 0.2, ])
enriched.motifs_RPC_Transitory <- FindMotifs(
  object = ATAC_merged,
  features = top.da.peak_RPC_Transitory,
  assay = "peaks"
)
head(enriched.motifs_RPC_Transitory)
# motif observed background percent.observed
# MA0471.3 MA0471.3     1466       6432         42.89058
# MA0732.2 MA0732.2     1509       6851         44.14862
# MA1713.2 MA1713.2     1741       7957         50.93622
# MA1959.2 MA1959.2     2166      10266         63.37039
# MA2328.1 MA2328.1     2242      10656         65.59391
# MA1516.2 MA1516.2     1720       8212         50.32183
# percent.background fold.enrichment pvalue motif.name
# MA0471.3            16.0800        2.667325      0       E2F6
# MA0732.2            17.1275        2.577646      0       EGR3
# MA1713.2            19.8925        2.560574      0     ZNF610
# MA1959.2            25.6650        2.469137      0       KLF7
# MA2328.1            26.6400        2.462234      0      ZBED4
# MA1516.2            20.5300        2.451136      0       KLF3
# p.adjust
# MA0471.3        0
# MA0732.2        0
# MA1713.2        0
# MA1959.2        0
# MA2328.1        0
# MA1516.2        0
write.csv(enriched.motifs_RPC_Transitory, "enriched_motifs_RPC_Transitory_270725.csv")
motifs_of_interest <- c("MA0700.3", "MA0069.1", "MA0718.2", "MA0726.2")
# MA0700.3 = LHX2 (former MA0700.2)
# MA0069.1 = PAX6
# MA0718.2 = RAX (former A0718.1)
# MA0726.2 = VSX2 (former MA0726.1)
# ID changed from JASPAR2020
enrichment_moi_RPC_Transitory <- enriched.motifs_RPC_Transitory[motifs_of_interest, , drop = FALSE]
# motif observed background percent.observed
# MA0700.3 MA0700.3     1124       5766        32.884728
# MA0069.1 MA0069.1      249       1546         7.284962
# MA0718.2 MA0718.2     1124       5766        32.884728
# MA0726.2 MA0726.2      401       1752        11.732007
# percent.background fold.enrichment        pvalue
# MA0700.3             14.415        2.281285 1.263442e-182
# MA0069.1              3.865        1.884854  6.202046e-23
# MA0718.2             14.415        2.281285 1.263442e-182
# MA0726.2              4.380        2.678540  2.153548e-78
# motif.name      p.adjust
# MA0700.3       LHX2 1.022110e-181
# MA0069.1       PAX6  6.755632e-23
# MA0718.2        RAX 1.022110e-181
# MA0726.2       VSX2  5.067172e-78
write.csv(enrichment_moi_RPC_Transitory, "260725_enrichment_moi_RPC_Transitory.csv")
DefaultAssay(ATAC_merged) <- "peaks"
MotifPlot(
  object = ATAC_merged,
  motifs = head(rownames(enriched.motifs_RPC_Transitory))
)
MotifPlot(
  object = ATAC_merged,
  motifs = c("LHX2", "PAX6", "RAX", "VSX2")
)
# save image as RPC_Transitory_enriched_motif_260725.tiff W800 H600

# footprinting
library(motifmatchr)
library(EnsDb.Hsapiens.v86)

library(Rsamtools)
file_names <- c("GSM5567518_d59_fragments.tsv.gz",
                "GSM5567519_d74_fragments.tsv.gz",
                "GSM5567520_d78_fragments.tsv.gz")
frag_list <- Fragments(ATAC_merged)
updated_frags <- lapply(seq_along(frag_list), function(i) {
  orig_frag <- frag_list[[i]]
  correct_path <- file_names[i] 
  if (!file.exists(paste0(correct_path, ".tbi"))) {
    indexTabix(file = correct_path, format = "bed")
  }
  CreateFragmentObject(
    path = correct_path,
    cells = GetFragmentData(orig_frag, slot = "cells")
  )
})
Fragments(ATAC_merged) <- NULL
Fragments(ATAC_merged) <- updated_frags

ATAC_merged_fp <- Footprint(
  object = ATAC_merged,
  motif.name = c("MA0700.3", "MA0069.1", "MA0718.2", "MA0726.2"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(ATAC_merged_fp, "260719_ATAC_RPC_footprint.rds")

sessionInfo()
# R version 4.5.2 (2025-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# LAPACK version 3.12.1
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.utf8 
# [2] LC_CTYPE=English_United Kingdom.utf8   
# [3] LC_MONETARY=English_United Kingdom.utf8
# [4] LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.utf8    
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils    
# [6] datasets  methods   base     
# 
# other attached packages:
#   [1] Rsamtools_2.26.0                 
# [2] EnsDb.Hsapiens.v86_2.99.0        
# [3] ensembldb_2.34.0                 
# [4] AnnotationFilter_1.34.0          
# [5] GenomicFeatures_1.62.0           
# [6] AnnotationDbi_1.72.0             
# [7] Biobase_2.70.0                   
# [8] motifmatchr_1.32.0               
# [9] ggseqlogo_0.2.2                  
# [10] RSQLite_3.53.2                   
# [11] BSgenome.Hsapiens.UCSC.hg38_1.4.5
# [12] BSgenome_1.78.0                  
# [13] rtracklayer_1.70.1               
# [14] BiocIO_1.20.0                    
# [15] Biostrings_2.78.0                
# [16] XVector_0.50.0                   
# [17] GenomeInfoDb_1.46.2              
# [18] TFBSTools_1.48.0                 
# [19] JASPAR2024_0.99.7                
# [20] BiocFileCache_3.0.0              
# [21] dbplyr_2.6.0                     
# [22] scales_1.4.0                     
# [23] GenomicRanges_1.62.1             
# [24] Seqinfo_1.0.0                    
# [25] IRanges_2.44.0                   
# [26] S4Vectors_0.48.1                 
# [27] BiocGenerics_0.56.0              
# [28] generics_0.1.4                   
# [29] Signac_1.17.1                    
# [30] ggplot2_4.0.2                    
# [31] patchwork_1.3.2                  
# [32] Seurat_5.4.0                     
# [33] SeuratObject_5.3.0               
# [34] sp_2.2-1                         
# [35] dplyr_1.2.0                      
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.23            splines_4.5.2              
# [3] later_1.4.8                 filelock_1.0.3             
# [5] bitops_1.0-9                tibble_3.3.1               
# [7] polyclip_1.10-7             XML_3.99-0.23              
# [9] DirichletMultinomial_1.52.0 fastDummies_1.7.6          
# [11] httr2_1.2.2                 lifecycle_1.0.5            
# [13] pwalign_1.6.0               globals_0.19.1             
# [15] lattice_0.22-9              MASS_7.3-65                
# [17] magrittr_2.0.4              plotly_4.12.0              
# [19] yaml_2.3.12                 httpuv_1.6.16              
# [21] otel_0.2.0                  sctransform_0.4.3          
# [23] spam_2.11-3                 spatstat.sparse_3.1-0      
# [25] reticulate_1.45.0           cowplot_1.2.0              
# [27] pbapply_1.7-4               DBI_1.3.0                  
# [29] RColorBrewer_1.1-3          abind_1.4-8                
# [31] Rtsne_0.17                  purrr_1.2.1                
# [33] RCurl_1.98-1.19             rappdirs_0.3.4             
# [35] ggrepel_0.9.7               irlba_2.3.7                
# [37] listenv_0.10.1              spatstat.utils_3.2-1       
# [39] seqLogo_1.76.0              goftest_1.2-3              
# [41] RSpectra_0.16-2             spatstat.random_3.4-4      
# [43] fitdistrplus_1.2-6          parallelly_1.46.1          
# [45] codetools_0.2-20            DelayedArray_0.36.1        
# [47] RcppRoll_0.3.2              tidyselect_1.2.1           
# [49] UCSC.utils_1.6.1            farver_2.1.2               
# [51] matrixStats_1.5.0           spatstat.explore_3.7-0     
# [53] GenomicAlignments_1.46.0    jsonlite_2.0.0             
# [55] progressr_0.19.0            ggridges_0.5.7             
# [57] survival_3.8-6              tools_4.5.2                
# [59] TFMPvalue_1.0.0             ica_1.0-3                  
# [61] Rcpp_1.1.1                  glue_1.8.0                 
# [63] gridExtra_2.3               SparseArray_1.10.10        
# [65] MatrixGenerics_1.22.0       withr_3.0.3                
# [67] fastmap_1.2.0               caTools_1.18.3             
# [69] digest_0.6.39               R6_2.6.1                   
# [71] mime_0.13                   scattermore_1.2            
# [73] gtools_3.9.5                tensor_1.5.1               
# [75] dichromat_2.0-0.1           spatstat.data_3.1-9        
# [77] cigarillo_1.0.0             tidyr_1.3.2                
# [79] data.table_1.18.2.1         httr_1.4.8                 
# [81] htmlwidgets_1.6.4           S4Arrays_1.10.1            
# [83] uwot_0.2.4                  pkgconfig_2.0.3            
# [85] gtable_0.3.6                blob_1.3.0                 
# [87] lmtest_0.9-40               S7_0.2.1                   
# [89] htmltools_0.5.9             dotCall64_1.2              
# [91] ProtGenerics_1.42.0         png_0.1-8                  
# [93] spatstat.univar_3.1-6       rstudioapi_0.19.0          
# [95] reshape2_1.4.5              rjson_0.2.23               
# [97] nlme_3.1-168                curl_7.1.0                 
# [99] cachem_1.1.0                zoo_1.8-15                 
# [101] stringr_1.6.0               KernSmooth_2.23-26         
# [103] parallel_4.5.2              miniUI_0.1.2               
# [105] restfulr_0.0.17             pillar_1.11.1              
# [107] grid_4.5.2                  vctrs_0.7.1                
# [109] RANN_2.6.2                  promises_1.5.0             
# [111] xtable_1.8-8                cluster_2.1.8.2            
# [113] cli_3.6.6                   compiler_4.5.2             
# [115] rlang_1.2.0                 crayon_1.5.3               
# [117] future.apply_1.20.2         labeling_0.4.3             
# [119] plyr_1.8.9                  stringi_1.8.7              
# [121] viridisLite_0.4.3           deldir_2.0-4               
# [123] BiocParallel_1.44.0         lazyeval_0.2.2             
# [125] spatstat.geom_3.7-0         Matrix_1.7-4               
# [127] RcppHNSW_0.6.0              sparseMatrixStats_1.22.0   
# [129] bit64_4.8.2                 future_1.70.0              
# [131] KEGGREST_1.50.0             shiny_1.13.0               
# [133] SummarizedExperiment_1.40.0 ROCR_1.0-12                
# [135] igraph_2.2.2                memoise_2.0.1              
# [137] fastmatch_1.1-8             bit_4.6.0

# see 260319_retina_chromatin_d59_d74_d78_merged_no10.R for the next step

