# Watabe et al. F1000Research ver5
# This file is Watabe_F1000_ver5_new_No4.R for result output
# Use 260319_retina_chromatin_d59_d74_d78_merged_No0.R for query scRNAseq import
# Use 260319_retina_chromatin_d59_d74_d78_merged_No1.R for query scRNAseq QC
# Use 260319_retina_chromatin_d59_d74_d78_merged_No2.R for scATACseq QC
# Use 260319_retina_chromatin_d59_d74_d78_merged_No3.R for scATACseq gene activity calculation
# Use Watabe_F1000_ver5_new_No1.R for scRNAseq processing
# Use Watabe_F1000_ver5_new_No2.R for scATACseq processing
# Use Watabe_F1000_ver5_new_No3.R for scATACseq processing

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

# Reference scRNAseq
reference <- readRDS("Kriukov_et_al_scRNAseq_ref.rds")
reference
# An object of class Seurat 
# 25275 features across 108838 samples within 1 assay 
# Active assay: RNA (25275 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 3 dimensional reductions calculated: pca, umap, tsne
Idents(object = reference) <- reference@meta.data$EK_PB_annov1

DotPlot_genes <- c("SOX2", "MKI67", "RLBP1",
                   "RCVRN", "OTX2", "PDE6H", "NRL",
                   "RHO", "TFAP2A", "TFAP2B",
                   "ONECUT2", "ATOH7", "POU4F2", "RBPMS",
                   "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "HES1",
                   "RAX", "PAX6", "VSX2", "LHX2")

levels(reference)
order <- c("Retinal progenitors", "Transitory", "Retinal ganglion cells",
           "Bipolar cells", "Horizontal cells", "Amacrine cells",
           "Cone precursors", "Cones", "Rods",
           "Muller glia", "Astrocytes", "Microglia"
           )
levels(reference) <- order

Figure1A <- DimPlot(reference, reduction = "umap")
# save image as Figure1A_260720.tiff W800 H600
# save image as Figure1A_260725.tiff W600 H450

# From Kriukov et al. Dev. Biol. 2025 Figure 2B
markers <- c("SOX2", "MKI67", "RLBP1", "PAX6",
             "RCVRN", "OTX2", "PDE6H", "NRL",
             "RHO", "VSX2", "TFAP2A", "TFAP2B",
             "ONECUT2", "ATOH7", "POU4F2", "RBPMS")

# Figure 1B: Feature Plot for RPC markers
Figure1B_SOX2 <- FeaturePlot(reference, features = "SOX2")
# save image as Figure1B_SOX2_260720.tiff W600 H600
# save image as Figure1B_SOX2_260725.tiff W450 H450
Figure1B_MKI67 <- FeaturePlot(reference, features = "MKI67")
# save image as Figure1B_MKI67_260720.tiff W600 H600
# save image as Figure1B_MKI67_260725.tiff W450 H450

Figure1C <- DotPlot(reference, 
                    features = DotPlot_genes) + RotatedAxis()
# save image as Figure1C_260725.tiff W1200 H450

Undiff <- WhichCells(reference, idents = c("Retinal progenitors", "Transitory"))
RPC_ref <-WhichCells(reference, idents = "Retinal progenitors")
Transitory_ref <-WhichCells(reference, idents = "Transitory")

Figure1D <- DimPlot(reference, reduction = "umap", 
                    cells.highlight = RPC_ref,
                    cols.highlight = "blue",
                    cols = "grey",
                    raster = FALSE,
                    )
# save image as Figure1D_RPC_260725.tiff W470 H350
Figure1E <- DimPlot(reference, reduction = "umap", 
                    cells.highlight = Transitory_ref,
                    cols.highlight = "blue",
                    cols = "grey",
                    raster = FALSE,
)
# save image as Figure1E_Transitory_260725.tiff W470 H350

Figure1F <- FeaturePlot(reference, features = "NOTCH1")
# save image as Figure1F_NOTCH1_260725.tiff W450 H450
Figure1G <- FeaturePlot(reference, features = "RAX")
# save image as Figure1G_RAX_260725.tiff W450 H450

all_markers_reference_260725 <- FindAllMarkers(object = reference,
                                               only.pos = TRUE)
write.csv(all_markers_reference_260725, "all_markers_reference_260725.csv")
# Additional Table 2

Additional_Table_1 <- table(reference@meta.data$EK_PB_annov1)
write.csv(Additional_Table_1, "Additional_Table_1.csv") # 260729

##### NOT RUN #####
# # Figure 1C: Dot Plot for marker genes
# Figure1C <- DotPlot(reference, features = markers)
# # save image as Figure1C_260720.tiff W2000 H600
# # Figure 1D: Feature Plot for NOTCH1-4
# Figure1D_NOTCH1 <- FeaturePlot(reference, features = "NOTCH1")
# # save image as Figure1D_NOTCH1_260720.tiff W600 H600
# Figure1D_NOTCH2 <- FeaturePlot(reference, features = "NOTCH2")
# # save image as Figure1D_NOTCH2_260720.tiff W600 H600
# Figure1D_NOTCH3 <- FeaturePlot(reference, features = "NOTCH3")
# # save image as Figure1D_NOTCH3_260720.tiff W600 H600
# Figure1D_NOTCH4 <- FeaturePlot(reference, features = "NOTCH4")
# # save image as Figure1D_NOTCH4_260720.tiff W600 H600
# # Figure 1E: Feature Plot for PAX6, LHX2, VSX, RAX
# Figure1E_PAX6 <- FeaturePlot(reference, features = "PAX6")
# # save image as Figure1E_PAX6_260720.tiff W600 H600
# Figure1E_LHX2 <- FeaturePlot(reference, features = "LHX2")
# # save image as Figure1E_LHX2_260720.tiff W600 H600
# Figure1E_VSX2 <- FeaturePlot(reference, features = "VSX2")
# # save image as Figure1E_VSX2_260720.tiff W600 H600
# Figure1E_RAX <- FeaturePlot(reference, features = "RAX")
# # save image as Figure1E_RAX_260720.tiff W600 H600
##### NOT RUN #####

# Query scRNAseq
query <- readRDS("Thomas_et_al_scRNAseq_query.rds")
query
# An object of class Seurat 
# 26242 features across 26388 samples within 1 assay 
# Active assay: RNA (26242 features, 3000 variable features)
# 7 layers present: counts.1, counts.2, counts.3, data.1, data.2, data.3, scale.data
# 3 dimensional reductions calculated: pca, harmony, umap

DotPlot_genes <- c("SOX2", "MKI67", "RLBP1",
                   "RCVRN", "OTX2", "PDE6H", "NRL",
                   "RHO", "TFAP2A", "TFAP2B",
                   "ONECUT2", "ATOH7", "POU4F2", "RBPMS",
                   "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "HES1",
                   "RAX", "PAX6", "VSX2", "LHX2")

levels(query)
order <- c("Retinal progenitors", "Transitory", "Retinal ganglion cells",
           "Bipolar cells", "Horizontal cells", "Amacrine cells",
           "Cone precursors", "Cones", "Rods",
           "Muller glia", "Astrocytes", "Microglia"
)
levels(query) <- order

Figure2A <- DimPlot(query, reduction = "umap")
# save image as Figure2A_260720.tiff W800 H600
# save image as Figure2A_260725.tiff W600 H450

# From Kriukov et al. Dev. Biol. 2025 Figure 2B
# markers <- c("SOX2", "MKI67", "RLBP1", "PAX6",
#              "RCVRN", "OTX2", "PDE6H", "NRL",
#              "RHO", "VSX2", "TFAP2A", "TFAP2B",
#              "ONECUT2", "ATOH7", "POU4F2", "RBPMS")

Figure2B_SOX2 <- FeaturePlot(query, features = "SOX2")
# save image as Figure2B_SOX2_260720.tiff W600 H600
# save image as Figure2B_SOX2_260725.tiff W450 H450
Figure2B_MKI67 <- FeaturePlot(query, features = "MKI67")
# save image as Figure2B_MKI67_260720.tiff W600 H600
# save image as Figure2B_MKI67_260725.tiff W450 H450

Figure2C <- DotPlot(query, 
                    features = DotPlot_genes) + RotatedAxis()
# save image as Figure2C_260725.tiff W1200 H450

RPC_query <-WhichCells(query, idents = "Retinal progenitors")
Transitory_query <-WhichCells(query, idents = "Transitory")

Figure2D <- DimPlot(query, reduction = "umap", 
                    cells.highlight = RPC_query,
                    cols.highlight = "blue",
                    cols = "grey",
                    raster = FALSE,
)
# save image as Figure2D_RPC_260725.tiff W470 H350
Figure2E <- DimPlot(query, reduction = "umap", 
                    cells.highlight = Transitory_query,
                    cols.highlight = "blue",
                    cols = "grey",
                    raster = FALSE,
)
# save image as Figure2E_Transitory_260725.tiff W470 H350

Figure2F <- FeaturePlot(query, features = "NOTCH1")
# save image as Figure2F_NOTCH1_260725.tiff W450 H450
Figure2G <- FeaturePlot(query, features = "RAX")
# save image as Figure2G_RAX_260725.tiff W450 H450

query <- JoinLayers(query)
all_markers_query_260725 <- FindAllMarkers(object = query,
                                           only.pos = TRUE)
write.csv(all_markers_query_260725, "all_markers_query_260725.csv")
# Additional Table 4

Additional_Table_3 <- table(query@active.ident)
write.csv(Additional_Table_3, "Additional_Table_3.csv") # 260729


##### NOT RUN #####
# Figure 2C: Dot Plot for marker genes
# Figure2C <- DotPlot(query, features = markers)
# # save image as Figure2C_260720.tiff W2000 H600
# Figure2D_NOTCH1 <- FeaturePlot(query, features = "NOTCH1")
# # save image as Figure2D_NOTCH1_260720.tiff W600 H600
# Figure2D_NOTCH2 <- FeaturePlot(query, features = "NOTCH2")
# # save image as Figure2D_NOTCH2_260720.tiff W600 H600
# Figure2D_NOTCH3 <- FeaturePlot(query, features = "NOTCH3")
# # save image as Figure2D_NOTCH3_260720.tiff W600 H600
# Figure2D_NOTCH4 <- FeaturePlot(query, features = "NOTCH4")
# # save image as Figure2D_NOTCH4_260720.tiff W600 H600
# # Figure 2E: Feature Plot for PAX6, LHX2, VSX, RAX
# Figure2E_PAX6 <- FeaturePlot(query, features = "PAX6")
# # save image as Figure2E_PAX6_260720.tiff W600 H600
# Figure2E_LHX2 <- FeaturePlot(query, features = "LHX2")
# # save image as Figure2E_LHX2_260720.tiff W600 H600
# Figure2E_VSX2 <- FeaturePlot(query, features = "VSX2")
# # save image as Figure2E_VSX2_260720.tiff W600 H600
# Figure2E_RAX <- FeaturePlot(query, features = "RAX")
# # save image as Figure2E_RAX_260720.tiff W600 H600
##### NOT RUN #####

# scATACseq
ATAC_merged <- readRDS("260719_ATAC_RPC_footprint.rds")
levels(ATAC_merged)
order <- c("Retinal progenitors", "Transitory", "Retinal ganglion cells",
           "Bipolar cells", "Horizontal cells", "Amacrine cells",
           "Cone precursors", "Cones", "Rods",
           "Muller glia", "Astrocytes", "Microglia"
)
levels(ATAC_merged) <- order
Figure3A <- DimPlot(ATAC_merged, reduction = "umap")
# save image as Figure3A_260727.tiff W600 H450
RPC_ATAC <-WhichCells(ATAC_merged, idents = "Retinal progenitors")
Figure3B <- DimPlot(ATAC_merged, reduction = "umap", 
                    cells.highlight = RPC_ATAC,
                    cols.highlight = "blue",
                    cols = "grey",
                    raster = FALSE,
)
# save image as Figure3B_RPC_260727.tiff W470 H350
Transitory_ATAC <-WhichCells(ATAC_merged, idents = "Transitory")
Figure3C <- DimPlot(ATAC_merged, reduction = "umap", 
                    cells.highlight = Transitory_ATAC,
                    cols.highlight = "blue",
                    cols = "grey",
                    raster = FALSE,
)
# save image as Figure3C_Transitory_260727.tiff W470 H350

ATAC_merged_plot <- ATAC_merged
ATAC_merged_plot$dataset <- sub("_.*", "", colnames(ATAC_merged_plot))
DefaultAssay(ATAC_merged_plot) <- 'peaks'
Idents(object = ATAC_merged_plot) <- "predicted.id"
ATAC_merged_plot <- SortIdents(ATAC_merged_plot)
peaks <- granges(ATAC_merged_plot[["peaks"]])
gene_anno <- Annotation(ATAC_merged_plot)
ATAC_merged_plot <- RegionStats(ATAC_merged_plot, genome = BSgenome.Hsapiens.UCSC.hg38)
ATAC_merged_plot <- LinkPeaks(
  object = ATAC_merged_plot,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4")
)
# An object of class Seurat 
# 230592 features across 22231 samples within 2 assays 
# Active assay: peaks (210972 features, 210972 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 2 dimensional reductions calculated: lsi, umap
links <- Links(ATAC_merged_plot)
idents.plot <- order
Idents(ATAC_merged_plot) <- factor(
  ATAC_merged_plot$predicted.id,
  levels = idents.plot
)
p_merged_N1 <- CoveragePlot(
  object = ATAC_merged_plot,
  region = "NOTCH1",
  features = "NOTCH1",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 100000,
  extend.downstream = 100000,
  region.highlight = GRanges(
    seqnames = "chr9",
    ranges = IRanges(
      start = 136541826,
      end = 136542945
    )
  ),
  ymax = 56
)
# save image as NOTCH1_CoveragePlot_260722.tiff W2000 H1000 ymax100
# save image as NOTCH1_CoveragePlot_260725.tiff W2000 H1000 ymax50
# save image as Figure3D_260727.tiff W2000 H1000 ymax50
# save image as Figure3D_260729.tiff W2000 H1000 ymax56 (if no indication 0-170)

##### not run ######
# Figure3E <- CoveragePlot(
#   object = ATAC_merged_plot,
#   region = "chr9-136542636-136542839",
#   expression.assay = "RNA",
#   idents = idents.plot,
#   extend.upstream = 100,
#   extend.downstream = 100,
#   ymax = 50
# )
# save image as Figure3E_260727.tiff W500 H700
##### not run ######

Figure3E <- CoveragePlot(
  object = ATAC_merged_plot,
  region = "chr9-136542126-136542839",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 300,
  extend.downstream = 300,
  region.highlight = GRanges(
    seqnames = "chr9",
    ranges = IRanges(
      start = 136542636,
      end = 136542645
    )
  ),
  ymax = 56
)
# save image as Figure3E_260727.tiff W1000 H700 ymax50
# save image as Figure3E_260729.tiff W1000 H700 ymax56

Figure3F <- PlotFootprint(
  ATAC_merged,
  features = "MA0718.2",
  label.top = 3,
  repel = TRUE
) +
  plot_annotation(
    title = "RAX",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
# save image as Figure3F_260727.tiff W700 H700

notch1_links <- links[links$gene == "NOTCH1"]
as.data.frame(notch1_links)
write.csv(as.data.frame(notch1_links),
          "NOTCH1_peak_zscore_RPC_260722.csv",
          row.names = FALSE)

# Additional Table 6

p_merged_N2 <- CoveragePlot(
  object = ATAC_merged_plot,
  region = "NOTCH2",
  features = "NOTCH2",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 100000,
  extend.downstream = 100000,
  ymax = 50
)
# save image as NOTCH2_CoveragePlot_260722.tiff W2000 H1000 ymax100
# save image as NOTCH2_CoveragePlot_260725.tiff W2000 H1000 ymax50

notch2_links <- links[links$gene == "NOTCH2"]
as.data.frame(notch2_links)
write.csv(as.data.frame(notch2_links),
          "NOTCH2_peak_zscore_RPC_260722.csv",
          row.names = FALSE)

p_merged_N3 <- CoveragePlot(
  object = ATAC_merged_plot,
  region = "NOTCH3",
  features = "NOTCH3",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 100000,
  extend.downstream = 100000,
  ymax = 50
)
# save image as NOTCH3_CoveragePlot_260722.tiff W2000 H1000 ymax100
# save image as NOTCH3_CoveragePlot_260725.tiff W2000 H1000 ymax50
notch3_links <- links[links$gene == "NOTCH3"]
as.data.frame(notch3_links)
write.csv(as.data.frame(notch3_links),
          "NOTCH3_peak_zscore_RPC_260722.csv",
          row.names = FALSE)

p_merged_N4 <- CoveragePlot(
  object = ATAC_merged_plot,
  region = "NOTCH4",
  features = "NOTCH4",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 100000,
  extend.downstream = 100000,
  ymax = 50
)
# save image as NOTCH4_CoveragePlot_260722.tiff W2000 H1000 ymax100
# save image as NOTCH4_CoveragePlot_260725.tiff W2000 H1000 ymax50
notch4_links <- links[links$gene == "NOTCH4"]
as.data.frame(notch4_links)
write.csv(as.data.frame(notch4_links),
          "NOTCH4_peak_zscore_RPC_260722.csv",
          row.names = FALSE)

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
#   [1] ggseqlogo_0.2.2                  
# [2] RSQLite_3.53.2                   
# [3] BSgenome.Hsapiens.UCSC.hg38_1.4.5
# [4] BSgenome_1.78.0                  
# [5] rtracklayer_1.70.1               
# [6] BiocIO_1.20.0                    
# [7] Biostrings_2.78.0                
# [8] XVector_0.50.0                   
# [9] GenomeInfoDb_1.46.2              
# [10] JASPAR2024_0.99.7                
# [11] BiocFileCache_3.0.0              
# [12] dbplyr_2.6.0                     
# [13] scales_1.4.0                     
# [14] GenomicRanges_1.62.1             
# [15] Seqinfo_1.0.0                    
# [16] IRanges_2.44.0                   
# [17] S4Vectors_0.48.1                 
# [18] BiocGenerics_0.56.0              
# [19] generics_0.1.4                   
# [20] Signac_1.17.1                    
# [21] ggplot2_4.0.2                    
# [22] patchwork_1.3.2                  
# [23] Seurat_5.4.0                     
# [24] SeuratObject_5.3.0               
# [25] sp_2.2-1                         
# [26] dplyr_1.2.0                      
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
# [33] RCurl_1.98-1.19             tweenr_2.0.3               
# [35] rappdirs_0.3.4              ggrepel_0.9.7              
# [37] irlba_2.3.7                 listenv_0.10.1             
# [39] spatstat.utils_3.2-1        seqLogo_1.76.0             
# [41] goftest_1.2-3               RSpectra_0.16-2            
# [43] spatstat.random_3.4-4       fitdistrplus_1.2-6         
# [45] parallelly_1.46.1           codetools_0.2-20           
# [47] DelayedArray_0.36.1         RcppRoll_0.3.2             
# [49] ggforce_0.5.0               tidyselect_1.2.1           
# [51] UCSC.utils_1.6.1            farver_2.1.2               
# [53] matrixStats_1.5.0           spatstat.explore_3.7-0     
# [55] GenomicAlignments_1.46.0    jsonlite_2.0.0             
# [57] progressr_0.19.0            ggridges_0.5.7             
# [59] survival_3.8-6              tools_4.5.2                
# [61] TFMPvalue_1.0.0             ica_1.0-3                  
# [63] Rcpp_1.1.1                  glue_1.8.0                 
# [65] gridExtra_2.3               SparseArray_1.10.10        
# [67] MatrixGenerics_1.22.0       withr_3.0.3                
# [69] fastmap_1.2.0               caTools_1.18.3             
# [71] digest_0.6.39               R6_2.6.1                   
# [73] mime_0.13                   scattermore_1.2            
# [75] gtools_3.9.5                tensor_1.5.1               
# [77] dichromat_2.0-0.1           spatstat.data_3.1-9        
# [79] cigarillo_1.0.0             tidyr_1.3.2                
# [81] data.table_1.18.2.1         httr_1.4.8                 
# [83] htmlwidgets_1.6.4           S4Arrays_1.10.1            
# [85] TFBSTools_1.48.0            uwot_0.2.4                 
# [87] pkgconfig_2.0.3             gtable_0.3.6               
# [89] blob_1.3.0                  lmtest_0.9-40              
# [91] S7_0.2.1                    htmltools_0.5.9            
# [93] dotCall64_1.2               Biobase_2.70.0             
# [95] png_0.1-8                   spatstat.univar_3.1-6      
# [97] rstudioapi_0.19.0           reshape2_1.4.5             
# [99] rjson_0.2.23                nlme_3.1-168               
# [101] curl_7.1.0                  cachem_1.1.0               
# [103] zoo_1.8-15                  stringr_1.6.0              
# [105] KernSmooth_2.23-26          parallel_4.5.2             
# [107] miniUI_0.1.2                restfulr_0.0.17            
# [109] pillar_1.11.1               grid_4.5.2                 
# [111] vctrs_0.7.1                 RANN_2.6.2                 
# [113] promises_1.5.0              xtable_1.8-8               
# [115] cluster_2.1.8.2             cli_3.6.6                  
# [117] compiler_4.5.2              Rsamtools_2.26.0           
# [119] rlang_1.2.0                 crayon_1.5.3               
# [121] future.apply_1.20.2         labeling_0.4.3             
# [123] plyr_1.8.9                  stringi_1.8.7              
# [125] viridisLite_0.4.3           deldir_2.0-4               
# [127] BiocParallel_1.44.0         lazyeval_0.2.2             
# [129] spatstat.geom_3.7-0         Matrix_1.7-4               
# [131] RcppHNSW_0.6.0              sparseMatrixStats_1.22.0   
# [133] bit64_4.8.2                 future_1.70.0              
# [135] shiny_1.13.0                SummarizedExperiment_1.40.0
# [137] ROCR_1.0-12                 igraph_2.2.2               
# [139] memoise_2.0.1               fastmatch_1.1-8            
# [141] bit_4.6.0

# hereafter copy-pasted 260319_retina_chromatin_d59_d74_d78_merged_N14.R for the next step

DefaultAssay(ATAC_merged) <- "peaks"

notch1_gr <- LookupGeneCoords(
  object = ATAC_merged,
  gene = "NOTCH1"
)
notch1_ext <- Extend(
  notch1_gr,
  upstream = 100000,
  downstream = 100000
)
genome(notch1_ext) <- "hg38"
strand(notch1_ext) <- "*"
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, notch1_ext)
names(seq) <- paste0(seqnames(notch1_ext), ":", start(notch1_ext), "-", end(notch1_ext)) # This is 260725
# names(seq) <- "NOTCH1" : This is 260722
writeXStringSet(
  DNAStringSet(seq),
  "NOTCH1_100kb_260725.fa"
)

# FIMO (https://meme-suite.org/meme/tools/fimo)
# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
# MA0069.1	PAX6	chr9	136611415	136611428	+	14.3	6.22e-06	1	TTCACACTTGAATT
# MA0069.1	PAX6	chr9	136472750	136472763	+	13.7667	9.13e-06	1	TTCACGCTTAAGCG
# MA0069.1	PAX6	chr9	136447054	136447067	-	11.9444	3.03e-05	1	TTTCTGCATGAGTG
# MA0069.1	PAX6	chr9	136421968	136421981	-	11.9333	3.04e-05	1	CTCACGGATGAGCG
# MA0069.1	PAX6	chr9	136494138	136494151	+	11.4	4.22e-05	1	CTCACGCATGGTGT
# MA0069.1	PAX6	chr9	136447319	136447332	+	11.3778	4.28e-05	1	GTCCCGCTTGGGTT
# MA0069.1	PAX6	chr9	136433775	136433788	-	10.9889	5.42e-05	1	TCCACGCCTGGGTG
# MA0069.1	PAX6	chr9	136603196	136603209	+	10.6222	6.75e-05	1	TTCATGCAACATTT
# MA0069.1	PAX6	chr9	136640563	136640576	+	10.5778	6.91e-05	1	TTCCTGCTTGGCTG
# MA0069.1	PAX6	chr9	136457718	136457731	+	10.4889	7.26e-05	1	TTCACTCATCATCG
# MA0069.1	PAX6	chr9	136473557	136473570	-	10.3111	8.02e-05	1	CACACGCATCAGCT
# MA0069.1	PAX6	chr9	136593277	136593290	-	10.1778	8.66e-05	1	TTCAAGGATGAGGG
# MA0069.1	PAX6	chr9	136527567	136527580	-	9.93333	9.9e-05	1	CTCAAGCCTCACTT
# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 21 2025 at 19:28:56
# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.
# fimo --oc . --verbosity 1 --bgfile --nrdb-- --thresh 1.0E-4 MA0069.1.meme NOTCH1_100kb_260725.fa

# Motif: MA0700.3 LHX2
# There were 0 motif occurences with a p-value less than 0.0001.

# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
# MA0718.1	RAX	chr9	136425119	136425128	+	12.1186	9.62e-06	1	ACCAATTAGC
# MA0718.1	RAX	chr9	136542636	136542645	-	11.8136	1.76e-05	1	GCCAATTAAG
# MA0718.1	RAX	chr9	136637708	136637717	-	11.8136	1.76e-05	1	GCCAATTAAG
# MA0718.1	RAX	chr9	136534411	136534420	+	11.7288	2.31e-05	1	GCCAATTAGA
# MA0718.1	RAX	chr9	136456273	136456282	-	11.5593	3.18e-05	1	GCTAATTAAA
# MA0718.1	RAX	chr9	136560911	136560920	+	11.2712	5.26e-05	1	GTCAATTAAA
# MA0718.1	RAX	chr9	136463520	136463529	+	11.0678	6.71e-05	1	GTCAATTAAG
# MA0718.1	RAX	chr9	136459278	136459287	+	10.8305	9.27e-05	1	GCCAATTATT
# MA0718.1	RAX	chr9	136520246	136520255	-	10.8305	9.27e-05	1	GGCAATTAGT
# MA0718.1	RAX	chr9	136460984	136460993	+	10.8136	9.91e-05	1	GTTAATTAAA
# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 21 2025 at 19:28:56
# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.
# fimo --oc . --verbosity 1 --bgfile --nrdb-- --thresh 1.0E-4 MA0718.1.meme NOTCH1_100kb_260725.fa

# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence
# MA0726.2	VSX2	chr9	136412272	136412278	-	11.6693	9.66e-05	1	CTAATTA
# MA0726.2	VSX2	chr9	136456275	136456281	-	11.6693	9.66e-05	1	CTAATTA
# MA0726.2	VSX2	chr9	136595525	136595531	+	11.6693	9.66e-05	1	CTAATTA
# MA0726.2	VSX2	chr9	136595526	136595532	-	11.6693	9.66e-05	1	CTAATTA
# MA0726.2	VSX2	chr9	136606981	136606987	+	11.6693	9.66e-05	1	CTAATTA
# FIMO (Find Individual Motif Occurrences): Version 5.5.9 compiled on Nov 21 2025 at 19:28:56
# The format of this file is described at https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results.
# fimo --oc . --verbosity 1 --bgfile --nrdb-- --thresh 1.0E-4 MA0726.2.meme NOTCH1_100kb_260725.fa

# Additional Table 7

notch3_gr <- LookupGeneCoords(
  object = ATAC_merged,
  gene = "NOTCH3"
)
notch3_ext <- Extend(
  notch3_gr,
  upstream = 100000,
  downstream = 100000
)
genome(notch3_ext) <- "hg38"
strand(notch3_ext) <- "*"
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, notch3_ext)
names(seq) <- paste0(seqnames(notch3_ext), ":", start(notch3_ext), "-", end(notch3_ext)) # This is 260725
# names(seq) <- "NOTCH3" : This is 260722
writeXStringSet(
  DNAStringSet(seq),
  "NOTCH3_100kb_260725.fa"
)

# FIMO (https://meme-suite.org/meme/tools/fimo)


# library(Signac)
# library(Seurat)
# library(JASPAR2024)
# library(TFBSTools)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(ggplot2)
# library(RSQLite)
# library(ggseqlogo)
# library(motifmatchr)
# library(EnsDb.Hsapiens.v86)
# 
# ATAC_merged_fp <- readRDS("260307_retina_chromatin_ATAC_merged_No6_sub3.rds")

# Footprint analysis
# 260722 Top5
# 260723 Top3
# 260723_2 Top4

p1 <- PlotFootprint(
  ATAC_merged,
  features = "MA0700.3",
  label.top = 3,
  repel = TRUE
) +
  plot_annotation(
    title = "LHX2",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_LHX2_260723.tiff",
  plot = p1,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning message:
#   Removed 9108 rows containing missing values or values outside the scale range (`geom_label_repel()`). 

p2 <- PlotFootprint(
  ATAC_merged,
  features = "MA0069.1",
  label.top = 3,
  repel = TRUE
) +
  plot_annotation(
    title = "PAX6",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_PAX6_260723.tiff",
  plot = p2,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning messages:
# 1: Removed 2 rows containing missing values or values outside the scale range (`geom_line()`). 
# 2: Removed 9252 rows containing missing values or values outside the scale range (`geom_label_repel()`).

p3 <- PlotFootprint(
  ATAC_merged,
  features = "MA0718.2",
  label.top = 3,
  repel = TRUE
) +
  plot_annotation(
    title = "RAX",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_RAX_260723.tiff",
  plot = p3,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning message:
#   Removed 9108 rows containing missing values or values outside the scale range (`geom_label_repel()`).

p4 <- PlotFootprint(
  ATAC_merged,
  features = "MA0726.2",
  label.top = 3,
  repel = TRUE
) +
  plot_annotation(
    title = "VSX2",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_VSX2_260723.tiff",
  plot = p4,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning message:
#   Removed 9126 rows containing missing values or values outside the scale range (`geom_label_repel()`). 

p1 <- PlotFootprint(
  ATAC_merged,
  features = "MA0700.3",
  label.top = 4,
  repel = TRUE
) +
  plot_annotation(
    title = "LHX2",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_LHX2_260723_2.tiff",
  plot = p1,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning message:
#   Removed 9108 rows containing missing values or values outside the scale range (`geom_label_repel()`). 

p2 <- PlotFootprint(
  ATAC_merged,
  features = "MA0069.1",
  label.top = 4,
  repel = TRUE
) +
  plot_annotation(
    title = "PAX6",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_PAX6_260723_2.tiff",
  plot = p2,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning messages:
# 1: Removed 2 rows containing missing values or values outside the scale range (`geom_line()`). 
# 2: Removed 9252 rows containing missing values or values outside the scale range (`geom_label_repel()`).

p3 <- PlotFootprint(
  ATAC_merged,
  features = "MA0718.2",
  label.top = 4,
  repel = TRUE
) +
  plot_annotation(
    title = "RAX",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_RAX_260723_2.tiff",
  plot = p3,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning message:
#   Removed 9108 rows containing missing values or values outside the scale range (`geom_label_repel()`).

p3 <- PlotFootprint(
  ATAC_merged,
  features = "MA0718.2",
  label.top = 3,
  repel = TRUE
) +
  plot_annotation(
    title = "RAX",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )

ggsave(
  filename = "Footprint_RAX_260723.tiff",
  plot = p3,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

p4 <- PlotFootprint(
  ATAC_merged,
  features = "MA0726.2",
  label.top = 4,
  repel = TRUE
) +
  plot_annotation(
    title = "VSX2",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_VSX2_260723_2.tiff",
  plot = p4,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

p1 <- PlotFootprint(
  ATAC_merged,
  features = "MA0700.3",
  label.top = 5,
  repel = TRUE
) +
  plot_annotation(
    title = "LHX2",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_LHX2_260722.tiff",
  plot = p1,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning message:
#   Removed 9108 rows containing missing values or values outside the scale range (`geom_label_repel()`). 

p2 <- PlotFootprint(
  ATAC_merged,
  features = "MA0069.1",
  label.top = 5,
  repel = TRUE
) +
  plot_annotation(
    title = "PAX6",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_PAX6_260722.tiff",
  plot = p2,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning messages:
# 1: Removed 2 rows containing missing values or values outside the scale range (`geom_line()`). 
# 2: Removed 9252 rows containing missing values or values outside the scale range (`geom_label_repel()`).

p3 <- PlotFootprint(
  ATAC_merged,
  features = "MA0718.2",
  label.top = 5,
  repel = TRUE
) +
  plot_annotation(
    title = "RAX",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_RAX_260722.tiff",
  plot = p3,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# Warning message:
#   Removed 9108 rows containing missing values or values outside the scale range (`geom_label_repel()`).

p4 <- PlotFootprint(
  ATAC_merged,
  features = "MA0726.2",
  label.top = 5,
  repel = TRUE
) +
  plot_annotation(
    title = "VSX2",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.45)
    )
  )
ggsave(
  filename = "Footprint_VSX2_260722.tiff",
  plot = p4,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
 

motifplot <- MotifPlot(
  object = ATAC_merged,
  motifs = c("LHX2", "PAX6", "RAX", "VSX2")
)
ggsave(
  filename = "MotifPlot_260723.tiff",
  plot = motifplot,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

table(ATAC_merged$celltype)
# Retinal ganglion cells         Amacrine cells 
# 6373                  10536 
# Transitory    Retinal progenitors 
# 1217                   1377 
# Cones       Horizontal cells 
# 97                    729 
# Rods        Cone precursors 
# 47                   1465 
# Astrocytes            Muller glia 
# 36                    222 
# Bipolar cells              Microglia 
# 125                      7 
ATAC_merged$dataset <- sub("_.*", "", colnames(ATAC_merged))
table(ATAC_merged$celltype, ATAC_merged$dataset)
# d59  d74  d78
# Retinal ganglion cells 2393 1633 2347
# Amacrine cells         2634 2831 5071
# Transitory              480  217  520
# Retinal progenitors     645  230  502
# Cones                     2   21   74
# Horizontal cells        330  112  287
# Rods                      2    9   36
# Cone precursors         329  353  783
# Astrocytes               33    0    3
# Muller glia              43   27  152
# Bipolar cells            24   28   73
# Microglia                 6    0    1

# Additional Table 5

all(
  rowSums(table(ATAC_merged$celltype, ATAC_merged$dataset)) ==
    table(ATAC_merged$celltype)
)
# [1] TRUE

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
#   [1] TFBSTools_1.48.0                 
# [2] ggseqlogo_0.2.2                  
# [3] RSQLite_3.53.2                   
# [4] BSgenome.Hsapiens.UCSC.hg38_1.4.5
# [5] BSgenome_1.78.0                  
# [6] rtracklayer_1.70.1               
# [7] BiocIO_1.20.0                    
# [8] Biostrings_2.78.0                
# [9] XVector_0.50.0                   
# [10] GenomeInfoDb_1.46.2              
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
# [33] RCurl_1.98-1.19             tweenr_2.0.3               
# [35] rappdirs_0.3.4              ggrepel_0.9.7              
# [37] irlba_2.3.7                 listenv_0.10.1             
# [39] spatstat.utils_3.2-1        seqLogo_1.76.0             
# [41] goftest_1.2-3               RSpectra_0.16-2            
# [43] spatstat.random_3.4-4       fitdistrplus_1.2-6         
# [45] parallelly_1.46.1           codetools_0.2-20           
# [47] DelayedArray_0.36.1         RcppRoll_0.3.2             
# [49] ggforce_0.5.0               tidyselect_1.2.1           
# [51] UCSC.utils_1.6.1            farver_2.1.2               
# [53] matrixStats_1.5.0           spatstat.explore_3.7-0     
# [55] GenomicAlignments_1.46.0    jsonlite_2.0.0             
# [57] progressr_0.19.0            ggridges_0.5.7             
# [59] survival_3.8-6              systemfonts_1.3.2          
# [61] tools_4.5.2                 ragg_1.5.2                 
# [63] TFMPvalue_1.0.0             ica_1.0-3                  
# [65] Rcpp_1.1.1                  glue_1.8.0                 
# [67] gridExtra_2.3               SparseArray_1.10.10        
# [69] MatrixGenerics_1.22.0       withr_3.0.3                
# [71] fastmap_1.2.0               caTools_1.18.3             
# [73] digest_0.6.39               R6_2.6.1                   
# [75] mime_0.13                   textshaping_1.0.5          
# [77] scattermore_1.2             gtools_3.9.5               
# [79] tensor_1.5.1                dichromat_2.0-0.1          
# [81] spatstat.data_3.1-9         cigarillo_1.0.0            
# [83] tidyr_1.3.2                 data.table_1.18.2.1        
# [85] httr_1.4.8                  htmlwidgets_1.6.4          
# [87] S4Arrays_1.10.1             uwot_0.2.4                 
# [89] pkgconfig_2.0.3             gtable_0.3.6               
# [91] blob_1.3.0                  lmtest_0.9-40              
# [93] S7_0.2.1                    htmltools_0.5.9            
# [95] dotCall64_1.2               Biobase_2.70.0             
# [97] png_0.1-8                   spatstat.univar_3.1-6      
# [99] rstudioapi_0.19.0           reshape2_1.4.5             
# [101] rjson_0.2.23                nlme_3.1-168               
# [103] curl_7.1.0                  cachem_1.1.0               
# [105] zoo_1.8-15                  stringr_1.6.0              
# [107] KernSmooth_2.23-26          parallel_4.5.2             
# [109] miniUI_0.1.2                restfulr_0.0.17            
# [111] pillar_1.11.1               grid_4.5.2                 
# [113] vctrs_0.7.1                 RANN_2.6.2                 
# [115] promises_1.5.0              xtable_1.8-8               
# [117] cluster_2.1.8.2             cli_3.6.6                  
# [119] compiler_4.5.2              Rsamtools_2.26.0           
# [121] rlang_1.2.0                 crayon_1.5.3               
# [123] future.apply_1.20.2         labeling_0.4.3             
# [125] plyr_1.8.9                  stringi_1.8.7              
# [127] viridisLite_0.4.3           deldir_2.0-4               
# [129] BiocParallel_1.44.0         lazyeval_0.2.2             
# [131] spatstat.geom_3.7-0         Matrix_1.7-4               
# [133] RcppHNSW_0.6.0              sparseMatrixStats_1.22.0   
# [135] bit64_4.8.2                 future_1.70.0              
# [137] shiny_1.13.0                SummarizedExperiment_1.40.0
# [139] ROCR_1.0-12                 igraph_2.2.2               
# [141] memoise_2.0.1               fastmatch_1.1-8            
# [143] bit_4.6.0

#######################################

# RNA_merged <- readRDS("260301_retina_chromatin_RNA_merged_No4.rds")
# table(RNA_merged$celltype)
# Early RPC             RPC1             RGC4             RPC3 
# 4487             4371             2600             2159 
# RPC2 Horizontal cell2             Cone             RGC3 
# 2118             1914             1704             1362 
# RGC1  PR/BC Precursor             RGC2             RPC5 
# 1227             1053              731              693 
# Amacrine cell             RPC4              Rod             RGC5 
# 635              463              450              190 
# Horizontal cell1      Muller glia 
# 150               81

# table(RNA_merged$celltype, RNA_merged$dataset)
# day 59 day 74 day 78
# Early RPC          1316   1293   1878
# RPC1               1301   1801   1269
# RGC4                746   1232    622
# RPC3                577    972    610
# RPC2                474    925    719
# Horizontal cell2    489    765    660
# Cone                334    761    609
# RGC3                446    588    328
# RGC1                604    459    164
# PR/BC Precursor     250    452    351
# RGC2                332    267    132
# RPC5                199    295    199
# Amacrine cell        55    242    338
# RPC4                133    201    129
# Rod                  34    185    231
# RGC5                 65    101     24
# Horizontal cell1     32     73     45
# Muller glia          69      7      5

# all(
#   rowSums(table(RNA_merged$celltype, RNA_merged$dataset)) ==
#     table(RNA_merged$celltype)
# )
# [1] TRUE