# GSE183684

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(harmony)

RNA_merged <- readRDS("260228_retina_chromatin_RNA_merged.rds")

# RNA_merged
# An object of class Seurat 
# 26242 features across 27819 samples within 1 assay 
# Active assay: RNA (26242 features, 0 variable features)
# 3 layers present: counts.1, counts.2, counts.3

RNA_merged[["percent.mt"]] <- PercentageFeatureSet(RNA_merged, pattern = "^MT-")

VlnPlot(RNA_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

RNA_merged <- subset(RNA_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

RNA_merged <- NormalizeData(RNA_merged, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_merged <- FindVariableFeatures(RNA_merged, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(RNA_merged)
RNA_merged <- ScaleData(RNA_merged, features = all.genes)

RNA_merged <- RunPCA(RNA_merged, features = VariableFeatures(object = RNA_merged), verbose = FALSE)

RNA_merged_harmony <- RunHarmony(RNA_merged, group.by.vars = "dataset")

RNA_merged_harmony <- RunUMAP(RNA_merged_harmony, reduction = "harmony", dims = 1:30)

RNA_merged_harmony <- FindNeighbors(RNA_merged_harmony, reduction = "harmony", dims = 1:30)

# RNA_merged_harmony_08 <- FindClusters(RNA_merged_harmony, resolution = 0.8)
# dim08 <- DimPlot(RNA_merged_harmony_08, reduction = 'umap', group.by = c("dataset", "ident"), ncol = 2)
# # Number of communities: 23
# 
# RNA_merged_harmony_06 <- FindClusters(RNA_merged_harmony, resolution = 0.6)
# dim06 <- DimPlot(RNA_merged_harmony_06, reduction = 'umap', group.by = c("dataset", "ident"), ncol = 2)
# # Number of communities: 21

# RNA_merged_harmony_051 <- FindClusters(RNA_merged_harmony, resolution = 0.51)
# dim051 <- DimPlot(RNA_merged_harmony_051, reduction = 'umap', group.by = c("dataset", "ident"), ncol = 2)
# # Number of communities: 18

# RNA_merged_harmony_04 <- FindClusters(RNA_merged_harmony, resolution = 0.4)
# dim04 <- DimPlot(RNA_merged_harmony_04, reduction = 'umap', group.by = c("dataset", "ident"), ncol = 2)
# # Number of communities: 17
# 
# dim08 | dim06 | dim051 | dim04

RNA_merged_harmony <- FindClusters(RNA_merged_harmony, resolution = 0.51)

# p1 <- DimPlot(RNA_merged_harmony, reduction = 'umap', group.by = c("dataset", "ident"), ncol = 2)
# 
# RNA_merged_umap <- RunUMAP(RNA_merged, dims = 1:30)
# 
# RNA_merged_umap <- FindNeighbors(RNA_merged_umap, dims = 1:30) %>% FindClusters()
# 
# p2 <- DimPlot(RNA_merged_umap, reduction = 'umap', group.by = c("dataset", "ident"), ncol = 2)
# 
# p1 | p2

# https://github.com/TCherryLab/Retina-Multiomic-Analysis/blob/main/Code/human_scRNASeq_seurat.R
# markers <- c('RCVRN', 'RHO', 'CRX', 'ARR3', 'GNAT2', 'VSX2', 'LHX4', 'TRPM1', 'GRM6', 'SLC1A3', 'RLBP1', 'PAX6', 'LHX1', 'ONECUT2', 'TFAP2B', 'GAD1', 'SLC6A9', 'RBPMS', 'NEFM', 'GFAP', 'CD74', 'P2RY12', 'BEST1', 'RPE65', 'SFRP2')
# markers2 <- c('KCNQ5', 'PRDM1', 'DCT', 'SLC35F3', 'NOL4', 'MPPED2', 'WDR72', 'NLK', 'PDE1C', 'PEX5L')

# https://ars.els-cdn.com/content/image/1-s2.0-S1534580722001204-mmc3.xlsx
marker_gene <- c('PROM1','CRX','RCVRN','OTX2','RHO','NR2E3','GNAT1','NRL','GADD45G','NEUROD1',
                 'RXRG','DCT','PRDM1','RAX2','CRABP2','PRPH2','SAG',
                 'ARR3','GNAT2','THRB','OPN1SW','PDE6H','HOTAIRM1',
                 'SLC1A3','SLN','RLBP1','SOX2','NFIA','CRYM','CLU','LINC00461',
                 'VSX1','VSX2','GRM6','PRKCA','LHX4','PROX1','PCP4','PCP2','TRPM1','PRDM8',
                 'ONECUT1','ONECUT2','ONECUT3','TFAP2B','LHX1','TFAP2A','ESRRB',
                 'SLC6A9','GAD1','SLC32A1','GAD2','SLC18A3','LHX9','MEIS2','TFAP2C',
                 'POU4F2','RBPMS','NEFM','GAP43','POU4F1','ELAVL4','POU6F2','ISL1',
                 'NHLH2','EBF1','EBF3','MYC',
                 'BEST1','RPE65','TIMP3',
                 'C8ORF46','ATOH7',
                 'VIM','SFRP2','MKI67','UBE2C','FGF19','CCND1','ID3',
                 'DLX1','DLX2',
                 'ASCL1','SOX4',
                 'OLIG2','NEUROG2','BTG2')

DotPlot(RNA_merged_harmony, features = marker_gene)
# Warning message:
#   The following requested variables were not found: C8ORF46, SLC18A3

new.cluster.ids <- c("Early RPC", "RPC1", "RGC4", 
                     "RPC3", "RPC2", "Horizontal cell2", 
                     "Cone", "RGC3",
                     "RGC1", "PR/BC Precursor", "RGC2", "RPC5", 
                     "Amacrine cell", "RPC4", "Rod", "RGC5", "Horizontal cell1",
                     "Muller glia")
names(new.cluster.ids) <- levels(RNA_merged_harmony)
RNA_merged_harmony <- RenameIdents(RNA_merged_harmony, new.cluster.ids)
RNA_merged_harmony$celltype <- Idents(RNA_merged_harmony)

DimPlot(RNA_merged_harmony, reduction = "umap", label = TRUE, repel = TRUE)

saveRDS(RNA_merged_harmony, "260301_retina_chromatin_merged_No1.rds")

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] future_1.69.0      harmony_1.2.4      Rcpp_1.0.14        ggplot2_4.0.0     
# [5] patchwork_1.3.2    Seurat_5.3.0       SeuratObject_5.0.2 sp_2.2-0          
# [9] dplyr_1.1.4       
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3         
# [4] rlang_1.1.5            magrittr_2.0.3         RcppAnnoy_0.0.22      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7        
# [10] compiler_4.4.1         spatstat.geom_3.3-5    png_0.1-8             
# [13] vctrs_0.6.5            reshape2_1.4.4         stringr_1.6.0         
# [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
# [19] promises_1.5.0         ggbeeswarm_0.7.3       purrr_1.0.4           
# [22] jsonlite_1.9.0         goftest_1.2-3          later_1.4.1           
# [25] spatstat.utils_3.1-2   irlba_2.3.5.1          parallel_4.4.1        
# [28] cluster_2.1.6          R6_2.6.1               ica_1.0-3             
# [31] stringi_1.8.4          RColorBrewer_1.1-3     spatstat.data_3.1-9   
# [34] reticulate_1.41.0      parallelly_1.45.1      spatstat.univar_3.1-1 
# [37] lmtest_0.9-40          scattermore_1.2        tensor_1.5.1          
# [40] future.apply_1.20.2    zoo_1.8-13             sctransform_0.4.1     
# [43] httpuv_1.6.15          Matrix_1.7-3           splines_4.4.1         
# [46] igraph_2.1.4           tidyselect_1.2.1       rstudioapi_0.18.0     
# [49] dichromat_2.0-0.1      abind_1.4-8            spatstat.random_3.3-2 
# [52] codetools_0.2-20       miniUI_0.1.2           spatstat.explore_3.3-4
# [55] listenv_0.10.0         lattice_0.22-6         tibble_3.2.1          
# [58] plyr_1.8.9             withr_3.0.2            shiny_1.13.0          
# [61] S7_0.2.0               ROCR_1.0-12            ggrastr_1.0.2         
# [64] Rtsne_0.17             fastDummies_1.7.5      survival_3.6-4        
# [67] polyclip_1.10-7        fitdistrplus_1.2-6     pillar_1.11.1         
# [70] KernSmooth_2.23-24     plotly_4.12.0          generics_0.1.4        
# [73] RcppHNSW_0.6.0         scales_1.4.0           globals_0.19.0        
# [76] xtable_1.8-8           RhpcBLASctl_0.23-42    glue_1.8.0            
# [79] lazyeval_0.2.2         tools_4.4.1            data.table_1.17.0     
# [82] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2         
# [85] cowplot_1.2.0          grid_4.4.1             tidyr_1.3.1           
# [88] nlme_3.1-164           beeswarm_0.4.0         vipor_0.4.7           
# [91] cli_3.6.4              spatstat.sparse_3.1-0  spam_2.11-1           
# [94] viridisLite_0.4.3      uwot_0.2.3             gtable_0.3.6          
# [97] digest_0.6.37          progressr_0.18.0       ggrepel_0.9.6         
# [100] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1     
# [103] lifecycle_1.0.5        httr_1.4.8             mime_0.12             
# [106] MASS_7.3-60.2