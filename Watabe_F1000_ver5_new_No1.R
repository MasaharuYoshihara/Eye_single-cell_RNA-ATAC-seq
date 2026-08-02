# Watabe et al. F1000Research ver5
# This file is Watabe_F1000_ver5_new_No1.R for scRNAseq processing

# Original data
# Reference scRNAseq from Kriukov et al. Dev. Biol. 2025
# Query scRNAseq from Thomas et al. Dev Cell. 2022
# scATACseq from Thomas et al. Dev Cell. 2022

# Milestones
# Processed Reference scRNAseq
# Code A
# RDS file: Kriukov_et_al_scRNAseq_ref.rds

# Processed Query scRNAseq
# Code B
# RDS file: Thomas_et_al_scRNAseq_query.rds

# This file (Watabe_F1000_ver5_new_No1.R) is for the codes A and B.
# This process finished on 20260709.

# scATACseq analysis starts from Watabe_F1000_ver5_new_No2.R

##### Code A #####

library(anndata)
library(Seurat)
library(Matrix)

# 1. Load the h5ad file into R
adata <- read_h5ad("860a9839-5d24-4073-9a67-6ad570f41da1.h5ad")

# 2. Extract and transpose the sparse expression matrix
# AnnData maps rows=cells, cols=genes. Seurat maps rows=genes, cols=cells.
counts_matrix <- Matrix::t(adata$X)

# 3. Assign cell and gene names to the matrix
colnames(counts_matrix) <- adata$obs_names
rownames(counts_matrix) <- adata$var_names

# 4. Create the Seurat object using the counts matrix and metadata
seurat_obj <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = adata$obs
)

# 5. Optional: Add dimensional reductions (e.g., UMAP) if present
if ("X_umap" %in% names(adata$obsm)) {
  umap_coords <- adata$obsm$X_umap
  rownames(umap_coords) <- adata$obs_names
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  seurat_obj[["umap"]] <- CreateDimReducObject(
    embeddings = umap_coords, 
    key = "UMAP_", 
    assay = "RNA"
  )
}

# 6. Save the structured object as an RDS file
saveRDS(seurat_obj, file = "Kriukov_et_al_scRNAseq.rds")

# An object of class Seurat 
# 32055 features across 108838 samples within 1 assay 
# Active assay: RNA (32055 features, 0 variable features)
# 1 layer present: counts
# 1 dimensional reduction calculated: umap

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
#   [1] stats     graphics  grDevices utils     datasets 
# [6] methods   base     
# 
# other attached packages:
#   [1] Matrix_1.7-4       Seurat_5.4.0      
# [3] SeuratObject_5.3.0 sp_2.2-1          
# [5] anndata_0.8.0      remotes_2.5.0     
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4         
# [3] gridExtra_2.3          rlang_1.2.0           
# [5] magrittr_2.0.4         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0     
# [9] ggridges_0.5.7         compiler_4.5.2        
# [11] spatstat.geom_3.7-0    png_0.1-8             
# [13] vctrs_0.7.1            reshape2_1.4.5        
# [15] stringr_1.6.0          pkgconfig_2.0.3       
# [17] fastmap_1.2.0          promises_1.5.0        
# [19] purrr_1.2.1            jsonlite_2.0.0        
# [21] goftest_1.2-3          later_1.4.8           
# [23] spatstat.utils_3.2-1   irlba_2.3.7           
# [25] parallel_4.5.2         cluster_2.1.8.2       
# [27] R6_2.6.1               ica_1.0-3             
# [29] stringi_1.8.7          RColorBrewer_1.1-3    
# [31] spatstat.data_3.1-9    reticulate_1.45.0     
# [33] parallelly_1.46.1      spatstat.univar_3.1-6 
# [35] lmtest_0.9-40          scattermore_1.2       
# [37] Rcpp_1.1.1             assertthat_0.2.1      
# [39] tensor_1.5.1           future.apply_1.20.2   
# [41] zoo_1.8-15             sctransform_0.4.3     
# [43] httpuv_1.6.16          splines_4.5.2         
# [45] igraph_2.2.2           tidyselect_1.2.1      
# [47] rstudioapi_0.18.0      abind_1.4-8           
# [49] spatstat.random_3.4-4  codetools_0.2-20      
# [51] miniUI_0.1.2           spatstat.explore_3.7-0
# [53] curl_7.0.0             listenv_0.10.1        
# [55] lattice_0.22-9         tibble_3.3.1          
# [57] plyr_1.8.9             withr_3.0.2           
# [59] shiny_1.13.0           S7_0.2.1              
# [61] ROCR_1.0-12            Rtsne_0.17            
# [63] future_1.69.0          fastDummies_1.7.5     
# [65] survival_3.8-6         polyclip_1.10-7       
# [67] fitdistrplus_1.2-6     pillar_1.11.1         
# [69] KernSmooth_2.23-26     plotly_4.12.0         
# [71] generics_0.1.4         rprojroot_2.1.1       
# [73] RcppHNSW_0.6.0         ggplot2_4.0.2         
# [75] scales_1.4.0           globals_0.19.0        
# [77] xtable_1.8-8           glue_1.8.0            
# [79] lazyeval_0.2.2         tools_4.5.2           
# [81] data.table_1.18.2.1    RSpectra_0.16-2       
# [83] RANN_2.6.2             dotCall64_1.2         
# [85] cowplot_1.2.0          grid_4.5.2            
# [87] tidyr_1.3.2            nlme_3.1-168          
# [89] patchwork_1.3.2        cli_3.6.6             
# [91] rappdirs_0.3.4         spatstat.sparse_3.1-0 
# [93] spam_2.11-3            viridisLite_0.4.3     
# [95] dplyr_1.2.0            uwot_0.2.4            
# [97] gtable_0.3.6           digest_0.6.39         
# [99] progressr_0.18.0       ggrepel_0.9.7         
# [101] htmlwidgets_1.6.4      farver_2.1.2          
# [103] htmltools_0.5.9        lifecycle_1.0.5       
# [105] here_1.0.2             httr_1.4.8            
# [107] mime_0.13              MASS_7.3-65

library(Seurat)

ref <- readRDS("Kriukov_et_al_scRNAseq.rds")
# ref
# An object of class Seurat 
# 32055 features across 108838 samples within 1 assay 
# Active assay: RNA (32055 features, 0 variable features)
# 1 layer present: counts
# 1 dimensional reduction calculated: umap

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast") # Alternatives: "uswest", "asia"
ensembl_ids <- rownames(ref)
gene_mapping <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name"),
  values = ensembl_ids,
  mart = ensembl
)
print(gene_mapping)

gene_mapping <- gene_mapping[gene_mapping$external_gene_name != "", ]
gene_mapping <- gene_mapping[!duplicated(gene_mapping$ensembl_gene_id), ]
gene_vector <- gene_mapping$external_gene_name
names(gene_vector) <- gene_mapping$ensembl_gene_id

genes_to_keep <- intersect(rownames(ref), names(gene_vector))
ref <- subset(ref, features = genes_to_keep)
counts_matrix <- GetAssayData(ref, layer = "counts")
new_gene_names <- gene_vector[rownames(counts_matrix)]
rownames(counts_matrix) <- new_gene_names
any(duplicated(rownames(counts_matrix)))
rownames(counts_matrix) <- make.unique(rownames(counts_matrix))
ref_rename <- CreateSeuratObject(counts = counts_matrix, meta.data = ref@meta.data)
# ref_rename
# An object of class Seurat 
# 25275 features across 108838 samples within 1 assay 
# Active assay: RNA (25275 features, 0 variable features)
# 1 layer present: counts

saveRDS(ref_rename, "Kriukov_et_al_scRNAseq_rename.rds")

sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin24.4.0
# Running under: macOS Tahoe 26.5.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /opt/homebrew/Cellar/r/4.5.1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] biomaRt_2.66.2      BiocManager_1.30.27 Seurat_5.5.1        SeuratObject_5.4.0 
# [5] sp_2.2-1           
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.19.0      jsonlite_2.0.0        
# [4] magrittr_2.0.5         spatstat.utils_3.2-3   farver_2.1.2          
# [7] vctrs_0.7.3            ROCR_1.0-12            memoise_2.0.1         
# [10] spatstat.explore_3.8-1 htmltools_0.5.9        progress_1.2.3        
# [13] curl_7.1.0             sctransform_0.4.3      parallelly_1.47.0     
# [16] KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3             
# [19] plyr_1.8.9             httr2_1.2.3            plotly_4.12.0         
# [22] zoo_1.8-15             cachem_1.1.0           igraph_2.3.3          
# [25] mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3       
# [28] Matrix_1.7-5           R6_2.6.1               fastmap_1.2.0         
# [31] fitdistrplus_1.2-6     future_1.70.0          shiny_1.14.0          
# [34] digest_0.6.39          S4Vectors_0.48.1       patchwork_1.3.2       
# [37] AnnotationDbi_1.72.0   tensor_1.5.1           RSpectra_0.16-2       
# [40] irlba_2.3.7            RSQLite_3.53.2         filelock_1.0.3        
# [43] progressr_0.19.0       spatstat.sparse_3.2-0  httr_1.4.8            
# [46] polyclip_1.10-7        abind_1.4-8            compiler_4.5.1        
# [49] withr_3.0.3            bit64_4.8.2            S7_0.2.2              
# [52] DBI_1.3.0              fastDummies_1.7.6      MASS_7.3-65           
# [55] rappdirs_0.3.4         tools_4.5.1            lmtest_0.9-40         
# [58] otel_0.2.0             httpuv_1.6.17          future.apply_1.20.2   
# [61] goftest_1.2-3          glue_1.8.1             nlme_3.1-169          
# [64] promises_1.5.0         grid_4.5.1             Rtsne_0.17            
# [67] cluster_2.1.8.2        reshape2_1.4.5         generics_0.1.4        
# [70] gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
# [73] data.table_1.18.4      hms_1.1.4              xml2_1.6.0            
# [76] XVector_0.50.0         BiocGenerics_0.56.0    spatstat.geom_3.8-1   
# [79] RcppAnnoy_0.0.23       ggrepel_0.9.8          RANN_2.6.2            
# [82] pillar_1.11.1          stringr_1.6.0          spam_2.11-4           
# [85] RcppHNSW_0.7.0         later_1.4.8            splines_4.5.1         
# [88] dplyr_1.2.1            BiocFileCache_3.0.0    lattice_0.22-9        
# [91] survival_3.8-6         bit_4.6.0              deldir_2.0-4          
# [94] tidyselect_1.2.1       Biostrings_2.78.0      miniUI_0.1.2          
# [97] pbapply_1.7-4          gridExtra_2.3.1        Seqinfo_1.0.0         
# [100] IRanges_2.44.0         scattermore_1.2        stats4_4.5.1          
# [103] Biobase_2.70.0         matrixStats_1.5.0      stringi_1.8.7         
# [106] lazyeval_0.2.3         codetools_0.2-20       tibble_3.3.1          
# [109] cli_3.6.6              uwot_0.2.4             xtable_1.8-8          
# [112] reticulate_1.46.0      Rcpp_1.1.1-1.1         globals_0.19.1        
# [115] spatstat.random_3.5-0  dbplyr_2.6.0           png_0.1-9             
# [118] spatstat.univar_3.2-0  parallel_4.5.1         ggplot2_4.0.3         
# [121] blob_1.3.0             prettyunits_1.2.0      dotCall64_1.2         
# [124] listenv_1.0.0          viridisLite_0.4.3      scales_1.4.0          
# [127] ggridges_0.5.7         purrr_1.2.2            crayon_1.5.3          
# [130] rlang_1.2.0            KEGGREST_1.50.0        cowplot_1.2.0

library(Seurat)

Seurat <- readRDS("Kriukov_et_al_scRNAseq_rename.rds")

# From https://github.com/mcrewcow/fetal_retina_Kriukov/blob/main/main.R
Seurat <- NormalizeData(Seurat)
Seurat <- FindVariableFeatures(Seurat, selection.method = "vst", nfeatures = 3000)
Seurat <- ScaleData(Seurat) #could be replaced with SCTransform
Seurat <- RunPCA(Seurat, npcs = 100)
Seurat <- FindNeighbors(Seurat, dims = 1:100)
Seurat <- FindClusters(Seurat, resolution = 1)
Seurat <- RunUMAP(Seurat, dims = 1:100)
Seurat <- RunTSNE(Seurat,  dims.use = 1:100)
DimPlot(object = Seurat, reduction = "umap")

print(Seurat@meta.data$cell_type)
# 9 Levels

saveRDS(Seurat, "Kriukov_et_al_scRNAseq_ref.rds")

sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin24.4.0
# Running under: macOS Tahoe 26.5.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /opt/homebrew/Cellar/r/4.5.1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] future_1.70.0      Seurat_5.5.1       SeuratObject_5.4.0 sp_2.2-1          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3.1       
# [4] rlang_1.2.0            magrittr_2.0.5         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7        
# [10] compiler_4.5.1         spatstat.geom_3.8-1    png_0.1-9             
# [13] vctrs_0.7.3            reshape2_1.4.5         stringr_1.6.0         
# [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
# [19] promises_1.5.0         purrr_1.2.2            jsonlite_2.0.0        
# [22] goftest_1.2-3          later_1.4.8            spatstat.utils_3.2-3  
# [25] irlba_2.3.7            parallel_4.5.1         cluster_2.1.8.2       
# [28] R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [31] RColorBrewer_1.1-3     spatstat.data_3.1-9    reticulate_1.46.0     
# [34] parallelly_1.47.0      spatstat.univar_3.2-0  lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.1.1-1.1         tensor_1.5.1          
# [40] future.apply_1.20.2    zoo_1.8-15             sctransform_0.4.3     
# [43] httpuv_1.6.17          Matrix_1.7-5           splines_4.5.1         
# [46] igraph_2.3.3           tidyselect_1.2.1       rstudioapi_0.19.0     
# [49] abind_1.4-8            spatstat.random_3.5-0  codetools_0.2-20      
# [52] miniUI_0.1.2           spatstat.explore_3.8-1 listenv_1.0.0         
# [55] lattice_0.22-9         tibble_3.3.1           plyr_1.8.9            
# [58] withr_3.0.3            shiny_1.14.0           S7_0.2.2              
# [61] ROCR_1.0-12            Rtsne_0.17             fastDummies_1.7.6     
# [64] survival_3.8-6         polyclip_1.10-7        fitdistrplus_1.2-6    
# [67] pillar_1.11.1          KernSmooth_2.23-26     plotly_4.12.0         
# [70] generics_0.1.4         RcppHNSW_0.7.0         ggplot2_4.0.3         
# [73] scales_1.4.0           globals_0.19.1         xtable_1.8-8          
# [76] glue_1.8.1             lazyeval_0.2.3         tools_4.5.1           
# [79] data.table_1.18.4      RSpectra_0.16-2        RANN_2.6.2            
# [82] dotCall64_1.2          cowplot_1.2.0          grid_4.5.1            
# [85] tidyr_1.3.2            nlme_3.1-169           patchwork_1.3.2       
# [88] cli_3.6.6              spatstat.sparse_3.2-0  spam_2.11-4           
# [91] viridisLite_0.4.3      dplyr_1.2.1            uwot_0.2.4            
# [94] gtable_0.3.6           digest_0.6.39          progressr_0.19.0      
# [97] ggrepel_0.9.8          htmlwidgets_1.6.4      farver_2.1.2          
# [100] htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8            
# [103] mime_0.13              MASS_7.3-65

library(Seurat)
library(ggplot2)

ref <- readRDS("Kriukov_et_al_scRNAseq_ref.rds")
# ref
# An object of class Seurat 
# 25275 features across 108838 samples within 1 assay 
# Active assay: RNA (25275 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 3 dimensional reductions calculated: pca, umap, tsne

print(ref@meta.data$EK_PB_annov1)
# [1] Retinal ganglion cells Retinal ganglion cells Transitory            
# [4] Retinal ganglion cells Transitory             Transitory
# 12 Levels: Muller glia Rods Retinal ganglion cells ... Microglia

Idents(object = ref) <- ref@meta.data$EK_PB_annov1

DimPlot(object = ref, reduction = "umap")
# save image as "DimPlot_ref_scRNAseq_20260629.tiff" W800 H600

FeaturePlot(ref, features = "PAX6")
# save image as "FeaturePlot_PAX6_ref_scRNAseq_20260629.tiff" W650 H600

FeaturePlot(ref, features = "NOTCH1")
# save image as "FeaturePlot_NOTCH1_ref_scRNAseq_20260629.tiff" W650 H600

DotPlot(ref, features = c("PAX6", "NOTCH1"))
# save image as "DotPlot_PAX6_NOTCH1_ref_scRNAseq_20260629.tiff" W650 H600

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggplot2_4.0.2      Seurat_5.4.0       SeuratObject_5.3.0
# [4] sp_2.2-1          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4          gridExtra_2.3         
# [4] rlang_1.2.0            magrittr_2.0.4         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0      ggridges_0.5.7        
# [10] compiler_4.5.2         spatstat.geom_3.7-0    png_0.1-8             
# [13] vctrs_0.7.1            reshape2_1.4.5         stringr_1.6.0         
# [16] pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3        
# [19] promises_1.5.0         purrr_1.2.1            jsonlite_2.0.0        
# [22] goftest_1.2-3          later_1.4.8            spatstat.utils_3.2-1  
# [25] irlba_2.3.7            parallel_4.5.2         cluster_2.1.8.2       
# [28] R6_2.6.1               ica_1.0-3              stringi_1.8.7         
# [31] RColorBrewer_1.1-3     spatstat.data_3.1-9    reticulate_1.45.0     
# [34] parallelly_1.46.1      spatstat.univar_3.1-6  lmtest_0.9-40         
# [37] scattermore_1.2        Rcpp_1.1.1             tensor_1.5.1          
# [40] future.apply_1.20.2    zoo_1.8-15             sctransform_0.4.3     
# [43] httpuv_1.6.16          Matrix_1.7-4           splines_4.5.2         
# [46] igraph_2.2.2           tidyselect_1.2.1       rstudioapi_0.19.0     
# [49] abind_1.4-8            spatstat.random_3.4-4  codetools_0.2-20      
# [52] miniUI_0.1.2           spatstat.explore_3.7-0 listenv_0.10.1        
# [55] lattice_0.22-9         tibble_3.3.1           plyr_1.8.9            
# [58] withr_3.0.3            shiny_1.13.0           S7_0.2.1              
# [61] ROCR_1.0-12            Rtsne_0.17             future_1.70.0         
# [64] fastDummies_1.7.6      survival_3.8-6         polyclip_1.10-7       
# [67] fitdistrplus_1.2-6     pillar_1.11.1          KernSmooth_2.23-26    
# [70] plotly_4.12.0          generics_0.1.4         RcppHNSW_0.6.0        
# [73] scales_1.4.0           globals_0.19.1         xtable_1.8-8          
# [76] glue_1.8.0             lazyeval_0.2.2         tools_4.5.2           
# [79] data.table_1.18.2.1    RSpectra_0.16-2        RANN_2.6.2            
# [82] dotCall64_1.2          cowplot_1.2.0          grid_4.5.2            
# [85] tidyr_1.3.2            nlme_3.1-168           patchwork_1.3.2       
# [88] cli_3.6.6              spatstat.sparse_3.1-0  spam_2.11-3           
# [91] viridisLite_0.4.3      dplyr_1.2.0            uwot_0.2.4            
# [94] gtable_0.3.6           digest_0.6.39          progressr_0.19.0      
# [97] ggrepel_0.9.7          htmlwidgets_1.6.4      farver_2.1.2          
# [100] htmltools_0.5.9        lifecycle_1.0.5        httr_1.4.8            
# [103] mime_0.13              MASS_7.3-65

##### Code B #####
library(Seurat)
library(harmony)

ref <- readRDS("Kriukov_et_al_scRNAseq_ref.rds")
Idents(object = ref) <- ref@meta.data$EK_PB_annov1
DimPlot(object = ref, reduction = "umap")

d59_RNA_count <- ReadMtx(mtx = "GSM5567526_d59_matrix.mtx.gz",
                         features = "GSM5567526_d59_features.tsv.gz",
                         cells = "GSM5567526_d59_barcodes.tsv.gz")
d74_RNA_count <- ReadMtx(mtx = "GSM5567527_d74_matrix.mtx.gz",
                         features = "GSM5567527_d74_features.tsv.gz",
                         cells = "GSM5567527_d74_barcodes.tsv.gz")
d78_RNA_count <- ReadMtx(mtx = "GSM5567528_d78_matrix.mtx.gz",
                         features = "GSM5567528_d78_features.tsv.gz",
                         cells = "GSM5567528_d78_barcodes.tsv.gz")

d59_RNA <- CreateSeuratObject(counts = d59_RNA_count,
                              min.cells = 3, min.features = 200)
d74_RNA <- CreateSeuratObject(counts = d74_RNA_count,
                              min.cells = 3, min.features = 200)
d78_RNA <- CreateSeuratObject(counts = d78_RNA_count,
                              min.cells = 3, min.features = 200)

d59_RNA$Sample <- 'd59'
d74_RNA$Sample <- 'd74'
d78_RNA$Sample <- 'd78'

RNA_merged <- merge(d59_RNA, y = c(d74_RNA, d78_RNA),
                    add.cell.ids = c("d59", "d74", "d78"), project = "query_scRNAseq")
head(RNA_merged$Sample)
# d59_AAACCCACAAAGCTAA-1 d59_AAACCCACACTATCGA-1 
# "d59"                  "d59" 
# d59_AAACCCAGTACCCAGC-1 d59_AAACCCAGTCTACAGT-1 
# "d59"                  "d59" 
# d59_AAACCCAGTCTTGGTA-1 d59_AAACCCAGTGACTAAA-1 
# "d59"                  "d59" 

RNA_merged[["percent.mt"]] <- PercentageFeatureSet(RNA_merged, pattern = "^MT-")
VlnPlot(RNA_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RNA_merged <- subset(RNA_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
RNA_merged <- NormalizeData(RNA_merged, normalization.method = "LogNormalize", scale.factor = 10000)
RNA_merged <- FindVariableFeatures(RNA_merged, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(RNA_merged)
RNA_merged <- ScaleData(RNA_merged, features = all.genes)
RNA_merged <- RunPCA(RNA_merged, features = VariableFeatures(object = RNA_merged), verbose = FALSE)
RNA_merged_harmony <- RunHarmony(RNA_merged, group.by.vars = "Sample")
Idents(RNA_merged) <- RNA_merged$Sample
Idents(RNA_merged)
DimPlot(RNA_merged)
# save as DimPlot_RNAmerged_260709.tiff W650 H600
Idents(RNA_merged_harmony) <- RNA_merged_harmony$Sample
Idents(RNA_merged_harmony)
DimPlot(RNA_merged_harmony)
# save as DimPlot_RNAmerged_harmony_260709.tiff W650 H600
RNA_merged_harmony <- RunUMAP(RNA_merged_harmony, reduction = "harmony", dims = 1:30)
RNA_merged_harmony <- FindNeighbors(RNA_merged_harmony, reduction = "harmony", dims = 1:30)
RNA_merged_harmony <- FindClusters(RNA_merged_harmony, resolution = 0.50)
DimPlot(RNA_merged_harmony, reduction = "harmony")
# Number of communities: 20
# save as DimPlot_RNAmerged_harmony_clustering_260709.tiff W650 H600
DimPlot(RNA_merged_harmony, reduction = "umap")
# save as DimPlot_RNAmerged_harmony_umap_260709.tiff W650 H600

common_features <- intersect(VariableFeatures(ref), VariableFeatures(RNA_merged_harmony))
anchors <- FindTransferAnchors(
  reference = ref,
  query = RNA_merged_harmony,
  reference.reduction = "pca",
  dims = 1:30,
  features = common_features
)
pred <- TransferData(
  anchorset = anchors,
  refdata = ref@meta.data$EK_PB_annov1
)
RNA_merged_harmony <- AddMetaData(RNA_merged_harmony, pred)

head(RNA_merged_harmony@meta.data$predicted.id)

Idents(object = RNA_merged_harmony) <- RNA_merged_harmony@meta.data$predicted.id
DimPlot(RNA_merged_harmony, reduction = "umap")
# save image as "DimPlot_query_scRNAseq_260709.tiff" W800 H600
FeaturePlot(RNA_merged_harmony, features = "NOTCH1")
# save image as "FeaturePlot_NOTCH1_query_scRNAseq_260709.tiff" W600 H600
FeaturePlot(RNA_merged_harmony, features = "PAX6")
# save image as "FeaturePlot_PAX6_query_scRNAseq_260709.tiff" W600 H600
DotPlot(RNA_merged_harmony, features = c("PAX6", "NOTCH1"))
# save image as "DotPlot_PAX6_NOTCH1_query_scRNAseq_260709.tiff" W600 H600

saveRDS(RNA_merged_harmony, "Thomas_et_al_scRNAseq_query.rds")
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
#   [1] stats     graphics  grDevices utils     datasets 
# [6] methods   base     
# 
# other attached packages:
#   [1] future_1.70.0      harmony_1.2.4      Rcpp_1.1.1        
# [4] Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-4         
# [3] gridExtra_2.3          rlang_1.2.0           
# [5] magrittr_2.0.4         RcppAnnoy_0.0.23      
# [7] otel_0.2.0             matrixStats_1.5.0     
# [9] ggridges_0.5.7         compiler_4.5.2        
# [11] spatstat.geom_3.7-0    png_0.1-8             
# [13] vctrs_0.7.1            reshape2_1.4.5        
# [15] stringr_1.6.0          pkgconfig_2.0.3       
# [17] fastmap_1.2.0          labeling_0.4.3        
# [19] promises_1.5.0         ggbeeswarm_0.7.3      
# [21] purrr_1.2.1            jsonlite_2.0.0        
# [23] goftest_1.2-3          later_1.4.8           
# [25] spatstat.utils_3.2-1   irlba_2.3.7           
# [27] parallel_4.5.2         cluster_2.1.8.2       
# [29] R6_2.6.1               ica_1.0-3             
# [31] stringi_1.8.7          RColorBrewer_1.1-3    
# [33] spatstat.data_3.1-9    reticulate_1.45.0     
# [35] parallelly_1.46.1      spatstat.univar_3.1-6 
# [37] lmtest_0.9-40          scattermore_1.2       
# [39] tensor_1.5.1           future.apply_1.20.2   
# [41] zoo_1.8-15             sctransform_0.4.3     
# [43] httpuv_1.6.16          Matrix_1.7-4          
# [45] splines_4.5.2          igraph_2.2.2          
# [47] tidyselect_1.2.1       rstudioapi_0.19.0     
# [49] abind_1.4-8            spatstat.random_3.4-4 
# [51] codetools_0.2-20       miniUI_0.1.2          
# [53] spatstat.explore_3.7-0 listenv_0.10.1        
# [55] lattice_0.22-9         tibble_3.3.1          
# [57] plyr_1.8.9             withr_3.0.3           
# [59] shiny_1.13.0           S7_0.2.1              
# [61] ROCR_1.0-12            ggrastr_1.0.2         
# [63] Rtsne_0.17             fastDummies_1.7.6     
# [65] survival_3.8-6         polyclip_1.10-7       
# [67] fitdistrplus_1.2-6     pillar_1.11.1         
# [69] KernSmooth_2.23-26     plotly_4.12.0         
# [71] generics_0.1.4         RcppHNSW_0.6.0        
# [73] ggplot2_4.0.2          scales_1.4.0          
# [75] globals_0.19.1         xtable_1.8-8          
# [77] RhpcBLASctl_0.23-42    glue_1.8.0            
# [79] lazyeval_0.2.2         tools_4.5.2           
# [81] data.table_1.18.2.1    RSpectra_0.16-2       
# [83] RANN_2.6.2             dotCall64_1.2         
# [85] cowplot_1.2.0          grid_4.5.2            
# [87] tidyr_1.3.2            nlme_3.1-168          
# [89] patchwork_1.3.2        beeswarm_0.4.0        
# [91] vipor_0.4.7            cli_3.6.6             
# [93] spatstat.sparse_3.1-0  spam_2.11-3           
# [95] viridisLite_0.4.3      dplyr_1.2.0           
# [97] uwot_0.2.4             gtable_0.3.6          
# [99] digest_0.6.39          progressr_0.19.0      
# [101] ggrepel_0.9.7          htmlwidgets_1.6.4     
# [103] farver_2.1.2           htmltools_0.5.9       
# [105] lifecycle_1.0.5        httr_1.4.8            
# [107] mime_0.13              MASS_7.3-65
