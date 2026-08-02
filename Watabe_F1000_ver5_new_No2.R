# Watabe et al. F1000Research ver5
# This file is Watabe_F1000_ver5_new_No2.R for scATACseq processing
# Use Watabe_F1000_ver5_new_No1.R for scRNAseq processing

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Signac)
library(GenomicRanges)
library(scales)
set.seed(123)

ATAC_merged <- readRDS("260301_retina_chromatin_ATAC_merged_No3.rds")
RNA_merged <- readRDS("Thomas_et_al_scRNAseq_query.rds")
head(RNA_merged@active.ident)
# d59_AAACCCACAAAGCTAA-1 d59_AAACCCACACTATCGA-1 
# Retinal ganglion cells         Amacrine cells 
# d59_AAACCCAGTACCCAGC-1 d59_AAACCCAGTCTACAGT-1 
# Amacrine cells Retinal ganglion cells 
# d59_AAACCCAGTCTTGGTA-1 d59_AAACCCAGTGACTAAA-1 
# Amacrine cells             Transitory 
# 12 Levels: Retinal ganglion cells ... Microglia
cell_counts <- table(RNA_merged@active.ident)
print(cell_counts)
# Retinal ganglion cells         Amacrine cells 
# 7933                   8870 
# Transitory    Retinal progenitors 
# 3074                   2632 
# Cones       Horizontal cells 
# 421                   1115 
# Rods        Cone precursors 
# 232                   1216 
# Astrocytes            Muller glia 
# 44                    534 
# Bipolar cells              Microglia 
# 297                     20
RNA_merged <- UpdateSeuratObject(RNA_merged)
transfer.anchors <- FindTransferAnchors(
  reference = RNA_merged,
  query = ATAC_merged,
  reduction = 'cca'
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = RNA_merged@active.ident,
  weight.reduction = ATAC_merged[['lsi']],
  dims = 2:30
)
ATAC_merged <- AddMetaData(object = ATAC_merged, metadata = predicted.labels)

# For remaining R codes, see 260319_retina_chromatin_d59_d74_d78_merged_No4.R
# celltypes <- levels(RNA_merged)
# color_vector <- hue_pal()(length(celltypes))
# rna_colors <- setNames(color_vector, celltypes)
# show_col(rna_colors)
# ATAC_merged$celltype <- factor(ATAC_merged$predicted.id, levels = names(rna_colors))

saveRDS(object = ATAC_merged, file = "260715_ATAC_annotated.rds")
celltypes <- levels(RNA_merged)
color_vector <- hue_pal()(length(celltypes))
rna_colors <- setNames(color_vector, celltypes)
show_col(rna_colors)
ATAC_merged$celltype <- factor(ATAC_merged$predicted.id, levels = names(rna_colors))
DimPlot(ATAC_merged,group.by = 'celltype')
# save image as DimPlot_ATAC_annotated_260715.tiff W800 H600
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
#   [1] future_1.70.0        scales_1.4.0        
# [3] GenomicRanges_1.62.1 Seqinfo_1.0.0       
# [5] IRanges_2.44.0       S4Vectors_0.48.1    
# [7] BiocGenerics_0.56.0  generics_0.1.4      
# [9] Signac_1.17.1        ggplot2_4.0.2       
# [11] patchwork_1.3.2      Seurat_5.4.0        
# [13] SeuratObject_5.3.0   sp_2.2-1            
# [15] dplyr_1.2.0         
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3       rstudioapi_0.19.0       
# [3] jsonlite_2.0.0           magrittr_2.0.4          
# [5] spatstat.utils_3.2-1     farver_2.1.2            
# [7] vctrs_0.7.1              ROCR_1.0-12             
# [9] spatstat.explore_3.7-0   Rsamtools_2.26.0        
# [11] RcppRoll_0.3.2           htmltools_0.5.9         
# [13] sctransform_0.4.3        parallelly_1.46.1       
# [15] KernSmooth_2.23-26       htmlwidgets_1.6.4       
# [17] ica_1.0-3                plyr_1.8.9              
# [19] plotly_4.12.0            zoo_1.8-15              
# [21] igraph_2.2.2             mime_0.13               
# [23] lifecycle_1.0.5          pkgconfig_2.0.3         
# [25] Matrix_1.7-4             R6_2.6.1                
# [27] fastmap_1.2.0            MatrixGenerics_1.22.0   
# [29] fitdistrplus_1.2-6       shiny_1.13.0            
# [31] digest_0.6.39            tensor_1.5.1            
# [33] RSpectra_0.16-2          irlba_2.3.7             
# [35] labeling_0.4.3           progressr_0.19.0        
# [37] spatstat.sparse_3.1-0    httr_1.4.8              
# [39] polyclip_1.10-7          abind_1.4-8             
# [41] compiler_4.5.2           withr_3.0.3             
# [43] S7_0.2.1                 BiocParallel_1.44.0     
# [45] fastDummies_1.7.6        MASS_7.3-65             
# [47] tools_4.5.2              lmtest_0.9-40           
# [49] otel_0.2.0               httpuv_1.6.16           
# [51] future.apply_1.20.2      goftest_1.2-3           
# [53] glue_1.8.0               nlme_3.1-168            
# [55] promises_1.5.0           grid_4.5.2              
# [57] Rtsne_0.17               cluster_2.1.8.2         
# [59] reshape2_1.4.5           gtable_0.3.6            
# [61] spatstat.data_3.1-9      tidyr_1.3.2             
# [63] data.table_1.18.2.1      XVector_0.50.0          
# [65] spatstat.geom_3.7-0      RcppAnnoy_0.0.23        
# [67] ggrepel_0.9.7            RANN_2.6.2              
# [69] pillar_1.11.1            stringr_1.6.0           
# [71] spam_2.11-3              RcppHNSW_0.6.0          
# [73] later_1.4.8              splines_4.5.2           
# [75] lattice_0.22-9           survival_3.8-6          
# [77] deldir_2.0-4             tidyselect_1.2.1        
# [79] Biostrings_2.78.0        miniUI_0.1.2            
# [81] pbapply_1.7-4            gridExtra_2.3           
# [83] scattermore_1.2          matrixStats_1.5.0       
# [85] stringi_1.8.7            UCSC.utils_1.6.1        
# [87] lazyeval_0.2.2           codetools_0.2-20        
# [89] tibble_3.3.1             cli_3.6.6               
# [91] uwot_0.2.4               xtable_1.8-8            
# [93] reticulate_1.45.0        dichromat_2.0-0.1       
# [95] Rcpp_1.1.1               GenomeInfoDb_1.46.2     
# [97] globals_0.19.1           spatstat.random_3.4-4   
# [99] png_0.1-8                spatstat.univar_3.1-6   
# [101] parallel_4.5.2           dotCall64_1.2           
# [103] sparseMatrixStats_1.22.0 bitops_1.0-9            
# [105] listenv_0.10.1           viridisLite_0.4.3       
# [107] ggridges_0.5.7           purrr_1.2.1             
# [109] crayon_1.5.3             rlang_1.2.0             
# [111] cowplot_1.2.0            fastmatch_1.1-8

# Memo
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

# This file (Watabe_F1000_ver5_new_No2.R) is for the scATACseq analysis.

# pre-processing

# library(Rsamtools)
# raw_file_d59 <- "C:/Users/ana_emb/Documents/notch-group/watabe/download_20260706/tsv/GSM5567518_d59_fragments.tsv"
# bgzip_file_d59 <- bgzip(raw_file_d59, dest = paste0(raw_file_d59, ".gz"))
# indexTabix(file = bgzip_file_d59, format = "bed")
# raw_file_d74 <- "C:/Users/ana_emb/Documents/notch-group/watabe/download_20260706/tsv/GSM5567519_d74_fragments.tsv"
# bgzip_file_d74 <- bgzip(raw_file_d74, dest = paste0(raw_file_d74, ".gz"))
# indexTabix(file = bgzip_file_d74, format = "bed")
# raw_file_d78 <- "C:/Users/ana_emb/Documents/notch-group/watabe/download_20260706/tsv/GSM5567520_d78_fragments.tsv"
# bgzip_file_d78 <- bgzip(raw_file_d78, dest = paste0(raw_file_d78, ".gz"))
# indexTabix(file = bgzip_file_d78, format = "bed")
# Tabiz files in "C:\\Users\\ana_emb\\Documents\\Notch-Group\\Watabe\\Download_20260706\\tsv"

# library(JASPAR2024)
# library(TFBSTools)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(ggplot2)
# library(RSQLite)
# library(ggseqlogo)
# library(motifmatchr)
# library(EnsDb.Hsapiens.v86)
# library(RSQLite)
# jaspar_db <- JASPAR2024()
# db_conn <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar_db))
# opts <- list(
#   species = "Homo sapiens",
#   collection = "CORE",
#   all_versions = FALSE
# )
# pfm <- getMatrixSet(x = db_conn, opts = opts)
# pwm <- toPWM(pfm)
# RSQLite::dbDisconnect(db_conn)
# saveRDS(pwm, "pwm_for_Footprint_20260708.rds")

# frag_path_d59 <- "C:/Users/ana_emb/Documents/notch-group/watabe/download_20260706/tsv/GSM5567518_d59_fragments.tsv.gz"
# frag_path_d74 <- "C:/Users/ana_emb/Documents/notch-group/watabe/download_20260706/tsv/GSM5567519_d74_fragments.tsv.gz"
# frag_path_d78 <- "C:/Users/ana_emb/Documents/notch-group/watabe/download_20260706/tsv/GSM5567520_d78_fragments.tsv.gz"

# From 260320_retina_chromatin_d59_before_merge.R
# cutoff <- 1000
# macs2_path <- "/home/notch/anaconda3/bin/macs2"
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
# 
# ATAC_total_counts_d59 <- CountFragments(frag_path_d59)
# ATAC_barcodes_d59 <- ATAC_total_counts_d59[ATAC_total_counts_d59$frequency_count > cutoff, ]$CB
# frags_d59 <- CreateFragmentObject(path = frag_path_d59, cells = ATAC_barcodes_d59)
# peaks_d59 <- CallPeaks(frags_d59, macs2.path = macs2_path)
# ATAC_counts_d59 <- FeatureMatrix(fragments = frags_d59,
#                                  features = peaks_d59,
#                                  cells = ATAC_barcodes_d59)
# chrom_assay_d59 <- CreateChromatinAssay(counts = ATAC_counts_d59,
#                                         sep = c(":", "-"),
#                                         fragments = frag_path_d59,
#                                         annotation = annotation)
# ATAC_SeuratObject_d59 <- CreateSeuratObject(counts = chrom_assay_d59,
#                                             assay = "peaks")
# saveRDS(ATAC_SeuratObject_d59, file = "241108_retina_chromatin_d59_No1.rds")

# d59 -> d74 "241108_retina_chromatin_d74_No1.rds"
# d59 -> d78 "241108_retina_chromatin_d78_No1.rds"

# From 260319_retina_chromatin_d59_d74_d78_merged_No0.R
# ATAC_d59 <- readRDS("241108_retina_chromatin_d59_No1.rds")
# ATAC_d74 <- readRDS("241108_retina_chromatin_d74_No1.rds")
# ATAC_d78 <- readRDS("241108_retina_chromatin_d78_No1.rds")
# 
# ATAC_d59@assays$peaks@fragments <- list()
# ATAC_d74@assays$peaks@fragments <- list()
# ATAC_d78@assays$peaks@fragments <- list()
# 
# frag_d59 <- CreateFragmentObject(
#   path = "C:/Users/ana_emb/Documents/Watabe/GSM5567518_d59_fragments.tsv.gz",
#   cells = colnames(ATAC_d59)
# )
# frag_d74 <- CreateFragmentObject(
#   path = "C:/Users/ana_emb/Documents/Watabe/GSM5567519_d74_fragments.tsv.gz",
#   cells = colnames(ATAC_d74)
# )
# frag_d78 <- CreateFragmentObject(
#   path = "C:/Users/ana_emb/Documents/Watabe/GSM5567520_d78_fragments.tsv.gz",
#   cells = colnames(ATAC_d78)
# )
# 
# ATAC_d59@assays$peaks@fragments <- list(frag_d59)
# ATAC_d74@assays$peaks@fragments <- list(frag_d74)
# ATAC_d78@assays$peaks@fragments <- list(frag_d78)
# 
# ATAC_merged <- merge(
#   x = ATAC_d59,
#   y = list(ATAC_d74, ATAC_d78),
#   add.cell.ids = c("d59", "d74", "d78")
# )

# saveRDS(ATAC_merged, "260228_retina_chromatin_ATAC_merged.rds")

# From 260319_retina_chromatin_d59_d74_d78_merged_No2.R
# ATAC_merged <- readRDS("260228_retina_chromatin_ATAC_merged.rds")
# 
# str(ATAC_merged)
# head(ATAC_merged)
# 
# ATAC_merged[['peaks']]
# ChromatinAssay data with 210972 features for 27440 cells
# Variable features: 0 
# Genome: 
#   Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 3 

# grange <- granges(ATAC_merged)
# 
# peaks.keep <- seqnames(granges(ATAC_merged)) %in% standardChromosomes(granges(ATAC_merged))
# ATAC_merged <- ATAC_merged[as.vector(peaks.keep), ]

# compute nucleosome signal score per cell
# ATAC_merged <- NucleosomeSignal(object = ATAC_merged)

# compute TSS enrichment score per cell
# ATAC_merged <- TSSEnrichment(object = ATAC_merged)

# not run
# peak_ranges should be a set of genomic ranges spanning the set of peaks to be quantified per cell
# peak_matrix <- FeatureMatrix(
#   fragments = Fragments(ATAC_merged),
#   features = grange
# )
# 
# # saveRDS(ATAC_merged, "260301_retina_chromatin_merged_No2_sub1.rds")
# # saveRDS(peak_matrix, "260301_retina_chromatin_merged_No2_sub2.rds")
# ATAC_merged <- readRDS("260301_retina_chromatin_merged_No2_sub1.rds")
# 
# # # not run
# # total_fragments <- CountFragments('GSM5567518_d59_fragments.tsv.gz')
# # rownames(total_fragments) <- total_fragments$CB
# # ATAC_merged$fragments <- total_fragments[colnames(ATAC_merged), "frequency_count"]
# 
# Fragments(ATAC_merged)
# # [[1]]
# # A Fragment object for 8071 cells
# # 
# # [[2]]
# # A Fragment object for 7095 cells
# # 
# # [[3]]
# # A Fragment object for 12274 cells
# 
# frag_list <- Fragments(ATAC_merged)
# 
# sample_ids <- c("d59", "d74", "d78")
# 
# frag_counts_list <- lapply(seq_along(frag_list), function(i) {
#   
#   counts <- CountFragments(frag_list[[i]]@path)
#   
#   # Add correct biological prefix
#   counts$CB <- paste0(sample_ids[i], "_", counts$CB)
#   
#   return(counts)
# })
# 
# total_fragments <- do.call(rbind, frag_counts_list)
# 
# rownames(total_fragments) <- total_fragments$CB
# 
# ATAC_merged$fragments <- 
#   total_fragments[colnames(ATAC_merged), "frequency_count"]
# 
# sum(is.na(ATAC_merged$fragments))
# # [1] 0
# summary(ATAC_merged$fragments)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# # 1001    6335   12298   15726   21206  311659 
# 
# ATAC_merged <- FRiP(
#   object = ATAC_merged,
#   assay = "peaks",
#   total.fragments = "fragments"
# )
# 
# # saveRDS(ATAC_merged, "260301_retina_chromatin_merged_No2_sub3.rds")
# 
# # add blacklist ratio
# ATAC_merged$blacklist_ratio <- FractionCountsInRegion(
#   object = ATAC_merged, 
#   assay = 'peaks',
#   regions = blacklist_hg38_unified
# )
# 
# DensityScatter(ATAC_merged, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
# 
# ATAC_merged$nucleosome_group <- ifelse(ATAC_merged$nucleosome_signal > 1.5, 'NS > 1.5', 'NS < 1.5')
# FragmentHistogram(object = ATAC_merged, group.by = 'nucleosome_group')
# # Warning messages:
# #   1: Removed 575 rows containing non-finite outside the scale range
# # (`stat_bin()`). 
# # 2: Removed 4 rows containing missing values or values outside the scale range
# # (`geom_bar()`). 
# 
# VlnPlot(
#   object = ATAC_merged,
#   features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'FRiP'),
#   pt.size = 0.1,
#   ncol = 5
# )
# 
# ATAC_merged_QC <- subset(
#   x = ATAC_merged,
#   subset = nCount_peaks > 2000 &
#     nCount_peaks < 30000 &
#     FRiP > 0.2 &
#     blacklist_ratio < 0.01 &
#     nucleosome_signal < 1.5 &
#     TSS.enrichment > 2.5
# )
# ATAC_merged_QC
# # An object of class Seurat 
# # 210972 features across 22231 samples within 1 assay 
# # Active assay: peaks (210972 features, 0 variable features)
# # 2 layers present: counts, data
# 
# # saveRDS(object = ATAC_merged_QC, file = "260301_retina_chromatin_merged_No2.rds")

# From 260319_retina_chromatin_d59_d74_d78_merged_No3.R
# ATAC_merged <- readRDS("260301_retina_chromatin_merged_No2.rds")
# 
# set.seed(1)
# 
# ATAC_merged <- RunTFIDF(ATAC_merged)
# ATAC_merged <- FindTopFeatures(ATAC_merged, min.cutoff = 'q0')
# ATAC_merged <- RunSVD(ATAC_merged)
# 
# DepthCor(ATAC_merged)
# 
# ATAC_merged <- RunUMAP(object = ATAC_merged, reduction = 'lsi', dims = 2:30)
# ATAC_merged <- FindNeighbors(object = ATAC_merged, reduction = 'lsi', dims = 2:30)
# ATAC_merged <- FindClusters(object = ATAC_merged, verbose = FALSE, algorithm = 3)
# DimPlot(object = ATAC_merged, label = TRUE) + NoLegend()
# 
# gene.activities <- GeneActivity(ATAC_merged)
# 
# # add the gene activity matrix to the Seurat object as a new assay and normalize it
# ATAC_merged[['RNA']] <- CreateAssayObject(counts = gene.activities)
# ATAC_merged <- NormalizeData(
#   object = ATAC_merged,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(ATAC_merged$nCount_RNA)
# )
# 
# DefaultAssay(ATAC_merged) <- 'RNA'
# 
# saveRDS(ATAC_merged, "260301_retina_chromatin_ATAC_merged_No3.rds")
