# Watabe et al. F1000Research ver5
# This file is Watabe_F1000_ver5_new_No5.R for peak search
# Use 260319_retina_chromatin_d59_d74_d78_merged_No0.R for query scRNAseq import
# Use 260319_retina_chromatin_d59_d74_d78_merged_No1.R for query scRNAseq QC
# Use 260319_retina_chromatin_d59_d74_d78_merged_No2.R for scATACseq QC
# Use 260319_retina_chromatin_d59_d74_d78_merged_No3.R for scATACseq gene activity calculation
# Use Watabe_F1000_ver5_new_No1.R for scRNAseq processing
# Use Watabe_F1000_ver5_new_No2.R for scATACseq processing
# Use Watabe_F1000_ver5_new_No3.R for scATACseq processing

library(tidyr)
library(dplyr)

N1_PAX6 <- read.csv("FIMO_NOTCH1_PAX6_260725.csv")
N1_RAX <- read.csv("FIMO_NOTCH1_RAX_260725.csv")
N1_VSX2 <- read.csv("FIMO_NOTCH1_VSX2_260725.csv")
N1_peaks <- read.csv("NOTCH1_peak_zscore_RPC_260722.csv")

head(N1_PAX6)
# motif_id motif_alt_id sequence_name     start      stop strand   score
# 1 MA0069.1         PAX6          chr9 136611415 136611428      + 14.3000
# 2 MA0069.1         PAX6          chr9 136472750 136472763      + 13.7667
# 3 MA0069.1         PAX6          chr9 136447054 136447067      - 11.9444
# 4 MA0069.1         PAX6          chr9 136421968 136421981      - 11.9333
# 5 MA0069.1         PAX6          chr9 136494138 136494151      + 11.4000
# 6 MA0069.1         PAX6          chr9 136447319 136447332      + 11.3778
# p.value q.value matched_sequence
# 1 6.22e-06       1   TTCACACTTGAATT
# 2 9.13e-06       1   TTCACGCTTAAGCG
# 3 3.03e-05       1   TTTCTGCATGAGTG
# 4 3.04e-05       1   CTCACGGATGAGCG
# 5 4.22e-05       1   CTCACGCATGGTGT
# 6 4.28e-05       1   GTCCCGCTTGGGTT
N1_PAX6_start <- N1_PAX6[, 4]
N1_PAX6_end <- N1_PAX6[, 5]
N1_RAX_start <- N1_RAX[, 4]
N1_RAX_end <- N1_RAX[, 5]
N1_VSX2_start <- N1_VSX2[, 4]
N1_VSX2_end <- N1_VSX2[, 5]
TF_start <- c(N1_PAX6_start, N1_RAX_start, N1_VSX2_start)
TF_end <- c(N1_PAX6_end, N1_RAX_end, N1_VSX2_end)

N1_peaks_peak <- N1_peaks[, 8]
N1_peaks_peak <- as.data.frame(N1_peaks_peak)
N1_peaks_peak_clean <- N1_peaks_peak %>%
  separate(N1_peaks_peak, into = c("chr", "start", "end"), sep = "-")
# chr     start       end
# 1 chr9 136050144 136052719
# 2 chr9 136135336 136136100
# 3 chr9 136147805 136148796
# 4 chr9 136266989 136269273
# 5 chr9 136308974 136309488
# 6 chr9 136325871 136329045
peak_start <- N1_peaks_peak_clean[,2]
peak_start <- as.numeric(peak_start)
peak_end <- N1_peaks_peak_clean[,3]
peak_end <- as.numeric(peak_end)

included_matrix <- outer(TF_start, peak_start, ">=") & outer(TF_end, peak_end, "<=")
print(included_matrix)
write.csv(included_matrix, "included_matrix_260725.csv")
# TRUE: [15, 18] only RAX TFBS
print(TF_start[15]) # 136542636
print(TF_end[15]) # 136542645
print(peak_start[18]) # 136542126
print(peak_end[18]) # 136542839
length(N1_PAX6) # 10
length(N1_RAX) # 10
length(N1_VSX2) # 10
length(TF_start) # 28

sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin24.4.0
# Running under: macOS Tahoe 26.5.2
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
#   [1] dplyr_1.2.1 tidyr_1.3.2
# 
# loaded via a namespace (and not attached):
#   [1] R6_2.6.1          tidyselect_1.2.1  magrittr_2.0.5    glue_1.8.1       
# [5] tibble_3.3.1      pkgconfig_2.0.3   generics_0.1.4    lifecycle_1.0.5  
# [9] cli_3.6.6         vctrs_0.7.3       withr_3.0.3       compiler_4.5.1   
# [13] purrr_1.2.2       rstudioapi_0.19.0 tools_4.5.1       pillar_1.11.1    
# [17] rlang_1.2.0