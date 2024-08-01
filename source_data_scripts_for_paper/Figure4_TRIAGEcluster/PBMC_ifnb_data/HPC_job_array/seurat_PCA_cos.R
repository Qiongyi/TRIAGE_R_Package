set.seed(37)
library(Seurat)
library(ggplot2)
library(dplyr)
library(lsa)

args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[1]) 

immune.combined.sct <- readRDS("PBMC_ifnb_data/PBMC_Seurat_CCAint_data_multires_seuob.RDS")
seures <- colnames(immune.combined.sct@meta.data)[grepl("integrated_snn_res", colnames(immune.combined.sct@meta.data))]

clustnum <- as.character(unique(immune.combined.sct[[]][[seures[n]]])[order(unique(immune.combined.sct[[]][[seures[n]]]))])
dfls <- list()
for (i in clustnum){
	immune.combined.sct_1 <- immune.combined.sct[, which(immune.combined.sct[[]][, seures[n]]==i)]
    immune.combined.sct_1_df <- t(Embeddings(immune.combined.sct_1, "pca"))
    corr <- as.data.frame(lsa::cosine(immune.combined.sct_1_df)) 
    dfcorr <- data.frame(correlation = mean(corr[lower.tri(corr)]), cluster_peak = i, method = "Seurat") 
    dfls[[paste0("Seurat_cluster_", i)]] <- dfcorr
    dfls.all <- do.call(rbind, dfls)
    dfls.all$clust_num <- length(clustnum)
    dfls.all$resolution <- seures[n]
}

saveRDS(dfls.all, paste0("PBMC_ifnb_data/Seurat_PCA_cosine/res_", seures[n], "_PCA_withincluster_mean_cosine.RDS"))
