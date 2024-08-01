#############################################################
######         Figure 2. TRIAGE-prioritized genes      ######
#############################################################

### Extract go slim annotaton and gene description from BIomaRt

# if (!requireNamespace("biomaRt", quietly = TRUE)) {
#   install.packages("biomaRt")
# }
library(biomaRt)

# Use Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Set the working directory and make sure the file "Priority_epimap_rts.csv" is in the working directory
#setwd("C:/Users/uqqzhao/UQ/1projects/TRIAGE_R_package/Figure2_TRIAGE_prioritized_genes") 
file = read.csv("Priority_epimap_rts.csv", header=T, quote="")
genes <- file[which(file$Priority == "Y"),]

attributes_description <- c('hgnc_symbol', 'description')
annotations_description <- getBM(attributes = attributes_description, filters = 'hgnc_symbol', values = genes$Gene, mart = ensembl)
write.table(annotations_description, "anno_description.xls", sep="\t", row.names=FALSE, quote=FALSE)

attributes_go_slim <- c('hgnc_symbol', 'goslim_goa_accession', 'goslim_goa_description')
annotations_go_slim <- getBM(attributes = attributes_go_slim, filters = 'hgnc_symbol', values = genes$Gene, mart = ensembl)
write.table(annotations_go_slim, "anno_go_slim.xls", sep="\t", row.names=FALSE, quote=FALSE)


################################################################################
### Figure 2c. - functional category
library(ggplot2)
library(dplyr)

data <- read.table("Priority_genes_category.xls", header = FALSE, sep = "\t")
data$V1 <- factor(data$V1, levels = rev(data$V1))

num_categories <- nrow(data)
colors <- rainbow(num_categories)
colors[2] <- "gray"

#lollipop plot
ggplot(data, aes(x = V2, y = V1)) +
  geom_point(aes(color = factor(V1), size = V2)) +
  geom_segment(aes(x = 10, xend = V2, y = V1, yend = V1, color = factor(V1)), linewidth = 1) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5)) +
  xlim(10, 1000) +
  labs(x = "Number of genes", y = "Category") +
  geom_text(aes(label = V2), hjust = -0.9) +
  scale_size_continuous(range = c(6, 12))

ggsave("function_category_lollipop.pdf", width = 4, height = 4)

################################################################################
### Figure 2d. GO enrichment analysis
file = read.csv("Priority_epimap_rts.csv", header=T, quote="")
genes <- file[which(file$Priority == "Y"),]

library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library("org.Hs.eg.db")
organism="org.Hs.eg.db"

go <- enrichGO(gene = genes$Gene,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "MF",
               pAdjustMethod = "fdr",
               pvalueCutoff = 0.05)

write.table(go,"GO_MF_priority_genes.xls", row.names=FALSE, quote = FALSE, sep = "\t")

# Simplify the results to reduce redundancy
go_simplified <- simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)
write.table(go_simplified, "GO_MF_priority_genes_simplified.xls", row.names = FALSE, quote = FALSE, sep = "\t")

pdf("GO_MF_priority_genes.top20.dotplot.pdf", width = 5, height = 6)
dotplot(go_simplified, color = "p.adjust", showCategory = 20, orderBy = "x", font.size = 8) +
  scale_color_continuous(labels = scales::scientific_format(digits = 1)) +
  theme(axis.text.y = element_text(size = 6, hjust = 1),  # Increase the size of Y axis labels and adjust their position to the right
        plot.margin = margin(0, 0, 0, 0.5, "cm"))  # Increase the right margin to provide more space for Y axis labels
dev.off()




########################################################################################
###### Supplementary Figure 1. SuperPath and disease category enrichment analyses ######
########################################################################################
### for GeneCards  Superpath
library(readxl)
library(ggplot2)
library(dplyr)

data <- read_excel("GeneCards_gene_analysis.xlsx", sheet = "Pathways")

data$Name <- gsub("SuperPath: ", "", data$Name)

top_data <- data %>%
  arrange(desc(Score)) %>%
  top_n(20, Score)

top_data$Name <- factor(top_data$Name, levels = top_data$Name[order(top_data$Score)])

#head(top_data)

colnames(top_data)=c("Score", "Name", "Matched_info")
top_data$Score <- as.numeric(top_data$Score)
top_data$Matched_Genes <- as.numeric(sub(" \\(.*", "", top_data$Matched_info))

colors <- rainbow(20)
colors[length(colors) - 1] <- "gray"
colors[length(colors)] <- "lightgray"
colors <- rev(colors)

ggplot(top_data, aes(x = Score, y = Name)) +
  geom_point(aes(color = factor(Score), size = Matched_Genes)) +
  geom_segment(aes(x = 15, xend = Score, y = Name, yend = Name, color = factor(Score)), linewidth = 1) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5)) +
  xlim(15, 85) +
  labs(x = "Score", y = "SuperPath") +
  geom_text(aes(label = Matched_Genes), hjust = -0.9) +
  scale_size_continuous(range = c(6, 12))

ggsave("Genecards_pathway_lollipop.pdf", width = 8, height = 8)

################################################################################
### for GeneCards Diseases
library(readxl)
library(ggplot2)

data <- read_excel("GeneCards_gene_analysis.xlsx", sheet = "Diseases")

#head(data)

top_data <- data[data$Score > 40, ]
top_data <- top_data[order(-top_data$Score), ]
top_data$Name <- factor(top_data$Name, levels = top_data$Name[order(top_data$Score)])

#head(top_data)
#top_data[1,3]
colnames(top_data)=c("Score", "Name", "Matched_info")

top_data$Score <- as.numeric(top_data$Score)
top_data$Matched_Genes <- as.numeric(sub(" \\(.*", "", top_data$Matched_info))

colors <- rainbow(13)
colors[length(colors) - 1] <- "gray"
colors[length(colors)] <- "lightgray"
colors <- rev(colors)

ggplot(top_data, aes(x = Score, y = Name)) +
  geom_point(aes(color = factor(Score), size = Matched_Genes)) +
  geom_segment(aes(x = 35, xend = Score, y = Name, yend = Name, color = factor(Score)), linewidth = 1) +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.position = "none", 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5)) +
  xlim(35, 60) +
  labs(x = "Score", y = "Disease") +
  #geom_text(aes(label = paste(Matched_Genes, "genes")), hjust = -0.4) +
  geom_text(aes(label = Matched_Genes), hjust = -1.2) +
  scale_size_continuous(range = c(6, 12))

ggsave("Genecards_disease_lollipop.pdf", width = 6, height = 6)



#############################################################
######        Figure 3 & Table 1. TRIAGEgene           ######
######       In vivo mouse RNA-seq data analysis       ######
#############################################################

### In vivo mouse RNA-seq dataset is from this study: Quaife-Ryan, G. A. et al. Multicellular Transcriptional Analysis of Mammalian Heart Regeneration. Circulation 136, 1123-1139 (2017). 
### https://pubmed.ncbi.nlm.nih.gov/28733351/

# Download data from GEO 95755:  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse95755

#setwd("C:/Users/uqqzhao/UQ/1projects/TRIAGE_R_package/Figure3_TRIAGEgene/GSE95755") 
file = read.table("GSE95755_MultiCellularRNAseq_EdgeR_CPM.txt", header=T, quote="", sep="\t")

#head(file)
# P56 Myo
cpm <- file[,c("Symbol", "ShamP56_Myo_1", "ShamP56_Myo_2", "ShamP56_Myo_3", "ShamP56_Myo_4", "MIP56_Myo_1", "MIP56_Myo_2", "MIP56_Myo_3", "MIP56_Myo_4")]
rownames(cpm) <- file$Symbol
cpm <- cpm[,-1]
#head(cpm)

library(readxl)
deg <- read_excel("GSE95755_MultiCellularRNAseq_AllSignificantRegulatedGenes.xlsx", sheet = "MIP56vsShP56.Myo")
# dim(deg)
# [1] 658   9

cpm_deg <- cpm[deg$Symbol, ]
#head(cpm_deg)


################################################################################
### Figure 3a) Jaccard index heatmap showing sample similarities using the ‘plotJaccard’ function in the TRIAGE R package. 

########### Run TRIAGEgene #############
# See the installation guideline here: https://triage-r-package.readthedocs.io/en/latest/Installation.html
#install.packages("path/to/TRIAGE_1.1.3.tar.gz", repos = NULL, type = "source")
library(TRIAGE)
ds <- TRIAGEgene(cpm_deg, species = "Mouse")
#head(ds)
plotJaccard(ds, "Fig3a_Jaccard_heatmap.pdf")


################################################################################
### Top 100 genes with DS ranking
library(dplyr)
library(TRIAGE)

file = read.table("GSE95755_MultiCellularRNAseq_EdgeR_CPM.txt", header=T, quote="", sep="\t")
cpm <- file[,c("Symbol", "ShamP56_Myo_1", "ShamP56_Myo_2", "ShamP56_Myo_3", "ShamP56_Myo_4", "MIP56_Myo_1", "MIP56_Myo_2", "MIP56_Myo_3", "MIP56_Myo_4")]
rownames(cpm) <- file$Symbol
cpm <- cpm[,-1]
ds <- TRIAGEgene(cpm, species = "Mouse")

ds$MIP56_Myo_DS <- rowMeans(ds[, c("MIP56_Myo_1", "MIP56_Myo_2", "MIP56_Myo_3", "MIP56_Myo_4")])
top100_genes_ds <- ds %>%
  arrange(desc(MIP56_Myo_DS)) %>%
  dplyr::slice(1:100) %>%
  rownames()

write.table(top100_genes_ds, file = "top100_genes_ds.xls", row.names = FALSE, col.names = FALSE, quote = FALSE)

################################################################################
### Table 1. Top 10 genes ranked by TRIAGE-weighted values (discordance scores).
# top100_genes_ds[1:10]
#[1] "Gata4"   "Nkx2-5"  "Tbx5"    "Irx4"    "Ntn1"    "Scn4b"   "Rnf220"  "Metap1d" "Hand2"   "Tbx20"



# GO enrichment analysis for top 100 genes with DS ranking
library(ggplot2)
# packageVersion("ggplot2")
# [1] ‘3.4.0’

library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
organism="org.Mm.eg.db"

top100_genes <- top100_genes_ds
go <- enrichGO(gene = top100_genes,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "BP",
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               qvalueCutoff = 1)

write.table(go,"P56myo_DStop100_GOBP.xls", row.names=FALSE, quote = FALSE, sep = "\t")



################################################################################
### Top 100 genes selected based on adjusted p-value
library(readxl)
deg <- read_excel("GSE95755_MultiCellularRNAseq_AllSignificantRegulatedGenes.xlsx", sheet = "MIP56vsShP56.Myo")

# dim(deg)
# [1] 658   9

top100_genes_pvalue <- head(deg, 100)$Symbol
top100_genes <- top100_genes_pvalue

# GO enrichment analysis
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
organism="org.Mm.eg.db"

go <- enrichGO(gene = top100_genes,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "BP",
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               qvalueCutoff = 1)

write.table(go,"P56myo_top100_GOBP.xls", row.names=FALSE, quote = FALSE, sep = "\t")

################################################################################
### Figure 3b) Venn diagram showing the distribution of the top 100 DS genes, top 100 DE genes, and all DE genes.  
overlap_genes <- intersect(top100_genes_ds, deg$Symbol)
# > overlap_genes
#[1] "Inhbb" "Gdf6"  "Olfm1"


################################################################################
### Figure 3c) Comparisons of the top 20 enriched GO terms between the top 100 DS genes and the top 100 DE genes.
library(dplyr)
library(ggplot2)

go1 <- read.table("P56myo_DStop100_GOBP.xls", header=T, quote="", sep="\t")
go1_top <- go1[order(go1$p.adjust), ][1:20, "Description"]

go2 <- read.table("P56myo_top100_GOBP.xls", header=T, quote="", sep="\t")
go2_top <- go2[order(go2$p.adjust), ][1:20, "Description"]

combined_terms <- unique(c(go1_top, go2_top))

# length(combined_terms)
# [1] 40

go1_filtered <- go1[go1$Description %in% combined_terms, c("Description", "p.adjust", "GeneRatio")]
go1_filtered$Source <- "DS"

go2_filtered <- go2[go2$Description %in% combined_terms, c("Description", "p.adjust", "GeneRatio")]
go2_filtered$Source <- "Adjusted p-value"

merged_go <- rbind(go1_filtered, go2_filtered)

merged_go$GeneRatio <- sapply(merged_go$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})

merged_go <- as.data.frame(merged_go)
# dim(merged_go)
# [1] 71  4

filtered_go <- merged_go[merged_go$p.adjust <= 0.05, ]
# dim(filtered_go)
# [1] 44  4
min_size <- -log10(0.05)
max_size <- max(-log10(filtered_go$p.adjust))

# max_size
# [1] 19.50746

pdf("Fig3c_GO_terms_compare_top20.pdf",  width=8, height=10)
ggplot(filtered_go, aes(x = reorder(Description, p.adjust), y = GeneRatio, color = Source, size = -log10(p.adjust))) +
  geom_point(alpha = 0.7, stroke = 0, shape = 16) +
  scale_color_manual(values = c("DS" = "blue", "Adjusted p-value" = "orange")) +
  coord_flip() +
  theme_bw() +
  labs(title = "Top GO terms from two ranking methods",
       x = "",
       y = "Gene Ratio",
       color = "Source",
       size = "-log10(p.adjust)")+
  scale_size(range = c(1, 8), limits = c(min_size, max_size)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

dev.off()



#############################################################
######            Figure 4. TRIAGEcluster              ######
###### In vivo human single-cell RNA-seq data analysis ######
#############################################################

set.seed(37)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SeuratData)
library(patchwork)
library(cowplot)

# packageVersion("Seurat")
# [1] ‘4.3.0’

#setwd("C:/Users/uqqzhao/UQ/1projects/TRIAGE_R_package/Figure4_TRIAGEcluster/PBMC_ifnb_data")

################################################################################
### Seurat data integration
#load ifnb data
install.packages("http://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source")
library(ifnb.SeuratData)
data("ifnb")

# > ifnb
# An object of class Seurat 
# 14053 features across 13999 samples within 1 assay 
# Active assay: RNA (14053 features, 0 variable features)

ifnb.list <- SplitObject(ifnb, split.by = "stim")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)

# Multiple resolution clustering for the correlation analysis within clusters
res <- seq(0, 5, by = 0.1)
for (i in res){
  immune.combined <- FindClusters(immune.combined, resolution = i)
}

saveRDS(immune.combined, "PBMC_ifnb_data/PBMC_Seurat_CCAint_data_multires_seuob.RDS")

metadf <- function(seu){
  g.meta <- cbind(seu@meta.data, as.data.frame(Embeddings(seu[["umap"]])))
  g.meta$cell_name <- rownames(g.meta) 
  return(g.meta)
}

# Read data from the saved RDS file
immune.combined <- readRDS("PBMC_ifnb_data/PBMC_Seurat_CCAint_data_multires_seuob.RDS")

intdata <- as.data.frame(as.matrix(GetAssayData(immune.combined, assay = "integrated", slot = "data")))
intpca <- as.data.frame(Embeddings(immune.combined, "pca"))
all_genes <- as.data.frame(as.matrix(GetAssayData(immune.combined, assay = "RNA", slot = "data")))
intmeta <- metadf(immune.combined)
write.csv(intdata, "PBMC_ifnb_data/pbmc_ifnb_Seurat_CCAint_data_HVG.csv", quote=F)
write.csv(intpca, "PBMC_ifnb_data/pbmc_ifnb_Seurat_CCAint_pca.csv", quote=F)
write.csv(all_genes, "PBMC_ifnb_data/pbmc_ifnb_Seurat_CCAint_allgenes.csv", quote=F)
write.csv(intmeta, "PBMC_ifnb_data/pbmc_ifnb_metadata_withUMAP.csv", quote=F)

immune.combined <- FindClusters(immune.combined, resolution = 0.5)

#head(immune.combined)
#p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim", pt.size = 1.6, raster = T) + ggtitle("batch")
#p2 <- DimPlot(immune.combined, reduction = "umap", group.by = "seurat_annotations", pt.size = 1.6, raster = T) + ggtitle("annotation")
#plot_grid(p1, p2)


################################################################################
### Figure 4a
# Read data from the RDS file ("PBMC_Seurat_CCAint_data_multires_seuob.RDS") to obtain the exact coordinates of points on a UMAP plot. This is for reproducing Figure 4a.
# Note here: https://github.com/satijalab/seurat/issues/3254
# "The exact location of points on a UMAP can chance across different computers and OSs. 
# We do our best to minimize any randomness to the procedure by fixing the random seed, 
# but some fluctuation across systems is inevitable, and nothing to worry about."
################################################################################
seu <- readRDS("./PBMC_ifnb_data/PBMC_Seurat_CCAint_data_multires_seuob.RDS")
p1 <- DimPlot(seu, reduction = "umap", group.by = "stim", pt.size = 0.1) + ggtitle("Conditions")
p2 <- DimPlot(seu, reduction = "umap", group.by = "seurat_annotations", pt.size = 0.1) + ggtitle("Cell types")
ggsave("Figure4a1_umap.pdf", plot = p1, width = 8, height = 8)
ggsave("Figure4a2_umap.pdf", plot = p2, width = 8, height = 8)



################################################################################
### Run TRIAGEcluster
# See the installation guideline here: https://triage-r-package.readthedocs.io/en/latest/Installation.html
#install.packages("path/to/TRIAGE_1.1.3.tar.gz", repos = NULL, type = "source")
library(TRIAGE)

TRIAGEcluster(expr = "PBMC_ifnb_data/pbmc_ifnb_Seurat_CCAint_allgenes.csv",
              metadata = "PBMC_ifnb_data/pbmc_ifnb_metadata_withUMAP.csv",
              outdir= "PBMC_ifnb_data/TRIAGE_Cluster/",
              output_prefix = "TRIAGE_Cluster",
              cell_column = "cell_name",
              umap_column = "UMAP_",
              seed=37,
              bw = c(seq(0.1, 0.2, 0.01), seq(0.3, 5.0, 0.1)))


################################################################################
### Figure 4b
# Figure 4b is the UMAP generated in the above step with a bandwidth of 0.4 (i.e. bw=0.4)
# "PBMC_ifnb_data/TRIAGE_Cluster/TRIAGE_Cluster_bw0.40_labelledUMAP.pdf"
################################################################################


################################################################################
### Run TRIAGE 'byPeak' to extract peak-level average data
### Supplementary Data 1
# Load the TRIAGE R package
library(TRIAGE)

result <- byPeak(expr = "PBMC_ifnb_data/pbmc_ifnb_Seurat_CCAint_allgenes.csv",
                 peak = "PBMC_ifnb_data/TRIAGE_Cluster/TRIAGE_Cluster_bw0.40_metadata.csv",
                 peak_column="Peak",
                 cell_column="cell_name")
write.csv(result, "PBMC_ifnb_data/Peak_AvgExp.csv", quote=F)
# Supplementary Data 1 is "Peak_AvgExp.csv".


################################################################################
### Run TRIAGEgene based on peak-level average data
ds <- TRIAGEgene(result)
write.csv(ds, "PBMC_ifnb_data/Peak_DS.csv", quote=F)


################################################################################
### Optional analysis
exp <- read.csv("PBMC_ifnb_data/Peak_AvgExp.csv", row.names = 1)
# Extract top 10 genes for each peak
get_top_genes <- function(df, peaks, top_n = 10) {
  top_genes_list <- list()
  
  for (peak in peaks) {
    if (peak %in% colnames(df)) {
      # Sort genes by expression level for the current peak
      sorted_genes <- df[order(-df[[peak]]), ]
      top_genes <- head(sorted_genes, top_n)
      top_genes_list[[peak]] <- top_genes
    } else {
      warning(paste("Peak column", peak, "not found in data"))
    }
  }
  
  return(top_genes_list)
}

peaks_of_interest <- c("Peak0", "Peak1", "Peak2")

# Get top 10 genes for each peak
top_genes <- get_top_genes(exp, peaks_of_interest)



################################################################################
### Run TRIAGE 'byPeak' to extract cluster-level average data
### Supplementary Data 2
# read the metadata file: can be "TRIAGE_Cluster_bw0.40_metadata.csv" or any other metadata files with Seurat cluster information
result <- byPeak(expr = "PBMC_ifnb_data/pbmc_ifnb_Seurat_CCAint_allgenes.csv",
                 peak = "PBMC_ifnb_data/TRIAGE_Cluster/TRIAGE_Cluster_bw0.40_metadata.csv",
                 peak_column="integrated_snn_res.0.6",
                 cell_column="cell_name",
                 prefix = "Cluster")
write.csv(result, "PBMC_ifnb_data/Cluster_AvgExp.csv", quote=F)
# Supplementary Data 2 is "Cluster_AvgExp.csv".


################################################################################
### Run TRIAGEgene based on cluster-level average data
ds <- TRIAGEgene(result)
write.csv(ds, "PBMC_ifnb_data/Cluster_DS.csv", quote=F)



################################################################################
### Run correlation analysis
# Seurat within cluster spearman correlation and cosine similarity run on HPC as job array
# PBMC_ifnb_data/HPC_job_array/seurat_integrated_corr.R
# PBMC_ifnb_data/HPC_job_array/seurat_integrated_cos.R
# PBMC_ifnb_data/HPC_job_array/seurat_PCA_corr.R
# PBMC_ifnb_data/HPC_job_array/seurat_PCA_cos.R
# combine results as below, saved .csv file used in comparison with TRIAGE-Cluster using "peak_correlation_and_plot.py" 


################################################################################
# integrated high-dimensional data
files <- list.files("PBMC_ifnb_data/Seurat_integrated_correlation", ".RDS", full.names = T)
files1 = files[-1]
seuls <- lapply(files1, readRDS)
corall <- do.call(rbind, seuls)
colnames(corall)[colnames(corall) == 'correlation'] <- 'spearman_correlation'

files <- list.files("PBMC_ifnb_data/Seurat_integrated_cosine", ".RDS", full.names = T)
files1 = files[-1]
seuls <- lapply(files1, readRDS)
corall1 <- do.call(rbind, seuls)
colnames(corall1)[colnames(corall1) == 'correlation'] <- 'cosine_similarity'

# identical(rownames(corall), rownames(corall1)) #TRUE
corall1["spearman_correlation"]=corall["spearman_correlation"]
corall2=corall1[, c('cosine_similarity', 'spearman_correlation', 'cluster_peak', 'resolution', 'clust_num', 'method')]
write.csv(corall2, "PBMC_ifnb_data/pbmc_ifnb_PCA_Seurat_within_cluster_mean_correlation_similarity.csv", quote = F, row.names = T)

################################################################################
# integrated low-dimensional PCA data, same as above
files <- list.files("PBMC_ifnb_data/Seurat_PCA_correlation", ".RDS", full.names = T)
files1 = files[-1]
seuls <- lapply(files1, readRDS)
corall <- do.call(rbind, seuls)
colnames(corall)[colnames(corall) == 'correlation'] <- 'spearman_correlation'

files <- list.files("PBMC_ifnb_data/Seurat_PCA_cosine", ".RDS", full.names = T)
files1 = files[-1]
seuls <- lapply(files1, readRDS)
corall1 <- do.call(rbind, seuls)
colnames(corall1)[colnames(corall1) == 'correlation'] <- 'cosine_similarity'

#identical(rownames(corall), rownames(corall1)) #TRUE
corall1["spearman_correlation"]=corall["spearman_correlation"]
corall2=corall1[, c('cosine_similarity', 'spearman_correlation', 'cluster_peak', 'resolution', 'clust_num', 'method')]
write.csv(corall2, "PBMC_ifnb_data/pbmc_ifnb_PCA_Seurat_within_cluster_mean_correlation_similarity.csv", quote = F, row.names = T)


################################################################################
### Figure 4c
# Run the python script "peak_correlation_and_plot.py" to generate Figure 4c.




#############################################################
######            Figure 5. TRIAGEparser               ######
######      In vitro human RNA-seq data analysis       ######
#############################################################

################################################################################
### In vitro human RNA-seq dataset is from this study: Wu, Z. et al. Wnt dose escalation during the exit from pluripotency identifies tranilast as a regulator of cardiac mesoderm. Dev Cell 59, 705-722 e708 (2024). 
### https://pubmed.ncbi.nlm.nih.gov/38354738/


### select 1431 DE genes (>2 fold changes and padj<0.05) for the GO enrichment analysis
#setwd("C:/Users/uqqzhao/UQ/1projects/TRIAGE_R_package/Figure5_TRIAGEparser/GSE246079")
file = read.table("1T_vs_1_DESeq2_p0.05_231011.xls", header=T, quote="", sep="\t")

file <- file[abs(file$log2FoldChange) > 1 & file$padj < 0.05, ]
#dim(file)
#[1] 1431   22
genes <- file$Row.names
writeLines(genes, "DE_genes_1431_list.txt")

################################################################################
### GO enrichment analysis for all DE genes
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
organism="org.Hs.eg.db"
go <- enrichGO(gene = genes,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "ALL",
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               qvalueCutoff = 1)
write.table(go,"GO_enrichment_results_1T_vs_1_DEgenes1431.xls", row.names=FALSE, quote = FALSE, sep = "\t")


################################################################################
### Run TRIAGE R package: TRIAGEparser
# See the installation guideline here: https://triage-r-package.readthedocs.io/en/latest/Installation.html
#install.packages("path/to/TRIAGE_1.1.3.tar.gz", repos = NULL, type = "source")
library(TRIAGE)
TRIAGEparser("DE_genes_1431_list.txt", input_type = "list", outdir = "TRIAGEparser_output_DE1431")
plotGO(indir = "./TRIAGEparser_output_DE1431", outdir = "./TRIAGEparser_output_DE1431")

### Run TRIAGE R package: getClusterGenes
cluster1_genes <- getClusterGenes("./TRIAGEparser_output_DE1431/gene_clusters/output_gene_clusters.csv", "cluster1")
writeLines(cluster1_genes, "cluster1_genes.txt")
# length(cluster1_genes)
# [1] 262

cluster2_genes <- getClusterGenes("./TRIAGEparser_output_DE1431/gene_clusters/output_gene_clusters.csv", "cluster2")
writeLines(cluster2_genes, "cluster2_genes.txt")
length(cluster2_genes)
# [1] 276

cluster3_genes <- getClusterGenes("./TRIAGEparser_output_DE1431/gene_clusters/output_gene_clusters.csv", "cluster3")
writeLines(cluster3_genes, "cluster3_genes.txt")
length(cluster3_genes)
# [1] 557

cluster4_genes <- getClusterGenes("./TRIAGEparser_output_DE1431/gene_clusters/output_gene_clusters.csv", "cluster4")
writeLines(cluster4_genes, "cluster4_genes.txt")
length(cluster4_genes)
# [1] 7


################################################################################
### GO enrichment analysis for gene clusters
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
organism="org.Hs.eg.db"

### gene cluster 1
file = read.table("cluster1_genes.txt", header=F, quote="", sep="\t")
go <- enrichGO(gene = file$V1,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "ALL",
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               qvalueCutoff = 1)
write.table(go,"GO_enrichment_results_1T_vs_1_genecluster1.xls", row.names=FALSE, quote = FALSE, sep = "\t")

### gene cluster 2
file = read.table("cluster2_genes.txt", header=F, quote="", sep="\t")
go <- enrichGO(gene = file$V1,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "ALL",
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               qvalueCutoff = 1)
write.table(go,"GO_enrichment_results_1T_vs_1_genecluster2.xls", row.names=FALSE, quote = FALSE, sep = "\t")

### gene cluster 3
file = read.table("cluster3_genes.txt", header=F, quote="", sep="\t")
go <- enrichGO(gene = file$V1,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "ALL",
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               qvalueCutoff = 1)
write.table(go,"GO_enrichment_results_1T_vs_1_genecluster3.xls", row.names=FALSE, quote = FALSE, sep = "\t")

### skip gene cluster 4 as there are only 7 genes in cluster 4


################################################################################
### Figure 5a
### compare GO terms - keywords: c("Wnt", "embryonic organ", "embryonic heart", "mesoderm", "mesenchyme", "pattern specification")
library(dplyr)
key_words <- c("Wnt", "embryonic organ", "embryonic heart", "mesoderm", "mesenchyme", "pattern specification")

process_file <- function(file) {
  data <- read.delim(file, header = TRUE, sep = "\t")
  
  filtered_data <- data %>%
    filter(grepl(paste(key_words, collapse = "|"), Description, ignore.case = TRUE))
  return(filtered_data)
}

files <- list("GO_enrichment_results_1T_vs_1_DEgenes1431.xls", "GO_enrichment_results_1T_vs_1_genecluster1.xls", "GO_enrichment_results_1T_vs_1_genecluster2.xls", "GO_enrichment_results_1T_vs_1_genecluster3.xls")
processed_data <- lapply(files, process_file)

all_go_terms <- unique(unlist(lapply(processed_data, function(x) x$Description)))

result <- data.frame(Description = all_go_terms, stringsAsFactors = FALSE)

for (i in 1:length(files)) {
  file_data <- processed_data[[i]]
  p_adjust_col <- file_data %>% select(Description, p.adjust)
  colnames(p_adjust_col)[2] <- gsub(".xls", "", files[i])
  result <- left_join(result, p_adjust_col, by = c("Description" = "Description"))
}

#head(result)
result[is.na(result)] <- 1
dim(result)
# [1] 38  5
result <- result %>% filter_if(is.numeric, any_vars(. < 0.05))
dim(result)
# [1] 16  5
colnames(result)<-c("Description", "all_genes", "cluster_1", "cluster_2", "cluster_3")

# Sort by cluster_2 p.adjust
result <- result %>% arrange(cluster_2)

colnames(result)<-c("Description", "all genes", "cluster 1", "cluster 2", "cluster 3")
write.table(result, "GO_1T_vs_1_DEgenes1431_combined.keywords.xls", sep = "\t", row.names = FALSE, quote = FALSE)

rownames(result) <- result$Description
go <- result[,-1,drop = FALSE]
go[go > 0.05] <- NA
go <- -log10(go)

# Create custom color palette : adjusted p-value>0.05 will be marked as lightgrey
color_palette <- colorRampPalette(c("lightgrey", "lavender", "purple", "red"))(100)
breaks <- c(seq(0, 1.3, length.out = 50), seq(1.31, max(go, na.rm = TRUE), length.out = 51))

pdf("GO_1T_vs_1_DEgenes1431_combined.keywords.pdf", width = 6, height = 6)
pheatmap(go, cluster_cols = FALSE, cluster_rows = FALSE, color = color_palette, breaks = breaks, border_color = NA, fontsize_row = 8)
dev.off()

# Figure 5a is "GO_1T_vs_1_DEgenes1431_combined.keywords.pdf".



################################################################################
### Figure 5b - gene network
# using the cnetplot function from the R package ClusterProfiler
# plot gene networks: mesoderm development, Wnt signalling pathway, cell-cell signaling by wnt, embryonic heart tube morphogenesis
library(clusterProfiler)
library(org.Hs.eg.db)
organism="org.Hs.eg.db"
library(dplyr)

#setwd("C:/Users/uqqzhao/UQ/1projects/TRIAGE_R_package/Figure5_TRIAGEparser/GSE246079")
go_results <- read.delim("GO_enrichment_results_1T_vs_1_genecluster2.xls", header = TRUE, sep = "\t")
keywords <- c("mesoderm development", "Wnt signaling pathway", "embryonic heart tube morphogenesis")

filtered_go <- go_results[go_results$Description %in% keywords, ]
genes <- unique(unlist(strsplit(filtered_go$geneID, "/")))

de_genes <- read.delim("1T_vs_1_DESeq2_p0.05_231011.xls", header = TRUE, sep = "\t")
de_genes <- de_genes %>% filter(Row.names %in% genes)
fold_changes <- setNames(de_genes$log2FoldChange, de_genes$Row.names)

### gene cluster 2
file = read.table("cluster2_genes.txt", header=F, quote="", sep="\t")
go <- enrichGO(gene = file$V1,
               OrgDb = organism, 
               keyType = 'SYMBOL',
               ont = "ALL",
               pAdjustMethod = "fdr",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)

# cnetplot
library(ggplot2)
pdf("gene_cluster2_heart_mesoderm_network.pdf", width = 6, height = 8)
cnetplot(go, showCategory = keywords, circular = F, color.params = list(foldChange = fold_changes, edge = T)) + guides(edge_color = "none")
dev.off()

# Figure 5b is "gene_cluster2_heart_mesoderm_network.pdf"



#############################################################
######          Supplementary Figure 2                 ######
#############################################################
### Comparison of top enriched GO terms for the genes in each cluster as well as for all DE genes
library(dplyr)
process_file <- function(file) {
  data <- read.delim(file, header = TRUE, sep = "\t")
  # sort based on p.adjust and take top 10 GO terms
  top_10 <- data %>%
    arrange(p.adjust) %>%
    slice(1:10)
  
  return(top_10)
}

# Note: The below four files were generated in previous steps
files <- list("GO_enrichment_results_1T_vs_1_DEgenes1431.xls", "GO_enrichment_results_1T_vs_1_genecluster1.xls", "GO_enrichment_results_1T_vs_1_genecluster2.xls", "GO_enrichment_results_1T_vs_1_genecluster3.xls")
processed_data <- lapply(files, process_file)

all_go_terms <- unique(unlist(lapply(processed_data, function(x) x$Description)))
result <- data.frame(Description = all_go_terms, stringsAsFactors = FALSE)

for (i in 1:length(files)) {
  file_data <- processed_data[[i]]
  p_adjust_col <- file_data %>% select(Description, p.adjust)
  colnames(p_adjust_col)[2] <- gsub(".xls", "", files[i])
  result <- left_join(result, p_adjust_col, by = c("Description" = "Description"))
}

#head(result)
result[is.na(result)] <- 1
colnames(result)<-c("Description", "all genes", "cluster 1", "cluster 2", "cluster 3")
write.table(result, "GO_1T_vs_1_DEgenes1431_combined.xls", sep = "\t", row.names = FALSE, quote = FALSE)

rownames(result) <- result$Description
go <- result[,-1,drop = FALSE]
go <- -log10(go)
color_palette <- colorRampPalette(c("lightgrey", "purple", "red"))(n = 100)

library(pheatmap)
pdf("GO_1T_vs_1_DEgenes1431_combined.pdf", width = 6, height = 8)
pheatmap(go, cluster_cols = FALSE, cluster_rows = FALSE, color = color_palette, border_color = NA, fontsize_row = 8)
dev.off()

# Supplementary Figure 2 is "GO_1T_vs_1_DEgenes1431_combined.pdf".

