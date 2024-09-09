# To facilitate the reproducibility of our analyses and to ensure that users can fully utilize the TRIAGE R package in various scenarios, we have provided detailed analysis steps below.


#############################################################
######         Figure 1. TRIAGE-prioritized genes      ######
#############################################################

###### Source data are listed at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Figure1_TRIAGE_prioritized_genes/"


#################
### Figure 1a ### Define the TRIAGE-prioritized genes: elbow point detection based on RTS values in "human_rts_epimap.txt".
#################

## To define the TRIAGE-prioritized genes, follow the steps below to detect the elbow point based on RTS values in "human_rts_epimap.txt".

# 1. Ensure you place "human_rts_epimap.txt" in the same directory as "elbow_point_detection.py".

# 2. Install necessary Python modules:
# You need 'numpy', 'matplotlib', and 'pandas' to run this python script. You can install these modules using pip:
pip install numpy matplotlib pandas

# 3. Run the script:
python3 ./elbow_point_detection.py

######################
### Figure 1b & 1c ### Annotation and functional categorization of 993 TRIAGE-prioritized genes
######################
# 1) Download gene symbols and gene aliases for human genes from NCBI Gene 
#     - using keywords: (Homo sapiens genome) AND "Homo sapiens"[porgn:__txid9606] -> Send to File ("gene_result.txt")

# 2) Add gene symbols and aliases into TRIAGE-prioritized genes
anno_gene_alias.pl ./NCBI/gene_result.txt Priority_epimap_rts.csv Priority_epimap_rts_alias.xls

# 3) Download gene-disease associations from the DisGeNET platform ("gene_associations.tsv") and add disease information for TRIAGE-prioritized genes
anno_DisGeNET.pl ./DisGeNET/gene_associations.tsv Priority_epimap_rts_alias.xls Priority_epimap_rts_alias_DisGeNET.xls

# 4) Add coding and noncoding information for TRIAGE-prioritized genes using NCBI RefSeq files. We used two sets of NCBI RefSeq files to accommodate both the original version of RefSeq annotation and a recently updated version of RefSeq annotation.
anno_coding_noncoding_refseq.pl ./NCBI/hg19_refseq.20240523.txt ./NCBI/hg19.ncbiRefSeq.105.20190906.txt Priority_epimap_rts_alias_DisGeNET.xls Priority_epimap_rts_alias_DisGeNET_refseq.xls

# 5) Get numbers of coding & noncoding genes for Figure 2b venn diagram
awk -F"\t" '$9=="coding"' Priority_epimap_rts_alias_DisGeNET_refseq.xls |wc -l
703
awk -F"\t" '$9=="noncoding"' Priority_epimap_rts_alias_DisGeNET_refseq.xls |wc -l
213
awk -F"\t" '$9=="both"' Priority_epimap_rts_alias_DisGeNET_refseq.xls |wc -l
77

# 6) Extract go slim annotaton and gene description from BIomaRt using the below R codes (also listed in "TRIAGE_paper_R_source_code.R")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Set the working directory and make sure the file "Priority_epimap_rts.csv" is in the working directory
setwd("C:/Users/uqqzhao/UQ/1projects/TRIAGE_R_package")
file = read.csv("Priority_epimap_rts.csv", header=T, quote="")
genes <- file[which(file$Priority == "Y"),]
attributes_description <- c('hgnc_symbol', 'description')
annotations_description <- getBM(attributes = attributes_description, filters = 'hgnc_symbol', values = genes$Gene, mart = ensembl)
write.table(annotations_description, "anno_description.xls", sep="\t", row.names=FALSE, quote=FALSE)
attributes_go_slim <- c('hgnc_symbol', 'goslim_goa_accession', 'goslim_goa_description')
annotations_go_slim <- getBM(attributes = attributes_go_slim, filters = 'hgnc_symbol', values = genes$Gene, mart = ensembl)
write.table(annotations_go_slim, "anno_go_slim.xls", sep="\t", row.names=FALSE, quote=FALSE)

# 7) Since 43 out of 993 genes do not have descriptions from biomaRt, we manually searched for these genes using the GeneCards website (https://www.genecards.org/). The gene annotations for these 43 genes are recorded in "gene_annotation_GeneCards_manual.txt".

# 8) Add go slim annotation, gene description, GeneCards annotation 
anno_go_slim.pl anno_go_slim.xls anno_description.xls gene_annotation_GeneCards_manual.txt Priority_epimap_rts_alias_DisGeNET_refseq.xls Priority_epimap_rts_goslim.xls

# 9) Add NCBI gene summary using NCBI API
# Note this script take >30min to run - so suggest to using nohup or submit the job to computing clusters
nohup NCBI_summary.pl Priority_epimap_rts_goslim.xls Priority_epimap_rts_annotation.xls &

# 10) GWAS catalog annotation was downloaded from "https://fuma.ctglab.nl/downloadPage" -> "MSigDB and GWAS catalog gene-set files"
anno_FUMA.pl FUMA_GENE2FUNC_genesets_v156plus Priority_epimap_rts_annotation.xls Priority_genes_annotation.xls

# 11) Functional categories for 993 TRIAGE-prioritized genes
function_category.pl Priority_genes_annotation.xls Priority_genes_category.xls

# 12) Genereate Figure 1c & 1d using custom R codes - see R source codes in "TRIAGE_paper_R_source_code.R"



########################################################################################
###### Supplementary Figure 1. SuperPath and disease category enrichment analyses ######
########################################################################################

###### Source data are listed at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Supplementary_Figure1/"

### Gene pathway analysis using GeneCards - the GeneAnalytics tool
# Score: The entities in the table are listed in descending order of their matching scores.
# The SuperPath matching score is based on the binomial distribution and is classified by quality (high, medium or low), indicated by the color of the score bar (dark green, light green or beige, respectively). 
# Results generated by the GeneAnalytics tool were saved as "GeneCards_gene_analysis.xlsx"

# For generating Supplementary Figure 1, see R source codes at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/TRIAGE_paper_R_source_code.R"



#############################################################
######     Figure 2. Overview of the TRIAGE R package  ######
#############################################################
###### For Figure 2, part of the visuals were manually created with BioRender.com and integrated using Adobe Illustrator.


#############################################################
######              Figure 3.                          ######
######       In vivo mouse RNA-seq data analysis       ######
#############################################################

###### Source data are listed at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Figure3_mouse_RNAseq/GSE95755/"

# For the in vivo mouse RNA-seq data analysis, including TRIAGEgene analysis, and the generation of Figure 3, please refer to the R source code available at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/TRIAGE_paper_R_source_code.R"



#############################################################
######            Figure 4.                            ######
###### In vivo human single-cell RNA-seq data analysis ######
#############################################################

###### Source data are listed at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Figure4_human_scRNAseq/PBMC_ifnb_data/"

# For the in vivo human single-cell RNA-seq data analysis, including scRNA-seq data integration, TRIAGE analysis, and the generation of Figure 4, please refer to the R source code available at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/TRIAGE_paper_R_source_code.R". 
A custom Python script ("peak_correlation_and_plot.py" at https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Figure4_human_scRNAseq) was used to generate Figure 4c.



#############################################################
######          Supplementary Figure 2                 ######
#############################################################

###### Source data are listed at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Figure4_human_scRNAseq/PBMC_ifnb_data/"

# For the generation of Supplementary Figure 2, please refer to the R source code available at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/TRIAGE_paper_R_source_code.R" - at the "Supplementary Figure 2" section.



#############################################################
######            Figure 5.                            ######
######      In vitro human RNA-seq data analysis       ######
#############################################################

###### Source data are listed at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Figure5_human_RNAseq/GSE246079/"

# For the in vitro human RNA-seq data analysis, including TRIAGE analysis, comparisons of enriched GO terms, gene network analysis, and the generation of Figure 5, please refer to the R source code available at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/TRIAGE_paper_R_source_code.R".




#############################################################
######          Supplementary Figure 3                 ######
#############################################################

###### Source data are listed at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/Supplementary_Figure3/"

# For the generation of Supplementary Figure 3, please refer to the R source code available at "https://github.com/Qiongyi/TRIAGE_R_Package/source_data_scripts_for_paper/TRIAGE_paper_R_source_code.R".