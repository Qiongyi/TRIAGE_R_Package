TRIAGEgene
==========

Description
-----------
TRIAGEgene is one of core functions of the TRIAGE R package, aiming to predict the regulatory potential of genes and further identify genetic drivers of cell identity. This approach calculates repressive tendency scores (RTS) for each gene by analyzing broad H3K27me3 domains near the gene (2.5kb upstream plus the gene body), and then integrates the RTS metric with gene expression data to calculate a TRIAGE-weighted value for each gene. This value, also referred to as Discordance Score (DS) in previous literature, indicates genes’ regulatory potential. After this TRIAGEgene transformation, it becomes instrumental in identifying potential regulatory and cell identity genes by ranking them in descending order based on the TRIAGE-weighted values within each group, which can vary from specific conditions to individual samples, distinct cell types, clusters, or even single cells. 
For more details, see: `Shim et al., Cell Systems 2020, Conserved Epigenetic Regulatory Logic Infers Genes Governing Cell Identity <https://linkinghub.elsevier.com/retrieve/pii/S2405-4712(20)30419-1>`_.

**Note:** TRIAGEgene is adaptable to any type of data mapped to protein-coding and non-coding genes, including RNAseq, proteomics, ChIP-seq, and more.



Input and Output
----------------

Input: The function requires a matrix or data frame of normalized gene expression data. This data can be in formats like Counts Per Million (CPM), Fragments Per Kilobase of transcript per Million mapped reads (FPKM), or Transcripts Per Million (TPM).

Output: The output is a matrix or data frame of TRIAGE-weighted gene expression data. This TRIAGE transformation can help to identify regulatory genes and genes crucial for cell identity.


Parameters
----------

- `m`: Input matrix or data frame of normalized gene expression data. Acceptable types include CPM, FPKM, TPM, among others.

..

- `species`: (Optional) Specifies the species. Default is "Human". Other options include "C.intestinalis", "Chicken", "Guinea Pig", "Mouse", "Pig", "Zebrafish".

..

- `log`: (Optional) Determines whether to apply a log transformation to the input data. The default value is NULL, allowing the function to make a decision based on data characteristics. Generally, it is recommended to use natural log-transformed normalized gene expression data as the input for TRIAGEgene. This transformation often enhances the analysis accuracy and is preferable for most datasets.

..

- `data_source`: (Optional) Data source selection, either "epimap" (default) or "roadmap". Originally utilizing H3K27me3 data from 111 cell types in the NIH Epigenome Roadmap dataset (Roadmap Epigenomics, et al., 2015), the current release of TRIAGEgene has expanded its scope to include the EpiMap dataset, which offers a more comprehensive H3K27me3 signatures across 833 cell and tissue types (`Boix, et al., 2021 <https://www.nature.com/articles/s41586-020-03145-z>`_). Users can choose either the default EpiMap dataset or the original Roadmap dataset for RTS calculations, ensuring backward compatibility and data reproducibility.

Usage Examples
--------------

.. code-block:: R

    # Example 1: Human data, with auto log transformation decision
    result <- TRIAGEgene(input)

    # Example 2: Human data, 'roadmap' data source, auto log transformation
    result <- TRIAGEgene(input, data_source = "roadmap")

    # Example 3: Human data, forced log transformation
    result <- TRIAGEgene(input, log = TRUE)

    # Example 4: Mouse data, without log transformation
    result <- TRIAGEgene(input, species = "Mouse", log = FALSE)
