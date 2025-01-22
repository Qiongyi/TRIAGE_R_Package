TRIAGEgene
==========

Description
-----------
TRIAGEgene is one of core functions of the TRIAGE R package, aiming to predict the regulatory potential of genes and further identify genetic drivers of cell identity. This approach calculates repressive tendency scores (RTS) for each gene by analyzing broad H3K27me3 domains near the gene (2.5kb upstream plus the gene body), and then integrates the RTS metric with gene expression data to calculate a TRIAGE-weighted value for each gene. This value, also referred to as Discordance Score (DS) in previous literature, indicates genes’ regulatory potential. After this TRIAGEgene transformation, it becomes instrumental in identifying potential regulatory and cell identity genes by ranking them in descending order based on the TRIAGE-weighted values within each group, which can vary from specific conditions to individual samples, distinct cell types, clusters, or even single cells. 
For more details, see: `Zhao et al., Briefings in Bioinformatics 2025, TRIAGE: an R package for regulatory gene analysis <https://doi.org/10.1093/bib/bbaf004>`_ and `Shim et al., Cell Systems 2020, Conserved Epigenetic Regulatory Logic Infers Genes Governing Cell Identity <https://pmc.ncbi.nlm.nih.gov/articles/PMC7781436/>`_.

The conservation of H3K27me3 patterns across eukaryotes (`Arthur et al., Genome Research, 2014 <https://pmc.ncbi.nlm.nih.gov/articles/PMC4079967/>`_) allows for the application of RTS calculation for a gene across species by identifying orthologous genes, even though the RTS values were originally generated using human H3K27me3 data. Our original publication validated this approach (`Shim et al., Cell Systems, 2020 <https://pmc.ncbi.nlm.nih.gov/articles/PMC7781436/>`_). In this R package, we utilize Ensembl BioMart to retrieve orthologous genes between humans and other species, enabling the application of RTS values from human genes to corresponding orthologs in non-human datasets. Detailed online instructions for retrieving orthologous genes can also be found here: https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/.

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

..

- `pvalue`: (Optional) Logical value indicating whether to calculate p-values using a rank-based Z-Score method. Default is False. If set to TRUE, the function will calculate p-values based on the rank of each gene's DS value compared to comparable genes with similar expression values.

..

- `min_comparable_genes`: (Optional) Integer specifying the minimum number of comparable genes to use for calculating p-values when `pvalue = TRUE`. This parameter ensures that at least this number of genes is used in the rank-based Z-Score calculation. Default is 100. Only applicable when `pvalue` is set to TRUE.

..

- `percentile`: (Optional) Numeric value representing the proportion of input genes to be used as comparable genes with similar gene expression values for p-value calculation when `pvalue = TRUE`. Default is 0.1 (10% of the input genes). If the calculated proportion results in fewer genes than `min_comparable_genes`, the `min_comparable_genes` value will be applied. Only applicable when `pvalue` is set to TRUE.

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

    # Example 5: Mouse data, calculate p-values using a rank-based Z-Score method
    result <- TRIAGEgene(input, species = "Mouse", pvalue = TRUE)