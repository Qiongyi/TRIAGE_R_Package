TRIAGEparser
============

Description
-----------
TRIAGEparser is one of core functions of the TRIAGE R package, designed to evaluate groups of genes, such as the top 100 genes ranked by TRIAGE-weighted values or differentially expressed genes, to classify genes with distinct biological functions. It performs principal component analysis to extract orthogonal patterns of H3K27me3 depositions from consortium-level epigenomic data and uses Bayesian information criterion to optimally determine gene clusters. TRIAGEparser then assesses each gene cluster by searching the protein-protein interaction (PPI) networks from the STRING database and conducts Gene Ontology (GO) enrichment analysis for genes with direct PPI interactions. 
For more details, see: `Sun et al., Nucleic Acid Research 2023, Inferring cell diversity in single cell data using consortium-scale epigenetic data as a biological anchor for cell identity <https://academic.oup.com/nar/article/51/11/e62/7147502>`_.

**Note:** TRIAGEparser is adaptable to any type of data mapped to protein-coding and non-coding genes, including RNAseq, proteomics, ChIP-seq, and more.


Input and Output
----------------

Input: TRIAGEparser requires an input file, which can be provided in two formats:

As a *Gene List*: A list of genes, typically in a text file - each line contains one gene name. This format is suitable when you want to analyze a specific set of genes.

As a *Table*: A more comprehensive data table, either in .csv or tab/space-delimited .txt format. This format is ideal for analyzing gene expression data along with other associated data points.


Output: The output from TRIAGEparser are two folders, "gene_clusters" and "go".

In the "gene_clusters" folder, there are "\*_gene_clusters.csv" files listing the probabilities of each gene being assigned to different gene clusters. For analyses involving multiple samples/groups, outputs are stored in distinct files. 

In the "go" folder, there are "\*_go.txt" files listing significance values (i.e., false discovery rates) for all associated GO terms descriptions across PPI-significant clusters. For analyses involving multiple samples/groups, outputs are stored in distinct files. 



Parameters
----------

- `input`: The input file, which can be a .csv file or a tab/space-delimited .txt file.

..

- `input_type`: (Optional) Specifies the input type, either 'table' or 'list'. Default is 'list'.

..

- `outdir`: (Optional) The path to the output directory. Default is 'TRIAGEparser_output'.

..

- `H3K27me3_pc`: (Optional) The pre-calculated H3K27me3 principal components. Default is 'pca_roadmap'.

..

- `number_of_pca`: (Optional) Number of principal components to use. Default is 10.

..

- `number_of_gene`: (Optional) Number of top genes to use if the input type is a table. Default is 100. If the input type is a list, all genes in the list will be used.

..

- `no_iter`: (Optional) Number of iterations for determining the best number of clusters using Bayesian Information Criterion (BIC). Default is 100.

..

- `EM_tol`: (Optional) Convergence threshold for the Expectation-Maximization (EM) iterations in the GaussianMixture function. Default is 1e-3.

..

- `EM_max_iter`: (Optional) Maximum number of EM iterations for the GaussianMixture function. Default is 100.

..

- `go_analysis`: (Optional) Option to perform GO enrichment analysis. (1: Yes, 0: No). Default is 1.

..

- `verbose`: (Optional) Level of verbosity (options: 1 or 0). Default is 1.

..

- `max_cluster`: (Optional) Maximum number of clusters to consider. Default is 10.

..

- `gene_order`: (Optional) Direction to sort genes (options: 'ascending' or 'descending'). Default is 'descending'.

..

- `go_threshold`: (Optional) Threshold for GO term enrichment (False Discovery Rate). Default is 0.01.



Usage Examples
--------------

.. code-block:: R

    # Example 1: Using a tab-delimited table file "input.txt" as the input 
    # and "TRIAGEparser_output" as the output directory
    TRIAGEparser("input.txt", input_type = "table")

    # Example 2: Using "input.txt" - a gene list as the input, 
    # and specifying the output directory
    TRIAGEparser("input.txt", outdir = "path/to/results")

    # Example 3: Using a CSV file "input.csv" and specifying the output 
    # directory. Using top 200 genes for the TRIAGEparser analysis.
    TRIAGEparser("input.csv", 
            input_type = "table", 
            outdir = "path/to/results", 
            number_of_gene = 200)
