User-friendly Functions
=======================

The TRIAGE R package includes several user-friendly functions designed to enhance data visualization and analysis. These functions provide intuitive graphical representations of complex data, making it easier for researchers to interpret their results. This document covers three key functions: `plotJaccard`, `byPeak`, and `plotGO`.

plotJaccard
-----------

The `plotJaccard` function generates Jaccard similarity index heatmaps based on the output from the `TRIAGEgene` function, allowing for intuitive data comparisons.

**Parameters:**

- `ds`: The output matrix from the `TRIAGEgene` function.

..

- `output_file`: The desired file name for the output heatmap PDF.

..

- `top_no`: (Optional) The number of top TRIAGE ranked genes to consider for the Jaccard index calculation. Default is 100.

**Usage Example:**

.. code-block:: R

    ds <- TRIAGEgene(input_matrix)
    plotJaccard(ds, "Jaccard_heatmap.pdf")


byPeak
------

The `byPeak` function calculates the average gene expression or average TRIAGE-weighted values for each gene grouped by 'Peak'. It supports direct data frame input or reading from CSV/TXT files.

**Parameters:**

- `expr`: The gene expression data, either as a data frame or a path to a CSV/TXT file.

..

- `peak`: The metadata containing cell IDs and peak values, either as a data frame or a path to a CSV/TXT file.

..

- `cell_column`: (Optional) Name of the column in metadata representing cell IDs. Default is "Barcode".

..

- `peak_column`: (Optional) Name of the column in metadata representing peak values. Default is "Peak".

..

- `prefix`: (Optional) The prefix for the column names in the output file. Default is NULL.


**Usage Example:**

.. code-block:: R

    # Example 1: using .csv files as the input files
    result <- byPeak(expr = "path/to/expression.csv", 
                    peak = "path/to/metadata.csv")

    # Example 2: using data frame ('expr_df' and 'metadata_df') as the input, 
    # grouped by 'Clusters', and cell IDs are in the "cell_name" column
    result <- byPeak(expr = expr_df, 
                    peak = metadata_df, 
                    peak_column="Clusters",
                    cell_column="cell_name")


topGenes
--------

The `topGenes` function extracts the top `n` genes with the highest values from each column of a given data frame or matrix. This function is useful for identifying the top genes based on DS values or gene expression values for each TRIAGE peak, cluster, or group.

**Parameters:**

- `ds`: A data frame or matrix where rows represent genes and columns represent TRIAGE peaks, cell clusters, or groups. The values should be DS values or gene expression values.
- `top_no` (optional): The number of top genes to extract for each column. The default value is 10.

**Returns:**

- A matrix where each element contains the top `n` genes in the format "GeneSymbol (DS value)" for each column of the input data.

**Version Added:**

- v1.1.4

**Usage Example:**

.. code-block:: R

    # Example: Extract the top 5 genes for each TRIAGE peak
    ds <- TRIAGEgene(input_matrix)
    top_genes <- topGenes(ds, top_no = 5)


getClusterGenes
---------------

The `getClusterGenes` function extracts genes assigned to a specific cluster based on the highest probability from the output of the `TRIAGEparser`. This is particularly useful for users conducting downstream analyses for individual gene clusters.

**Parameters:**

- `input_file`: The path to the CSV file containing genes and their cluster probabilities.

..

- `cluster_name`: The name of the cluster to extract genes from.

**Version Added:**

- v1.1.3

**Usage Example:**

.. code-block:: R


    # Example: Extract genes assigned to cluster1 from the TRIAGEparser output
    cluster1_genes <- getClusterGenes("TRIAGEparser_output/gene_clusters/output_gene_clusters.csv", "cluster1")



plotGO
------

The `plotGO` function creates GO enrichment heatmaps from the output of the `TRIAGEparser`. It visualizes the GO enrichment analysis results for specific groups or IDs.

**Parameters:**

- `indir`: The path to the output directory from `TRIAGEparser`.

..

- `outdir`: The directory where the generated heatmap PDF files will be saved.

..

- `id`: (Optional) Parameter to specify a particular group or ID for heatmap generation. Default is NULL (generates heatmaps for all groups/IDs).

..

- `color_palette`: (Optional) Parameter for custom heatmap color palette. Default is a gradient from light grey to red.

..

- `top_terms`: (Optional) The number of top GO terms for each gene cluster to include in the heatmap. Default is 10.

..

- `fdr`: (Optional) The FDR threshold for the heatmap visulization of TRIAGEparser results. Default is 0.01.

..

- `width`: (Optional) The width of the output PDF heatmap. Default is NULL, which uses default behavior of pdf().

..

- `height`: (Optional) The height of the output PDF heatmap. Default is NULL, which uses default behavior of pdf().


**Usage Example:**

.. code-block:: R

    # Example 1: Generate heatmaps for all groups/IDs
    plotGO(indir = "path/to/TRIAGEparser_output", 
        outdir = "path/to/heatmap_output")

    # Example 2: Generate heatmap for a specific group “Peak1”, 
    # with the PDF size 6X7
    plotGO(indir = "path/to/TRIAGEparser_output", 
        outdir = "path/to/heatmap_output", 
        id = "Peak1",
        width=6, height=7)

    # Example 3: Generate heatmap for two groups: “Peak0” and "Peak1", 
    # with the PDF size 6X7
    plotGO(indir = "path/to/TRIAGEparser_output", 
        outdir = "path/to/heatmap_output", 
        id = c("Peak0", "Peak1"),
        width=6, height=7)