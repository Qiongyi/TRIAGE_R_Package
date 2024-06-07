User-friendly Functions
=======================

The TRIAGE R package includes several user-friendly functions designed to enhance data visualization and analysis. These functions provide intuitive graphical representations of complex data, making it easier for researchers to interpret their results. This document covers three key functions: `plotJaccard`, `byPeak`, and `plotGO`.

plotJaccard Function
--------------------

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


byPeak Function
---------------

The `byPeak` function calculates the average gene expression or average TRIAGE-weighted values for each gene grouped by 'Peak'. It supports direct data frame input or reading from CSV/TXT files.

**Parameters:**

- `expr`: The gene expression data, either as a data frame or a path to a CSV/TXT file.

..

- `peak`: The metadata containing cell IDs and peak values, either as a data frame or a path to a CSV/TXT file.

..

- `cell_column`: (Optional) Name of the column in metadata representing cell IDs. Default is "Barcode".

..

- `peak_column`: (Optional) Name of the column in metadata representing peak values. Default is "Peak".

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


plotGO Function
---------------

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

- `top_terms`: (Optional) The number of top GO terms to include in the heatmap. Default is 10.
..

- `width`: (Optional) The width of the output PDF heatmap. Default is NULL, which uses default behavior of pdf().
..

- `height`: (Optional) The height of the output PDF heatmap. Default is NULL, which uses default behavior of pdf().


**Usage Example:**

.. code-block:: R

    # Example 1: Generate heatmaps for all groups/IDs
    plotGO(indir = "path/to/TRIAGEparser_output", 
        outdir = "path/to/heatmap_output")

    # Example 2: Generate heatmap for a specific group “Peak01”, with the PDF size 6X7
    plotGO(indir = "path/to/TRIAGEparser_output", 
        outdir = "path/to/heatmap_output", 
        id = "Peak01",
        width=6, height=7)
