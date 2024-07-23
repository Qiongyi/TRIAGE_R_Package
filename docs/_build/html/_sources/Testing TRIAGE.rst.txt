Testing TRIAGE
==============

The TRIAGE R package offers a comprehensive suite of tools for analyzing transcriptomic data. This document provides a guide to testing the key functionalities of TRIAGE, including `TRIAGEgene`, `TRIAGEcluster`, and `TRIAGEparser`, along with their associated visualization and analysis functions. These tests are designed to demonstrate the capabilities of each function and ensure their correct operation.

Test TRIAGEgene + plotJaccard()
----------------------------------

`TRIAGEgene` is used for gene-level analysis and generating TRIAGE-weighted gene expression data.

**# Test 1: Run TRIAGEgene on Demo Human Data**

Objective: To test `TRIAGEgene` using human data and generate a Jaccard Index Heatmap for visualization.

**Steps:**

1. Read the input file (tab delimited .txt file).
2. Run `TRIAGEgene` (Auto-selection for log transformation is enabled).
3. Generate Jaccard Index Heatmap based on `TRIAGEgene` output (using top 100 genes by default).

.. code-block:: R

    library(TRIAGE)
    # Read input file
    input_file <- system.file("extdata", "TRIAGEgene_demo_Human.txt", package = "TRIAGE")
    demo <- read.table(input_file, header = TRUE, sep = "\t", quote = "", row.names = 1)

    # Run TRIAGEgene
    ds <- TRIAGEgene(demo)

    # Generate Jaccard Index Heatmap
    setwd("/path/to/working/directory")
    if (!dir.exists("tests")) {
      dir.create("tests")
    }
    plotJaccard(ds, "tests/Jaccard_heatmap_Human_test1.pdf")


**# Test 2: Run TRIAGEgene on Demo Mouse Data**

Objective: To test `TRIAGEgene` using mouse data and generate a Jaccard Index Heatmap for visualization.

**Steps:**

1. Read the input file (CSV format).
2. Run `TRIAGEgene` with species specified as "Mouse". Auto-selection for log transformation is enabled.
3. Generate a Jaccard Index Heatmap using the top 100 genes.

.. code-block:: R

    library(TRIAGE)
    # Read input file (CSV)
    input_file <- system.file("extdata", "TRIAGEgene_demo_Mouse.csv", package = "TRIAGE")
    demo <- read.csv(input_file, row.names = 1)

    # Run TRIAGEgene for Mouse data
    ds <- TRIAGEgene(demo, species = "Mouse")

    # Generate Jaccard Index Heatmap
    plotJaccard(ds, "tests/Jaccard_heatmap_Mouse_test2.pdf", top_no = 100)


**# Test 3: Run TRIAGEgene on Mouse Data with Matrix Input**

Objective: To evaluate the functionality of `TRIAGEgene` using mouse data in matrix format and generate a Jaccard Index Heatmap for visualization.

**Steps:**

1. Read the input file and convert it to a matrix (CSV format).
2. Run `TRIAGEgene` with matrix input, specifying "Mouse" as the species.
3. Generate a Jaccard Index Heatmap using the top 88 genes.

.. code-block:: R

    library(TRIAGE)
    # 1) Read input file (CSV) and convert to matrix
    input_file <- system.file("extdata", "TRIAGEgene_demo_Mouse.csv", package = "TRIAGE")
    demo <- read.csv(input_file, row.names = 1)
    demo_matrix <- as.matrix(demo)

    # 2) Run TRIAGEgene with matrix input for Mouse data
    ds <- TRIAGEgene(demo_matrix, species = "Mouse")

    # 3) Generate Jaccard Index Heatmap
    plotJaccard(ds, "tests/Jaccard_heatmap_Mouse_test3.pdf", top_no = 88)


Test TRIAGEcluster + byPeak()
------------------------------

`TRIAGEcluster` is used for refining cell clustering in scRNA-seq data.

**# Test 4: Run TRIAGEcluster and TRIAGEgene on Human Data**

Objective: To use `TRIAGEcluster` for cell clustering, `byPeak()` for analyzing average expression data by peak, and `TRIAGEgene` for generating TRIAGE-weighted expression data (DS).

**Steps:**

1. Run `TRIAGEcluster` for Cell Clustering, using CSV files for expression data and metadata.
2. Select a suitable Bandwidth based on UMAP reviews and Calculate Average Gene Expression by Peak using the `byPeak()` function.
3. Run `TRIAGEgene` to generate TRIAGE-weighted expression data.

.. code-block:: R

    library(TRIAGE)
    library(reticulate)
    setwd("/path/to/working/directory")
    
    # 1) Run TRIAGEcluster
    expr_file <- system.file("extdata", "TRIAGEcluster_demo_expr_human.csv", package = "TRIAGE")
    metadata_file <- system.file("extdata", "TRIAGEcluster_demo_metadata_human.csv", package = "TRIAGE")
    TRIAGEcluster(expr_file, metadata_file, outdir = "tests/test4", output_prefix = "demo")

    # 2) Select a suitable bandwidth and calculate average gene expression
    peak_file <- "tests/test4/demo_bw0.80_metadata.csv"
    avg_peak <- byPeak(expr_file, peak_file, cell_column = "Barcode", peak_column = "Peak")
    # Save the average gene expression result to a CSV file
    write.csv(avg_peak, file = "tests/test4/AverageByPeak.csv", row.names = TRUE, quote = FALSE)

    # 3) Run TRIAGEgene to generate TRIAGE-weighted expression data (DS)
    ds <- TRIAGEgene(avg_peak)
    # Save the average DS result to a CSV file
    write.csv(ds, file = "tests/test4/AverageByPeak_DS.csv", row.names = TRUE, quote = FALSE)
    # Save the average DS result to a tab-delimited text file
    write.table(ds, file = "tests/test4/AverageByPeak_DS.txt", sep = "\t", 
                row.names = TRUE, col.names = NA, quote = FALSE)



Test TRIAGEparser + plotGO() + getClusterGenes()
-------------------------------------------------

`TRIAGEparser` is a machine learning-based method for classifying genes with distinct biological functions.

**# Test 5: Run TRIAGEparser with "AverageByPeak_DS.csv"**

Objective: To demonstrate `TRIAGEparser` functionality using a CSV file with four peak clusters.

**Steps:**

1. Run `TRIAGEparser`.
2. Generate GO Heatmaps for All Groups.
3. Extract genes for cluster1 from the "Peak0_gene_clusters.csv" output of TRIAGEparser.

.. code-block:: R

    library(TRIAGE)
    library(reticulate)
    # 1) Run TRIAGEparser with "AverageByPeak_DS.csv" generated in Test 4
    input_file <- "tests/test4/AverageByPeak_DS.csv"
    TRIAGEparser(input_file, input_type = "table", outdir="tests/test5")

    # 2) Generate Heatmaps using plotGO()
    plotGO(indir="tests/test5", outdir="tests/test5")

    # 3) Extract genes for cluster1 from the "Peak0_gene_clusters.csv" using getClusterGenes()
    cluster1_genes <- getClusterGenes("tests/test5/gene_clusters/Peak0_gene_clusters.csv", "cluster1")

**# Test 6: Run TRIAGEparser with "AverageByPeak_DS.txt"**

Objective: To demonstrate `TRIAGEparser` functionality using a tab-delimited text file and generate a specific gene group heatmap.

**Steps:**

1. Run `TRIAGEparser` with tab-delimited text file input.
2. Generate GO Heatmap for the "Peak0" group.

.. code-block:: R

    library(TRIAGE)
    library(reticulate)
    # 1) Run TRIAGEparser with "AverageByPeak_DS.txt" generated in Test 4
    input_file <- "tests/test4/AverageByPeak_DS.txt"
    TRIAGEparser(input_file, input_type = "table", outdir="tests/test6")

    # 2) Generate heatmap for "Peak0" group
    plotGO(indir="tests/test6", outdir="tests/test6", id = "Peak0")


**# Test 7: Run TRIAGEparser with a Gene List**

Objective: To test `TRIAGEparser` using a gene list and visualize gene ontology enrichment.

**Steps:**

1. Run `TRIAGEparser` with a gene list file as input.
2. Generate Gene Ontology Heatmap.

.. code-block:: R

    # 1) Run TRIAGEparser with gene list file
    input_file <- system.file("extdata", "TRIAGEparser_demo_genelist.txt", package = "TRIAGE")
    TRIAGEparser(input_file, input_type = "list", outdir="tests/test7")

    # 2) Generate Gene Ontology Heatmap
    plotGO(indir="tests/test7", outdir="tests/test7")


These tests serve as a practical demonstration of how to apply the TRIAGE R package for analyzing and visualizing complex transcriptomic data. Researchers can adapt these procedures to their specific datasets, ensuring the effective use of TRIAGE in research projects.
