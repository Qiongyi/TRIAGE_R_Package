Testing TRIAGE
==============

The TRIAGE R package offers a comprehensive suite of tools for analyzing transcriptomic data. This document provides a guide to testing the key functionalities of TRIAGE, including `TRIAGEgene`, `TRIAGEcluster`, and `TRIAGEparser`, along with their associated visualization and analysis functions. These tests are designed to demonstrate the capabilities of each function and ensure their correct operation.

Testing TRIAGEgene + plotJaccard()
----------------------------------

`TRIAGEgene` is used for gene-level analysis and generating TRIAGE-weighted gene expression data.

### Test 1: Running TRIAGEgene on Demo Human Data

Objective: To test `TRIAGEgene` using human data and generate a Jaccard Index Heatmap for visualization.

**Steps:**

1. Read the input file (tab delimited .txt file).
2. Run `TRIAGEgene` (Auto-selection for log transformation is enabled).
3. Generate Jaccard Index Heatmap based on `TRIAGEgene` output (using top 100 genes by default).

.. code-block:: R

    library(Triage)
    # Read input file
    input_file <- system.file("extdata", "TRIAGEgene_demo_Human.txt", package = "Triage")
    demo <- read.table(input_file, header = TRUE, sep = "\t", quote = "", row.names = 1)

    # Run TRIAGEgene
    ds <- TRIAGEgene(demo)

    # Generate Jaccard Index Heatmap
    setwd("/path/to/working/directory")
    if (!dir.exists("tests")) {
      dir.create("tests")
    }
    plotJaccard(ds, "tests/Jaccard_heatmap_Human_demo.pdf")


### Test 2: Running TRIAGEgene on Demo Mouse Data

Objective: To test `TRIAGEgene` using mouse data and generate a Jaccard Index Heatmap for visualization.

**Steps:**

1. Read the input file (CSV format).
2. Run `TRIAGEgene` with species specified as "Mouse". Auto-selection for log transformation is enabled.
3. Generate a Jaccard Index Heatmap using the top 100 genes.

.. code-block:: R

    library(Triage)
    # Read input file (CSV)
    input_file <- system.file("extdata", "TRIAGEgene_demo_Mouse.csv", package = "Triage")
    demo <- read.csv(input_file, row.names = 1)

    # Run TRIAGEgene for Mouse data
    ds <- TRIAGEgene(demo, species = "Mouse")

    # Generate Jaccard Index Heatmap
    plotJaccard(ds, "tests/Jaccard_heatmap_Mouse_demo.pdf", top_no = 100)


### Test 3: Running TRIAGEgene on Mouse Data with Matrix Input

Objective: To evaluate the functionality of `TRIAGEgene` using mouse data in matrix format and generate a Jaccard Index Heatmap for visualization.

**Steps:**

1. Read the input file and convert it to a matrix (CSV format).
2. Run `TRIAGEgene` with matrix input, specifying "Mouse" as the species.
3. Generate a Jaccard Index Heatmap using the top 88 genes.

.. code-block:: R

    library(Triage)
    # Read input file (CSV) and convert to matrix
    input_file <- system.file("extdata", "TRIAGEgene_demo_Mouse.csv", package = "Triage")
    demo <- read.csv(input_file, row.names = 1)
    demo_matrix <- as.matrix(demo)

    # Run TRIAGEgene with matrix input for Mouse data
    ds <- TRIAGEgene(demo_matrix, species = "Mouse")

    # Generate Jaccard Index Heatmap
    plotJaccard(ds, "tests/Jaccard_heatmap_Mouse_demo2.pdf", top_no = 88)


Testing TRIAGEcluster + byPeak() + TRIAGEgene() + plotJaccard()
---------------------------------------------------------------

`TRIAGEcluster` is used for refining cell clustering in scRNA-seq data. 

### Test 4: Running TRIAGEcluster and TRIAGEgene on Human Data
# Data Source: 
# - This test uses a publicly available single-nucleus RNA sequencing (snRNA-seq) dataset from Kramann et al., Nature, 2022, titled "Spatial multi-omic map of human myocardial infarction".
# - Dataset URL: https://www.nature.com/articles/s41586-022-05060-x#data-availability
# - The demonstration uses control cells to showcase the TRIAGE analysis pipeline.

Objective: To use `TRIAGEcluster` for cell clustering and `TRIAGEgene` for analyzing average expression data by peak.

**Steps:**

1. Run `TRIAGEcluster` for Cell Clustering (using CSV files for expression data and metadata).
2. Select the Most Suitable Bandwidth (manual review of plots).
3. Calculate Average Gene Expression by Peak (using `byPeak` function).
4. Save Results and Generate Jaccard Index Heatmap.

.. code-block:: R

    library(Triage)
    library(reticulate)
    setwd("/path/to/working/directory")
    # Run TRIAGEcluster
    expr_file <- system.file("extdata", "TRIAGEcluster_demo_expr_human.csv", package = "Triage")
    metadata_file <- system.file("extdata", "TRIAGEcluster_demo_metadata_human.csv", package = "Triage")
    TRIAGEcluster(expr_file, metadata_file, outdir = "tests", output_prefix = "demo")

    # Select suitable bandwidth and calculate average expression
    peak_file <- "tests/demo_bw0.80_metadata.csv"
    result <- byPeak(expr_file, peak_file)

    # Save and generate heatmap
    write.csv(result, file = "tests/AverageByPeak.csv", row.names = TRUE, quote = FALSE)
    ds <- TRIAGEgene(result)
    plotJaccard(ds, "tests/Jaccard_heatmap_peak.pdf")


### Alternative Calculations 1: Average Gene Expression by Cluster

Objective: To calculate average gene expression based on cluster categories using `byPeak` function, followed by `TRIAGEgene` analysis and Jaccard index heatmap generation.

**Steps:**

1. Calculate Average Gene Expression by Cluster.
2. Save Results.
3. Run `TRIAGEgene` on the results.
4. Generate Jaccard Index Heatmap for cluster-based averages.

.. code-block:: R

    # Calculate average expression by cluster
    result2 <- byPeak(expr_file, peak_file, peak_col = "final_cluster")

    # Save results
    write.csv(result2, file = "tests/AverageByCluster.csv", row.names = TRUE, quote = FALSE)
    write.table(result2, file = "tests/AverageByCluster.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

    # Run TRIAGEgene and generate Jaccard index heatmap
    ds2 <- TRIAGEgene(result2)
    plotJaccard(ds2, "tests/Jaccard_heatmap_cluster.pdf")


### Alternative Calculations 2: Average Gene Expression by Cell Type

Objective: To calculate average gene expression based on cell type categories using `byPeak` function, followed by `TRIAGEgene` analysis and Jaccard index heatmap generation.

**Steps:**

1. Calculate Average Gene Expression by Cell Type.
2. Save Results.
3. Run `TRIAGEgene` on the results.
4. Generate Jaccard Index Heatmap for cell type-based averages.

.. code-block:: R

    # Calculate average expression by cell type
    result3 <- byPeak(expr_file, peak_file, peak_col = "cell_type")

    # Save results
    write.csv(result3, file = "tests/AverageByCelltype.csv", row.names = TRUE, quote = FALSE)
    write.table(result3, file = "tests/AverageByCelltype.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

    # Run TRIAGEgene and generate Jaccard index heatmap
    ds3 <- TRIAGEgene(result3)
    plotJaccard(ds3, "tests/Jaccard_heatmap_celltype.pdf")


Testing TRIAGEparser + plotGO()
-------------------------------

`TRIAGEparser` is a machine learning-based method for evaluating gene expression rank lists.

### Test 5: Running TRIAGEparser with "AverageByPeak.csv"

Objective: To demonstrate `TRIAGEparser` functionality using a CSV file with four peak clusters.

**Steps:**

1. Run `TRIAGEparser`.
2. Generate GO Heatmaps for All Groups.

.. code-block:: R

    library(Triage)
    library(reticulate)
    # Run TRIAGEparser
    input_file <- "tests/AverageByPeak.csv"
    TRIAGEparser(input_file, input_type = "table", outdir="tests/TRIAGEparser_test5")

    # Generate Heatmaps
    plotGO(indir="tests/TRIAGEparser_test5", outdir="tests/TRIAGEparser_test5")


### Test 6: Running TRIAGEparser with "AverageByPeak.txt"

Objective: To demonstrate `TRIAGEparser` functionality using a tab-delimited text file and generate a specific gene group heatmap.

**Steps:**

1. Run `TRIAGEparser` with tab-delimited text file input.
2. Generate GO Heatmap for the "Peak0" group.

.. code-block:: R

    library(Triage)
    library(reticulate)
    # Run TRIAGEparser
    input_file <- "tests/AverageByPeak.txt"
    TRIAGEparser(input_file, input_type = "table", outdir="tests/TRIAGEparser_test6")

    # Generate heatmap for "Peak0" group
    plotGO(indir="tests/TRIAGEparser_test6", outdir="tests/TRIAGEparser_test6", id = "Peak0")


### Test 7: Running TRIAGEparser with Gene List

Objective: To test `TRIAGEparser` using a gene list and visualize gene ontology enrichment.

**Steps:**

1. Run `TRIAGEparser` with a gene list file as input.
2. Generate Gene Ontology Heatmap.

.. code-block:: R

    # Run TRIAGEparser with gene list file
    input_file <- system.file("extdata", "TRIAGEparser_demo_genelist.txt", package = "Triage")
    TRIAGEparser(input_file, input_type = "list", outdir="tests/TRIAGEparser_test7")

    # Generate Gene Ontology Heatmap
    plotGO(indir="tests/TRIAGEparser_test7", outdir="tests/TRIAGEparser_test7")


These tests serve as a practical demonstration of how to apply the TRIAGE R package for analyzing and visualizing complex transcriptomic data. Researchers can adapt these procedures to their specific datasets, ensuring the effective use of TRIAGE in research projects.