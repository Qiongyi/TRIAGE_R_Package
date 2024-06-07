TRIAGEcluster
=============

Description
-----------
TRIAGEcluster is one of core functions of the TRIAGE R package, representing a notable advancement in scRNA-seq data analysis by integrating epigenomic signatures to demarcate cell diversity within heterogeneous data. Utilizing a combination of genes with high RTS with weighted kernel density estimation in a two-dimensional space, TRIAGEcluster can be leveraged to refine cell clusters and identify biologically distinct cell populations, which we refer to as “peaks” in our method. For more details, see: Sun et al., Nucleic Acid Research 2023, "Inferring cell diversity in single cell data using consortium-scale epigenetic data as a biological anchor for cell identity".
For more details, see: `Sun et al., Nucleic Acid Research 2023, Inferring cell diversity in single cell data using consortium-scale epigenetic data as a biological anchor for cell identity <https://academic.oup.com/nar/article/51/11/e62/7147502>`_.



Input and Output
----------------

Input: TRIAGEcluster requires two input files. 
    1. A normalized gene expression matrix file or a TRIAGE-weighted file. This matrix file can be in either .csv or tab-delimited .txt format.
    2. A metadata file for scRNA-seq data analysis. The metadata file needs to contain information about cell identifiers and UMAP coordinates and can also be in .csv or tab-delimited .txt format.

Output: TRIAGEcluster generates a total of 18 output files in the specified directory. These outputs include nine UMAP plots, each corresponding to a different bandwidth resolution ranging from 0.1 to 0.9. These plots visually represent the clustering of cell populations within the scRNA-seq data. Alongside the UMAP plots, there are nine corresponding metadata files. Each of these files details the characteristics of the identified cell populations, referred to as "Peaks," for each bandwidth resolution.


Parameters
----------
- `expr`: The input file for the normalized gene expression matrix or TRIAGE-weighted matrix. Accepts .csv or tab-delimited .txt files.

..

- `metadata`: The metadata file for scRNA-seq data analysis, also in .csv or tab-delimited .txt format.

..

- `outdir`: (Optional) Specifies the output directory for the analysis results. Default is "TRIAGEcluster_results".

..

- `output_prefix`: (Optional) The prefix for output files. Default is "TRIAGEcluster".

..

- `cell_column`: (Optional) Indicates the column name in the metadata file representing cell identifiers. Default is "Barcode".

..

- `umap_column`: (Optional) Specifies the prefix for UMAP coordinate columns. Default is "UMAP\_".

..

- `priority_rts`: (Optional) Specifies the path to the priority RTS gene list file. The default setting is "Priority_epimap_rts.csv". Note that this option is intended for advanced users who are capable of generating their own priority RTS gene lists. For most users, it is recommended to use the default setting to ensure optimal functionality and results.

..

- `min_cells_per_peak`: (Optional) Sets the minimum number of cells required per peak to be considered in the analysis. Default is 5.

..

- `seed`: (Optional) Sets the seed value for reproducibility. Default is NULL.

Usage Examples
--------------

TRIAGEcluster can be run using various combinations of parameters, as shown in the following examples:

.. code-block:: R

    # Example 1: Using a .csv TRIAGE-weighted (i.e. discordance score) matrix 
    # and a tab-delimited .txt metadata file.
    TRIAGEcluster("ds.csv", 
                "metadata.txt", 
                outdir = "TRIAGEcluster_results", 
                output_prefix = "project1")

    # Example 2: Using a tab-delimited .txt file for the discordance score matrix 
    # and a .csv file for the metadata. The cell identifiers are specified in the 
    # "Cells" column and the UMAP info are named in "UMAP1" and "UMAP2" columns 
    # in the metadata file. Set the seed value to 123.
    TRIAGEcluster(expr = "ds.txt", 
              metadata = "metadata.csv", 
              output_prefix = "project2", 
              cell_column = "Cells",  
              umap_column = "UMAP",
              seed=123)

    # Example 3: Using .csv files for both gene expression matrix and metadata.
    TRIAGEcluster("Expr_matrix.csv", 
                "metadata.csv", 
                output_prefix = "project3")

    # Example 4: Using tab-delimited .txt files for both gene expression matrix 
    # and metadata.
    TRIAGEcluster("Expr_matrix.txt", 
                "metadata.txt", 
                outdir = "results", 
                output_prefix = "project4")

