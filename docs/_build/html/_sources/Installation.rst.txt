Installation
============

This section covers the installation and setup of the TRIAGE R package. The TRIAGE R package requires both R and Python dependencies. Follow these steps to set up the environment and install the package.

.. _installation:


R Dependencies
--------------

Install the required R packages:

.. code-block:: R

    install.packages("reticulate")
    install.packages("data.table")
    install.packages("pheatmap")


**Additional Packages for** `compareGO`

To enable the `compareGO` function, the following R packages are required:

.. code-block:: R

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install(c("clusterProfiler", "enrichplot"))
    install.packages(c("dplyr", "ggplot2", "reshape2"))


These packages are only required for `compareGO`, which performs GO enrichment analysis on gene lists. Note it also requires an organism-specific annotation database (e.g., `org.Hs.eg.db` for human), which can be installed with:

.. code-block:: R

    BiocManager::install("org.Hs.eg.db")

Common options include `org.Hs.eg.db` for humans, `org.Mm.eg.db` for mice, and `org.Rn.eg.db` for rats. For additional species databases, see the [`Bioconductor Annotation Packages <https://bioconductor.org/packages/release/BiocViews.html#___OrgDb>`_].

Python Environment Setup
------------------------

To install the TRIAGE R package, you'll first need to set up the required Python modules within R using the reticulate package. There are two main options for this setup: using a Conda environment or a Python virtual environment. Details for each option are provided below.

**Option 1: Using a Conda environment**

To set up a Python environment with Conda for use in R, follow the steps below. This approach is distinct from the Python virtual environment (pyenv) option (Option 2) outlined later. The Conda option provides more comprehensive package management and may be preferable for complex projects.

# Bash commands for setting up the Conda environment "r-env":

.. code-block:: bash

    # List all existing Conda environments to see if any relevant environments have already been created
    conda env list

    # If you don’t already have a Conda environment for R, create a new one (e.g., "r-env")
    conda create -n r-env -c conda-forge r

    # Activate the "r-env" Conda environment
    conda activate r-env

    # Install R within the Conda environment
    # Note: This step is only needed if R was not already installed.
    conda install -c conda-forge r

It’s recommended to create a separate Conda environment for Python modules used by reticulate (e.g., "r-reticulate") to avoid potential conflicts with the "r-env" environment.

# Bash commands for setting up the Conda environment "r-reticulate":

.. code-block:: bash

    # Create a new Conda environment specifically for reticulate (e.g. r-reticulate)
    conda create -n r-reticulate python

    # Verify that Python is installed in this environment
    python --version

    # If Python is not installed, you can install it with the following command:
    conda install python

    # Install the required Python modules (e.g., pandas, scipy, matplotlib) for the TRIAGE R package
    conda install pandas scipy matplotlib requests scikit-learn seaborn


In future R sessions (when "r-env" is active), you can activate the "r-reticulate" environment within R using `use_condaenv`:

.. code-block:: R

    library(reticulate)
    use_condaenv("r-reticulate", required = TRUE)

    # You can verify that the modules were installed successfully. For example:
    py_module_available("pandas")
    py_module_available("matplotlib")


**Option 2: Using a Python virtual environment**

To use a Python virtual environment, you’ll need to create and install required Python modules within R using the reticulate package. This setup is suitable for projects that need isolated Python environments without the broader management features provided by Conda.

Start by installing the required Python version and modules:

.. code-block:: R

    library(reticulate)
    # Note: any Python with version >=3.9 works. Here we install Python version 3.9.5 as an example. 
    reticulate::install_python(version = '3.9.5')
    reticulate::py_install("pandas", envname = "r-reticulate")
    reticulate::py_install("scipy", envname = "r-reticulate")
    reticulate::py_install("matplotlib", envname = "r-reticulate")
    reticulate::py_install("requests", envname = "r-reticulate")
    reticulate::py_install("scikit-learn", envname = "r-reticulate")
    reticulate::py_install("seaborn", envname = "r-reticulate")


In future R sessions, you can activate the “r-reticulate” environment to load these modules:

.. code-block:: R

    library(reticulate)
    use_virtualenv("r-reticulate", required = TRUE)


Installing TRIAGE R Package
---------------------------

1. Download the TRIAGE R package (e.g., TRIAGE_1.1.5.tar.gz) from the GitHub repository: https://github.com/Qiongyi/TRIAGE_R_Package

2. Install the TRIAGE R package from the source file:

.. code-block:: R

    install.packages("path/to/TRIAGE_1.1.5.tar.gz", repos = NULL, type = "source")
    library(TRIAGE)
