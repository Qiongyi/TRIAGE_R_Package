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

Python Environment Setup
------------------------

Install Python modules within R using the reticulate package. It's recommended to manage Python dependencies in a separate environment.

.. code-block:: R

    library(reticulate)
    # Note: any Python with version >=3.9 works. Here we install Python version 3.9.5 as an example. 
    reticulate::install_python(version = '3.9.5')
    reticulate::py_install("pandas", envname = "r-reticulate")
    reticulate::py_install("numpy", envname = "r-reticulate")
    reticulate::py_install("scipy", envname = "r-reticulate")
    reticulate::py_install("matplotlib", envname = "r-reticulate")
    reticulate::py_install("requests", envname = "r-reticulate")
    reticulate::py_install("scikit-learn", envname = "r-reticulate")
    reticulate::py_install("seaborn", envname = "r-reticulate")

Knowledge Base
--------------

When using the reticulate package in R to install Python modules, it is important to understand how reticulate::py_install manages the installation of these modules. Below are key points to consider:

- **Intelligent Installation**: The reticulate::py_install function is designed to intelligently manage Python modules within your specified Python environment. It checks whether the required modules are already installed in the targeted environment.

- **No Redundant Installation**: If a module is already installed in the given environment, py_install will not attempt to reinstall it. This avoids redundant installation and saves time, especially when setting up the environment for the first time or in a new session.

- **Specified Environment**: You can specify the environment where you want the modules to be installed. This is done using the envname parameter in the py_install function. For instance, reticulate::py_install("pandas", envname = "r-reticulate") will install the pandas module in the r-reticulate environment.


In future R sessions, activate the r-reticulate environment:

.. code-block:: R

    library(reticulate)
    use_virtualenv("r-reticulate", required = TRUE)


Installing TRIAGE R Package
---------------------------

1. Download the TRIAGE R package (e.g., TRIAGE_1.1.1.tar.gz) from the GitHub repository: https://github.com/Qiongyi/TRIAGE_R_Package

2. Install the TRIAGE R package from the source file:

.. code-block:: R

    install.packages("path/to/TRIAGE_1.1.1.tar.gz", repos = NULL, type = "source")
    library(TRIAGE)
