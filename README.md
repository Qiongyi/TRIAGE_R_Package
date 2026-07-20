# TRIAGE R Package

> [!IMPORTANT]
> **TRIAGE Toolkit v2.0.0 (also referred to as TRIAGE R package v2) is now available.**  
> The updated toolkit is described in our 2026 *Current Protocols* paper:  
> [TRIAGE Toolkit: Streamlined Discovery of Regulatory Genes and Elements](https://doi.org/10.1002/cpz1.70413)


## TRIAGE Toolkit v2.0.0

TRIAGE is a streamlined computational toolkit for discovering and prioritizing regulatory genes and genomic elements involved in cell identity, development, and disease. The toolkit supports regulatory analysis using diverse input types, including gene expression matrices, gene lists, and genomic loci, and can be applied to both bulk and single-cell RNA-seq data.

The TRIAGE toolkit includes four main components:

- **TRIAGEgene**: prioritizes regulatory genes from gene expression data
- **TRIAGEcluster**: identifies biologically distinct cell populations from single-cell RNA-seq data
- **TRIAGEparser**: groups genes into functionally related regulatory modules and supports downstream functional interpretation
- **TRIAGEccs**: prioritizes regulatory genomic elements at single-base resolution


## Access and Licensing

TRIAGE Toolkit v2.0.0 is available through the [UniQuest Online Store](https://uniquest.store/product/triage2) under one of the following licences:

- **Academic Research and Teaching Licence**: free for academic research and teaching
- **General Use Licence**: for other uses, including commercial research

## Documentation

The TRIAGE Toolkit v2.0.0 documentation is available at:

[https://tinyurl.com/triage2doc](https://tinyurl.com/triage2doc)



# TRIAGE R Package v1.x

This GitHub repository provides the TRIAGE R package v1.x releases and associated source materials. With the exception of TRIAGEccs, the core functions available in TRIAGE Toolkit v2.0.0 are also available in the v1.x package and provide equivalent functionality.

Users who do not require the single-base-resolution regulatory element analysis provided by TRIAGEccs may continue to use the TRIAGE R package v1.x from this GitHub repository.

Important bug fixes and essential updates affecting the shared functions will also be incorporated into future v1.x releases where applicable.


# Overview
The TRIAGE R package is a powerful, integrated toolkit developed by our team to decipher the regulatory complexities of transcriptome dynamics and cell identity in development and disease. 

TRIAGE R package is particularly useful for researchers focusing on:
- Transcriptome analysis to uncover regulatory gene networks, applicable to both bulk RNA-seq and single-cell RNA-seq
- Discovery, characterization, and prioritization of regulatory genes
- Exploring key regulatory genes in developmental and disease contexts

# Features
- Integrated toolkit: We consolidate all TRIAGE methods into a single, cohesive package, simplifying usage for researchers with even basic R knowledge.
- User-friendly functions for regulatory gene analysis and visualization: The package provides a suite of streamlined functions that allow for seamless integration of regulatory gene analysis into standard workflows.
- Multi-species support: Extends its utility across a variety of biological models, making it versatile for different research contexts.
- Enhanced performance: Features significant improvements in computational efficiency and compatibility.

# Download the Latest Release
You can download the latest release of TRIAGE [here](https://github.com/Qiongyi/TRIAGE_R_Package/releases/latest/download/TRIAGE_1.1.7.tar.gz).

# Documentation
The complete TRIAGE R package manual is available via [readthedocs](https://triage-r-package.readthedocs.io/en/latest/index.html).

# Citing TRIAGE R Package
If you use or discuss the TRIAGE R package in your research, please cite our paper: 

**Qiongyi Zhao**, et al. "TRIAGE: an R package for regulatory gene analysis"
*Briefings in Bioinformatics*, Volume 26, Issue 1, January 2025, bbaf004.
[https://doi.org/10.1093/bib/bbaf004](https://doi.org/10.1093/bib/bbaf004)
