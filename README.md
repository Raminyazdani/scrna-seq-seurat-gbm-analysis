# Single Cell RNA Sequencing Project

**Project Type:** University Project  
**Primary Stack:** R

## Description

This is a Single Cell RNA Sequencing analysis project using Seurat and other Bioconductor packages. The project analyzes single-cell RNA-seq data with focus on data preprocessing, quality control, dimensionality reduction, clustering, and differential expression analysis.

## Tech Stack

- R (version 4.x recommended)
- Seurat
- GEOquery
- SingleCellExperiment
- dplyr
- Bioconductor packages

## Folder Structure

```
scRNA-seq_project/
├── Project.R                     # Main project analysis script
├── test.R                        # Test/helper script
├── exprMatrix__GBM_2.RData       # Expression matrix data (GBM dataset)
├── data/                         # Additional data files
├── scRNA-seq_project.docx        # Project documentation
├── .Rhistory                     # R history file
└── README.md                     # This file
```

## Setup / Installation

The project includes an automatic package installation function. When you run the script, it will:

1. Check for required packages
2. Install missing packages from CRAN or Bioconductor automatically

Required packages include:
- GEOquery (for downloading data from GEO)
- Seurat (for single-cell analysis)
- dplyr (for data manipulation)
- SingleCellExperiment
- Other Bioconductor dependencies

Manual installation:
```r
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c("GEOquery", "Seurat", "SingleCellExperiment"))
install.packages("dplyr")
```

## How to Run

1. Navigate to the project directory and open R:
```bash
cd scRNA-seq_project
R
```

2. Source the main script:
```r
source("Project.R")
```

The script will:
- Install any missing packages automatically
- Load the expression matrix data
- Perform single-cell RNA-seq analysis
- Generate visualizations and results

Alternatively, run from command line:
```bash
cd scRNA-seq_project
Rscript Project.R
```

## Inputs/Outputs

**Inputs:**
- `exprMatrix__GBM_2.RData` - Expression matrix for GBM (Glioblastoma) dataset
- Additional data files in `data/` directory
- GEO datasets (downloaded automatically if needed)

**Outputs:**
- Quality control plots
- UMAP/t-SNE dimensionality reduction visualizations
- Cluster identification and marker genes
- Differential expression analysis results
- Cell type annotations

## Notes

- The project uses GBM (Glioblastoma) single-cell RNA-seq data
- Automatic package installation is included in the script
- All paths are relative to the project directory
- The script clears the workspace at the beginning: `rm(list = ls())`
- For reproducibility, consider setting a random seed

## Troubleshooting

- If package installation fails, ensure BiocManager is up to date: `BiocManager::install(version = "3.18")`
- For Seurat issues: `install.packages("Seurat")` with the latest version
- If memory issues occur: `options(future.globals.maxSize = 8000 * 1024^2)`
- For GEOquery download issues, check internet connection and GEO availability
- Ensure you have sufficient disk space for downloading datasets
