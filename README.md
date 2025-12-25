# Single-Cell RNA-Seq Analysis of Glioblastoma with Seurat

**Primary Stack:** R | Seurat | Bioconductor

## Overview

A comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline for studying glioblastoma tumor heterogeneity. This pipeline processes raw expression data through quality control, normalization, dimensionality reduction, clustering, cell type annotation, and differential expression analysis to identify cell populations and their molecular signatures.

## What Problem Does This Solve?

Glioblastoma (GBM) is a highly heterogeneous brain tumor where understanding cellular composition at single-cell resolution is crucial for identifying therapeutic targets. This pipeline enables:

- Identification of distinct cell populations within tumor samples
- Characterization of cell type-specific gene expression patterns
- Discovery of marker genes for different cell types
- Analysis of differential expression between sample conditions
- Biological pathway enrichment to understand cellular functions

## Tech Stack

- **R** (version 4.x recommended)
- **Seurat** - Single-cell analysis toolkit
- **Bioconductor packages:**
  - GEOquery - Download data from Gene Expression Omnibus
  - SingleCellExperiment - Single-cell data container
  - scDblFinder - Doublet detection
  - SingleR - Automated cell type annotation
  - celldex - Reference datasets for annotation
- **Analysis packages:**
  - dplyr, tidyverse - Data manipulation
  - ggplot2, patchwork, cowplot - Visualization
  - clusterProfiler, enrichplot - Pathway enrichment
  - sctransform, glmGamPoi - Normalization

## Repository Structure

```
scrna-seq-seurat-gbm-analysis/
├── Project.R                  # Main analysis pipeline
├── exploratory_analysis.R     # Exploratory breast cancer analysis
├── data/                      # Data directory (created on first run)
│   └── GSE75688_*.txt        # Downloaded expression matrices
├── figures/                   # Generated plots (created on first run)
├── LARGE_FILES.md            # Documentation for external large files
├── .Rhistory                 # R session history
└── README.md                 # This file
```

**Note:** Large data files (>100MB) are documented in `LARGE_FILES.md` and stored externally.

## Setup / Installation

This pipeline includes automatic package installation. When you run the script, it will:

1. Check for required packages
2. Install missing packages from CRAN or Bioconductor automatically

### Prerequisites

- R (version 4.0 or higher recommended)
- Internet connection for downloading packages and data

### Manual Installation (Optional)

If you prefer to install packages manually before running:

```r
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
  "GEOquery", "Seurat", "SingleCellExperiment", 
  "scDblFinder", "SingleR", "celldex",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot"
))

# Install CRAN packages
install.packages(c(
  "dplyr", "tidyverse", "ggplot2", "patchwork", 
  "cowplot", "viridis", "knitr", "kableExtra"
))
```

## How to Run

### Main Analysis Pipeline (Project.R)

**From within R/RStudio:**

```r
# Set working directory to repository root
setwd("/path/to/scrna-seq-seurat-gbm-analysis")

# Source the main script
source("Project.R")
```

**From command line:**

```bash
# Navigate to repository root
cd scrna-seq-seurat-gbm-analysis

# Run with Rscript
Rscript Project.R
```

**Note:** The main script uses `rstudioapi` for automatic working directory detection when run from RStudio. When running from command line with `Rscript`, ensure you're in the repository root directory.

### What the Pipeline Does

The script will automatically:
1. Install any missing R packages
2. Create `data/` and `figures/` directories if they don't exist
3. Download GEO dataset GSE75688 (Breast Cancer single-cell data)
4. Process and filter expression matrices
5. Perform quality control and doublet detection
6. Run normalization and dimensionality reduction
7. Identify cell clusters and annotate cell types
8. Generate visualizations and save to `figures/`
9. Perform differential expression and pathway enrichment analysis

### Runtime

- First run: ~30-60 minutes (includes package installation and data download)
- Subsequent runs: ~15-30 minutes (data already cached in `data/` directory)

## Data Sources and Inputs

### Primary Dataset
- **GEO Accession:** GSE75688
- **Description:** Breast Cancer single-cell RNA-seq TPM expression data
- **Source:** Automatically downloaded from NCBI GEO on first run
- **Samples:** Focused on BC07 patient samples (primary and lymph node metastasis)

### Data Files
The pipeline automatically downloads and caches:
- `GSE75688_first_element.rds` - GEO metadata object
- `GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt` - Expression matrix

All downloaded files are stored in the `data/` directory (created automatically).

### Large External Files
Some optional large files are documented in `LARGE_FILES.md` and available on Google Drive:
- `exprMatrix__GBM_2.RData` (215 MB) - GBM expression matrix
- `data/GSE89567_IDH_A_processed_data.txt` (470 MB) - Alternative dataset

**These are not required for the main pipeline.**

## Outputs

The pipeline generates multiple outputs organized by analysis stage:

### Quality Control
- Violin plots of QC metrics (nCount_RNA, nFeature_RNA, percent.mt, etc.)
- Feature scatter plots showing relationships between metrics
- Distribution plots (histograms, density, boxplots) for all QC metrics
- Doublet detection results

### Dimensionality Reduction
- PCA elbow plot for optimal component selection
- PCA dimension plots colored by patient
- UMAP plots (by cluster and by patient)
- t-SNE plots (by cluster and by patient)

### Clustering & Cell Type Annotation
- Cluster assignments stored in Seurat object metadata
- Cell type annotations using SingleR with Human Primary Cell Atlas
- Heatmaps of top marker genes per cluster
- Bar plots showing cell type composition per sample

### Differential Expression
- Marker genes for each cluster
- Differential expression between BC07 primary and metastatic samples
- Top differentially expressed genes with fold changes

### Pathway Enrichment
- Gene Ontology (GO) enrichment for each cluster
- GSEA results for marker genes
- Dot plots visualizing enriched pathways
- Results from multiple databases (KEGG, Reactome, GO, WikiPathway, etc.)

### Storage
- Plots stored in `plot_dict` list object (in R session)
- Tables stored in `table_dict` list object (in R session)
- Intermediate RDS files cached in `data/` for faster re-runs

## Analysis Methods

### Pipeline Stages

1. **Data Acquisition**
   - Download from GEO using GEOquery
   - Load and parse expression matrices
   - Extract and format metadata

2. **Quality Control**
   - Calculate QC metrics (gene counts, UMI counts, mitochondrial %)
   - Identify doublets using scDblFinder
   - Remove outliers based on IQR and MAD thresholds
   - Filter low-quality cells

3. **Normalization & Feature Selection**
   - SCTransform normalization with glmGamPoi
   - Select 2,000 highly variable features
   - Scale and center data

4. **Dimensionality Reduction**
   - PCA with automatic optimal component selection
   - UMAP for 2D visualization
   - t-SNE for alternative visualization

5. **Clustering**
   - Graph-based clustering (Louvain algorithm)
   - Optimal resolution selection via modularity optimization
   - Cluster visualization in reduced dimensions

6. **Cell Type Annotation**
   - Automated annotation using SingleR
   - Reference: Human Primary Cell Atlas
   - Manual validation via marker genes

7. **Differential Expression**
   - FindAllMarkers for cluster-specific genes
   - FindMarkers for pairwise comparisons
   - Statistical testing with MAST

8. **Pathway Enrichment**
   - Gene Ontology enrichment analysis
   - GSEA for ranked gene lists
   - Multiple pathway databases queried

## Reproducibility Notes

- **Random Seed:** Not currently set - results may vary slightly between runs due to stochastic algorithms (PCA, clustering)
- **Software Versions:** Pipeline tested with R 4.x and Seurat 4.x/5.x
- **Data Caching:** Downloaded data is cached in `data/` directory to ensure consistency across runs
- **Deterministic Steps:** Data loading, QC filtering, and normalization are deterministic
- **Stochastic Steps:** PCA initialization, clustering, UMAP/t-SNE embeddings may vary

### For Full Reproducibility

Add this line at the start of your R session:
```r
set.seed(42)
```

## Troubleshooting

### Package Installation Issues

**Problem:** BiocManager installation fails
```r
# Update BiocManager to latest version
install.packages("BiocManager")
BiocManager::install(version = "3.18")
```

**Problem:** Seurat installation fails
```r
# Install from CRAN with dependencies
install.packages("Seurat", dependencies = TRUE)
```

### Memory Issues

**Problem:** R runs out of memory during analysis
```r
# Increase memory limit (adjust value as needed)
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB
```

### Data Download Issues

**Problem:** GEOquery cannot download dataset
- Check internet connection
- Verify GEO database is accessible: https://www.ncbi.nlm.nih.gov/geo/
- Try downloading manually and placing in `data/` directory

**Problem:** Insufficient disk space
- Ensure ~500MB free space for data downloads
- Check `data/` directory size

### Path and Working Directory Issues

**Problem:** Script cannot find files when run from command line
- Ensure you're running from repository root directory
- Use absolute paths if needed
- If using RStudio, the script auto-detects working directory

**Problem:** `rstudioapi` error when using Rscript
- This is expected - the script falls back to current directory
- Ensure you're in repository root: `cd scrna-seq-seurat-gbm-analysis`

### Runtime Issues

**Problem:** Script takes very long on first run
- First run includes package installation (~10-15 min)
- Data download from GEO (~5-10 min)
- Actual analysis (~15-30 min)
- Subsequent runs are faster due to caching

## License

This analysis pipeline is available for academic and research use. Please cite the appropriate tools and datasets:

- **Seurat:** Hao et al., Cell 2021
- **scDblFinder:** Germain et al., F1000Research 2021  
- **SingleR:** Aran et al., Nature Immunology 2019
- **Dataset GSE75688:** Chung et al., Nature Communications 2017

## Citation

If you use this pipeline in your research, please cite the original dataset and key tools as listed above.
