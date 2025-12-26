# Project Identity: Single-Cell RNA-Seq Analysis of Glioblastoma

## Professional Project Identity

### Display Title
**Single-Cell RNA-Seq Analysis of Glioblastoma with Seurat**

### Repo Slug
`scrna-seq-seurat-gbm-analysis`

### Tagline
End-to-end single-cell RNA sequencing analysis pipeline for glioblastoma (GBM) tumor characterization using Seurat and Bioconductor

### GitHub Description
A comprehensive single-cell RNA-seq analysis pipeline for glioblastoma samples, featuring quality control, dimensionality reduction, clustering, cell type annotation, and differential expression analysis using Seurat and Bioconductor packages.

### Primary Stack
- R (4.x)
- Seurat
- Bioconductor (GEOquery, SingleCellExperiment, scDblFinder)
- tidyverse, ggplot2, patchwork

### Topics/Keywords
- single-cell-rna-seq
- bioinformatics
- seurat
- glioblastoma
- cancer-genomics
- bioconductor
- dimensionality-reduction
- cell-type-annotation
- differential-expression
- quality-control
- r-language

### Problem & Approach
**Problem:** Glioblastoma (GBM) is a highly heterogeneous brain tumor. Understanding cellular composition and gene expression patterns at single-cell resolution is crucial for identifying therapeutic targets and understanding tumor biology.

**Approach:** This pipeline processes single-cell RNA-seq data from GEO (Gene Expression Omnibus), performs comprehensive quality control including doublet detection, applies normalization and dimensionality reduction (PCA, UMAP, t-SNE), identifies cell clusters, annotates cell types using reference databases, and performs differential expression analysis to identify marker genes and enriched biological pathways.

### Inputs & Outputs Overview

**Inputs:**
- GEO dataset accession (GSE75688 - Breast Cancer dataset used for demonstration)
- Raw TPM expression matrices
- Sample metadata

**Outputs:**
- Quality control metrics and visualizations (violin plots, scatter plots, histograms)
- Dimensionality reduction plots (PCA elbow plots, UMAP, t-SNE)
- Cell cluster assignments and visualizations
- Cell type annotations
- Marker gene lists
- Differential expression results
- Gene ontology enrichment analysis results
- Heatmaps of top marker genes

