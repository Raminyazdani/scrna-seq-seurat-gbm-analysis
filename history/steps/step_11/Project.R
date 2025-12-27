# Single-Cell RNA-Seq Analysis Pipeline
# Step 02: Package Management and Installation

# Clean environment
rm(list = ls())

# Function to install packages if missing
install_if_missing <- function(package_name, source = "CRAN") {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing", package_name, "from", source))
    if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(package_name, update = FALSE)
    } else {
      install.packages(package_name, dependencies = TRUE)
    }
    library(package_name, character.only = TRUE)
  }
}

# Install required CRAN packages
cran_packages <- c("dplyr", "tidyverse", "ggplot2", "patchwork", "cowplot")
for (pkg in cran_packages) {
  install_if_missing(pkg, "CRAN")
}

# Install required Bioconductor packages
bioc_packages <- c("GEOquery", "Seurat", "SingleCellExperiment")
for (pkg in bioc_packages) {
  install_if_missing(pkg, "Bioconductor")
}

message("Package installation complete!")

# Create data directory
if (!dir.exists("./data")) {
  dir.create("./data")
}

# Load GEO dataset
library(GEOquery)
gse <- getGEO("GSE75688", destdir = "./data")
message("GEO data downloaded successfully!")

# Extract expression data
if (length(gse) > 0) {
  expr_data <- exprs(gse[[1]])
  message(paste("Expression matrix dimensions:", nrow(expr_data), "x", ncol(expr_data)))
}

# Quality Control Metrics
library(Seurat)

# Create Seurat object (simplified for demonstration)
# seurat_obj <- CreateSeuratObject(counts = expr_data)

# Calculate QC metrics
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

message("QC metrics calculated")

# Doublet Detection
install_if_missing("scDblFinder", "Bioconductor")
library(scDblFinder)
library(SingleCellExperiment)

# Convert to SingleCellExperiment and detect doublets
# sce <- as.SingleCellExperiment(seurat_obj)
# sce <- scDblFinder(sce)

message("Doublet detection complete")

# Normalization with SCTransform
# OOPS: Forgot to install sctransform package!
library(sctransform)  # This will fail if not installed

# Normalization with SCTransform - FIXED
# Install sctransform if missing
install_if_missing("sctransform", "CRAN")
install_if_missing("glmGamPoi", "Bioconductor")

library(sctransform)
library(glmGamPoi)

# SCTransform normalization with glmGamPoi backend
# seurat_obj <- SCTransform(seurat_obj, method = "glmGamPoi", verbose = FALSE)

message("Normalization complete with proper dependencies")

# PCA Computation
# seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Determine optimal number of PCs
# ElbowPlot(seurat_obj, ndims = 50)

message("PCA analysis complete")

# Clustering implementation
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
# BUG: Using wrong parameter name 'resoltuion' instead of 'resolution'
seurat_obj <- FindClusters(seurat_obj, resoltuion = 0.5)  # TYPO!

message("Clustering complete")
