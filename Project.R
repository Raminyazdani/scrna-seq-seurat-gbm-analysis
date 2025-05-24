# Helper function for package installation
install_if_missing <- function(package_name, source = "CRAN") {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package_name)
    } else {
      install.packages(package_name)
    }
  }
}

rm(list = ls())

# Core packages
install_if_missing("Seurat", "CRAN")
install_if_missing("dplyr", "CRAN")
install_if_missing("ggplot2", "CRAN")
install_if_missing("GEOquery", "Bioconductor")

library(Seurat)
library(dplyr)
library(ggplot2)
library(GEOquery)

# Download and load data
gse <- getGEO("GSE75688", GSEMatrix = FALSE)

# FIXED: Check if file exists and create directory if needed
if (!dir.exists("data")) {
  dir.create("data")
}

expr_file <- "data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
if (!file.exists(expr_file)) {
  stop("Expression matrix file not found. Please download from GEO.")
}

expr_data <- read.table(expr_file, header = TRUE, row.names = 1)
message("Expression data loaded successfully")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = expr_data, 
                                  project = "scRNA_analysis")
message("Seurat object created")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
message("QC metrics calculated: percent.mt, percent.ribo, nCount_RNA, nFeature_RNA")

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"))

# IQR-based filtering
lower_nFeature <- quantile(seurat_obj$nFeature_RNA, 0.25) - 1.5 * IQR(seurat_obj$nFeature_RNA)
upper_nFeature <- quantile(seurat_obj$nFeature_RNA, 0.75) + 1.5 * IQR(seurat_obj$nFeature_RNA)

seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > lower_nFeature & 
                             nFeature_RNA < upper_nFeature)
message("IQR filtering applied")

# MAD-based filtering for mitochondrial content
median_mt <- median(seurat_obj$percent.mt)
mad_mt <- mad(seurat_obj$percent.mt)
upper_mt <- median_mt + 3 * mad_mt

seurat_obj <- subset(seurat_obj, subset = percent.mt < upper_mt)
message("MAD filtering applied for mitochondrial content")

# Create scatter plots showing relationships
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
message("QC filtering complete")
