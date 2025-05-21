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
