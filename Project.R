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

# Environment setup
rm(list = ls())

message("Package installation helper created")

# Install core packages
install_if_missing("Seurat", "CRAN")
install_if_missing("dplyr", "CRAN")  
install_if_missing("ggplot2", "CRAN")

library(Seurat)
library(dplyr)
library(ggplot2)

message("Core packages loaded successfully")
