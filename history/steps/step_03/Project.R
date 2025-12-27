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
