rm = (list = ls())
install_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        message(paste("Attempting to install", pkg, "from Bioconductor..."))
        BiocManager::install(pkg)
        library(pkg, character.only = TRUE)
      })
    }
  }
}

packages_list = c(
  "GEOquery", "Seurat", "dplyr", "SingleCellExperiment",
  "ggplot2", "tidyverse", "Matrix"
)

install_if_missing(packages_list)
rm(install_if_missing)

for (pkg in packages_list) {
  library(pkg, character.only = TRUE)
}

# Set working directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Create directories
dir.create("./figures", showWarnings = FALSE)
dir.create("./data", showWarnings = FALSE)

# Define paths
gse_file_path <- "./data/GSE75688_first_element.rds"
expr_matrix_path <- "./data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75688/suppl/GSE75688%5FGEO%5Fprocessed%5FBreast%5FCancer%5Fraw%5FTPM%5Fmatrix.txt.gz"

# Download GEO data
if (!file.exists(gse_file_path)) {
  gse_list <- getGEO("GSE75688", GSEMatrix = TRUE)
  gse <- gse_list[[1]]
  saveRDS(gse, file = gse_file_path)
} else {
  gse <- readRDS(gse_file_path)
}

# Download expression matrix
if (file.exists(expr_matrix_path)) {
  exprMatrix <- read.delim(expr_matrix_path, header = T)
  message("File loaded successfully.")
} else {
  download.file(url, destfile = paste0(expr_matrix_path, ".gz"), mode = "wb")
  message("File downloaded successfully.")
  gunzip(paste0(expr_matrix_path, ".gz"), destname = expr_matrix_path, remove = TRUE)
  exprMatrix <- read.delim(expr_matrix_path, header = T)
  message("File loaded successfully.")
}

# Filter data
keep_columns <- grepl("BC07", colnames(exprMatrix)) | !grepl("BC", colnames(exprMatrix))
exprMatrix <- exprMatrix[, keep_columns]

# Extract metadata
metadata <- pData(phenoData(gse))
colnames(metadata)[colnames(metadata) == "patient id:ch1"] <- "patient_name"
colnames(metadata)[colnames(metadata) == "title"] <- "patient_id"
metadata <- metadata[grepl("BC07", metadata$patient_id), ]

message("Data loading complete!")
