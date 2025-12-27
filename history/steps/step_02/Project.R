# Single-Cell RNA-Seq Analysis Pipeline
# Step 02: Core Package Installation

# Install core packages (Seurat, dplyr, ggplot2)
packages_list <- c("Seurat", "dplyr", "ggplot2")

install_if_missing <- function(pkg, source = "CRAN") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (source == "CRAN") {
      install.packages(pkg, dependencies = TRUE)
    } else if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg)
    }
  }
}

for (pkg in packages_list) {
  install_if_missing(pkg, "CRAN")
}

message("Core packages installed successfully")
