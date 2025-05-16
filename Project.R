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
