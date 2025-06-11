# ============================================================================
# Exploratory Single-Cell RNA-Seq Analysis - Breast Cancer Dataset
# ============================================================================
#
# Description:
#   Exploratory analysis script for breast cancer single-cell RNA-seq data
#   from GEO dataset GSE75688. This script demonstrates the analysis pipeline
#   with a focus on BC07 patient samples.
#
# Usage:
#   From RStudio: source("exploratory_analysis.R")
#   From command line: Rscript exploratory_analysis.R (ensure you're in repo root)
#
# Inputs:
#   - GEO dataset GSE75688 (downloaded automatically)
#   - Breast Cancer single-cell TPM expression matrix
#
# Outputs:
#   - QC plots and statistics
#   - Dimensionality reduction visualizations (PCA, UMAP, t-SNE)
#   - Cell cluster assignments
#   - All outputs stored in plot_dict and table_dict objects
#
# Note:
#   This is a simplified version compared to Project.R, focusing on core
#   analysis steps without extensive enrichment analysis.
#
# Repository: https://github.com/Raminyazdani/scrna-seq-seurat-gbm-analysis
# ============================================================================

rm = (list = ls())
install_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      # First try installing from CRAN
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        # If CRAN installation fails, try Bioconductor
        message(paste("Attempting to install", pkg, "from Bioconductor..."))
        BiocManager::install(pkg)
        library(pkg, character.only = TRUE)
      })
    }
  }
}



packages_list = c(
  "GEOquery",
  "Seurat",
  "dplyr",
  "SingleCellExperiment",
  "Seurat",
  "ggplot2",
  "sctransform",
  "BiocParallel",
  "scDblFinder",
  "cowplot",
  "RColorBrewer",
  "viridis",
  "knitr",
  "viridisLite",
  "patchwork",
  "tidyverse",
  "dplyr",
  "Matrix",
  "kableExtra",
  "harmony",
  "glmGamPoi"
)

install_if_missing(packages_list)
rm(install_if_missing)
# Load the packages

for (pkg in packages_list) {
  library(pkg, character.only = TRUE)
}
rm(list = ls())
plot_dict <- list()
table_dict <- list()

# Set working directory
# Note: When running from RStudio, this auto-detects the script location
# When running from command line with Rscript, ensure you're in the repo root
if (requireNamespace("rstudioapi", quietly = TRUE) && 
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # When running from command line, assume current directory is correct
  # User should run from repository root
  cat("Running from command line. Current directory:", getwd(), "\n")
  cat("Ensure you're in the repository root directory.\n")
}
getwd()

# Create directories for figures and data if they don't exist
dir.create("./figures", showWarnings = FALSE)
dir.create("./data", showWarnings = FALSE)

#############################################################
# Initialize pipeline inputs
#############################################################

# Define the path to save the GSE data
gse_file_path <- "./data/GSE75688_first_element.rds"

# Define the path and URL for the expression matrix
expr_matrix_path <- "./data/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75688/suppl/GSE75688%5FGEO%5Fprocessed%5FBreast%5FCancer%5Fraw%5FTPM%5Fmatrix.txt.gz"

#############################################################
# End of pipeline inputs initialization
#############################################################

# Check if the GSE file exists in the data directory
if (!file.exists(gse_file_path)) {
  # If the file doesn't exist, download the GSE dataset
  gse_list <- getGEO("GSE75688", GSEMatrix = TRUE)
  
  # Extract the first element from the list
  gse <- gse_list[[1]]
  
  # Save the first element of the GSE object as an RDS file
  saveRDS(gse, file = gse_file_path)
} else {
  # If the file exists, load the first element of the GSE object from the RDS file
  gse <- readRDS(gse_file_path)
}


# Check if the expression matrix file exists
if (file.exists(expr_matrix_path)) {
  # Load the file if it exists
  exprMatrix <- read.delim(expr_matrix_path, header = T)
  
  message("File loaded successfully.")
} else {
  # Download the file if it doesn't exist
  download.file(url,
                destfile = paste0(expr_matrix_path, ".gz"),
                mode = "wb")
  message("File downloaded successfully.")
  
  # Extract the file
  gunzip(paste0(expr_matrix_path, ".gz"),
         destname = expr_matrix_path,
         remove = TRUE)
  message("File extracted successfully.")
  
  # Load the file
  exprMatrix <- read.delim(expr_matrix_path, header = T)
  message("File loaded successfully.")
}
rm(url, gse_file_path)

# Filter columns related to BC07 and non-BC patients
keep_columns <- grepl("BC07", colnames(exprMatrix)) |
  !grepl("BC", colnames(exprMatrix))
exprMatrix <- exprMatrix[, keep_columns]

# Extract metadata from the phenoData of the GSE object
metadata <- pData(phenoData(gse))

# Rename columns for easier reference
colnames(metadata)[colnames(metadata) == "patient id:ch1"] <- "patient_name"
colnames(metadata)[colnames(metadata) == "title"] <- "patient_id"

# Filter metadata to only include data related to patient BC07
metadata <- metadata[grepl("BC07", metadata$patient_id), ]
rm(keep_columns, expr_matrix_path)  # Remove unused variables

# Create a backup of the original expression matrix
exprMatrix.org = exprMatrix

##########################################
# Select specific columns from the expression matrix
gene_data = exprMatrix %>% select(gene_id, gene_name, gene_type)

# Set row names for the expression matrix and return to the original expression matrix
exprMatrix = exprMatrix.org
# Step 1: Identify the rows with duplicate gene_name
duplicates <- exprMatrix %>%
  filter(gene_name %in% gene_name[duplicated(gene_name)])

# Step 2: Calculate the mean of non-zero numeric columns for the duplicates
merged_duplicates <- duplicates %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), ~ {
    mean_value <- mean(.x[.x != 0], na.rm = TRUE)
    ifelse(is.nan(mean_value), 0, mean_value)
  }), across(where( ~ !is.numeric(.)), ~ .x[1]))

# Step 3: Create a dataframe with rows that do not have duplicated gene_name
non_duplicates <- exprMatrix %>%
  filter(!gene_name %in% gene_name[duplicated(gene_name)])

# Step 4: Combine the merged duplicates with the non-duplicated dataframe
exprMatrix <- bind_rows(merged_duplicates, non_duplicates)
exprMatrix <- as.data.frame(exprMatrix)

rm(duplicates, merged_duplicates, non_duplicates)
rownames(exprMatrix) = exprMatrix$gene_name
exprMatrix = exprMatrix %>% select(-gene_name, -gene_id, -gene_type)

##########################################

# Create a SingleCellExperiment object with counts and logcounts matrices
sce <- SingleCellExperiment(
  list(
    counts = as.matrix(exprMatrix),
    logcounts = as.matrix(exprMatrix)
  ),
  colData = DataFrame(
    patient.id = metadata$patient_id,
    patient.name = metadata$patient_name
  ),
  metadata = list(study = "gbm_scrnaseq_analysis")
)

# Add gene symbols to rowData
rowData(sce)$Symbol = (rownames(sce$gene_name))

# Convert counts matrix to dgCMatrix format for efficient storage and computation
all.counts = as.matrix(counts(sce))
cell_metadata = as.data.frame(colData(sce))
rm(sce)
rownames(cell_metadata) <- 1:nrow(cell_metadata)
all.counts <- as(all.counts, "dgCMatrix")

# Create a Seurat object from the counts matrix and cell metadata
gbm = CreateSeuratObject(counts = all.counts, meta.data = cell_metadata)
rm(all.counts, cell_metadata)


## Calculate Quality Control (QC) metrics:
gbm[["percent.mt"]] <-  PercentageFeatureSet(gbm, pattern = "^MT-")  # Mitochondrial gene percentage

gbm[["percent.ERCC"]] <-  PercentageFeatureSet(gbm, pattern = "^ERCC")  # ERCC spike-in percentage
gbm[["percent.ribo"]] <-  PercentageFeatureSet(gbm, pattern = "^RPS |^RPL")  # Ribosomal gene percentage
gbm[["percent.hb"]] <-  PercentageFeatureSet(gbm, pattern = "^HB[^(P)]")  # Hemoglobin gene percentage

table_dict$QC_metrics = gbm@meta.data %>% select(1:ncol(gbm@meta.data)) %>% rename("nCount_RNA" =
                                                                                     "genes detected per cell") %>% rename("nFeature_RNA" = "Total UMI counts per cell")

# Extract metadata from the Seurat object
numeric_metadata <- gbm@meta.data[, sapply(gbm@meta.data, function(x)
  is.numeric(x) || is.integer(x))]

# Plot Violin plots for numeric metadata features
violin_plots <- VlnPlot(
  gbm,
  features = colnames(numeric_metadata),
  ncol = 3,
  group.by = "patient.name",
  pt.size = 0.001,
  combine = FALSE
)

# If you want to combine them into one plot, use CombinePlots (if multiple plots) or leave it as is
combined_violin_plot <- wrap_plots(violin_plots)

# Create a list and assign the combined plot to the key 'violin_plot'
plot_dict$violin_plot <- combined_violin_plot
rm(violin_plots, combined_violin_plot)
plot_dict$violin_plot

# Create scatter plots to explore relationships between different QC metrics
f1 <- FeatureScatter(gbm,
                     feature1 = "nCount_RNA",
                     feature2 = "nFeature_RNA",
                     group.by = "patient.name")
f2 <- FeatureScatter(gbm,
                     feature1 = "nCount_RNA",
                     feature2 = "percent.mt",
                     group.by = "patient.name")
f3 <- FeatureScatter(gbm,
                     feature1 = "percent.mt",
                     feature2 = "nFeature_RNA",
                     group.by = "patient.name")

plot_dict$FeaturScatterPlots <- wrap_plots(f1, f2, f3)
wrap_plots(f1, f2, f3)
rm(f1, f2, f3)
##########################################

# Function to calculate summary statistics
calculate_stats <- function(x) {
  return(
    c(
      min = min(x, na.rm = TRUE),
      q1 = quantile(x, 0.25, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      mean = mean(x, na.rm = TRUE),
      q3 = quantile(x, 0.75, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      iqr = IQR(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  )
}

# Apply the calculate_stats function to each numeric column in metadata
stats_summary <- t(sapply(numeric_metadata, calculate_stats))

# Convert stats summary to a data frame
stats_summary_df <- as.data.frame(stats_summary)
rm(stats_summary)

# Display the statistics summary as a formatted HTML table
kable(
  stats_summary_df,
  format = "html",
  digits = 2,
  row.names = TRUE
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = TRUE
  ) %>%
  scroll_box(width = "100%", height = "400px")
table_dict$initial_stats_summary_df <- stats_summary_df
rm(stats_summary_df)


# Initialize a list to store plots
plot_list <- list()

# Loop through each numeric metadata column to create histograms, density plots, and box plots
for (col in colnames(numeric_metadata)) {
  # Histogram
  p_hist <- ggplot(gbm@meta.data, aes_string(x = col)) +
    geom_histogram(bins = 30,
                   color = "black",
                   fill = "lightblue") +
    labs(title = paste("Histogram of", col),
         x = col,
         y = "Frequency") +
    theme_minimal(base_size = 10) +  # Smaller base font size
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      # Adjust title font size
      axis.title = element_text(size = 10),
      # Axis title size
      axis.text = element_text(size = 8)
    )  # Axis text size
  
  # Density Plot
  p_density <- ggplot(gbm@meta.data, aes_string(x = col)) +
    geom_density(fill = "lightblue", alpha = 0.5) +
    labs(title = paste("Density Plot of", col),
         x = col,
         y = "Density") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  # Box Plot
  p_box <- ggplot(gbm@meta.data, aes_string(x = "1", y = col)) +
    geom_boxplot(fill = "lightblue") +
    labs(title = paste("Box Plot of", col),
         x = "",
         y = col) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
  
  # Combine the three plots into one row
  combined_plot <- p_hist + p_density + p_box + plot_layout(ncol = 3)  # Arrange in a single row
  
  # Add the combined plot to the list
  plot_list[[col]] <- combined_plot
  rm(combined_plot)
}
rm(col)
# Combine all plots into one big picture
final_plot <- wrap_plots(plot_list, ncol = 1) +  # Combine into a single column layout
  plot_annotation(title = "Combined Plots for All Metrics",
                  theme = theme(plot.title = element_text(size = 16)))  # Add a title to the combined plot

plot_dict$Metrics_plot <- final_plot

# Display the final combined plot
print(final_plot)

rm(final_plot,
   p_box,
   p_density,
   p_hist,
   plot_list,
   numeric_metadata)

# Convert Seurat object to SingleCellExperiment object for further analysis
gbm.sce = as.SingleCellExperiment(gbm)
gbm.sce = scDblFinder(gbm.sce, samples = "orig.ident", clusters = F)  # Detect doublets

# Log-transform the counts data
logcounts(gbm.sce) <- log1p(counts(gbm.sce))  # log1p(x) is log(1 + x)

# Display the number of doublets i66622dentified
table_dict$scDblFinder.class <- table(gbm.sce$scDblFinder.class)
table_dict$scDblFinder.class

# Convert SingleCellExperiment object back to Seurat object
gbm.seurat <- as.Seurat(gbm.sce, counts = "counts", data = "logcounts")

# Apply Quality Control: Doublet removal
gbm <- subset(gbm.seurat, subset = scDblFinder.class == "singlet")
rm(gbm.seurat, gbm.sce)
# Identify numeric columns in the metadata
numeric_metadata <- gbm@meta.data[, sapply(gbm@meta.data, function(x)
  is.numeric(x) || is.integer(x))]

stats_summary <- t(sapply(numeric_metadata, calculate_stats))

# Convert stats_summary to a data frame if it isn't already
stats_summary_df <- as.data.frame(stats_summary)

kable(
  stats_summary_df,
  format = "html",
  digits = 2,
  row.names = TRUE
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = TRUE
  ) %>%
  scroll_box(width = "100%", height = "400px")
table_dict$After_scDblFinder_stats_summary_df <- stats_summary_df
rm(stats_summary, numeric_metadata)


# Define thresholds using the interquartile range (IQR)
thresh_iqr = 5
nCount_RNA_upper = stats_summary_df["nCount_RNA", "q3.75%"] + thresh_iqr *
  stats_summary_df["nCount_RNA", "iqr"]
nCount_RNA_lower = stats_summary_df["nCount_RNA", "q1.25%"] - thresh_iqr *
  stats_summary_df["nCount_RNA", "iqr"]
nFeature_RNA_upper = stats_summary_df["nFeature_RNA", "q3.75%"] + thresh_iqr *
  stats_summary_df["nFeature_RNA", "iqr"]
nFeature_RNA_lower = stats_summary_df["nFeature_RNA", "q1.25%"] - thresh_iqr *
  stats_summary_df["nFeature_RNA", "iqr"]
rm(stats_summary_df)


# Subset the Seurat object based on the calculated thresholds
gbm = subset(
  gbm,
  subset = nCount_RNA >= nCount_RNA_lower &
    nCount_RNA <= nCount_RNA_upper
  &
    nFeature_RNA >= nFeature_RNA_lower &
    nFeature_RNA <= nFeature_RNA_upper
)
rm(
  nCount_RNA_upper,
  nCount_RNA_lower,
  nFeature_RNA_upper,
  nFeature_RNA_lower,
  thresh_iqr
)

# Function to identify outliers based on median absolute deviation (MAD)
is_outlier <- function(scdata, metric, nmads) {
  M <- scdata@meta.data[, metric]
  outlier <- (M < median(M) - nmads * mad(M)) |
    (median(M) + nmads * mad(M) < M)
  return(outlier)
}

# Identify and mark outliers
gbm@meta.data["outlier"] = (is_outlier(gbm, "nCount_RNA", 5) |
                              is_outlier(gbm, "nFeature_RNA", 5))
gbm@meta.data["ribo_outlier"] = is_outlier(gbm, "percent.ribo", 3)
gbm@meta.data["hb_outlier"] = is_outlier(gbm, "percent.hb", 3)
gbm@meta.data["ERCC_outlier"] = is_outlier(gbm, "percent.ERCC", 3)
gbm@meta.data["mt_outlier"] = is_outlier(gbm, "percent.mt", 3)

gbm@meta.data["all_outliers"] <- gbm@meta.data["outlier"] |
  gbm@meta.data["ribo_outlier"] |
  gbm@meta.data["hb_outlier"] |
  gbm@meta.data["ERCC_outlier"] | gbm@meta.data["mt_outlier"]

# Display the number of outliers
table(gbm@meta.data["all_outliers"])

# Subset the Seurat object to exclude outliers
gbm <- subset(gbm, subset = outlier == FALSE)

#################################################

# Normalize, find variable features, and scale data
gbm <- SCTransform(
  gbm,
  method = "glmGamPoi",
  variable.features.n = 2000,
  verbose = FALSE
)

# finding optimal pca factor
optimal_pcs = 50
opts = vector()
opts = append(opts, 50)

while (TRUE) {
  gbm_test = RunPCA(gbm,
                    features = VariableFeatures(gbm),
                    npcs = optimal_pcs,
                    verbose = FALSE)
  dim(gbm_test[["pca"]])[2]
  
  stdev <- gbm_test[["pca"]]@stdev
  cum_variance <- cumsum(stdev ^ 2) / sum(stdev ^ 2)
  optimal_pcs <- which(cum_variance > 0.9)[1]
  opts <- append(opts, optimal_pcs)
  if (length(opts) > 3) {
    if (opts[length(opts) - 1] - opts[length(opts) - 2] > opts[length(opts)
                                                               ] - opts[length(opts)]-1) {
      break
    }
  }
}

npcs_n = ceiling(sum(opts)/length(opts))
rm(e,cum_variance,optimal_pcs,opts,stdev,gbm_test)


# running with optimal pca
gbm = RunPCA(gbm, features = VariableFeatures(gbm), npcs = npcs_n)

e = ElbowPlot(gbm, ndims = npcs_n)

plot_dict$elbowPCA <- e

rm(e)


all_dims = list()
for (i in 1:npcs_n) {
  for (j in i:npcs_n) {
    if (j <= i) {
      next
    }
    d = DimPlot(
      gbm ,
      reduction = "pca",
      dims = c(i, j),
      group.by = "patient.name",
      pt.size = 0.5,
    )+ theme_cowplot(font_size = 4) + guides(color = guide_legend(override.aes = list(size=0.5), ncol=1) )
    all_dims[[paste(i, "-", j)]] <- d
  }
}

plot_dict$all_dims <-all_dims
rm(i, j, d, npcs_n,all_dims)

# all_dims has all the dimensions dimplots  , visualizations

wrap_plots(plot_dict$all_dims$`1 - 2`,plot_dict$all_dims$`2 - 3`)
wrap_plots(plot_dict$all_dims)


# Apply Harmony to correct for batch effects based on patient name
#gbm <- RunHarmony(gbm, "patient.name")
#DimPlot(gbm,
#        reduction = "harmony",
#        dims = c(1, 2),
#        group.by = "patient.name")
# NOOOOOOOO NEEEEEEEEED

# Find neighbors and clusters
gbm <- FindNeighbors(gbm, dims = 1:20, reduction = "pca")
gbm <- FindClusters(gbm,
                    resolution = 0.5,
                    method = "igraph",
                    n.iter = 20)
gbm <- FindClusters(gbm,
                    resolution = 0.01,
                    method = "igraph",
                    n.iter = 20)
gbm <- FindClusters(gbm,
                    resolution = 1,
                    method = "igraph",
                    n.iter = 20)

# Run UMAP and t-SNE for visualization
gbm <- RunUMAP(gbm, reduction = "pca", dims = 1:20)
gbm <- RunTSNE(gbm,
               reduction = "pca",
               dims = 1:20,
               perplexity = 30)

# Plot UMAP and t-SNE results
p1_umap <- DimPlot(object = gbm, reduction = "umap")
p2_umap <- DimPlot(object = gbm,
                   reduction = "umap",
                   group.by = "patient.name")

# Combine UMAP plots into a grid
plot_grid(p1_umap, p2_umap, labels = c("A", "B"))

# Plot t-SNE results
p1_tsne <- DimPlot(object = gbm, reduction = "tsne")
p2_tsne <- DimPlot(object = gbm,
                   reduction = "tsne",
                   group.by = "patient.name")

# Combine t-SNE plots into a grid
plot_grid(p1_tsne, p2_tsne, labels = c("C", "D"))