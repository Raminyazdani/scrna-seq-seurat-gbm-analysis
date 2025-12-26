# [Package installation and data loading code from step_02...]
# [... plus QC metrics calculation ...]

# Calculate QC metrics
gbm[["percent.mt"]] <- PercentageFeatureSet(gbm, pattern = "^MT-")
gbm[["percent.ERCC"]] <- PercentageFeatureSet(gbm, pattern = "^ERCC")
gbm[["percent.ribo"]] <- PercentageFeatureSet(gbm, pattern = "^RPS |^RPL")
gbm[["percent.hb"]] <- PercentageFeatureSet(gbm, pattern = "^HB[^(P)]")

# Violin plots
violin_plots <- VlnPlot(gbm, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)

# [Additional QC visualizations...]
message("QC metrics calculated and visualized!")
