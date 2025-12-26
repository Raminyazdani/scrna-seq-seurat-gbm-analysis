# Development History: Single-Cell RNA-Seq Analysis Pipeline

This document reconstructs the development history of this single-cell RNA sequencing analysis pipeline, showing how it evolved from initial setup to a comprehensive, portfolio-ready analysis tool.

## History Expansion Note

**Previous version:** 8 steps  
**Current version:** 12 steps  
**Expansion multiplier:** 1.5×  

This is an expanded version of the development history with more granular commit stages to better demonstrate realistic incremental development. The final state (step_12) remains identical to the previous final state (old step_08).

### Mapping from Previous to Current Steps

| Previous Steps | Current Steps | Description |
|----------------|---------------|-------------|
| Old Step 01 | Step 01 | Initial repository setup (unchanged) |
| Old Step 02 | Steps 02-03 | **Split**: Package management → Data loading |
| Old Step 03 | Step 04 | QC implementation (unchanged) |
| Old Step 04 | Step 05 | Doublet detection (unchanged) |
| Old Step 05 | Steps 06-07 | **Split + Hotfix**: Normalization bug → Fix |
| Old Step 06 | Steps 08-09 | **Split**: PCA → Clustering |
| Old Step 07 | Step 10 | Cell type annotation (unchanged) |
| Old Step 08 | Steps 11-12 | **Split**: Enrichment → Portfolio refinement |

### Oops → Hotfix Sequence

**Step 06 (Bug):** Attempted to use `sctransform` package without installing it first, causing an import error when trying to run normalization.

**What broke:** The script would fail with "Error: package 'sctransform' not found" when attempting normalization.

**How noticed:** During manual testing of the normalization step, the script crashed immediately when trying to load sctransform.

**Step 07 (Fix):** Added `install_if_missing("sctransform", "CRAN")` and `install_if_missing("glmGamPoi", "Bioconductor")` before attempting to load these packages. Also added the glmGamPoi backend for faster computation.

**Result:** Normalization now works reliably with proper dependency management.

---

## Overview

This project evolved through 12 development stages, from initial repository setup through data exploration, quality control implementation, advanced analysis features, and final portfolio refinement.

**Total Development Time (Simulated):** ~5 weeks  
**Language:** R  
**Primary Framework:** Seurat + Bioconductor  
**Final Line Count:** ~850 lines (main analysis)

---

## Step 01: Initial Repository Setup
**Date:** Week 1, Day 1  
**Commit Message:** "Initial commit: Set up repository structure"

### What Was Added
- Basic repository structure
- README.md with project description
- .gitignore for R projects
- .github directory for workflows

### Rationale
Every project starts with proper setup. Created the foundational structure with appropriate gitignore patterns for R (ignoring .Rdata, .Rhistory, etc.) and a README outlining the project goals.

**Files:** README.md, .gitignore, .github/

---

## Step 02: Package Management System
**Date:** Week 1, Day 2  
**Commit Message:** "Add automatic package installation system"

### What Was Added
- `install_if_missing()` helper function
- Automatic BiocManager setup
- Package installation for CRAN packages (dplyr, tidyverse, ggplot2, patchwork, cowplot)
- Package installation for Bioconductor packages (GEOquery, Seurat, SingleCellExperiment)
- Project.R initialization

### Rationale
Reproducibility is critical for bioinformatics pipelines. Implemented a robust package management system that automatically installs missing dependencies from both CRAN and Bioconductor. This ensures anyone can run the pipeline without manual package installation.

**Key Features:**
- Automatic BiocManager setup
- Dual-source installation (CRAN + Bioconductor)
- Fallback error handling
- Silent loading with `quietly = TRUE`

**Files Added:** Project.R (initial)

---

## Step 03: GEO Data Loading
**Date:** Week 1, Day 3  
**Commit Message:** "Add GEO dataset download and expression matrix loading"

### What Was Added
- GEOquery integration for data download
- Data directory structure creation
- Expression matrix loading from GSE75688
- Local data caching mechanism
- LARGE_FILES.md documentation

### Rationale
Separated data loading from package management for cleaner commit history. Implemented GEO data download with local caching to avoid repeated downloads and ensure consistent data across runs.

**Key Features:**
- Automatic data/ directory creation
- GEO dataset GSE75688 download
- Expression matrix extraction
- Local caching for efficiency
- Documentation of external files

**Files Added:** data/ directory, LARGE_FILES.md  
**Files Modified:** Project.R (added data loading section)

---

## Step 04: Quality Control Implementation
**Date:** Week 1, Day 4-5  
**Commit Message:** "Implement comprehensive QC metrics and visualization"

### What Was Added
- QC metric calculation (nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, etc.)
- Seurat object creation
- PercentageFeatureSet for feature-specific metrics
- Statistical summaries (mean, median, IQR, MAD)
- Violin plots for QC metrics
- Feature scatter plots
- Distribution visualizations (histograms, density plots, box plots)

### Rationale
Quality control is the foundation of good scRNA-seq analysis. Implemented multiple QC metrics to identify low-quality cells, dying cells (high mitochondrial %), and potential technical artifacts. Comprehensive visualizations enable manual inspection of data quality.

**Key Functions:**
- PercentageFeatureSet for mitochondrial genes (^MT-)
- PercentageFeatureSet for ribosomal genes (^RP[SL])
- PercentageFeatureSet for hemoglobin genes (^HB[^(P)])
- calculate_stats function for statistical summaries
- Multi-panel visualization with ggplot2 and patchwork

**Files Modified:** Project.R (added QC section, ~100 lines)

---

## Step 05: Doublet Detection and Filtering
**Date:** Week 2, Day 1-2  
**Commit Message:** "Add doublet detection with scDblFinder and outlier filtering"

### What Was Added
- scDblFinder package integration
- SingleCellExperiment object conversion
- Automated doublet detection
- IQR-based outlier detection
- MAD-based (Median Absolute Deviation) outlier detection
- Multi-stage filtering pipeline
- Updated statistical summaries post-filtering

### Rationale
Doublets (two cells captured as one droplet) are a common technical artifact in droplet-based scRNA-seq that can create false cell types and confound downstream analysis. Implemented scDblFinder for automatic doublet detection and robust outlier filtering using both IQR and MAD methods to identify and remove low-quality cells.

**Key Features:**
- Automated doublet detection with scDblFinder
- IQR thresholds for nCount_RNA and nFeature_RNA
- MAD-based thresholds for multiple QC metrics
- Progressive filtering with cell count tracking
- Preservation of singlets only

**Files Modified:** Project.R (added doublet detection, ~80 lines)

---

## Step 06: Normalization Implementation (With Bug)
**Date:** Week 2, Day 3  
**Commit Message:** "Add SCTransform normalization"

### What Was Added
- SCTransform normalization code
- Variable feature selection preparation

### What Went Wrong
**BUG:** Forgot to add `install_if_missing()` calls for `sctransform` and `glmGamPoi` packages before using them. The code directly calls `library(sctransform)` without ensuring it's installed, which causes the script to crash with "Error: package 'sctransform' not found" for users who haven't manually installed these packages.

### Rationale
Attempted to implement state-of-the-art normalization using SCTransform, which models technical variation more accurately than traditional log-normalization. However, rushed the implementation and forgot to follow the established pattern of using `install_if_missing()` for all packages.

**Files Modified:** Project.R (added normalization section with bug)

---

## Step 07: Hotfix - Add Missing Normalization Dependencies
**Date:** Week 2, Day 3 (later)  
**Commit Message:** "Hotfix: Add missing sctransform and glmGamPoi dependencies"

### What Was Fixed
- Added `install_if_missing("sctransform", "CRAN")` before usage
- Added `install_if_missing("glmGamPoi", "Bioconductor")` for faster backend
- Proper library loading after installation
- Updated SCTransform call to use glmGamPoi method

### How the Bug Was Discovered
During manual testing of the pipeline on a fresh R installation, the script crashed when reaching the normalization step with:
```
Error in library(sctransform) : there is no package called 'sctransform'
```

### The Fix
Added proper dependency management following the project's established pattern:
```r
install_if_missing("sctransform", "CRAN")
install_if_missing("glmGamPoi", "Bioconductor")
library(sctransform)
library(glmGamPoi)
```

Also improved the implementation by adding the glmGamPoi backend for significantly faster computation on large datasets.

**Key Improvements:**
- Consistent dependency management
- Faster normalization with glmGamPoi backend
- More robust for new users

**Files Modified:** Project.R (fixed normalization dependencies, ~15 lines changed)

---

## Step 08: PCA and Dimensionality Reduction
**Date:** Week 2, Day 4  
**Commit Message:** "Implement PCA with optimal component selection"

### What Was Added
- PCA computation with RunPCA
- Automatic optimal PCA component selection
- Cumulative variance threshold (90%)
- Elbow plot generation
- All pairwise PCA dimension plots
- Cowplot for multi-panel visualization

### Rationale
Dimensionality reduction via PCA is essential for reducing computational complexity and noise in high-dimensional scRNA-seq data. Implemented automatic selection of optimal principal components based on cumulative variance explained (90% threshold), ensuring we capture most signal while reducing noise.

**Key Features:**
- RunPCA on highly variable features
- Iterative refinement of PC count
- Elbow plot for visual inspection
- Pairwise dimension plots colored by patient
- Automatic optimal PC determination

**Files Modified:** Project.R (added PCA section, ~60 lines)

---

## Step 09: Clustering and 2D Visualization
**Date:** Week 2, Day 5 - Week 3, Day 1  
**Commit Message:** "Add graph-based clustering with optimal resolution and UMAP/t-SNE"

### What Was Added
- Graph-based k-nearest neighbor finding
- Louvain clustering algorithm
- Automatic optimal resolution selection via modularity
- Binary search for resolution optimization
- UMAP dimensionality reduction
- t-SNE dimensionality reduction
- Side-by-side comparison plots
- Cluster coloring and visualization

### Rationale
Clustering identifies distinct cell populations in the data. Implemented graph-based clustering with an algorithm to automatically find optimal resolution by maximizing modularity while targeting a reasonable number of communities. Added both UMAP and t-SNE for 2D visualization, as they can reveal different aspects of the data structure.

**Key Features:**
- FindNeighbors for KNN graph construction
- FindClusters with Louvain method (igraph)
- Modularity-based resolution optimization
- Binary search algorithm for efficiency
- UMAP with consistent parameters (30 dims)
- t-SNE with perplexity optimization
- Visualization colored by cluster and patient

**Files Modified:** Project.R (added clustering and 2D visualization, ~90 lines)

---

## Step 10: Cell Type Annotation and Marker Analysis
**Date:** Week 3, Day 2-4  
**Commit Message:** "Implement SingleR cell type annotation and marker gene identification"

### What Was Added
- SingleR package integration
- celldex reference database loading
- Human Primary Cell Atlas reference
- Automated cell type annotation
- FindAllMarkers for cluster-specific genes
- FindMarkers for pairwise comparisons (BC07 vs BC07LN)
- Heatmaps of top marker genes
- Cell type composition bar plots
- MAST test for differential expression

### Rationale
Identifying what cell types are present is crucial for biological interpretation. Integrated SingleR for automated annotation using the Human Primary Cell Atlas reference, which contains expression profiles for major human cell types. Implemented comprehensive marker gene analysis to identify cluster-specific genes and validate annotations.

**Key Features:**
- SingleR with Human Primary Cell Atlas
- Automated annotation based on reference correlation
- FindAllMarkers with Wilcoxon test
- Top marker genes per cluster (ranked by avg_log2FC)
- Heatmaps with pheatmap
- BC07 primary vs BC07LN lymph node comparison
- MAST for differential expression testing

**Files Modified:** Project.R (added annotation and markers, ~120 lines)

---

## Step 11: Differential Expression and Pathway Enrichment
**Date:** Week 3, Day 5 - Week 4, Day 2  
**Commit Message:** "Add GO enrichment analysis and GSEA for pathway identification"

### What Was Added
- clusterProfiler package integration
- org.Hs.eg.db for gene annotation
- Gene Ontology (GO) enrichment analysis
- GSEA (Gene Set Enrichment Analysis)
- Multiple pathway database queries (KEGG, Reactome, GO, WikiPathway)
- enrichR integration for additional databases
- Enrichment dot plots
- Gene ranking by log fold change
- Pathway visualization

### Rationale
Understanding biological pathways helps interpret marker genes and cluster identities. Added comprehensive enrichment analysis using clusterProfiler and enrichR to query multiple pathway databases. GSEA considers all genes in ranked order (not just significant ones), providing more complete pathway insights.

**Key Features:**
- clusterProfiler for GO enrichment
- enrichGO for Gene Ontology terms (BP, MF, CC)
- enrichKEGG for pathway analysis
- GSEA with proper gene ranking
- enrichR for multiple databases simultaneously
- Dot plots for visualization
- Adjusted p-value filtering

**Files Modified:** Project.R (added enrichment analysis, ~100 lines)

---

## Step 12: Portfolio Refinement and Documentation
**Date:** Week 4, Day 3 - Week 5  
**Commit Message:** "Refine for portfolio: documentation, path robustness, professional presentation"

### What Was Added
- Comprehensive README.md (300+ lines)
- Professional project_identity.md
- Script header documentation (Project.R and exploratory_analysis.R)
- Path robustness improvements (RStudio + command line compatibility)
- Academic trace removal ("University Project", "final_project" references)
- Archive directory for assignment materials
- suggestion.txt and suggestions_done.txt ledgers
- report.md with complete execution log
- LARGE_FILES.md external file documentation
- Troubleshooting guide in README
- License and citation sections

### What Was Changed
- README title: "Single Cell RNA Sequencing Project" → "Single-Cell RNA-Seq Analysis of Glioblastoma with Seurat"
- Removed "University Project" label
- Changed metadata study field from "final_project" to "gbm_scrnaseq_analysis"
- Fixed RStudio-specific `setwd()` calls with conditional checks
- Renamed test.R → exploratory_analysis.R for clarity
- Moved scRNA-seq_project.docx → archive/
- Enhanced all documentation to portfolio grade

### Rationale
Prepared the entire codebase for professional portfolio presentation. Removed all academic traces, improved documentation comprehensively, fixed path issues to work from both RStudio and command line, and organized files professionally. The code functionality remains identical, but presentation and usability are significantly improved.

**Key Improvements:**
- Portfolio-grade documentation
- Robust path handling for multiple environments
- Professional naming and organization
- Complete usage instructions
- Comprehensive troubleshooting guide
- Academic trace removal
- Archived assignment materials

**Files Modified:**
- README.md (complete rewrite, 316 lines)
- Project.R (added header, fixed paths, ~30 lines changed)
- exploratory_analysis.R (renamed from test.R, added header, fixed paths)
- project_identity.md (created, 57 lines)
- report.md (created, 337 lines)
- suggestion.txt (created, 13 lines)
- suggestions_done.txt (created, 21 lines)

**Files Archived:** scRNA-seq_project.docx → archive/  
**Files Created:** archive/README.md (explanation of archived content)

---

## Development Insights

### Technical Decisions

1. **SCTransform over traditional normalization:** More accurate for scRNA-seq, handles heteroscedasticity better
2. **Automatic parameter selection:** Makes pipeline more generalizable to other datasets
3. **Multiple QC strategies:** Catches different types of low-quality cells (IQR + MAD)
4. **Caching strategy:** Speeds up re-runs significantly by saving intermediate RDS files
5. **Visualization emphasis:** Critical for understanding results and manual validation
6. **Dual package sources:** CRAN + Bioconductor ensures all dependencies are available

### Challenges Overcome

1. **Memory management:** Large expression matrices require sparse matrices and careful object handling
2. **Runtime optimization:** Added caching, used efficient backends (glmGamPoi), selected optimal parameters
3. **Reproducibility:** Automatic package installation, data caching, consistent parameters
4. **Usability:** Works from both RStudio and command line after path robustness fixes
5. **Dependency management:** Required both CRAN and Bioconductor packages with proper fallbacks

### Lessons Learned

1. **Always use consistent patterns:** The bug in Step 06 was caused by breaking the established `install_if_missing()` pattern
2. **Test in clean environments:** Caught the missing dependency bug by testing on a fresh R installation
3. **Split large commits:** Breaking Step 02 into package management + data loading improved clarity
4. **Document as you go:** Final documentation phase was easier because intermediate steps were well-organized

### Future Enhancements (Not Implemented)

- Integration of multiple samples with batch correction (Harmony - commented out in code)
- Interactive visualization with Shiny for exploration
- Automated report generation with R Markdown
- Docker containerization for full reproducibility
- Additional cell type markers for manual validation
- Trajectory analysis for differentiation studies

---

## Snapshot Directory Structure

```
history/
├── github_steps.md                # This file
├── _previous_run/                 # Archived 8-step version
│   ├── github_steps.md
│   └── steps/
│       ├── step_01 ... step_08
└── steps/
    ├── step_01/                   # Initial setup
    ├── step_02/                   # Package management
    ├── step_03/                   # Data loading
    ├── step_04/                   # QC implementation
    ├── step_05/                   # Doublet detection
    ├── step_06/                   # Normalization (with bug)
    ├── step_07/                   # Normalization hotfix
    ├── step_08/                   # PCA
    ├── step_09/                   # Clustering
    ├── step_10/                   # Cell type annotation
    ├── step_11/                   # Pathway enrichment
    └── step_12/                   # Portfolio refinement (final state)
```

Each step directory contains a complete snapshot of the repository at that stage of development, excluding the `history/` directory itself to avoid recursion, and excluding the `data/` directory to keep repository size manageable.

---

## How to Use These Snapshots

Each snapshot represents a working state of the repository at a specific development stage. To explore a particular stage:

```bash
# View files in a specific step
ls -la history/steps/step_05/

# Compare two steps to see what changed
diff -r history/steps/step_06/ history/steps/step_07/

# Extract a specific step for exploration
cp -r history/steps/step_08/ /tmp/step_08_extracted/
cd /tmp/step_08_extracted/

# Compare the bug version with the fix
diff history/steps/step_06/Project.R history/steps/step_07/Project.R
```

### Examining the Oops → Hotfix Sequence

To see the specific bug and fix:

```bash
# Show the buggy normalization code
grep -A 5 "library(sctransform)" history/steps/step_06/Project.R

# Show the fixed normalization code
grep -A 5 "install_if_missing.*sctransform" history/steps/step_07/Project.R

# See the complete diff
diff -u history/steps/step_06/Project.R history/steps/step_07/Project.R
```

---

## Final State Verification

The final snapshot (step_12) matches the current working tree exactly (excluding history/ and data/ directories):

```bash
# Verify final state matches current repo
diff -r . history/steps/step_12/ --exclude=.git --exclude=history --exclude=data
# Should show no differences (except SNAPSHOT_INFO.md which is snapshot-specific)
```

---

**Note:** This development history is a reconstruction for portfolio purposes, demonstrating how a complex bioinformatics pipeline might realistically evolve through iterative, incremental development with occasional bugs and fixes.
