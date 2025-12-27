# Development History: Single-Cell RNA-Seq Analysis Pipeline

This document reconstructs the development history of this single-cell RNA sequencing analysis pipeline, showing how it evolved from initial setup to a comprehensive, portfolio-ready analysis tool.

## History Expansion Note

**Previous version (v1):** 8 steps  
**Previous version (v2):** 12 steps  
**Current version (v3):** 18 steps  
**Expansion multiplier (v2→v3):** 1.5×  
**Total expansion (v1→v3):** 2.25×

This is the second expansion of the development history, now with even more granular commit stages to demonstrate realistic incremental development with multiple iterations, bug discovery, and fixes. The final state (step_18) remains functionally identical to previous final states, but the journey shows more realistic development patterns including:

- More granular feature implementation splits
- Two distinct oops→hotfix sequences  
- Iterative refinement of documentation
- Realistic bug discovery and resolution patterns

### Mapping from Previous v2 (12 steps) to Current v3 (18 steps)

| Previous v2 Steps | Current v3 Steps | Description |
|-------------------|------------------|-------------|
| Old Step 01 | Step 01 | Initial repository setup (unchanged) |
| Old Step 02 | Steps 02-03 | **Split**: Core packages → Bioconductor packages |
| Old Step 03 | Step 04 | Data loading (unchanged) |
| Old Step 04 | Steps 05-06 | **Split**: QC calculation → QC visualization |
| Old Step 05 | Step 07 | Doublet detection (unchanged) |
| Old Step 06 | Steps 08-09 | **Oops→Hotfix #1**: Normalization bug → Dependency fix |
| Old Step 07 | (merged) | Hotfix merged into step 09 |
| Old Step 08 | Step 10 | PCA (unchanged) |
| Old Step 09 | Steps 11-12 | **Oops→Hotfix #2**: Clustering typo → Parameter fix |
| Old Step 10 | Steps 13-14 | **Split**: Hotfix separation → Cell annotation |
| Old Step 11 | Steps 15-16 | **Split**: Marker genes → Pathway enrichment |
| Old Step 12 | Steps 17-18 | **Split**: Documentation → Final polish |

### Oops → Hotfix Sequences (2 total)

#### Sequence #1: Normalization Dependency Bug (Steps 08-09)

**Step 08 (Bug):** Added normalization code using `library(sctransform)` and `library(glmGamPoi)` without first calling `install_if_missing()` for these packages.

**What broke:** Script would fail with "Error: package 'sctransform' not found" or "Error: package 'glmGamPoi' not found" when attempting to run normalization on a fresh R installation.

**How noticed:** During testing on a clean R environment (without pre-installed packages), the script crashed immediately when trying to load sctransform. The error message was clear: package not available.

**Step 09 (Fix):** Added proper dependency installation:
- `install_if_missing("sctransform", "CRAN")`  
- `install_if_missing("glmGamPoi", "Bioconductor")`
- Placed these calls before the corresponding `library()` statements
- Followed the established pattern used for all other packages

**Result:** Normalization now works reliably on fresh installations with automatic dependency management.

**Lesson:** Consistency matters - breaking the established `install_if_missing()` pattern caused an avoidable bug.

#### Sequence #2: Clustering Parameter Typo (Steps 11-13)

**Step 11 (Bug):** Implemented clustering but made a typo in the `FindClusters()` call: used `resoltuion = 0.5` instead of `resolution = 0.5`.

**What broke:** The clustering function either failed with an "unknown parameter" error or silently ignored the typo and used the default resolution, leading to unexpected cluster counts.

**How noticed:** After running clustering, the number of clusters didn't match expectations. Reviewing the code revealed the typo in the parameter name.

**Step 12-13 (Fix):** 
- Step 12: Added proper UMAP/t-SNE visualizations and began optimization
- Step 13: Explicit hotfix commit correcting the 'resoltuion' typo to 'resolution'
- Verified that clusters now form properly with the intended resolution parameter

**Result:** Clustering operates correctly with proper parameterization.

**Lesson:** Typos in parameter names can cause silent failures or unexpected behavior. Code review and testing are essential.

---

## Overview

This project evolved through 18 development stages over simulated 6-7 weeks, from initial repository setup through data exploration, quality control implementation, advanced analysis features, bug fixes, and final portfolio refinement.

**Total Development Time (Simulated):** ~6-7 weeks  
**Language:** R  
**Primary Framework:** Seurat + Bioconductor  
**Final Line Count:** ~850 lines (main analysis)  
**Bugs Found & Fixed:** 2 (dependency management, parameter typo)

---

## Detailed Step-by-Step Development History

### Step 01: Initial Repository Setup
**Commit:** "Initial commit: repository structure"

**Changes:**
- Created repository with basic structure
- Added README.md with project description
- Created .gitignore for R projects
- Added .github/ directory with issue templates
- Established data/ directory placeholder

**Files:** README.md, .gitignore, .github/

---

### Step 02: Core Package Installation
**Commit:** "Add core package installation framework"

**Changes:**
- Created Project.R with package management
- Implemented install_if_missing() helper
- Added core packages: Seurat, dplyr, ggplot2
- Set up automatic installation

---

### Step 03: Bioconductor Packages
**Commit:** "Add Bioconductor support"

**Changes:**
- Extended install_if_missing() for Bioconductor
- Added GEOquery, SingleCellExperiment, scDblFinder
- Added BiocManager installation

---

### Step 04: Data Loading
**Commit:** "Implement GEO data download"

**Changes:**
- Added GEOquery functionality (GSE75688)
- Implemented expression matrix loading
- Created data caching
- Parsed sample metadata

---

### Step 05: QC Metrics Calculation
**Commit:** "Calculate QC metrics"

**Changes:**
- Created Seurat object
- Calculated percent.mt, percent.ribo
- Computed nCount_RNA, nFeature_RNA

---

### Step 06: QC Visualization
**Commit:** "Add QC visualization and filtering"

**Changes:**
- Created violin plots
- Implemented IQR/MAD filtering
- Generated scatter plots
- Applied quality filters

---

### Step 07: Doublet Detection
**Commit:** "Implement doublet detection"

**Changes:**
- Added scDblFinder algorithm
- Identified and removed doublets
- Generated visualizations

---

### Step 08: Normalization (WITH BUG)
**Commit:** "Add SCTransform normalization"

**Changes:**
- Added SCTransform code
- Configured glmGamPoi backend
- **BUG:** Missing install_if_missing() calls

**Status:** BROKEN

---

### Step 09: Normalization Hotfix
**Commit:** "HOTFIX: Add missing dependencies"

**Changes:**
- Added install_if_missing("sctransform")
- Added install_if_missing("glmGamPoi")
- Fixed installation pattern

**Status:** FIXED

---

### Step 10: PCA
**Commit:** "Implement PCA"

**Changes:**
- Ran PCA on normalized data
- Created elbow plot
- Selected 30 PCs
- Generated visualizations

---

### Step 11: Clustering (WITH BUG)
**Commit:** "Add clustering"

**Changes:**
- Implemented FindNeighbors
- Applied FindClusters
- **BUG:** Typo 'resoltuion' instead of 'resolution'

**Status:** BROKEN

---

### Step 12: Clustering Optimization
**Commit:** "Add UMAP/t-SNE visualizations"

**Changes:**
- Added UMAP reduction
- Added t-SNE reduction
- Created cluster visualizations
- Began parameter optimization

---

### Step 13: Clustering Hotfix
**Commit:** "HOTFIX: Fix parameter typo"

**Changes:**
- Corrected 'resoltuion' to 'resolution'
- Verified proper parameterization
- Tested cluster formation

**Status:** FIXED

---

### Step 14: Cell Type Annotation
**Commit:** "Implement cell type annotation"

**Changes:**
- Integrated SingleR
- Used Human Primary Cell Atlas
- Assigned cell types
- Created composition plots

---

### Step 15: Marker Genes
**Commit:** "Identify marker genes"

**Changes:**
- Implemented FindAllMarkers
- Identified top markers per cluster
- Calculated DE statistics
- Created ranked gene lists

---

### Step 16: Pathway Enrichment
**Commit:** "Add pathway enrichment"

**Changes:**
- Implemented clusterProfiler
- Queried GO, KEGG, Reactome
- Generated enrichment plots
- Performed GSEA

---

### Step 17: Documentation Setup
**Commit:** "Add portfolio documentation"

**Changes:**
- Created project_identity.md
- Created report.md
- Created suggestion.txt
- Created suggestions_done.txt
- Added file headers

---

### Step 18: Final Portfolio Refinement
**Commit:** "Final polish: complete documentation"

**Changes:**
- Completed all portfolio files
- Finalized README (316 lines)
- Renamed test.R to exploratory_analysis.R
- Fixed RStudio/CLI path handling
- Removed academic traces
- Added archive/ directory
- Enhanced all documentation sections

**Status:** FINAL - Portfolio ready

---

## Summary

This 18-step history demonstrates:
- ✓ Realistic incremental development
- ✓ Two oops→hotfix sequences
- ✓ Iterative documentation refinement
- ✓ 1.5× expansion from v2 (12 steps)
- ✓ Final state matches working tree
- ✓ No .git or history in snapshots
- ✓ All portfolio deliverables complete

---

## Verification Commands

```bash
# Count steps
ls -d history/steps/step_* | wc -l
# Result: 18

# Verify final state
diff -r . history/steps/step_18/ --exclude=.git --exclude=history --exclude=data --exclude=figures
# Result: No differences (except SNAPSHOT_INFO.md)

# Check for .git/history in snapshots
find history/steps/step_18 -name ".git" -o -name "history"
# Result: None found
```
