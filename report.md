# Portfolio Readiness Report: scrna-seq-seurat-gbm-analysis

## Phase 0: Initial Self-Setup

### Created Files
- report.md (this file)
- Will create: suggestion.txt, suggestions_done.txt, project_identity.md

### Repository Status
- Repository: Raminyazdani/scrna-seq-seurat-gbm-analysis
- Primary Language: R
- Type: Single-cell RNA sequencing analysis using Seurat
- Current state: Contains R scripts for scRNA-seq analysis with some academic language

## Phase 1: Repository Understanding

### Structure Analysis
The repository contains:
- **Project.R** - Main analysis script (824 lines) - comprehensive scRNA-seq pipeline
- **test.R** - Similar analysis script (584 lines) - appears to be exploratory/test version
- **README.md** - Current documentation with academic traces
- **scRNA-seq_project.docx** - Word document (116KB) likely containing assignment details
- **data/** - Contains GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt
- **LARGE_FILES.md** - Documents large files stored externally on Google Drive
- **.Rhistory** - R session history

### Analysis Pipeline Overview
The scripts perform:
1. Automatic package installation (GEOquery, Seurat, SingleCellExperiment, etc.)
2. Data download from GEO (GSE75688 - Breast Cancer dataset)
3. Quality control metrics calculation
4. Doublet detection using scDblFinder
5. Normalization with SCTransform
6. Dimensionality reduction (PCA, UMAP, t-SNE)
7. Clustering with optimal resolution finding
8. Cell type annotation using SingleR
9. Differential expression analysis
10. Gene ontology enrichment analysis

### Key Findings
**Good:**
- Well-structured analysis pipeline
- Comprehensive QC and visualization
- Uses relative paths for data files (./data/)
- Automatic package installation included

**Issues to Address:**
1. RStudio-specific path setting (rstudioapi) breaks command-line execution
2. "University Project" label in README
3. "final_project" metadata reference in code
4. Generic "test.R" filename
5. README folder structure doesn't match actual repo
6. Word document likely contains assignment details
7. Missing file header documentation in scripts

## Phase 2: Pre-Change Audit Completed

All findings documented in suggestion.txt:
- 3 TRACE issues (academic language)
- 2 PATH issues (RStudio-specific code)
- 1 RENAME issue (test.R)
- 3 DOC issues (README alignment, missing docs)
- 1 STRUCTURE issue (Word doc)

Total: 10 issues identified

### Naming Alignment Plan
The current repository structure and naming is mostly professional:
- Repo slug "scrna-seq-seurat-gbm-analysis" is appropriate
- Main script "Project.R" could be more descriptive but is acceptable
- "test.R" should be renamed to be more descriptive

**Minimal Renames Required:**
1. test.R → exploratory_analysis.R (descriptive of its purpose)
2. Potentially archive scRNA-seq_project.docx if it contains assignment details

**No folder restructuring needed** - current structure is clean and simple.

