# Snapshot Info: Step 09

**Commit:** HOTFIX - Add missing dependency installation

**Changes:**
- Added install_if_missing("sctransform", "CRAN")
- Added install_if_missing("glmGamPoi", "Bioconductor")  
- Fixed normalization to work on fresh installations

**Status:** FIXED - Normalization now works reliably
**Fix for:** Step 08 bug - missing sctransform dependencies
