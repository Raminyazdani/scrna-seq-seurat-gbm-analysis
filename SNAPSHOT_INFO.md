# Snapshot Info: Step 08

**Commit:** Add normalization (WITH BUG)

**Changes:**
- Added SCTransform normalization
- Configured glmGamPoi backend
- **BUG:** Used library(sctransform) and library(glmGamPoi) without install_if_missing()

**Status:** BROKEN - Will fail on fresh R installations
**Bug:** Missing dependency installation for sctransform and glmGamPoi packages
