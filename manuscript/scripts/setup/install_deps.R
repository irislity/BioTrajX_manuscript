# =============================================================================
# install_deps.R
#
# Installs the packages this pipeline needs that CANNOT be captured by
# `renv::restore()` alone: BioTrajX itself, plus the Bioconductor- and
# GitHub-only TI-method packages. renv snapshots record exact source hashes
# from a real, already-installed environment — they can't be hand-written —
# so this script has to actually install these once, after which you should
# re-run `renv::snapshot()` to lock them in for next time.
#
# Order:
#   1. renv::restore()      -- CRAN packages from renv.lock
#   2. Rscript manuscript/scripts/setup/install_deps.R   -- this script
#   3. renv::snapshot()     -- lock the newly-installed packages too
#
# Usage (from repo root):
#   Rscript manuscript/scripts/setup/install_deps.R
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

# ── Bioconductor ─────────────────────────────────────────────────────────────
BiocManager::install(c(
  "SingleCellExperiment",
  "slingshot",
  "destiny",
  "TSCAN",
  "fgsea",
  "AnnotationDbi",
  "org.Hs.eg.db",  # get_markers_msigdb(species = "Homo sapiens") -- S8, S9
  "org.Mm.eg.db"   # get_markers_msigdb(species = "Mus musculus") -- S8, S10
), update = FALSE, ask = FALSE)

# destiny's 'smoother' dependency is also archived off CRAN (same situation
# as SCORPIUS below) -- BiocManager::install() silently skips it rather than
# erroring, so it has to be installed explicitly.
remotes::install_version("smoother", version = "1.3")

# ── CRAN ──────────────────────────────────────────────────────────────────────
install.packages(c(
  "ggh4x",    # facet_grid2(..., axes = "all") for the GAM-facet supplementary panels
  "mgcv",     # GAM backend for geom_smooth(method = "gam"); usually ships with R
              # itself as a "recommended" package -- only needed here if your R
              # build doesn't already have it.
  "msigdbr"   # get_markers_msigdb() -- used by every real-data figure script
              # (S8, S9, S10, S11)
))
# SCORPIUS has been archived off CRAN entirely (no longer in the live repo,
# source or binary) -- install_version() pulls it from the CRAN Archive.
# A plain install.packages("SCORPIUS") 404s with a confusing "error code 56"
# rather than a clear "package not available" message.
remotes::install_version("SCORPIUS", version = "1.0.9")

# ── GitHub-only ───────────────────────────────────────────────────────────────
remotes::install_github("irislity/BioTrajX")
remotes::install_github("cole-trapnell-lab/monocle3")
remotes::install_github("carmonalab/STACAS")     # ProjecTILs dependency, not on CRAN
remotes::install_github("carmonalab/ProjecTILs")
# CytoTRACE (original, Gulati et al. 2020): the official site
# (https://cytotrace.stanford.edu/) is a Shiny app whose download button only
# appears after a license-agreement click-through in the UI, but the tarball
# itself sits at a static, directly-fetchable URL underneath (no auth gate on
# the file itself) -- so a plain download + local install works fine.
# (gunsagargulati/CytoTRACE on GitHub also mirrors an older version, 0.3.1,
# if this URL ever moves.)
tf <- tempfile(fileext = ".tar.gz")
download.file("https://cytotrace.stanford.edu/CytoTRACE_0.3.3.tar.gz", tf, mode = "wb")
install.packages(tf, repos = NULL, type = "source")

message("\nDone. Now run renv::snapshot() to lock these into renv.lock.")
