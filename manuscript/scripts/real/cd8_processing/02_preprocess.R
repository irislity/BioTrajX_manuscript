# =============================================================================
# 02_preprocess.R — QC, normalization, and CD8+ T cell subsetting
#
# Usage:
#   Rscript manuscript/scripts/real/cd8_processing/02_preprocess.R
# =============================================================================

library(Seurat)

OUT_DIR <- "output"

SO.GSE164690.PBL   <- readRDS(file.path(OUT_DIR, "SO.GSE164690.PBL.rds"))
SO.GSE164690.cd45p <- readRDS(file.path(OUT_DIR, "SO.GSE164690.cd45p.rds"))

# ── QC ────────────────────────────────────────────────────────────────────────
so_list <- list(
  PBL   = SO.GSE164690.PBL,
  cd45p = SO.GSE164690.cd45p
)

so_list <- lapply(so_list, PercentageFeatureSet, pattern = "^MT-", col.name = "percent.mt")

VlnPlot(so_list$PBL,   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(so_list$cd45p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ── Filter ────────────────────────────────────────────────────────────────────
# Thresholds chosen from pre-filter QC plots above
filter_so <- function(so) {
  subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
                      nCount_RNA < 30000 & percent.mt < 10)
}

so_list <- lapply(so_list, filter_so)

VlnPlot(so_list$PBL,   features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(so_list$cd45p, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ── Normalize ─────────────────────────────────────────────────────────────────
so_list <- lapply(so_list, NormalizeData)
so_list <- lapply(so_list, FindVariableFeatures, selection.method = "vst", nfeatures = 2000)
so_list <- lapply(so_list, ScaleData)

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(so_list$PBL,   file.path(OUT_DIR, "GSE164690.PBL.norm.rds"))
saveRDS(so_list$cd45p, file.path(OUT_DIR, "GSE164690.cd45p.norm.rds"))
