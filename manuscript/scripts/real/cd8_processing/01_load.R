# =============================================================================
# 01_load.R — Load GSE164690 raw data
#
# GSE164690: CD8+ T cell exhaustion dataset (PBL + CD45+ tumor-infiltrate
# samples, head & neck cancer). Raw 10x matrices from GEO accession GSE164690.
#
# Usage:
#   Rscript manuscript/scripts/real/cd8_processing/01_load.R
# =============================================================================

library(Seurat)

# Set DATA_DIR to the folder containing the extracted GSE164690_RAW files
DATA_DIR <- "data/GSE164690_RAW"
OUT_DIR  <- "output"
dir.create(OUT_DIR, showWarnings = FALSE)

# ── Load PBL samples ─────────────────────────────────────────────────────────
pbl_prefixes <- c(
  "GSM5017021_HN01_PBL", "GSM5017024_HN02_PBL", "GSM5017026_HN03_PBL",
  "GSM5017028_HN04_PBL", "GSM5017030_HN05_PBL", "GSM5017033_HN06_PBL",
  "GSM5017036_HN07_PBL", "GSM5017039_HN08_PBL", "GSM5017042_HN09_PBL",
  "GSM5017048_HN11_PBL", "GSM5017051_HN12_PBL", "GSM5017054_HN13_PBL",
  "GSM5017057_HN14_PBL", "GSM5017060_HN15_PBL", "GSM5017063_HN16_PBL",
  "GSM5017066_HN17_PBL", "GSM5017069_HN18_PBL"
)

read_mtx <- function(prefix, dir) {
  ReadMtx(
    mtx      = file.path(dir, paste0(prefix, "_matrix.mtx.gz")),
    cells    = file.path(dir, paste0(prefix, "_barcodes.tsv.gz")),
    features = file.path(dir, paste0(prefix, "_features.tsv.gz"))
  )
}

pbl_matrices <- lapply(setNames(pbl_prefixes, pbl_prefixes), read_mtx, dir = DATA_DIR)
SO.GSE164690.PBL <- CreateSeuratObject(pbl_matrices, assay = "RNA")
SO.GSE164690.PBL$orig.ident <- "GSE164690.PBL"

saveRDS(SO.GSE164690.PBL, file.path(OUT_DIR, "SO.GSE164690.PBL.rds"))

# ── Load CD45+ tumor-infiltrate samples ─────────────────────────────────────
cd45p_prefixes <- c(
  "GSM5017022_HN01_CD45p", "GSM5017025_HN02_CD45p", "GSM5017027_HN03_CD45p",
  "GSM5017029_HN04_CD45p", "GSM5017031_HN05_CD45p", "GSM5017034_HN06_CD45p",
  "GSM5017037_HN07_CD45p", "GSM5017040_HN08_CD45p", "GSM5017043_HN09_CD45p",
  # HN10 excluded: no matched PBL sample
  "GSM5017049_HN11_CD45p", "GSM5017052_HN12_CD45p",
  "GSM5017055_HN13_CD45p", "GSM5017058_HN14_CD45p", "GSM5017061_HN15_CD45p",
  "GSM5017064_HN16_CD45p", "GSM5017067_HN17_CD45p", "GSM5017070_HN18_CD45p"
)

cd45p_matrices <- lapply(setNames(cd45p_prefixes, cd45p_prefixes), read_mtx, dir = DATA_DIR)
SO.GSE164690.cd45p <- CreateSeuratObject(cd45p_matrices, assay = "RNA")
SO.GSE164690.cd45p$orig.ident <- "GSE164690.CD45P"

saveRDS(SO.GSE164690.cd45p, file.path(OUT_DIR, "SO.GSE164690.cd45p.rds"))
