# =============================================================================
# 03_ProjecTILs_CytoTRACE.R — ProjecTILs classification + DEG analysis
#
# Usage:
#   Rscript manuscript/scripts/real/cd8_processing/03_ProjecTILs_CytoTRACE.R
# =============================================================================

library(Seurat)
library(ProjecTILs)

OUT_DIR <- "output"

pbl   <- readRDS(file.path(OUT_DIR, "GSE164690.PBL.norm.rds"))
cd45p <- readRDS(file.path(OUT_DIR, "GSE164690.cd45p.norm.rds"))

GSE164690.all <- merge(cd45p, y = pbl)

# ── ProjecTILs classification ────────────────────────────────────────────────
GSE164690.all <- JoinLayers(GSE164690.all)

GSE164690.all.projectils <- ProjecTILs.classifier(
  query = GSE164690.all,
  ref   = CD8T_human_ref_v1
)

GSE164690.all.projectils <- subset(
  GSE164690.all.projectils,
  functional.cluster != "NA"
)

DimPlot(GSE164690.all.projectils,
        group.by = "functional.cluster", label = TRUE, repel = TRUE)

# ── DEG analysis ──────────────────────────────────────────────────────────────
fdr_threshold <- 0.01

deg <- FindMarkers(GSE164690.all.projectils, group.by = "functional.cluster")
deg$p_val_adj.fdr <- p.adjust(deg$p_val, method = "fdr")
deg_fdr <- subset(deg, p_val_adj.fdr <= fdr_threshold)

write.csv(deg_fdr, file.path(OUT_DIR, "projectils_deg_fdr0.01.csv"))

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(GSE164690.all.projectils,
        file.path(OUT_DIR, "GSE164690.all.projectils.rds"))
