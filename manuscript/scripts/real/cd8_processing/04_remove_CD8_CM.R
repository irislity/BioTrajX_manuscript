# =============================================================================
# 04_remove_CD8_CM.R — drop CD8.CM cells from the CD8T dataset
#
# Usage:
#   Rscript manuscript/scripts/real/cd8_processing/04_remove_CD8_CM.R
# =============================================================================

library(Seurat)

repo_root <- here::here()
data_path <- file.path(repo_root, "data", "cd8t.rds")

seurat <- readRDS(data_path)
message("Before:")
print(table(seurat$functional.cluster))

seurat <- subset(seurat, functional.cluster != "CD8.CM")
seurat$functional.cluster <- droplevels(factor(seurat$functional.cluster))

message("\nAfter:")
print(table(seurat$functional.cluster))

saveRDS(seurat, data_path)
message(sprintf("\nSaved: %s (%d cells)", data_path, ncol(seurat)))
