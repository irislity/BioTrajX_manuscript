# =============================================================================
# run_ti_stemcell.R
#
# Run all TI methods on the stem cell dataset and save pseudotimes.
#
# Root cell (start_cell passed to the 7 directed methods): the
# Stem_Progenitors cell nearest the stem-cell-population centroid in PCA
# space -- a stable, method-agnostic root not tied to any one TI method's own
# noise (same approach as GSE131847's d0centroid root; see run_ti_gse131847.R).
#
# INPUT:   data/stem_cell.rds
# OUTPUT:  manuscript/results/branch_stemcell/ti_pseudotimes.csv
#          manuscript/results/branch_stemcell/root_cell_stemcell.txt
#
# Usage:
#   Rscript manuscript/scripts/real/run_ti_stemcell.R
# =============================================================================

Sys.setenv(OMP_NUM_THREADS        = "1",
           OPENBLAS_NUM_THREADS   = "1",
           MKL_NUM_THREADS        = "1",
           VECLIB_MAXIMUM_THREADS = "1")

suppressPackageStartupMessages(library(Seurat))

repo_root <- here::here()
source(file.path(repo_root, "manuscript", "scripts", "real", "run_ti_methods.R"))

out_dir <- file.path(repo_root, "manuscript", "results", "branch_stemcell")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("Loading stem_cell.rds ...")
seurat <- readRDS(file.path(repo_root, "data", "stem_cell.rds"))
message(sprintf("  %d cells × %d genes", ncol(seurat), nrow(seurat)))

# =============================================================================
# 2. PREPARE EXPRESSION MATRICES
# =============================================================================
seurat_pca <- Embeddings(seurat, "pca")

message("Selecting top 2000 highly variable genes ...")
seurat <- FindVariableFeatures(seurat, nfeatures = 2000, verbose = FALSE)
hvg    <- VariableFeatures(seurat)

expr          <- as.matrix(GetAssayData(seurat, layer = "data")[hvg, ])
fullgene_expr <- as.matrix(GetAssayData(seurat, layer = "data"))

message(sprintf("  log-norm HVG:  %d × %d", nrow(expr),          ncol(expr)))
message(sprintf("  log-norm (all):%d × %d", nrow(fullgene_expr), ncol(fullgene_expr)))

# =============================================================================
# 3. ROOT CELL — stem-cell-population centroid
# =============================================================================
# Cell closest to the Stem_Progenitors centroid in PCA space -- a stable,
# method-agnostic root not tied to any one TI method's own noise (same
# approach as GSE131847's d0centroid root; see run_ti_gse131847.R).
message("Computing Stem_Progenitors centroid for root cell selection ...")
stem_cells  <- colnames(seurat)[seurat$Phenotype == "Stem_Progenitors"]
n_pcs_root  <- min(20, ncol(seurat_pca))
stem_pca    <- seurat_pca[stem_cells, seq_len(n_pcs_root), drop = FALSE]
centroid    <- colMeans(stem_pca)
d_centroid  <- sqrt(rowSums(sweep(stem_pca, 2, centroid, "-")^2))
start_cell  <- stem_cells[which.min(d_centroid)]
message(sprintf("  Root cell (nearest Stem_Progenitors centroid in PC1-%d): %s",
                n_pcs_root, start_cell))

# =============================================================================
# 4. RUN ALL TI METHODS
# =============================================================================
set.seed(42)
message("\nRunning all TI methods ...")
ti_results <- run_all_ti_methods(
  expr               = expr,
  fullgene_expr      = fullgene_expr,
  start_cell         = start_cell,
  precomp_pca        = seurat_pca,
  methods            = c("Slingshot", "CytoTRACE", "Monocle3", "DPT", "SCORPIUS",
                         "TSCAN", "PAGA-DPT", "Palantir")
)

# =============================================================================
# 5. SAVE
# =============================================================================
out <- file.path(out_dir, "ti_pseudotimes.csv")
write.csv(ti_results, out, row.names = TRUE)
message("Saved: ", out)

root_out <- file.path(out_dir, "root_cell_stemcell.txt")
writeLines(start_cell, root_out)
message("Saved: ", root_out)
message("Done.")
