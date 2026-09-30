# =============================================================================
# run_ti_cd8t.R
#
# Run all TI methods on the CD8 T-cell dataset and save pseudotimes.
#
# Root cell (start_cell passed to the 7 directed methods): the CD8.NaiveLike
# cell nearest the naive-population centroid in PCA space -- a stable,
# method-agnostic root not tied to any one TI method's own noise (same
# approach as GSE131847's d0centroid root; see run_ti_gse131847.R).
#
# INPUT:   data/cd8t.rds
# OUTPUT:  manuscript/results/linear_cd8t/ti_pseudotime_cd8t.csv
#          manuscript/results/linear_cd8t/root_cell_cd8t.txt
#
# Usage:
#   Rscript manuscript/scripts/real/run_ti_cd8t.R
# =============================================================================

Sys.setenv(OMP_NUM_THREADS        = "1",
           OPENBLAS_NUM_THREADS   = "1",
           MKL_NUM_THREADS        = "1",
           VECLIB_MAXIMUM_THREADS = "1")

suppressPackageStartupMessages(library(Seurat))

repo_root <- here::here()
source(file.path(repo_root, "manuscript", "scripts", "real", "run_ti_methods.R"))

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("Loading cd8t.rds ...")
seurat <- readRDS(file.path(repo_root, "data", "cd8t.rds"))
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
# 3. ROOT CELL — naive-population centroid
# =============================================================================
# Naive cell closest to the CD8.NaiveLike centroid in PCA space -- a stable,
# method-agnostic root not tied to any one TI method's own noise (same
# approach as GSE131847's d0centroid root; see run_ti_gse131847.R).
message("Computing CD8.NaiveLike centroid for root cell selection ...")
naive_cells <- colnames(seurat)[seurat$functional.cluster == "CD8.NaiveLike"]
n_pcs_root  <- min(20, ncol(seurat_pca))
naive_pca   <- seurat_pca[naive_cells, seq_len(n_pcs_root), drop = FALSE]
centroid    <- colMeans(naive_pca)
d_centroid  <- sqrt(rowSums(sweep(naive_pca, 2, centroid, "-")^2))
start_cell  <- naive_cells[which.min(d_centroid)]
message(sprintf("  Root cell (nearest CD8.NaiveLike centroid in PC1-%d): %s",
                n_pcs_root, start_cell))

# =============================================================================
# 4. RUN TI METHODS
# =============================================================================
set.seed(42)
message("\nRunning TI methods ...")
ti_results <- run_all_ti_methods(
  expr          = expr,
  fullgene_expr = fullgene_expr,
  start_cell    = start_cell,
  precomp_pca   = seurat_pca,
  methods       = c("Slingshot", "CytoTRACE", "Monocle3", "DPT", "SCORPIUS",
                    "TSCAN", "PAGA-DPT", "Palantir")
)

# =============================================================================
# 5. SAVE
# =============================================================================
out_dir <- file.path(repo_root, "manuscript", "results", "linear_cd8t")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "ti_pseudotime_cd8t.csv")
write.csv(ti_results, out, row.names = TRUE)
message("\nSaved: ", out)

root_out <- file.path(out_dir, "root_cell_cd8t.txt")
writeLines(start_cell, root_out)
message("Saved: ", root_out)
message("Done.")
