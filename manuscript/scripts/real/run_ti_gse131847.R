# =============================================================================
# run_ti_gse131847.R
#
# Run all TI methods on the GSE131847 LCMV time-course CD8 T-cell dataset
# (naive, d3, d4, d5, d6, d7, d10, d14, d21, d32, d60, d90 post-infection)
# and save pseudotimes.
#
# Root cell (start_cell passed to the 7 directed methods): one of four
# modes, selected via a command-line argument --
#   d0centroid (default) -- naive cell nearest the D0/naive centroid in PCA
#                            space (step 3 below)
#   cytonaive             -- naive cell with the min CytoTRACE score
#                            (original approach: start_cell <-
#                            naive_cells[which.min(cytotrace_pt[naive_cells])])
#   cytoglobal             -- cell with the min CytoTRACE score genome-wide,
#                            with no naive restriction (CytoTRACE's own
#                            unrestricted root)
#   dualcyto               -- run each of the 7 directed methods TWICE, once
#                            rooted at cytonaive and once at cytoglobal,
#                            saved as separate columns (suffixed
#                            "_cytonaive" / "_cytoglobal"), alongside a
#                            single root-independent CytoTRACE column --
#                            14 + 1 = 15 trajectories total.
#
# INPUT:   data/GSE131847_seu.rds
# OUTPUT:  manuscript/results/linear_gse131847_<mode>/ti_pseudotime_gse131847.csv
#          manuscript/results/linear_gse131847_<mode>/root_cell_gse131847.txt
#          (dualcyto: root_cells_gse131847.csv instead, with one row per root)
#
# Usage:
#   Rscript manuscript/scripts/real/run_ti_gse131847.R [d0centroid|cytonaive|cytoglobal|dualcyto]
# =============================================================================

args      <- commandArgs(trailingOnly = TRUE)
root_mode <- if (length(args) >= 1) args[1] else "d0centroid"
stopifnot(root_mode %in% c("d0centroid", "cytonaive", "cytoglobal", "dualcyto"))
message("Root mode: ", root_mode)

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
message("Loading GSE131847_seu.rds ...")
seurat <- readRDS(file.path(repo_root, "data", "GSE131847_seu.rds"))
message(sprintf("  %d cells x %d genes", ncol(seurat), nrow(seurat)))
message("  cell_type (day) distribution:")
print(table(seurat$cell_type))

DefaultAssay(seurat) <- "SCT"

# =============================================================================
# 2. PREPARE EXPRESSION MATRICES
# =============================================================================
seurat_pca <- Embeddings(seurat, "pca")

message("Selecting top 2000 highly variable genes ...")
seurat <- FindVariableFeatures(seurat, nfeatures = 2000, verbose = FALSE)
hvg    <- VariableFeatures(seurat)

expr          <- as.matrix(GetAssayData(seurat, assay = "SCT", layer = "data")[hvg, ])
fullgene_expr <- as.matrix(GetAssayData(seurat, assay = "SCT", layer = "data"))

message(sprintf("  log-norm HVG:  %d x %d", nrow(expr),          ncol(expr)))
message(sprintf("  log-norm (all):%d x %d", nrow(fullgene_expr), ncol(fullgene_expr)))

directed_methods <- c("Slingshot", "Monocle3", "DPT", "SCORPIUS", "TSCAN",
                      "PAGA-DPT", "Palantir")
out_dir <- file.path(repo_root, "manuscript", "results", paste0("linear_gse131847_", root_mode))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (root_mode != "dualcyto") {
  # ===========================================================================
  # 3. ROOT CELL -- single root, shared by all 7 directed methods
  # ===========================================================================
  # Ground-truth day-of-infection is known here, so anchor the root in the
  # naive population, using one of three schemes for comparison:
  if (root_mode == "d0centroid") {
    # Naive cell closest to the D0 centroid in PCA space -- a stable,
    # method-agnostic root not tied to any one TI method's own noise.
    message("Computing D0 (naive) centroid for root cell selection ...")
    naive_cells <- colnames(seurat)[seurat$cell_type == "naive"]
    n_pcs_root  <- min(20, ncol(seurat_pca))
    naive_pca   <- seurat_pca[naive_cells, seq_len(n_pcs_root), drop = FALSE]
    centroid    <- colMeans(naive_pca)
    d_centroid  <- sqrt(rowSums(sweep(naive_pca, 2, centroid, "-")^2))
    start_cell  <- naive_cells[which.min(d_centroid)]
    message(sprintf("  Root cell (nearest D0/naive centroid in PC1-%d): %s",
                    n_pcs_root, start_cell))
  } else if (root_mode == "cytonaive") {
    # Original approach: most-primitive naive cell by CytoTRACE score.
    message("Computing CytoTRACE score for root cell selection (restricted to naive cells) ...")
    cytotrace_pt <- run_cytotrace(expr, fullgene_expr)
    naive_cells  <- colnames(seurat)[seurat$cell_type == "naive"]
    naive_cells  <- intersect(naive_cells, names(cytotrace_pt))
    # run_cytotrace() inverts: 0 = primitive, 1 = differentiated -- so root = min
    start_cell   <- naive_cells[which.min(cytotrace_pt[naive_cells])]
    message(sprintf("  Root cell (min CytoTRACE score among naive cells): %s", start_cell))
  } else {
    # CytoTRACE's own unrestricted root -- most-primitive cell genome-wide,
    # regardless of true cell_type label.
    message("Computing CytoTRACE score for root cell selection (genome-wide, no naive restriction) ...")
    cytotrace_pt <- run_cytotrace(expr, fullgene_expr)
    # run_cytotrace() inverts: 0 = primitive, 1 = differentiated -- so root = min
    start_cell   <- names(cytotrace_pt)[which.min(cytotrace_pt)]
    message(sprintf("  Root cell (min CytoTRACE score genome-wide): %s (cell_type = %s)",
                    start_cell, as.character(seurat$cell_type[start_cell])))
  }

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
    methods       = c(directed_methods, "CytoTRACE")
  )

  # =============================================================================
  # 5. SAVE
  # =============================================================================
  out <- file.path(out_dir, "ti_pseudotime_gse131847.csv")
  write.csv(ti_results, out, row.names = TRUE)
  message("\nSaved: ", out)

  root_out <- file.path(out_dir, "root_cell_gse131847.txt")
  writeLines(start_cell, root_out)
  message("Saved: ", root_out)

} else {
  # ===========================================================================
  # 3+4. dualcyto -- run all 7 directed methods twice, once per CytoTRACE root
  # ===========================================================================
  message("Computing CytoTRACE score (used for both roots + its own column) ...")
  cytotrace_pt <- run_cytotrace(expr, fullgene_expr)

  naive_cells      <- colnames(seurat)[seurat$cell_type == "naive"]
  naive_cells      <- intersect(naive_cells, names(cytotrace_pt))
  root_cytonaive   <- naive_cells[which.min(cytotrace_pt[naive_cells])]
  root_cytoglobal  <- names(cytotrace_pt)[which.min(cytotrace_pt)]
  message(sprintf("  Root (cytonaive):  %s", root_cytonaive))
  message(sprintf("  Root (cytoglobal): %s (cell_type = %s)",
                  root_cytoglobal, as.character(seurat$cell_type[root_cytoglobal])))

  message("\nRunning 7 directed TI methods rooted at cytonaive ...")
  set.seed(42)
  ti_cytonaive <- run_all_ti_methods(
    expr = expr, fullgene_expr = fullgene_expr,
    start_cell = root_cytonaive, precomp_pca = seurat_pca,
    methods = directed_methods
  )
  colnames(ti_cytonaive) <- paste0(colnames(ti_cytonaive), "_cytonaive")

  message("\nRunning 7 directed TI methods rooted at cytoglobal ...")
  set.seed(42)
  ti_cytoglobal <- run_all_ti_methods(
    expr = expr, fullgene_expr = fullgene_expr,
    start_cell = root_cytoglobal, precomp_pca = seurat_pca,
    methods = directed_methods
  )
  colnames(ti_cytoglobal) <- paste0(colnames(ti_cytoglobal), "_cytoglobal")

  ti_results <- cbind(ti_cytonaive[colnames(seurat), , drop = FALSE],
                      ti_cytoglobal[colnames(seurat), , drop = FALSE],
                      CytoTRACE = cytotrace_pt[colnames(seurat)])
  rownames(ti_results) <- colnames(seurat)

  out <- file.path(out_dir, "ti_pseudotime_gse131847.csv")
  write.csv(ti_results, out, row.names = TRUE)
  message("\nSaved: ", out)

  root_out <- file.path(out_dir, "root_cells_gse131847.csv")
  write.csv(data.frame(root_type = c("cytonaive", "cytoglobal"),
                       cell      = c(root_cytonaive, root_cytoglobal)),
            root_out, row.names = FALSE)
  message("Saved: ", root_out)
}

message("Done.")
