# =============================================================================
# pipeline.R
#
# General-purpose pipeline: run all TI methods on any dataset, then evaluate
# each with BioTrajX DOE metrics (D, O, E).
#
# ── QUICK START ───────────────────────────────────────────────────────────────
# Edit the CONFIG block below, then:
#   Rscript manuscript/pipeline.R
#
# ── PROGRAMMATIC USE ──────────────────────────────────────────────────────────
# source("manuscript/pipeline.R")
# run_pipeline(
#   expr           = my_expr_matrix,   # genes × cells log-norm
#   naive_markers  = c("TCF7", ...),
#   term_markers   = c("GZMB", ...),
#   out_dir        = "results/my_dataset"
# )
#
# ── INPUTS ────────────────────────────────────────────────────────────────────
# expr          genes × cells log-normalised matrix (matrix / dgCMatrix / path to .rds)
# counts        genes × cells raw count matrix (optional; enables CytoTRACE2)
# naive_markers character vector of naive/progenitor marker genes
# term_markers  character vector of terminal/differentiated marker genes
# start_cell    name of the root cell (NULL = lowest pseudotime cell in first cluster)
# species       "human" or "mouse" (used by CytoTRACE2)
# methods       TI methods to attempt; any subset of the list below
# E_method      BioTrajX endpoint metric: "gmm", "clusters", or "combined"
#
# ── OUTPUTS (written to out_dir/) ─────────────────────────────────────────────
# ti_pseudotimes.csv      pseudotime [0,1] per cell per TI method
# doe_summary.csv         D / O / E / DOE_score per method, ranked
# multi_doe_bar.pdf       grouped bar chart
# multi_doe_heatmap.pdf   score heatmap
# =============================================================================

# ── CONFIG ────────────────────────────────────────────────────────────────────
# Edit this block when running via Rscript; ignored when sourced.

CONFIG <- list(
  expr           = "data/stemcell_biotrajx/expr.rds",
  counts         = NULL,         # path to counts .rds, or NULL to skip CytoTRACE2
  naive_markers  = NULL,         # NULL = load from markers_rds
  term_markers   = NULL,
  markers_rds    = "data/stemcell_biotrajx/marker_genes.rds",
  start_cell     = NULL,         # NULL = auto (cell with lowest mean expression rank)
  species        = "human",
  out_dir        = "manuscript/plots/pipeline_out",
  dataset_name   = "stemcell",
  E_method       = "gmm",
  methods        = c("Slingshot", "CytoTRACE", "CytoTRACE2",
                     "DPT", "SCORPIUS",
                     "PAGA-DPT", "CellRank")
)

# =============================================================================
# PIPELINE FUNCTION
# =============================================================================

#' Run the full TI → BioTrajX DOE pipeline on any dataset
#'
#' @param expr           genes × cells log-norm matrix, dgCMatrix, or path to .rds
#' @param naive_markers  character vector of naive/progenitor marker genes
#' @param term_markers   character vector of terminal marker genes
#' @param out_dir        directory to write all outputs (created if absent)
#' @param counts         genes × cells raw count matrix or path to .rds (for CytoTRACE2)
#' @param start_cell     root cell name; NULL = cell with minimum summed expression
#' @param species        "human" or "mouse" (CytoTRACE2)
#' @param methods        TI method names to attempt
#' @param E_method       BioTrajX E metric: "gmm", "clusters", or "combined"
#' @param dataset_name   label used in plot titles
#' @param verbose        print progress messages
#' @return list with $ti (pseudotime data.frame) and $doe (multi_doe_results)
run_pipeline <- function(
    expr,
    naive_markers,
    term_markers,
    out_dir,
    counts         = NULL,
    start_cell     = NULL,
    species        = "human",
    methods        = c("Slingshot", "CytoTRACE", "CytoTRACE2",
                       "DPT", "SCORPIUS",
                       "PAGA-DPT", "CellRank"),
    E_method       = "gmm",
    dataset_name   = "dataset",
    verbose        = TRUE
) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Load expression matrix ──────────────────────────────────────────────────
  if (is.character(expr)) {
    if (verbose) message("Loading expr from ", expr)
    expr <- readRDS(expr)
  }
  expr <- as.matrix(expr)
  cells <- colnames(expr)
  if (verbose) message(sprintf("Expression: %d genes × %d cells", nrow(expr), ncol(expr)))

  # ── Load counts (optional) ──────────────────────────────────────────────────
  if (is.character(counts)) {
    if (verbose) message("Loading counts from ", counts)
    counts <- as.matrix(readRDS(counts))
  } else if (!is.null(counts)) {
    counts <- as.matrix(counts)
  }

  # ── Auto start_cell — cell with the lowest total expression (most naive) ────
  if (is.null(start_cell)) {
    col_sums  <- colSums(expr)
    start_cell <- cells[which.min(col_sums)]
    if (verbose) message(sprintf("Auto start_cell: %s (colSum = %.1f)",
                                 start_cell, col_sums[start_cell]))
  }

  # ── Step 1: Run all TI methods ──────────────────────────────────────────────
  if (verbose) message("\n── Step 1: TI methods ───────────────────────────────────────────────")
  ti_df <- run_all_ti_methods(
    expr        = expr,
    start_cell  = start_cell,
    counts_expr = counts,
    methods     = methods
  )

  write.csv(ti_df, file.path(out_dir, "ti_pseudotimes.csv"), row.names = TRUE)
  if (verbose) message("Saved ti_pseudotimes.csv")

  # ── Drop methods that returned all-NA ──────────────────────────────────────
  has_data <- sapply(ti_df, function(pt) sum(!is.na(pt)) > 10)
  if (verbose) {
    message(sprintf("Valid methods:   %s", paste(names(ti_df)[has_data],  collapse = ", ")))
    message(sprintf("Skipped (all NA): %s", paste(names(ti_df)[!has_data], collapse = ", ")))
  }
  pt_list <- as.list(ti_df[, has_data, drop = FALSE])

  if (length(pt_list) == 0)
    stop("No TI methods produced valid pseudotime. Cannot run BioTrajX.")

  # ── Step 2: BioTrajX multi-DOE ─────────────────────────────────────────────
  if (verbose) message("\n── Step 2: BioTrajX DOE ─────────────────────────────────────────────")
  res_multi <- compute_multi_doe_linear(
    expr_or_seurat   = expr,
    pseudotime_list  = pt_list,
    naive_markers    = naive_markers,
    terminal_markers = term_markers,
    E_method         = E_method,
    plot_E           = FALSE,
    verbose          = verbose
  )

  # ── Summary table ───────────────────────────────────────────────────────────
  summary_df <- res_multi$comparison_summary[, c("trajectory", "D_naive", "D_term",
                                                   "O", "E_naive", "E_term",
                                                   "E_comp", "DOE_score")]
  write.csv(summary_df, file.path(out_dir, "doe_summary.csv"), row.names = FALSE)

  if (verbose) {
    message("\nResults:")
    print(summary_df[, c("trajectory","D_naive","D_term","O","E_naive","E_term","DOE_score")])
    message(sprintf("\nBest: %s  (DOE = %.3f)", res_multi$best_trajectory,
                    summary_df$DOE_score[summary_df$trajectory == res_multi$best_trajectory]))
  }

  # ── Plots ───────────────────────────────────────────────────────────────────
  n <- length(pt_list)

  p_bar <- plot(res_multi, type = "bar") +
    ggplot2::labs(title = sprintf("%s — DOE metrics by TI method", dataset_name))
  ggplot2::ggsave(file.path(out_dir, "multi_doe_bar.pdf"),
                  p_bar, width = max(8, n * 1.2), height = 5)

  p_heat <- plot(res_multi, type = "heatmap") +
    ggplot2::labs(title = sprintf("%s — DOE heatmap", dataset_name))
  ggplot2::ggsave(file.path(out_dir, "multi_doe_heatmap.pdf"),
                  p_heat, width = 8, height = max(4, n * 0.6))

  if (verbose) message(sprintf("Saved plots to %s", out_dir))

  invisible(list(ti = ti_df, doe = res_multi))
}

# =============================================================================
# ENTRY POINT  (Rscript manuscript/pipeline.R)
# =============================================================================

if (!interactive() && sys.nframe() == 0L) {

  suppressPackageStartupMessages({
    library(ggplot2)
    library(reshape2)
  })

  # ── Locate repo root ──────────────────────────────────────────────────────
  repo_root <- tryCatch(
    normalizePath(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "..")),
    error = function(e) getwd()
  )

  # ── Load BioTrajX ─────────────────────────────────────────────────────────
  library(BioTrajX)

  # ── Source TI methods ─────────────────────────────────────────────────────
  source(file.path(repo_root, "manuscript", "run_ti_methods.R"))

  # ── Resolve CONFIG paths relative to repo root ────────────────────────────
  .abs <- function(p) if (!is.null(p)) file.path(repo_root, p) else NULL

  # ── Load markers ─────────────────────────────────────────────────────────
  naive <- CONFIG$naive_markers
  term  <- CONFIG$term_markers
  if (is.null(naive) || is.null(term)) {
    mrk   <- readRDS(.abs(CONFIG$markers_rds))
    naive <- naive %||% mrk$naive_markers
    term  <- term  %||% mrk$term_markers
  }

  run_pipeline(
    expr          = .abs(CONFIG$expr),
    counts        = .abs(CONFIG$counts),
    naive_markers = naive,
    term_markers  = term,
    start_cell    = CONFIG$start_cell,
    species       = CONFIG$species,
    out_dir       = .abs(CONFIG$out_dir),
    dataset_name  = CONFIG$dataset_name,
    E_method      = CONFIG$E_method,
    methods       = CONFIG$methods
  )
}
