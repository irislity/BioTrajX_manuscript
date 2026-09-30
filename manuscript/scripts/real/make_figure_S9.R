# =============================================================================
# make_figure_S9.R
#
# Supplementary figure panels for:
#   S9 — CD8+ T cell exhaustion (linear trajectory)
#
# Panel structure:
#   S9_a_umap.pdf          GT UMAP + per-method pseudotime UMAPs (Spectral)
#                          (methods ordered by DOE, high to low)
#   S9_b_module_trends.pdf Module score trends, faceted by TI method
#                          (methods ordered by DOE, high to low)
#   S9_c_doe_heatmap.pdf   BioTrajX plot(res, type = "heatmap")
#   S9_d_gam_tf_trends.pdf GAM fits of Expr ~ Pseudotime per Method x TF
#                          (early: TCF7, LEF1, FOXP1; late: TOX, BATF, EOMES;
#                          methods ordered by DOE, high to low)
#   S9_e_violin_celltype.pdf Violin of Pseudotime by cell type, faceted by
#                          method (ordered by DOE, high to low)
#
# Prerequisites:
#   data/cd8t.rds
#   manuscript/results/linear_cd8t/ti_pseudotime_cd8t.csv  (from run_ti_cd8t.R)
#
# Usage:
#   Rscript manuscript/scripts/real/make_figure_S9.R [marker_mode]
#     marker_mode  msigdb (default) | slimr_pctit
#                  -- msigdb: KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP (early) /
#                     JIANG_MELANOMA_TRM2_CD8 (terminal), top_n = 30
#                  -- slimr_pctit: SlimR::Markers_list_PCTIT "CD8+ Tn" (early) /
#                     "CD8+ GZMK+ Tex" (terminal), top_n = NULL (no truncation)
#                  Output goes to figures/real/S9 (msigdb) or
#                  figures/real/S9_slimr_pctit.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(ggh4x)
})

repo_root <- here::here()

library(BioTrajX)

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

args        <- commandArgs(trailingOnly = TRUE)
marker_mode <- if (length(args) >= 1) args[1] else "msigdb"
stopifnot(marker_mode %in% c("msigdb", "slimr_pctit"))
message("marker_mode = ", marker_mode)

out_S9 <- file.path(repo_root, "manuscript", "figures", "real",
                    if (marker_mode == "msigdb") "S9" else paste0("S9_", marker_mode))
dir.create(out_S9, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SHARED HELPERS
# =============================================================================

#' UMAP coloured by pseudotime (Spectral palette, S2 style)
plot_pt_umap <- function(df, title = NULL) {
  ggplot(df, aes(x = UMAP1, y = UMAP2, colour = Pseudotime)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_distiller(palette = "Spectral", direction = -1,
                           na.value = "grey80", name = "Pseudotime") +
    theme_classic(base_size = 10) +
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = 9, face = "bold")) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2")
}

#' UMAP coloured by discrete cell type
plot_gt_umap <- function(umap_df, label_col, title = "Cell type") {
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = .data[[label_col]])) +
    geom_point(size = 0.5, alpha = 0.7) +
    theme_classic(base_size = 10) +
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = 9, face = "bold")) +
    labs(title = title, colour = NULL, x = "UMAP 1", y = "UMAP 2")
}

#' Combine GT + per-method UMAP plots with patchwork and save
save_umap_grid <- function(gt_plot, method_plots, ncol = 3,
                            out_path, width = 14, height = NULL) {
  nrow <- ceiling((length(method_plots) + 1) / ncol)
  if (is.null(height)) height <- nrow * 2.8 + 0.5
  if (has_patchwork) {
    combined <- Reduce(`+`, c(list(gt_plot), method_plots)) +
      plot_layout(ncol = ncol, guides = "collect") &
      theme(legend.position = "right")
    ggsave(out_path, combined, width = width, height = height)
  } else {
    message("  (patchwork unavailable — saving GT UMAP only)")
    ggsave(out_path, gt_plot, width = 7, height = 5)
  }
}

#' Save a base-R plotting call to PDF
save_base_pdf <- function(path, width, height, plot_fn) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  tryCatch(plot_fn(), error = function(e) message("  plot error: ", e$message))
}

#' Methods with a DOE score, ordered high to low
.doe_order <- function(res, methods) {
  if (is.null(res)) return(methods)
  cs <- res$comparison_summary
  ord <- cs$trajectory[order(-cs$DOE_score)]
  ord[ord %in% methods]
}

# =============================================================================
# FIGURE S9 — CD8T LINEAR TRAJECTORY
# =============================================================================
message("\n", strrep("=", 70))
message("FIGURE S9 — CD8T linear trajectory")
message(strrep("=", 70))

# ── Load Seurat + pseudotimes ─────────────────────────────────────────────────
message("Loading cd8t.rds ...")
obj_s9 <- tryCatch(
  readRDS(file.path(repo_root, "data", "cd8t.rds")),
  error = function(e) { message("  ERROR: ", e$message); NULL }
)

csv_s9 <- file.path(repo_root, "manuscript", "results", "linear_cd8t", "ti_pseudotime_cd8t.csv")
if (!file.exists(csv_s9))
  stop("ti_pseudotime_cd8t.csv not found. Run run_ti_cd8t.R first.")
ti_df_s9 <- read.csv(csv_s9, row.names = 1, check.names = FALSE)
message(sprintf("  Pseudotimes: %d cells × %d methods", nrow(ti_df_s9), ncol(ti_df_s9)))

# ── Marker genes ──────────────────────────────────────────────────────────────
if (marker_mode == "msigdb") {
  # Naive pole:    KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP (C7)
  #   Curated Kaech naive vs day-8 effector signature; yields highest DOE
  #   scores (top DOE 0.689) and strongest D_early recovery (0.62-0.68).
  # Terminal pole: JIANG_MELANOMA_TRM2_CD8 (C2)
  message("  Fetching MSigDB markers for S9 ...")
  ms_s9 <- get_markers_msigdb(
    early      = "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP",
    terminal   = "JIANG_MELANOMA_TRM2_CD8",
    collection = NULL,
    species    = "Homo sapiens"
  )
  ms_s9 <- filter_markers(ms_s9, obj_s9, top_n = 30, min_detection = 0.10)
} else {
  # SlimR curated pan-cancer T-cell (PCTIT) atlas signatures.
  message("  Fetching SlimR PCTIT markers for S9 ...")
  pctit_s9 <- SlimR::Markers_list_PCTIT
  ms_s9 <- BioTrajX:::.marker_set(
    early    = pctit_s9[["CD8+ Tn"]]$Markers,
    terminal = pctit_s9[["CD8+ GZMK+ Tex"]]$Markers,
    source   = "SlimR_PCTIT",
    metadata = list(early_set = "CD8+ Tn", terminal_set = "CD8+ GZMK+ Tex")
  )
  ms_s9 <- filter_markers(ms_s9, obj_s9, top_n = NULL, min_detection = 0.10)
}
early_genes_s9 <- ms_s9$early
term_genes_s9  <- ms_s9$terminal
message(sprintf("  early markers: %d, terminal markers: %d",
                length(early_genes_s9), length(term_genes_s9)))

# ── Compute multi-DOE (keep result object for heatmap plot) ───────────────────
message("  Computing multi-DOE for CD8T ...")
pt_list_s9 <- as.list(ti_df_s9)
pt_list_s9 <- pt_list_s9[sapply(pt_list_s9, function(x) sum(!is.na(x)) > 10)]

res_s9 <- tryCatch(
  compute_multi_doe_linear(
    expr_or_seurat   = obj_s9,
    pseudotime_list  = pt_list_s9,
    early_markers    = early_genes_s9,
    terminal_markers = term_genes_s9,
    E_method         = "gmm",
    plot_E           = FALSE,
    verbose          = TRUE
  ),
  error = function(e) { message("  ERROR: ", e$message); NULL }
)
if (!is.null(res_s9))
  write.csv(res_s9$comparison_summary,
            file.path(repo_root, "manuscript", "results", "linear_cd8t",
                      if (marker_mode == "msigdb") "doe_scores_cd8t.csv"
                      else sprintf("doe_scores_cd8t_%s.csv", marker_mode)),
            row.names = FALSE)

# ── Panel a: UMAP grid ────────────────────────────────────────────────────────
message("\n[S9-a] UMAP grid")
tryCatch({
  if (is.null(obj_s9)) stop("Seurat object not loaded")
  umap_df <- as.data.frame(Embeddings(obj_s9, "umap"))
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$CellType  <- obj_s9$functional.cluster[rownames(umap_df)]
  celltype_order_s9a <- c("CD8.NaiveLike", "CD8.CM", "CD8.EM", "CD8.TPEX", "CD8.TEX")
  umap_df$CellType  <- factor(umap_df$CellType,
                              levels = intersect(celltype_order_s9a, unique(umap_df$CellType)))

  p_gt <- plot_gt_umap(umap_df, "CellType", title = "Cell type") +
    theme(plot.title = element_text(size = 16, face = "bold"))

  method_order_s9 <- .doe_order(res_s9, colnames(ti_df_s9))
  shared <- intersect(rownames(umap_df), rownames(ti_df_s9))

  # Root cell: identical across all methods — the CD8.NaiveLike centroid
  # cell (see run_ti_cd8t.R), passed as start_cell to every TI method.
  # Persisted to root_cell_cd8t.txt since it can no longer be read off any
  # one method's pseudotime column (unlike CytoTRACE, it isn't a raw score).
  root_cell_path_s9 <- file.path(repo_root, "manuscript", "results", "linear_cd8t",
                                 "root_cell_cd8t.txt")
  root_cell_s9 <- if (file.exists(root_cell_path_s9))
    readLines(root_cell_path_s9, n = 1) else NA_character_
  root_cell_s9 <- if (!is.na(root_cell_s9) && root_cell_s9 %in% shared)
    root_cell_s9 else NA_character_
  root_coord_s9 <- if (!is.na(root_cell_s9))
    umap_df[root_cell_s9, c("UMAP1", "UMAP2")] else NULL

  method_plots <- lapply(method_order_s9, function(m) {
    doe_val <- if (!is.null(res_s9))
      res_s9$comparison_summary$DOE_score[res_s9$comparison_summary$trajectory == m]
    else NA
    p <- plot_pt_umap(
      data.frame(UMAP1      = umap_df[shared, "UMAP1"],
                 UMAP2      = umap_df[shared, "UMAP2"],
                 Pseudotime = ti_df_s9[shared, m]),
      title = m
    )
    if (!is.null(root_coord_s9))
      p <- p + geom_point(data = root_coord_s9, aes(x = UMAP1, y = UMAP2),
                          colour = "black", shape = 17, size = 4, inherit.aes = FALSE)
    p + labs(subtitle = if (!is.na(doe_val) && length(doe_val) == 1)
               sprintf("DOE = %.3f", doe_val) else NULL) +
      theme(plot.subtitle = element_text(size = 7, colour = "grey40"))
  })

  save_umap_grid(p_gt, method_plots, ncol = 3,
                 out_path = file.path(out_S9, "S9_a_umap.pdf"), width = 14)
  message("  Saved S9_a_umap.pdf")
}, error = function(e) message("  SKIPPED S9-a: ", e$message))

# ── Panel b: Module score trends per method ────────────────────────────────────
message("\n[S9-b] Module score trends")
tryCatch({
  if (is.null(obj_s9)) stop("Seurat object not loaded")
  obj_s9 <- AddModuleScore(obj_s9, features = list(early_genes_s9), name = "early_mod")
  obj_s9 <- AddModuleScore(obj_s9, features = list(term_genes_s9),  name = "term_mod")

  scores <- list(
    Early    = setNames(obj_s9$early_mod1, colnames(obj_s9)),
    Terminal = setNames(obj_s9$term_mod1,  colnames(obj_s9))
  )

  cells_all <- Reduce(intersect, lapply(scores, names))
  cells_all <- intersect(cells_all, rownames(ti_df_s9))

  trend_long <- do.call(rbind, lapply(colnames(ti_df_s9), function(m) {
    pt   <- ti_df_s9[cells_all, m]
    keep <- !is.na(pt)
    do.call(rbind, lapply(names(scores), function(mod) {
      data.frame(Pseudotime = pt[keep], Score = scores[[mod]][cells_all][keep],
                 Module = mod, Method = m, stringsAsFactors = FALSE)
    }))
  }))
  trend_long$Module <- factor(trend_long$Module, levels = c("Early", "Terminal"))
  trend_long$Method <- factor(trend_long$Method, levels = .doe_order(res_s9, colnames(ti_df_s9)))

  p_b9 <- ggplot(trend_long, aes(x = Pseudotime, y = Score, colour = Module)) +
    geom_point(size = 0.2, alpha = 0.15) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.9, span = 0.4) +
    scale_colour_manual(values = c(Early = "#4E9AF1", Terminal = "#2D7A3E")) +
    facet_wrap(~ Method, ncol = 4, scales = "free_x", axes = "all") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 15),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16),
          panel.grid.minor = element_blank()) +
    labs(x = "Pseudotime [0, 1]", y = "Module score", colour = NULL)

  n_m <- length(unique(trend_long$Method))
  ggsave(file.path(out_S9, "S9_b_module_trends.pdf"), p_b9,
         width = 14, height = ceiling(n_m / 4) * 3 + 1)
  message("  Saved S9_b_module_trends.pdf")
}, error = function(e) message("  SKIPPED S9-b: ", e$message))

# ── Panel c: BioTrajX DOE heatmap (style/font size matching S11 panel e) ─────
message("\n[S9-c] DOE heatmap (BioTrajX)")
tryCatch({
  if (is.null(res_s9)) stop("res_s9 is NULL")
  if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
    sysfonts::font_add("Arial",
      regular    = file.path(repo_root, "manuscript", "fonts", "Arial.ttf"),
      bold       = file.path(repo_root, "manuscript", "fonts", "Arial Bold.ttf"),
      italic     = file.path(repo_root, "manuscript", "fonts", "Arial Italic.ttf"),
      bolditalic = file.path(repo_root, "manuscript", "fonts", "Arial Bold Italic.ttf"))
    showtext::showtext_auto()
    showtext::showtext_opts(dpi = 300)
  }
  p_c9 <- plot(res_s9, type = "heatmap") +
    labs(title = NULL) +
    theme(panel.grid = element_blank(),
          text        = element_text(family = "Arial"),
          axis.text   = element_text(size = 15, family = "Arial"),
          axis.title  = element_text(size = 16, family = "Arial"),
          legend.text  = element_text(size = 14, family = "Arial"),
          legend.title = element_text(size = 15, family = "Arial"))
  p_c9$layers[[2]]$aes_params$size <- 5.1
  p_c9$layers[[2]]$aes_params$family <- "Arial"
  ggsave(file.path(out_S9, "S9_c_doe_heatmap.pdf"), p_c9, width = 10, height = 5)
  if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto(FALSE)
  message("  Saved S9_c_doe_heatmap.pdf")
}, error = function(e) message("  SKIPPED S9-c: ", e$message))

# =============================================================================
# Panel d: GAM fits of Expr ~ Pseudotime, per Method x TF
#
# Early/naive-like: TCF7, LEF1, FOXP1.  Late/exhaustion-associated: TOX,
# BATF, EOMES. Methods ordered by DOE score, high to low.
# =============================================================================
message("\n[S9-d] GAM fits x Pseudotime, per Method x TF")
tryCatch({
  if (is.null(obj_s9)) stop("Seurat object not loaded")
  if (is.null(res_s9)) stop("res_s9 is NULL")

  naive_genes    <- c("TCF7", "LEF1", "FOXP1")
  terminal_genes <- c("TOX", "BATF", "EOMES")
  tf_genes  <- c(naive_genes, terminal_genes)
  expr_mat  <- GetAssayData(obj_s9, layer = "data")
  genes_present <- intersect(tf_genes, rownames(expr_mat))
  if (length(genes_present) == 0) stop("none of the requested TFs found in obj_s9")

  method_order_s9d <- .doe_order(res_s9, colnames(ti_df_s9))

  df_d <- do.call(rbind, lapply(method_order_s9d, function(m) {
    cells <- intersect(rownames(ti_df_s9), colnames(obj_s9))
    pt    <- ti_df_s9[cells, m]; keep <- !is.na(pt)
    cells <- cells[keep]; pt <- pt[keep]
    do.call(rbind, lapply(genes_present, function(g) {
      data.frame(Method = m, Gene = g, Pseudotime = pt,
                 Expr = as.numeric(expr_mat[g, cells]),
                 GeneGroup = if (g %in% naive_genes) "Early" else "Terminal",
                 stringsAsFactors = FALSE)
    }))
  }))
  df_d$Method    <- factor(df_d$Method,    levels = method_order_s9d)
  df_d$Gene      <- factor(df_d$Gene,      levels = genes_present)
  df_d$GeneGroup <- factor(df_d$GeneGroup, levels = c("Early", "Terminal"))

  p_d9 <- ggplot(df_d, aes(x = Pseudotime, y = Expr)) +
    geom_point(size = 0.15, alpha = 0.08, colour = "grey50") +
    geom_smooth(aes(colour = GeneGroup, fill = GeneGroup),
                method = "gam", formula = y ~ s(x, bs = "cs"),
                se = TRUE, linewidth = 0.9) +
    scale_colour_manual(values = c(Early = "#4E9AF1", Terminal = "#2D7A3E"), guide = "none") +
    scale_fill_manual(values = c(Early = "#4E9AF1", Terminal = "#2D7A3E"), guide = "none") +
    ggh4x::facet_grid2(rows = vars(Method), cols = vars(GeneGroup, Gene),
                       scales = "free_y", independent = "y", axes = "all",
                       strip = ggh4x::strip_nested()) +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold", size = 11),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          panel.border = element_rect(colour = "grey20", fill = NA),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "Pseudotime [0, 1]", y = "Expression", title = NULL)

  ggsave(file.path(out_S9, "S9_d_gam_tf_trends.pdf"), p_d9,
         width = length(genes_present) * 2.2 + 1,
         height = length(method_order_s9d) * 1.8 + 1)
  message("  Saved S9_d_gam_tf_trends.pdf")
}, error = function(e) message("  SKIPPED S9-d: ", e$message))

# =============================================================================
# Panel e: violin plot of Pseudotime by cell type, faceted by method
# =============================================================================
message("\n[S9-e] Pseudotime by cell type, violin")
tryCatch({
  if (is.null(obj_s9)) stop("Seurat object not loaded")
  if (is.null(res_s9)) stop("res_s9 is NULL")

  cell_type <- obj_s9$functional.cluster
  method_order_s9e <- .doe_order(res_s9, colnames(ti_df_s9))
  celltype_order_s9e <- c("CD8.NaiveLike", "CD8.CM", "CD8.EM", "CD8.TPEX", "CD8.TEX")

  df_e <- do.call(rbind, lapply(method_order_s9e, function(m) {
    cells <- intersect(rownames(ti_df_s9), colnames(obj_s9))
    pt    <- ti_df_s9[cells, m]; keep <- !is.na(pt)
    cells <- cells[keep]
    data.frame(Method = m, CellType = cell_type[cells], Pseudotime = pt[keep],
               stringsAsFactors = FALSE)
  }))
  df_e$Method   <- factor(df_e$Method,   levels = method_order_s9e)
  df_e$CellType <- factor(df_e$CellType,
                          levels = intersect(celltype_order_s9e, unique(df_e$CellType)))

  doe_labels_e <- data.frame(
    Method = method_order_s9e,
    label  = sprintf("DOE = %.3f",
                     res_s9$comparison_summary$DOE_score[
                       match(method_order_s9e, res_s9$comparison_summary$trajectory)]),
    stringsAsFactors = FALSE
  )
  doe_labels_e$Method <- factor(doe_labels_e$Method, levels = method_order_s9e)

  p_e9 <- ggplot(df_e, aes(x = CellType, y = Pseudotime)) +
    geom_violin(aes(fill = CellType), scale = "width", colour = NA,
                alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.size = 0.4, outlier.alpha = 0.3,
                 fill = "white", alpha = 0.7, linewidth = 0.3) +
    geom_text(data = doe_labels_e, aes(label = label), x = -Inf, y = Inf,
              hjust = -0.05, vjust = 1.3, size = 3.2, colour = "grey20",
              inherit.aes = FALSE) +
    facet_wrap(~ Method, ncol = 4, axes = "all") +
    scale_fill_hue(guide = "none") +
    theme_classic(base_size = 14) +
    theme(strip.text = element_text(face = "bold", size = 14),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    labs(x = "Cell type", y = "Pseudotime", title = NULL)

  n_m <- length(method_order_s9e)
  ggsave(file.path(out_S9, "S9_e_violin_celltype.pdf"), p_e9,
         width = 14, height = ceiling(n_m / 4) * 3.4 + 1)
  message("  Saved S9_e_violin_celltype.pdf")
}, error = function(e) message("  SKIPPED S9-e: ", e$message))

message("\nFigure S9 complete. Output: ", out_S9)
