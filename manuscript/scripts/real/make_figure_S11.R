# =============================================================================
# make_figure_S11.R
#
# Demonstrates BioTrajX as a root selection aid for Monocle3.
# Candidate root cells are selected first (before any Monocle3 run), then
# Monocle3 is run independently from each root — exactly as a user would do
# when exploring different starting points. BioTrajX DOE scores are used to
# automatically identify the biologically correct root.
#
# Candidate roots:
#   CD8.NaiveLike — root 1: same CytoTRACE selection as S8 (which.min on
#                            inverted run_cytotrace() score across full dataset)
#   CD8.NaiveLike — root 2: closest to PCA centroid of cluster
#   CD8.CM / CD8.EM / CD8.TPEX / CD8.TEX — PCA centroid cell per cluster
#
# Panels:
#   S11_a_umap.pdf        UMAP coloured by pseudotime for each candidate root
#   S11_b_gene_trends.pdf Naive/exhausted module score trends vs pseudotime per root
#   S11_c_doe_bar.pdf     DOE score bar chart (NaiveLike blue, other orange)
#   S11_d_radar.pdf       BioTrajX radar plot (D / O / E) across all roots
#   S11_e_doe_heatmap.pdf BioTrajX DOE heatmap (D / O / E / DOE) across all roots
#
# Prerequisites:
#   data/cd8t.rds
#
# Usage:
#   Rscript manuscript/scripts/real/make_figure_S11.R [marker_mode]
#     marker_mode  msigdb (default) | slimr_pctit
#                  -- msigdb: KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP (early) /
#                     JIANG_MELANOMA_TRM2_CD8 (terminal), top_n = 30
#                  -- slimr_pctit: SlimR::Markers_list_PCTIT "CD8+ Tn" (early) /
#                     "CD8+ GZMK+ Tex" (terminal), top_n = NULL (no truncation)
#                  Output goes to figures/real/S11 (msigdb) or
#                  figures/real/S11_slimr_pctit.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
})

repo_root <- here::here()
library(BioTrajX)

# Source run_monocle3() and helpers
source(file.path(repo_root, "manuscript", "scripts", "real", "run_ti_methods.R"))

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

args        <- commandArgs(trailingOnly = TRUE)
marker_mode <- if (length(args) >= 1) args[1] else "msigdb"
stopifnot(marker_mode %in% c("msigdb", "slimr_pctit"))
message("marker_mode = ", marker_mode)

out_S11 <- file.path(repo_root, "manuscript", "figures", "real",
                     if (marker_mode == "msigdb") "S11" else paste0("S11_", marker_mode))
dir.create(out_S11, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA AND MARKERS
# =============================================================================
message("Loading cd8t.rds ...")
seurat <- readRDS(file.path(repo_root, "data", "cd8t.rds"))

# Subset to top 2000 HVGs — identical to S8 (run_ti_cd8t.R)
message("Selecting top 2000 HVGs ...")
seurat <- FindVariableFeatures(seurat, nfeatures = 2000, verbose = FALSE)
hvg    <- VariableFeatures(seurat)

expr_mat      <- as.matrix(GetAssayData(seurat, layer = "data")[hvg, ])
fullgene_expr <- as.matrix(GetAssayData(seurat, layer = "data"))

if (marker_mode == "msigdb") {
  message("Fetching MSigDB markers (same as S8) ...")
  ms <- get_markers_msigdb(
    early      = "KAECH_NAIVE_VS_DAY8_EFF_CD8_TCELL_UP",
    terminal   = "JIANG_MELANOMA_TRM2_CD8",
    collection = NULL,
    species    = "Homo sapiens"
  )
  ms <- filter_markers(ms, seurat, top_n = 30, min_detection = 0.10)
} else {
  message("Fetching SlimR PCTIT markers for S11 ...")
  pctit_s11 <- SlimR::Markers_list_PCTIT
  ms <- BioTrajX:::.marker_set(
    early    = pctit_s11[["CD8+ Tn"]]$Markers,
    terminal = pctit_s11[["CD8+ GZMK+ Tex"]]$Markers,
    source   = "SlimR_PCTIT",
    metadata = list(early_set = "CD8+ Tn", terminal_set = "CD8+ GZMK+ Tex")
  )
  ms <- filter_markers(ms, seurat, top_n = NULL, min_detection = 0.10)
}
naive_genes  <- ms$early
tex_genes    <- ms$terminal
message(sprintf("  Markers: %d naive, %d exhausted", length(naive_genes), length(tex_genes)))

# Compute naive module score now so it can be used for root selection below
message("Computing naive module score ...")
seurat <- AddModuleScore(seurat, features = list(naive_genes), name = "naive_mod")

# =============================================================================
# 2. SELECT CANDIDATE ROOT CELLS  (before running Monocle3)
# =============================================================================
message("\nSelecting candidate root cells ...")
meta  <- seurat@meta.data
cells <- colnames(seurat)

# PCA coords for centroid-based selection
pca_coords <- Embeddings(seurat, "pca")[, 1:20]

# Helper: cell closest to cluster PCA centroid
centroid_cell <- function(cluster_name) {
  idx   <- which(meta$functional.cluster == cluster_name)
  ctr   <- colMeans(pca_coords[idx, , drop = FALSE])
  dists <- rowSums(sweep(pca_coords[idx, , drop = FALSE], 2, ctr)^2)
  cells[idx[which.min(dists)]]
}

# NaiveLike root 1: same CytoTRACE selection as S8
# run_cytotrace() inverts the score (0 = primitive, 1 = differentiated) so root = which.min
message("Computing CytoTRACE score for root selection (same as S8) ...")
cytotrace_pt <- run_cytotrace(expr_mat, fullgene_expr)
nl_root1     <- names(which.min(cytotrace_pt))

# NaiveLike root 2: PCA centroid
nl_root2 <- centroid_cell("CD8.NaiveLike")

# NaiveLike root 3: cell with maximum naive module score (global which.max)
nl_root3 <- names(which.max(setNames(seurat$naive_mod1, colnames(seurat))))
message(sprintf("  NaiveLike (naive markers) root: %s  (naive mod score = %.3f)",
                nl_root3, max(seurat$naive_mod1)))

# One centroid cell per other cluster
roots <- c(
  "NaiveLike (CytoTRACE)"     = nl_root1,
  "NaiveLike (centroid)"      = nl_root2,
  "NaiveLike (naive markers)" = nl_root3,
  "CD8.CM"                    = centroid_cell("CD8.CM"),
  "CD8.EM"                    = centroid_cell("CD8.EM"),
  "CD8.TPEX"                  = centroid_cell("CD8.TPEX"),
  "CD8.TEX"                   = centroid_cell("CD8.TEX")
)

message("\nCandidate roots:")
for (nm in names(roots)) message(sprintf("  %-30s %s", nm, roots[nm]))

# Cluster colour palette — ggplot2's default hue scale, but ordered along the
# naive-to-exhausted trajectory (NaiveLike first) rather than alphabetically,
# matching S9 panel A's cell-type ordering.
cluster_levels <- c("CD8.NaiveLike", "CD8.CM", "CD8.EM", "CD8.TPEX", "CD8.TEX")
cluster_levels <- intersect(cluster_levels, unique(meta$functional.cluster))
cluster_pal    <- setNames(scales::hue_pal()(length(cluster_levels)), cluster_levels)

# Map each root to its cluster colour
root_cluster <- c(
  "NaiveLike (CytoTRACE)"     = "CD8.NaiveLike",
  "NaiveLike (centroid)"      = "CD8.NaiveLike",
  "NaiveLike (naive markers)" = "CD8.NaiveLike",
  "CD8.CM"                    = "CD8.CM",
  "CD8.EM"                    = "CD8.EM",
  "CD8.TPEX"                  = "CD8.TPEX",
  "CD8.TEX"                   = "CD8.TEX"
)
root_pal <- setNames(cluster_pal[root_cluster], names(root_cluster))

# =============================================================================
# 3. RUN MONOCLE3 INDEPENDENTLY FROM EACH ROOT
# =============================================================================
message("\nRunning Monocle3 independently from each candidate root ...")

pt_list <- lapply(names(roots), function(nm) {
  rc <- roots[[nm]]
  message(sprintf("\n  [%s] root = %s", nm, rc))
  tryCatch(
    run_monocle3(expr_mat, start_cell = rc),
    error = function(e) {
      message("    failed: ", e$message)
      setNames(rep(NA_real_, length(cells)), cells)
    }
  )
})
names(pt_list) <- names(roots)

# =============================================================================
# 4. COMPUTE MULTI-DOE
# =============================================================================
message("\nComputing BioTrajX DOE for each root ...")
res_s10 <- tryCatch(
  compute_multi_doe_linear(
    expr_or_seurat   = seurat,
    pseudotime_list  = pt_list,
    early_markers    = naive_genes,
    terminal_markers = tex_genes,
    E_method         = "gmm",
    plot_E           = FALSE,
    verbose          = TRUE
  ),
  error = function(e) { message("ERROR: ", e$message); NULL }
)

if (!is.null(res_s10)) {
  message("\n── DOE scores per root ──")
  print(res_s10$comparison_summary[, c("trajectory", "DOE_score")],
        row.names = FALSE, digits = 3)
  write.csv(res_s10$comparison_summary,
            file.path(out_S11, "doe_scores_by_root.csv"), row.names = FALSE)
}

# =============================================================================
# 5. PANEL A — UMAP grid (pseudotime per root)
# =============================================================================
message("\n[S11-a] UMAP grid")
tryCatch({
  umap_coords <- as.data.frame(Embeddings(seurat, "umap"))
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$cluster   <- factor(meta$functional.cluster[match(rownames(umap_coords), cells)],
                                  levels = cluster_levels)

  # Font sizes chosen so that, after this panel is scaled down (~0.79x) inside
  # the assembled Figure S11 (and further inside Figure 1's reuse), the
  # smallest text (legend/subtitle) stays >=7pt, axis titles >=8pt, and panel
  # headings land in [9, 11]pt on the final printed page.
  BASE_A <- 11.63; HEAD_A <- 12.64; SUB_A <- 10.11

  # Ground-truth cluster UMAP
  p_gt <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, colour = cluster)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_manual(values = cluster_pal, name = NULL) +
    theme_classic(base_size = BASE_A) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 16, face = "bold")) +
    labs(title = "Cell type", x = "UMAP 1", y = "UMAP 2")

  # Per-root pseudotime UMAPs; mark the root cell with a star
  pt_plots <- lapply(names(roots), function(nm) {
    pt  <- pt_list[[nm]]
    df  <- data.frame(UMAP1      = umap_coords$UMAP1,
                      UMAP2      = umap_coords$UMAP2,
                      Pseudotime = pt[rownames(umap_coords)])
    rc_coord <- umap_coords[roots[[nm]], , drop = FALSE]
    doe_val  <- if (!is.null(res_s10))
      res_s10$comparison_summary$DOE_score[
        res_s10$comparison_summary$trajectory == nm] else NA

    ggplot(df, aes(x = UMAP1, y = UMAP2, colour = Pseudotime)) +
      geom_point(size = 0.4, alpha = 0.6) +
      geom_point(data = rc_coord, aes(x = UMAP1, y = UMAP2),
                 colour = "black", shape = 17, size = 4, inherit.aes = FALSE) +
      scale_colour_distiller(palette = "Spectral", direction = -1,
                             na.value = "grey80", name = "Pseudotime") +
      theme_classic(base_size = BASE_A) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            legend.position = "right",
            legend.key.height = unit(0.4, "cm"),
            plot.title    = element_text(size = HEAD_A, face = "bold"),
            plot.subtitle = element_text(size = SUB_A, colour = "grey40")) +
      labs(title    = nm,
           subtitle = if (!is.na(doe_val)) sprintf("DOE = %.3f", doe_val) else "",
           x = NULL, y = NULL)
  })

  if (has_patchwork) {
    all_plots <- c(list(p_gt), pt_plots)
    combined  <- Reduce(`+`, all_plots) +
      plot_layout(ncol = 4, guides = "collect") &
      theme(legend.position = "right")
    nrow_val  <- ceiling(length(all_plots) / 4)
    ggsave(file.path(out_S11, "S11_a_umap.pdf"), combined,
           width = 14, height = nrow_val * 2.8 + 0.5)
  } else {
    ggsave(file.path(out_S11, "S11_a_umap.pdf"), p_gt, width = 7, height = 5)
  }
  message("  Saved S11_a_umap.pdf")
}, error = function(e) message("  SKIPPED S11-a: ", e$message))

# =============================================================================
# 6. PANEL B — Gene expression trends (module scores vs pseudotime per root)
# =============================================================================
message("\n[S11-b] Gene expression trends")
tryCatch({
  seurat <- AddModuleScore(seurat, features = list(tex_genes), name = "tex_mod")

  scores <- list(
    Naive     = setNames(seurat$naive_mod1, colnames(seurat)),
    Exhausted = setNames(seurat$tex_mod1,   colnames(seurat))
  )
  cells_sc <- Reduce(intersect, lapply(scores, names))

  trend_long <- do.call(rbind, lapply(names(pt_list), function(nm) {
    pt   <- pt_list[[nm]][cells_sc]
    keep <- !is.na(pt)
    do.call(rbind, lapply(names(scores), function(mod) {
      data.frame(Pseudotime = pt[keep], Score = scores[[mod]][cells_sc][keep],
                 Module = mod, Root = nm, stringsAsFactors = FALSE)
    }))
  }))
  trend_long$Module <- factor(trend_long$Module, levels = c("Naive", "Exhausted"))
  trend_long$Root   <- factor(trend_long$Root, levels = names(roots))

  # Font sizes target: tick/legend text >=7pt, axis titles >=8pt, facet
  # heading (Root name) in [9,11]pt after the ~0.92x scale-down in S11.
  BASE_B <- 9.97; HEAD_B <- 10.84

  p_b <- ggplot(trend_long, aes(x = Pseudotime, y = Score, colour = Module)) +
    geom_point(size = 0.15, alpha = 0.10) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.9, span = 0.4) +
    scale_colour_manual(values = c(Naive = "#56B4E9", Exhausted = "#D55E00")) +
    facet_wrap(~ Root, ncol = 3, scales = "free_x") +
    theme_minimal(base_size = BASE_B) +
    theme(legend.position  = "bottom",
          strip.text       = element_text(face = "bold", size = HEAD_B),
          panel.grid.minor = element_blank()) +
    labs(x = "Pseudotime [0, 1]", y = "Module score", colour = NULL)

  n_roots <- length(unique(trend_long$Root))
  ggsave(file.path(out_S11, "S11_b_gene_trends.pdf"), p_b,
         width = 12, height = ceiling(n_roots / 3) * 3 + 1)
  message("  Saved S11_b_gene_trends.pdf")
}, error = function(e) message("  SKIPPED S11-b: ", e$message))

# =============================================================================
# 7. PANEL C — DOE bar chart (DOE score, NaiveLike blue vs other orange)
# =============================================================================
message("\n[S11-c] DOE bar chart")
tryCatch({
  if (is.null(res_s10)) stop("res_s10 is NULL")
  df_bar <- res_s10$comparison_summary[, c("trajectory", "DOE_score")]
  df_bar$is_naive   <- grepl("NaiveLike", df_bar$trajectory)
  df_bar$trajectory <- factor(df_bar$trajectory,
                               levels = df_bar$trajectory[order(df_bar$DOE_score)])

  # This panel is used near its native size (scale ~1.0) in the assembled
  # figure, so base_size=12 already clears the axis/tick floors; only the
  # plot title needs an explicit cap so it doesn't exceed the 11pt heading
  # ceiling once printed at full size.
  p_bar <- ggplot(df_bar, aes(x = trajectory, y = DOE_score, fill = is_naive)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#56B4E9", "FALSE" = "#D55E00"),
                      labels = c("TRUE" = "CD8.NaiveLike", "FALSE" = "Other cluster"),
                      name   = "Root cluster") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(title = "BioTrajX DOE score by Monocle3 root cell",
         x = "Root candidate", y = "DOE score [0, 1]")

  ggsave(file.path(out_S11, "S11_c_doe_bar.pdf"), p_bar, width = 8, height = 5)
  message("  Saved S11_c_doe_bar.pdf")
}, error = function(e) message("  SKIPPED S11-c: ", e$message))

# =============================================================================
# 8. PANEL D — BioTrajX radar plot
# =============================================================================
message("\n[S11-d] Radar plot")
tryCatch({
  if (is.null(res_s10)) stop("res_s10 is NULL")
  if (!requireNamespace("fmsb", quietly = TRUE)) stop("fmsb package required for radar plots")

  # Replicates plot.multi_doe_results(..., type = "radar") locally (rather
  # than calling the package function) so we can shrink the title via
  # par(cex.main), since fmsb::radarchart's `title` argument otherwise
  # renders at the device default (~13-14pt), which is too large for this
  # panel's [9, 11]pt heading target at its ~1.0x assembled scale. The axis
  # (calcex), vertex (vlcex), and legend cex values below are unchanged from
  # the package defaults, since those already read >=7pt/>=8pt at this scale.
  metrics    <- c("D_early", "D_term", "O", "E_early", "E_term", "DOE_score")
  plot_data  <- res_s10$comparison_summary[, c("trajectory", metrics), drop = FALSE]
  radar_data <- plot_data[, metrics, drop = FALSE]
  radar_data[is.na(radar_data)] <- 0
  for (metric in metrics)
    radar_data[[metric]] <- pmax(0, pmin(1, radar_data[[metric]]))
  radar_data <- rbind(rep(1, length(metrics)), rep(0, length(metrics)), radar_data)
  rownames(radar_data) <- c("Max", "Min", plot_data$trajectory)
  n_trajectories <- nrow(plot_data)
  colors      <- rainbow(n_trajectories, alpha = 0.3)
  line_colors <- rainbow(n_trajectories, alpha = 0.8)

  grDevices::pdf(file.path(out_S11, "S11_d_radar.pdf"), width = 10, height = 6, pointsize = 12)
  layout(matrix(c(1, 2), ncol = 2), widths = c(3, 1))
  old_par <- par(cex.main = 10 / 12)  # target ~10pt effective title at scale ~1.0
  fmsb::radarchart(radar_data, axistype = 1, pcol = line_colors, pfcol = colors,
                   plwd = 2, plty = 1, cglcol = "grey80", cglty = 1, axislabcol = "grey40",
                   caxislabels = c("0", "0.25", "0.5", "0.75", "1"), calcex = 0.72,
                   cglwd = 0.6, vlcex = 0.9, title = "DOE Metrics")
  par(old_par)
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = plot_data$trajectory, col = line_colors, lty = 1, lwd = 2,
        cex = 0.8, bty = "n")
  layout(1)
  grDevices::dev.off()
  message("  Saved S11_d_radar.pdf")
}, error = function(e) message("  SKIPPED S11-d: ", e$message))

# =============================================================================
# 9. PANEL E — DOE heatmap
# =============================================================================
message("\n[S11-e] DOE heatmap")
tryCatch({
  if (is.null(res_s10)) stop("res_s10 is NULL")
  if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
    sysfonts::font_add("Arial",
      regular    = file.path(repo_root, "manuscript", "fonts", "Arial.ttf"),
      bold       = file.path(repo_root, "manuscript", "fonts", "Arial Bold.ttf"),
      italic     = file.path(repo_root, "manuscript", "fonts", "Arial Italic.ttf"),
      bolditalic = file.path(repo_root, "manuscript", "fonts", "Arial Bold Italic.ttf"))
    showtext::showtext_auto()
    showtext::showtext_opts(dpi = 300)
  }
  # No in-plot title — the assembled main-figure panel carries its own
  # narrative title/subtitle above the panel instead. Sizes are ~15% larger
  # than the previous pass, and everything renders in Arial to match the
  # rest of Figure 1.
  p_hm <- plot(res_s10, type = "heatmap") +
    labs(title = NULL, y = "Root population / initialization") +
    theme(panel.grid = element_blank(),
          text        = element_text(family = "Arial"),
          axis.text   = element_text(size = 15, family = "Arial"),
          axis.title  = element_text(size = 16, family = "Arial"),
          legend.text  = element_text(size = 14, family = "Arial"),
          legend.title = element_text(size = 15, family = "Arial"))
  # Enlarge the in-cell score labels (2nd layer is the geom_text from
  # plot.multi_doe_results()) to match the larger axis/legend text above.
  p_hm$layers[[2]]$aes_params$size <- 5.1
  p_hm$layers[[2]]$aes_params$family <- "Arial"
  ggsave(file.path(out_S11, "S11_e_doe_heatmap.pdf"), p_hm, width = 11.5, height = 5.75, dpi = 300)
  if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto(FALSE)
  message("  Saved S11_e_doe_heatmap.pdf")
}, error = function(e) message("  SKIPPED S11-e: ", e$message))

message("\nFigure S11 complete. Output: ", out_S11)
