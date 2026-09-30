# =============================================================================
# make_figure_S10.R
#
# Supplementary figure panels for:
#   S10 — Stem cell differentiation (branched trajectory)
#
# Panel structure:
#   S10_a_umap.pdf                  GT UMAP + per-method pseudotime UMAPs
#   S10_b_doe_heatmap.pdf           BioTrajX plot(res, scope="branch", type="heatmap"),
#                                    overall aggregate DOE column appended
#   S10_c_branch_umap_pseudotime.pdf Per-branch UMAP grid, per-method pseudotime
#                                    (methods ordered by branch DOE, high to low)
#   S10_d_module_trends.pdf         Module score trends, faceted by method x branch
#                                    (methods ordered by aggregate DOE, high to low)
#   S10_e_gam_stem_to_ery.pdf,      GAM fits of Expr ~ Pseudotime per Method x TF,
#   S10_e_gam_stem_to_b.pdf         one facet plot per branch (methods ordered by
#                                    branch DOE, high to low)
#   S10_f_violin_stem_to_ery.pdf,   Violin of Pseudotime by cell type, one per
#   S10_f_violin_stem_to_b.pdf      branch, faceted by method (ordered by branch DOE)
#
# Prerequisites:
#   data/stem_cell.rds
#   manuscript/results/branch_stemcell/ti_pseudotimes.csv
#     (from run_biotrajx_stemcell.R)
#
# Usage:
#   Rscript manuscript/scripts/real/make_figure_S10.R
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

out_S10 <- file.path(repo_root, "manuscript", "figures", "real", "S10")
dir.create(out_S10, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# SHARED HELPERS
# =============================================================================

#' UMAP coloured by pseudotime (Spectral palette)
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

# =============================================================================
# FIGURE S10 — STEM CELL BRANCHED TRAJECTORY
# =============================================================================
message("\n", strrep("=", 70))
message("FIGURE S10 — Stem cell branched trajectory")
message(strrep("=", 70))

# ── Load Seurat + pseudotimes ─────────────────────────────────────────────────
message("Loading stem_cell.rds ...")
obj_s10 <- tryCatch(
  readRDS(file.path(repo_root, "data", "stem_cell.rds")),
  error = function(e) { message("  ERROR: ", e$message); NULL }
)

csv_s10 <- file.path(repo_root, "manuscript", "results", "branch_stemcell", "ti_pseudotimes.csv")
if (!file.exists(csv_s10))
  stop("ti_pseudotimes.csv not found. Run run_ti_stemcell.R first.")
ti_df_s10 <- read.csv(csv_s10, row.names = 1, check.names = FALSE)
message(sprintf("  Pseudotimes: %d cells × %d methods", nrow(ti_df_s10), ncol(ti_df_s10)))

# ── Marker genes (MSigDB — independently curated, not derived from this dataset) ──
# Early/HSC (both branches): HAY_BONE_MARROW_CD34_POS_HSC (C8)
#   Human Cell Atlas bone marrow CD34+ HSC signature; ortholog-mapped to mouse.
#   Spearman cor with Slingshot PT: Ery branch = -0.735, B branch = -0.465.
# Ery terminal: WP_ERYTHROPOIESIS (C2 / WikiPathways)
#   Curated erythropoiesis pathway; cor with Ery-branch PT = +0.659.
# B terminal: HADDAD_B_LYMPHOCYTE_PROGENITOR (C2)
#   Published B lymphocyte progenitor signature (Haddad et al.); cor = +0.700.
message("  Fetching MSigDB markers for S10 (Mus musculus) ...")
ms_ery_s10 <- get_markers_msigdb(
  early      = "HAY_BONE_MARROW_CD34_POS_HSC",
  terminal   = "WP_ERYTHROPOIESIS",
  collection = NULL,
  species    = "Mus musculus"
)
ms_ery_s10 <- filter_markers(ms_ery_s10, obj_s10, top_n = 30, min_detection = 0.10)

ms_b_s10 <- get_markers_msigdb(
  early      = "HAY_BONE_MARROW_CD34_POS_HSC",
  terminal   = "HADDAD_B_LYMPHOCYTE_PROGENITOR",
  collection = NULL,
  species    = "Mus musculus"
)
ms_b_s10 <- filter_markers(ms_b_s10, obj_s10, top_n = 30, min_detection = 0.10)

stem_progenitor_genes <- ms_ery_s10$early   # same HSC set for both branches
erythrocyte_genes     <- ms_ery_s10$terminal
bcell_genes           <- ms_b_s10$terminal

early_markers_list <- list(
  Stem_to_Ery = stem_progenitor_genes,
  Stem_to_B   = stem_progenitor_genes
)
terminal_markers_list <- list(
  Stem_to_Ery = erythrocyte_genes,
  Stem_to_B   = bcell_genes
)

branch_filters <- list(
  Stem_to_Ery = list(
    include   = c("Stem_Progenitors", "Erythroid_progenitors_Erythroblasts", "Erythrocytes"),
    min_cells = 50
  ),
  Stem_to_B = list(
    include   = c("Stem_Progenitors", "Immature_B", "Mature_B"),
    min_cells = 50
  )
)

#' Methods with a valid DOE score for branch `br`, ordered high to low
.branch_doe_order <- function(res, br, methods) {
  doe_br <- sapply(methods, function(m) {
    b <- res$results[[m]]$branches[[br]]
    if (is.null(b)) NA_real_ else b$DOE_score
  })
  sort(doe_br[!is.na(doe_br)], decreasing = TRUE)
}

# ── Compute branched multi-DOE ────────────────────────────────────────────────
message("  Computing multi-DOE branched for S10 ...")
pt_list_s10 <- as.list(ti_df_s10)
pt_list_s10 <- pt_list_s10[sapply(pt_list_s10, function(x) sum(!is.na(x)) > 10)]

res_s10 <- tryCatch(
  compute_multi_doe_branched(
    expr_or_seurat        = obj_s10,
    pseudotime_list       = pt_list_s10,
    early_markers_list    = early_markers_list,
    terminal_markers_list = terminal_markers_list,
    cluster_labels        = "Phenotype",
    branch_filters        = branch_filters,
    E_method              = "gmm",
    plot_E                = FALSE,
    verbose               = TRUE
  ),
  error = function(e) { message("  ERROR: ", e$message); NULL }
)

# Save branch summary CSV
if (!is.null(res_s10)) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  all_br <- sort(unique(unlist(lapply(res_s10$results, function(x) names(x$branches)))))
  doe_rows <- unlist(lapply(all_br, function(br) lapply(names(res_s10$results), function(m) {
    b <- res_s10$results[[m]]$branches[[br]]
    if (is.null(b)) return(NULL)
    data.frame(Branch=br, Method=m,
               D_early=b$D$D_early%||%NA, D_term=b$D$D_term%||%NA,
               O=b$O$O%||%NA,
               E_early=b$E$E_early%||%NA, E_term=b$E$E_term%||%NA,
               DOE=b$DOE_score%||%NA, stringsAsFactors=FALSE)
  })), recursive = FALSE)
  branch_dir <- file.path(repo_root, "manuscript", "results", "branch_stemcell")
  write.csv(do.call(rbind, Filter(Negate(is.null), doe_rows)),
            file.path(branch_dir, "doe_summary_branch.csv"),
            row.names = FALSE)
  message("  Saved doe_summary_branch.csv")
  if (!is.null(res_s10$comparison_overall))
    write.csv(res_s10$comparison_overall,
              file.path(branch_dir, "doe_summary_overall.csv"),
              row.names = FALSE)
  message("  Saved doe_summary_overall.csv")
}

# ── Panel a: UMAP grid ────────────────────────────────────────────────────────
message("\n[S10-a] UMAP grid")
tryCatch({
  if (is.null(obj_s10)) stop("Seurat object not loaded")
  umap_df_s10 <- as.data.frame(Embeddings(obj_s10, "umap"))
  colnames(umap_df_s10) <- c("UMAP1", "UMAP2")
  umap_df_s10$Phenotype <- obj_s10$Phenotype[rownames(umap_df_s10)]
  phenotype_levels_s10  <- unique(as.character(umap_df_s10$Phenotype))
  phenotype_levels_s10  <- c("Stem_Progenitors",
                             sort(setdiff(phenotype_levels_s10, "Stem_Progenitors")))
  umap_df_s10$Phenotype <- factor(umap_df_s10$Phenotype, levels = phenotype_levels_s10)

  p_gt_s10 <- plot_gt_umap(umap_df_s10, "Phenotype", title = "Cell type") +
    theme(plot.title = element_text(size = 16, face = "bold"))

  # Order methods by aggregate DOE score, high to low
  method_order_s10 <- colnames(ti_df_s10)
  if (!is.null(res_s10) && !is.null(res_s10$comparison_overall)) {
    ov <- res_s10$comparison_overall
    method_order_s10 <- ov$trajectory[order(-ov$aggregate_DOE)]
    method_order_s10 <- method_order_s10[method_order_s10 %in% colnames(ti_df_s10)]
  }

  shared_s10 <- intersect(rownames(umap_df_s10), rownames(ti_df_s10))

  # Root cell: identical across all methods — the Stem_Progenitors centroid
  # cell (see run_ti_stemcell.R), passed as start_cell to every TI method.
  # Persisted to root_cell_stemcell.txt since it can no longer be read off
  # any one method's pseudotime column (unlike CytoTRACE, it isn't a raw
  # score).
  root_cell_path_s10 <- file.path(repo_root, "manuscript", "results", "branch_stemcell",
                                  "root_cell_stemcell.txt")
  root_cell_s10 <- if (file.exists(root_cell_path_s10))
    readLines(root_cell_path_s10, n = 1) else NA_character_
  root_cell_s10 <- if (!is.na(root_cell_s10) && root_cell_s10 %in% shared_s10)
    root_cell_s10 else NA_character_
  root_coord_s10 <- if (!is.na(root_cell_s10))
    umap_df_s10[root_cell_s10, c("UMAP1", "UMAP2")] else NULL

  method_plots_s10 <- lapply(method_order_s10, function(m) {
    doe_val <- if (!is.null(res_s10) && !is.null(res_s10$comparison_overall))
      res_s10$comparison_overall$aggregate_DOE[res_s10$comparison_overall$trajectory == m]
    else NA
    p <- plot_pt_umap(
      data.frame(UMAP1      = umap_df_s10[shared_s10, "UMAP1"],
                 UMAP2      = umap_df_s10[shared_s10, "UMAP2"],
                 Pseudotime = ti_df_s10[shared_s10, m]),
      title = m
    )
    if (!is.null(root_coord_s10))
      p <- p + geom_point(data = root_coord_s10, aes(x = UMAP1, y = UMAP2),
                          colour = "black", shape = 17, size = 4, inherit.aes = FALSE)
    p + labs(subtitle = if (!is.na(doe_val) && length(doe_val) == 1)
               sprintf("DOE = %.3f", doe_val) else NULL) +
      theme(plot.subtitle = element_text(size = 7, colour = "grey40"))
  })

  save_umap_grid(p_gt_s10, method_plots_s10, ncol = 3,
                 out_path = file.path(out_S10, "S10_a_umap.pdf"), width = 12)
  message("  Saved S10_a_umap.pdf")
}, error = function(e) message("  SKIPPED S10-a: ", e$message))

# ── Panel c: per-branch UMAP grid, methods ranked by branch DOE ───────────────
message("\n[S10-c] Branch pseudotime UMAPs (ranked by DOE)")
tryCatch({
  if (is.null(obj_s10)) stop("Seurat object not loaded")
  if (is.null(res_s10)) stop("res_s10 is NULL")

  branch_labels_s10 <- c(Stem_to_Ery = "Ery branch", Stem_to_B = "B branch")

  row_plots_s10 <- lapply(names(branch_labels_s10), function(br) {
    br_cells <- colnames(obj_s10)[obj_s10$Phenotype %in% branch_filters[[br]]$include]
    umap_br  <- umap_df_s10[intersect(rownames(umap_df_s10), br_cells), , drop = FALSE]
    umap_br$Phenotype <- droplevels(factor(umap_br$Phenotype))

    gt_plot <- plot_gt_umap(umap_br, "Phenotype", title = branch_labels_s10[[br]]) +
      theme(plot.title = element_text(size = 16, face = "bold"))

    doe_br <- .branch_doe_order(res_s10, br, colnames(ti_df_s10))

    shared_br <- intersect(rownames(umap_br), rownames(ti_df_s10))
    method_names_br <- names(doe_br)
    last_br_flag    <- br == names(branch_labels_s10)[length(branch_labels_s10)]
    last_m_flag     <- method_names_br[length(method_names_br)]
    method_plots_br <- lapply(method_names_br, function(m) {
      pt_raw  <- ti_df_s10[shared_br, m]
      pt_rng  <- range(pt_raw, na.rm = TRUE)
      pt_norm <- (pt_raw - pt_rng[1]) / (pt_rng[2] - pt_rng[1])
      p <- plot_pt_umap(
        data.frame(UMAP1      = umap_br[shared_br, "UMAP1"],
                   UMAP2      = umap_br[shared_br, "UMAP2"],
                   Pseudotime = pt_norm),
        title = m
      ) +
        scale_colour_distiller(palette = "Spectral", direction = -1,
                               na.value = "grey80", name = "Pseudotime",
                               limits = c(0, 1)) +
        labs(subtitle = sprintf("DOE = %.3f", doe_br[[m]])) +
        theme(plot.subtitle = element_text(size = 7, colour = "grey40"))
      if (!last_br_flag || m != last_m_flag)
        p <- p + theme(legend.position = "none")
      p
    })

    list(gt = gt_plot, methods = method_plots_br)
  })

  if (has_patchwork) {
    row_combined <- lapply(row_plots_s10, function(rp) {
      Reduce(`+`, c(list(rp$gt), rp$methods)) +
        plot_layout(ncol = 1 + length(rp$methods))
    })
    p_c10 <- wrap_plots(row_combined, ncol = 1) +
      plot_layout(guides = "collect") &
      theme(legend.position = "right")
    ggsave(file.path(out_S10, "S10_c_branch_umap_pseudotime.pdf"), p_c10,
           width = 16, height = 6)
    message("  Saved S10_c_branch_umap_pseudotime.pdf")
  } else {
    message("  (patchwork unavailable — skipping S10-c)")
  }
}, error = function(e) message("  SKIPPED S10-c: ", e$message))

# ── Panel b: BioTrajX DOE heatmap (branch scope + overall DOE column) ─────────
message("\n[S10-b] DOE heatmap (BioTrajX, branch scope)")
tryCatch({
  if (is.null(res_s10)) stop("res_s10 is NULL")
  p_b10 <- plot(res_s10, scope = "branch", type = "heatmap", branch_mode = "facet") +
    labs(title = NULL) +
    theme(panel.grid = element_blank())
  ggsave(file.path(out_S10, "S10_b_doe_heatmap.pdf"), p_b10, width = 11, height = 7)
  message("  Saved S10_b_doe_heatmap.pdf")
}, error = function(e) message("  SKIPPED S10-b: ", e$message))

# ── Panel d: Module score trends per method, one row per branch ───────────────
# Arranged like panel c: each branch is its own row, with methods ranked by
# that branch's own DOE score (high to low) rather than a shared aggregate
# order, so each row (including the B-cell branch) carries its own method-name
# strip labels in that row's DOE ranking.
message("\n[S10-d] Module score trends (ranked by branch DOE)")
tryCatch({
  if (is.null(obj_s10)) stop("Seurat object not loaded")
  if (is.null(res_s10)) stop("res_s10 is NULL")
  obj_s10 <- AddModuleScore(obj_s10, features = list(stem_progenitor_genes), name = "early_mod")
  obj_s10 <- AddModuleScore(obj_s10, features = list(erythrocyte_genes),     name = "ery_mod")
  obj_s10 <- AddModuleScore(obj_s10, features = list(bcell_genes),           name = "bcell_mod")

  scores_s10 <- list(
    Early       = setNames(obj_s10$early_mod1, colnames(obj_s10)),
    Erythrocyte = setNames(obj_s10$ery_mod1,   colnames(obj_s10)),
    `B cell`    = setNames(obj_s10$bcell_mod1, colnames(obj_s10))
  )
  module_pal_s10 <- c(Early = "#4E9AF1", Erythrocyte = "#2D7A3E", `B cell` = "#9B59B6")

  row_plots_d10 <- lapply(names(branch_labels_s10), function(br) {
    br_cells <- colnames(obj_s10)[obj_s10$Phenotype %in% branch_filters[[br]]$include]
    mods     <- if (br == "Stem_to_Ery") c("Early", "Erythrocyte") else c("Early", "B cell")

    doe_br          <- .branch_doe_order(res_s10, br, colnames(ti_df_s10))
    method_order_br <- names(doe_br)

    df_br <- do.call(rbind, lapply(method_order_br, function(m) {
      shared <- intersect(br_cells, rownames(ti_df_s10))
      pt <- ti_df_s10[shared, m]; keep <- !is.na(pt)
      do.call(rbind, lapply(mods, function(mod) {
        data.frame(Pseudotime = pt[keep], Score = scores_s10[[mod]][shared][keep],
                   Module = mod, Method = m, stringsAsFactors = FALSE)
      }))
    }))
    df_br$Method <- factor(df_br$Method, levels = method_order_br)
    df_br$Module <- factor(df_br$Module, levels = mods)

    ggplot(df_br, aes(x = Pseudotime, y = Score, colour = Module)) +
      geom_point(size = 0.2, alpha = 0.12) +
      geom_smooth(method = "loess", se = TRUE, linewidth = 0.9, span = 0.5) +
      scale_colour_manual(values = module_pal_s10[mods]) +
      facet_wrap(~ Method, nrow = 1, scales = "free", axes = "all") +
      theme_minimal(base_size = 9) +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold", size = 8),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      labs(x = "Pseudotime [0, 1]", y = "Module score", colour = NULL,
           title = branch_labels_s10[[br]])
  })

  if (has_patchwork) {
    p_d10 <- wrap_plots(row_plots_d10, ncol = 1) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    n_cols_d10 <- max(vapply(row_plots_d10,
                             function(p) nlevels(p$data$Method), integer(1)))
    ggsave(file.path(out_S10, "S10_d_module_trends.pdf"), p_d10,
           width = n_cols_d10 * 2.2 + 1, height = 7)
    message("  Saved S10_d_module_trends.pdf")
  } else {
    message("  (patchwork unavailable — skipping S10-d)")
  }
}, error = function(e) message("  SKIPPED S10-d: ", e$message))

# =============================================================================
# Panel e: GAM fits of Expr ~ Pseudotime, per Method x TF, one plot per branch
# =============================================================================
gam_genes <- list(
  Stem_to_Ery = list(early = c("Hlf", "Meis1"), late = c("Klf1", "Nfe2")),
  Stem_to_B   = list(early = c("Hlf", "Meis1"), late = c("Bach2", "Ikzf3"))
)

for (br in names(gam_genes)) {
  message(sprintf("\n[S10-e] GAM fits x Pseudotime (%s)", br))
  tryCatch({
    if (is.null(obj_s10)) stop("Seurat object not loaded")
    if (is.null(res_s10)) stop("res_s10 is NULL")

    early_genes_br <- gam_genes[[br]]$early
    late_genes_br  <- gam_genes[[br]]$late
    genes <- c(early_genes_br, late_genes_br)
    expr_mat <- GetAssayData(obj_s10, layer = "data")
    genes_present <- intersect(genes, rownames(expr_mat))
    if (length(genes_present) == 0) stop("none of the requested TFs found in obj_s10")

    br_cells    <- colnames(obj_s10)[obj_s10$Phenotype %in% branch_filters[[br]]$include]
    doe_br      <- .branch_doe_order(res_s10, br, colnames(ti_df_s10))
    method_order_br <- names(doe_br)

    df <- do.call(rbind, lapply(method_order_br, function(m) {
      cells <- intersect(br_cells, rownames(ti_df_s10))
      pt    <- ti_df_s10[cells, m]; keep <- !is.na(pt)
      cells <- cells[keep]; pt <- pt[keep]
      do.call(rbind, lapply(genes_present, function(g) {
        data.frame(Method = m, Gene = g, Pseudotime = pt,
                   Expr = as.numeric(expr_mat[g, cells]),
                   GeneGroup = if (g %in% early_genes_br) "Early" else "Terminal",
                   stringsAsFactors = FALSE)
      }))
    }))
    df$Method    <- factor(df$Method,    levels = method_order_br)
    df$Gene      <- factor(df$Gene,      levels = genes_present)
    df$GeneGroup <- factor(df$GeneGroup, levels = c("Early", "Terminal"))

    # Terminal colour matches the module-score colour for this branch's
    # terminal population in panel D (Erythrocyte = green, B cell = purple).
    terminal_colour_br <- if (br == "Stem_to_B") "#9B59B6" else "#2D7A3E"

    p_e10 <- ggplot(df, aes(x = Pseudotime, y = Expr)) +
      geom_point(size = 0.2, alpha = 0.12, colour = "grey50") +
      geom_smooth(aes(colour = GeneGroup, fill = GeneGroup),
                  method = "gam", formula = y ~ s(x, bs = "cs"),
                  se = TRUE, linewidth = 0.9) +
      scale_colour_manual(values = c(Early = "#4E9AF1", Terminal = terminal_colour_br), guide = "none") +
      scale_fill_manual(values = c(Early = "#4E9AF1", Terminal = terminal_colour_br), guide = "none") +
      ggh4x::facet_grid2(rows = vars(Method), cols = vars(GeneGroup, Gene),
                         scales = "free_y", independent = "y", axes = "all",
                         strip = ggh4x::strip_nested()) +
      theme_minimal(base_size = 9) +
      theme(strip.text = element_text(face = "bold", size = 11),
            strip.background = element_rect(fill = "grey85", colour = "grey20"),
            panel.border = element_rect(colour = "grey20", fill = NA),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = "Pseudotime [0, 1]", y = "Expression", title = gsub("_", " ", br))

    out_path <- file.path(out_S10, sprintf("S10_e_gam_%s.pdf", tolower(br)))
    ggsave(out_path, p_e10,
           width = length(genes_present) * 2.2 + 1,
           height = length(method_order_br) * 1.8 + 1)
    message("  Saved ", basename(out_path))
  }, error = function(e) message(sprintf("  SKIPPED S10-e (%s): ", br), e$message))
}

# =============================================================================
# Panel f: violin plot of Pseudotime by cell type, one plot per branch
# =============================================================================
for (br in names(branch_filters)) {
  message(sprintf("\n[S10-f] Pseudotime by cell type, violin (%s)", br))
  tryCatch({
    if (is.null(obj_s10)) stop("Seurat object not loaded")
    if (is.null(res_s10)) stop("res_s10 is NULL")

    br_cells  <- colnames(obj_s10)[obj_s10$Phenotype %in% branch_filters[[br]]$include]
    phenotype <- droplevels(factor(obj_s10$Phenotype[br_cells]))
    names(phenotype) <- br_cells

    doe_br          <- .branch_doe_order(res_s10, br, colnames(ti_df_s10))
    method_order_br <- names(doe_br)

    df <- do.call(rbind, lapply(method_order_br, function(m) {
      cells <- intersect(br_cells, rownames(ti_df_s10))
      pt    <- ti_df_s10[cells, m]; keep <- !is.na(pt)
      cells <- cells[keep]
      data.frame(Method = m, CellType = phenotype[cells], Pseudotime = pt[keep],
                 stringsAsFactors = FALSE)
    }))
    df$Method   <- factor(df$Method,   levels = method_order_br)
    df$CellType <- factor(df$CellType, levels = branch_filters[[br]]$include)

    doe_labels <- data.frame(Method = names(doe_br),
                             label  = sprintf("DOE = %.3f", doe_br),
                             stringsAsFactors = FALSE)
    doe_labels$Method <- factor(doe_labels$Method, levels = method_order_br)

    p_f10 <- ggplot(df, aes(x = CellType, y = Pseudotime)) +
      geom_violin(aes(fill = CellType), scale = "width", colour = NA,
                  alpha = 0.6, trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.size = 0.4, outlier.alpha = 0.3,
                   fill = "white", alpha = 0.7, linewidth = 0.3) +
      geom_text(data = doe_labels, aes(label = label), x = -Inf, y = Inf,
                hjust = -0.05, vjust = 1.3, size = 3.2, colour = "grey20",
                inherit.aes = FALSE) +
      facet_wrap(~ Method, ncol = 3, axes = "all") +
      scale_fill_hue(guide = "none") +
      scale_x_discrete(labels = function(x)
        gsub("Erythroid_progenitors_Erythroblasts", "Erythroid\nprogenitors\nErythroblasts", x)) +
      theme_classic(base_size = 14) +
      theme(strip.text = element_text(face = "bold", size = 14),
            strip.background = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 15),
            plot.title = element_text(hjust = 0.5)) +
      labs(x = "Cell type", y = "Pseudotime", title = gsub("_", " ", br))

    n_m <- length(method_order_br)
    out_path <- file.path(out_S10, sprintf("S10_f_violin_%s.pdf", tolower(br)))
    ggsave(out_path, p_f10, width = 10.5, height = ceiling(n_m / 3) * 3.2 + 1)
    message("  Saved ", basename(out_path))
  }, error = function(e) message(sprintf("  SKIPPED S10-f (%s): ", br), e$message))
}

message("\nFigure S10 complete. Output: ", out_S10)
