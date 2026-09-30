# =============================================================================
# make_figure_S8.R
#
# Supplementary figure panels for:
#   S8 — GSE131847 LCMV time-course CD8 T-cell dataset (naive -> d90
#        post-infection), including a validation step unique to this
#        dataset: since the true day of infection is known for every cell
#        (cell_type), we can directly check whether a method's DOE score
#        (label-free) tracks how well its pseudotime recovers the known
#        day ordering (ground-truth check) — mouse LCMV-specific
#        naive/effector marker genes are from GSE41867.
#
# Root comparison: each of the 7 directed TI methods is run under two roots --
#   NCR (Naive_Centroid_Root)  naive cell nearest the D0/naive centroid in PCA
#                              space (linear_gse131847_d0centroid)
#   CR  (CytoTRACE_Root)       cell with the lowest CytoTRACE score genome-wide
#                              (linear_gse131847_cytoglobal)
# plus 1 root-independent CytoTRACE column (run_cytotrace() never receives a
# start_cell) = 15 trajectory columns total, suffixed "_NCR"/"_CR". Both
# halves are read straight off their own results dirs -- no TI methods are
# rerun here, they're merged in R.
#
# Panel structure:
#   S8_a_umap.pdf              UMAP grid: ground-truth day + per-method-per-root
#                              pseudotime, root marked in NCR/CR colour
#   S8_b_module_trends.pdf     early/terminal module score vs pseudotime, per method
#   S8_c_doe_heatmap.pdf       BioTrajX DOE heatmap
#   S8_d_day_corr_vs_doe.pdf   per-method (pseudotime vs. true day) correlation vs. DOE score
#   S8_e_pseudotime_vs_day.pdf per-method pseudotime vs. true day scatter + Spearman rho
#
# Prerequisites:
#   data/GSE131847_seu.rds
#   manuscript/results/linear_gse131847_d0centroid/ti_pseudotime_gse131847.csv (from run_ti_gse131847.R d0centroid)
#   manuscript/results/linear_gse131847_cytoglobal/ti_pseudotime_gse131847.csv (from run_ti_gse131847.R cytoglobal)
#
# Usage:
#   Rscript manuscript/scripts/real/make_figure_S8.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
})

repo_root <- here::here()

library(BioTrajX)

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

out_dir     <- file.path(repo_root, "manuscript", "figures", "real", "S8")
results_dir <- file.path(repo_root, "manuscript", "results", "linear_gse131847_ncr_cr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# NCR = Naive_Centroid_Root, CR = CytoTRACE_Root -- the two roots compared.
# Shared colours/labels/shapes so panel a's root markers and panel d's shape
# legend stay in sync. CytoTRACE itself is root-independent (run_cytotrace()
# never receives start_cell) and is folded into CR for labelling purposes,
# since CR IS the cell it naturally ranks lowest genome-wide -- diamond marks
# that combined CR/CytoTRACE category, triangle marks NCR.
NCR_LABEL <- "Naive-centroid root (NCR)"
CR_LABEL  <- "CytoTRACE root (CR)"
NCR_COLOR <- "blue"
CR_COLOR  <- "red"
ROOT_DIRS <- c(NCR = "linear_gse131847_d0centroid", CR = "linear_gse131847_cytoglobal")

BASE_A <- 11.63; HEAD_A <- 12.64; SUB_A <- 10.11

plot_pt_umap <- function(df, title = NULL) {
  ggplot(df, aes(x = UMAP1, y = UMAP2, colour = Pseudotime)) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_colour_distiller(palette = "Spectral", direction = -1,
                           na.value = "grey80", name = "Pseudotime") +
    theme_classic(base_size = BASE_A) +
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = HEAD_A, face = "bold")) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2")
}

plot_gt_umap <- function(umap_df, label_col, title = "Cell type") {
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = .data[[label_col]])) +
    geom_point(size = 0.5, alpha = 0.7) +
    theme_classic(base_size = BASE_A) +
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = HEAD_A, face = "bold")) +
    labs(title = title, colour = NULL, x = "UMAP 1", y = "UMAP 2")
}

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
    message("  (patchwork unavailable -- saving GT UMAP only)")
    ggsave(out_path, gt_plot, width = 7, height = 5)
  }
}

#' Methods with a DOE score, ordered high to low
.doe_order <- function(res, methods) {
  if (is.null(res)) return(methods)
  cs <- res$comparison_summary
  ord <- cs$trajectory[order(-cs$DOE_score)]
  ord[ord %in% methods]
}

message("\n", strrep("=", 70))
message("GSE131847 (LCMV time course) -- DOE analysis")
message(strrep("=", 70))

# ── Load Seurat + pseudotimes ─────────────────────────────────────────────────
message("Loading GSE131847_seu.rds ...")
obj <- readRDS(file.path(repo_root, "data", "GSE131847_seu.rds"))
DefaultAssay(obj) <- "SCT"

ti_methods_8 <- c("PAGA-DPT", "CytoTRACE", "Monocle3", "DPT", "Slingshot", "SCORPIUS",
                  "TSCAN", "Palantir")
directed_8   <- setdiff(ti_methods_8, "CytoTRACE")

# Dual-root data assembled from two already-computed, single-root runs -- no
# TI methods are rerun. CytoTRACE itself is root-independent (run_cytotrace()
# never receives start_cell), so its score is identical across both source
# runs -- taken from the first one listed (NCR).
root_csvs <- lapply(ROOT_DIRS, function(d) {
  p <- file.path(repo_root, "manuscript", "results", d, "ti_pseudotime_gse131847.csv")
  if (!file.exists(p))
    stop("Missing ti_pseudotime_gse131847.csv in ", d, " -- run run_ti_gse131847.R first.")
  read.csv(p, row.names = 1, check.names = FALSE)
})
stopifnot(identical(rownames(root_csvs[[1]]), rownames(root_csvs[[2]])))

directed_parts <- Map(function(csv, suffix) {
  part <- csv[, intersect(colnames(csv), directed_8), drop = FALSE]
  colnames(part) <- paste0(colnames(part), "_", suffix)
  part
}, root_csvs, names(ROOT_DIRS))

# unname() on directed_parts is essential: cbind()-ing a NAMED data.frame
# argument prefixes every one of its columns with that name (e.g.
# "NCR.Slingshot_NCR") -- it only names the column directly for a plain
# vector argument, which is why the CytoTRACE element below stays named.
ti_df <- do.call(cbind, c(unname(directed_parts), list(CytoTRACE = root_csvs[[1]][["CytoTRACE"]])))
rownames(ti_df) <- rownames(root_csvs[[1]])

csv_path <- file.path(results_dir, "ti_pseudotime_gse131847.csv")
write.csv(ti_df, csv_path, row.names = TRUE)
message("  Saved combined NCR/CR pseudotimes: ", csv_path)

root_cells <- vapply(ROOT_DIRS, function(d)
  readLines(file.path(repo_root, "manuscript", "results", d, "root_cell_gse131847.txt"), n = 1),
  character(1))
roots_out <- file.path(results_dir, "root_cells_gse131847.csv")
write.csv(data.frame(root_type = names(ROOT_DIRS), cell = unname(root_cells)),
          roots_out, row.names = FALSE)
message("  Saved: ", roots_out)
message(sprintf("  Pseudotimes: %d cells x %d methods", nrow(ti_df), ncol(ti_df)))

# ── Marker genes (MSigDB, mouse LCMV: naive vs. day8 LCMV effector) ─────────
# Naive pole:    GSE41867_NAIVE_VS_DAY8_LCMV_EFFECTOR_CD8_TCELL_UP
#   Genes higher in naive vs. day-8 LCMV-infected effector CD8 T cells.
# Terminal pole: GSE41867_NAIVE_VS_DAY8_LCMV_EFFECTOR_CD8_TCELL_DN
#   Genes higher in day-8 LCMV-infected effector CD8 T cells vs. naive.
# Generic (not Armstrong/Clone13-specific) pairing, chosen to make no
# assumption about whether this time course resolves (memory) or persists
# (exhaustion) past day 8.
message("  Fetching MSigDB (mouse) markers for GSE131847 ...")
ms <- get_markers_msigdb(
  early      = "GSE41867_NAIVE_VS_DAY8_LCMV_EFFECTOR_CD8_TCELL_UP",
  terminal   = "GSE41867_NAIVE_VS_DAY8_LCMV_EFFECTOR_CD8_TCELL_DN",
  collection = NULL,
  species    = "Mus musculus"
)
ms           <- filter_markers(ms, obj, top_n = 30, min_detection = 0.10)
early_genes  <- ms$early
term_genes   <- ms$terminal
message(sprintf("  early (naive) markers: %d, terminal (LCMV effector) markers: %d",
                length(early_genes), length(term_genes)))

# ── Compute multi-DOE ────────────────────────────────────────────────────────
message("  Computing multi-DOE for GSE131847 ...")
pt_list <- as.list(ti_df)
pt_list <- pt_list[sapply(pt_list, function(x) sum(!is.na(x)) > 10)]

res <- tryCatch(
  compute_multi_doe_linear(
    expr_or_seurat   = obj,
    pseudotime_list  = pt_list,
    early_markers    = early_genes,
    terminal_markers = term_genes,
    E_method         = "gmm",
    plot_E           = FALSE,
    verbose          = TRUE
  ),
  error = function(e) { message("  ERROR: ", e$message); NULL }
)
if (!is.null(res))
  write.csv(res$comparison_summary,
            file.path(results_dir, "doe_scores_gse131847.csv"),
            row.names = FALSE)

# ── Panel a: UMAP grid ────────────────────────────────────────────────────────
message("\n[a] UMAP grid")
tryCatch({
  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$CellType  <- obj$cell_type[rownames(umap_df)]

  p_gt <- plot_gt_umap(umap_df, "CellType", title = "Day post-infection") +
    theme(plot.title = element_text(size = 16, face = "bold"))

  method_order <- .doe_order(res, colnames(ti_df))
  shared <- intersect(rownames(umap_df), rownames(ti_df))

  # Two exact roots, one per directed-method run -- read back verbatim rather
  # than approximated. The single CytoTRACE column is folded into CR (that
  # IS the cell it naturally ranks lowest genome-wide), so it shares CR's
  # colour/label -- diamond marks CR/CytoTRACE, triangle marks NCR.
  root_types_a   <- c("NCR", "CR")
  root_label_map <- c(NCR = NCR_LABEL, CR = CR_LABEL)
  root_color_map <- c(NCR = NCR_COLOR, CR = CR_COLOR)
  root_shape_map <- c(NCR = 24, CR = 23)  # triangle / diamond

  roots_df    <- read.csv(file.path(results_dir, "root_cells_gse131847.csv"),
                          stringsAsFactors = FALSE)
  root_coords <- setNames(lapply(names(ROOT_DIRS), function(rt) {
    cell <- roots_df$cell[roots_df$root_type == rt]
    if (length(cell) == 1 && cell %in% rownames(umap_df))
      umap_df[cell, c("UMAP1", "UMAP2")] else NULL
  }), names(ROOT_DIRS))

  # Which root type does column m belong to? ("CytoTRACE" -> CR)
  root_type_of <- function(m) {
    if (m == "CytoTRACE") return("CR")
    if (grepl("_NCR$", m)) return("NCR")
    if (grepl("_CR$",  m)) return("CR")
    NA_character_
  }
  root_label_for <- function(m) { rt <- root_type_of(m); if (is.na(rt)) NA_character_ else root_label_map[[rt]] }
  # Prettify facet titles: "Slingshot_NCR" -> "Slingshot [NCR]"
  display_title <- function(m) {
    if (m == "CytoTRACE") return("CytoTRACE [CR]")
    m <- sub("_NCR$", " [NCR]", m)
    m <- sub("_CR$",  " [CR]",  m)
    m
  }
  ncol_a <- 4
  root_fill_values  <- setNames(root_color_map, root_label_map)
  root_shape_values <- setNames(root_shape_map, root_label_map)

  method_plots <- lapply(method_order, function(m) {
    doe_val <- if (!is.null(res))
      res$comparison_summary$DOE_score[res$comparison_summary$trajectory == m]
    else NA
    p <- plot_pt_umap(
      data.frame(UMAP1      = umap_df[shared, "UMAP1"],
                 UMAP2      = umap_df[shared, "UMAP2"],
                 Pseudotime = ti_df[shared, m]),
      title = display_title(m)
    )
    # Every panel gets a row for EACH root type (one real geom_point per
    # level), so every panel's built-in legend already has all keys with
    # correct fill colours/shapes -- the inactive roots are just drawn
    # fully transparent. (A phantom break via scale `limits` alone leaves
    # an empty glyph for a level absent from a panel's own data -- ggplot
    # needs a real row.)
    root_pts <- do.call(rbind, lapply(root_types_a, function(rt) {
      if (is.null(root_coords[[rt]])) return(NULL)
      data.frame(root_coords[[rt]], Root = root_label_map[[rt]])
    }))
    root_pts$Active <- root_pts$Root == root_label_for(m)
    p <- p + geom_point(data = root_pts,
                        aes(x = UMAP1, y = UMAP2, fill = Root, shape = Root, alpha = Active),
                        colour = "black", size = 4, stroke = 0.6,
                        inherit.aes = FALSE) +
      scale_fill_manual(name = "Root", values = root_fill_values,
                        guide = guide_legend(override.aes = list(alpha = 1),
                                              title.theme = element_text(size = 11),
                                              label.theme = element_text(size = 9))) +
      scale_shape_manual(name = "Root", values = root_shape_values,
                         guide = guide_legend(override.aes = list(alpha = 1),
                                               title.theme = element_text(size = 11),
                                               label.theme = element_text(size = 9))) +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0), guide = "none")
    p + labs(subtitle = if (!is.na(doe_val) && length(doe_val) == 1)
               sprintf("DOE = %.3f", doe_val) else NULL) +
      theme(plot.subtitle = element_text(size = SUB_A, colour = "grey40"))
  })

  save_umap_grid(p_gt, method_plots, ncol = ncol_a,
                 out_path = file.path(out_dir, "S8_a_umap.pdf"), width = 14)
  message("  Saved S8_a_umap.pdf")
}, error = function(e) message("  SKIPPED a: ", e$message))

# ── Panel b: Module score trends per method ────────────────────────────────────
message("\n[b] Module score trends")
tryCatch({
  obj <- AddModuleScore(obj, features = list(early_genes), name = "early_mod")
  obj <- AddModuleScore(obj, features = list(term_genes),  name = "term_mod")

  scores <- list(
    Early    = setNames(obj$early_mod1, colnames(obj)),
    Terminal = setNames(obj$term_mod1,  colnames(obj))
  )

  cells_all <- Reduce(intersect, lapply(scores, names))
  cells_all <- intersect(cells_all, rownames(ti_df))

  trend_long <- do.call(rbind, lapply(colnames(ti_df), function(m) {
    pt   <- ti_df[cells_all, m]
    keep <- !is.na(pt)
    do.call(rbind, lapply(names(scores), function(mod) {
      data.frame(Pseudotime = pt[keep], Score = scores[[mod]][cells_all][keep],
                 Module = mod, Method = m, stringsAsFactors = FALSE)
    }))
  }))
  trend_long$Module <- factor(trend_long$Module, levels = c("Early", "Terminal"))
  trend_long$Method <- factor(trend_long$Method, levels = .doe_order(res, colnames(ti_df)))

  # 15 method-root columns -- 3x5 (5 cols) rather than the default 4-wide wrap.
  ncol_b <- 5

  p_b <- ggplot(trend_long, aes(x = Pseudotime, y = Score, colour = Module)) +
    geom_point(size = 0.2, alpha = 0.15) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 0.9, span = 0.4) +
    scale_colour_manual(values = c(Early = "#4E9AF1", Terminal = "#2D7A3E")) +
    facet_wrap(~ Method, ncol = ncol_b, scales = "free_x", axes = "all") +
    theme_minimal(base_size = 11.63) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 12.64),
          panel.grid.minor = element_blank()) +
    labs(x = "Pseudotime [0, 1]", y = "Module score", colour = NULL)

  n_m <- length(unique(trend_long$Method))
  ggsave(file.path(out_dir, "S8_b_module_trends.pdf"), p_b,
         width = ncol_b * 3.5, height = ceiling(n_m / ncol_b) * 3 + 1)
  message("  Saved S8_b_module_trends.pdf")
}, error = function(e) message("  SKIPPED b: ", e$message))

# ── Panel c: BioTrajX DOE heatmap (style/font size matching S11 panel e) ─────
message("\n[c] DOE heatmap (BioTrajX)")
tryCatch({
  if (is.null(res)) stop("res is NULL")
  if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
    sysfonts::font_add("Arial",
      regular    = file.path(repo_root, "manuscript", "fonts", "Arial.ttf"),
      bold       = file.path(repo_root, "manuscript", "fonts", "Arial Bold.ttf"),
      italic     = file.path(repo_root, "manuscript", "fonts", "Arial Italic.ttf"),
      bolditalic = file.path(repo_root, "manuscript", "fonts", "Arial Bold Italic.ttf"))
    showtext::showtext_auto()
    showtext::showtext_opts(dpi = 300)
  }
  p_c <- plot(res, type = "heatmap") +
    labs(title = NULL) +
    theme(panel.grid = element_blank(),
          text        = element_text(family = "Arial"),
          axis.text   = element_text(size = 15, family = "Arial"),
          axis.title  = element_text(size = 16, family = "Arial"),
          legend.text  = element_text(size = 14, family = "Arial"),
          legend.title = element_text(size = 15, family = "Arial"))
  p_c$layers[[2]]$aes_params$size <- 5.1
  p_c$layers[[2]]$aes_params$family <- "Arial"
  ggsave(file.path(out_dir, "S8_c_doe_heatmap.pdf"), p_c, width = 11.5, height = 7.5, dpi = 300)
  if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto(FALSE)
  message("  Saved S8_c_doe_heatmap.pdf")
}, error = function(e) message("  SKIPPED c: ", e$message))

# =============================================================================
# Panel d: ground-truth day recovery vs. DOE score
#
# Unlike the CD8T dataset (S8), here every cell has a known day of infection
# (cell_type). This lets us directly test whether the label-free DOE score
# tracks how well a method's pseudotime recovers the *actual* temporal
# ordering -- rather than relying on a proxy like literature TF recovery.
# =============================================================================
message("\n[d] Pseudotime-vs-true-day recovery, compared against DOE score")
tryCatch({
  if (is.null(res)) stop("res is NULL")

  day_lookup <- c(naive = 0, d3 = 3, d4 = 4, d5 = 5, d6 = 6, d7 = 7,
                  d10 = 10, d14 = 14, d21 = 21, d32 = 32, d60 = 60, d90 = 90)
  true_day <- day_lookup[as.character(obj$cell_type)]
  names(true_day) <- colnames(obj)

  doe_lookup <- setNames(res$comparison_summary$DOE_score,
                         res$comparison_summary$trajectory)

  day_corr <- do.call(rbind, lapply(colnames(ti_df), function(m) {
    pt    <- ti_df[[m]]
    names(pt) <- rownames(ti_df)
    cells <- intersect(names(pt)[!is.na(pt)], names(true_day))
    if (length(cells) < 10) return(NULL)

    sp <- suppressWarnings(cor.test(pt[cells], true_day[cells], method = "spearman"))
    pe <- suppressWarnings(cor.test(pt[cells], true_day[cells], method = "pearson"))

    data.frame(
      Method       = m,
      DOE_score    = unname(doe_lookup[m]),
      n_cells      = length(cells),
      spearman_rho = unname(sp$estimate),
      spearman_p   = sp$p.value,
      pearson_r    = unname(pe$estimate),
      pearson_p    = pe$p.value,
      stringsAsFactors = FALSE
    )
  }))
  day_corr <- day_corr[order(-day_corr$DOE_score), ]

  # Root type per method, for panel d's shape legend (explains what NCR/CR
  # mean directly in the plot rather than relying on the facet titles alone).
  # CytoTRACE is folded into CR (that IS the cell it naturally ranks lowest
  # genome-wide), so it shares CR's diamond shape rather than getting its own.
  day_corr$RootType <- ifelse(grepl("_NCR$", day_corr$Method), NCR_LABEL, CR_LABEL)
  day_corr$RootType <- factor(day_corr$RootType, levels = c(NCR_LABEL, CR_LABEL))
  # Point labels drop the "_NCR"/"_CR" suffix -- shape already encodes root
  # type, so the label just needs the method name (CytoTRACE has no suffix).
  day_corr$MethodLabel <- sub("_(NCR|CR)$", "", day_corr$Method)

  message("\nPer-method: pseudotime vs. true day-of-infection correlation, and DOE score")
  print(day_corr[, c("Method", "DOE_score", "spearman_rho", "spearman_p")])

  write.csv(day_corr, file.path(results_dir, "day_recovery_vs_doe_gse131847.csv"),
            row.names = FALSE)
  message("  Saved day_recovery_vs_doe_gse131847.csv")

  # Correlate the day-recovery metric (Spearman rho vs. true day) itself
  # against DOE score, across methods -- analogous to S8 panel d, but using
  # ground-truth day recovery instead of literature TF recovery.
  if (nrow(day_corr) >= 3) {
    sp_meta <- suppressWarnings(cor.test(day_corr$DOE_score, day_corr$spearman_rho,
                                         method = "spearman"))
    pe_meta <- cor.test(day_corr$DOE_score, day_corr$spearman_rho, method = "pearson")

    fmt_p <- function(p) if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
    message(sprintf("\nDOE score vs. day-recovery (rho): r=%.3f (%s), rho=%.3f (%s)",
                    unname(pe_meta$estimate), fmt_p(pe_meta$p.value),
                    unname(sp_meta$estimate), fmt_p(sp_meta$p.value)))

    score_range <- range(day_corr$DOE_score)
    ann_x <- score_range[1] + 0.2 * diff(score_range)

    # showtext (not the base pdf() device) renders the "ρ" glyph -- the base
    # device's mbcsToSbcs font-encoding path can't represent it and silently
    # drops/garbles the character even though the source string is valid UTF-8.
    if (requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
      sysfonts::font_add("Arial",
        regular    = file.path(repo_root, "manuscript", "fonts", "Arial.ttf"),
        bold       = file.path(repo_root, "manuscript", "fonts", "Arial Bold.ttf"),
        italic     = file.path(repo_root, "manuscript", "fonts", "Arial Italic.ttf"),
        bolditalic = file.path(repo_root, "manuscript", "fonts", "Arial Bold Italic.ttf"))
      showtext::showtext_auto()
      showtext::showtext_opts(dpi = 300)
    }

    root_shape_values <- setNames(c(17, 18), c(NCR_LABEL, CR_LABEL))  # triangle / diamond

    p_d <- ggplot(day_corr, aes(x = DOE_score, y = spearman_rho)) +
      geom_smooth(method = "lm", se = TRUE, colour = "#AAAAAA",
                 fill = "#DDDDDD", linewidth = 0.8) +
      geom_point(aes(colour = DOE_score, shape = RootType), size = 5) +
      geom_text_repel(aes(label = MethodLabel), size = 4.3, max.overlaps = 20,
                      box.padding = 0.5, point.padding = 0.3, force = 3,
                      family = "Arial") +
      scale_colour_viridis_c(name = "DOE score", option = "plasma", direction = -1) +
      scale_shape_manual(name = "Root", values = root_shape_values,
                         guide = guide_legend(
                           title.theme = element_text(size = 13, family = "Arial"),
                           label.theme = element_text(size = 11, family = "Arial"))) +
      annotate("text", x = ann_x, y = max(day_corr$spearman_rho) * 0.98,
               size = 5, colour = "grey30", family = "Arial",
               label = sprintf("r = %.2f (%s)",
                               unname(pe_meta$estimate), fmt_p(pe_meta$p.value))) +
      theme_classic(base_size = 18) +
      theme(text = element_text(family = "Arial"),
            legend.position = "right",
            axis.text  = element_text(size = 16, family = "Arial"),
            axis.title = element_text(size = 18, family = "Arial"),
            legend.text  = element_text(size = 14, family = "Arial"),
            legend.title = element_text(size = 16, family = "Arial")) +
      labs(x = "DOE score",
           y = "Spearman ρ (pseudotime vs. infection day)")

    ggsave(file.path(out_dir, "S8_d_day_corr_vs_doe.pdf"), p_d,
           width = 8.5, height = 6, dpi = 300)
    if (requireNamespace("showtext", quietly = TRUE)) showtext::showtext_auto(FALSE)
    message("  Saved S8_d_day_corr_vs_doe.pdf")
  }
}, error = function(e) message("  SKIPPED d: ", e$message))

# =============================================================================
# Panel e: pseudotime vs. true day-of-infection, per method
#
# Violin + boxplot per true day (categorical x-axis), one panel per method,
# ordered by Spearman rho (pseudotime vs. true day) -- the per-cell
# distribution underlying panel d's summary.
# =============================================================================
message("\n[e] Pseudotime vs. true day-of-infection")
tryCatch({
  day_lookup <- c(naive = 0, d3 = 3, d4 = 4, d5 = 5, d6 = 6, d7 = 7,
                  d10 = 10, d14 = 14, d21 = 21, d32 = 32, d60 = 60, d90 = 90)
  day_levels <- names(day_lookup)
  true_day   <- day_lookup[as.character(obj$cell_type)]
  day_label  <- factor(as.character(obj$cell_type), levels = day_levels)
  names(true_day)  <- colnames(obj)
  names(day_label) <- colnames(obj)

  day_long <- do.call(rbind, lapply(colnames(ti_df), function(m) {
    pt    <- ti_df[[m]]
    names(pt) <- rownames(ti_df)
    cells <- intersect(names(pt)[!is.na(pt)], names(true_day))
    if (length(cells) < 10) return(NULL)
    data.frame(Method = m, Day = day_label[cells], TrueDay = true_day[cells],
               Pseudotime = pt[cells], stringsAsFactors = FALSE)
  }))
  day_long$Day <- factor(day_long$Day, levels = day_levels)

  doe_lookup <- if (!is.null(res))
    setNames(res$comparison_summary$DOE_score, res$comparison_summary$trajectory)
  else NULL

  rho_labels <- do.call(rbind, lapply(split(day_long, day_long$Method), function(d) {
    sp  <- suppressWarnings(cor.test(d$Pseudotime, d$TrueDay, method = "spearman"))
    m   <- unique(d$Method)
    doe <- if (!is.null(doe_lookup) && m %in% names(doe_lookup)) doe_lookup[[m]] else NA
    data.frame(Method = m,
               rho    = unname(sp$estimate),
               DOE    = doe,
               label  = if (!is.na(doe))
                 sprintf("rho = %.2f     DOE = %.2f", unname(sp$estimate), doe)
               else sprintf("rho = %.2f", unname(sp$estimate)),
               stringsAsFactors = FALSE)
  }))
  # Order facets by DOE score, high to low (falls back to rho if DOE unavailable)
  method_order <- if (!is.null(doe_lookup))
    rho_labels$Method[order(-rho_labels$DOE)]
  else rho_labels$Method[order(-rho_labels$rho)]
  day_long$Method    <- factor(day_long$Method,    levels = method_order)
  rho_labels$Method  <- factor(rho_labels$Method,  levels = method_order)

  # 15 method-root columns -- 3x5 (5 cols) rather than the default 3-wide wrap.
  ncol_e <- 5

  p_e <- ggplot(day_long, aes(x = Day, y = Pseudotime)) +
    geom_violin(aes(fill = Day), scale = "width", colour = NA, alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.size = 0.4, outlier.alpha = 0.3,
                 fill = "white", alpha = 0.7, linewidth = 0.3) +
    geom_text(data = rho_labels, aes(label = label), x = -Inf, y = Inf,
              hjust = -0.025, vjust = 1.1, size = 3.4, colour = "grey20",
              inherit.aes = FALSE) +
    facet_wrap(~ Method, ncol = ncol_e, axes = "all") +
    scale_fill_hue(guide = "none") +
    theme_classic(base_size = 11.63) +
    theme(strip.text = element_text(face = "bold", size = 12.64),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Day post-infection (true)", y = "Pseudotime")

  n_m <- length(levels(day_long$Method))
  ggsave(file.path(out_dir, "S8_e_pseudotime_vs_day.pdf"), p_e,
         width = ncol_e * 3.5, height = ceiling(n_m / ncol_e) * 3.2 + 1)
  message("  Saved S8_e_pseudotime_vs_day.pdf")
}, error = function(e) message("  SKIPPED e: ", e$message))

message("\nGSE131847 DOE analysis complete. Output: ", out_dir)
