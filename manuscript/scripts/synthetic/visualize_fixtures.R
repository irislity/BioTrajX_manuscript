# =============================================================================
# visualize_fixtures.R
#
# Visualization script for BioTrajX test fixture scenarios
#
# CEILING / CORRECTNESS fixtures  (Gaussian noise — known expected values)
#   - Gaussian linear      D≈1, O≈1, E≈1
#   - Gaussian reversed    D≈0, O≈1, E≈0
#   - Endpoint             D≈0.8, O≈0.85, E≈1  (step-function ZINB; E ground truth)
#
# STRESS / ROBUSTNESS fixtures  (ZINB noise — directional expectations)
#   - ZINB linear          baseline
#   - ZINB reversed        D↓, E↓, O unchanged (orientation-invariant)
#   - Flat / no signal     all metrics near floor
#   - Branched trajectory
#   - High dropout  (pi0 = 0.9)
#   - Low signal    (slope = 0.5)
#   - Small dataset (20 × 30)
#   - High overdispersion (size = 0.2)
#
# NEGATIVE CONTROL fixtures  (shuffled pseudotime — signal destroyed)
#   - ZINB linear, shuffled pseudotime    D/O/E → floor
#   - Branched ZINB, shuffled pseudotime  D/O/E → floor (per branch)
#
# Per-scenario plots (a–e):
#   a  UMAP embedding coloured by pseudotime
#   b  Gene expression trends along pseudotime
#   c  D metric base-R plot
#   d  Order consistency (O) — isotonic fits
#   e  Endpoints validity (E) — precision@k with GMM labels
#   f  [endpoint fixture only] GMM labels vs ground-truth overlay
#
# Summary plots:
#   00_metric_heatmap.pdf   — metric heatmap across all non-branched scenarios
#   00_metric_summary.csv   — numeric table of all D / O / E values
#   00_branched_summary.csv — per-branch metrics for branched scenarios
#
# Usage (from repo root):
#   Rscript manuscript/visualize_fixtures.R
#   # or interactively:
#   source("manuscript/visualize_fixtures.R"); run_all_visualizations()
# =============================================================================

library(ggplot2)
library(reshape2)

# ── Locate repo root & source fixtures ───────────────────────────────────────
.repo_root <- tryCatch({
  # When sourced from RStudio, go up one level from /manuscript/
  src <- rstudioapi::getSourceEditorContext()$path
  normalizePath(file.path(dirname(src), ".."))
}, error = function(e) {
  # When run via Rscript or getwd()-based sourcing, use working directory
  getwd()
})

source(file.path(.repo_root, "manuscript", "scripts", "synthetic", "helper-fixtures.R"))

# ── Load BioTrajX ──────────────────────────────────────────────────────────
library(BioTrajX)

# =============================================================================
# HELPER: GMM LABELS FROM FIXTURE
# =============================================================================

#' Compute GMM-based naive / terminal cell labels from a linear fixture
#'
#' @param f     A fixture list (output from make_*_fixture())
#' @param plot  Logical; whether to emit base-R density plots (default FALSE)
#' @return      An "endpoints_validity" object from metrics_e()
compute_gmm_labels <- function(f, plot = FALSE) {
  naive_scores <- colMeans(f$expr[f$early_markers, , drop = FALSE])
  term_scores  <- colMeans(f$expr[f$term_markers,  , drop = FALSE])
  suppressWarnings(
    metrics_e(f$pseudotime,
              early_marker_scores    = naive_scores,
              terminal_marker_scores = term_scores,
              method = "gmm",
              plot   = plot)
  )
}

#' Ground-truth overlay for the endpoint fixture
#'
#' Scatter of per-cell marker score vs pseudotime, with points coloured by
#' whether the cell is a TRUE endpoint (ground truth) vs GMM-labelled endpoint.
#' Shows agreement and any mismatches between the two.
#'
#' @param f     An endpoint fixture (must have $true_early_endpoint /
#'              $true_terminal_endpoint fields)
#' @param E_res An endpoints_validity object from metrics_e()
plot_endpoint_ground_truth <- function(f, E_res) {
  nm <- f$early_markers; tm <- f$term_markers
  pt <- f$pseudotime

  naive_scores <- colMeans(f$expr[nm, , drop = FALSE])
  term_scores  <- colMeans(f$expr[tm, , drop = FALSE])
  cells        <- names(pt)

  gmm_naive <- cells[!is.na(E_res$early_labels)    & E_res$early_labels    == 1]
  gmm_term  <- cells[!is.na(E_res$terminal_labels) & E_res$terminal_labels == 1]

  classify <- function(cell, true_ep, gmm_ep) {
    is_true <- cell %in% true_ep
    is_gmm  <- cell %in% gmm_ep
    dplyr_like <- ifelse(is_true & is_gmm,  "True positive",
                  ifelse(is_true & !is_gmm, "False negative",
                  ifelse(!is_true & is_gmm, "False positive",
                                            "True negative")))
    factor(dplyr_like,
           levels = c("True positive", "False negative",
                      "False positive", "True negative"))
  }

  df <- rbind(
    data.frame(
      Pseudotime = pt,
      Score      = naive_scores,
      Program    = "Naive",
      Label      = classify(cells, f$true_early_endpoint, gmm_naive)
    ),
    data.frame(
      Pseudotime = pt,
      Score      = term_scores,
      Program    = "Terminal",
      Label      = classify(cells, f$true_terminal_endpoint, gmm_term)
    )
  )

  pal <- c("True positive"  = "#2ca02c",
           "False negative" = "#d62728",
           "False positive" = "#ff7f0e",
           "True negative"  = "grey80")

  ggplot(df, aes(x = Pseudotime, y = Score, colour = Label)) +
    geom_point(size = 1.8, alpha = 0.85) +
    scale_colour_manual(values = pal, drop = FALSE) +
    facet_wrap(~ Program, ncol = 2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      title    = "Endpoint fixture: GMM labels vs ground truth",
      subtitle = sprintf("E_early = %.3f  |  E_term = %.3f  |  E_comp = %.3f",
                         E_res$E_early, E_res$E_term,
                         ifelse(is.na(E_res$E_comp), NA, E_res$E_comp)),
      x = "Pseudotime", y = "Mean marker score", colour = NULL
    )
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' UMAP embedding of fixture cells, coloured by pseudotime
#' @param f     Fixture list with $expr, $pseudotime, and optionally $cluster_labels
#' @param title Plot title (NULL = no title)
plot_umap <- function(f, title = NULL) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    message("  [a] UMAP skipped — install uwot: install.packages('uwot')")
    return(NULL)
  }
  mat <- t(f$expr)   # cells × genes
  set.seed(42L)
  umap_res <- uwot::umap(mat,
                         n_neighbors = min(15L, nrow(mat) - 1L),
                         min_dist    = 0.1,
                         n_threads   = 1L,
                         verbose     = FALSE)
  cells <- colnames(f$expr)
  pt    <- f$pseudotime[cells]
  df    <- data.frame(UMAP1 = umap_res[, 1], UMAP2 = umap_res[, 2],
                      Pseudotime = pt)

  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = Pseudotime)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_colour_distiller(palette = "Spectral", direction = 1) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right") +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", colour = "Pseudotime")

  if (!is.null(f$cluster_labels)) {
    df$Branch <- f$cluster_labels[cells]
    ctr <- stats::aggregate(cbind(UMAP1, UMAP2) ~ Branch, data = df, FUN = mean)
    p <- p + ggplot2::geom_label(
      data        = ctr,
      mapping     = aes(x = UMAP1, y = UMAP2, label = Branch),
      colour      = "black", fill = "white",
      size        = 3.5, fontface = "bold",
      inherit.aes = FALSE
    )
  }
  p
}

#' Gene expression trends along pseudotime
#' @param f       Fixture list
#' @param n_genes Number of genes per type to show
#' @param title   Plot title
plot_gene_trends <- function(f, n_genes = 3, title = "Gene Expression Trends") {
  pt   <- f$pseudotime
  expr <- f$expr

  nm <- intersect(f$early_markers, rownames(expr))
  tm <- intersect(f$term_markers,  rownames(expr))
  noise_genes <- setdiff(rownames(expr), c(nm, tm))

  genes <- c(
    head(nm,          min(n_genes, length(nm))),
    head(tm,          min(n_genes, length(tm))),
    head(noise_genes, min(n_genes, length(noise_genes)))
  )
  if (length(genes) == 0) return(NULL)

  df_list <- lapply(genes, function(g) {
    data.frame(
      Pseudotime = pt,
      Expression = expr[g, ],
      Gene       = g,
      GeneType   = ifelse(g %in% nm, "Naive",
                          ifelse(g %in% tm, "Terminal", "Noise"))
    )
  })
  df          <- do.call(rbind, df_list)
  df$GeneType <- factor(df$GeneType, levels = c("Naive", "Terminal", "Noise"))

  ggplot(df, aes(x = Pseudotime, y = Expression, color = Gene)) +
    geom_point(alpha = 0.3, size = 0.7) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, span = 0.75) +
    facet_wrap(~ GeneType, ncol = 3) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(title = title, x = "Pseudotime", y = "Expression")
}

#' Heatmap of D / O / E metric values across scenarios
#' @param summary_df  Data frame with columns Scenario, D_naive, D_term, O, E_naive, E_term, E_comp
#' @param title       Optional plot title
plot_metric_heatmap <- function(summary_df, title = NULL) {
  metric_cols <- c("D_early", "D_term", "O", "E_early", "E_term", "DOE_score")
  keep <- intersect(metric_cols, names(summary_df))
  # Carry O_orientation as an id variable if present
  id_vars <- if ("O_orientation" %in% names(summary_df)) c("Scenario", "O_orientation") else "Scenario"
  plot_df <- summary_df[, c(id_vars, keep), drop = FALSE]
  df_long <- melt(plot_df, id.vars = id_vars,
                  variable.name = "Metric", value.name = "Value")
  df_long$Metric   <- factor(df_long$Metric, levels = metric_cols)
  df_long$Scenario <- factor(df_long$Scenario,
                             levels = rev(unique(summary_df$Scenario)))
  # Build cell label: O gets orientation on second line
  if ("O_orientation" %in% names(df_long)) {
    df_long$cell_label <- ifelse(
      is.na(df_long$Value), "\u2013",
      ifelse(
        df_long$Metric == "O" & !is.na(df_long$O_orientation),
        sprintf("%.2f\n(%s)", df_long$Value, df_long$O_orientation),
        sprintf("%.2f", df_long$Value)
      )
    )
  } else {
    df_long$cell_label <- ifelse(is.na(df_long$Value), "\u2013", sprintf("%.2f", df_long$Value))
  }

  ggplot(df_long, aes(x = Metric, y = Scenario, fill = Value)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = cell_label),
              size = 3, colour = "black", lineheight = 0.9) +
    scale_fill_gradient2(
      low      = "#d95f4f",
      mid      = "#f3e96b",
      high     = "#5fbf7a",
      midpoint = 0.5,
      limits   = c(0, 1),
      na.value = "grey90",
      name     = "Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid  = element_blank()) +
    labs(title = title, x = NULL, y = NULL)
}

# =============================================================================
# METRIC WRAPPERS (base-R plots → PNG/PDF via png()/pdf() / dev.off())
# =============================================================================

#' Save a base-R plotting call to PNG or PDF
save_base_plot <- function(plot_fn, path, width = 8, height = 5, res = 150, ...) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "pdf") {
    grDevices::pdf(path, width = width, height = height)
  } else {
    grDevices::png(path, width = width, height = height, units = "in", res = res)
  }
  tryCatch(plot_fn(...), finally = grDevices::dev.off())
  invisible(path)
}

# =============================================================================
# CROSS-SCENARIO COMPARISON
# =============================================================================

#' Compute a summary data frame of D / O / E metrics for a list of scenarios
compute_metric_summary <- function(scenarios) {
  rows <- lapply(scenarios, function(sc) {
    f <- sc$fixture
    if (is.null(f$early_markers) || is.null(f$term_markers)) return(NULL)

    nm <- intersect(f$early_markers, rownames(f$expr))
    tm <- intersect(f$term_markers,  rownames(f$expr))
    if (length(nm) == 0 || length(tm) == 0) return(NULL)

    # D metric
    D_res <- tryCatch(
      metrics_d(f$expr, nm, tm, f$pseudotime),
      error = function(e) NULL
    )

    # O metric
    O_res <- tryCatch(
      metrics_o(f$expr, nm, tm, f$pseudotime),
      error = function(e) NULL
    )

    # E metric (GMM)
    E_res <- tryCatch(
      compute_gmm_labels(f),
      error = function(e) NULL
    )

    D_early <- if (!is.null(D_res)) D_res$D_early else NA
    D_term  <- if (!is.null(D_res)) D_res$D_term  else NA
    O       <- if (!is.null(O_res)) O_res$O       else NA
    E_early <- if (!is.null(E_res)) E_res$E_early else NA
    E_term  <- if (!is.null(E_res)) E_res$E_term  else NA
    E_comp  <- if (!is.null(E_res)) E_res$E_comp  else NA
    D_comp    <- mean(c(D_early, D_term), na.rm = TRUE)
    if (is.nan(D_comp)) D_comp <- NA_real_
    comp_vec  <- c(D_comp, O, E_comp)
    DOE_score <- if (all(is.na(comp_vec))) NA_real_ else mean(comp_vec, na.rm = TRUE)

    data.frame(
      Scenario      = sc$label,
      D_early       = D_early,
      D_term        = D_term,
      O             = O,
      O_orientation = if (!is.null(O_res)) O_res$orientation else NA_character_,
      E_early       = E_early,
      E_term        = E_term,
      DOE_score     = DOE_score,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

#' Compute per-branch D / O / E metrics for branched fixtures
#'
#' For each branch, subsets cells and markers, then computes all three metrics
#' on that branch's cells only.
#'
#' @param scenarios Full scenario list; only type=="branched" entries are used
#' @return A data frame with one row per branch per branched scenario
compute_branched_metric_summary <- function(scenarios) {
  branched <- Filter(function(sc) sc$type %in% c("branched", "branched_endpoint"), scenarios)

  rows <- lapply(branched, function(sc) {
    f    <- sc$fixture
    # Iterate over trajectory pairs (AB, AC) defined in branch_filters
    trajs <- names(f$branch_filters)

    traj_rows <- lapply(trajs, function(traj) {
      include_labels <- f$branch_filters[[traj]]$include
      cells <- names(f$cluster_labels)[f$cluster_labels %in% include_labels]
      if (length(cells) == 0) return(NULL)

      nm <- intersect(f$early_markers_list[[traj]], rownames(f$expr))
      tm <- intersect(f$term_markers_list[[traj]],  rownames(f$expr))
      if (length(nm) == 0 || length(tm) == 0) return(NULL)

      expr_tr <- f$expr[, cells, drop = FALSE]
      pt_tr   <- f$pseudotime[cells]

      # D metric
      D_res <- tryCatch(
        metrics_d(expr_tr, nm, tm, pt_tr),
        error = function(e) NULL
      )

      # O metric
      O_res <- tryCatch(
        metrics_o(expr_tr, nm, tm, pt_tr),
        error = function(e) NULL
      )

      # E metric (GMM) — on trajectory cells only
      naive_scores <- colMeans(expr_tr[nm, , drop = FALSE])
      term_scores  <- colMeans(expr_tr[tm, , drop = FALSE])
      E_res <- tryCatch(
        suppressWarnings(
          metrics_e(pt_tr,
                    early_marker_scores    = naive_scores,
                    terminal_marker_scores = term_scores,
                    method = "gmm", plot = FALSE)
        ),
        error = function(e) NULL
      )

      D_early <- if (!is.null(D_res)) D_res$D_early else NA
      D_term  <- if (!is.null(D_res)) D_res$D_term  else NA
      O_val   <- if (!is.null(O_res)) O_res$O       else NA
      E_early <- if (!is.null(E_res)) E_res$E_early else NA
      E_term  <- if (!is.null(E_res)) E_res$E_term  else NA
      E_comp  <- if (!is.null(E_res)) E_res$E_comp  else NA
      D_comp    <- mean(c(D_early, D_term), na.rm = TRUE)
      if (is.nan(D_comp)) D_comp <- NA_real_
      comp_vec  <- c(D_comp, O_val, E_comp)
      DOE_score <- if (all(is.na(comp_vec))) NA_real_ else mean(comp_vec, na.rm = TRUE)

      data.frame(
        Scenario      = sc$label,
        Branch        = traj,
        n_cells       = length(cells),
        D_early       = D_early,
        D_term        = D_term,
        O             = O_val,
        O_orientation = if (!is.null(O_res)) O_res$orientation else NA_character_,
        E_early       = E_early,
        E_term        = E_term,
        DOE_score     = DOE_score,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, Filter(Negate(is.null), traj_rows))
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

#' Run all visualizations and save to output_dir
#'
#' @param output_dir  Directory to save plots (created if absent)
#' @param format      "png" or "pdf"
#' @param width_heatmap,height_heatmap  Dimensions for heatmap panels
#' @param width_std,height_std          Dimensions for standard panels
run_all_visualizations <- function(output_dir   = file.path(.repo_root, "manuscript/plots"),
                                   format       = "pdf",
                                   width_heatmap = 11,
                                   height_heatmap = 6,
                                   width_std    = 8,
                                   height_std   = 5) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", output_dir)

  # ── Define scenarios ───────────────────────────────────────────────────────
  scenarios <- list(
    # ── Ceiling / correctness (Gaussian noise — known expected values) ──────────
    list(name = "01_gaussian_linear",       label = "Gaussian Linear",
         fixture = make_gaussian_fixture(n_cells = 1000),
         type = "linear",   group = "ceiling"),
    list(name = "02_gaussian_reversed",     label = "Gaussian Reversed",
         fixture = make_gaussian_reversed_fixture(n_cells = 1000),
         type = "linear",   group = "ceiling"),
    list(name = "03_stepwise_repression",   label = "Stepwise Repression",
         fixture = make_stepwise_repression_fixture(n_cells = 1000),
         type = "endpoint", group = "ceiling"),
    list(name = "04_stepwise_activation",   label = "Stepwise Activation",
         fixture = make_stepwise_activation_fixture(n_cells = 1000),
         type = "endpoint", group = "ceiling"),
    # ── ZINB stress / robustness (linear) ──────────────────────────────────────
    list(name = "05_zinb_linear",           label = "ZINB Linear",
         fixture = make_linear_fixture(n_cells = 1000),
         type = "linear",   group = "stress"),
    list(name = "06_zinb_reversed",         label = "ZINB Reversed",
         fixture = make_reversed_fixture(n_cells = 1000),
         type = "linear",   group = "stress"),
    list(name = "07_flat",                  label = "Flat / No Signal",
         fixture = make_flat_fixture(n_cells = 1000),
         type = "linear",   group = "stress"),
    list(name = "08_high_dropout",          label = "High Dropout (pi0=0.9)",
         fixture = make_high_dropout_fixture(n_cells = 1000),
         type = "linear",   group = "stress"),
    list(name = "09_low_signal",            label = "Low Signal (slope=0.5)",
         fixture = make_low_signal_fixture(n_cells = 1000),
         type = "linear",   group = "stress"),
    # ── Branched trajectories ───────────────────────────────────────────────────
    list(name = "10_gaussian_branched",     label = "Gaussian Branched",
         fixture = make_gaussian_branched_fixture(n_cells_per_branch = 500),
         type = "branched", group = "branched"),
    list(name = "11_gaussian_branched_reverse", label = "Gaussian Branched Reversed",
         fixture = make_gaussian_branched_reversed_fixture(n_cells_per_branch = 500),
         type = "branched", group = "branched"),
    # ── Branched stepwise endpoint fixtures (ceiling / floor for E metric) ──────
    list(name = "12_branched_stepwise_repression", label = "Branched Stepwise Repression",
         fixture = make_branched_stepwise_fixture(n_cells_per_branch = 500),
         type = "branched_endpoint", group = "branched"),
    list(name = "13_branched_stepwise_activation", label = "Branched Stepwise Activation",
         fixture = make_branched_stepwise_reversed_fixture(n_cells_per_branch = 500),
         type = "branched_endpoint", group = "branched"),
    # ── Branched stress / robustness fixtures ───────────────────────────────────
    list(name = "14_branched_ZINB",         label = "Branched ZINB",
         fixture = make_branched_zinb_fixture(n_cells_per_branch = 500),
         type = "branched", group = "branched"),
    list(name = "15_branched_low_signal",   label = "Branched Low Signal",
         fixture = make_branched_low_signal_fixture(n_cells_per_branch = 500),
         type = "branched", group = "branched"),
    # ── Negative controls (shuffled pseudotime — signal destroyed) ──────────────
    list(name = "16_zinb_linear_shuffled", label = "ZINB Linear (Shuffled Pseudotime)",
         fixture = make_shuffled_fixture(n_cells = 1000),
         type = "linear",   group = "negative_control"),
    list(name = "17_branched_zinb_shuffled", label = "Branched ZINB (Shuffled Pseudotime)",
         fixture = make_branched_shuffled_fixture(n_cells_per_branch = 500),
         type = "branched", group = "negative_control")
  )

  # ── Per-scenario plots ─────────────────────────────────────────────────────
  for (sc in scenarios) {
    f   <- sc$fixture
    nm  <- sc$name
    lb  <- sc$label
    typ <- sc$type
    message(sprintf("\n[%s]  %s", nm, lb))

    # --- a. UMAP (always all cells) ---------------------------------------------
    message("  [a] UMAP")
    f_umap <- f
    if (typ %in% c("branched", "branched_endpoint")) {
      f_umap$early_markers <- f$early_markers_list$AB
      f_umap$term_markers  <- f$term_markers_list$AB
    }
    p_umap <- plot_umap(f_umap, title = NULL)
    if (!is.null(p_umap)) {
      ggsave(
        file.path(output_dir, paste0(nm, "_a_umap.pdf")),
        p_umap, width = width_std, height = height_std
      )
    }

    if (typ %in% c("branched", "branched_endpoint")) {
      # ── Branched: generate panels b–f separately for each branch (AB, AC) ────
      branches <- names(f$branch_filters)

      for (br in branches) {
        include_labels <- f$branch_filters[[br]]$include
        br_cells <- names(f$cluster_labels)[f$cluster_labels %in% include_labels]

        f_br <- f
        f_br$expr           <- f$expr[, br_cells, drop = FALSE]
        f_br$pseudotime     <- f$pseudotime[br_cells]
        f_br$cluster_labels <- f$cluster_labels[br_cells]
        f_br$early_markers  <- f$early_markers_list[[br]]
        f_br$term_markers   <- f$term_markers_list[[br]]

        if (typ == "branched_endpoint") {
          f_br$true_early_endpoint    <- f$true_early_endpoint
          fate <- substr(br, nchar(br), nchar(br))   # "B" from "AB", "C" from "AC"
          f_br$true_terminal_endpoint <- f[[paste0("true_terminal_endpoint_", fate)]]
        }

        nm_br <- intersect(f_br$early_markers, rownames(f_br$expr))
        tm_br <- intersect(f_br$term_markers,  rownames(f_br$expr))

        # b. Gene trends (branch cells + branch markers — no cross-branch overlap)
        message(sprintf("  [b] Gene trends (%s)", br))
        p_trends <- plot_gene_trends(f_br, n_genes = 3, title = NULL)
        if (!is.null(p_trends)) {
          ggsave(
            file.path(output_dir, paste0(nm, "_b_gene_trends_", br, ".pdf")),
            p_trends, width = width_heatmap, height = 4
          )
        }

        if (length(nm_br) == 0 || length(tm_br) == 0) next

        # c. D metric
        message(sprintf("  [c] D metric (%s)", br))
        D_res <- tryCatch(
          metrics_d(f_br$expr, nm_br, tm_br, f_br$pseudotime),
          error = function(e) { message("    D failed: ", e$message); NULL }
        )
        if (!is.null(D_res)) {
          .pt <- f_br$pseudotime; .D <- D_res
          save_base_plot(
            function() plot_metrics_d(.pt, .D, main = ""),
            path   = file.path(output_dir, paste0(nm, "_c_D_metric_", br, ".pdf")),
            width  = width_std, height = height_std
          )
        }

        # d. O metric
        message(sprintf("  [d] O metric (%s)", br))
        O_res <- tryCatch(
          metrics_o(f_br$expr, nm_br, tm_br, f_br$pseudotime),
          error = function(e) { message("    O failed: ", e$message); NULL }
        )
        if (!is.null(O_res) && !is.na(O_res$O)) {
          .expr <- f_br$expr; .pt <- f_br$pseudotime; .O <- O_res
          .nm <- nm_br; .tm <- tm_br
          save_base_plot(
            function() plot_metrics_o(.expr, .pt, .O,
                                      early_markers = .nm, terminal_markers = .tm,
                                      main = ""),
            path   = file.path(output_dir, paste0(nm, "_d_O_metric_", br, ".pdf")),
            width  = width_std, height = height_std
          )
        }

        # e. E metric
        message(sprintf("  [e] E metric (%s)", br))
        E_res <- tryCatch(
          compute_gmm_labels(f_br),
          error = function(e) { message("    E failed: ", e$message); NULL }
        )
        if (!is.null(E_res)) {
          .E <- E_res
          save_base_plot(
            function() plot_metrics_e(.E, overlay = TRUE, main = ""),
            path   = file.path(output_dir, paste0(nm, "_e_E_metric_", br, ".pdf")),
            width  = width_std, height = height_std
          )
        }

        # f. Ground-truth overlay (branched_endpoint only)
        if (typ == "branched_endpoint" &&
            !is.null(E_res) && !is.null(f_br$true_early_endpoint)) {
          message(sprintf("  [f] Ground-truth overlay (%s)", br))
          p_gt <- tryCatch(
            plot_endpoint_ground_truth(f_br, E_res),
            error = function(e) { message("    GT plot failed: ", e$message); NULL }
          )
          if (!is.null(p_gt)) {
            ggsave(
              file.path(output_dir, paste0(nm, "_f_ground_truth_", br, ".pdf")),
              p_gt, width = width_std + 2, height = height_std
            )
          }
        }
      }  # end branch loop

    } else {
      # ── Non-branched: single set of panels ────────────────────────────────────
      nm_genes    <- intersect(f$early_markers, rownames(f$expr))
      tm_genes    <- intersect(f$term_markers,  rownames(f$expr))
      has_markers <- length(nm_genes) > 0 && length(tm_genes) > 0

      # b. Gene expression trends
      message("  [b] Gene expression trends")
      p_trends <- plot_gene_trends(f, n_genes = 3, title = NULL)
      if (!is.null(p_trends)) {
        ggsave(
          file.path(output_dir, paste0(nm, "_b_gene_trends.pdf")),
          p_trends, width = width_heatmap, height = 4
        )
      }

      if (!has_markers) {
        message("  (skipping metric plots — no marker overlap)")
        next
      }

      # c. D metric
      message("  [c] D metric")
      D_res <- tryCatch(
        metrics_d(f$expr, nm_genes, tm_genes, f$pseudotime),
        error = function(e) { message("    D failed: ", e$message); NULL }
      )
      if (!is.null(D_res)) {
        save_base_plot(
          function() plot_metrics_d(f$pseudotime, D_res, main = ""),
          path   = file.path(output_dir, paste0(nm, "_c_D_metric.pdf")),
          width  = width_std, height = height_std
        )
      }

      # d. O metric
      message("  [d] O metric")
      O_res <- tryCatch(
        metrics_o(f$expr, nm_genes, tm_genes, f$pseudotime),
        error = function(e) { message("    O failed: ", e$message); NULL }
      )
      if (!is.null(O_res) && !is.na(O_res$O)) {
        save_base_plot(
          function() plot_metrics_o(f$expr, f$pseudotime, O_res,
                                    early_markers    = nm_genes,
                                    terminal_markers = tm_genes,
                                    main = ""),
          path   = file.path(output_dir, paste0(nm, "_d_O_metric.pdf")),
          width  = width_std, height = height_std
        )
      }

      # e. E metric
      message("  [e] E metric (GMM)")
      E_res <- tryCatch(
        compute_gmm_labels(f),
        error = function(e) { message("    E failed: ", e$message); NULL }
      )
      if (!is.null(E_res)) {
        save_base_plot(
          function() plot_metrics_e(E_res, overlay = TRUE, main = ""),
          path   = file.path(output_dir, paste0(nm, "_e_E_metric.pdf")),
          width  = width_std, height = height_std
        )
      }

      # f. Ground-truth overlay (endpoint fixtures only)
      if (typ == "endpoint" &&
          !is.null(E_res) && !is.null(f$true_early_endpoint)) {
        message("  [f] Ground-truth overlay")
        p_gt <- tryCatch(
          plot_endpoint_ground_truth(f, E_res),
          error = function(e) { message("    GT plot failed: ", e$message); NULL }
        )
        if (!is.null(p_gt)) {
          ggsave(
            file.path(output_dir, paste0(nm, "_f_ground_truth.pdf")),
            p_gt, width = width_std + 2, height = height_std
          )
        }
      }
    }  # end non-branched
  }

  # ── Summary: metric heatmap (all non-branched scenarios) ───────────────────
  summarisable <- Filter(function(sc) sc$type %in% c("linear", "endpoint"), scenarios)

  summary_df <- tryCatch(
    compute_metric_summary(summarisable),
    error = function(e) { message("  Summary failed: ", e$message); NULL }
  )

  if (!is.null(summary_df) && nrow(summary_df) > 0) {
    # Add group column
    group_map  <- setNames(sapply(summarisable, `[[`, "group"),
                           sapply(summarisable, `[[`, "label"))
    summary_df$Group <- group_map[summary_df$Scenario]

    # Metric heatmap (replaces all bar-chart comparison figures)
    message("\n[summary] Metric heatmap")
    p_heat <- plot_metric_heatmap(summary_df)
    ggsave(
      file.path(output_dir, "00_metric_heatmap.pdf"),
      p_heat,
      width  = 8,
      height = max(4, nrow(summary_df) * 0.55 + 1.5)
    )

    # Summary CSV
    write.csv(summary_df[, setdiff(names(summary_df), "Group")],
              file.path(output_dir, "00_metric_summary.csv"),
              row.names = FALSE)
    message("  Saved 00_metric_summary.csv")
  }

  # ── Branched summary CSV ────────────────────────────────────────────────────
  message("\n[summary] Branched metric summary (per branch)")
  branched_df <- tryCatch(
    compute_branched_metric_summary(scenarios),
    error = function(e) { message("  Branched summary failed: ", e$message); NULL }
  )
  if (!is.null(branched_df) && nrow(branched_df) > 0) {
    write.csv(branched_df,
              file.path(output_dir, "00_branched_summary.csv"),
              row.names = FALSE)
    message("  Saved 00_branched_summary.csv")
    print(branched_df)

    # Branched metric heatmap (AB branch only)
    message("\n[summary] Branched metric heatmap")
    branched_heat_df <- branched_df[branched_df$Branch == "AB", , drop = FALSE]
    branched_heat_df$Scenario <- paste0(branched_heat_df$Scenario, " (", branched_heat_df$Branch, ")")
    branched_heat_df$Branch   <- NULL
    branched_heat_df$n_cells  <- NULL
    p_branched_heat <- plot_metric_heatmap(branched_heat_df)
    ggsave(
      file.path(output_dir, "00_branched_metric_heatmap.pdf"),
      p_branched_heat,
      width  = 8,
      height = max(4, nrow(branched_heat_df) * 0.55 + 1.5)
    )
    message("  Saved 00_branched_metric_heatmap.pdf")
  }

  message(sprintf("\nDone! %d files written to: %s",
                  length(list.files(output_dir)), output_dir))
  invisible(output_dir)
}

# =============================================================================
# ENTRY POINT
# =============================================================================

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) run_all_visualizations(output_dir = args[1])
  else run_all_visualizations()
}
