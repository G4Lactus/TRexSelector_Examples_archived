#!/usr/bin/env Rscript
# ==============================================================================
# trex_spca_plt_utils.R
# ==============================================================================
#
# Plotting utilities for the T-Rex SPCA demonstration suite — ggplot2 port of
# the Python plotter of the same name (trex_spca_plt_utils.py), kept visually
# aligned with it: tab10 method colours, hollow per-method markers, ONE legend
# outside the panels, PNG + PDF output into the sibling
# simulation_results/plots/ folder.
#
# Reads the tidy result CSVs written by the demo:
#
#   method,metric,snr_db,value  — SNR sweep (metrics FDR/TPR/PEV/PEVsig/
#                                 PEVsigmix). Renders the TPR|FDR sweep figure
#                                 and, when PEVsig is present, the extra
#                                 two-normalization PEV figure (<stem>_pev).
#   pc,metric,snr_db,value      — union-support FDR per PC (metrics
#                                 FDRunion/k). Renders the two-panel
#                                 union-FDR figure.
#
# The sweep axis is the signal-to-noise ratio in decibels and is plotted as
# such (linear in dB, integer-dB major ticks, one minor tick per dB).
#
# Usage:
#   Rscript trex_spca_plt_utils.R <results.csv> [--outdir DIR] [--title T]
#                                 [--pev-title T] [--tfdr F] [--formats png,pdf]
#                                 [--stem S]
# ==============================================================================

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

DEFAULT_TFDR <- 0.1

# matplotlib "tab10" — the colour convention shared by every suite's figures.
TAB10 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

# Hollow per-method markers (o, s, ^, D, v, P) — overlapping curves stay
# distinguishable. ggplot2 shapes 21-25 honour fill; 3 is the plus.
METHOD_SHAPES <- c(21L, 22L, 24L, 23L, 25L, 3L)

# Baselines first, then the T-Rex SPCA variants, so the legend reads top-down.
METHOD_ORDER <- c("OrdPCA", "OraclePCA",
                  "TRexSPCA-EN-Act", "TRexSPCA-EN-Thr",
                  "TRexSPCA-ENaug-Act", "TRexSPCA-ENaug-Thr")


# ------------------------------------------------------------------------------
# Shared helpers (house conventions)
# ------------------------------------------------------------------------------

.order_methods <- function(present) {
  c(METHOD_ORDER[METHOD_ORDER %in% present],
    present[!present %in% METHOD_ORDER])
}

#' Integer-dB x-axis: major ticks on whole decibels at the coarsest step that
#' keeps <= 8 major intervals, one minor tick per dB, 4% padding.
.db_x_scale <- function(values) {
  lo <- min(values)
  hi <- max(values)
  span <- hi - lo
  pad <- if (span > 0) span * 0.04 else 0.5
  steps <- c(1, 2, 5, 10, 20, 25, 50, 100)
  step <- steps[which(span <= 8 * steps)[1]]
  if (is.na(step)) step <- max(1, round(span / 8))
  scale_x_continuous(
    breaks = seq(ceiling(lo / step) * step, hi, by = step),
    minor_breaks = seq(ceiling(lo), floor(hi), by = 1),
    limits = c(lo - pad, hi + pad),
    expand = expansion(0)
  )
}

#' The shared panel look: light major/minor grid, bold panel title.
.panel_theme <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25,
                                      linetype = "dashed"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_rect(colour = "grey55", linewidth = 0.3)
    )
}

#' Method colour/fill/shape scales over a fixed ordered set.
.method_scales <- function(methods, group_label = "Method") {
  cols <- setNames(TAB10[seq_along(methods)], methods)
  shps <- setNames(METHOD_SHAPES[(seq_along(methods) - 1L) %%
                                   length(METHOD_SHAPES) + 1L], methods)
  list(
    scale_color_manual(name = group_label, values = cols),
    scale_shape_manual(name = group_label, values = shps)
  )
}

save_figure <- function(plot, outdir, stem, formats = c("png", "pdf"),
                        width = 12.4, height = 4.8) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (ext in formats) {
    path <- file.path(outdir, paste0(stem, ".", ext))
    ggsave(path, plot = plot, width = width, height = height,
           dpi = 300, bg = "white")
    cat(sprintf("Wrote: %s\n", path))
  }
}

#' `.../simulation_results/data/x.csv` -> `.../simulation_results/plots`.
default_plots_dir <- function(csv_path) {
  parent <- dirname(csv_path)
  if (basename(parent) == "data") file.path(dirname(parent), "plots") else file.path(parent, "plots")
}

load_csv <- function(path) {
  if (!file.exists(path)) stop("Input CSV not found: ", path)
  df <- utils::read.csv(path, strip.white = TRUE, stringsAsFactors = FALSE)
  names(df) <- trimws(names(df))
  df
}


# ------------------------------------------------------------------------------
# Figure: TPR / FDR vs. SNR (dB)
# ------------------------------------------------------------------------------

.sweep_panel <- function(df, methods, metric, panel_title, ylim = NULL,
                         tfdr = 0, ybreaks_step = 0.1,
                         group_label = "Method", percent = FALSE) {
  sub <- df[df$metric == metric & df$method %in% methods, ]
  sub$method <- factor(sub$method, levels = methods)
  if (percent) sub$value <- 100 * sub$value

  p <- ggplot(sub, aes(x = snr_db, y = value,
                       color = method, shape = method)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2.4, fill = NA, stroke = 0.9)

  if (tfdr > 0) {
    p <- p + geom_hline(
      aes(yintercept = tfdr, linetype = sprintf("tFDR = %g", tfdr)),
      color = "black", linewidth = 0.8
    ) +
      scale_linetype_manual(name = NULL, values = "dotted",
                            guide = guide_legend(order = 1))
  }

  p <- p + .method_scales(methods, group_label) +
    guides(color = guide_legend(order = 2), shape = guide_legend(order = 2)) +
    labs(title = panel_title, x = "SNR (dB)", y = NULL) +
    .db_x_scale(sort(unique(sub$snr_db))) +
    scale_y_continuous(breaks = scales::breaks_width(ybreaks_step),
                       minor_breaks = scales::breaks_width(ybreaks_step / 2)) +
    .panel_theme()
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

plot_spca_sweep <- function(df, title, tfdr = DEFAULT_TFDR) {
  methods <- .order_methods(unique(df$method))
  p_tpr <- .sweep_panel(df, methods, "TPR", "TPR (true positive rate)",
                        ylim = c(0, 1.05))
  p_fdr <- .sweep_panel(df, methods, "FDR", "FDR (false discovery rate)",
                        tfdr = tfdr)
  (p_tpr | p_fdr) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(
                      size = 14, face = "bold", hjust = 0.5
                    )))
}


# ------------------------------------------------------------------------------
# Extra figure: the two PEV normalizations side by side
# ------------------------------------------------------------------------------

plot_spca_pev_normalizations <- function(df, title) {
  methods <- .order_methods(unique(df$method))

  p_total <- .sweep_panel(df, methods, "PEV",
                          "Normalized by the total variance of X",
                          ybreaks_step = 10, percent = TRUE) +
    labs(y = "Cumulative PEV (%)")

  p_sig <- .sweep_panel(df, methods, "PEVsig",
                        "Normalized by Signal + Mixed EV (Definition 1)",
                        ybreaks_step = 10, percent = TRUE) +
    labs(y = "Cumulative PEV (%)")

  # Everything above 100% is Null EV — variance harvested from null loadings
  # and reported as if it were explained signal.
  ref_layers <- list(
    geom_hline(aes(yintercept = 100,
                   linetype = "Signal + Mixed EV = 100%"),
               color = "black", linewidth = 0.8)
  )

  # The paper's "Ordinary PCA (Sig + Mix)" dashed reference curve.
  ref_sub <- df[df$method == "OrdPCA" & df$metric == "PEVsigmix", ]
  lt_values <- c("Signal + Mixed EV = 100%" = "dotted")
  if (nrow(ref_sub) > 0) {
    ref_sub$value <- 100 * ref_sub$value
    ref_layers <- c(ref_layers, list(
      geom_line(data = ref_sub,
                aes(x = snr_db, y = value,
                    linetype = "OrdPCA (Sig + Mix part)"),
                color = "grey35", linewidth = 0.7, inherit.aes = FALSE)
    ))
    lt_values <- c(lt_values, "OrdPCA (Sig + Mix part)" = "dashed")
  }
  p_sig <- p_sig + ref_layers +
    scale_linetype_manual(name = NULL, values = lt_values,
                          guide = guide_legend(order = 1))

  (p_total | p_sig) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(
                      size = 14, face = "bold", hjust = 0.5
                    )))
}


# ------------------------------------------------------------------------------
# Figure: union-support FDR per PC (demo 01, section 2)
# ------------------------------------------------------------------------------

plot_union_fdr <- function(df, title, tfdr = DEFAULT_TFDR) {
  df$pc_lbl <- sprintf("PC%g", df$pc)
  pcs <- sprintf("PC%g", sort(unique(df$pc)))

  panel <- function(metric, panel_title, ylab, ylim = NULL, tfdr = 0) {
    sub <- df[df$metric == metric, ]
    sub$pc_lbl <- factor(sub$pc_lbl, levels = pcs)
    p <- ggplot(sub, aes(x = snr_db, y = value,
                         color = pc_lbl, shape = pc_lbl)) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 2.4, fill = NA, stroke = 0.9)
    if (tfdr > 0) {
      p <- p + geom_hline(
        aes(yintercept = tfdr, linetype = sprintf("tFDR = %g", tfdr)),
        color = "black", linewidth = 0.8
      ) +
        scale_linetype_manual(name = NULL, values = "dotted",
                              guide = guide_legend(order = 1))
    }
    p <- p + .method_scales(pcs, "Component") +
      guides(color = guide_legend(order = 2),
             shape = guide_legend(order = 2)) +
      labs(title = panel_title, x = "SNR (dB)", y = ylab) +
      .db_x_scale(sort(unique(sub$snr_db))) +
      .panel_theme()
    if (!is.null(ylim)) {
      p <- p + coord_cartesian(ylim = ylim) +
        scale_y_continuous(breaks = scales::breaks_width(0.1),
                           minor_breaks = scales::breaks_width(0.05))
    }
    p
  }

  p_fdr <- panel("FDRunion", "Union-support FDR", "FDR (vs. union support)",
                 ylim = c(0, 1), tfdr = tfdr)
  p_k <- panel("k", "Selected support size", "mean |A_m|")

  (p_fdr | p_k) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(
                      size = 14, face = "bold", hjust = 0.5
                    )))
}


# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------

.parse_args <- function(argv) {
  args <- list(csv = character(0), outdir = NULL, title = NULL,
               pev_title = NULL, tfdr = DEFAULT_TFDR,
               formats = c("png", "pdf"), stem = NULL)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    flag <- function() {
      i <<- i + 1L
      argv[[i]]
    }
    switch(a,
           "--outdir"    = {
             args$outdir <- flag()
           },
           "--title"     = {
             args$title <- flag()
           },
           "--pev-title" = {
             args$pev_title <- flag()
           },
           "--tfdr"      = {
             args$tfdr <- as.numeric(flag())
           },
           "--formats"   = {
             args$formats <- strsplit(flag(), ",")[[1]]
           },
           "--stem"      = {
             args$stem <- flag()
           },
           {
             args$csv <- c(args$csv, a)
           })
    i <- i + 1L
  }
  if (length(args$csv) != 1L)
    stop("Usage: trex_spca_plt_utils.R <results.csv> [--outdir DIR] ",
         "[--title T] [--pev-title T] [--tfdr F] [--formats png,pdf] ",
         "[--stem S]")
  args
}

main <- function() {
  args <- .parse_args(commandArgs(trailingOnly = TRUE))
  csv_path <- args$csv[[1]]
  df <- load_csv(csv_path)

  outdir <- if (is.null(args$outdir)) default_plots_dir(csv_path) else args$outdir
  stem <- if (is.null(args$stem)) sub("[.]csv$", "", basename(csv_path)) else args$stem

  # Union-FDR CSV (pc,metric,snr_db,value) gets its own two-panel figure.
  if ("pc" %in% names(df)) {
    fig <- plot_union_fdr(
      df, if (is.null(args$title)) "T-Rex SPCA: union-support FDR per PC" else args$title, args$tfdr
    )
    save_figure(fig, outdir, stem, args$formats)
    return(invisible(0L))
  }

  missing <- setdiff(c("method", "metric", "snr_db", "value"), names(df))
  if (length(missing) > 0)
    stop("CSV is missing columns: ", paste(missing, collapse = ", "))

  fig <- plot_spca_sweep(
    df, if (is.null(args$title)) stem else args$title, args$tfdr
  )
  save_figure(fig, outdir, stem, args$formats)

  # Extra figure — only if the demo reported the second normalization.
  if ("PEVsig" %in% df$metric) {
    fig2 <- plot_spca_pev_normalizations(
      df, if (is.null(args$pev_title))
        "T-Rex SPCA: cumulative PEV under both normalizations" else args$pev_title
    )
    save_figure(fig2, outdir, paste0(stem, "_pev"), args$formats)
  }
  invisible(0L)
}

if (sys.nframe() == 0L) main()
