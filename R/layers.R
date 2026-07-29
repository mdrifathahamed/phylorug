# layers.R
# Layer functions for the phylorug layered plotting system.
# Each function returns a phylorug_layer object containing a draw() closure
# that receives the rendering context (ctx) at plot time.


# ── theme_phylorug ───────────────────────────────────────────────────────────

#' Set global theme
#'
#' Controls plot-wide styling: font family, base size multiplier, and margins.
#' All auto-scaled sizes (tips, labels, dots, legend) are multiplied by
#' `base_size`, so `base_size = 1.5` makes everything 50 percent larger.
#'
#' @param base_size Numeric multiplier applied to all auto-scaled elements.
#'   Default `1`.
#' @param family Font family. Default `"Helvetica"`.
#' @param mar Numeric vector of length 4. Plot margins
#'   (bottom, left, top, right). Default `c(0.5, 0.5, 0.5, 2.5)`.
#'
#' @returns A `phylorug_layer` object.
#' @export
theme_phylorug <- function(base_size = 1,
                           family    = "Helvetica",
                           mar       = c(0.5, 0.5, 0.5, 2.5)) {

  structure(
    list(
      type   = "theme",
      draw   = function(ctx) invisible(NULL),  # applied before drawing
      params = list(
        base_size = base_size,
        family    = family,
        mar       = mar
      )
    ),
    class = "phylorug_layer"
  )
}


# ── rug_layer ────────────────────────────────────────────────────────────────

#' Add rug cells at variable nodes
#'
#' Draws the coloured grid at each internal node where at least one
#' comparison tree disagrees with the backbone topology.
#'
#' @param cell_scale Numeric. Cell height as a fraction of median tip
#'   spacing. Default `0.45`.
#' @param n_rows,n_cols Grid shape. `NULL` for automatic.
#' @param x_offset,y_offset Shift the grid from the node position.
#' @param rug_position `"inside"` (default) or `"outside"`.
#' @param border_col Cell border colour. Default `"grey40"`.
#' @param border_lwd Cell border line width. Default `0.2`.
#'
#' @returns A `phylorug_layer` object.
#' @export
rug_layer <- function(cell_scale   = 0.45,
                      n_rows       = NULL,
                      n_cols       = NULL,
                      x_offset     = 0,
                      y_offset     = 0,
                      rug_position = c("inside", "outside"),
                      border_col   = "grey40",
                      border_lwd   = 0.2) {

  rug_position <- match.arg(rug_position)

  draw_fn <- function(ctx) {
    cols <- if (is.null(n_cols)) choose_grid(ctx$n_tree)$n_cols else n_cols
    rows <- if (is.null(n_rows)) ceiling(ctx$n_tree / cols)    else n_rows

    cell_h <- ctx$dy * cell_scale
    cell_w <- cell_h * (diff(ctx$last_pp$x.lim) / ctx$pin[1]) /
      (diff(ctx$last_pp$y.lim) / ctx$pin[2])

    variable <- !ctx$unanimous
    if (!any(variable)) return(invisible(NULL))

    plot_node_rug(
      presence     = ctx$presence[variable, , drop = FALSE],
      support      = if (is.null(ctx$support)) NULL else
        ctx$support[variable, , drop = FALSE],
      support_type = ctx$support_type,
      thresholds   = ctx$thresholds,
      cell_h       = cell_h,
      cell_w       = cell_w,
      n_cols       = cols,
      x_offset     = x_offset,
      y_offset     = y_offset,
      rug_position = rug_position,
      last_pp      = ctx$last_pp
    )
  }

  structure(list(type = "rug", draw = draw_fn), class = "phylorug_layer")
}


# ── dot_layer ────────────────────────────────────────────────────────────────

#' Add dots at unanimous nodes
#'
#' Draws a dot at each internal node where every comparison tree
#' recovers the same clade.
#'
#' @param col Dot colour. Default `"black"`.
#' @param cex Dot size. Default `0.45`. Auto-scales with tree size.
#' @param pch Point character. Default `16` (solid circle). Use `21` for
#'   a filled circle with separate border (set `fill` for inner colour).
#' @param fill Fill colour when using `pch = 21`. Default `"black"`.
#' @param auto_scale Logical. Scale dot size with tree size. Default `TRUE`.
#'
#' @returns A `phylorug_layer` object.
#' @export
dot_layer <- function(col        = "black",
                      cex        = 0.45,
                      pch        = 16,
                      fill       = "black",
                      auto_scale = TRUE) {

  draw_fn <- function(ctx) {
    if (!any(ctx$unanimous)) return(invisible(NULL))

    uni_rows <- which(ctx$unanimous)
    node_ids <- as.integer(rownames(ctx$presence)[uni_rows])
    valid    <- node_ids > ctx$ntip & node_ids <= length(ctx$last_pp$xx)

    if (any(valid)) {
      ids <- node_ids[valid]

      dot_cex <- if (auto_scale) {
        max(0.15, min(cex, 1.2 - 0.003 * ctx$ntip))
      } else {
        cex
      }

      graphics::points(
        ctx$last_pp$xx[ids],
        ctx$last_pp$yy[ids],
        pch = pch,
        cex = dot_cex,
        col = col,
        bg  = fill
      )
    }
  }

  structure(list(type = "dot", draw = draw_fn), class = "phylorug_layer")
}


# ── support_labels ───────────────────────────────────────────────────────────

#' Add backbone support labels
#'
#' Draws the backbone tree's own node support values (e.g., bootstrap)
#' beside each non-unanimous node.
#'
#' @param col Label colour. Default `"red"`.
#' @param cex Label size. `NULL` for auto-scaling with tree size.
#' @param font Font style: 1 = normal, 2 = bold, 3 = italic,
#'   4 = bold-italic. Default `1`.
#' @param adj Numeric vector of length 2. Horizontal and vertical
#'   adjustment for label position. Default `c(1.1, 1.4)`.
#' @param skip_unanimous Logical. Hide labels at unanimous nodes.
#'   Default `TRUE`.
#' @param round Integer. Number of decimal places. Default `2`.
#'
#' @returns A `phylorug_layer` object.
#' @export
support_labels <- function(col            = "red",
                           cex            = NULL,
                           font           = 1,
                           adj            = c(1.1, 1.4),
                           skip_unanimous = TRUE,
                           round          = 2) {

  draw_fn <- function(ctx) {
    labs <- ctx$backbone$node.label
    if (is.null(labs)) return(invisible(NULL))

    show_idx <- which(!is.na(labs) & nzchar(labs))
    if (length(show_idx) == 0) return(invisible(NULL))

    node_ids <- show_idx + ctx$ntip

    if (skip_unanimous) {
      uni_ids  <- as.integer(rownames(ctx$presence)[ctx$unanimous])
      keep     <- !(node_ids %in% uni_ids)
      show_idx <- show_idx[keep]
      node_ids <- node_ids[keep]
    }
    if (length(show_idx) == 0) return(invisible(NULL))

    num <- suppressWarnings(as.numeric(labs[show_idx]))
    txt <- ifelse(is.na(num), labs[show_idx],
                  format(round(num, round), trim = TRUE))

    label_cex <- if (is.null(cex)) {
      max(0.20, min(0.65, 0.85 - 0.002 * ctx$ntip))
    } else {
      cex
    }

    ape::nodelabels(
      text  = txt,
      node  = node_ids,
      frame = "none",
      cex   = label_cex,
      col   = col,
      font  = font,
      adj   = adj
    )
  }

  structure(list(type = "support_labels", draw = draw_fn),
            class = "phylorug_layer")
}


# ── legend_layer ─────────────────────────────────────────────────────────────

#' Add position and threshold legends
#'
#' Draws the position-numbering legend (which cell = which analysis) and,
#' in support mode, the threshold-shading legend.
#'
#' @param size Numeric multiplier for overall legend size. `1` = default,
#'   `1.5` = 50 percent larger, `0.7` = 30 percent smaller. Scales both
#'   text and cells proportionally.
#' @param text_cex Legend text size. `NULL` for auto-scaling (modified by
#'   `size`).
#' @param cell_size Legend cell size in inches. `NULL` for auto-scaling
#'   (modified by `size`).
#' @param font Font style: 1 = normal, 2 = bold, 3 = italic. Default `1`.
#' @param max_chars Maximum characters before truncating analysis names.
#'   Default `35`.
#'
#' @returns A `phylorug_layer` object.
#' @export
legend_layer <- function(size      = 1,
                         text_cex  = NULL,
                         cell_size = NULL,
                         font      = 1,
                         max_chars = 35L) {

  draw_fn <- function(ctx) {
    tree_names <- colnames(ctx$presence)
    n_tree     <- length(tree_names)
    n_cols_leg <- choose_grid(n_tree)$n_cols

    usr      <- graphics::par("usr")
    din      <- graphics::par("din")
    margin_x <- (usr[2] - usr[1]) * 0.015
    margin_y <- (usr[4] - usr[3]) * 0.015

    y0_top  <- usr[4] - margin_y
    x0_left <- usr[1] + margin_x

    # Legend sizing, scaled by size multiplier
    w_ref        <- din[1] * 0.85 + din[2] * 0.15
    leg_cell_in  <- if (is.null(cell_size)) w_ref * 0.018 * size else cell_size
    th_sq_in     <- if (is.null(cell_size)) w_ref * 0.008 * size else cell_size * 0.45
    pos_text_cex <- if (is.null(text_cex)) {
      min(0.55, max(0.20, w_ref * 0.045)) * size
    } else {
      text_cex
    }
    th_text_cex <- if (is.null(text_cex)) {
      min(0.50, max(0.18, w_ref * 0.041)) * size
    } else {
      text_cex
    }

    leg_cell_h <- graphics::yinch(leg_cell_in)
    leg_cell_w <- graphics::xinch(leg_cell_in)
    th_sq_h    <- graphics::yinch(th_sq_in)
    th_sq_w    <- graphics::xinch(th_sq_in)

    # Truncate long names
    display_names <- ifelse(
      nchar(tree_names) > max_chars,
      paste0(substr(tree_names, 1, max_chars - 1), "\u2026"),
      tree_names
    )

    old_family <- graphics::par("family")
    graphics::par(family = "sans")
    on.exit(graphics::par(family = old_family), add = TRUE)

    # Position legend — always topleft
    draw_position_legend(
      display_names, n_cols_leg,
      cell_w   = leg_cell_w,
      cell_h   = leg_cell_h,
      x0       = x0_left,
      y0       = y0_top,
      text_cex = pos_text_cex
    )

    # Threshold legend — always topright, support mode only
    if (ctx$mode == "support") {
      longest_label <- paste0(
        ">=80-97.9 (SH-aLRT) or >=95-97 (UFBoot2) ",
        "or >=0.95-0.98 (ASTRAL)"
      )
      th_width_in <- graphics::strwidth(
        longest_label, units = "inches", cex = th_text_cex
      ) + th_sq_in + 0.2

      x0_right <- usr[2] - margin_x - graphics::xinch(th_width_in)

      draw_threshold_legend(
        x0       = x0_right,
        y0       = y0_top,
        sq_h     = th_sq_h,
        sq_w     = th_sq_w,
        text_cex = th_text_cex
      )
    }
  }

  structure(list(type = "legend", draw = draw_fn), class = "phylorug_layer")
}
