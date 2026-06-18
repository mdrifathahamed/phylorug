#' Draws the support-value grid(rugs) at every node of a plotted tree
#'
#' A phylogenetic tree must already be drawn on the active graphics device.
#' For each internal node, this function paints a small grid of coloured
#' rectangles, one per comparison tree, where the colour shows how strongly
#' that tree supports the node's clade. Together these grids are the rug, a
#' compact way to read agreement and conflict across analyses directly on the
#' tree. Cells are coloured from the normalised support value so that scales
#' from different methods are comparable; an absent clade is drawn white.
#' Optionally the raw support value is printed inside each cell for an exact
#' read when the figure is zoomed.
#'
#' Users do not call this directly; \code{plot_phylorug()} calls it after
#' drawing the tree and working out the cell geometry.
#'
#' @param tree The reference tree of class \code{"phylo"}, already plotted on
#'   the active device.
#' @param rug_mt The support matrix from \code{node_presence_matrix()}: first
#'   column \code{node_id}, remaining columns the comparison trees, holding raw
#'   support values.
#' @param cell_h,cell_w Numeric. Height and width of one cell, in the tree's
#'   plotting coordinates.
#' @param x_offset,y_offset Numeric. Shift the whole grid away from the node,
#'   as a fraction of the tree's width and height.
#' @param map_to_color A function taking a normalised value (0-1) and
#'   \code{pal_info}, returning a colour string.
#' @param pal_info Palette configuration passed to \code{map_to_color}.
#' @param n_cols Integer. Columns in each node's grid. Default \code{2}.
#' @param adaptive Logical. Reserved for adaptive cell sizing on dense trees.
#'   Default \code{TRUE}.
#' @param fill_fraction Numeric. Fraction of available space a grid fills,
#'   used by adaptive sizing. Default \code{0.4}.
#' @param show_values Logical. If \code{TRUE}, print the raw support value
#'   inside each cell. Default \code{FALSE}.
#'
#' @return Invisibly \code{NULL}; called for its drawing side effect.
#'
#' @keywords internal
plot_node_rug <- function(tree, rug_mt,
                          cell_h,
                          cell_w,
                          x_offset      = -0.0095,
                          y_offset      = 0.0023,
                          map_to_color,
                          pal_info,
                          n_cols        = 2,
                          adaptive      = TRUE,
                          fill_fraction = 0.4,
                          show_values   = FALSE) {
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylogenetic tree of class \"phylo\".",
         call. = FALSE
    )
  }
  if (!is.matrix(rug_mt) && !is.data.frame(rug_mt)) {
    stop("`rug_mt` must be a matrix or data frame.", call. = FALSE)
  }

  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  dx_offset <- max(last_pp$xx) * x_offset
  dy_offset <- max(last_pp$yy) * y_offset

  n_nodes <- nrow(rug_mt)
  n_cells <- ncol(rug_mt) - 1

  n_rows  <- ceiling(n_cells / n_cols)
  total_w <- n_cols * cell_w
  total_h <- n_rows * cell_h

  # Normalise each tree column once, for colour. Raw values are kept for text.
  raw_mat  <- rug_mt[, -1, drop = FALSE]
  norm_mat <- apply(raw_mat, 2, normalize_support)

  for (i in seq_len(n_nodes)) {
    node_id  <- rug_mt[i, 1]
    raw_vals <- as.numeric(raw_mat[i, ])
    nrm_vals <- as.numeric(norm_mat[i, ])

    x_center <- last_pp$xx[node_id] + dx_offset
    y_center <- last_pp$yy[node_id] + dy_offset

    x0 <- x_center - total_w / 2
    y0 <- y_center + total_h / 2

    for (k in seq_along(nrm_vals)) {
      nrm <- nrm_vals[k]
      raw <- raw_vals[k]

      row_idx <- ceiling(k / n_cols)
      col_idx <- ((k - 1) %% n_cols) + 1

      xleft   <- x0 + (col_idx - 1) * cell_w
      xright  <- xleft + cell_w
      ytop    <- y0 - (row_idx - 1) * cell_h
      ybottom <- ytop - cell_h

      # Absent clade -> white; otherwise colour from the normalised value
      col <- if (is.na(nrm)) "white" else map_to_color(nrm, pal_info)

      graphics::rect(xleft, ybottom, xright, ytop,
                     col    = col,
                     border = "black",
                     lwd    = 0.4
      )

      # Optional raw value printed in the cell (blank for absent clades)
      if (show_values && !is.na(raw)) {
        graphics::text(
          x      = (xleft + xright) / 2,
          y      = (ytop + ybottom) / 2,
          labels = format(raw, digits = 2),
          cex    = cell_h * 4,
          col    = "black"
        )
      }
    }
  }
  invisible(NULL)
}
