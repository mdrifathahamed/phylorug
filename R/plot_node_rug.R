#' Draw the support cells at every node of a plotted tree
#'
#' For each internal node, this function paints a small grid of rectangles,
#' one per analysis, where the shade shows how strongly that analysis supports
#' the node's clade. By default the cells are greyscale and a darker cell means
#' stronger support, so each analysis is identified by its position in the grid.
#' In colour mode each analysis keeps its own hue instead. An absent clade is
#' drawn white. Together these grids are the rug. When values are shown, the
#' support number is printed inside each cell in a contrasting colour.
#'
#' Users do not call this directly; \code{plot_phylorug()} calls it after
#' drawing the tree and working out the cell geometry.
#'
#' @param rug_mt The support matrix from \code{node_presence_matrix()}: first
#'   column \code{node_id}, remaining columns the analyses, holding raw values.
#' @param hues Character vector of base colours, one per analysis, in column
#'   order. Ignored when \code{bw = TRUE}.
#' @param cell_h,cell_w Numeric. Height and width of one cell, in the tree's
#'   plotting coordinates.
#' @param n_cols Integer. Columns in each node's grid.
#' @param x_offset,y_offset Numeric. Shift the whole grid away from the node,
#'   as a fraction of the tree's width and height.
#' @param rug_position One of \code{"outside"} (default) or \code{"inside"}.
#'   \code{"outside"} centres each rug on its node; \code{"inside"} shifts it
#'   left toward the root and up, into the angle where the branch splits.
#' @param show_values Logical. If \code{TRUE}, print the support value
#'   inside each cell. Ignored in presence mode. Default \code{FALSE}.
#' @param bw Logical. If \code{TRUE} (default), cells are greyscale and the
#'   support value sets how dark each cell is, so analysis identity comes from
#'   the cell's position in the grid rather than its colour. If \code{FALSE},
#'   each analysis is given its own hue. Default \code{TRUE}.
#' @param last_pp The stored \code{plot.phylo} coordinates, passed in by
#'   \code{plot_phylorug()}. If \code{NULL}, they are fetched from the active
#'   device.
#'
#' @return Invisibly \code{NULL}; called for its drawing side effect.
#'
#' @keywords internal
plot_node_rug <- function(rug_mt,
                          hues,
                          cell_h,
                          cell_w,
                          n_cols,
                          x_offset     = -0.0095,
                          y_offset     = 0.0023,
                          rug_position = c("outside", "inside"),
                          show_values  = FALSE,
                          bw           = TRUE,
                          last_pp      = NULL) {
  if (!is.matrix(rug_mt) && !is.data.frame(rug_mt)) {
    stop("`rug_mt` must be a matrix or data frame.", call. = FALSE)
  }
  rug_position <- match.arg(rug_position)

  if (is.null(last_pp)) {
    last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  }

  dx_offset <- max(last_pp$xx) * x_offset
  dy_offset <- max(last_pp$yy) * y_offset

  n_nodes <- nrow(rug_mt)
  n_an    <- ncol(rug_mt) - 1
  n_rows  <- ceiling(n_an / n_cols)
  total_w <- n_cols * cell_w
  total_h <- n_rows * cell_h

  # Raw values for the text; normalised (per column) for the colour
  raw_mat  <- rug_mt[, -1, drop = FALSE]
  norm_mat <- apply(raw_mat, 2, normalize_support)

  # Presence mode: every present cell is exactly 1, so no values to print
  non_na      <- raw_mat[!is.na(raw_mat)]
  is_presence <- length(non_na) > 0 && all(non_na == 1)
  if (is_presence) {
    show_values <- FALSE
  }

  for (i in seq_len(n_nodes)) {
    node_id  <- rug_mt[i, 1]
    raw_vals <- as.numeric(raw_mat[i, ])
    nrm_vals <- as.numeric(norm_mat[i, ])

    # Skip this node entirely if every analysis is absent (all-white rug)
    if (all(is.na(raw_vals))) next

    x_center <- last_pp$xx[node_id] + dx_offset
    y_center <- last_pp$yy[node_id] + dy_offset

    # outside: centred on the node (the default).
    # inside: tucked into the crook, just left of the branch with a small gap,
    # lifted up above the node.
    if (rug_position == "inside") {
      gap_x <- cell_w * 0.5
      gap_y <- cell_h * 0.23
      x0    <- x_center - total_w - gap_x
      y0    <- y_center + total_h + gap_y
    } else {
      x0 <- x_center
      y0 <- y_center + total_h / 2
    }

    for (k in seq_len(n_an)) {
      raw <- raw_vals[k]
      nrm <- nrm_vals[k]

      row_idx <- ceiling(k / n_cols)
      col_idx <- ((k - 1) %% n_cols) + 1

      xleft   <- x0 + (col_idx - 1) * cell_w
      xright  <- xleft + cell_w
      ytop    <- y0 - (row_idx - 1) * cell_h
      ybottom <- ytop - cell_h

      fill <- cell_color(nrm, hues[k], bw = bw)
      graphics::rect(xleft, ybottom, xright, ytop,
                     col    = fill,
                     border = "grey40",
                     lwd    = 0.2
      )
      if (show_values && !is.na(raw)) {
        graphics::text(
          x      = (xleft + xright) / 2,
          y      = (ytop + ybottom) / 2,
          labels = format(nrm, digits = 2),
          cex    = cell_h * 0.35,
          font   = 2,
          col    = contrast_text_color(fill)
        )
      }
    }
  }
  invisible(NULL)
}
