# gg_rug_cells.R
# Builds the rug grid as plain data frames (rectangles + pattern overlays)
# instead of drawing directly with graphics::rect()/points()/segments().
# Reuses resolve_tier(), resolve_cell(), bin_support(), bin_fill() from
# plot_node_rug.R untouched -- those are pure functions with no base
# graphics calls, so they work identically here.

#' Build rug cell rectangles + pattern overlays for ggplot2
#'
#' @param presence Presence matrix (node ids as rownames).
#' @param support Support matrix, same shape, or `NULL`.
#' @param support_type Named vector mapping tree name to support measure.
#' @param thresholds Optional threshold overrides.
#' @param node_xy Data frame with columns `node`, `x`, `y`, from the
#'   `nodes` element of `gg_tree_layout()`.
#' @param cell_h,cell_w Cell height/width in data (tree) coordinates.
#' @param n_cols Integer. Grid columns.
#' @param rug_position `"outside"` or `"inside"`.
#' @param x_offset,y_offset Fractional shift of the whole grid.
#' @param max_x,max_y Tree extent, for offset scaling (matches
#'   `plot_node_rug()`'s `max(last_pp$xx) * x_offset` convention).
#'
#' @returns A list with `rects` (data.frame: xmin,xmax,ymin,ymax,fill),
#'   `cross` (data.frame of X-pattern segments for not-computed cells), and
#'   `dots` (data.frame of point positions for low-support hatch cells).
#' @noRd
gg_rug_cells <- function(presence,
                         support      = NULL,
                         support_type = NULL,
                         thresholds   = NULL,
                         node_xy,
                         cell_h,
                         cell_w,
                         n_cols,
                         rug_position = c("outside", "inside"),
                         x_offset     = 0,
                         y_offset     = 0,
                         max_x        = 1,
                         max_y        = 1) {

  rug_position <- match.arg(rug_position)
  tier <- resolve_tier(support, support_type)

  node_ids <- as.integer(rownames(presence))
  n_nodes  <- nrow(presence)
  n_tree   <- ncol(presence)
  n_rows   <- ceiling(n_tree / n_cols)
  total_w  <- n_cols * cell_w
  total_h  <- n_rows * cell_h
  tree_names <- colnames(presence)

  dx_offset <- max_x * x_offset
  dy_offset <- max_y * y_offset

  rects_list <- vector("list", n_nodes * n_tree)
  cross_list <- vector("list", n_nodes * n_tree)
  dots_list  <- vector("list", n_nodes * n_tree)
  ri <- 0L; ci <- 0L; di <- 0L

  for (i in seq_len(n_nodes)) {
    node_id <- node_ids[i]
    xy_row  <- node_xy[node_xy$node == node_id, ]
    if (nrow(xy_row) == 0) next  # node not in current layout, skip

    x_center <- xy_row$x + dx_offset
    y_center <- xy_row$y + dy_offset

    if (rug_position == "inside") {
      gap_x <- cell_w * 0.25
      gap_y <- cell_h * 0.23
      x0 <- x_center - total_w - gap_x
      y0 <- y_center + total_h + gap_y
    } else {
      x0 <- x_center + cell_w * 0.25
      y0 <- y_center + total_h / 2
    }

    p_vals <- as.numeric(presence[i, ])
    s_vals <- if (is.null(support)) rep(NA_real_, n_tree) else as.numeric(support[i, ])

    for (k in seq_len(n_tree)) {
      row_idx <- ceiling(k / n_cols)
      col_idx <- ((k - 1) %% n_cols) + 1

      xleft   <- x0 + (col_idx - 1) * cell_w
      xright  <- xleft + cell_w
      ytop    <- y0 - (row_idx - 1) * cell_h
      ybottom <- ytop - cell_h

      st <- if (is.null(support_type)) NA_character_ else support_type[[tree_names[k]]]

      cell <- resolve_cell(
        p = p_vals[k], s = s_vals[k],
        support_type = st, thresholds = thresholds, tier = tier
      )

      ri <- ri + 1L
      rects_list[[ri]] <- data.frame(
        xmin = xleft, xmax = xright, ymin = ybottom, ymax = ytop,
        fill = cell$fill, border = cell$border
      )

      if (identical(cell$pattern, "cross")) {
        ci <- ci + 1L
        cross_list[[ci]] <- data.frame(
          x = c(xleft, xleft), xend = c(xright, xright),
          y = c(ybottom, ytop), yend = c(ytop, ybottom)
        )
      } else if (identical(cell$pattern, "dots")) {
        n_per_side <- 4L
        off <- 1 / (2 * n_per_side)
        fracs <- seq(off, 1 - off, length.out = n_per_side)
        grid  <- expand.grid(fx = fracs, fy = fracs)
        di <- di + 1L
        dots_list[[di]] <- data.frame(
          x = xleft + grid$fx * (xright - xleft),
          y = ybottom + grid$fy * (ytop - ybottom)
        )
      }
    }
  }

  list(
    rects = do.call(rbind, rects_list[seq_len(ri)]),
    cross = if (ci > 0) do.call(rbind, cross_list[seq_len(ci)]) else NULL,
    dots  = if (di > 0) do.call(rbind, dots_list[seq_len(di)]) else NULL
  )
}
