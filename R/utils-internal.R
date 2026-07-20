# utils-internal.R
# Internal helper functions for phylorug. None are exported.


# choose_grid -----------------------------------------------------------------
choose_grid <- function(n_cells) {
  n_cols <- ceiling(sqrt(n_cells))
  n_rows <- ceiling(n_cells / n_cols)
  list(n_rows = n_rows, n_cols = n_cols)
}


# draw_position_legend --------------------------------------------------------
draw_position_legend <- function(analyses, n_cols, cell_w, cell_h,
                                 x0 = NULL, y0 = NULL) {
  n_an   <- length(analyses)
  n_rows <- ceiling(n_an / n_cols)

  usr <- graphics::par("usr")
  if (is.null(x0)) x0 <- usr[1] + 0.01 * (usr[2] - usr[1])
  if (is.null(y0)) y0 <- usr[4] - 0.01 * (usr[4] - usr[3])

  for (k in seq_len(n_an)) {
    row_idx <- ceiling(k / n_cols)
    col_idx <- ((k - 1) %% n_cols) + 1
    xleft   <- x0 + (col_idx - 1) * cell_w
    xright  <- xleft + cell_w
    ytop    <- y0 - (row_idx - 1) * cell_h
    ybottom <- ytop - cell_h
    graphics::rect(xleft, ybottom, xright, ytop,
                   col = "white", border = "black", lwd = 0.5)
    graphics::text((xleft + xright) / 2, (ytop + ybottom) / 2,
                   labels = k, cex = 0.5)
  }

  grid_right <- x0 + n_cols * cell_w
  key_x      <- grid_right + 0.5 * cell_w
  key_lines  <- paste0(seq_len(n_an), " \u2013 ", analyses)
  graphics::text(key_x, y0, labels = paste(key_lines, collapse = "\n"),
                 adj = c(0, 1), cex = 0.55, family = "sans")

  grid_bottom <- y0 - n_rows * cell_h
  invisible(grid_bottom)
}


# draw_threshold_legend -------------------------------------------------------
draw_threshold_legend <- function(x0, y0, cell_size = NULL) {
  x_inch <- graphics::xinch(1)
  y_inch <- graphics::yinch(1)

  sq   <- if (is.null(cell_size)) y_inch * 0.08 else cell_size
  sq_x <- x_inch * 0.08
  gap  <- sq * 0.4
  text_gap <- sq_x * 0.3

  rows <- list(
    list(fill = "#000000", pattern = "none",
         label = "\u226598 (SH-aLRT/UFBoot2) or \u22650.99 (ASTRAL)"),
    list(fill = "#5F5E5A", pattern = "none",
         label = "80\u201397.9 (SH-aLRT) or 95\u201397 (UFBoot2) or 0.95\u20130.98 (ASTRAL)"),
    list(fill = "#B4B2A9", pattern = "none",
         label = "50\u201379.9 (SH-aLRT) or 50\u201394 (UFBoot2) or 0.5\u20130.95 (ASTRAL)"),
    list(fill = "white", pattern = "dots",
         label = "<50 (SH-aLRT/UFBoot2) or <0.5 (ASTRAL)"),
    list(fill = "white", pattern = "none",
         label = "monophyly not supported"),
    list(fill = "white", pattern = "cross",
         label = "not computed")
  )

  for (i in seq_along(rows)) {
    r      <- rows[[i]]
    y_top  <- y0 - (i - 1) * (sq + gap)
    y_bot  <- y_top - sq
    xleft  <- x0
    xright <- x0 + sq_x

    graphics::rect(xleft, y_bot, xright, y_top,
                   col = r$fill, border = "black", lwd = 0.4)

    if (identical(r$pattern, "dots")) {
      n_per_side <- 4L
      w   <- xright - xleft
      h   <- y_top - y_bot
      off <- 1 / (2 * n_per_side)
      fracs <- seq(off, 1 - off, length.out = n_per_side)
      grid  <- expand.grid(fx = fracs, fy = fracs)
      graphics::points(xleft + grid$fx * w, y_bot + grid$fy * h,
                       pch = 16, cex = 0.2, col = "black")
    } else if (identical(r$pattern, "cross")) {
      graphics::segments(xleft, y_bot, xright, y_top,
                         col = "grey30", lwd = 0.5)
      graphics::segments(xleft, y_top, xright, y_bot,
                         col = "grey30", lwd = 0.5)
    }

    graphics::text(xright + text_gap, (y_top + y_bot) / 2,
                   labels = r$label, adj = c(0, 0.5),
                   cex = 0.45, family = "sans")
  }

  invisible(NULL)
}
# auto_canvas -----------------------------------------------------------------
# Compute PDF width and height from the number of tips.
# Targets ~0.18 inches per tip vertically, minimum 6 inches.
# Width is 60% of height, minimum 8 inches.
auto_canvas <- function(ntip) {
  height <- max(6, ntip * 0.18)
  width  <- max(8, height * 0.6)
  list(width = width, height = height)
}
