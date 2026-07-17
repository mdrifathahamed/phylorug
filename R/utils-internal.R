
# normalize_support()
# Rescale a vector of support values to the 0-1 range.
# Bootstrap-style values (0-100) are divided by 100; values already on a
# 0-1 scale are returned unchanged. Detection is per vector: if any value
# exceeds 1, the whole vector is treated as 0-100. All-NA vectors pass
# through untouched. Called at colour time, once per tree column.
normalize_support <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  if (any(x > 1, na.rm = TRUE)) {
    x / 100
  } else {
    x
  }
}

# choose_grid()
# Pick a roughly square grid that holds a given number of cells.
# Columns are the ceiling of the square root; rows are then whatever
# is needed to fit all cells. Keeps the rug compact rather than a long
# thin strip. Used by plot_phylorug() when the user does not set a grid.
choose_grid <- function(n_cells) {
  n_cols <- ceiling(sqrt(n_cells))
  n_rows <- ceiling(n_cells / n_cols)
  list(n_rows = n_rows, n_cols = n_cols)
}


# contrast_text_color()
# Pick black or white text so a value printed on a coloured cell stays
# readable. Uses perceptual luminance: light fills get black text, dark
# fills get white. Used by plot_node_rug() when show_values is TRUE.
contrast_text_color <- function(fill) {
  rgb_val <- grDevices::col2rgb(fill)
  lum <- (0.299 * rgb_val[1] + 0.587 * rgb_val[2] + 0.114 * rgb_val[3]) / 255
  if (lum > 0.6) "black" else "white"
}
# draw_support_gradient()
# Draw one small greyscale bar showing that cell darkness tracks support:
# light at 0, dark at 1. Drawn in the tree's own coordinate space, just below
# the swatch legend. Called by plot_phylorug() when gradient_legend = TRUE.
draw_support_gradient <- function(above, label = "Support value", n_steps = 60) {
  span_x <- graphics::par("usr")[2] - graphics::par("usr")[1]
  line_h <- graphics::strheight("M", cex = 0.8)

  x0    <- above$rect$left
  bar_w <- 0.045 * span_x
  bar_h <- line_h * 0.6

  # bar sits above the legend's top edge, leaving room for its "0/1" line and
  #label
  y_bar_top <- above$rect$top + line_h * 2.2

  # label above the bar
  graphics::text(x0, y_bar_top + line_h * 0.4, label, adj = 0, cex = 0.6)

  ramp  <- grDevices::colorRampPalette(c("grey90", "grey10"))(n_steps)
  seg_w <- bar_w / n_steps
  for (s in seq_len(n_steps)) {
    graphics::rect(x0 + (s - 1) * seg_w, y_bar_top - bar_h,
                   x0 + s * seg_w, y_bar_top, col = ramp[s], border = NA)
  }
  graphics::rect(x0, y_bar_top - bar_h, x0 + bar_w, y_bar_top,
                 border = "black", lwd = 0.3)

  # 0 and 1 just under the bar (sits between bar and legend)
  graphics::text(x0, y_bar_top - bar_h - line_h * 0.4, "0", adj = 0, cex = 0.55)
  graphics::text(x0 + bar_w, y_bar_top - bar_h - line_h * 0.4, "1", adj = 1, cex = 0.55)
  invisible(NULL)
}

# cell_color()
# Pick a cell's fill from its (normalised) support value.
# bw = TRUE : greyscale ramp, white (low) to black (high). Identity comes
#             from the cell's position, not its colour.
# bw = FALSE: white-to-hue ramp, so each analysis keeps its own colour.
# An absent value (NA) is always white, so an absent clade reads as blank.
cell_color <- function(value, hue, bw = FALSE) {
  if (is.na(value)) {
    return("white")
  }
  v <- max(0, min(1, value))
  if (bw) {
    grey_level <- 1 - v
    grDevices::rgb(grey_level, grey_level, grey_level)
  } else {
    rgb_hue <- grDevices::col2rgb(hue) / 255
    mixed <- (1 - v) * c(1, 1, 1) + v * rgb_hue[, 1]
    grDevices::rgb(mixed[1], mixed[2], mixed[3])
  }
}
# get_palette()
# Return n distinct colours for colour mode. Uses the colour-blind-safe
# Okabe-Ito palette first (up to 8), then falls back to the larger
# Polychrome set so colours stay distinguishable. Categorical, not a ramp.
get_palette <- function(n) {
  okabe_ito <- grDevices::palette.colors(palette = "Okabe-Ito")
  okabe_ito <- okabe_ito[okabe_ito != "#000000"]
  if (n <= length(okabe_ito)) {
    unname(okabe_ito[seq_len(n)])
  } else {
    unname(grDevices::palette.colors(n, palette = "Polychrome 36"))
  }
}
# draw_position_legend()
# The legend for black-and-white mode. Identity comes from a cell's position
# in the grid, not its colour, so this draws a small reference grid with each
# slot numbered 1..n_an, and beside it a key listing which analysis each
# number is. Laid out like the rugs (same n_cols, same fill order) so a reader
# can match a cell's position on the tree to a number here.
#
# cell_w and cell_h should be the rug's own cell width and height (from
# plot_phylorug), so the legend boxes are the same shape as the rug cells.
# Both are needed because on a tall canvas the x and y scales differ, so a
# single size would collapse the boxes into a thin line.
draw_position_legend <- function(analyses, n_cols, cell_w, cell_h,
                                 x0 = NULL, y0 = NULL) {
  n_an   <- length(analyses)
  n_rows <- ceiling(n_an / n_cols)

  usr <- graphics::par("usr")
  if (is.null(x0)) x0 <- usr[1] + 0.02 * (usr[2] - usr[1])
  if (is.null(y0)) y0 <- usr[4] - 0.02 * (usr[4] - usr[3])

  for (k in seq_len(n_an)) {
    row_idx <- ceiling(k / n_cols)
    col_idx <- ((k - 1) %% n_cols) + 1
    xleft   <- x0 + (col_idx - 1) * cell_w
    xright  <- xleft + cell_w
    ytop    <- y0 - (row_idx - 1) * cell_h
    ybottom <- ytop - cell_h
    graphics::rect(xleft, ybottom, xright, ytop,
                   col = "white", border = "black", lwd = 0.5
    )
    graphics::text((xleft + xright) / 2, (ytop + ybottom) / 2,
                   labels = k, cex = 0.5
    )
  }

  grid_right <- x0 + n_cols * cell_w
  key_x      <- grid_right + 0.5 * cell_w
  key_lines  <- paste0(seq_len(n_an), " - ", analyses)
  graphics::text(
    x = key_x, y = y0,
    labels = paste(key_lines, collapse = "\n"),
    adj = c(0, 1), cex = 0.6
  )

  invisible(NULL)
}
