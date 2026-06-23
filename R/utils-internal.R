# detect Format automatically _used in read_trees()
detect_format <- function(path) {
  lines <- readLines(path, n = 5, warn = FALSE)

  # Empty file check
  if (length(lines) == 0 || all(!nzchar(trimws(lines)))) {
    stop("File appears to be empty: ", basename(path),
      call. = FALSE
    )
  }

  # NEXUS -- always starts with #NEXUS
  if (grepl("^\\s*#NEXUS", lines[1], ignore.case = TRUE)) {
    return("nexus")
  }

  if (any(grepl("\\[&&NHX:", lines, ignore.case = TRUE))) {
    return("newick")
  }

  # Default - Newick
  "newick"
}
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

# cell_color()
# Choose a cell's fill colour. The analysis's base hue sets the colour;
# the normalised support value (0-1) sets how deep it is, from near-white
# at 0 to the full hue at 1. An absent clade (NA) is white. Used by
# plot_node_rug() so each analysis keeps its own hue while intensity
# tracks support strength.
cell_color <- function(value, hue) {
  if (is.na(value)) {
    return("white")
  }
  value <- max(0, min(1, value))
  ramp <- grDevices::colorRampPalette(c("white", hue))
  ramp(101)[round(value * 100) + 1]
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

  # bar sits above the legend's top edge, leaving room for its "0/1" line and label
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
# get_palette()
# Return n distinct colours for n analyses. Uses a curated 12-colour set
# first (each chosen to be easy to tell apart), and falls back to R's
# Polychrome palette for larger n so colours stay distinguishable.
# Categorical, not a gradient: adjacent colours should look different.
get_palette <- function(n) {
  curated <- c(
    "#D85A30", "#185FA5", "#1D9E75", "#993556",
    "#854F0B", "#534AB7", "#3B6D11", "#A32D2D",
    "#2C8C9C", "#B5651D", "#6B4C9A", "#557A2F"
  )
  if (n <= length(curated)) {
    curated[seq_len(n)]
  } else {
    grDevices::palette.colors(n, palette = "Polychrome 36")
  }
}
