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

