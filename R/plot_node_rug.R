#' Draw the rug cells at every node of a plotted tree
#'
#' For each internal backbone node, paints a small grid of rectangles, one per
#' comparison tree, showing whether that tree recovered the node's clade and, in
#' support mode, how strongly. Together these grids are the rug.
#'
#' The appearance is set by the tier, which is resolved from the arguments
#' rather than named directly:
#' \itemize{
#'   \item Presence (Tier 1): `support` is `NULL`. Black = present,
#'     white = absent, grey = partial recovery in a pool.
#'   \item Not-computed (Tier 2): `support` is given but `support_type` is not.
#'     Adds a hatched cell for clades that could not be evaluated.
#'   \item Four-state (Tier 3): both `support` and `support_type` are given.
#'     Present cells are shaded by binned support; absent cells are
#'     marked distinctly; not-computed cells are hatched.
#' }
#'
#' Users do not call this directly; [plot_phylorug()] calls it after drawing the
#' tree and working out the cell geometry.
#'
#' @param npm Node presence or support matrix from [node_presence_matrix()].
#'   one row per internal backbone node (node numbers as rownames), one column
#'   per comparison tree. Cells are `1`, `0`, or a pool proportion.
#' @param support Support matrix of the same shape, or `NULL`. When `NULL`, the
#'  rug is drawn in presence mode.
#' @param support_type Named character vector mapping each comparison tree to
#'   its support measure `"ufboot"`, `"sh_alrt"`, `"lpp"`, `"posterior"`), or
#'   `NULL`. Required for binned shading.
#' @param thresholds Optional list overriding the built-in bin thresholds, keyed
#'   by support type.`NULL` uses the literature defaults.
#' @param cell_h,cell_w Numeric. Height and width of one cell, in the tree's
#'   plotting coordinates.
#' @param n_cols Integer. Columns in each node's grid.
#' @param x_offset,y_offset Numeric. Shift the whole grid away from the node, as
#'   a fraction of the tree's width and height.
#' @param rug_position One of `"outside"` (default) or `"inside"`.
#' @param last_pp The stored [ape::plot.phylo()] coordinates from
#'   [plot_phylorug()]. If `NULL`, fetched from the active device.
#'
#' @return Invisibly `NULL`; called for its drawing side effect.
#'
#' @keywords internal
plot_node_rug <- function(npm,
                          support      = NULL,
                          support_type = NULL,
                          thresholds   = NULL,
                          cell_h,
                          cell_w,
                          n_cols,
                          x_offset     = 0,
                          y_offset     = 0,
                          rug_position = c("outside", "inside"),
                          last_pp      = NULL) {

  if (!is.matrix(npm)) {
    stop("`npm` must be a matrix.", call. = FALSE)
  }
  if (!is.null(support) && !identical(dim(support), dim(npm))) {
    stop("`support` must have the same dimensions as `npm`.", call. = FALSE)
  }
  rug_position <- match.arg(rug_position)

  if (is.null(last_pp)) {
    last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  }

  # Tier is resolved from what was supplied.
  tier <- resolve_tier(support, support_type)

  # Node ids come from the rownames now, not a first column.
  node_ids <- as.integer(rownames(npm))

  dx_offset <- max(last_pp$xx) * x_offset
  dy_offset <- max(last_pp$yy) * y_offset

  n_nodes <- nrow(npm)
  n_tree  <- ncol(npm)
  n_rows  <- ceiling(n_tree / n_cols)
  total_w <- n_cols * cell_w
  total_h <- n_rows * cell_h

  tree_names <- colnames(npm)

  for (i in seq_len(n_nodes)) {
    node_id <- node_ids[i]
    p_vals  <- as.numeric(npm[i, ])
    s_vals  <- if (is.null(support)) rep(NA_real_, n_tree) else as.numeric(support[i, ])

    x_center <- last_pp$xx[node_id] + dx_offset
    y_center <- last_pp$yy[node_id] + dy_offset

    # outside: centred on the node.
    # inside: tucked into the crook, left of the branch and lifted up.
    if (rug_position == "inside") {
      gap_x <- cell_w * 0.25
      gap_y <- cell_h * 0.23
      x0    <- x_center - total_w - gap_x
      y0    <- y_center + total_h + gap_y
    } else {
      x0 <- x_center + cell_w * 0.25
      y0 <- y_center + total_h / 2
    }

    for (k in seq_len(n_tree)) {
      row_idx <- ceiling(k / n_cols)
      col_idx <- ((k - 1) %% n_cols) + 1

      xleft   <- x0 + (col_idx - 1) * cell_w
      xright  <- xleft + cell_w
      ytop    <- y0 - (row_idx - 1) * cell_h
      ybottom <- ytop - cell_h

      st <- if (is.null(support_type)) NA_character_ else support_type[[tree_names[k]]]

      cell <- resolve_cell(
        p            = p_vals[k],
        s            = s_vals[k],
        support_type = st,
        thresholds   = thresholds,
        tier         = tier
      )

      draw_cell(xleft, ybottom, xright, ytop, cell)
    }
  }

  invisible(NULL)
}


#' Resolve which tier to draw from the supplied arguments
#'
#' @noRd
resolve_tier <- function(support, support_type) {
  if (is.null(support)) {
    return(1L)                      # presence only
  }
  if (is.null(support_type)) {
    return(2L)                      # presence + not-computed, no binning
  }
  3L                                # full four-state with binned support
}


#' Decide one cell's fill and pattern
#'
#' Internal. Returns a list describing how to draw one cell: \code{fill} colour,
#' \code{hatch} logical, and \code{border}. This is the single place where the
#' four-state grammar lives; every tier routes through it.
#'
#' State priority is the same across tiers:
#'   not evaluable (presence NA)  ->  absent (presence 0)  ->  present.
#' What each state looks like depends on the tier.
#'
#' @noRd
resolve_cell <- function(p, s, support_type, thresholds, tier) {

  # Not evaluable: cross.
  if (is.na(p)) {
    return(list(fill = "white", pattern = "cross", border = "grey40"))
  }

  # Absent: white, no pattern.
  if (p == 0) {
    return(list(fill = "white", pattern = "none", border = "grey40"))
  }

  # Tier 1 and 2: greyscale by presence proportion.
  if (tier < 3L) {
    grey <- 1 - p
    fill <- grDevices::rgb(grey, grey, grey)
    return(list(fill = fill, pattern = "none", border = "grey40"))
  }

  # Tier 3, present but no support value.
  if (is.na(s)) {
    return(list(fill = "black", pattern = "none", border = "grey40"))
  }

  bin  <- bin_support(s, support_type, thresholds)
  spec <- bin_fill(bin)
  list(fill = spec$fill, pattern = spec$pattern, border = "grey40")
}

#' Draw one resolved cell
#'
#' @noRd
draw_cell <- function(xleft, ybottom, xright, ytop, cell) {
  graphics::rect(
    xleft, ybottom, xright, ytop,
    col    = cell$fill,
    border = cell$border,
    lwd    = 0.2
  )

  if (identical(cell$pattern, "cross")) {
    # Single X across the cell.
    graphics::segments(xleft, ybottom, xright, ytop,
                       col = "grey30", lwd = 0.5)
    graphics::segments(xleft, ytop, xright, ybottom,
                       col = "grey30", lwd = 0.5)
  }

  invisible(NULL)
}

#' Map a support value to an integer bin, by support type
#'
#' Internal. Returns an integer tier: 4 = very high, 3 = high, 2 = moderate,
#' 1 = low/rejected, NA passes through. Thresholds are per measure, following
#' each measure's own literature; no normalization across measures.
#'
#' STUB: the default thresholds below are the working values from the plan and
#' must be checked against the primary literature before release.
#'
#' @noRd
bin_support <- function(value, support_type, thresholds = NULL) {
  if (is.na(value)) {
    return(NA_integer_)
  }

  th <- if (!is.null(thresholds) && !is.null(thresholds[[support_type]])) {
    thresholds[[support_type]]
  } else {
    default_thresholds(support_type)
  }

  # th is c(very_high, high, moderate): the lower bound of each tier.
  if (value >= th[["very_high"]]) return(4L)
  if (value >= th[["high"]])      return(3L)
  if (value >= th[["moderate"]])  return(2L)
  1L
}


#' Built-in bin thresholds per support measure
#'
#' Internal. Lower bound of each tier. STUB values from the development plan,
#' pending literature validation (Minh et al. 2013 for UFBoot2, Guindon et al.
#' 2010 for SH-aLRT, Sayyari and Mirarab 2016 for ASTRAL LPP).
#'
#' @noRd
default_thresholds <- function(support_type) {
  switch(
    support_type %||% "ufboot",
    ufboot    = c(very_high = 95,   high = 95,   moderate = 50),
    sh_alrt   = c(very_high = 98,   high = 80,   moderate = 50),
    lpp       = c(very_high = 0.99, high = 0.95, moderate = 0.5),
    posterior = c(very_high = 0.99, high = 0.95, moderate = 0.5),
    # Unknown type: fall back to a 0-100 percentage scale.
    c(very_high = 95, high = 80, moderate = 50)
  )
}


#' Fill colour for an integer support bin
#'
#' Internal. Greyscale for the four support tiers, red for the lowest. The
#' absent (yellow) and not-computed (hatch) states are handled in
#' \code{resolve_cell()}, not here.
#'
#' Fill and pattern for an integer support bin
#'
#' @noRd
bin_fill <- function(bin) {
  if (is.na(bin)) {
    return(list(fill = "black", pattern = "none"))
  }
  switch(
    as.character(bin),
    "4" = list(fill = "#000000", pattern = "none"),   # very high
    "3" = list(fill = "#5F5E5A", pattern = "none"),   # high
    "2" = list(fill = "#B4B2A9", pattern = "none"),   # moderate
    "1" = list(fill = "#E8C547", pattern = "none"),   # low (<50)
    list(fill = "black", pattern = "none")
  )
}

#' Default value for NULL
#'
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
