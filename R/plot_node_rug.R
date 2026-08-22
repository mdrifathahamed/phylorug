#' Draw the rug at every internal node of a plotted tree
#'
#' @description
#' For each internal node of the backbone tree, paints a small grid of
#' rectangles, one per comparison tree, showing whether that tree recovered the
#' node's clade and, in support mode, how strongly. Together these grids are the
#' rug.
#'
#' The appearance is set by a tier, resolved from which arguments are supplied
#' rather than named directly. ([plot_phylorug()] selects the tier for the user
#' through its `mode` argument.)
#' \itemize{
#'   \item Tier 1, presence: `support` is `NULL`. Cells are black for present
#'     and white for absent. Grey cells appear only when `pool_threshold = 0`
#'     was used in [node_presence_matrix()], indicating partial recovery
#'     across a pool of equally optimal trees.
#'   \item Tier 2, support: both `support` and `support_type` are supplied.
#'     Recovered cells are shaded by binned support strength, from black (very
#'     high) through greys to yellow (low). A cell is white when the tree does
#'     not recover the clade at all, and red when the tree recovers the clade
#'     but carries no support value for it (an unscored node).
#' }
#'
#' Users do not call this directly; [plot_phylorug()] calls it after drawing the
#' tree and working out the cell geometry.
#'
#' @param npm The presence matrix from [node_presence_matrix()]. One row per
#'   internal backbone node (node numbers as rownames), one column per
#'   comparison tree. Cells are `1` (recovered), `0` (not recovered), or
#'   a proportion between 0 and 1 when `pool_threshold = 0` was used.
#'
#' @param support Support matrix of the same shape as `npm`, or `NULL`. When
#'   `NULL`, the rug is drawn in presence mode.
#'
#' @param support_type Named character vector mapping each comparison tree to
#'   its support measure (e.g. `"ufboot"`, `"sh_alrt"`, `"lpp"`,
#'   `"jackknife"`, `"bremer_ratio"`), or `NULL`. When `NULL`, values >1
#'   are auto-normalized to 0-1 and binned against universal thresholds.
#'
#' @param thresholds Optional list overriding the built-in bin thresholds, keyed
#'   by support type.`NULL` uses the literature defaults.
#'
#' @param cell_h,cell_w Numeric. Height and width of one cell, in the tree's
#'   plotting coordinates.
#'
#' @param n_cols Integer. Columns in each node's grid.
#'
#' @param x_offset,y_offset Numeric. Shift the whole grid away from the node, as
#'   a fraction of the tree's width and height.
#'
#' @param rug_position One of `"outside"` or `"inside"`(default).
#'
#' @param last_pp Plot coordinates from [ape::plot.phylo()] giving the x/y
#'   position of every node and tip, used to place each rug grid at its node.
#'   [plot_phylorug()] always supplies this. If called directly with `NULL`,
#'   the coordinates are fetched from the most recently drawn tree.
#'
#' @return Returns nothing; it draws the rug cells directly onto the tree.
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
                          rug_position = c("inside", "outside"),
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
    s_vals  <- if (is.null(support)) {
      rep(NA_real_, n_tree)
    } else {
      as.numeric(support[i, ])
    }

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

      st <- if (is.null(support_type)) {
        NA_character_
      } else {
        support_type[[tree_names[k]]]
      }

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
    return(1L)   # presence
  }
  2L             # support, fully specified (support_type guaranteed by caller)
}


#' Decide one cell's fill and pattern
#'
#' Internal. Returns a list describing how to draw one cell: `fill` colour,
#' `pattern`, and `border`. This is the single place where the cell colour
#' grammar lives; both tiers route through it.
#'
#' State priority is the same across tiers:
#'  not computed (presence NA)  ->  not recovered (presence 0)  ->  present.
#'
#' @noRd
resolve_cell <- function(p, s, support_type, thresholds, tier) {
  # Not computed: red.
  if (is.na(p)) {
    return(list(fill = "#D64545", pattern = "none", border = "grey40"))
  }
  # Not recovered: white.
  if (p == 0) {
    return(list(fill = "white", pattern = "none", border = "grey40"))
  }
  # Tier 1 (presence): greyscale by recovery proportion.
  if (tier == 1L) {
    grey <- 1 - p
    return(list(fill = grDevices::rgb(grey, grey, grey),
                pattern = "none", border = "grey40"))
  }
  # Tier 2 (support): recovered but no support value -> treat as not-computed.
  if (is.na(s)) {
    return(list(fill = "#D64545", pattern = "none", border = "grey40"))
  }
  # Tier 2 (support): bin the value.
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
  invisible(NULL)
}

#' Map a support value to an integer bin
#'
#' Internal. Returns an integer bin: 4 = very high, 3 = high, 2 = moderate,
#' 1 = low. NA passes through. When `support_type` is a recognized string,
#' metric-specific thresholds are applied. When `support_type` is NULL,
#' values >1 are divided by 100 and binned against universal 0-1 thresholds.
#' Custom thresholds override both, keyed by metric name or `"universal"`.
#'
#' @noRd
bin_support <- function(value, support_type = NULL, thresholds = NULL) {
  if (is.na(value)) {
    return(NA_integer_)
  }

  th <- if (!is.null(thresholds) && !is.null(support_type) &&
            !is.null(thresholds[[support_type]])) {
    # User custom thresholds for a named metric
    thresholds[[support_type]]
  } else if (!is.null(thresholds) &&
             (is.null(support_type) || is.na(support_type)) &&
             !is.null(thresholds[["universal"]])) {
    # User custom universal thresholds
    thresholds[["universal"]]
  } else {
    # Named metric defaults or universal defaults
    default_thresholds(support_type)
  }

  # Auto-normalize when thresholds are on 0-1 scale but value is > 1
  if (value > 1 && all(th <= 1)) {
    value <- value / 100
  }

  if (value >= th[["very_high"]]) return(4L)
  if (value >= th[["high"]])      return(3L)
  if (value >= th[["moderate"]])  return(2L)
  1L
}
#' Built-in bin thresholds per support measure
#'
#' Internal. Lower bound of each tier. When `support_type` is NULL,
#' returns universal 0-1 thresholds for use with auto-normalized values.
#' When a recognized string is passed, returns metric-specific thresholds
#' on each metric's native scale. Unrecognized strings fall through to
#' universal 0-1 thresholds.
#'
#' @noRd
default_thresholds <- function(support_type) {
  if (is.null(support_type) || is.na(support_type)) {
    return(c(very_high = 0.95, high = 0.80, moderate = 0.50))
  }
  switch(
    support_type,
    ufboot         = c(very_high = 95,   high = 80,   moderate = 50),
    sh_alrt        = c(very_high = 80,   high = 70,   moderate = 50),
    lpp            = c(very_high = 0.95, high = 0.90, moderate = 0.50),
    posterior      = c(very_high = 0.99, high = 0.95, moderate = 0.50),
    jackknife      = c(very_high = 95,   high = 80,   moderate = 50),
    bootstrap      = c(very_high = 95,   high = 80,   moderate = 50),
    bremer_ratio   = c(very_high = 0.95, high = 0.80, moderate = 0.50),
    transfer_boot  = c(very_high = 95,   high = 80,   moderate = 50),
    # Unrecognized string: universal 0-1 thresholds
    c(very_high = 0.95, high = 0.80, moderate = 0.50)
  )
}

#' Fill and pattern for an integer support bin
#'
#' Internal. Greyscale for the top three support tiers, yellow for the lowest.
#' The not-recovered (white) and not-computed (red) states are handled in
#' `resolve_cell()`, not here.
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
