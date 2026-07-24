# utils-internal.R

choose_grid <- function(n_cells) {
  n_cols <- ceiling(sqrt(n_cells))
  n_rows <- ceiling(n_cells / n_cols)
  list(n_rows = n_rows, n_cols = n_cols)
}

draw_position_legend <- function(analyses, n_cols, cell_w, cell_h,
                                 x0, y0, text_cex = 0.5) {
  n_an   <- length(analyses)
  n_rows <- ceiling(n_an / n_cols)

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
                   labels = k, cex = text_cex)
  }

  grid_right <- x0 + n_cols * cell_w
  key_x      <- grid_right + cell_w * 0.5
  key_lines  <- paste0(seq_len(n_an), " \u2013 ", analyses)
  graphics::text(key_x, y0,
                 labels = paste(key_lines, collapse = "\n"),
                 adj = c(0, 1), cex = text_cex, family = "sans")
  invisible(y0 - n_rows * cell_h)
}

draw_threshold_legend <- function(x0, y0, sq_h, sq_w, text_cex = 0.5) {
  gap      <- sq_h * 0.4
  text_gap <- sq_w * 0.3

  rows <- list(
    list(fill = "#000000", pattern = "none",
         label = ">=98 (SH-aLRT/UFBoot2) or >=0.99 (ASTRAL)"),
    list(fill = "#5F5E5A", pattern = "none",
         label = ">=80-97.9 (SH-aLRT) or >=95-97 (UFBoot2) or >=0.95-0.98 (ASTRAL)"),
    list(fill = "#B4B2A9", pattern = "none",
         label = ">=50-79.9 (SH-aLRT) or >=50-94 (UFBoot2) or >=0.5-0.95 (ASTRAL)"),
    list(fill = "#E8C547", pattern = "none",
         label = "<50 (SH-aLRT/UFBoot2) or <0.5 (ASTRAL)"),
    list(fill = "white",   pattern = "none",
         label = "monophyly not supported"),
    list(fill = "white",   pattern = "cross",
         label = "not computed")
  )

  for (i in seq_along(rows)) {
    r      <- rows[[i]]
    y_top  <- y0 - (i - 1) * (sq_h + gap)
    y_bot  <- y_top - sq_h
    xleft  <- x0
    xright <- x0 + sq_w
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
                       pch = 16, cex = text_cex * 0.4, col = "black")
    } else if (identical(r$pattern, "cross")) {
      graphics::segments(xleft, y_bot, xright, y_top, col = "grey30", lwd = 0.5)
      graphics::segments(xleft, y_top, xright, y_bot, col = "grey30", lwd = 0.5)
    }
    graphics::text(xright + text_gap, (y_top + y_bot) / 2,
                   labels = r$label, adj = c(0, 0.5),
                   cex = text_cex, family = "sans")
  }
  invisible(NULL)
}

# auto_canvas -----------------------------------------------------------------
#' Compute canvas dimensions from actual tree properties
#'
#' Height is set by per_tip. Width is computed from three real measurements:
#'
#' 1. TREE DEPTH — if branch lengths exist, the max root-to-tip distance
#'    determines how much horizontal room the phylogram needs. Deeper trees
#'    get more space (log-scaled, so it doesn't explode). If no branch
#'    lengths, falls back to a cladogram estimate from tip count.
#'
#' 2. LABEL WIDTH — longest tip label × character width at taxa_cex.
#'
#' 3. LEGEND BAND — scales with n_tree and mode.
#'
#' This prevents unnecessary stretching when branches are short, and ensures
#' deep phylograms get enough room to show branch-length variation.
#'
#' @param backbone The backbone tree.
#' @param ntip Integer.
#' @param n_tree Integer.
#' @param mode Character.
#' @param has_legend Logical.
#' @param per_tip Numeric.
#' @param taxa_cex Numeric.
#' @noRd
auto_canvas <- function(backbone,
                        ntip,
                        n_tree     = 4L,
                        mode       = "presence",
                        has_legend = TRUE,
                        per_tip    = 0.15,
                        taxa_cex   = 0.6) {

  # --- Height ---
  height <- max(6, ntip * per_tip)

  # --- Width component 1: Tree depth (inches) ---
  has_bl <- !is.null(backbone$edge.length) &&
    length(backbone$edge.length) > 0 &&
    any(backbone$edge.length > 0, na.rm = TRUE)

  if (has_bl) {
    # Phylogram: use actual root-to-tip distance.
    # node.depth.edgelength() returns distances from root for each node/tip.
    # Max value = deepest root-to-tip path = total tree span on x-axis.
    depths    <- ape::node.depth.edgelength(backbone)
    max_depth <- max(depths)

    # Convert tree depth to inches via log scale.
    # Log prevents tiny-depth trees from being too narrow and deep trees
    # from being absurdly wide. Range: 2–8 inches for tree portion.
    #   depth 0.01 → ~2.0″   depth 0.1 → ~4.5″
    #   depth 0.5  → ~6.8″   depth 2.0 → ~8.0″
    tree_width <- max(2, min(8, 2 * log2(1 + max_depth * 50)))
  } else {
    # Cladogram: no branch lengths. Width from tip count (log-scaled).
    # More tips = more branching levels = slightly more horizontal room.
    #   10 tips → ~3.3″   50 tips → ~5.6″   300 tips → ~6.0″
    tree_width <- max(2, min(6, 1.5 * log2(ntip)))
  }

  # --- Width component 2: Label width (inches) ---
  # Helvetica average character width ≈ 0.065 inches at cex=1
  max_nchar   <- max(nchar(backbone$tip.label))
  label_width <- max_nchar * 0.065 * taxa_cex

  # --- Width component 3: Rug space ---
  rug_width <- 0.5

  # --- Width component 4: Legend band ---
  leg_width <- 0
  if (has_legend) {
    leg_width <- max(2.5, n_tree * 0.25 + 1.5)
    if (mode == "support") leg_width <- leg_width + 3.5
  }

  # --- Total width ---
  width <- max(7, tree_width + label_width + rug_width + leg_width)

  list(width = width, height = height)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
