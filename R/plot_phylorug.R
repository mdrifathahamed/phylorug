#' Draw a phylorug: a backbone tree with node rugs
#'
#' Draws a backbone tree and overlays, at each internal node, a small grid of
#' cells (a rug) summarising how several comparison trees treat that node's
#' clade. The appearance is set by \code{mode}:
#'
#' \itemize{
#'   \item \code{"presence"} (default): each cell is black (clade recovered),
#'     white (rejected), or grey (recovered in part of a pool). This is the
#'     simplest, most robust view and needs no support values.
#'   \item \code{"support"}: present cells are shaded by binned support, using
#'     thresholds appropriate to each comparison tree's support measure. Absent
#'     cells are marked distinctly and clades that could not be evaluated are
#'     hatched. Requires \code{support_type}.
#' }
#'
#' Nodes where every comparison tree recovered the clade can be marked with a
#' single dot instead of a full grid, keeping the figure clean and drawing the
#' eye to conflict.
#'
#' @param backbone The reference tree of class \code{"phylo"}.
#'
#' @param npm The list returned by \code{\link{node_presence_matrix}}, with a
#'   \code{presence} element and one or more \code{support_1}, \code{support_2},
#'   \code{support_3} elements.
#'
#' @param mode One of \code{"presence"} (default) or \code{"support"}.
#'
#' @param support_col Integer 1-3. Which support matrix to shade by, when
#'   \code{mode = "support"}. Selects \code{npm$support_1},
#'   \code{npm$support_2}, or \code{npm$support_3}. Default \code{1}.
#'
#' @param support_type Named character vector mapping each comparison tree to
#'   its support measure (\code{"ufboot"}, \code{"sh_alrt"}, \code{"lpp"},
#'   \code{"posterior"}). Required when \code{mode = "support"}.
#'
#' @param thresholds Optional list overriding the built-in bin thresholds, keyed
#'   by support type. \code{NULL} uses the literature defaults.
#'
#' @param n_rows,n_cols Optional grid shape at each node. Leave \code{NULL}
#'   (default) to choose a roughly square grid automatically.
#'
#' @param legend Logical. If \code{TRUE} (default), draw the legend.
#'
#' @param legend_position Where to place the legend: \code{"topleft"}
#'   (default), \code{"bottomleft"}, \code{"topright"}, or
#'   \code{"bottomright"}.
#'
#' @param show_support Logical. If \code{TRUE} (default), display the
#'   backbone tree's own node support labels (in red) on the tree.
#'
#' @param cell_scale Numeric. Multiplier on the automatic cell height. Default
#'   \code{0.3}.
#'
#' @param x_offset,y_offset Numeric. Shift the whole grid away from the node, as
#'   a fraction of the tree's width and height. Default \code{0}.
#'
#' @param rug_position One of \code{"outside"} (default) or \code{"inside"}.
#'
#' @param dot_unanimous Logical. If \code{TRUE} (default), nodes where every
#'   comparison tree recovered the clade are marked with a single dot instead of
#'   a full grid.
#'
#' @param dot_col Colour of the unanimous-node dot. Default \code{"black"}.
#'
#' @param dot_cex Size of the unanimous-node dot. Default \code{0.6}.
#'
#' @param ... Passed to \code{ape::plot.phylo()}.
#'
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' backbone <- trees[["iqtree"]]
#' others   <- trees[names(trees) != "iqtree"]
#' npm      <- node_presence_matrix(backbone, others, support_col = c(1, 2))
#'
#' # Presence mode (default)
#' plot_phylorug(backbone, npm)
#'
#' # Support mode, shading by the first support column
#' plot_phylorug(
#'   backbone, npm,
#'   mode         = "support",
#'   support_col  = 1,
#'   support_type = c(astral = "lpp", raxml = "ufboot")
#' )
#' }
plot_phylorug <- function(backbone, npm,
                          mode          = c("presence", "support"),
                          support_col   = 1,
                          support_type  = NULL,
                          thresholds    = NULL,
                          n_rows        = NULL,
                          n_cols        = NULL,
                          legend        = TRUE,
                          legend_position = c("topleft", "bottomleft",
                                              "topright", "bottomright"),
                          show_support  = TRUE,
                          cell_scale    = 0.3,
                          x_offset      = 0,
                          y_offset      = 0,
                          rug_position  = c("outside", "inside"),
                          dot_unanimous = TRUE,
                          dot_col       = "black",
                          dot_cex       = 0.6,
                          ...) {

  # --- Input validation ------------------------------------------------------
  if (!inherits(backbone, "phylo")) {
    stop("`backbone` must be a phylogenetic tree of class \"phylo\".",
         call. = FALSE)
  }
  if (!is.list(npm) || !("presence" %in% names(npm))) {
    stop("`npm` must be the list returned by `node_presence_matrix()`, with a ",
         "`presence` element.", call. = FALSE)
  }

  mode            <- match.arg(mode)
  rug_position    <- match.arg(rug_position)
  legend_position <- match.arg(legend_position)

  presence <- npm$presence

  if (ncol(presence) < 1L) {
    stop("`npm` has no comparison trees to plot.", call. = FALSE)
  }

  # --- Support matrix selection -----------------------------------------------
  support <- NULL
  if (mode == "support") {
    if (is.null(support_type)) {
      stop("`mode = \"support\"` requires `support_type`, naming each ",
           "comparison tree's support measure (e.g. c(astral = \"lpp\")).",
           call. = FALSE)
    }
    support_name <- paste0("support_", support_col)
    if (!support_name %in% names(npm)) {
      stop("`npm` has no `", support_name, "`. Rebuild the matrix with a ",
           "matching `support_col`, or choose a lower `support_col`.",
           call. = FALSE)
    }
    support <- npm[[support_name]]
  }

  tree_names <- colnames(presence)
  n_tree     <- length(tree_names)

  # --- Grid shape -------------------------------------------------------------
  if (is.null(n_cols)) n_cols <- choose_grid(n_tree)$n_cols
  if (is.null(n_rows)) n_rows <- ceiling(n_tree / n_cols)
  if (n_rows * n_cols < n_tree) {
    stop("A grid of ", n_rows, " by ", n_cols, " cannot hold ", n_tree,
         " comparison trees. Increase `n_rows` or `n_cols`.", call. = FALSE)
  }

  # --- Auto y.lim for legend space -------------------------------------------
  dots <- list(...)
  ntip <- ape::Ntip(backbone)
  if (is.null(dots$y.lim) && legend) {
    # Scale legend space with cell_scale so larger rugs get more room.
    # Base: position grid (~3 rows) + threshold key (~6 rows) + padding.
    base_rows   <- if (mode == "support") 9 else 3
    legend_rows <- ceiling(base_rows * (1 + cell_scale))
    if (grepl("top", legend_position)) {
      dots$y.lim <- c(0, ntip + legend_rows)
    } else {
      dots$y.lim <- c(-legend_rows, ntip + 1)
    }
  }

  do.call(ape::plot.phylo, c(list(x = backbone), dots))
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  # Compute early — needed by backbone support label filtering below.
  unanimous <- apply(presence, 1, function(p) all(p == 1))

  # --- Backbone support labels ------------------------------------------------
  if (show_support && !is.null(backbone$node.label)) {
    labs <- backbone$node.label
    show <- which(!is.na(labs) & nzchar(labs))
    if (length(show) > 0) {
      node_ids <- show + ntip
      # Hide labels at unanimous nodes (already marked with a dot)
      if (dot_unanimous) {
        unanimous_ids <- as.integer(rownames(presence)[unanimous])
        keep <- !(node_ids %in% unanimous_ids)
        show     <- show[keep]
        node_ids <- node_ids[keep]
      }
      if (length(show) > 0) {
        ape::nodelabels(
          text  = labs[show],
          node  = node_ids,
          frame = "none",
          cex   = 0.35,
          col   = "red",
          adj   = c(1.1, 1.4)
        )
      }
    }
  }

  # --- Cell geometry ----------------------------------------------------------
  yy_tip     <- sort(last_pp$yy[seq_len(ntip)])
  dy         <- stats::median(diff(yy_tip))
  pin        <- graphics::par("pin")
  x_per_inch <- diff(last_pp$x.lim) / pin[1]
  y_per_inch <- diff(last_pp$y.lim) / pin[2]
  cell_h     <- dy * cell_scale
  cell_w     <- cell_h * (x_per_inch / y_per_inch)


  if (dot_unanimous && any(unanimous)) {
    ids <- as.integer(rownames(presence)[unanimous])
    graphics::points(
      last_pp$xx[ids], last_pp$yy[ids],
      pch = 16, cex = dot_cex, col = dot_col
    )
  }

  variable <- if (dot_unanimous) !unanimous else rep(TRUE, nrow(presence))
  if (any(variable)) {
    plot_node_rug(
      presence     = presence[variable, , drop = FALSE],
      support      = if (is.null(support)) NULL else support[variable, , drop = FALSE],
      support_type = support_type,
      thresholds   = thresholds,
      cell_h       = cell_h,
      cell_w       = cell_w,
      n_cols       = n_cols,
      x_offset     = x_offset,
      y_offset     = y_offset,
      rug_position = rug_position,
      last_pp      = last_pp
    )
  }

  # --- Legend -----------------------------------------------------------------
  if (legend) {
    usr    <- graphics::par("usr")
    y_inch <- graphics::yinch(1)
    margin <- 0.02

    # Horizontal anchor
    if (grepl("left", legend_position)) {
      leg_x0 <- usr[1] + margin * (usr[2] - usr[1])
    } else {
      leg_x0 <- usr[2] - margin * (usr[2] - usr[1]) - graphics::xinch(2.5)
    }

    # Vertical anchor
    if (grepl("top", legend_position)) {
      leg_y0 <- usr[4] - margin * (usr[4] - usr[3])
    } else {
      leg_y0 <- usr[3] + margin * (usr[4] - usr[3]) + y_inch * 1.5
    }

    legend_bottom <- draw_position_legend(
      tree_names, n_cols,
      cell_w = cell_w * 3, cell_h = cell_h * 3,
      x0 = leg_x0, y0 = leg_y0
    )

    if (mode == "support") {
      draw_threshold_legend(
        x0 = leg_x0,
        y0 = legend_bottom - y_inch * 0.14
      )
    }
  }

  invisible(NULL)
}
