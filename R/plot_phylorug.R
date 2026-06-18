#' Draw a phylorug: a tree with node-support cells
#'
#' Draws a backbone tree and overlays, at each internal node, a small grid of
#' coloured cells summarising how several analyses support that node's clade.
#' Each analysis is given its own hue, so a cell's colour tells you which
#' analysis it is; the support value sets how deep that hue is, so a darker
#' cell means stronger support. An absent clade is drawn white. A legend maps
#' each hue to its analysis, and the support value can be printed inside each
#' cell for an exact read.
#'
#' This is the main user-facing plotting function. It draws the tree, works out
#' the cell geometry, chooses a grid that fits the analyses, and draws the rug.
#'
#' @param backbone The reference tree of class \code{"phylo"}.
#' @param rug_mt The support matrix from \code{node_presence_matrix()}: first
#'   column \code{node_id}, remaining columns the analyses.
#' @param hues Character vector of base colours, one per analysis, in the same
#'   order as the matrix columns. If \code{NULL} (default), a built-in palette
#'   is used.
#' @param n_rows,n_cols Optional grid shape. Leave \code{NULL} (default) to
#'   choose automatically a roughly square grid that fits the analyses. Set one
#'   or both to fix it; if only one is set, the other is chosen to fit.
#' @param show_values Logical. If \code{TRUE}, print each support value inside
#'   its cell. Has no effect on a presence matrix. Default \code{FALSE}.
#' @param legend Logical. If \code{TRUE} (default), draw a legend mapping hues
#'   to analyses.
#' @param cell_scale Numeric. Multiplier on the automatic cell height, for
#'   tuning density on crowded trees. Default \code{0.3}.
#' @param ... Passed to \code{ape::plot.phylo()}.
#'
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' backbone <- ape::rtree(6)
#' trees <- list(
#'   IQTREE = ape::rtree(6),
#'   ASTRAL = ape::rtree(6),
#'   MrBayes = ape::rtree(6)
#' )
#' mat <- node_presence_matrix(backbone, trees)
#' plot_phylorug(backbone, mat)
#' }
plot_phylorug <- function(backbone, rug_mt,
                          hues        = NULL,
                          n_rows      = NULL,
                          n_cols      = NULL,
                          show_values = FALSE,
                          legend      = TRUE,
                          cell_scale  = 0.3,
                          ...) {
  if (!inherits(backbone, "phylo")) {
    stop("`backbone` must be a phylogenetic tree of class \"phylo\".",
         call. = FALSE
    )
  }
  if (!is.matrix(rug_mt) && !is.data.frame(rug_mt)) {
    stop("`rug_mt` must be a matrix or data frame.", call. = FALSE)
  }
  if (ncol(rug_mt) < 2) {
    stop("`rug_mt` must have a node_id column and at least one analysis.",
         call. = FALSE
    )
  }

  analyses <- colnames(rug_mt)[-1]
  n_an     <- length(analyses)

  # Default hue palette, recycled if there are more analyses than colours
  if (is.null(hues)) {
    base_hues <- c("#D85A30", "#185FA5", "#1D9E75", "#993556",
                   "#854F0B", "#534AB7", "#3B6D11", "#A32D2D")
    hues <- base_hues[((seq_len(n_an) - 1) %% length(base_hues)) + 1]
  }
  if (length(hues) != n_an) {
    stop("`hues` must have one colour per analysis (", n_an, ").",
         call. = FALSE
    )
  }

  # Resolve the grid: automatic unless the user fixed it
  if (is.null(n_cols)) n_cols <- choose_grid(n_an)$n_cols
  if (is.null(n_rows)) n_rows <- ceiling(n_an / n_cols)
  if (n_rows * n_cols < n_an) {
    stop("A grid of ", n_rows, " by ", n_cols, " cannot hold ", n_an,
         " analyses. Increase `n_rows` or `n_cols`.", call. = FALSE
    )
  }

  # Draw the tree, then read back where each node landed
  ape::plot.phylo(backbone, ...)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  # Cell size from tip spacing, corrected for canvas aspect ratio
  ntip   <- ape::Ntip(backbone)
  yy_tip <- sort(last_pp$yy[seq_len(ntip)])
  dy     <- stats::median(diff(yy_tip))

  pin        <- graphics::par("pin")
  x_per_inch <- diff(last_pp$x.lim) / pin[1]
  y_per_inch <- diff(last_pp$y.lim) / pin[2]

  cell_h <- dy * cell_scale
  cell_w <- cell_h * (x_per_inch / y_per_inch)
  # Split nodes into unanimous (all analyses present) and variable
  cells     <- rug_mt[, -1, drop = FALSE]
  all_present <- apply(cells, 1, function(x) all(!is.na(x)))

  unanimous_mt <- rug_mt[all_present, , drop = FALSE]
  variable_mt  <- rug_mt[!all_present, , drop = FALSE]

  # Unanimous nodes: a single dot
  if (dot_unanimous && nrow(unanimous_mt) > 0) {
    node_ids <- unanimous_mt[, 1]
    graphics::points(
      last_pp$xx[node_ids], last_pp$yy[node_ids],
      pch = 16, cex = dot_cex, col = dot_col
    )
  }

  # Variable nodes: the full rug
  if (nrow(variable_mt) > 0) {
    plot_node_rug(
      rug_mt      = variable_mt,
      hues        = hues,
      cell_h      = cell_h,
      cell_w      = cell_w,
      n_cols      = n_cols,
      x_offset    = x_offset,
      y_offset    = y_offset,
      show_values = show_values,
      last_pp     = last_pp
    )
  }

  # Legend mapping hue to analysis
  if (legend) {
    graphics::legend(
      "topleft",
      legend = analyses,
      fill   = hues,
      bty    = "n",
      cex    = 0.8
    )
  }

  invisible(NULL)
}
