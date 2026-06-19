#' Map each comparison tree to a row and column in the node grid
#'
#' Each backbone node carries a clade, and every analysis either recovered that
#' clade or did not, with some level of support. To show all analyses at once,
#' phylorug draws a small grid of cells at each node, one cell per analysis. The
#' grid is what makes the comparison readable. Instead of laying the analyses
#' out as a single long strip beside each node, which would run off the page
#' when there are several, they wrap into a compact block of rows and columns.
#'
#' This function builds that block. Given the support matrix and the chosen
#' grid shape, it assigns every analysis a fixed row and column, filling left to
#' right and top to bottom, and returns those positions for every node. It does
#' not draw anything and does not touch colour, it only decides where each cell
#' goes, so the drawing step has an explicit map to follow.
#'
#' \code{plot_phylorug()} calls this after choosing a grid shape that fits the
#' number of analyses. It is internal, but documented here because the cell
#' arrangement in the figure comes directly from these positions.
#'
#' @param rug_mt The support matrix from \code{node_presence_matrix()}. It is a
#'   matrix whose first column is \code{node_id} and whose remaining columns are
#'   the comparison trees. Each cell holds that tree's value for that node's
#'   clade, either a binary presence flag (\code{1} or \code{NA}) or a support
#'   value (\code{0}-\code{1} or \code{NA}), depending on how the matrix was
#'   built.
#'
#' @param n_rows,n_cols Optional. The grid shape at each node. Leave both
#'   \code{NULL} (default) to size the grid automatically to the number of
#'   trees, choosing a roughly square shape. Set one or both to fix the grid
#'   yourself; if you set only one, the other is derived to fit. There must be
#'   at least as many cells as trees.
#'
#' @return A long data frame, one row per node-and-tree combination. For a tree
#'   with 50 nodes and 4 analyses, that is 200 rows. Each row places one
#'   analysis's cell in the grid at one node, using five columns: \code{node_id}
#'   the backbone node the cell belongs to; \code{tree_name} which analysis the
#'   cell is for; \code{row} and \code{col} where the cell sits in that node's
#'   grid (\code{row} 1 is the top, \code{col} 1 is the left); and
#'   \code{cell_index} the analysis's order in the sequence (1, 2, 3 ...). The
#'   drawing function reads this table to know where to paint each cell.
#'
#' @keywords internal
rug_layout_map <- function(rug_mt,
                           n_rows = NULL,
                           n_cols = NULL) {
  if (!is.matrix(rug_mt) && !is.data.frame(rug_mt)) {
    stop(
      "`rug_mt` must be a matrix or data frame.",
      call. = FALSE
    )
  }
  if (ncol(rug_mt) < 2) {
    stop(
      "`rug_mt` must have at least two columns.",
      call. = FALSE
    )
  }
  tree_names <- colnames(rug_mt)[-1]
  n_cells    <- length(tree_names)

  # Size the grid to the number of trees when not given. choose_grid() picks a
  # roughly square shape; if only one dimension is fixed, derive the other.
  if (is.null(n_cols)) n_cols <- choose_grid(n_cells)$n_cols
  if (is.null(n_rows)) n_rows <- ceiling(n_cells / n_cols)

  if (n_cells > n_rows * n_cols) {
    stop(
      "There are ", n_cells, " trees but only ", n_rows * n_cols,
      " cells in a ", n_rows, " by ", n_cols,
      " grid. Increase `n_rows` or `n_cols`.",
      call. = FALSE
    )
  }
  layout_df <- do.call(rbind, lapply(seq_len(nrow(rug_mt)), function(i) {
    node_id <- rug_mt[i, 1]
    do.call(rbind, lapply(seq_along(tree_names), function(k) {
      data.frame(
        node_id          = node_id,
        tree_name        = tree_names[k],
        row              = ceiling(k / n_cols),
        col              = ((k - 1) %% n_cols) + 1,
        cell_index       = k,
        stringsAsFactors = FALSE
      )
    }))
  }))
  layout_df
}
