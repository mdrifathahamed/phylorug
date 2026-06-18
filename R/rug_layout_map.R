#' Makes a layout for the comparison trees as a grid of cells at each node
#'
#' When several analyses are compared on one tree, each backbone node ends up
#' carrying a small block of cells, one for every comparison tree. Before
#' those cells can be drawn, this function decides their arrangement. Given the
#' support matrix and a grid shape, it returns the row and column each tree
#' occupies at every node, filling the grid from left to right and top to
#' bottom. It computes positions only.
#'
#' Users do not normally call this themselves. \code{plot_phylorug()} calls it,
#' having first chosen a grid shape that fits the number of trees. It is
#' documented here so the arrangement of cells in the final figure can be traced
#' and, if needed, adjusted.
#'
#' @param rug_mt The support matrix from \code{node_presence_matrix()}: a matrix
#'   or data frame whose first column is \code{node_id} and whose remaining
#'   columns are the comparison trees.
#'
#' @param n_rows Integer. The number of rows in the grid at each node. Together
#'   with \code{n_cols} this sets how many cells are available; there must be at
#'   least as many cells as trees. Default \code{2}.
#'
#' @param n_cols Integer. The number of columns in the grid at each node.
#'   Default \code{3}.
#'
#' @return A long-format data frame with one row for every node-and-tree pair
#'   and five columns: \code{node_id} (the backbone node), \code{tree_name} (the
#'   comparison tree), \code{row} (its grid row, 1 = top), \code{col} (its grid
#'   column, 1 = left), and \code{cell_index} (its position in sequence). The
#'   plotting functions use these positions to place each cell.
#'
#' @keywords internal
rug_layout_map <- function(rug_mt,
                           n_rows = 2,
                           n_cols = 3) {
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
