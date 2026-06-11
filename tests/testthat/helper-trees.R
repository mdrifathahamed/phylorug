make_tree_list <- function(n_trees = 3, n_tips = 5) {
  tips  <- paste0("Sp_", seq_len(n_tips))
  trees <- replicate(
    n_trees,
    ape::rtree(n_tips, tip.label = tips),
    simplify = FALSE
  )
  names(trees) <- paste0("tree_", seq_len(n_trees))
  trees
}

make_compressed_multiphylo <- function(n_trees = 3, n_tips = 5) {
  tips  <- paste0("Sp_", seq_len(n_tips))
  trees <- replicate(
    n_trees,
    ape::rtree(n_tips, tip.label = tips),
    simplify = FALSE
  )
  class(trees) <- "multiPhylo"
  attr(trees, "TipLabel") <- tips
  trees
}
