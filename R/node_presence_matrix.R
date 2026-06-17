#' Construct a node presence/support matrix for multiple phylogenetic trees
#'
#' Compares a list of phylogenetic trees against a reference tree topology set
#' as a backbone. For each internal node of the backbone, the function checks
#' whether the same clade (the same set of tip labels) appears in other trees in
#' comparison and records either binary presence value or the extracts the
#' support values.
#'
#' @param backbone A phylogenetic tree of class \code{"phylo"}. Each of its
#'   internal nodes becomes one row of the output matrix.
#'
#' @param trees A named list of \code{"phylo"} objects, or a \code{"multiPhylo"}
#'   object. Each tree becomes one column of the output matrix.
#'
#' @param use_support Logical. If \code{FALSE} (default), cells record binary
#'   presence: \code{1} if the backbone clade is found in the comparison tree,
#'   \code{NA} if it is not. If \code{TRUE}, cells record the extracted support
#'   value for matching clades, and \code{NA} for clades that are absent.
#'
#' @param support_col Integer, \code{1} or \code{2}. Selects which support value
#'   to read. Column 1 is the primary support metric for all standard methods
#'   (bootstrap, posterior probability, or ASTRAL local posterior) and is the
#'   right choice in almost all cases. Column 2 holds the optional second value
#'   that tools such as ASTRAL write (for example quartet support); methods that
#'   write only one value have \code{NA} there. Used only when
#'   \code{use_support = TRUE}. Defaults to \code{1}.
#'
#' @param round_support Integer or \code{NULL}. If supplied, support values are
#'   rounded to this many decimal places. Defaults to \code{NULL}.
#'
#' @return A numeric matrix. The first column, \code{node_id}, holds the
#'   backbone node numbers. Each following column is one comparison tree, named
#'   after the list element (or \code{tree_1}, \code{tree_2} and so on if the
#'   list has no names). A cell is \code{1} or a support value where the clade
#'   was found, and \code{NA} where it was not.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' shared_backbone <- ape::rtree(5, tip.label = c("A", "B", "C", "D", "E"))
#' alternative_trees <- list(
#'   UCE_partition  = ape::rtree(5, tip.label = c("A", "B", "C", "D", "E")),
#'   Mito_partition = ape::rtree(5, tip.label = c("A", "B", "C", "D", "E"))
#' )
#'
#' presence <- node_presence_matrix(
#'   backbone    = shared_backbone,
#'   trees       = alternative_trees,
#'   use_support = FALSE
#' )
#'
#' support <- node_presence_matrix(
#'   backbone    = shared_backbone,
#'   trees       = alternative_trees,
#'   use_support = TRUE,
#'   support_col = 1
#' )
#' }
node_presence_matrix <- function(backbone,
                                 trees,
                                 use_support   = FALSE,
                                 support_col   = 1,
                                 round_support = NULL) {
  if (!inherits(backbone, "phylo")) {
    stop(
      "`backbone` must be a phylogenetic tree of class \"phylo\".",
      call. = FALSE
    )
  }
  if (!inherits(trees, c("list", "multiPhylo"))) {
    stop(
      "`trees` must be a list of phylo objects or a multiPhylo object.",
      call. = FALSE
    )
  }
  if (!support_col %in% c(1, 2)) {
    stop(
      "`support_col` must be 1 or 2.",
      call. = FALSE
    )
  }
  # Get the sorted tip labels (the clade) below a given node
  get_clade <- function(tree, node) {
    tips <- phangorn::Descendants(tree, node, type = "tips")[[1]]
    sort(tree$tip.label[tips])
  }
  ntip     <- ape::Ntip(backbone)
  nnodes   <- backbone$Nnode
  bb_nodes <- (ntip + 1):(ntip + nnodes)
  # Clade at each backbone node
  bb_clades <- lapply(bb_nodes, get_clade, tree = backbone)
  # Build one column per comparison tree
  node_matrix <- vapply(seq_along(trees), function(i) {
    tr        <- trees[[i]]
    tr_nodes  <- (ape::Ntip(tr) + 1):(ape::Ntip(tr) + tr$Nnode)
    tr_clades <- lapply(tr_nodes, get_clade, tree = tr)
    if (use_support) {
      supp <- node_support(tr, round = round_support)
      if (nrow(supp) != length(tr_nodes)) {
        stop(
          "Number of node labels in tree ", i,
          " does not match the number of internal nodes.",
          call. = FALSE
        )
      }
      support_vec <- as.numeric(supp[[support_col]])
    }
    # For each backbone clade, record presence or support
    vapply(bb_clades, function(cl) {
      idx <- which(vapply(
        tr_clades,
        function(cl2) identical(cl2, cl),
        logical(1)
      ))
      if (length(idx) == 0) {
        NA_real_
      } else if (!use_support) {
        1
      } else {
        support_vec[idx[1]]
      }
    }, numeric(1))
  }, numeric(length(bb_clades)))
  # Ensure matrix shape (one tree would drop to a vector)
  node_matrix <- as.matrix(node_matrix)
  if (!is.null(names(trees))) {
    colnames(node_matrix) <- names(trees)
  } else {
    colnames(node_matrix) <- paste0("tree_", seq_along(trees))
  }
  # Prepend the backbone node ids
  cbind(node_id = bb_nodes, node_matrix)
}
