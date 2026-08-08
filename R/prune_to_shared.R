#' Trim the backbone and its comparison trees to their shared set of taxa
#'
#' First finds the taxa that the backbone and all comparison trees share, then
#' trims everything else from the backbone and from every comparison tree.
#'
#' @details
#' A clade cannot be compared across trees that were not run on the same
#' taxon. Where a comparison tree is missing a backbone taxon, every
#' backbone clade containing that taxon becomes unevaluable in that tree, and
#' [node_presence_matrix()] refuses to proceed. Trimming to the shared taxa
#' resolves this by making every tree consist the same taxon.
#'
#' The cost is that the questions become narrower. Any clade containing a
#' dropped taxon no longer exists on the backbone and cannot be reported, even
#' where most comparison trees recovered it. Check the report from
#' [check_taxa()] before trimming: if only one or two taxa are involved the loss
#' is minimal, but if several comparison trees each lack a different taxon, the
#' shared set can shrink quickly.
#'
#' @param backbone A phylogenetic tree of class `"phylo"`.
#'
#' @param trees A named list of `"phylo"` or `"multiPhylo"` objects,
#'    as returned by [read_trees()].
#'
#' @param verbose Logical. If `TRUE` (default), reports how many taxa were
#'   dropped and names them.
#'
#' @return A list with two elements: `backbone`, the pruned backbone tree, and
#'   `trees`, the pruned comparison trees, retaining their names and their
#'   original single-tree or pooled structure. The list carries a `"dropped"`
#'   attribute naming the taxa that were removed.
#'
#' @seealso [check_taxa()] to see which taxa are missing and from where, before
#'   deciding to prune.
#'
#' @export
#'
#' @examples
#' # Build a backbone and comparison trees that do NOT all share the same taxa,
#' # so there is something to trim. (In practice these come from read_trees().)
#' backbone <- ape::read.tree(
#'   text = "((((A,B),(C,D)),((E,F),(G,H))),(I,J));"
#' )
#' others <- list(
#'   # All ten taxa, same as the backbone:
#'   tree_1 = ape::read.tree(text = "((((A,B),(C,D)),((E,F),(G,H))),(I,J));"),
#'   # Missing J:
#'   tree_2 = ape::read.tree(text = "((((A,B),(C,D)),((E,F),(G,H))),I);"),
#'   # Missing I:
#'   tree_3 = ape::read.tree(text = "((((A,B),(C,D)),((E,F),(G,H))),J);")
#' )
#'
#' # See which taxa are missing, and from which trees, before trimming:
#' check_taxa(backbone, others)
#'
#' # Trim the backbone and every comparison tree to their shared taxa:
#' shared <- prune_to_shared(backbone, others)
#'
#' # Which taxa were dropped:
#' attr(shared, "dropped")
#'
#' # The pruned trees are ready for node_presence_matrix():
#' shared$backbone
#' shared$trees
prune_to_shared <- function(backbone, trees, verbose = TRUE) {

  if (!inherits(backbone, "phylo") && !inherits(backbone, "multiPhylo")) {
    stop(
      "`backbone` must be a `phylo` object, not ",
      paste(class(backbone), collapse = "/"), ".",
      call. = FALSE
    )
  }
  if (!inherits(trees, "list")) {
    stop(
      "`trees` must be a list of `phylo` and/or `multiPhylo` objects, one ",
      "element per analysis, as returned by `read_trees()`.",
      call. = FALSE
    )
  }
  if (length(trees) == 0L) {
    stop("`trees` is empty.", call. = FALSE)
  }
  bb_taxa <- tree_taxa(backbone)
  shared <- Reduce(
    intersect,
    lapply(trees, tree_taxa),
    accumulate = FALSE,
    init        = bb_taxa
  )
  if (length(shared) < 3L) {
    stop(
      "Only ", length(shared), " taxa are shared by the backbone and every ",
      "comparison tree. A tree needs at least three. The comparison trees are ",
      "too different to be compared after pruning.",
      call. = FALSE
    )
  }
  dropped <- setdiff(bb_taxa, shared)
  if (verbose) {
    if (length(dropped) == 0L) {
      message(
        "All ", length(bb_taxa), " taxa are shared. Nothing was pruned."
      )
    } else {
      message(
        "Dropped ", length(dropped), " of ", length(bb_taxa), " taxa, ",
        "leaving ", length(shared), ": ",
        paste(utils::head(dropped, 10L), collapse = ", "),
        if (length(dropped) > 10L) ", ..." else "",
        ". Any clade containing a dropped taxon can no longer be reported."
      )
    }
  }
  out <- list(
    backbone = prune_one(backbone, shared),
    trees    = stats::setNames(
      lapply(trees, prune_one, keep = shared),
      names(trees)
    )
  )
  attr(out, "dropped") <- dropped
  out
}

#' Sorted tip labels of one tree or pool
#'
#' Internal. A pool must share one taxon set across its trees; the first tree's
#' labels stand for the pool.
#'
#' @noRd
tree_taxa <- function(x) {
  sort(as_pool(x)[[1L]]$tip.label)
}

#' Prune one comparison tree to a set of taxa
#'
#' Internal. Handles a single tree or a pool, returning the same structure it
#' was given. Tips not in \code{keep} are dropped from every tree.
#'
#' @noRd
prune_one <- function(x, keep) {
  if (inherits(x, "phylo")) {
    drop <- setdiff(x$tip.label, keep)
    if (length(drop) == 0L) {
      return(x)
    }
    return(ape::drop.tip(x, drop))
  }
  pool <- as_pool(x)
  out  <- lapply(pool, function(tr) {
    drop <- setdiff(tr$tip.label, keep)
    if (length(drop) == 0L) tr else ape::drop.tip(tr, drop)
  })
  class(out) <- "multiPhylo"
  out
}
