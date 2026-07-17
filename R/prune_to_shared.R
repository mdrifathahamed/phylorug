#' Reduce a backbone and its comparison trees to their shared taxa
#'
#' Drops from the backbone, and from every comparison tree, any taxon that is
#' not present in all of them. The result is a set of trees that can be compared
#' clade by clade by \code{\link{node_presence_matrix}}.
#'
#' @details
#' A clade cannot be compared across trees that were not run on the same
#' terminals. Where a comparison tree is missing a backbone taxon, every
#' backbone clade containing that taxon becomes unevaluable in that tree, and
#' \code{\link{node_presence_matrix}} refuses to proceed. Pruning to the shared
#' taxa resolves this by making every tree ask the same question.
#'
#' The cost is that the question narrows. Any clade containing a dropped taxon
#' no longer exists on the backbone and cannot be reported, even where most
#' comparison trees recovered it. Inspect the report from
#' \code{\link{check_taxa}} before pruning: if only one or two taxa are involved
#' the loss is slight, but if several comparison trees each lack a different
#' taxon the shared set can shrink quickly.
#'
#' @param backbone A phylogenetic tree of class \code{"phylo"}.
#'
#' @param trees A named list of \code{"phylo"} or \code{"multiPhylo"} objects,
#'   one per comparison, as returned by \code{\link{read_trees}}.
#'
#' @param verbose Logical. If \code{TRUE} (default), reports how many taxa were
#'   dropped and names them.
#'
#' @return A list with two elements: \code{backbone}, the pruned backbone tree,
#'   and \code{trees}, the pruned comparison trees, retaining their names and
#'   their original single-tree or pooled structure. The list carries a
#'   \code{"dropped"} attribute naming the taxa that were removed.
#'
#' @seealso \code{\link{check_taxa}} to see which taxa are missing and from
#'   where, before deciding to prune.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' trees    <- read_trees("path/to/your/trees")
#' backbone <- trees[["iqtree"]]
#' others   <- trees[names(trees) != "iqtree"]
#'
#' # See what is missing, and from which comparison trees
#' ok <- check_taxa(backbone, others)
#' attr(ok, "diagnostics")
#'
#' # Reduce everything to the shared taxa
#' shared <- prune_to_shared(backbone, others)
#' attr(shared, "dropped")
#'
#' m <- node_presence_matrix(shared$backbone, shared$trees)
#' }
prune_to_shared <- function(backbone, trees, verbose = TRUE) {

  if (!inherits(backbone, "phylo") && !inherits(backbone, "multiPhylo")) {
    stop(
      "`backbone` must be a `phylo` object, not ",
      paste(class(backbone), collapse = "/"), ".",
      call. = FALSE
    )
  }
  if (!inherits(trees, c("list", "multiPhylo"))) {
    stop(
      "`trees` must be a list of trees, as returned by `read_trees()`.",
      call. = FALSE
    )
  }
  if (length(trees) == 0L) {
    stop("`trees` is empty.", call. = FALSE)
  }

  bb_taxa <- tree_taxa(backbone)

  # Taxa shared by the backbone and every comparison tree.
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
