#' Check tip label consistency between backbone and a set of trees
#'
#' Compares each comparison tree against the backbone to confirm that all
#' analyses were run on the same set of taxa. A clade cannot be evaluated in an
#' analysis that lacks its taxa, so this check should be run before a phylorug
#' is built.
#'
#' @details
#' This function reports; it does not modify the trees. If the backbone is
#' present in \code{trees} it is silently removed before comparison. Three
#' outcomes are
#'
#' \describe{
#'   \item{identical}{The comparison tree shares the backbone's taxa exactly.}
#'   \item{superset}{The comparison tree contains every backbone taxon, plus
#'     some others.
#'   \item{superset}{The comparison trees contains every backbone taxon, plus
#'     some extra. Every backbone clade can still be evaluated, so no action is
#'     needed; the extra taxa are ignored.}
#'   \item{missing}{The comparison trees lacks one or more backbone taxa. Any
#'     backbone clade containing a missing taxon \emph{might produce false
#'     absence} in that tree. Scoring such a clade as absent would report a
#'     rejection where no question was ever put. Use [prune_to_shared()] to
#'     reduce the backbone and the comparison trees to their common taxa.}
#' }
#'
#' Where a comparison tree is a pool of several equally optimal trees, every
#' tree in the pool is checked. A pool whose trees disagree about their own taxa
#' is a data problem, and is an error.
#'
#' @param backbone The reference tree, as a \code{phylo} object. The rug is
#'   drawn on this tree, and its clades are searched for in each comparison
#'   tree.
#'
#' @param trees A named list of \code{phylo} or \code{multiPhylo} objects, one
#'   per comparison, as returned by [read_trees()].If the backbone is present in the
#'   list it is removed automatically before comparison, so you can pass the
#'   full list from [read_trees()] without subsetting first.
#'
#' @param verbose Logical. If \code{TRUE} (default), reports the outcome and
#'   names any mismatched comparison trees.
#'
#' @return A single logical value: \code{TRUE} if every comparison tree shares
#'   the backbone's taxa exactly. The value carries a \code{"diagnostics"}
#'   attribute, a data frame with one row per comparison tree giving its status
#'   (\code{"identical"}, \code{"superset"} or \code{"missing"}) and the taxa
#'   missing from, or extra to, the backbone. Diagnostics are attached whatever
#'   \code{verbose} is set to.
#'
#' @seealso \code{\link{read_trees}} to load the trees;
#'   \code{\link{prune_to_shared}} to reduce them to a common taxon set.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' trees    <- read_trees("path/to/your/trees")
#' backbone <- trees[["iqtree"]]
#'
#' ok <- check_taxa(backbone, trees)
#'
#'
#' # The full report, whether or not anything went wrong
#' report <- attr(ok, "diagnostics")
#' report
#'
#' # Which comparison trees cannot evaluate every backbone clade?
#' report[report$status == "missing", ]
#' }
check_taxa <- function(backbone, trees, verbose = TRUE) {

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
    stop("`trees` is empty. Supply at least one comparison tree.",
         call. = FALSE)
  }
  if (any(vapply(trees, is.null, logical(1)))) {
    stop(
      "`trees` contains NULL elements. Ensure every tree was read successfully.",
      call. = FALSE
    )
  }

  # A compressed multiPhylo stores tip labels once, as an attribute, rather than
  # per tree. Taxon sets cannot be compared until it is decompressed.
  if (inherits(trees, "multiPhylo") && !is.null(attr(trees, "TipLabel"))) {
    stop(
      "`trees` is a compressed multiPhylo object (tip labels stored as an ",
      "attribute, not per tree). Decompress it first with ",
      "`ape::uncompressTipLabel(trees)`.",
      call. = FALSE
    )
  }

  # Remove the backbone from `trees` if it was passed in, identified by
  # object identity. This lets the user pass the full read_trees() output
  # without subsetting first.
  is_backbone <- vapply(trees, identical, logical(1L), backbone)
  if (any(is_backbone)) {
    trees <- trees[!is_backbone]
  }
  if (length(trees) == 0L) {
    stop(
      "`trees` contains only the backbone. Supply at least one comparison tree.",
      call. = FALSE
    )
  }

  bb_taxa <- pool_taxa(backbone, name = "backbone")
  nm <- names(trees)
  if (is.null(nm)) {
    nm <- paste0("tree_", seq_along(trees))
  }

  status <- character(length(trees))
  miss   <- character(length(trees))
  extra  <- character(length(trees))
  n_taxa <- integer(length(trees))

  for (i in seq_along(trees)) {
    taxa <- pool_taxa(trees[[i]], name = nm[i])

    m <- setdiff(bb_taxa, taxa)   # in backbone, absent from this comparison
    e <- setdiff(taxa, bb_taxa)   # in this comparison, absent from backbone

    status[i] <- if (length(m) == 0L && length(e) == 0L) {
      "identical"
    } else if (length(m) == 0L) {
      "superset"
    } else {
      "missing"
    }

    miss[i]   <- paste(m, collapse = ", ")
    extra[i]  <- paste(e, collapse = ", ")
    n_taxa[i] <- length(taxa)
  }

  diagnostics <- data.frame(
    comparison = nm,
    status     = status,
    n_taxa     = n_taxa,
    missing    = miss,
    extra      = extra,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  ok <- all(status == "identical")

  if (verbose) {
    if (ok) {
      message(
        "All ", length(trees), " comparison trees share the same ",
        length(bb_taxa), " taxa as the backbone."
      )
    } else {
      is_superset <- status == "superset"
      is_missing  <- status == "missing"

      if (any(is_superset)) {
        message(
          sum(is_superset), " comparison tree(s) contain taxa not in the ",
          "backbone: ", paste(nm[is_superset], collapse = ", "),
          ". Every backbone clade can still be evaluated, so no action is ",
          "needed; the extra taxa are ignored."
        )
      }
      if (any(is_missing)) {
        message(
          sum(is_missing), " comparison tree(s) are MISSING backbone taxa: ",
          paste(nm[is_missing], collapse = ", "),
          ". Any backbone clade containing a missing taxon cannot be evaluated ",
          "in those trees. Use `prune_to_shared()` to reduce the backbone and ",
          "the comparison trees to their common taxa. See ",
          "`attr(result, \"diagnostics\")` for the taxa involved."
        )
      }
    }
  }

  attr(ok, "diagnostics") <- diagnostics
  ok
}


#' Tip labels of one comparison tree
#'
#' Internal. Returns the sorted tip labels of a comparison, whether it holds one
#' tree or a pool of several. Every tree in a pool must carry the same taxa: a
#' pool whose trees disagree about their own taxon set is a data problem, not a
#' phylogenetic result, and is therefore an error rather than a mismatch.
#'
#' @noRd
pool_taxa <- function(x, name = "comparison") {
  pool <- as_pool(x)

  taxa <- lapply(pool, function(tr) sort(tr$tip.label))

  if (length(taxa) > 1L) {
    same <- vapply(taxa, identical, logical(1), taxa[[1L]])
    if (!all(same)) {
      stop(
        "The trees within \"", name, "\" do not share the same taxa. A pool ",
        "must be the equally optimal trees from one search, all run on the ",
        "same data.",
        call. = FALSE
      )
    }
  }

  taxa[[1L]]
}
