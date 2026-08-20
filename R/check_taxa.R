#' Diagnose taxon consistency between backbone and comparison trees
#'
#' For each comparison tree, checks whether its tip labels (taxon) match the
#' backbone's exactly, contain extra labels, or are missing some. A tree
#' missing a tip label of the backbone cannot be scored for any clade
#' containing that tip. We recommend running this diagnostic before building
#' the matrix, especially when combining trees from different studies or
#' pipelines. However, it is not a mandatory step [node_presence_matrix()]
#' enforces taxon matching internally and will error with a clear message if any
#' comparison tree is missing backbone taxa. If you are certain all trees share
#' the same taxon set, you can proceed directly to [node_presence_matrix()].
#'
#' @details
#' This function reports; it does not modify the trees. If the backbone is
#' present in comparison `trees` it is silently removed before comparison. Three
#' outcomes are
#'
#' \describe{
#'   \item{identical}{All comparison trees share the backbone's taxa exactly.}
#'
#'   \item{superset}{One or more comparison trees contain every backbone taxon,
#'     plus some extra. Every backbone clade can still be evaluated, so no
#'     action is needed; the extra taxa can be ignored.}
#'
#'   \item{missing}{One or more comparison trees lack one or more backbone
#'     taxa. Any backbone clade containing a missing taxon might produce false
#'     absence in that tree. Scoring such a clade as absent would report a
#'     rejection where no question was ever put. Use [prune_to_shared()] to
#'     reduce the backbone and the comparison trees to their common taxa.}
#' }
#'
#' Where a comparison tree is a pool of several equally optimal trees, every
#' tree in the pool is checked. A pool whose trees disagree about their own taxa
#' is a data problem, and is an error.
#'
#' @param backbone The reference tree, as a `"phylo"` object. The rug is
#'   drawn on this tree, and its clades are searched for in each comparison
#'   tree.
#'
#' @param trees A named list of `"phylo"` or `"multiPhylo"` objects, one per
#'   comparison, as returned by [read_trees()].If the backbone is present in the
#'   list it is removed automatically before comparison, so you can pass the
#'   full list from [read_trees()] without subsetting first.
#'
#' @param verbose Logical. If `TRUE` (default), reports the outcome and
#'   names any mismatched comparison trees.
#'
#' @return A single logical value: `TRUE` if every comparison tree shares the
#'   backbone's taxa exactly. The value carries a `"diagnostics"` attribute, a
#'   data frame with one row per comparison tree giving its status
#'   (`"identical"`, `"superset"` or `"missing"`) and the taxa missing from, or
#'   extra to, the backbone. Diagnostics are attached whatever `verbose` is set
#'   to.
#'
#' @seealso [read_trees()] to load the trees and [prune_to_shared()] to reduce
#'   them to a common taxon set.
#'
#' @export
#'
#' @examples
#' # sample_trees is a named list of 5 phylo objects shipped with phylorug,
#' # the same structure you get from read_trees(). Pick one analysis as the
#' # backbone, and the rest as comparison trees:
#' backbone <- sample_trees[["70p_uce"]]
#' others   <- sample_trees[names(sample_trees) != "70p_uce"]
#'
#' # Diagnose taxon consistency (verbose = TRUE by default reports the outcome):
#' result <- check_taxa(backbone, others)
#'
#' # If any tree is missing or has extra taxa, inspect the full report:
#' attr(result, "diagnostics")

check_taxa <- function(backbone, trees, verbose = TRUE) {

  if (!inherits(backbone, "phylo") && !inherits(backbone, "multiPhylo")) {
    stop(
      "`backbone` must be a `phylo` object, not ",
      paste(class(backbone), collapse = "/"), ".",
      call. = FALSE
    )
  }
  if (!inherits(trees, "list")) {
    stop(
      "`trees` must be a list of `phylo` and/or `multiPhylo` ",
      "objects, one element per analysis, as returned by ",
      "`read_trees()`. A single `multiPhylo` is a pool of ",
      "trees from ONE analysis and should be wrapped in a list.",
      call. = FALSE
    )
  }
  if (length(trees) == 0L) {
    stop("`trees` is empty. Supply at least one comparison tree.",
         call. = FALSE)
  }
  if (any(vapply(trees, is.null, logical(1)))) {
    stop(
      "`trees` contains NULL elements. Ensure every tree was ",
      "read successfully.",
      call. = FALSE
    )
  }
  valid_element <- vapply(trees, function(x) {
    inherits(x, "phylo") || inherits(x, "multiPhylo")
  }, logical(1))
  if (!all(valid_element)) {
    stop(
      "All elements of `trees` must be `phylo` or `multiPhylo` objects. ",
      "Invalid elements at positions: ",
      paste(which(!valid_element), collapse = ", "), ".",
      call. = FALSE
    )
  }
  is_compressed <- vapply(trees, function(x) {
    inherits(x, "multiPhylo") && !is.null(attr(x, "TipLabel"))
  }, logical(1))
  if (any(is_compressed)) {
    stop(
      "Element(s) at position(s) ",
      paste(which(is_compressed), collapse = ", "),
      " are compressed multiPhylo objects (tip labels stored as an ",
      "attribute, not per tree). Decompress with ",
      "`ape::uncompressTipLabel()` first.",
      call. = FALSE
    )
  }
  is_backbone <- vapply(trees, identical, logical(1L), backbone)
  if (any(is_backbone)) {
    trees <- trees[!is_backbone]
  }
  if (length(trees) == 0L) {
    stop(
      "`trees` contains only the backbone. Supply at least ",
      "one comparison tree.",
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

    m <- setdiff(bb_taxa, taxa)
    e <- setdiff(taxa, bb_taxa)

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
          sum(is_missing), " comparison tree(s) are MISSING ",
          "backbone taxa: ",
          paste(nm[is_missing], collapse = ", "),
          ". Any backbone clade containing a missing taxon ",
          "cannot be evaluated in those trees. Use ",
          "`prune_to_shared()` to reduce the backbone and ",
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
