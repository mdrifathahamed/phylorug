#' Check tip label consistency between backbone and a set of trees
#'
#' Compares each comparison tree against the backbone to confirm that all
#' analyses were run on the same set of taxa. A clade cannot be evaluated in an
#' analysis that lacks its taxa, so this check should be run before a rug is
#' built.
#'
#' @param backbone The reference tree, a \code{phylo} object. The rug is
#'   drawn on this tree, and clades are searched with the reference  with this
#'   backbone and compared with other trees.
#'
#' @param trees A named list of \code{phylo} or \code{multiPhylo} objects
#'   returned by \code{\link{read_trees}}, one per tree.
#'
#' @param verbose Logical. If \code{TRUE} (default), reports the outcome and
#'   names any mismatched analyses.
#'
#' @return A single logical value: \code{TRUE} if every analysis shares the
#'   backbone's taxa exactly. The value carries a \code{"diagnostics"}
#'   attribute, a data frame with one row per tree giving its status
#'   (\code{"identical"}, \code{"superset"} or \code{"missing"}) and the taxa
#'   missing from, or extra to, the backbone. Diagnostics are attached whatever
#'   \code{verbose} is set to.
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
#' # The full report, whether or not anything went wrong
#' report <- attr(ok, "diagnostics")
#' report
#'
#' # Which analyses cannot evaluate every backbone clade?
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
    stop("`trees` is empty. Supply at least one analysis to compare.",
         call. = FALSE)
  }
  if (any(vapply(trees, is.null, logical(1)))) {
    stop(
      "`trees` contains NULL elements. Ensure every tree was read successfully.",
      call. = FALSE
    )
  }
  if (inherits(trees, "multiPhylo") && !is.null(attr(trees, "TipLabel"))) {
    stop(
      "`trees` is a compressed multiPhylo object (tip labels stored as an ",
      "attribute, not per tree). Decompress it first with ",
      "`ape::uncompressTipLabel(trees)`.",
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

    m <- setdiff(bb_taxa, taxa)   # in backbone, absent from this analysis
    e <- setdiff(taxa, bb_taxa)   # in this analysis, absent from backbone

    status[i] <- if (length(m) == 0L && length(e) == 0L) {
      "identical"
    } else if (length(m) == 0L) {
      # Every backbone taxon is present, so every backbone clade can still be
      # evaluated. The extra taxa are irrelevant to the comparison.
      "superset"
    } else {
      # A backbone taxon is absent. Clades containing it cannot be evaluated
      # here, and must not be scored as absent.
      "missing"
    }

    miss[i]   <- paste(m, collapse = ", ")
    extra[i]  <- paste(e, collapse = ", ")
    n_taxa[i] <- length(taxa)
  }

  diagnostics <- data.frame(
    analysis = nm,
    status   = status,
    n_taxa   = n_taxa,
    missing  = miss,
    extra    = extra,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )

  ok <- all(status == "identical")

  if (verbose) {
    if (ok) {
      message(
        "All ", length(trees), " analyses share the same ", length(bb_taxa),
        " taxa as the backbone."
      )
    } else {
      is_superset <- status == "superset"
      is_missing  <- status == "missing"

      if (any(is_superset)) {
        message(
          sum(is_superset), " analysis(es) contain taxa not in the backbone: ",
          paste(nm[is_superset], collapse = ", "),
          ". Every backbone clade can still be evaluated, so no action is ",
          "needed; the extra taxa are ignored."
        )
      }
      if (any(is_missing)) {
        message(
          sum(is_missing), " analysis(es) are MISSING backbone taxa: ",
          paste(nm[is_missing], collapse = ", "),
          ". Any backbone clade containing a missing taxon cannot be evaluated ",
          "in those analyses, and must not be scored as absent. Either prune ",
          "the backbone to the shared taxa, which narrows the question being ",
          "asked, or mark the affected cells as not evaluable. See ",
          "`attr(result, \"diagnostics\")` for the taxa involved."
        )
      }
    }
  }

  attr(ok, "diagnostics") <- diagnostics
  ok
}


#' Tip labels of one analysis
#'
#' Internal. Returns the sorted tip labels of an analysis, whether it holds one
#' tree or a pool of several. Every tree in a pool must carry the same taxa: a
#' pool whose trees disagree about their own taxon set is a data problem, not a
#' phylogenetic result, and is therefore an error rather than a mismatch.
#'
#' @noRd
pool_taxa <- function(x, name = "analysis") {
  pool <- as_pool(x)

  taxa <- lapply(pool, function(tr) sort(tr$tip.label))

  if (length(taxa) > 1L) {
    same <- vapply(taxa, identical, logical(1), taxa[[1L]])
    if (!all(same)) {
      stop(
        "The trees within analysis \"", name, "\" do not share the same taxa. ",
        "A pool must be the equally optimal trees from one search, all run on ",
        "the same data.",
        call. = FALSE
      )
    }
  }

  taxa[[1L]]
}
