#' Tabulate clade recovery and support across trees
#'
#' Compares a set of phylogenetic trees against a reference topology, the
#' `"backbone"`. For each internal node of the backbone, the function asks
#' whether the same clade (the same set of tips) appears in each tree, and
#' records either its presence or its support value.
#'
#' @details
#' Every tree must include the complete set of backbone taxa. Because an
#' tree lacking a backbone taxon cannot assess the presence of a clade
#' containing that taxon,scoring such a clade as 'absent' would introduce a
#' false negative.The function therefore enforces strict taxon matching. Use
#' [check_taxa()] to diagnose the discrepancies and [prune_to_shared()] to
#' harmonize the backbone and comparing tree taxa.
#'
#' The function always computes both a presence matrix and a support matrix
#' in a single pass. The plotting function [plot_phylorug()] decide which to
#' use based on the visualisation context.
#'
#' In presence matrix (`use_support = FALSE`),For multiple equally most
#' parsimonious trees (MPTs) the cell records the clade frequency among MPTs
#' and, for a single-tree this is 1 for presence  or 0 for absence.
#'
#' In support matrix `use_support = TRUE`, the cell records the support value
#' from the tree that recovered the clade, such as bootstrap or posterior
#' probability. For MPTs, support values are averaged across the trees that
#' recovered the clade; trees that did not recover it contribute nothing. For a
#' single tree this is the support value where the clade is present, or `NA`
#' where it is absent, since a clade that is not in the tree has no node to
#' carry a support value.
#'
#' @inheritParams check_taxa
#'
#' @param support_col Integer or vector of integers (max 3). Specifies which
#'   value(s) to extract from multi-metric node labels. For example, if an
#'   IQ-TREE tree stores SH-aLRT and UFBoot as "80/95", passing `1` extracts
#'   only the first metric into a `support_1` matrix. Passing `c(1, 2)`
#'   efficiently extracts both simultaneously into `support_1` and `support_2`
#'   matrices, allowing you to easily switch between them during plotting
#'   without recalculating. Default is `1` .
#'
#' @return A named list with one row per internal backbone node and one column
#'   per comparison tree:
#'   \describe{
#'     \item{presence}{Clade presence: `1` where recovered, `0` where absent, or
#'      the proportion of pool trees recovering the clade.}
#'     \item{support_1, support_2, ...}{One matrix per value in `support_col`.
#'      Raw support values where the clade was recovered,`NA` where absent.
#'      Named in the order requested, so `support_col = c(1, 2)` produces
#'      `support_1` and `support_2`.}
#'   }
#' @export
#'
#' @examples
#' # phylorug ships `sample_trees`: a named list of 5 phylo objects (beetles),
#' # the same structure you get from read_trees(). Just type `sample_trees`.
#' # Pick one analysis as the backbone, and the rest as comparison trees:
#' backbone <- sample_trees[["70p_uce"]]
#' others   <- sample_trees[names(sample_trees) != "70p_uce"]
#'
#' # It is good practice to check the taxa line up before building the matrix:
#' check_taxa(backbone, others)
#'
#' # For a presence/absence matrix, just pass the backbone and comparison trees
#' # (support_col defaults to 1):
#' npmatrix <- node_presence_matrix(backbone, others)
#' npmatrix$presence     # clade presence/absence
#' npmatrix$support_1    # support values (from column 1 of the node labels)
#'
#' # To extract several support metrics at once (e.g. SH-aLRT and UFBoot2 stored
#' # as "80/95"), pass their column positions to support_col:
#' rugmt <- node_presence_matrix(backbone, others, support_col = c(1, 2))
#' rugmt$support_1       # first metric
#' rugmt$support_2       # second metric
node_presence_matrix <- function(backbone,
                                 trees,
                                 support_col = 1) {

  if (!inherits(backbone, "phylo")) {
    stop(
      "`backbone` must be a phylogenetic tree of class \"phylo\".",
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
  if (!is.numeric(support_col) || any(!support_col %in% seq_len(3))) {
    stop(paste0("`support_col` must be an integer or integer vector ",
                "with values 1, 2, or 3."), call. = FALSE)
  }
  support_col <- as.integer(support_col)
  nm <- names(trees)
  if (is.null(nm)) {
    nm <- paste0("tree_", seq_along(trees))
  }

  # Backbone check if its rooted or not
  if (!ape::is.rooted(backbone)) {
    stop(
      "`backbone` is unrooted. Clade comparison requires a rooted tree: the ",
      "descendants of a node are undefined without a root. Root it with ",
      "`ape::root()`.",
      call. = FALSE
    )
  }

  # Comparison trees check if its rooted or not
  unrooted <- vapply(seq_along(trees), function(j) {
    pool <- as_pool(trees[[j]])
    !all(vapply(pool, ape::is.rooted, logical(1)))
  }, logical(1))

  if (any(unrooted)) {
    stop(
      length(which(unrooted)), " comparison tree(s) contain unrooted trees: ",
      paste(nm[unrooted], collapse = ", "),
      ". Clade comparison requires rooted trees. Root them with ",
      "`ape::root(tree, outgroup = \"...\", resolve.root = TRUE)`.",
      call. = FALSE
    )
  }
  # every trees must carry the same  set of taxa
  bb_taxa <- sort(backbone$tip.label)
  assert_shared_taxa(bb_taxa, trees, nm)
  #setup the empty matrix
  ntip     <- ape::Ntip(backbone)
  bb_nodes <- (ntip + 1L):(ntip + backbone$Nnode)
  bb_keys  <- clade_keys(backbone)

  # Presence matrix: 1 = clade recovered, 0 = absent, proportion for pools
  presence_matrix <- matrix(
    NA_real_,
    nrow     = length(bb_nodes),
    ncol     = length(trees),
    dimnames = list(as.character(bb_nodes), nm)
  )

  # One support matrix per requested column
  support_matrices <- lapply(seq_along(support_col), function(s) {
    matrix(
      NA_real_,
      nrow     = length(bb_nodes),
      ncol     = length(trees),
      dimnames = list(as.character(bb_nodes), nm)
    )
  })

  #Double loop begins -to fill each cell of the matrix
  #outer loop: j iterates over trees(columns)
  #inner loop: i iterates over backbone nodes(rows)
  # Double loop begins - to fill each cell of the matrix
  for (j in seq_along(trees)) {
    pool      <- as_pool(trees[[j]])
    pool_keys <- lapply(pool, clade_keys)

    for (i in seq_along(bb_nodes)) {
      key  <- bb_keys[i]
      hits <- vapply(pool_keys, function(k) key %in% k, logical(1))

      # 1. Fill presence
      presence_matrix[i, j] <- mean(hits)

      # 2. Fill support (This must be INSIDE the i loop)
      if (any(hits)) {
        for (s in seq_along(support_col)) {
          col <- support_col[s]
          vals <- vapply(which(hits), function(w) {
            supp <- node_support(pool[[w]])
            idx  <- match(key, pool_keys[[w]]) # Use current key
            if (is.na(idx) || idx > nrow(supp) || col > ncol(supp)) {
              return(NA_real_)
            }
            as.numeric(supp[[col]][idx])
          }, numeric(1))

          support_matrices[[s]][i, j] <-
            if (all(is.na(vals))) NA_real_
            else mean(vals, na.rm = TRUE)
        }
      } else {
        # If no hits, explicitly set support to NA
        for (s in seq_along(support_col)) {
          support_matrices[[s]][i, j] <- NA_real_
        }
      }
    }
  }
  #Attach metadata and return
  pool_sizes <- vapply(trees, pool_size, integer(1))

  # Build named list: presence + support_1, support_2, ...
  support_named <- stats::setNames(
    support_matrices,
    paste0("support_", seq_along(support_col))
  )

  result <- c(list(presence = presence_matrix), support_named)

  attr(result, "node_id")    <- bb_nodes
  attr(result, "pool_sizes") <- pool_sizes

  result
}
#' Refuse to proceed unless every trees carries the backbone's taxa
#'
#' @noRd
assert_shared_taxa <- function(bb_taxa, trees, nm) {
  missing_taxa <- lapply(seq_along(trees), function(j) {
    taxa <- sort(as_pool(trees[[j]])[[1L]]$tip.label)
    setdiff(bb_taxa, taxa)
  })

  bad <- which(lengths(missing_taxa) > 0L)
  if (length(bad) == 0L) {
    return(invisible(TRUE))
  }

  detail <- vapply(bad, function(j) {
    tx <- missing_taxa[[j]]
    paste0(
      "  ", nm[j], ": missing ", length(tx), " taxa (",
      paste(utils::head(tx, 5L), collapse = ", "),
      if (length(tx) > 5L) ", ..." else "", ")"
    )
  }, character(1))

  stop(
    length(bad), " tree(s) do not contain every backbone taxon:\n",
    paste(detail, collapse = "\n"),
    "\nA clade containing a missing taxon cannot be evaluated in that ",
    "tree, so phylorug will not proceed. Run `check_taxa()` for the full ",
    "report, then `prune_to_shared()` to reduce the backbone and the ",
    "comparison trees to their common taxa.",
    call. = FALSE
  )
}

#' Canonical clade keys for one tree
#'
#' @noRd
clade_keys <- function(tree) {
  ntip  <- ape::Ntip(tree)
  nodes <- (ntip + 1L):(ntip + tree$Nnode)

  vapply(nodes, function(nd) {
    tips <- phangorn::Descendants(tree, nd, type = "tips")[[1L]]
    paste(sort(tree$tip.label[tips]), collapse = "|")
  }, character(1))
}
