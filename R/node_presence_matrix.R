#' Tabulate clade recovery and support across trees
#'
#' Compares a set of phylogenetic trees against a reference topology, the
#' `"backbone"`. For each internal node of the backbone, the function asks
#' whether the same clade appears in each tree, and records either its presence
#' or its support value.
#'
#' @details
#' This is the core function of the phylorug pipeline. It can be called
#' directly after [read_trees()], [check_taxa()] is a useful diagnosis but not
#' a prerequisite, as [node_presence_matrix()] requires perfectly matched taxon
#' sets and the pipline can not proceed if any mismatches in taxa set are
#' detected between the backbone and comparison trees.
#'
#' Every tree must include the complete set of backbone taxa. Because a tree
#' lacking a backbone taxon cannot assess the presence of a clade containing
#' that taxon, scoring such a clade as 'absent' would introduce a false
#' negative. The function therefore enforces strict taxon matching. Use
#' [check_taxa()] to diagnose the discrepancies and [prune_to_shared()] to
#' harmonize the backbone and comparing tree taxa.
#'
#' The function always computes both a presence matrix and a support matrix
#' in a single pass. The plotting function [plot_phylorug()] decide which to
#' use based on the visualisation context.
#'
#' The `presence` matrix records clade recovery: `1` where the bipartition
#' is found, `0` where absent. For multiple equally most parsimonious trees
#' (MPTs), the behaviour depends on `pool_threshold`: the default (`1.0`,
#' strict consensus) requires the clade to appear in every pool tree;
#' `0.5` applies majority rule; `0` records the raw proportion.
#'
#' Each `support_*` matrix records the raw support value from the node label
#' of the tree that recovered the clade, such as bootstrap percentage or
#' posterior probability. Support extraction is only meaningful for single
#' trees (`phylo`), not for pools of equally optimal trees (`multiPhylo`).
#' Averaging node labels across MPTs is not scientifically valid: support
#' values come from resampling analyses (bootstrap, jackknife), which
#' produce a separate consensus tree that should be passed to phylorug as
#' a single `phylo` comparison tree. If a pool is supplied, only presence
#' mode is meaningful; support cells for that column will be `NA`.
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
#' @param support_type Optional named character vector mapping comparison trees
#'   to their support metrics (e.g. `"ufboot"`, `"sh_alrt"`, `"lpp"`,
#'   `"jackknife"`). Stored as an attribute of the returned list and used
#'   automatically by [plot_phylorug()] when `mode = "support"`. If omitted,
#'   [plot_phylorug()] will auto-normalize support values and apply universal
#'   thresholds.
#'
#' @param pool_threshold Numeric between 0 and 1. Controls how pools of equally
#'   optimal trees (e.g. from TNT or PAUP*) are scored. Default `1.0`
#'   (strict consensus): a clade must appear in every pool tree to be scored as
#'   present. Set to `0.5` for majority rule (>50% of pool trees). Set to `0` to
#'   record the raw proportion (gradient cells in the rug). See Simmons &
#'   Freudenstein (2011) for why strict consensus is the recommended default for
#'   parsimony analyses. This parameter is only applicable when evaluating
#'   parsimony-based `multiPhylo` objects (e.g., Most Parsimonious Trees
#'   generated via TNT or similar software).
#'
#' @return A named list with one row per internal backbone node and one column
#'   per comparison tree:
#'   \describe{
#'     \item{presence}{Clade presence: `1` where recovered, `0` where absent.
#'      For pools, the value depends on `pool_threshold`: strict consensus
#'      (default) gives `1` or `0`; `pool_threshold = 0` gives the raw
#'      proportion.}
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
#' # Optional: diagnose taxon overlap before building the matrix.
#' # This is not required — node_presence_matrix() enforces matching internally.
#' check_taxa(backbone, others)
#'
#' # --- Presence/absence matrix -----------------------------------------------
#' # Just pass the backbone and comparison trees (support_col defaults to 1).
#' # For pools of MPTs, pool_threshold defaults to 1.0 (strict consensus):
#' npmatrix <- node_presence_matrix(backbone, others)
#' npmatrix$presence     # clade presence/absence (1/0)
#' npmatrix$support_1    # support values (from column 1 of the node labels)
#'
#' # --- Multiple support metrics at once --------------------------------------
#' # IQ-TREE stores SH-aLRT and UFBoot2 as "80/95". Pass their column positions
#' # to support_col to extract both in a single pass:
#' rugmt <- node_presence_matrix(backbone, others, support_col = c(1, 2))
#' rugmt$support_1       # first metric (SH-aLRT)
#' rugmt$support_2       # second metric (UFBoot2)
#'
#' # --- With support_type (stored in npm, used by plot_phylorug) --------------
#' # Declaring support_type lets plot_phylorug() apply metric-specific
#' # thresholds automatically. If omitted, universal 0-1 thresholds are used.
#' support_type <- c(
#'   "70p_ASTRAL_partition_entropy" = "lpp",
#'   "70p_ASTRAL_uce"               = "lpp",
#'   "70p_ghost"                    = "ufboot",
#'   "70p_partition_entropy"        = "ufboot"
#' )
#' npm_typed <- node_presence_matrix(backbone, others,
#'                                   support_type = support_type)
#' attr(npm_typed, "support_type")   # stored as an attribute
node_presence_matrix <- function(backbone,
                                 trees,
                                 support_col    = 1,
                                 support_type   = NULL,
                                 pool_threshold = 1.0) {

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
  if (!is.numeric(pool_threshold) || length(pool_threshold) != 1L ||
      pool_threshold < 0 || pool_threshold > 1) {
    stop(
      "`pool_threshold` must be a single number between 0 and 1.",
      call. = FALSE
    )
  }
  nm <- names(trees)
  if (is.null(nm)) {
    nm <- paste0("tree_", seq_along(trees))
  }
  if (!ape::is.rooted(backbone)) {
    stop(
      "`backbone` is unrooted. Clade comparison requires a rooted tree: the ",
      "descendants of a node are undefined without a root. Root it with ",
      "`ape::root()`.",
      call. = FALSE
    )
  }
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

  bb_taxa <- sort(backbone$tip.label)
  assert_shared_taxa(bb_taxa, trees, nm)

  ntip     <- ape::Ntip(backbone)
  bb_nodes <- (ntip + 1L):(ntip + backbone$Nnode)
  bb_keys  <- clade_keys(backbone)

  presence_matrix <- matrix(
    NA_real_,
    nrow     = length(bb_nodes),
    ncol     = length(trees),
    dimnames = list(as.character(bb_nodes), nm)
  )

  support_matrices <- lapply(seq_along(support_col), function(s) {
    matrix(
      NA_real_,
      nrow     = length(bb_nodes),
      ncol     = length(trees),
      dimnames = list(as.character(bb_nodes), nm)
    )
  })

  for (j in seq_along(trees)) {
    pool      <- as_pool(trees[[j]])
    pool_keys <- lapply(pool, clade_keys)

    for (i in seq_along(bb_nodes)) {
      key  <- bb_keys[i]
      hits <- vapply(pool_keys, function(k) key %in% k, logical(1))

      freq <- mean(hits)
      presence_matrix[i, j] <- if (pool_threshold == 0) {
        freq
      } else {
        if (freq >= pool_threshold) 1 else 0
      }

      if (any(hits) && length(pool) == 1L) {
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

        for (s in seq_along(support_col)) {
          support_matrices[[s]][i, j] <- NA_real_
        }
      }
    }
  }
  pool_sizes <- vapply(trees, pool_size, integer(1))

  support_named <- stats::setNames(
    support_matrices,
    paste0("support_", seq_along(support_col))
  )
  result <- c(list(presence = presence_matrix), support_named)
  attr(result, "node_id")      <- bb_nodes
  attr(result, "pool_sizes")   <- pool_sizes
  attr(result, "support_type")   <- support_type
  attr(result, "pool_threshold") <- pool_threshold
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
