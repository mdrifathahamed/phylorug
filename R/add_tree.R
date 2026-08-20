#' Add a comparison tree to an existing node presence matrix
#'
#' @description
#' Adds a new comparison tree to an existing node presence matrix without
#' rerunning [node_presence_matrix()] from scratch. The function computes
#' clade presence and support for the new tree only and attaches the results as
#' a new column in every matrix of the npm. All existing columns remain
#' unchanged.
#'
#' @details
#' Three taxon scenarios are handled automatically:
#' \describe{
#'   \item{Identical}{The new tree shares all backbone taxa exactly. Every
#'     clade is evaluated normally.}
#'   \item{Superset}{The new tree contains all backbone taxa plus extras.
#'     The extra taxa are pruned internally so that clade keys match the
#'     backbone's. Without this pruning, a backbone clade `(A, B)`
#'     would fail to match a new-tree node grouping `(A, B, X)`,
#'     producing a false absence. A message reports the pruning. The
#'     user's original tree object is not modified.}
#'   \item{Missing}{The new tree lacks one or more backbone taxa. Any
#'     backbone clade containing a missing taxon cannot be evaluated in
#'     this tree and is set to `NA` in both the presence and support
#'     matrices. These appear as red ("not computed") cells in the rug
#'     plot. A message lists the missing taxa and suggests rebuilding
#'     the npm from scratch with [prune_to_shared()] if the red cells
#'     are undesirable. Unlike [node_presence_matrix()], which errors
#'     on missing taxa, `add_tree()` is deliberately permissive:
#'     when adding a tree after the initial analysis, different taxon
#'     sampling is a legitimate use case, not a data error.}
#' }
#'
#' The new tree must be rooted. If it is a `multiPhylo` object (a pool
#' of equally optimal trees from one search), clade presence is recorded as
#' the proportion of pool trees recovering each clade, and support values
#' are averaged across trees that recovered the clade.
#'
#' If `support_type` is supplied, it is appended to the
#' `"support_type"` attribute stored in the npm object, keyed by the
#' new tree's name. This allows [plot_phylorug()] to apply the correct
#' metric-specific thresholds without the user re-declaring the full vector.
#'
#' @param npm A node presence matrix list, as returned by
#'   [node_presence_matrix()] or a previous call to `add_tree()`.
#'
#' @param backbone The backbone tree used to create `npm`. Must be the same
#'   `phylo` object (same topology and tip labels) that was passed to
#'   [node_presence_matrix()].
#'
#' @param new_tree A `phylo` or `multiPhylo` object representing the new
#'   comparison analysis to add. Must be rooted.
#'
#' @param name Character string. Column name for the new tree in the npm
#'   matrices. Required. Must be a single non-empty string and must not
#'   duplicate an existing column name.
#'
#' @param support_col Integer or vector of integers (max 3). Which support
#'   value(s) to extract from the new tree's node labels, matching the
#'   `support_col` used when the original `npm` was built. Default is `1`. For
#'   example, if the new tree is from IQ-TREE with compound labels `"80/95"`,
#'   passing `c(1, 2)` extracts both metrics.
#'
#' @param support_type Optional character string naming the support metric
#'   of the new tree (e.g. `"ufboot"`, `"lpp"`, `"jackknife"`). If supplied, it
#'   is appended to the npm's `"support_type"` attribute so that
#'   [plot_phylorug()] can apply metric-specific thresholds automatically. If
#'   `NULL`, the new tree's support values will be auto-normalized (values >1
#'   divided by 100) and binned against universal thresholds at plot time. For
#'   this reason, raw Bremer support (an unbounded integer) is not supported;
#'   use the Bremer ratio `"bremer_ratio"` instead, which is bounded between
#'   0 and 1.
#'
#' @return The updated npm list with the new tree appended as an additional
#'   column in the `presence` matrix and each `support_*` matrix. Attributes
#'   `"node_id"` and `"pool_sizes"` are updated, and `"support_type"` is
#'   extended if the argument was supplied.
#'
#' @seealso [node_presence_matrix()] to build the initial npm, [plot_phylorug()]
#'   to visualise the result, and [check_taxa()] to diagnose taxon overlap
#'   before adding.
#'
#' @export
#'
#' @examples
#' # Build an npm with 3 of the 5 shipped trees:
#' backbone <- sample_trees[["70p_uce"]]
#' others   <- sample_trees[c("70p_partition_entropy", "70p_ghost")]
#' npm      <- node_presence_matrix(backbone, others)
#' ncol(npm$presence)   # 2 comparison trees
#'
#' # Add a fourth tree without re-running the pipeline:
#' npm <- add_tree(npm, backbone,
#'                 sample_trees[["70p_ASTRAL_uce"]],
#'                 name = "ASTRAL_uce")
#' ncol(npm$presence)   # now 3
#'
#' # Add the fifth tree with a declared support type:
#' npm <- add_tree(npm, backbone,
#'                 sample_trees[["70p_ASTRAL_partition_entropy"]],
#'                 name = "ASTRAL_part_entropy",
#'                 support_type = "lpp")
#' ncol(npm$presence)   # now 4
#' attr(npm, "support_type")
#'
#' # To add a tree from an external file (not run):
#' #   new_tree <- ape::read.tree("path/to/new_analysis.treefile")
#' #   new_tree <- ape::root(new_tree, outgroup = "outgroup_sp",
#' #                         resolve.root = TRUE)
#' #   new_tree <- ape::drop.tip(new_tree, "outgroup_sp")
#' #   npm <- add_tree(npm, backbone, new_tree, name = "new_analysis",
#' #                   support_type = "ufboot")
add_tree <- function(npm, backbone, new_tree,
                     name,
                     support_col  = 1,
                     support_type = NULL) {
  # --- 1. Validate inputs ----------------------------------------------------
  if (!is.list(npm) || !("presence" %in% names(npm))) {
    stop("`npm` must be a valid node presence matrix list.", call. = FALSE)
  }
  if (!inherits(backbone, "phylo")) {
    stop("`backbone` must be a phylo object.", call. = FALSE)
  }
  if (!inherits(new_tree, "phylo") && !inherits(new_tree, "multiPhylo")) {
    stop(
      "`new_tree` must be a `phylo` or `multiPhylo` object.",
      call. = FALSE
    )
  }
  if (!ape::is.rooted(backbone)) {
    stop("`backbone` is unrooted.", call. = FALSE)
  }

  pool <- as_pool(new_tree)
  if (!all(vapply(pool, ape::is.rooted, logical(1)))) {
    stop("`new_tree` contains unrooted trees. Root with `ape::root()`.",
         call. = FALSE)
  }

  if (!is.numeric(support_col) || any(!support_col %in% seq_len(3))) {
    stop(
      "`support_col` must be an integer or integer vector with values 1, 2, ",
      "or 3.",
      call. = FALSE
    )
  }
  support_col <- as.integer(support_col)
  # --- 2. Name ---------------------------------------------------------------
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop(
      "`name` must be a single non-empty character string, e.g. ",
      "`name = \"ASTRAL_uce\"`.",
      call. = FALSE
    )
  }
  if (name %in% colnames(npm$presence)) {
    stop(
      "Column name \"", name, "\" already exists in npm. ",
      "Use a unique name.",
      call. = FALSE
    )
  }
  # --- 3. Taxon comparison ---------------------------------------------------
  bb_taxa  <- sort(backbone$tip.label)
  new_taxa <- sort(pool[[1L]]$tip.label)
  missing  <- setdiff(bb_taxa, new_taxa)
  extra    <- setdiff(new_taxa, bb_taxa)

  # Superset: prune extra taxa from the internal copy so clade keys match.
  # Without this, a backbone clade (A, B) would not match a new-tree node
  # grouping (A, B, X), producing a false absence.
  if (length(extra) > 0L) {
    message(
      "`new_tree` contains ", length(extra), " taxa not in the backbone: ",
      paste(utils::head(extra, 5L), collapse = ", "),
      if (length(extra) > 5L) ", ..." else "",
      ". These are pruned internally to match the backbone. ",
      "The original tree object is not modified."
    )
    pool <- lapply(pool, function(tr) ape::drop.tip(tr, extra))
  }

  # Missing: warn and mark affected clades as NA.
  if (length(missing) > 0L) {
    message(
      "`new_tree` is missing ", length(missing), " backbone taxa: ",
      paste(utils::head(missing, 5L), collapse = ", "),
      if (length(missing) > 5L) ", ..." else "",
      ". Clades containing these taxa are marked as not computed (NA) ",
      "and will appear as red cells in the rug. To avoid this, rebuild ",
      "the npm from scratch with all trees pruned to shared taxa using ",
      "`prune_to_shared()`."
    )
  }

  # --- 4. Build clade keys ---------------------------------------------------
  bb_keys   <- clade_keys(backbone)
  pool_keys <- lapply(pool, clade_keys)

  ntip     <- ape::Ntip(backbone)
  bb_nodes <- (ntip + 1L):(ntip + backbone$Nnode)

  # Descendant tip labels for each backbone node
  bb_tip_sets <- lapply(bb_nodes, function(nd) {
    tips <- phangorn::Descendants(backbone, nd, type = "tips")[[1L]]
    sort(backbone$tip.label[tips])
  })

  # --- 5. Compute new column -------------------------------------------------
  new_presence <- rep(NA_real_, length(bb_nodes))
  new_supports <- lapply(seq_along(support_col), function(s) {
    rep(NA_real_, length(bb_nodes))
  })

  for (i in seq_along(bb_nodes)) {
    clade_tips <- bb_tip_sets[[i]]

    # If any descendant tip is missing from new_tree, this clade cannot be
    # evaluated: leave as NA (not computed → red cell in the rug).
    if (length(missing) > 0L && any(clade_tips %in% missing)) {
      next
    }

    # Clade can be evaluated: match key normally
    key  <- bb_keys[i]
    hits <- vapply(pool_keys, function(k) key %in% k, logical(1))
    new_presence[i] <- mean(hits)

    if (any(hits)) {
      for (s in seq_along(support_col)) {
        col <- support_col[s]
        vals <- vapply(which(hits), function(w) {
          supp <- node_support(pool[[w]])
          idx  <- match(key, pool_keys[[w]])
          if (is.na(idx) || idx > nrow(supp) || col > ncol(supp)) {
            return(NA_real_)
          }
          as.numeric(supp[[col]][idx])
        }, numeric(1))

        new_supports[[s]][i] <-
          if (all(is.na(vals))) NA_real_
        else mean(vals, na.rm = TRUE)
      }
    }
  }

  # --- 6. Append to npm ------------------------------------------------------
  new_p_mat <- matrix(new_presence, ncol = 1,
                      dimnames = list(as.character(bb_nodes), name))
  npm$presence <- cbind(npm$presence, new_p_mat)

  for (s in seq_along(support_col)) {
    sname <- paste0("support_", s)
    new_s_mat <- matrix(new_supports[[s]], ncol = 1,
                        dimnames = list(as.character(bb_nodes), name))
    if (sname %in% names(npm)) {
      npm[[sname]] <- cbind(npm[[sname]], new_s_mat)
    } else {
      npm[[sname]] <- new_s_mat
    }
  }

  # --- 7. Update attributes --------------------------------------------------
  ps <- attr(npm, "pool_sizes")
  ps <- c(ps, stats::setNames(pool_size(new_tree), name))
  attr(npm, "pool_sizes") <- ps

  if (!is.null(support_type)) {
    existing_st <- attr(npm, "support_type")
    new_st      <- stats::setNames(support_type, name)
    attr(npm, "support_type") <- c(existing_st, new_st)
  }

  npm
}
