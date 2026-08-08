#' Parse support values from internal nodes of a phylogenetic tree
#'
#' @description
#' Extracts the node labels of one tree and returns their support values as a
#' data frame: one row per internal node, three columns (`support_1`,
#' `support_2`, `support_3`). How a label fills those columns:
#'
#' * A single value (such as an IQ-TREE bootstrap) fills `support_1`; the other
#'   two are `NA`.
#' * A compound label split by `sep` (default `/`) fills one column per value:
#'   * Two for IQ-TREE SH-aLRT/UFBoot2 or ASTRAL pp1/pp2.
#'   * Three for IQ-TREE SH-aLRT/UFBoot2/aBayes or ASTRAL pp1/pp2/pp3.
#' * Empty or non-numeric pieces become `NA`.
#' * A tree with no node labels returns a data frame of all `NA`.
#'
#' @param tree A phylogenetic tree of class `"phylo"`, with node labels stored
#'   in `tree$node.label`.
#'
#' @param sep A single character string giving the delimiter that separates
#'   compound support values. Defaults to `"/"`.
#'
#' @param digits Integer or `NULL`. If supplied, support values are rounded to
#'   this many decimal places. Defaults to `NULL`.
#'
#' @return A data frame with three columns, `support_1`, `support_2` and
#'   `support_3`, and one row per internal node of the tree. Columns not present
#'   in a given label are `NA`.
#'
#' @keywords internal
node_support <- function(tree,
                         sep   = "/",
                         digits = NULL) {
  if (!inherits(tree, "phylo")) {
    stop(
      "`tree` must be a phylogenetic tree of class \"phylo\".",
      call. = FALSE
    )
  }
  if (!is.null(digits) && (!is.numeric(digits) || length(digits) != 1)) {
    stop(
      "`digits` must be a single integer or NULL.",
      call. = FALSE
    )
  }
  n_node <- ape::Nnode(tree)
  n_col  <- 3L
  lbl <- tree$node.label
 # No labels: return all NA, one row per internal node
  if (is.null(lbl)) {
    out <- as.data.frame(matrix(NA_real_, nrow = n_node, ncol = n_col))
    names(out) <- paste0("support_", seq_len(n_col))
    return(out)
  }
  # Pad short label vectors so there is one entry per internal node
  if (length(lbl) < n_node) {
    lbl <- c(lbl, rep(NA_character_, n_node - length(lbl)))
  }
  # Split each label on the separator
  spl <- strsplit(lbl, sep, fixed = TRUE)
  #Coerce each split label to exactly n_col numbers, padding with NA
  supp_mat <- t(
    vapply(
      spl,
      function(x) {
        x <- x[nzchar(x) & !is.na(x)]        # drop empty and NA pieces
        vals <- suppressWarnings(as.numeric(x))
        length(vals) <- n_col                # truncate or pad with NA
        vals
      },
      numeric(n_col)
    )
  )
  #Build the output data frame
  result <- as.data.frame(supp_mat)
  names(result) <- paste0("support_", seq_len(n_col))
  #Round if requested
  if (!is.null(digits)) {
    result[] <- lapply(result, function(v) round(v, digits))
  }
  result
}
