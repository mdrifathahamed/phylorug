#' Parse support values from internal nodes of a phylogenetic tree.
#'
#' Extracts up to three node support values from the node labels of a
#' phylogenetic tree. A single value (such as an IQ-TREE bootstrap) fills the
#' first column and leaves the rest as NA. A compound label separated by `sep`
#' fills as many columns as it carries: two for IQ-TREE SH-aLRT/UFBoot2 or
#' ASTRAL pp1/pp2, three for IQ-TREE SH-aLRT/UFBoot2/aBayes or ASTRAL
#' pp1/pp2/pp3. Empty or non-numeric labels become NA. Trees without node labels
#' return a data frame of NAs, one row per internal node.
#'
#' @param tree A phylogenetic tree of class `"phylo"`, with node labels stored
#'   in tree$node.label.
#'
#' @param sep A single character string giving the delimiter that separates
#'   compound support values. Defaults to `"/"`.
#'
#' @param round Integer or `NULL`. If supplied, support values are rounded to this
#'   many decimal places. Defaults to `NULL`.
#'
#' @return A data frame with three columns, support_1, support_2 and support_3,
#'   and one row per internal node of the tree. Columns not present in a given
#'   label are `NA`.
#'
#' @noRd
node_support <- function(tree,
                         sep   = "/",
                         round = NULL) {

  # 1. Validate the tree
  if (!inherits(tree, "phylo")) {
    stop(
      "`tree` must be a phylogenetic tree of class \"phylo\".",
      call. = FALSE
    )
  }

  # 2. Validate the round argument
  if (!is.null(round) && (!is.numeric(round) || length(round) != 1)) {
    stop(
      "`round` must be a single integer or NULL.",
      call. = FALSE
    )
  }

  n_node <- ape::Nnode(tree)
  n_col  <- 3L   # up to three support values: bootstrap / UFBoot2 / aBayes

  # 3. Pull out the node labels
  lbl <- tree$node.label

  # 4. No labels: return all NA, one row per internal node
  if (is.null(lbl)) {
    out <- as.data.frame(matrix(NA_real_, nrow = n_node, ncol = n_col))
    names(out) <- paste0("support_", seq_len(n_col))
    return(out)
  }

  # 5. Pad short label vectors so there is one entry per internal node
  if (length(lbl) < n_node) {
    lbl <- c(lbl, rep(NA_character_, n_node - length(lbl)))
  }

  # 6. Split each label on the separator
  spl <- strsplit(lbl, sep, fixed = TRUE)

  # 7. Coerce each split label to exactly n_col numbers, padding with NA
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

  # 8. Build the output data frame
  result <- as.data.frame(supp_mat)
  names(result) <- paste0("support_", seq_len(n_col))

  # 9. Round if requested
  if (!is.null(round)) {
    result[] <- lapply(result, function(v) round(v, round))
  }

  result
}
