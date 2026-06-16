#' Parse node support values from internal nodes of a phylogenetic tree.
#'
#' Extracts one or two node support values from the node labels of a
#' phylogenetic tree. A single value (such as an IQ-TREE bootstrap) fills the
#' first column and leaves the second as NA. A compound label separated by
#' \code{sep} (such as an ASTRAL local posterior probability and quartet support
#' value) fills both columns. Empty or non-numeric labels become NA. Trees
#' without node labels return a data frame of NAs, one row per internal node.
#'
#' @param tree A phylogenetic tree of class "phylo", with node labels stored
#'   in tree$node.label.
#'
#' @param sep A single character string giving the delimiter that separates
#'   compound support values. Defaults to "/".
#'
#' @param round Integer or NULL. If supplied, support values are rounded to
#'   this many decimal places. Defaults to NULL.
#'
#' @return A data frame with two columns, support_1 and support_2, and one
#'   row per internal node of the tree.
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
  # 3. Pull out the node labels
  lbl <- tree$node.label
  # 4. If there are no labels, return all NA, one row per internal node
  if (is.null(lbl)) {
    return(data.frame(
      support_1 = rep(NA_real_, ape::Nnode(tree)),
      support_2 = rep(NA_real_, ape::Nnode(tree))
    ))
  }
  # 5. If labels are shorter than the node count, pad with NA
  if (length(lbl) < ape::Nnode(tree)) {
    lbl <- c(lbl, rep(NA_character_, ape::Nnode(tree) - length(lbl)))
  }
  # 6. Split each label on the separator
  spl <- strsplit(lbl, sep, fixed = TRUE)
  # 7. Turn each split label into exactly two numbers
  supp_mat <- t(
    vapply(
      spl,
      function(x) {
        x <- x[nzchar(x) & !is.na(x)]
        if (length(x) == 0) {
          c(NA_real_, NA_real_)
        } else if (length(x) == 1) {
          c(suppressWarnings(as.numeric(x[[1]])), NA_real_)
        } else {
          suppressWarnings(as.numeric(x[1:2]))
        }
      },
      numeric(2)
    )
  )
  # 8. Build the output data frame
  result <- data.frame(
    support_1 = supp_mat[, 1],
    support_2 = supp_mat[, 2]
  )
  # 9. Round if requested
  if (!is.null(round)) {
    result$support_1 <- base::round(result$support_1, round)
    result$support_2 <- base::round(result$support_2, round)
  }
  result
}
