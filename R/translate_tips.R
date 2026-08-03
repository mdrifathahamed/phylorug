#' Relabel tips across a list of trees using a lookup table
#'
#' Translates tip labels of each tree in a list using a lookup table supplied
#' as a data frame. This function is optional in the phylorug workflow , use it
#' only if your tip labels are in  specimen codes, accession numbers, or any
#' other identifiers that need converting to a different format, such as
#' scientific species names.
#' If your tip labels are already in the correc format, skip this step and
#' proceed directly to [node_presence_matrix()]. Only tips with a matching entry
#' in `from_col` re replaced with the corresponding value in `to_col`. Unmatched
#' tips are left unchanged.This function modifies only the tip.label field of
#' each tree. Tree topology, branch lengths, and node labels are unaffected.
#'
#' @param trees A named list of `"phylo"` objects, returned by [read_trees()].
#'
#' @param data A data frame containing the label translation lookup table. Can
#'   be a standard `data.frame` or a `tibble`.
#'
#' @param from_col A character string specifying the column name in `data`
#'   that holds the current tip labels of the trees. Defaults to `"from"`.
#'
#' @param to_col A character string specifying the column name in `data`
#'   that holds the replacement labels. Defaults to `to`.
#'
#' @param verbose Logical. If `TRUE` (default), reports how many tip
#'   labels were translated and how many were left unchanged for each tree.
#'
#' @return A named list of `"phylo"` objects identical in structure to `trees`,
#'   with matching tip labels replaced according to `data`.Tree topology and
#'   edge lengths remain unchanged.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Load trees
#' trees <- read_trees("path/to/trees")
#'
#' # Translation table mapping specimen codes to scientific names
#' dict <- data.frame(
#'   from = c("NicorbUCE", "NicvesUCE"),
#'   to   = c("Nicrophorus_orbus", "Nicrophorus_vespillo")
#' )
#'
#' # Translate tip labels
#' trees <- translate_tips(trees,
#'                         data     = dict,
#'                         from_col = "from",
#'                         to_col   = "to")
#' }
translate_tips <- function(trees,
                           data,
                           from_col = "from",
                           to_col = "to",
                           verbose = TRUE) {
  # Validate inputs
  if (!inherits(trees, c("list", "multiPhylo"))) {
    stop(
      "`trees` must be a list of phylo objects or a multiPhylo object.",
      call. = FALSE
    )
  }
  if (!inherits(data, "data.frame")) {
    stop(
      "`data` must be a data frame.",
      call. = FALSE
    )
  }
  if (!all(c(from_col, to_col) %in% colnames(data))) {
    stop(
      "Columns \"", from_col, "\" and \"", to_col,
      "\" not found in `data`.",
      call. = FALSE
    )
  }
  if (anyDuplicated(data[[from_col]])) {
    stop(
      "Values in `from_col` must be unique.",
      call. = FALSE
    )
  }
  # Build translation dictionary
  trans_dict <- stats::setNames(
    as.character(data[[to_col]]),
    as.character(data[[from_col]])
  )
  # Translate each tree in the list
  trees <- lapply(trees, function(tr) {
    match_idx <- match(tr$tip.label, names(trans_dict))
    matched   <- !is.na(match_idx)

    tr$tip.label[matched] <- trans_dict[match_idx[matched]]

    if (verbose) {
      message("Tips translated : ", sum(matched))
      message("Tips unchanged  : ", sum(!matched))
    }
    tr
  })
  trees
}
