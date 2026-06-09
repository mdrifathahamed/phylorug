#' Translate tip labels across a list of phylogenetic trees
#'
#' Renames the tip labels of each tree in a list using a lookup table supplied
#' as a data frame. This function is optional in the phylorug workflow — use it
#' only if your tip labels are specimen codes, accession numbers, or any other
#' identifiers that need converting to a different format, such as scientific
#' species names. If your tip labels are already in the correct format, skip
#' this step and proceed directly to \code{\link{node_presence_matrix}}.
#'
#' Only tips with a matching entry in \code{from_col} are renamed. Unmatched
#' tips are left unchanged. Both \code{from_col} and \code{to_col} must be
#' character columns, not factors.
#'
#' @param trees A named list of \code{"phylo"} objects, as returned by
#'   \code{\link{read_trees}}.
#'
#' @param data A data frame containing the label translation lookup table. Can
#'   be a standard \code{data.frame} or a \code{tibble}.
#'
#' @param from_col A character string specifying the column name in \code{data}
#'   that holds the current tip labels of the trees. Defaults to \code{"from"}.
#'
#' @param to_col A character string specifying the column name in \code{data}
#'   that holds the replacement labels. Defaults to \code{"to"}.
#'
#' @param verbose Logical. If \code{TRUE} (default), reports how many tip
#'   labels were translated and how many were left unchanged for each tree.
#'
#' @return A named list of \code{"phylo"} objects identical in structure to
#'   \code{trees}, with matching tip labels replaced according to \code{data}.
#'   Tree topology and edge lengths remain unchanged.
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
                           to_col   = "to",
                           verbose  = TRUE) {
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
  # Build translation dictionary
  # as.character() protects against factor columns
  trans_dict <- stats::setNames(
    as.character(data[[to_col]]),
    as.character(data[[from_col]])
  )
  # Translate each tree in the list
  trees <- lapply(trees, function(tr) {
    tips_to_translate <- tr$tip.label %in% names(trans_dict)
    tr$tip.label[tips_to_translate] <-
      trans_dict[tr$tip.label[tips_to_translate]]
    if (verbose) {
      message("Tips translated : ", sum(tips_to_translate))
      message("Tips unchanged  : ", sum(!tips_to_translate))
    }
    tr
  })
  trees
}
