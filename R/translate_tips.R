#' Relabel tips across a list of trees using a lookup table
#'
#' @description
#' Translates tip labels of each tree in a list using a lookup table supplied
#' as a data frame. This function is optional in the phylorug workflow, use it
#' only if your tip labels are specimen codes, accession numbers, or any other
#' identifiers that need converting to a different format, such as scientific
#' species names.
#'
#' If your tip labels are already in the correct format, skip this step and
#' proceed directly to [node_presence_matrix()]. Only tips with a matching entry
#' in `from_col` are replaced with the corresponding value in `to_col`.
#' Unmatched tips are left unchanged. This function modifies only the
#' `tip.label` field of each tree. Tree topology, branch lengths, and node
#' labels are unaffected.
#'
#' @param trees A named list of `"phylo"` and/or `"multiPhylo"` objects, one
#'   element per analysis, as returned by [read_trees()]. A `"multiPhylo"`
#'   element is a pool of tied-optimal trees (POY/TNT/PAUP*) treated as a
#'   single analysis.
#'
#' @param data A data frame containing the label translation lookup table. Can
#'   be a standard `data.frame` or a `tibble`.
#'
#' @param from_col A character string naming the column in `data` that holds the
#'   current tip labels (for example `"specimen_code"`). Required.
#'
#' @param to_col A character string naming the column in `data` that holds the
#'   replacement labels (for example `"scientific_name"`). Required.
#'
#' @param verbose Logical. If `TRUE` (default), reports how many tip
#'   labels were translated and how many were left unchanged for each tree.
#'
#' @return A named list of `"phylo"` objects identical in structure to `trees`,
#'   with matching tip labels replaced according to `data`. Tree topology and
#'   edge lengths remain unchanged.
#'
#' @export
#'
#' @examples
#' # Example trees and a lookup table ship inside phylorug, so we first read
#' # them into the environment, then tell translate_tips() which column holds
#' # the current labels and which holds the replacements:
#' dir  <- system.file("extdata", "beetles_50p", package = "phylorug")
#' file <- system.file("extdata", "beetles_50p", "biogeo.csv", package = "phylorug")
#'
#' trees <- read_trees(dir)
#' dict  <- utils::read.csv(file)
#'
#' # Relabel tips from specimen codes to species names. The message reports how
#' # many tips were translated and how many were left unchanged in each tree:
#' translated <- translate_tips(trees, dict,
#'                              from_col = "specimen_code",
#'                              to_col   = "species_name")
translate_tips <- function(trees,
                           data,
                           from_col,
                           to_col,
                           verbose = TRUE) {
  if (missing(from_col) || missing(to_col)) {
    stop("`from_col` and `to_col` are required: name the columns in `data` ",
         "that hold the current and replacement labels.", call. = FALSE)
  }
  if (!inherits(trees, "list")) {
    stop("`trees` must be a list of `phylo` and/or `multiPhylo` objects, ",
         "one element per analysis (as returned by `read_trees()`). A pool of ",
         "tied-optimal trees (POY/TNT/PAUP*) is a `multiPhylo` element *inside* ",
         "this list, representing a single analysis.",
         call. = FALSE)
  }
  if (length(trees) == 0L) {
    stop("`trees` is empty; nothing to translate.", call. = FALSE)
  }
  # Each element must be a phylo or multiPhylo (pool)
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
  if (!inherits(data, "data.frame")) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!all(c(from_col, to_col) %in% colnames(data))) {
    stop("Columns \"", from_col, "\" and \"", to_col,
         "\" not found in `data`.", call. = FALSE)
  }
  if (anyDuplicated(data[[from_col]])) {
    stop("Values in `from_col` must be unique.", call. = FALSE)
  }
  if (anyNA(data[[from_col]])) {
    stop("`from_col` contains NA values. Remove or fill them.", call. = FALSE)
  }
  if (anyNA(data[[to_col]])) {
    stop("`to_col` contains NA values. Remove or fill them.", call. = FALSE)
  }
  # Preserve attributes that lapply would strip
  original_attrs <- attributes(trees)

  # Build translation dictionary
  trans_dict <- stats::setNames(
    as.character(data[[to_col]]),
    as.character(data[[from_col]])
  )

  # Translate tips in a single phylo
  translate_one <- function(tr) {
    match_idx <- match(tr$tip.label, names(trans_dict))
    matched   <- !is.na(match_idx)
    tr$tip.label[matched] <- trans_dict[match_idx[matched]]
    list(tree = tr, n_matched = sum(matched), n_unmatched = sum(!matched))
  }

  # Translate each analysis-- unwrap pools, translate every tree, rewrap
  trees <- lapply(names(trees), function(nm) {
    element <- trees[[nm]]

    if (inherits(element, "phylo")) {
      result <- translate_one(element)
      if (verbose) {
        message(nm, ": ", result$n_matched, " tips translated, ",
                result$n_unmatched, " unchanged")
      }
      return(result$tree)
    }

    # multiPhylo pool-- translate each tree inside
    n_trees  <- length(element)
    results  <- lapply(element, translate_one)
    element[] <- lapply(results, `[[`, "tree")

    if (verbose) {
      total_matched   <- sum(vapply(results, `[[`, integer(1), "n_matched"))
      total_unmatched <- sum(vapply(results, `[[`, integer(1), "n_unmatched"))
      message(nm, " (", n_trees, " trees): ",
              total_matched, " tips translated, ",
              total_unmatched, " unchanged (totals across pool)")
    }
    element
  })
  names(trees) <- original_attrs$names

  # Restore attributes lost by lapply (pool_sizes, class, etc.)
  for (a in names(original_attrs)) {
    if (is.null(attr(trees, a))) {
      attr(trees, a) <- original_attrs[[a]]
    }
  }

  trees
}
