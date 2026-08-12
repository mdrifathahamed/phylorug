#' @keywords internal
#' @description
#' Generates rug plot visualizations directly on tree nodes to illustrate
#' clade recovery and support from multiple phylogenetic and phylogenomic
#' analyses. See `vignette("phylorug")` for a detailed guide.
#'
#' @section Core functions:
#' \describe{
#'   \item{\code{\link{read_trees}}}{Read tree files from disk}
#'   \item{\code{\link{translate_tips}}}{Rename specimen codes to species names}
#'   \item{\code{\link{check_taxa}}}{Diagnose taxon consistency across trees}
#'   \item{\code{\link{prune_to_shared}}}{Prune trees to shared taxa}
#'   \item{\code{\link{node_presence_matrix}}}{Make presence/support matrix }
#'   \item{\code{\link{plot_phylorug}}}{Draw the rug plot on a reference tree}
#' }
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
