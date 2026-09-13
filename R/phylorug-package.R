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
#'   \item{\code{\link{node_presence_matrix}}}{Build the presence/support matrix}
#'   \item{\code{\link{plot_phylorug}}}{Draw the rug plot on a reference tree}
#'   \item{\code{\link{add_tree}}}{Add a tree to an existing matrix}
#' }
#'
#' @references
#' Wheeler, W.C. (1995). Sequence alignment, parameter sensitivity, and
#' the phylogenetic analysis of molecular data. \emph{Systematic Biology},
#' 44(3), 321-331. \doi{10.1093/sysbio/44.3.321}
#'
#' Giribet, G. (2003). Stability in phylogenetic formulations and its
#' relationship to nodal support. \emph{Systematic Biology}, 52(4),
#' 554-564. \doi{10.1080/10635150390223730}
#'
#' Sanders, J.G. (2010). Program note: Cladescan, a program for automated
#' phylogenetic sensitivity analysis. \emph{Cladistics}, 26(1), 114-116.
#' \doi{10.1111/j.1096-0031.2009.00280.x}
#'
#' Machado, D.J. (2015). YBYRA facilitates comparison of large
#' phylogenetic trees. \emph{BMC Bioinformatics}, 16, 204.
#' \doi{10.1186/s12859-015-0642-9}
#'
#' Wolfe, J.M. et al. (2019). A phylogenomic framework, evolutionary
#' timeline and genomic resources for comparative studies of decapod
#' crustaceans. \emph{Proceedings of the Royal Society B}, 286(1901),
#' 20190079. \doi{10.1098/rspb.2019.0079}
#'
#' Fu, Y. et al. (2025). Phylogenomic insights into the higher level
#' relationships within Culicomorpha (Diptera). \emph{Insect Systematics
#' and Diversity}, 9(6), ixaf056. \doi{10.1093/isd/ixaf056}
#'
#' Montanaro, G., Lopes, F. et al. (2026). Phylogenomics resolves a
#' 200-year-old puzzle: a revised tribal classification of Afro-Eurasian
#' dung beetles (Coleoptera: Scarabaeinae). \emph{bioRxiv}.
#' \doi{10.64898/2026.07.22.740134}
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
