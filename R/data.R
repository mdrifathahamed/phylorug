#' Sample beetle phylogenies (15-taxon subset)
#'
#' A named list (the same structure returned by [read_trees()]) of five
#' phylogenetic trees from a dung beetle phylogenomic study, Lopes et al.
#' (2024), pruned to 15 taxa for compact demonstration of the `phylorug`
#' workflow. The trees are already rooted, outgroup-removed, and tip labels
#' translated from museum codes to species names. Users can select a backbone
#' tree and comparison trees, then pass them directly to [check_taxa()],
#' [node_presence_matrix()], and [plot_phylorug()].
#'
#' The five trees represent independent phylogenomic analyses of the same
#' set of taxa using different inference methods and data types:
#' \describe{
#'   \item{70p_uce}{IQ-TREE maximum likelihood analysis of ultraconserved
#'     element (UCE) data. Recommended as the backbone tree.}
#'   \item{70p_partition_entropy}{IQ-TREE maximum likelihood analysis of
#'     partitioned sequence data.}
#'   \item{70p_ghost}{IQ-TREE GHOST heterotachous model on partitioned data.}
#'   \item{70p_ASTRAL_uce}{ASTRAL coalescent analysis of UCE gene trees.}
#'   \item{70p_ASTRAL_partition_entropy}{ASTRAL coalescent analysis of
#'     partitioned gene trees.}
#' }
#'
#' @format A named list of 5 objects of class `"phylo"`, each with 15 tips. Node
#'   labels contain support values: SH-aLRT/UFBoot2 for IQ-TREE trees and local
#'   posterior probability for ASTRAL trees.
#'
#' @source Lopes, F., Gunter, N., Gillett, C. P. D. T., et al. (2024).
#'   From museum drawer to tree: Historical DNA phylogenomics clarifies the
#'   systematics of rare dung beetles (Coleoptera: Scarabaeinae) from museum
#'   collections. \emph{PLOS ONE}, 19(12), e0309596.
#'   \doi{10.1371/journal.pone.0309596}
#'
#' @examples
#' # Select backbone and comparison trees
#' backbone <- sample_trees[["70p_uce"]]
#' others   <- sample_trees[names(sample_trees) != "70p_uce"]
#'
#' # Validate taxa
#' check_taxa(backbone, others)
#'
#' # Build the node presence matrix
#' npm <- node_presence_matrix(backbone, others, support_col = c(1, 2))
#'
#' # Plot in presence mode
#' plot_phylorug(backbone, npm)
"sample_trees"
