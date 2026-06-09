#' Check tip label consistency across a set of phylogenetic trees
#'
#' Compares a list of phylogenetic trees against a reference tree
#' (the first tree in the list) to ensure that all trees contain the exact same
#' set of taxa (tip labels). If discrepancies are found, detailed diagnostic
#' reports listing missing and extra taxa are provided.
#'
#' @param trees A list containing multiple \code{phylo} objects, or a
#'   uncompressed multiPhylo object containing the trees to be compared.
#'
#' @param verbose Logical. If \code{TRUE} (default), detailed status messages
#'   and diagnostic reports for mismatched trees are printed to the console.
#'
#'
#' @return A single logical value. \code{TRUE} if all trees contain
#'   identical taxa. \code{FALSE} if any tree differs, in which case
#'   the pipeline should not proceed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Assuming 'my_trees' is a named list of phylo objects loaded into R:
#' trees <- read_trees("path/to/your/trees")
#' result <- check_taxa(trees)
#'
#' # Check if it's safe to proceed:
#' if (result) {
#'   message("Safe to proceed.")
#' }
#' }
check_taxa <- function(trees, verbose = TRUE) {
  # Validate input

  if (!inherits(trees, c("list", "multiPhylo"))) {
    stop(
      "`trees` must be a list of phylo objects or a multiPhylo object.",
      call. = FALSE
    )
  }
  if (inherits(trees, "multiPhylo") && !is.null(attr(trees, "TipLabel"))) {
    stop(
      "`trees` appears to be a compressed multiPhylo object ",
      "(tip labels stored as an attribute, not per-tree). ",
      "Decompress it first with `ape::uncompressTipLabel(trees)`.",
      call. = FALSE
    )
  }
  if (length(trees) == 0) {
    stop(
      "`trees` is empty. Supply at least two trees.",
      call. = FALSE
    )
  }

  if (any(vapply(trees, is.null, logical(1)))) {
    stop(
      "`trees` contains NULL elements. ",
      "Ensure all tree files or objects were loaded successfully.",
      call. = FALSE
    )
  }

  # Extract and sort tip labels from each tree
  taxa_list <- lapply(trees, function(tr) sort(tr$tip.label))

  # Use first tree as reference
  ref <- taxa_list[[1]]

  # Compare all trees to reference
  same <- vapply(taxa_list, function(x) identical(x, ref), logical(1))

  # Report results if verbose
  if (verbose) {
    if (all(same)) {
      message("\u2705 All trees contain the same set of taxa.")
    } else {
      bad <- which(!same)

      # 1. Quietly compile data logs for programmatic backup
      full_diagnostics <- lapply(bad, function(i) {
        list(
          tree_index = i,
          missing    = setdiff(ref, taxa_list[[i]]),
          extra      = setdiff(taxa_list[[i]], ref)
        )
      })
      names(full_diagnostics) <- paste0("tree_", bad)
      options(phylorug_diagnostics = full_diagnostics)

      # 2. Strict, elegant, ONE-line status report
      message(
        "\u274c Taxa mismatch detected in ", length(bad), " tree(s). ",
        "Run `getOption('phylorug_diagnostics')` to view the full, detailed report."
      )
    }
  }

  # Return single logical
  all(same)
}
