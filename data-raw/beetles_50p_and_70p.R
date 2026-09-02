# data-raw/beetles_50p_and_70p.R
#
# Prunes the full (unpublished) beetle trees down to a package-friendly
# taxon subset, and records where the data came from.
#
# Both folders originally contained the full, unpublished phylogenomic
# tree files from Montanaro, Lopes et al. (2026) (see Source below), a
# study of dung beetle systematics (Coleoptera: Scarabaeinae). The
# "50p" and "70p" refer to the minimum taxon-occupancy threshold used to
# filter the sequence alignment before tree inference.
#
# Because the full trees are not yet published, a random subset of taxa
# is pruned before shipping them with the package: 50 taxa from the 50p
# dataset, 70 taxa from the 70p dataset (plus the 2 outgroup taxa in
# each, so 52 and 72 tips respectively). The exact same taxon set is
# kept across all 5 trees within a folder, which node_presence_matrix()
# requires. Support values, topology, and branch lengths for the
# retained taxa are otherwise unmodified from the original inference.
#
# Each folder holds five trees from independent analyses:
#   *_uce.tre                      IQ-TREE ML on UCE loci
#   *_partition_entropy.tre        IQ-TREE ML on partitioned alignment
#   *_ghost.tre                    IQ-TREE GHOST heterotachous model
#   *_ASTRAL_uce.tre               ASTRAL coalescent on UCE gene trees
#   *_ASTRAL_partition_entropy.tre ASTRAL coalescent on partitioned genes
#
# The file biogeo.csv (in beetles_50p/) is a small demonstration lookup
# table that maps short specimen codes to full species names, used by
# the translate_tips() examples. It serves both the 50p and 70p datasets.
#
# All trees include two outgroup taxa (NicorbUCE, NicvesUCE) that
# users must root on and remove before analysis. Tip labels are
# specimen codes that users translate with translate_tips() and the
# biogeo.csv lookup table.
#
# Source
#   Montanaro, G., Lopes, F., Gunter, N.L., Scholtz, C., Davis, A.L.,
#   Losacco, F., Rossini, M., Gillett, C.P.D.T., Saxton, N.A., Stone, R.L.,
#   Daniel, G.M., and Tarasov, S. (2026). Phylogenomics resolves a
#   200-year-old puzzle: a revised tribal classification of Afro-Eurasian
#   dung beetles (Coleoptera: Scarabaeinae). bioRxiv.
#   https://doi.org/10.64898/2026.07.22.740134
#
# HOW TO REPRODUCE THE PRUNING
#   1. Edit the four paths below so they point to the right folders on
#      your machine.
#   2. Run the script down to STEP 3. Two tree plots will appear --
#      check they look reasonable (not a flat comb, all tips distinct).
#   3. If they look good, uncomment the two write_pruned(...) lines at
#      the bottom and run the script again -- that's what actually
#      writes the files into the package.
#   4. Not happy with the tree shape? Change `seed` (in the two
#      prune_dataset() calls below) to a different number and repeat
#      from step 2 -- that draws a different random set of taxa.

library(ape)

# --- STEP 1: EDIT THESE FOUR PATHS to match your machine -------------------
input_dir_50p  <- "C:/Users/1/Desktop/phylorug_analysis/data/bettles/50p"
input_dir_70p  <- "C:/Users/1/Desktop/phylorug_analysis/data/bettles/70p"
output_dir_50p <- "C:/Users/1/Desktop/phylorug/inst/extdata/beetles_50p"
output_dir_70p <- "C:/Users/1/Desktop/phylorug/inst/extdata/beetles_70p"
# -----------------------------------------------------------------------------

# These two taxa are outgroups in every tree. They are always kept, and
# don't count toward the 50 / 70 target.
outgroup <- c("NicorbUCE", "NicvesUCE")

# Reads every .tre file in `input_dir`, picks ONE random set of `n_target`
# ingroup taxa (using `seed` so the draw is reproducible), and prunes every
# tree down to that same set + the outgroup. All trees in a dataset must
# end up with identical taxa, because node_presence_matrix() requires it.
prune_dataset <- function(input_dir, n_target, seed = 1) {

  files <- list.files(input_dir, pattern = "\\.tre$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No .tre files found in ", input_dir)
  }

  trees <- lapply(files, read.tree)
  names(trees) <- basename(files)

  # Check all trees already share the same taxa, before we prune anything.
  tip_sets <- lapply(trees, function(tr) sort(tr$tip.label))
  if (length(unique(tip_sets)) != 1L) {
    stop("Trees in ", input_dir, " do not currently share identical taxa -- ",
         "fix that before pruning.")
  }

  ingroup <- setdiff(trees[[1]]$tip.label, outgroup)

  set.seed(seed)
  keep <- c(sample(ingroup, n_target), outgroup)

  pruned <- lapply(trees, function(tr) keep.tip(tr, keep))

  # Check every pruned tree ended up with the same taxa and the right size.
  pruned_sets <- lapply(pruned, function(tr) sort(tr$tip.label))
  stopifnot(length(unique(pruned_sets)) == 1L)
  stopifnot(ape::Ntip(pruned[[1]]) == n_target + length(outgroup))

  list(trees = pruned, filenames = basename(files))
}

# Writes a prune_dataset() result out to `output_dir`, one file per tree,
# reusing the original filenames.
write_pruned <- function(result, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  mapply(
    function(tr, fname) write.tree(tr, file = file.path(output_dir, fname)),
    result$trees, result$filenames
  )
  invisible(NULL)
}

# --- STEP 2: run the pruning ---------------------------------------------
result_50p <- prune_dataset(input_dir_50p, n_target = 50, seed = 1)
result_70p <- prune_dataset(input_dir_70p, n_target = 70, seed = 1)

# --- STEP 3: look at the trees before writing anything ---------------------
plot(result_50p$trees[[1]], cex = 0.6, main = "50p (pruned, 52 tips)")
plot(result_70p$trees[[1]], cex = 0.6, main = "70p (pruned, 72 tips)")

# --- STEP 4:  write the files -------
write_pruned(result_50p, output_dir_50p)
write_pruned(result_70p, output_dir_70p)
