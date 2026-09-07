# data-raw/beetles_50p_and_70p.R
#
# Prunes the full beetle trees to the 43 taxa of Onthophagina sensu
# lato for use in the phylorug package vignettes and paper figures.
#
# The full 316-taxon tree from Montanaro, Lopes et al. (2026) is too
# large for a package vignette, so we prune it to a smaller subset
# for demonstration purposes. We chose the 43 taxa of Onthophagina
# sensu lato so that the subset contains a complete subtribe rather
# than an arbitrary selection of tips.
#
# The "50p" and "70p" refer to the minimum taxon-occupancy threshold
# used to filter the UCE sequence alignment before tree inference.
# The same 43 Onthophagina ingroup taxa are present in both datasets
# and across all 5 trees within each folder, plus 2 outgroup taxa
# (NicorbUCE, NicvesUCE), giving 45 tips total per tree.
#
# Each folder holds five trees from independent analyses:
#   *_uce.tre                      IQ-TREE ML on UCE loci
#   *_partition_entropy.tre        IQ-TREE ML on partitioned alignment
#   *_ghost.tre                    IQ-TREE GHOST heterotachous model
#   *_ASTRAL_uce.tre               ASTRAL coalescent on UCE gene trees
#   *_ASTRAL_partition_entropy.tre ASTRAL coalescent on partitioned genes
#
# The file biogeo.csv (in beetles_50p/) maps all 316 specimen codes
# to species names and serves both datasets. Users translate tip
# labels with translate_tips() using this lookup table.
#
# All trees include two outgroup taxa (NicorbUCE, NicvesUCE) that
# users must root on and remove before analysis.
#
# Support values, topology, and branch lengths for the retained taxa
# are unmodified from the original inference.
#
# A prototype version of phylorug was used in the source study itself
# to generate the rug plots in Montanaro, Lopes et al. (2026).
#
# Source
#   Montanaro, G., Lopes, F., Gunter, N.L., Scholtz, C., Davis, A.L.,
#   Losacco, F., Rossini, M., Gillett, C.P.D.T., Saxton, N.A.,
#   Stone, R.L., Daniel, G.M., and Tarasov, S. (2026). Phylogenomics
#   resolves a 200-year-old puzzle: a revised tribal classification of
#   Afro-Eurasian dung beetles (Coleoptera: Scarabaeinae). bioRxiv.
#   https://doi.org/10.64898/2026.07.22.740134
#
# HOW TO REPRODUCE
#   1. Edit the paths below to point to the full unpublished trees
#      on your machine and the package output directories.
#   2. Run the script. It will identify 43 Onthophagina taxa from
#      biogeo.csv, verify they exist in both tree sets, prune all
#      trees, and plot the results for inspection.
#   3. If the plots look good, uncomment the write_pruned() calls
#      at the bottom and run again to write files into the package.

library(ape)

# --- EDIT THESE PATHS -------------------------------------------------------
input_dir_50p  <- "C:/Users/1/Desktop/phylorug_analysis/data/bettles/50p"
input_dir_70p  <- "C:/Users/1/Desktop/phylorug_analysis/data/bettles/70p"
output_dir_50p <- "C:/Users/1/Desktop/phylorug/inst/extdata/beetles_50p"
output_dir_70p <- "C:/Users/1/Desktop/phylorug/inst/extdata/beetles_70p"
biogeo_path    <- "C:/Users/1/Desktop/phylorug/inst/extdata/beetles_50p/biogeo.csv"
# ----------------------------------------------------------------------------

outgroup <- c("NicorbUCE", "NicvesUCE")

# --- STEP 1: Identify Onthophagina taxa from biogeo.csv ---------------------

biogeo <- read.csv(biogeo_path)

onthophagina_genera <- paste0(
  "Kurtops|Hamonthophagus|Onthophagus|Digitonthophagus|",
  "Phalops|Amietina|Caccobius|Eusaproceius|Proagoderus|",
  "Parascatonomus|Macronthophagus|Cleptocaccobius|Euonthophagus"
)

onthophagina_rows <- biogeo[grep(onthophagina_genera,
                                 biogeo$species_name,
                                 ignore.case = TRUE), ]

keep_codes <- onthophagina_rows$specimen_code
cat("Onthophagina taxa identified:", length(keep_codes), "\n")

# --- STEP 2: Verify codes exist in both tree sets --------------------------

files_70p    <- list.files(input_dir_70p, pattern = "\\.tre$", full.names = TRUE)
ref_tree_70p <- read.tree(files_70p[1])
files_50p    <- list.files(input_dir_50p, pattern = "\\.tre$", full.names = TRUE)
ref_tree_50p <- read.tree(files_50p[1])

keep_70p <- intersect(keep_codes, ref_tree_70p$tip.label)
keep_50p <- intersect(keep_codes, ref_tree_50p$tip.label)

cat("70p:", length(keep_70p), "found,",
    length(setdiff(keep_codes, ref_tree_70p$tip.label)), "missing\n")
cat("50p:", length(keep_50p), "found,",
    length(setdiff(keep_codes, ref_tree_50p$tip.label)), "missing\n")

# --- STEP 3: Prune all trees to Onthophagina + outgroup --------------------

prune_to_clade <- function(input_dir, keep_tips, outgroup) {
  files <- list.files(input_dir, pattern = "\\.tre$", full.names = TRUE)
  trees <- lapply(files, read.tree)
  names(trees) <- basename(files)

  available <- intersect(keep_tips, trees[[1]]$tip.label)
  keep_all  <- c(available, outgroup)

  pruned <- lapply(trees, function(tr) keep.tip(tr, keep_all))

  pruned_sets <- lapply(pruned, function(tr) sort(tr$tip.label))
  stopifnot(length(unique(pruned_sets)) == 1L)
  cat("Pruned to", Ntip(pruned[[1]]), "tips (",
      length(available), "ingroup +", length(outgroup), "outgroup)\n")

  list(trees = pruned, filenames = basename(files))
}

result_70p <- prune_to_clade(input_dir_70p, keep_70p, outgroup)
result_50p <- prune_to_clade(input_dir_50p, keep_50p, outgroup)

# --- STEP 4: Visual inspection ----------------------------------------------

par(mfrow = c(1, 2))
plot(result_70p$trees[[1]], cex = 0.5,
     main = paste0("70p Onthophagina (", Ntip(result_70p$trees[[1]]), " tips)"))
plot(result_50p$trees[[1]], cex = 0.5,
     main = paste0("50p Onthophagina (", Ntip(result_50p$trees[[1]]), " tips)"))
par(mfrow = c(1, 1))

# --- STEP 5: Write files (uncomment when satisfied) -------------------------

write_pruned <- function(result, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  mapply(
    function(tr, fname) write.tree(tr, file = file.path(output_dir, fname)),
    result$trees, result$filenames
  )
  invisible(NULL)
}

# write_pruned(result_70p, output_dir_70p)
# write_pruned(result_50p, output_dir_50p)
# cat("\nFiles written to:\n", output_dir_70p, "\n", output_dir_50p, "\n")
