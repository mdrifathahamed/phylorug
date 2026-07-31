# data-raw/beetles_50p_and_70p.R
#
# This script does not create any R objects. It records where the raw
# tree files in inst/extdata/beetles_50p/ and inst/extdata/beetles_70p/
# came from.
#
# Both folders contain unmodified phylogenomic tree files from Lopes
# et al. (2024), a study of dung beetle systematics (Coleoptera:
# Scarabaeinae and Aphodiinae). The "50p" and "70p" refer to the
# minimum taxon-occupancy threshold used to filter the sequence
# alignment before tree inference.
#
# Each folder holds five trees from independent analyses:
#   *_uce.tre                      IQ-TREE ML on UCE loci
#   *_partition_entropy.tre        IQ-TREE ML on partitioned alignment
#   *_ghost.tre                    IQ-TREE GHOST heterotachous model
#   *_ASTRAL_uce.tre               ASTRAL coalescent on UCE gene trees
#   *_ASTRAL_partition_entropy.tre ASTRAL coalescent on partitioned genes
#
# The file biogeo.xlsx (in beetles_50p/) is a lookup table that maps
# short specimen codes to full species names. It serves both the 50p
# and 70p datasets.
#
# All trees include two outgroup taxa (NicorbUCE, NicvesUCE) that
# users must root on and remove before analysis. Tip labels are
# specimen codes that users translate with translate_tips() and the
# biogeo.xlsx lookup table.
#
# Reference
#   Lopes et al. (2024) From museum drawer to tree: Historical DNA phylogenomics
#   clarifies the systematics of rare dung beetles (Coleoptera: Scarabaeinae)
#   from museum collections. PLOS ONE 19(12): e0309596.
#   https://doi.org/10.1371/journal.pone.0309596
