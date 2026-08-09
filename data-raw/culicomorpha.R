# data-raw/culicomorpha.R
#
# Documents the origin of inst/extdata/culicomorpha/
#
# Ten phylogenomic trees from a whole-genome study of higher-level
# relationships within Culicomorpha (Diptera). 46 species from 8
# Culicomorpha families plus 3 outgroup taxa (2 Psychodidae:
# Phlebotomus_chinensis, Clogmia_albipunctata; 1 Scatopsidae:
# Coboldia_fuscipes). Analysed across three amino-acid matrices
# (Matrix1, Matrix1-kpi, Matrix1-smart) and one nucleotide matrix
# (Matrix2), using ten inference strategies: IQ-TREE partitioned
# site-homogeneous (LG/GTR), PMSF site-heterogeneous (C60),
# LG+C20+F+R, and ASTRAL coalescent.
#
# Source
#   Fu, Y., Du, S., Fang, X., Xu, Z. & Wang, X. (2025).
#   Phylogenomic insights into the higher-level relationships within
#   Culicomorpha (Diptera) revealed by whole-genome sequencing.
#   Insect Systematics and Diversity, 9(6), ixaf056.
#   https://doi.org/10.1093/isd/ixaf056
#
#   Data from Dryad: https://doi.org/10.5061/dryad.3bk3j9kvs
#
# Files copied unmodified into inst/extdata/culicomorpha/:
#   Matrix1_LG_C20_F_R.treefile
#   Matrix1-kpi_ASTRAL.treefile
#   Matrix1-kpi_partitioning.treefile
#   Matrix1-kpi_PMSF(H1.guide).treefile
#   Matrix1-kpi_PMSF(H2.guide).treefile
#   Matrix1-smart_ASTRAL.tre
#   Matrix1-smart_partitioning.treefile
#   Matrix1-smart_PMSF(H1.guide).treefile
#   Matrix1-smart_PMSF(H2.guide).treefile
#   Matrix2_partitioning.treefile
#
# No processing applied. Tree files are stored as deposited by the
# authors. Users root on the three outgroup taxa in the vignette:
#
#   trees <- lapply(trees, function(tr) {
#     tr <- ape::root(tr, outgroup = c("Coboldia_fuscipes",
#                                       "Phlebotomus_chinensis",
#                                       "Clogmia_albipunctata"),
#                     resolve.root = TRUE)
#     ape::drop.tip(tr, c("Coboldia_fuscipes",
#                          "Phlebotomus_chinensis",
#                          "Clogmia_albipunctata"))
#   })
