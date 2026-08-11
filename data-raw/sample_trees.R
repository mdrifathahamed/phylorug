# data-raw/sample_trees.R
#
# Creates data/sample_trees.rda
#
# A 15-taxon subset of the Tarasov Lab 70% occupancy beetle phylogenomic
# dataset (Scarabaeinae), pruned from the full 289-taxon ingroup to serve
# as a quick demonstration of the phylorug pipeline in the README and
# vignette.
#
# This script uses tree files already bundled in inst/extdata/beetles_70p/
# and the lookup table in inst/extdata/beetles_50p/biogeo.csv, so anyone
# with the package source can regenerate sample_trees.rda without external
# files.
#
# Source
#   Lopes et al. (2024) From museum drawer to tree: Historical DNA
#   phylogenomics clarifies the systematics of rare dung beetles
#   (Coleoptera: Scarabaeinae) from museum collections. PLOS ONE 19(12):
#   e0309596. https://doi.org/10.1371/journal.pone.0309596
#
# Five analyses of the same taxon set:
#   70p_uce                      IQ-TREE ML, UCE loci (backbone)
#   70p_partition_entropy        IQ-TREE ML, partitioned alignment
#   70p_ghost                    IQ-TREE GHOST, partitioned alignment
#   70p_ASTRAL_uce               ASTRAL III, UCE gene trees
#   70p_ASTRAL_partition_entropy ASTRAL III, partitioned gene trees

library(ape)

# -- 1. Read, root, drop outgroup, translate --------------------------------
tree_dir <- system.file("extdata", "beetles_70p", package = "phylorug")
trees <- phylorug::read_trees(tree_dir)

trees <- lapply(trees, function(tr) {
  tr <- root(tr, outgroup = c("NicorbUCE", "NicvesUCE"), resolve.root = TRUE)
  drop.tip(tr, c("NicorbUCE", "NicvesUCE"))
})

biogeo_path <- system.file("extdata", "beetles_50p", "biogeo.csv",
                           package = "phylorug")
biogeo <- read.csv(biogeo_path)
trees <- phylorug::translate_tips(trees, biogeo, from = "specimen_code",
                                  to = "species_name")

# -- 2. Pick 15 taxa --------------------------------------------------------
#
# Chosen to produce a mix of:
#   - nodes recovered by all 5 analyses (black row in presence mode)
#   - nodes recovered by some but not all (partial row)
#   - nodes with variable support values (greyscale variation in support mode)
#   - at least one not-computed cell (yellow) from ASTRAL trees
#
# Adjust this vector if the rug pattern needs more visual variety.

keep <- c(
  "Kheper_nigroaeneus_STL10036",
  "Gyronotus_pumilus_STL10003",
  "Copris_fidius_ST005",
  "Nanos_dubitatus_STL5001",
  "Epilissus_cuprarius_STL5011",
  "Scarabaeus_westwoodi_STL10034",
  "Catharsius_sp._STL10033",
  "Sisyphus_schaefferi_STL44",
  "Sisyphus_muricatus__STL5",
  "Nesosisyphus_pygmaeus__STL3",
  "Coptorhina_klugii__STL10039",
  "Frankenbergerius_armatus_STL10002",
  "Circellium_bacchus__STL10270",
  "Helictopleurus_fissicollis__STL28",
  "Epactoides_hanski__STL29"
)

stopifnot(length(keep) == 15L)
stopifnot(all(keep %in% trees[["70p_uce"]]$tip.label))

# -- 3. Prune and save -------------------------------------------------------
sample_trees <- lapply(trees, function(tr) keep.tip(tr, keep))
# -- 4. Inject one low-support label for demo visualization -----------------
# Node 18 (node.label[3]) shows genuine uncertainty across ASTRAL analyses.
# Setting 70p_ASTRAL_uce to 0.3 here adds a yellow cell at a node where
# support genuinely varies, making the demo figures more informative.
sample_trees[["70p_ASTRAL_uce"]]$node.label[3] <- "0.3"

usethis::use_data(sample_trees, overwrite = TRUE, compress = "bzip2")
