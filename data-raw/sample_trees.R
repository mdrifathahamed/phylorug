# data-raw/sample_trees.R
#
# Creates data/sample_trees.rda
#
# A 20-taxon subset of the Tarasov Lab 70% occupancy beetle phylogenomic
# dataset (Scarabaeinae + Aphodiinae), pruned from the full 289-taxon
# ingroup to serve as a quick demonstration of the phylorug pipeline.
#
# Source
#   Lopes et al. (2024) From museum drawer to tree: Historical DNA phylogenomics
#    clarifies the systematics of rare dung beetles (Coleoptera: Scarabaeinae)
#    from museum collections. PLOS ONE 19(12): e0309596.
#    https://doi.org/10.1371/journal.pone.0309596
#
# Five analyses of the same taxon set:
#   70p_uce                      IQ-TREE ML, UCE loci (backbone)
#   70p_partition_entropy        IQ-TREE ML, partitioned alignment
#   70p_ghost                    IQ-TREE GHOST, partitioned alignment
#   70p_ASTRAL_uce               ASTRAL III, UCE gene trees
#   70p_ASTRAL_partition_entropy ASTRAL III, partitioned gene trees
#
# The raw tree files and lookup table live outside the package and are
# not shipped to users.

library(ape)
library(readxl)

# -- 1. Read, root, drop outgroup, translate --------------------------------

trees <- phylorug::read_trees(
  "C:/Users/1/Desktop/phylorug_analysis/data/bettles/70p"
)

trees <- lapply(trees, function(tr) {
  tr <- root(tr, outgroup = c("NicorbUCE", "NicvesUCE"), resolve.root = TRUE)
  drop.tip(tr, c("NicorbUCE", "NicvesUCE"))
})

biogeo <- read_excel(
  "C:/Users/1/Desktop/phylorug_analysis/data/biogeo.xlsx",
  sheet = "BioGeo"
)
trees <- phylorug::translate_tips(trees, biogeo, from = "from", to = "to")

# -- 2. Pick 20 taxa -------------------------------------------------------
#
# Chosen by inspecting which backbone nodes are recovered by some but not
# all comparison trees, then selecting taxa that sit across a mix of
# unanimous, variable, and absent clades.

keep <- c(
  "Sisyphus_schaefferi_STL44",
  "Sisyphus_muricatus _STL5",
  "Nesosisyphus_pygmaeus _STL3",
  "Scarabaeus_westwoodi_STL10034",
  "Copris_fidius_ST005",
  "Catharsius_sp._STL10033",
  "Kheper_nigroaeneus_STL10036",
  "Circellium_bacchus _STL10270",
  "Helictopleurus_fissicollis _STL28",
  "Epactoides_hanski _STL29",
  "Nanos_dubitatus_STL5001",
  "Epilissus_cuprarius_STL5011",
  "Onthophagus_taurus_ST002",
  "Onthophagus_probus_STL10015",
  "Digitonthophagus_gazella _STL10213",
  "Frankenbergerius_armatus_STL10002",
  "Sarophorus_costatus_ST007",
  "Coptorhina_klugii _STL10039",
  "Gyronotus_pumilus_STL10003",
  "Bohepilissus_sp1 _STL10279"
)

stopifnot(length(keep) == 20L)
stopifnot(all(keep %in% trees[["70p_uce"]]$tip.label))

# -- 3. Prune and save ------------------------------------------------------

sample_trees <- lapply(trees, function(tr) keep.tip(tr, keep))

usethis::use_data(sample_trees, overwrite = TRUE, compress = "bzip2")
