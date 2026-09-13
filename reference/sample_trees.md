# Sample beetle phylogenies (15-taxon subset)

A named list (the same structure returned by
[`read_trees()`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md))
of five phylogenetic trees from a dung beetle phylogenomic study,
Montanaro, Lopes et al. (2026), pruned to 15 taxa for compact
demonstration of the `phylorug` workflow. The trees are already rooted,
outgroup-removed, and tip labels translated from museum codes to species
names. Users can select a backbone tree and comparison trees, then pass
them directly to
[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md),
[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
and
[`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md).

## Usage

``` r
sample_trees
```

## Format

A named list of 5 objects of class `"phylo"`, each with 15 tips. Node
labels contain support values: SH-aLRT/UFBoot2 for IQ-TREE trees and
local posterior probability for ASTRAL trees.

## Source

Montanaro, G., Lopes, F., Gunter, N.L., et al. (2026). Phylogenomics
resolves a 200-year-old puzzle: a revised tribal classification of
Afro-Eurasian dung beetles (Coleoptera: Scarabaeinae). *bioRxiv*.
[doi:10.64898/2026.07.22.740134](https://doi.org/10.64898/2026.07.22.740134)

## Details

The five trees represent independent phylogenomic analyses of the same
set of taxa using different inference methods and data types:

- 70p_uce:

  IQ-TREE maximum likelihood analysis of ultraconserved element (UCE)
  data. Recommended as the backbone tree.

- 70p_partition_entropy:

  IQ-TREE maximum likelihood analysis of partitioned sequence data.

- 70p_ghost:

  IQ-TREE GHOST heterotachous model on partitioned data.

- 70p_ASTRAL_uce:

  ASTRAL coalescent analysis of UCE gene trees.

- 70p_ASTRAL_partition_entropy:

  ASTRAL coalescent analysis of partitioned gene trees.

## Examples

``` r
# Core pipeline:
backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[names(sample_trees) != "70p_uce"]

# Optional: diagnose taxon overlap
check_taxa(backbone, others)
#> All 4 comparison trees share the same 15 taxa as the backbone.
#> [1] TRUE
#> attr(,"diagnostics")
#>                     comparison    status n_taxa missing extra
#> 1 70p_ASTRAL_partition_entropy identical     15              
#> 2               70p_ASTRAL_uce identical     15              
#> 3                    70p_ghost identical     15              
#> 4        70p_partition_entropy identical     15              

npm <- node_presence_matrix(backbone, others, support_col = c(1, 2))
tmp <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm, file = tmp)
unlink(tmp)
```
