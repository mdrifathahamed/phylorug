# Diagnose taxon consistency between backbone and comparison trees

For each comparison tree, checks whether its tip labels (taxon) match
the backbone's exactly, contain extra labels, or are missing some. A
tree missing a tip label of the backbone cannot be scored for any clade
containing that tip. We recommend running this diagnostic before
building the matrix, especially when combining trees from different
studies or pipelines. However, it is not a mandatory step
[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
enforces taxon matching internally and will error with a clear message
if any comparison tree is missing backbone taxa. If you are certain all
trees share the same taxon set, you can proceed directly to
[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md).

## Usage

``` r
check_taxa(backbone, trees, verbose = TRUE)
```

## Arguments

- backbone:

  The reference tree, as a `"phylo"` object. The rug is drawn on this
  tree, and its clades are searched for in each comparison tree.

- trees:

  A named list of `"phylo"` or `"multiPhylo"` objects, one per
  comparison, as returned by
  [`read_trees()`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md).If
  the backbone is present in the list it is removed automatically before
  comparison, so you can pass the full list from
  [`read_trees()`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md)
  without subsetting first.

- verbose:

  Logical. If `TRUE` (default), reports the outcome and names any
  mismatched comparison trees.

## Value

A single logical value: `TRUE` if every comparison tree shares the
backbone's taxa exactly. The value carries a `"diagnostics"` attribute,
a data frame with one row per comparison tree giving its status
(`"identical"`, `"superset"` or `"missing"`) and the taxa missing from,
or extra to, the backbone. Diagnostics are attached whatever `verbose`
is set to.

## Details

This function reports; it does not modify the trees. If the backbone is
present in comparison `trees` it is silently removed before comparison.
Three outcomes are

- identical:

  All comparison trees share the backbone's taxa exactly.

- superset:

  One or more comparison trees contain every backbone taxon, plus some
  extra. Every backbone clade can still be evaluated, so no action is
  needed; the extra taxa can be ignored.

- missing:

  One or more comparison trees lack one or more backbone taxa. Any
  backbone clade containing a missing taxon might produce false absence
  in that tree. Scoring such a clade as absent would report a rejection
  where no question was ever put. Use
  [`prune_to_shared()`](https://mdrifathahamed.github.io/phylorug/reference/prune_to_shared.md)
  to reduce the backbone and the comparison trees to their common taxa.

Where a comparison tree is a pool of several equally optimal trees,
every tree in the pool is checked. A pool whose trees disagree about
their own taxa is a data problem, and is an error.

## See also

[`read_trees()`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md)
to load the trees and
[`prune_to_shared()`](https://mdrifathahamed.github.io/phylorug/reference/prune_to_shared.md)
to reduce them to a common taxon set.

## Examples

``` r
# sample_trees is a named list of 5 phylo objects shipped with phylorug,
# the same structure you get from read_trees(). Pick one analysis as the
# backbone, and the rest as comparison trees:
backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[names(sample_trees) != "70p_uce"]

# Diagnose taxon consistency (verbose = TRUE by default reports the outcome):
result <- check_taxa(backbone, others)
#> All 4 comparison trees share the same 15 taxa as the backbone.

# If any tree is missing or has extra taxa, inspect the full report:
attr(result, "diagnostics")
#>                     comparison    status n_taxa missing extra
#> 1 70p_ASTRAL_partition_entropy identical     15              
#> 2               70p_ASTRAL_uce identical     15              
#> 3                    70p_ghost identical     15              
#> 4        70p_partition_entropy identical     15              
```
