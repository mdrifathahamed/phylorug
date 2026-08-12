# Trim the backbone and its comparison trees to their shared set of taxa

First finds the taxa that the backbone and all comparison trees share,
then trims everything else from the backbone and from every comparison
tree.

## Usage

``` r
prune_to_shared(backbone, trees, verbose = TRUE)
```

## Arguments

- backbone:

  A phylogenetic tree of class `"phylo"`.

- trees:

  A named list of `"phylo"` or `"multiPhylo"` objects, as returned by
  [`read_trees()`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md).

- verbose:

  Logical. If `TRUE` (default), reports how many taxa were dropped and
  names them.

## Value

A list with two elements: `backbone`, the pruned backbone tree, and
`trees`, the pruned comparison trees, retaining their names and their
original single-tree or pooled structure. The list carries a `"dropped"`
attribute naming the taxa that were removed.

## Details

A clade cannot be compared across trees that were not run on the same
taxon. Where a comparison tree is missing a backbone taxon, every
backbone clade containing that taxon becomes unevaluable in that tree,
and
[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
refuses to proceed. Trimming to the shared taxa resolves this by making
every tree consist the same taxon.

The cost is that the questions become narrower. Any clade containing a
dropped taxon no longer exists on the backbone and cannot be reported,
even where most comparison trees recovered it. Check the report from
[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md)
before trimming: if only one or two taxa are involved the loss is
minimal, but if several comparison trees each lack a different taxon,
the shared set can shrink quickly.

## See also

[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md)
to see which taxa are missing and from where, before deciding to prune.

## Examples

``` r
# Build a backbone and comparison trees that do NOT all share the same taxa,
# so there is something to trim. (In practice these come from read_trees().)
backbone <- ape::read.tree(
  text = "((((A,B),(C,D)),((E,F),(G,H))),(I,J));"
)
others <- list(
  # All ten taxa, same as the backbone:
  tree_1 = ape::read.tree(text = "((((A,B),(C,D)),((E,F),(G,H))),(I,J));"),
  # Missing J:
  tree_2 = ape::read.tree(text = "((((A,B),(C,D)),((E,F),(G,H))),I);"),
  # Missing I:
  tree_3 = ape::read.tree(text = "((((A,B),(C,D)),((E,F),(G,H))),J);")
)

# See which taxa are missing, and from which trees, before trimming:
check_taxa(backbone, others)
#> 2 comparison tree(s) are MISSING backbone taxa: tree_2, tree_3. Any backbone clade containing a missing taxon cannot be evaluated in those trees. Use `prune_to_shared()` to reduce the backbone and the comparison trees to their common taxa. See `attr(result, "diagnostics")` for the taxa involved.
#> [1] FALSE
#> attr(,"diagnostics")
#>   comparison  status n_taxa missing extra
#> 1     tree_2 missing      9       J      
#> 2     tree_3 missing      9       I      

# Trim the backbone and every comparison tree to their shared taxa:
shared <- prune_to_shared(backbone, others)
#> Dropped 2 of 10 taxa, leaving 8: I, J. Any clade containing a dropped taxon can no longer be reported.

# Which taxa were dropped:
attr(shared, "dropped")
#> [1] "I" "J"

# The pruned trees are ready for node_presence_matrix():
shared$backbone
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   A, B, C, D, E, F, ...
#> 
#> Rooted; no branch length.
shared$trees
#> $tree_1
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   A, B, C, D, E, F, ...
#> 
#> Rooted; no branch length.
#> 
#> $tree_2
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   A, B, C, D, E, F, ...
#> 
#> Rooted; no branch length.
#> 
#> $tree_3
#> 
#> Phylogenetic tree with 8 tips and 7 internal nodes.
#> 
#> Tip labels:
#>   A, B, C, D, E, F, ...
#> 
#> Rooted; no branch length.
#> 
```
