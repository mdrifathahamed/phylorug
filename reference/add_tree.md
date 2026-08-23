# Add a comparison tree to an existing node presence matrix

Adds a new comparison tree to an existing node presence matrix without
rerunning whole pipeline from scratch. The function computes clade
recovery(presence) and support for the new tree only and attaches the
results as a new column in every matrix of the npm. All existing columns
remain unchanged.

## Usage

``` r
add_tree(
  npm,
  backbone,
  new_tree,
  name,
  support_col = 1,
  support_type = NULL,
  pool_threshold = 1
)
```

## Arguments

- npm:

  A node presence matrix list, as returned by
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
  or a previous call to `add_tree()`.

- backbone:

  The backbone tree used to create `npm`. Must be the same `phylo`
  object (same topology and tip labels) that was passed to
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md).

- new_tree:

  A `phylo` or `multiPhylo` object representing the new comparison
  analysis to add. Must be rooted.

- name:

  Character string. Column name for the new tree in the npm matrices.
  Required. Must be a single non-empty string and must not duplicate an
  existing column name.

- support_col:

  Integer or vector of integers (max 3). Which support value(s) to
  extract from the new tree's node labels, matching the `support_col`
  used when the original `npm` was built. Default is `1`. For example,
  if the new tree is from IQ-TREE with compound labels `"80/95"`,
  passing `c(1, 2)` extracts both metrics.

- support_type:

  Optional character string naming the support metric of the new tree
  (e.g. `"ufboot"`, `"lpp"`, `"jackknife"`). If supplied, it is appended
  to the npm's `"support_type"` attribute so that
  [`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
  can apply metric-specific thresholds automatically. If `NULL`, the new
  tree's support values will be auto-normalized (values \>1 divided
  by 100) and binned against universal thresholds at plot time. For this
  reason, raw Bremer support (an unbounded integer) is not supported;
  use the`"bremer_ratio"` instead, which is bounded between 0 and 1.

- pool_threshold:

  Numeric between 0 and 1. Controls how pools are scored. Default `1.0`
  (strict consensus). Should match the value used when the original npm
  was built. See
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
  for details. This parameter is only applicable when evaluating
  parsimony-based `multiPhylo` objects (e.g., Most Parsimonious Trees
  generated via TNT or similar software).

## Value

The updated npm list with the new tree appended as an additional column
in the `presence` matrix and each `support_*` matrix. Attributes
`"node_id"` and `"pool_sizes"` are updated, and `"support_type"` is
extended if the argument was supplied.

## Details

Three taxon scenarios are handled automatically:

- Identical:

  The new tree shares all backbone taxa exactly. Every clade is
  evaluated normally.

- Superset:

  The new tree contains all backbone taxa plus extras. The extra taxa
  are pruned internally so that clade keys match the backbone's. Without
  this pruning, a backbone clade `(A, B)` would fail to match a new-tree
  node grouping `(A, B, X)`, producing a false absence. A message
  reports the pruning. The user's original tree object is not modified.

- Missing:

  The new tree lacks one or more backbone taxa. Any backbone clade
  containing a missing taxon cannot be evaluated in this tree and is set
  to `NA` in both the presence and support matrices. These appear as red
  ("not computed") cells in the rug plot. A message lists the missing
  taxa and suggests rebuilding the npm from scratch with
  [`prune_to_shared()`](https://mdrifathahamed.github.io/phylorug/reference/prune_to_shared.md)
  if the red cells are undesirable. Unlike
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md),
  which errors on missing taxa, `add_tree()` is deliberately permissive:
  when adding a tree after the initial analysis, different taxon
  sampling is a legitimate use case, not a data error.

The new tree must be rooted. If it is a `multiPhylo` object (a pool of
equally optimal trees from one search), clade presence is determined by
`pool_threshold`: the default (`1.0`, strict consensus) requires the
clade in every pool tree; `0.5` applies majority rule; `0` records the
raw proportion. Support extraction is skipped for pools because
averaging node labels across equally optimal trees is not scientifically
valid; support cells for that column will be `NA`. To include support
values from a parsimony analysis, compute the consensus tree with
bootstrap or jackknife support upstream (e.g. in TNT) and pass that
consensus as a single `phylo` object.

If `support_type` is supplied, it is appended to the `"support_type"`
attribute stored in the npm object, keyed by the new tree's name. This
allows
[`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
to apply the correct metric-specific thresholds without the user
re-declaring the full vector.

## See also

[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
to build the initial npm,
[`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
to visualise the result, and
[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md)
to diagnose taxon overlap before adding.

## Examples

``` r
# Build an npm with 3 of the 5 shipped trees:
backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[c("70p_partition_entropy", "70p_ghost")]
npm      <- node_presence_matrix(backbone, others)
ncol(npm$presence)   # 2 comparison trees
#> [1] 2

# Add a fourth tree without re-running the pipeline:
npm <- add_tree(npm, backbone,
                sample_trees[["70p_ASTRAL_uce"]],
                name = "ASTRAL_uce")
ncol(npm$presence)   # now 3
#> [1] 3

# Add the fifth tree with a declared support type:
npm <- add_tree(npm, backbone,
                sample_trees[["70p_ASTRAL_partition_entropy"]],
                name = "ASTRAL_part_entropy",
                support_type = "lpp")
ncol(npm$presence)   # now 4
#> [1] 4
attr(npm, "support_type")
#> ASTRAL_part_entropy 
#>               "lpp" 

# To add a tree from an external file (not run):
#   new_tree <- ape::read.tree("path/to/new_analysis.treefile")
#   new_tree <- ape::root(new_tree, outgroup = "outgroup_sp",
#                         resolve.root = TRUE)
#   new_tree <- ape::drop.tip(new_tree, "outgroup_sp")
#   npm <- add_tree(npm, backbone, new_tree, name = "new_analysis",
#                   support_type = "ufboot")
```
