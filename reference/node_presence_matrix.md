# Tabulate clade recovery and support across trees

Compares a set of phylogenetic trees against a reference topology, the
`"backbone"`. For each internal node of the backbone, the function asks
whether the same clade (the same set of tips) appears in each tree, and
records either its presence or its support value.

## Usage

``` r
node_presence_matrix(backbone, trees, support_col = 1, support_type = NULL)
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

- support_col:

  Integer or vector of integers (max 3). Specifies which value(s) to
  extract from multi-metric node labels. For example, if an IQ-TREE tree
  stores SH-aLRT and UFBoot as "80/95", passing `1` extracts only the
  first metric into a `support_1` matrix. Passing `c(1, 2)` efficiently
  extracts both simultaneously into `support_1` and `support_2`
  matrices, allowing you to easily switch between them during plotting
  without recalculating. Default is `1` .

- support_type:

  Optional named character vector mapping comparison trees to their
  support metrics (e.g. `"ufboot"`, `"sh_alrt"`, `"lpp"`,
  `"jackknife"`). Stored as an attribute of the returned list and used
  automatically by
  [`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
  when `mode = "support"`. If omitted,
  [`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
  will auto-normalize support values and apply universal thresholds.

## Value

A named list with one row per internal backbone node and one column per
comparison tree:

- presence:

  Clade presence: `1` where recovered, `0` where absent, or the
  proportion of pool trees recovering the clade.

- support_1, support_2, ...:

  One matrix per value in `support_col`. Raw support values where the
  clade was recovered,`NA` where absent. Named in the order requested,
  so `support_col = c(1, 2)` produces `support_1` and `support_2`.

## Details

This is the core function of the phylorug pipeline. It can be called
directly after
[`read_trees()`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md)
—
[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md)
is a useful diagnostic but not a prerequisite, as this function enforces
taxon matching internally. Every tree must include the complete set of
backbone taxa. Because an tree lacking a backbone taxon cannot assess
the presence of a clade containing that taxon,scoring such a clade as
'absent' would introduce a false negative.The function therefore
enforces strict taxon matching. Use
[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md)
to diagnose the discrepancies and
[`prune_to_shared()`](https://mdrifathahamed.github.io/phylorug/reference/prune_to_shared.md)
to harmonize the backbone and comparing tree taxa.

The function always computes both a presence matrix and a support matrix
in a single pass. The plotting function
[`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
decide which to use based on the visualisation context.

In presence matrix (`use_support = FALSE`),For multiple equally most
parsimonious trees (MPTs) the cell records the clade frequency among
MPTs and, for a single-tree this is 1 for presence or 0 for absence.

In support matrix `use_support = TRUE`, the cell records the support
value from the tree that recovered the clade, such as bootstrap or
posterior probability. For MPTs, support values are averaged across the
trees that recovered the clade; trees that did not recover it contribute
nothing. For a single tree this is the support value where the clade is
present, or `NA` where it is absent, since a clade that is not in the
tree has no node to carry a support value.

## Examples

``` r
# phylorug ships `sample_trees`: a named list of 5 phylo objects (beetles),
# the same structure you get from read_trees(). Just type `sample_trees`.
# Pick one analysis as the backbone, and the rest as comparison trees:
backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[names(sample_trees) != "70p_uce"]

# Optional: diagnose taxon overlap before building the matrix.
# This is not required — node_presence_matrix() enforces matching internally.
check_taxa(backbone, others)
#> All 4 comparison trees share the same 15 taxa as the backbone.
#> [1] TRUE
#> attr(,"diagnostics")
#>                     comparison    status n_taxa missing extra
#> 1 70p_ASTRAL_partition_entropy identical     15              
#> 2               70p_ASTRAL_uce identical     15              
#> 3                    70p_ghost identical     15              
#> 4        70p_partition_entropy identical     15              

# --- Presence/absence matrix -----------------------------------------------
# Just pass the backbone and comparison trees (support_col defaults to 1):
npmatrix <- node_presence_matrix(backbone, others)
npmatrix$presence     # clade presence/absence (1/0 or pool proportion)
#>    70p_ASTRAL_partition_entropy 70p_ASTRAL_uce 70p_ghost 70p_partition_entropy
#> 16                            1              1         1                     1
#> 17                            1              1         0                     1
#> 18                            1              1         0                     1
#> 19                            0              0         0                     1
#> 20                            0              0         0                     1
#> 21                            0              0         0                     0
#> 22                            0              0         0                     1
#> 23                            0              0         1                     1
#> 24                            1              1         1                     1
#> 25                            1              1         1                     1
#> 26                            1              1         1                     1
#> 27                            1              1         1                     1
#> 28                            1              1         1                     1
#> 29                            1              1         1                     1
npmatrix$support_1    # support values (from column 1 of the node labels)
#>    70p_ASTRAL_partition_entropy 70p_ASTRAL_uce 70p_ghost 70p_partition_entropy
#> 16                         0.00           0.01       100                   100
#> 17                         1.00           1.00        NA                   100
#> 18                         1.00           0.30        NA                   100
#> 19                           NA             NA        NA                   100
#> 20                           NA             NA        NA                   100
#> 21                           NA             NA        NA                    NA
#> 22                           NA             NA        NA                   100
#> 23                           NA             NA       100                   100
#> 24                         1.00           1.00       100                   100
#> 25                         1.00           1.00       100                   100
#> 26                         1.00           1.00       100                   100
#> 27                         1.00           1.00       100                   100
#> 28                         1.00           1.00       100                   100
#> 29                         0.91           1.00       100                   100

# --- Multiple support metrics at once --------------------------------------
# IQ-TREE stores SH-aLRT and UFBoot2 as "80/95". Pass their column positions
# to support_col to extract both in a single pass:
rugmt <- node_presence_matrix(backbone, others, support_col = c(1, 2))
rugmt$support_1       # first metric (SH-aLRT)
#>    70p_ASTRAL_partition_entropy 70p_ASTRAL_uce 70p_ghost 70p_partition_entropy
#> 16                         0.00           0.01       100                   100
#> 17                         1.00           1.00        NA                   100
#> 18                         1.00           0.30        NA                   100
#> 19                           NA             NA        NA                   100
#> 20                           NA             NA        NA                   100
#> 21                           NA             NA        NA                    NA
#> 22                           NA             NA        NA                   100
#> 23                           NA             NA       100                   100
#> 24                         1.00           1.00       100                   100
#> 25                         1.00           1.00       100                   100
#> 26                         1.00           1.00       100                   100
#> 27                         1.00           1.00       100                   100
#> 28                         1.00           1.00       100                   100
#> 29                         0.91           1.00       100                   100
rugmt$support_2       # second metric (UFBoot2)
#>    70p_ASTRAL_partition_entropy 70p_ASTRAL_uce 70p_ghost 70p_partition_entropy
#> 16                           NA             NA       100                   100
#> 17                           NA             NA        NA                   100
#> 18                           NA             NA        NA                   100
#> 19                           NA             NA        NA                   100
#> 20                           NA             NA        NA                   100
#> 21                           NA             NA        NA                    NA
#> 22                           NA             NA        NA                   100
#> 23                           NA             NA       100                   100
#> 24                           NA             NA       100                   100
#> 25                           NA             NA       100                   100
#> 26                           NA             NA       100                   100
#> 27                           NA             NA       100                   100
#> 28                           NA             NA       100                   100
#> 29                           NA             NA       100                   100

# --- With support_type (stored in npm, used by plot_phylorug) --------------
# Declaring support_type lets plot_phylorug() apply metric-specific
# thresholds automatically. If omitted, universal 0-1 thresholds are used.
support_type <- c(
  "70p_ASTRAL_partition_entropy" = "lpp",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ghost"                    = "ufboot",
  "70p_partition_entropy"        = "ufboot"
)
npm_typed <- node_presence_matrix(backbone, others,
                                  support_type = support_type)
attr(npm_typed, "support_type")   # stored as an attribute
#> 70p_ASTRAL_partition_entropy               70p_ASTRAL_uce 
#>                        "lpp"                        "lpp" 
#>                    70p_ghost        70p_partition_entropy 
#>                     "ufboot"                     "ufboot" 
```
