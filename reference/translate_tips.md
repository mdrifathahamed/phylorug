# Relabel tips across a list of trees using a lookup table

Translates tip labels of each tree in a list using a lookup table
supplied as a data frame. This function is optional in the phylorug
workflow, use it only if your tip labels are specimen codes, accession
numbers, or any other identifiers that need converting to a different
format, such as scientific species names.

If your tip labels are already in the correct format, skip this step and
proceed directly to
[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md).
Only tips with a matching entry in `from_col` are replaced with the
corresponding value in `to_col`. Unmatched tips are left unchanged. This
function modifies only the `tip.label` field of each tree. Tree
topology, branch lengths, and node labels are unaffected.

## Usage

``` r
translate_tips(trees, data, from_col, to_col, verbose = TRUE)
```

## Arguments

- trees:

  A named list of `"phylo"` and/or `"multiPhylo"` objects, one element
  per analysis, as returned by
  [`read_trees()`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md).
  A `"multiPhylo"` element is a pool of tied-optimal trees
  (POY/TNT/PAUP\*) treated as a single analysis.

- data:

  A data frame containing the label translation lookup table. Can be a
  standard `data.frame` or a `tibble`.

- from_col:

  A character string naming the column in `data` that holds the current
  tip labels (for example `"specimen_code"`). Required.

- to_col:

  A character string naming the column in `data` that holds the
  replacement labels (for example `"scientific_name"`). Required.

- verbose:

  Logical. If `TRUE` (default), reports how many tip labels were
  translated and how many were left unchanged for each tree.

## Value

A named list of `"phylo"` objects identical in structure to `trees`,
with matching tip labels replaced according to `data`. Tree topology and
edge lengths remain unchanged.

## Examples

``` r
# Example trees and a lookup table ship inside phylorug, so we first read
# them into the environment, then tell translate_tips() which column holds
# the current labels and which holds the replacements:
dir  <- system.file("extdata", "beetles_50p", package = "phylorug")
file <- system.file("extdata", "beetles_50p", "biogeo.csv",
                    package = "phylorug")

trees <- read_trees(dir)
#> Read 5 analyses (5 trees) from: /home/runner/work/_temp/Library/phylorug/extdata/beetles_50p
dict  <- utils::read.csv(file)

# Relabel tips from specimen codes to species names. The message reports how
# many tips were translated and how many were left unchanged in each tree:
translated <- translate_tips(trees, dict,
                             from_col = "specimen_code",
                             to_col   = "species_name")
#> 50p_ASTRAL_partition_entropy: 52 tips translated, 0 unchanged
#> 50p_ASTRAL_uce: 52 tips translated, 0 unchanged
#> 50p_ghost: 52 tips translated, 0 unchanged
#> 50p_partition_entropy: 52 tips translated, 0 unchanged
#> 50p_uce: 52 tips translated, 0 unchanged
```
