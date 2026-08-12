# Parse support values from internal nodes of a phylogenetic tree

Extracts the node labels of one tree and returns their support values as
a data frame: one row per internal node, three columns (`support_1`,
`support_2`, `support_3`). How a label fills those columns:

- A single value (such as an IQ-TREE bootstrap) fills `support_1`; the
  other two are `NA`.

- A compound label split by `sep` (default `/`) fills one column per
  value:

  - Two for IQ-TREE SH-aLRT/UFBoot2 or ASTRAL pp1/pp2.

  - Three for IQ-TREE SH-aLRT/UFBoot2/aBayes or ASTRAL pp1/pp2/pp3.

- Empty or non-numeric pieces become `NA`.

- A tree with no node labels returns a data frame of all `NA`.

## Usage

``` r
node_support(tree, sep = "/", digits = NULL)
```

## Arguments

- tree:

  A phylogenetic tree of class `"phylo"`, with node labels stored in
  `tree$node.label`.

- sep:

  A single character string giving the delimiter that separates compound
  support values. Defaults to `"/"`.

- digits:

  Integer or `NULL`. If supplied, support values are rounded to this
  many decimal places. Defaults to `NULL`.

## Value

A data frame with three columns, `support_1`, `support_2` and
`support_3`, and one row per internal node of the tree. Columns not
present in a given label are `NA`.
