# Import phylogenetic trees from a directory

Scans a directory for phylogenetic tree files, parses each one, and
returns them as a named list. Format is detected from file content, so a
directory may mix Newick and NEXUS files.

## Usage

``` r
read_trees(
  dir = ".",
  ext = c("tre", "tree", "treefile", "nwk", "newick", "nex", "nexus", "contree"),
  format = c("auto", "newick", "nexus"),
  verbose = TRUE
)
```

## Arguments

- dir:

  Path to the directory containing the tree files. Defaults to the
  working directory.

- ext:

  Character vector of file extensions to look for, without the leading
  dot (e.g. `"tre"`, not `".tre"`). Case-insensitive. The default covers
  common tree extensions (`tre`, `tree`, `treefile`, `nwk`, `nex`,
  `nexus`, `contree`). Pass your own to narrow or widen the search.

- format:

  The parsing strategy to use. Defaults to `"auto"`, which detects the
  format of each file individually, so a directory may mix formats. Use
  `"newick"` or `"nexus"` to force one parser for all files.

- verbose:

  Logical. If `TRUE` (default), reports trees read, files skipped, and
  pooled analyses. Errors are always reported.

## Value

A named list, one element per tree file, named after the file names with
extensions removed. Single-tree files yield `"phylo"` objects;
multi-tree files yield `"multiPhylo"`. Carries a `"pool_sizes"`
attribute giving the number of trees per analysis.

## Details

Each file is one analysis. A file holding a single tree is returned as a
`"phylo"` ; a file holding several equally optimal trees from one search
(as from POY, TNT, or PAUP\*) is returned as a `"multiPhylo"` and is
scored as a pool. How many pool trees must recover a clade before it
counts as present is controlled by the `pool_threshold` argument in
[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md).

Do not supply posterior samples, bootstrap replicates, or sets of gene
trees. These are distributions rather than analyses and must be
summarised in the upstream analysis before use (e.g. into a consensus
tree). Files holding more than 100 trees are rejected with an error.

Support values stored as numeric node labels in the tree file (e.g.
`((A,B)95,(C,D)100);`) are read normally. BEAST-style bracket
annotations (`[&posterior=0.98]`) are discarded by ape; the topology is
unaffected, and a message is emitted. Support values are imported
exactly as written in the tree file; `read_trees()` does not recompute
or verify them.

Files matching `ext` that contain no tree, such as a NEXUS character
matrix, are skipped with a message. A tree file that fails to parse is
an error.

## Examples

``` r
# phylorug ships example  raw tree files in its `extdata` directory.
# Point read_trees() at one of those folders:
dir <- system.file("extdata", "beetles_70p", package = "phylorug")

# With only the directory, read_trees() detects each file's format from its
# contents (NEXUS or Newick, the two formats this version supports) and
# matches the default extension filter (tre, tree, treefile, nwk, ...).

trees <- read_trees(dir)
#> Read 5 analyses (5 trees) from: /home/runner/work/_temp/Library/phylorug/extdata/beetles_70p

# The result is a named list, one element per analysis, ready to pass to
# [check_taxa()] or [node_presence_matrix()]:

names(trees)
#> [1] "70p_ASTRAL_partition_entropy" "70p_ASTRAL_uce"              
#> [3] "70p_ghost"                    "70p_partition_entropy"       
#> [5] "70p_uce"                     

# If a folder holds many files and you want only some, narrow `ext`
# to one extension, or to a set of them:
trees <- read_trees(dir, ext = "tre")
#> Read 5 analyses (5 trees) from: /home/runner/work/_temp/Library/phylorug/extdata/beetles_70p
trees <- read_trees(dir, ext = c("tre", "tree", "treefile"))
#> Read 5 analyses (5 trees) from: /home/runner/work/_temp/Library/phylorug/extdata/beetles_70p

# For full control you can also force a single parser. Note this applies
# one format to every matched file, so a file of a different format would
# not be read.
trees <- read_trees(dir, format = "newick")
#> Read 5 analyses (5 trees) from: /home/runner/work/_temp/Library/phylorug/extdata/beetles_70p
```
