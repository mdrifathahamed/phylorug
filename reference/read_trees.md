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

  Character vector of file extensions to search for, without the leading
  dot. Matching is case-insensitive. You may supply your own, for
  example `ext = "treefile"`. If you do not, the default filters the
  directory to the common tree file extensions, so that alignments, log
  files and configuration files sitting beside the trees are not read.
  Files with no extension, such as classic RAxML output, cannot be
  matched.

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
scored as a pool. How a pool's clade recovery is summarized (as a
continuous proportion, or binarised at a threshold) is decided later, by
[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md),
not here.

Do not supply posterior samples, bootstrap replicates, or sets of gene
trees. These are distributions rather than analyses and must be
summarised before use. A file holding more than 100 trees is an error.

Support values written as Newick node labels are read normally.
BEAST-style bracket annotations (`[&posterior=0.98]`) are discarded by
ape; the topology is unaffected, and a message is emitted.

Support values are imported exactly as written in the tree file;
`read_trees()` does not recompute or verify them. When a file holds
several tied-optimal trees, be aware that upstream programs differ in
how they summarize support across such trees: TNT and POY4 default to
the more conservative strict-consensus approach, whereas PAUP\* and
PHYLIP default to the frequency-within-replicates approach, which
Simmons and Freudenstein (2011) showed can inflate apparent support for
unsupported clades. This choice is made by the upstream software before
the file reaches.

Files matching `ext` that contain no tree, such as a NEXUS character
matrix, are skipped with a message. A tree file that fails to parse is
an error.

## References

Simmons, M.P. & Freudenstein, J.V. (2011). Spurious 99% bootstrap and
jackknife support for unsupported clades. *Molecular Phylogenetics and
Evolution*, 61(1), 177-191.
[doi:10.1016/j.ympev.2011.06.003](https://doi.org/10.1016/j.ympev.2011.06.003)

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
