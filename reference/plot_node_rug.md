# Draw the rug at every internal node of a plotted tree

For each internal node of the backbone tree, paints a small grid of
rectangles, one per comparison tree, showing whether that tree recovered
the node's clade and, in support mode, how strongly. Together these
grids are the rug.

The appearance is set by a tier, resolved from which arguments are
supplied rather than named directly.
([`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
selects the tier for the user through its `mode` argument.)

- Tier 1, presence: `support` is `NULL`. Cells are black for present,
  white for absent, and grey for partial recovery in a pool.

- Tier 2, support: both `support` and `support_type` are supplied.
  Recovered cells are shaded by binned support strength, from black
  (very high) through greys to yellow (low). A cell is white when the
  tree does not recover the clade at all, and red when the tree recovers
  the clade but carries no support value for it (an unscored node).

Users do not call this directly;
[`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
calls it after drawing the tree and working out the cell geometry.

## Usage

``` r
plot_node_rug(
  npm,
  support = NULL,
  support_type = NULL,
  thresholds = NULL,
  cell_h,
  cell_w,
  n_cols,
  x_offset = 0,
  y_offset = 0,
  rug_position = c("inside", "outside"),
  last_pp = NULL
)
```

## Arguments

- npm:

  The presence matrix from
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md).
  One row per internal backbone node (node numbers as rownames), one
  column per comparison tree. Cells are `1` (recovered), `0` (not
  recovered), or a proportion (pool recovery).

- support:

  Support matrix of the same shape as `npm`, or `NULL`. When `NULL`, the
  rug is drawn in presence mode.

- support_type:

  Named character vector mapping each comparison tree to its support
  measure (`"ufboot"`, `"sh_alrt"`, `"lpp"`, `"posterior"`), or `NULL`.
  Required for binned shading.

- thresholds:

  Optional list overriding the built-in bin thresholds, keyed by support
  type.`NULL` uses the literature defaults.

- cell_h, cell_w:

  Numeric. Height and width of one cell, in the tree's plotting
  coordinates.

- n_cols:

  Integer. Columns in each node's grid.

- x_offset, y_offset:

  Numeric. Shift the whole grid away from the node, as a fraction of the
  tree's width and height.

- rug_position:

  One of `"outside"` or `"inside"`(default).

- last_pp:

  Plot coordinates from
  [`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html)
  giving the x/y position of every node and tip, used to place each rug
  grid at its node.
  [`plot_phylorug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md)
  always supplies this. If called directly with `NULL`, the coordinates
  are fetched from the most recently drawn tree.

  \#' @return Returns nothing; it draws the rug cells directly onto the
  tree.
