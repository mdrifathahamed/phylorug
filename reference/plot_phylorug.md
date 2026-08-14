# Draw a phylorug: a backbone tree with node rugs

`plot_phylorug()` overlays clade stability grids (rugs) on a backbone
phylogeny, comparing how multiple analyses treat each internal node. In
presence mode, cells are black (recovered) or white (absent). In support
mode, cells are shaded by binned support strength. The function handles
canvas sizing, legend placement, and font scaling automatically.

## Usage

``` r
plot_phylorug(
  backbone,
  npm,
  file = NULL,
  width = NULL,
  height = NULL,
  mode = c("presence", "support"),
  support_idx = 1,
  support_type = NULL,
  thresholds = NULL,
  n_rows = NULL,
  n_cols = NULL,
  include_backbone = FALSE,
  legend = TRUE,
  show_support = FALSE,
  cell_scale = 0.45,
  x_offset = 0,
  y_offset = 0,
  rug_position = c("inside", "outside"),
  dot_identical = TRUE,
  dot_col = "black",
  dot_cex = NULL,
  support_label_cex = NULL,
  support_label_col = "red",
  rug_on_identical = FALSE,
  hide_unsupported = FALSE,
  ...
)
```

## Arguments

- backbone:

  A `phylo` object representing the backbone tree.

- npm:

  A named list containing the `presence` and `support` matrices, exactly
  as returned by
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md).

- file:

  Optional character string specifying the output file name or full
  directory path. If left as `NULL` (the default), the plot renders in
  the active graphics device (e.g., RStudio) for quick drafting. To
  avoid aspect ratio distortion caused by GUI exports and to generate
  perfectly scaled, publication-ready figures, provide a file path here
  (must end in `.pdf`, `.png`, or `.jpg`) to utilize the package's
  internal scaling engine.

- width, height:

  Numeric. Optional canvas dimensions in inches, applied only when
  exporting to a `file`.If left as `NULL` (the default), the package's
  internal engine dynamically calculates the optimal canvas dimensions
  based on the tree size and legend layout. Providing values here
  overrides the automatic scaling, which is useful for meeting strict
  journal dimension requirements.

- mode:

  One of `"presence"` (default) or `"support"`.

- support_idx:

  Integer (1, 2, or 3). Default is `1`. Specifies which single support
  matrix from the `npm` list to visualize when `mode = "support"`. This
  corresponds directly to your extraction order in
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md).
  For example, if you generated the data using `support_col = c(1, 2)`,
  passing `2` here tells the plotting engine to physically shade the
  grid cells using the second metric (stored in your list as
  `support_2`).

- support_type:

  Named character vector mapping comparison trees to their support
  metrics (e.g., `"ufboot"`, `"sh_alrt"`, `"lpp"`). Required only when
  `mode = "support"`. Names must match your comparison tree names.

- thresholds:

  Optional list overriding built-in bin thresholds. While default
  thresholds are provided based on common literature, it is highly
  recommended to define your own custom thresholds to suit your specific
  analytical framework.

- n_rows, n_cols:

  Integer. Grid shape for the rug at each node. If `NULL` (the default),
  a roughly square grid is chosen automatically.

- include_backbone:

  Logical. If `TRUE`, the backbone tree occupies cell 1 of every rug,
  showing its own presence (always `1`) or its own support value.
  Default `FALSE`, since the backbone defines the topology and trivially
  recovers every clade. Useful when the figure will be edited in a
  vector editor and every analysis must appear in the grid.

- legend:

  Logical. Draw the legend. Default TRUE.

- show_support:

  Logical. If `TRUE`, backbone node support labels are drawn beside each
  node. Default `FALSE`.

- cell_scale:

  Numeric multiplier on cell height. Default 0.45.

- x_offset, y_offset:

  Numeric grid shift. Default 0.

- rug_position:

  One of `"inside"` (default) or `"outside"`. Controls where the node
  rug grid is placed relative to the backbone node: `"inside"` tucks the
  grid into the crook above-left (toward the root), while `"outside"`
  places the grid to the right of the node (toward the tips).

- dot_identical:

  Logical. Default `TRUE`. If `TRUE`, draws a small dot on backbone
  nodes where every comparison tree is identical in recovering the
  clade.

- dot_col, dot_cex:

  Colour and size of the identical-clade dot. Defaults are `"black"` and
  NULL (auto-scales).

- support_label_cex:

  Numeric or `NULL`. Size of the backbone support labels. Default `NULL`
  auto-scales with tree size.

- support_label_col:

  Colour of the backbone support labels. Default `"red"`.

- rug_on_identical:

  Logical. Default `FALSE`. When it and `dot_identical`both are `TRUE`,
  unanimous nodes receive both the dot and a support rug, making
  per-tree support strength visible even at universally recovered
  clades.

- hide_unsupported:

  Logical. Default `FALSE`. When `TRUE`, nodes where no comparison tree
  recovers the clade are left bare, no rug is drawn. The absence of both
  a dot and a rug signals that the clade is unique to the backbone
  topology.

- ...:

  Additional arguments passed to
  [`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html),
  such as `cex`, `edge.width`, `font`, or `label.offset`. These override
  the automatic scaling when provided.

## Value

Invisibly, the file path if a file was written, or `NULL` if plotted
directly to the active graphics device (not recommended).

## See also

[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
to build the input data,
[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md)
to verify taxon sets, and
[`plot_node_rug()`](https://mdrifathahamed.github.io/phylorug/reference/plot_node_rug.md)
which handles the cell-level drawing(not used by the user).

## Examples

``` r
# Build the plotting input from real trees:
backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[names(sample_trees) != "70p_uce"]
npm <- node_presence_matrix(backbone, others, support_col = 1)

# --- Presence mode --------------------------------------------------------
# Each cell shows whether an analysis recovered the backbone clade.
# Writing to a file uses the internal scaling engine for a clean figure.
tmp <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm, file = tmp)
unlink(tmp)

# --- Support mode ---------------------------------------------------------
# Cells are shaded by support strength. `support_type` tells plot_phylorug
# how to read each tree's values: ASTRAL trees carry local posterior
# probability ("lpp"), IQ-TREE trees carry UFBoot2 ("ufboot").
support_type <- c(
  "70p_ASTRAL_partition_entropy" = "lpp",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ghost"                    = "ufboot",
  "70p_partition_entropy"        = "ufboot"
)
tmp2 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm,
              file         = tmp2,
              mode         = "support",
              support_idx  = 1,
              support_type = support_type)
unlink(tmp2)

# --- Some optional controls -----------------------------------------------
# include_backbone = TRUE adds the backbone as its own cell in every rug;
# rug_position = "outside" places grids toward the tips instead of the crook;
# hide_unsupported = TRUE leaves clades no analysis recovered bare.
tmp3 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm,
              file             = tmp3,
              include_backbone = TRUE,
              rug_position     = "outside")
unlink(tmp3)
```
