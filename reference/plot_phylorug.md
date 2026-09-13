# Draw a phylorug: a backbone tree with rug plots on every internal nodes

`plot_phylorug()` overlays clade stability grids (rugs) on a backbone
phylogeny, compares how multiple analyses treat each internal node. In
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
  thresholds = NULL,
  n_rows = NULL,
  n_cols = NULL,
  include_backbone = FALSE,
  nodes = NULL,
  legend = TRUE,
  show_support = FALSE,
  show_support_idx = NULL,
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
  hide_unsupported = TRUE,
  ...
)
```

## Arguments

- backbone:

  A `phylo` object representing the backbone tree.

- npm:

  A named list containing the `presence` and `support` matrices,
  returned by
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
  based on the tree size and legend layout.

- mode:

  One of `"presence"` (default) or `"support"`.

- support_idx:

  Integer (1, 2, or 3). Default is `1`. Specifies which single support
  matrix from the `npm` list to visualize. For example, if you generated
  the `npm` using `support_col = c(1, 2)`, passing `2` here tells the
  plotting engine to physically shade the grid cells using the second
  metric (stored in your list as `support_2`).

- thresholds:

  Optional list overriding built-in bin thresholds. Made by metric name
  for named metrics, or `"universal"` when `support_type` is not
  declared. For example:
  `thresholds = list(ufboot = c(very_high = 97, high = 85, moderate = 50))`
  or
  `thresholds = list(universal = c(very_high = 0.90, high = 0.70, moderate = 0.50))`.
  If `NULL` (default) and `support_type` was declared in
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md),
  literature-based thresholds for each metric are used. If neither
  `support_type` nor `thresholds` is provided, universal thresholds
  (0.95 / 0.80 / 0.50) are applied to all trees.

- n_rows, n_cols:

  Integer. Grid shape for the rug at each node. If `NULL` (the default),
  a roughly square grid is chosen automatically.

- include_backbone:

  Logical or character. If `FALSE` (default), the backbone does not
  occupy a rug cell. If `TRUE`, the backbone is added as cell 1 using
  universal thresholds in support mode (with a message). If a character
  string naming a support metric (`include_backbone = "ufboot"`), the
  backbone is added as cell 1 and binned against that metric's
  thresholds.

- nodes:

  Optional. Restricts the plot to specific backbone nodes. Three ways to
  specify:

  - **By taxa (recommended):** a list of character vectors, each with
    \>= 2 tip labels. phylorug finds the MRCA of each set. Example:
    `nodes = list(c("sp_A", "sp_B"), c("sp_C", "sp_D"))`.

  - **By node ID (advanced):** a numeric vector of ape internal node
    IDs, as shown in `rownames(npm$presence)`. To see which taxa belong
    to a node, use `ape::extract.clade(backbone, node_id)$tip.label`.

  - **By node label:** a character vector matching
    `backbone$node.label`.

  If `NULL` (default), all internal nodes are plotted. When supplied,
  only selected nodes get dots, rugs, or support labels.

- legend:

  Logical. Draw the legend. Default TRUE.

- show_support:

  Logical. If `TRUE`, backbone node support labels are drawn beside each
  node. Default `FALSE`.

- show_support_idx:

  Integer or NULL. Which value from compound node labels (e.g. "80/95")
  to display when `show_support = TRUE`. Default `NULL` displays the
  full compound label as-is. Set to `1` or `2` to display only a single
  metric.

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

  Logical. Default `FALSE`. When this and `dot_identical` are both
  `TRUE`, unanimous nodes receive both the dot and a support rug, making
  per-tree support strength visible even at universally recovered
  clades.

- hide_unsupported:

  Logical. Default `TRUE`. Nodes where no comparison tree recovers the
  clade are left bare, the absence of both a dot and a rug signals that
  the clade is unique to the backbone topology. Set to `FALSE` to draw
  all-white rugs at these nodes, which can be useful when you want every
  internal node to carry a visible grid for annotation or figure
  editing.

- ...:

  Additional arguments passed to
  [`ape::plot.phylo()`](https://rdrr.io/pkg/ape/man/plot.phylo.html),
  such as `cex`, `edge.width`, `font`, or `label.offset`. These override
  the automatic scaling when provided.

## Value

Invisibly, the file path if a file was written, or `NULL` if plotted
directly to the active graphics device.

## See also

[`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md)
to build the input data,
[`check_taxa()`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md)
to verify taxon sets, and
[`add_tree()`](https://mdrifathahamed.github.io/phylorug/reference/add_tree.md)
to append trees to an existing matrix.

## Examples

``` r
# Build the plotting input from real trees:
backbone <- sample_trees[["70p_uce"]]
others   <- sample_trees[names(sample_trees) != "70p_uce"]
support_type <- c(
  "70p_ASTRAL_partition_entropy" = "lpp",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ghost"                    = "ufboot",
  "70p_partition_entropy"        = "ufboot"
)
npm_st <- node_presence_matrix(backbone, others,
                               support_col = 1,
                               support_type = support_type)

# --- Presence mode --------------------------------------------------------
# Each cell shows whether an analysis recovered the backbone clade.
tmp <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm_st, file = tmp)
unlink(tmp)

# --- Support mode ---------------------------------------------------------
tmp2 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm_st,
              file         = tmp2,
              mode         = "support",
              support_idx  = 1)
unlink(tmp2)

# --- Some optional controls -----------------------------------------------
tmp3 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm_st,
              file             = tmp3,
              include_backbone = TRUE,
              rug_position     = "outside")
unlink(tmp3)
# --- Restrict to specific nodes --------------------------------------------
# By taxa (recommended): select a clade by the tips that define it.
tmp4 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm_st,
              file  = tmp4,
              nodes = list(backbone$tip.label[1:2]))
unlink(tmp4)

# By node ID, if you already know it (see rownames(npm_st$presence)):
tmp5 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm_st,
              file  = tmp5,
              nodes = as.integer(rownames(npm_st$presence))[1:2])
unlink(tmp5)
```
