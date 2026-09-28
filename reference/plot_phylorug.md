# Draw a phylorug: a backbone tree with rug plots on every internal nodes

`plot_phylorug()` overlays rug plots on a backbone phylogeny, compares
how multiple analyses treat each internal node. In presence mode, cells
are black (recovered) or white (absent). In support mode, cells are
shaded by binned support strength. The function handles canvas sizing,
legend placement, and font scaling automatically.

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
  show_support = TRUE,
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
  dot_on_very_high = TRUE,
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

  Optional list overriding the built-in bin thresholds used in support
  mode. Support values are binned into four levels: very high, high,
  moderate, and low. Each entry supplies the lower bounds for the top
  three levels as a named numeric vector with exactly three names:
  `very_high`, `high`, and `moderate`. Any value below `moderate` is
  automatically binned as low, so no fourth cutpoint is needed. Key each
  entry by metric name when `support_type` was declared in
  [`node_presence_matrix()`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md),
  or by `"universal"` when it was not. Examples:
  `thresholds = list(ufboot = c(very_high = 97, high = 85, moderate = 50))`
  means UFBoot2 values \>= 97 are very high, 85–96 are high, 50–84 are
  moderate, and \< 50 are low.
  `thresholds = list(universal = c(very_high = 0.90, high = 0.70, moderate = 0.50))`
  applies a single scale to all trees. If `NULL` (default) and
  `support_type` was declared, default thresholds for each metric are
  used. If neither `support_type` nor `thresholds` is provided,
  universal thresholds (0.95 / 0.80 / 0.50) are applied.

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

  Logical. Default `TRUE`. Controls both legends: the position legend
  (top-left, numbered grid mapping each cell to its analysis) and, in
  support mode, the threshold legend (top-right, colour key for the four
  support bins). Set to `FALSE` to suppress both, for example when
  annotating the figure manually in Inkscape or Illustrator.

- show_support:

  Logical. If `TRUE` (default), backbone node support labels are drawn
  beside each node, and the postion legend includes a colored line
  naming what those numbers are.

- show_support_idx:

  Optional. Controls what appears in the label on the tree AND what the
  legend calls it. Accepts several shapes:

  - `NULL` (default): the raw `node.label` is drawn on the tree and the
    legend says `<value> (backbone support)` without naming metrics.

  - An integer (`1` or `2`): only that slot of a compound label (e.g.
    "80/95") is drawn, and the legend shows just that number.

  - A named integer, e.g. `c("1" = "sh_alrt")`: same slot behaviour, and
    the legend adds the metric name (e.g. "SH-aLRT").

  - Two named integers, e.g. `c("1" = "sh_alrt", "2" = "ufboot")`: the
    full compound is drawn and the legend names both metrics.

  - A plain string, e.g. `"UFBoot2"`: legal only when the tree has
    single-value labels (no `/`). On compound labels, this is ambiguous
    (which slot?) and the function falls back to the raw-label behaviour
    with a message explaining the shape it accepts.

  Recognised metric keys for named-integer form: `"ufboot"`,
  `"sh_alrt"`, `"lpp"`, `"posterior"`, `"jackknife"`, `"bootstrap"`,
  `"bremer_ratio"`, `"transfer_boot"`. Any other key is shown verbatim
  in the legend.

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

  Logical. Default `TRUE`. **Presence mode only.** draws a small dot on
  backbone nodes where every comparison tree is identical in recovering
  the clade. Has no effect in support mode; see `dot_on_very_high` for
  the support-mode equivalent.

- dot_col, dot_cex:

  Colour and size of the identical-clade dot. Defaults are `"black"` and
  NULL (auto-scales). Shared by both `dot_identical` (presence mode) and
  `dot_on_very_high` (support mode).

- support_label_cex:

  Numeric or `NULL`. Size of the backbone support labels. Default `NULL`
  auto-scales with tree size.

- support_label_col:

  Colour of the backbone support labels. Default `"red"`.

- dot_on_very_high:

  Logical. Default `TRUE`. **Support mode only.** Controls what a dot
  means there. When `TRUE`, a dot is drawn only where a clade is
  recovered by every comparison analysis AND every one of those analyses
  independently rates it very-high support. Anything less (not recovered
  everywhere, or recovered everywhere but not uniformly very-high) is
  drawn as a rug instead. Set to `FALSE` to turn off dots entirely in
  support mode: every node is then drawn as a rug, however strong or
  weak its support, the most conservative and transparent option. Has no
  effect in presence mode; see `dot_identical` for the presence-mode
  equivalent.

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

## Details

Each internal node is drawn as either a dot (compact summary) or a rug
plot (full per-analysis grid). What a dot means depends on the mode. In
presence mode, a dot means every comparison analysis recovered the
clade; control this with `dot_identical`, setting it to `FALSE` to show
rug plots on every node, even unanimous ones. In support mode, a dot
means every comparison analysis recovered the clade and every one
independently rates it very-high support, control this with
`dot_on_very_high`, setting it to `FALSE` to show rug plots on every
node regardless of support level. Dot appearance (`dot_col`, `dot_cex`)
is shared across both modes.

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
# --- Restrict to specific nodes -------------------------------------------
# By taxa: name two or more tips whose MRCA is the node you want.
# Useful when you know which species define a clade of interest.
tmp4 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm_st,
              file  = tmp4,
              nodes = list(c("Sisyphus_schaefferi_STL44",
                             "Sisyphus_muricatus__STL5")))
unlink(tmp4)

# By node ID: useful when exploring interactively after inspecting
# rownames(npm_st$presence) or ape::nodelabels() output.
tmp5 <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm_st,
              file  = tmp5,
              nodes = as.integer(rownames(npm_st$presence))[1:2])
unlink(tmp5)
```
