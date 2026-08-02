# gg_tree_layout.R
# Extracts node/tip/edge coordinates from ape's layout engine, without
# reimplementing tree geometry. We let ape::plot.phylo() do the actual
# layout math (it already handles phylogram/cladogram/unrooted spacing
# correctly), but with plot = FALSE so nothing is drawn, then harvest the
# coordinates into plain data frames that ggplot2 can consume.

#' Compute tree layout for ggplot2 rendering
#'
#' @param backbone A `phylo` object.
#' @param ... Passed to `ape::plot.phylo()` (e.g. `type`, `use.edge.length`).
#'
#' @returns A list with `edges` (horizontal + vertical segments), `tips`
#'   (label positions), `nodes` (internal node positions), `ntip`, and the
#'   raw `last_pp` object from ape.
#' @noRd
gg_tree_layout <- function(backbone, ...) {

  if (!inherits(backbone, "phylo"))
    stop("`backbone` must be a phylo object.", call. = FALSE)

  ntip <- ape::Ntip(backbone)

  # Open a throwaway device so plot.phylo has somewhere to compute par(),
  # even with plot = FALSE. Nothing is ever shown to the user.
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)

  ape::plot.phylo(backbone, plot = FALSE, ...)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  xx <- last_pp$xx
  yy <- last_pp$yy
  edge <- backbone$edge  # matrix: col 1 = parent, col 2 = child

  # Rectangular (phylogram/cladogram) elbow edges: one horizontal segment
  # (the branch itself) and one vertical connector per edge. This matches
  # what ape draws internally for type = "phylogram" / "cladogram".
  edges_h <- data.frame(
    x    = xx[edge[, 1]],
    xend = xx[edge[, 2]],
    y    = yy[edge[, 2]],
    yend = yy[edge[, 2]]
  )
  edges_v <- data.frame(
    x    = xx[edge[, 1]],
    xend = xx[edge[, 1]],
    y    = yy[edge[, 1]],
    yend = yy[edge[, 2]]
  )

  tips <- data.frame(
    node  = seq_len(ntip),
    label = backbone$tip.label,
    x     = xx[seq_len(ntip)],
    y     = yy[seq_len(ntip)]
  )

  n_internal <- backbone$Nnode
  node_ids   <- (ntip + 1L):(ntip + n_internal)
  nodes <- data.frame(
    node = node_ids,
    x    = xx[node_ids],
    y    = yy[node_ids]
  )

  list(
    edges_h = edges_h,
    edges_v = edges_v,
    tips    = tips,
    nodes   = nodes,
    ntip    = ntip,
    last_pp = last_pp
  )
}
