# phylorug_plot.R
# Layered plotting system for phylorug using base graphics + ape
#
# Usage:
#   phylorug_plot(backbone, npm) +
#     rug_layer() +
#     dot_layer() +
#     support_labels() +
#     legend_layer()

#' Create a layered phylorug plot object
#'
#' Builds a composable plot object that can be extended with `+` to add
#' visual layers (rugs, dots, labels, legends). Nothing is drawn until
#' the object is printed or explicitly plotted.
#'
#' @param backbone A `phylo` object.
#' @param npm A named list from [node_presence_matrix()] containing at least
#'   a `presence` matrix.
#' @param mode `"presence"` or `"support"`.
#' @param support_idx Integer. Which support matrix to visualize. Default `1`.
#' @param support_type Named character vector mapping trees to support types.
#' @param thresholds Optional list overriding bin thresholds.
#' @param include_backbone Logical. Include backbone as cell 1. Default `FALSE`.
#' @param file Output file path (`.pdf`, `.png`, `.jpg`), or `NULL` for
#'   on-screen rendering.
#' @param width,height Canvas dimensions in inches. `NULL` for auto.
#' @param ... Additional arguments passed to [ape::plot.phylo()].
#'
#' @returns An S3 object of class `"phylorug_plot"`. Add layers with `+`.
#'
#' @export
#'
#' @examples
#' backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
#' npm <- list(
#'   presence = matrix(c(1,0,1,1, 1,1,0,1), nrow = 4, ncol = 2,
#'     dimnames = list(as.character(6:9), c("tree_1", "tree_2")))
#' )
#'
#' # Build and render
#' p <- phylorug_plot(backbone, npm) +
#'   rug_layer() +
#'   dot_layer()
#' \dontrun{
#' p   # prints to screen
#' }
phylorug_plot <- function(backbone, npm,
                          mode             = c("presence", "support"),
                          support_idx      = 1,
                          support_type     = NULL,
                          thresholds       = NULL,
                          include_backbone = FALSE,
                          file             = NULL,
                          width            = NULL,
                          height           = NULL,
                          ...) {

  if (!inherits(backbone, "phylo"))
    stop("`backbone` must be a phylo object.", call. = FALSE)
  if (!is.list(npm) || !("presence" %in% names(npm)))
    stop("`npm` must contain a `presence` matrix.", call. = FALSE)

  mode <- match.arg(mode)

  presence <- npm$presence
  if (ncol(presence) < 1L)
    stop("`npm` has no comparison trees.", call. = FALSE)

  ntip <- ape::Ntip(backbone)

  # Build support matrix if needed
  support <- NULL
  if (mode == "support") {
    support <- npm[[paste0("support_", support_idx)]]
  }

  # Include backbone as first column
  if (include_backbone) {
    bb_p <- matrix(1, nrow = nrow(presence), ncol = 1,
                   dimnames = list(rownames(presence), "backbone"))
    presence <- cbind(bb_p, presence)
    if (mode == "support") {
      bb_s <- node_support(backbone)[[paste0("support_", support_idx)]]
      bb_col <- matrix(bb_s[as.integer(rownames(presence)) - ntip],
                       ncol = 1,
                       dimnames = list(rownames(presence), "backbone"))
      support <- cbind(bb_col, support)
    }
  }

  unanimous <- apply(presence, 1, function(row) all(row == 1, na.rm = FALSE))
  tier      <- resolve_tier(support, support_type)

  obj <- structure(
    list(
      backbone     = backbone,
      presence     = presence,
      support      = support,
      support_type = support_type,
      thresholds   = thresholds,
      mode         = mode,
      tier         = tier,
      ntip         = ntip,
      unanimous    = unanimous,
      tree_args    = list(...),
      file         = file,
      width        = width,
      height       = height,
      layers       = list()
    ),
    class = "phylorug_plot"
  )
  obj
}


#' Add a layer to a phylorug plot
#'
#' @param e1 A `phylorug_plot` object.
#' @param e2 A `phylorug_layer` object from [rug_layer()], [dot_layer()],
#'   [support_labels()], or [legend_layer()].
#'
#' @returns The modified `phylorug_plot` object.
#'
#' @export
"+.phylorug_plot" <- function(e1, e2) {
  if (!inherits(e2, "phylorug_layer"))
    stop("Can only add phylorug_layer objects with `+`.", call. = FALSE)
  e1$layers <- c(e1$layers, list(e2))
  e1
}


#' @export
print.phylorug_plot <- function(x, ...) {
  render_phylorug(x)
  invisible(x)
}


#' @export
plot.phylorug_plot <- function(x, ...) {
  render_phylorug(x)
  invisible(x)
}


# --- Internal rendering engine ---------------------------------------------

#' Render a phylorug_plot object
#'
#' Draws the tree and all layers in order. This is the engine behind
#' `print()` and `plot()`.
#'
#' @param obj A `phylorug_plot` object.
#' @noRd
render_phylorug <- function(obj) {

  backbone  <- obj$backbone
  presence  <- obj$presence
  ntip      <- obj$ntip
  n_tree    <- ncol(presence)

  # --- Scaling (from current plot_phylorug) --------------------------------
  per_tip  <- max(0.08, min(0.22, 0.35 - 0.001 * ntip))
  taxa_cex <- max(0.25, min(1.0, 1.4 - 0.004 * ntip))

  dots <- obj$tree_args
  if (is.null(dots$cex))          dots$cex          <- taxa_cex
  if (is.null(dots$edge.width))   dots$edge.width   <- max(0.3, 1.5 - 0.004 * ntip)
  if (is.null(dots$label.offset)) dots$label.offset  <- 0.001
  if (is.null(dots$show.tip.label)) dots$show.tip.label <- TRUE
  if (is.null(dots$no.margin))    dots$no.margin    <- TRUE
  if (is.null(dots$mar))          dots$mar          <- c(0.5, 0.5, 0.5, 2.5)
  dots$family <- "Helvetica"

  # --- Device routing ------------------------------------------------------
  if (!is.null(obj$file)) {
    canvas <- auto_canvas(backbone = backbone, ntip = ntip,
                          n_tree = n_tree, mode = obj$mode,
                          has_legend = has_layer(obj, "legend"),
                          per_tip = per_tip, taxa_cex = taxa_cex)
    w <- if (is.null(obj$width))  canvas$width  else obj$width
    h <- if (is.null(obj$height)) canvas$height else obj$height

    ext <- tolower(tools::file_ext(obj$file))
    if (ext == "pdf") {
      grDevices::pdf(obj$file, width = w, height = h)
    } else if (ext == "png") {
      grDevices::png(obj$file, width = w, height = h, units = "in", res = 300)
    } else if (ext %in% c("jpg", "jpeg")) {
      grDevices::jpeg(obj$file, width = w, height = h, units = "in", res = 300)
    } else {
      stop("Unsupported file extension. Use .pdf, .png, or .jpg", call. = FALSE)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # --- y.lim padding for legend --------------------------------------------
  if (is.null(dots$y.lim) && has_layer(obj, "legend")) {
    top_pad <- ceiling(ntip * 0.12)
    dots$y.lim <- c(-1.0, ntip + top_pad)
  }

  # --- Draw tree -----------------------------------------------------------
  do.call(ape::plot.phylo, c(list(x = backbone), dots))
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  # --- Cell geometry -------------------------------------------------------
  yy_tip <- sort(last_pp$yy[seq_len(ntip)])
  dy     <- stats::median(diff(yy_tip))
  pin    <- graphics::par("pin")

  # Shared context passed to every layer
  ctx <- list(
    backbone     = backbone,
    presence     = presence,
    support      = obj$support,
    support_type = obj$support_type,
    thresholds   = obj$thresholds,
    mode         = obj$mode,
    tier         = obj$tier,
    ntip         = ntip,
    n_tree       = n_tree,
    unanimous    = obj$unanimous,
    last_pp      = last_pp,
    dy           = dy,
    pin          = pin
  )

  # --- Execute layers in order ---------------------------------------------
  for (layer in obj$layers) {
    layer$draw(ctx)
  }

  invisible(NULL)
}


#' Check if a layer type exists in the plot
#' @noRd
has_layer <- function(obj, type) {
  any(vapply(obj$layers, function(l) identical(l$type, type), logical(1)))
}
