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
#' @param layout Tree layout: `"phylogram"` (default, with branch lengths),
#'   `"cladogram"` (equal spacing), `"fan"` (circular), or `"unrooted"`.
#' @param file Output file path (`.pdf`, `.png`, `.jpg`), or `NULL` for
#'   on-screen rendering.
#' @param width,height Canvas dimensions in inches. `NULL` for auto.
#' @param ... Additional arguments passed to [ape::plot.phylo()], such as
#'   `cex`, `font`, `edge.width`, `edge.color`, `tip.color`,
#'   `label.offset`, or `show.tip.label`.
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
                          layout           = c("phylogram", "cladogram",
                                               "fan", "unrooted"),
                          file             = NULL,
                          width            = NULL,
                          height           = NULL,
                          ...) {

  if (!inherits(backbone, "phylo"))
    stop("`backbone` must be a phylo object.", call. = FALSE)
  if (!is.list(npm) || !("presence" %in% names(npm)))
    stop("`npm` must contain a `presence` matrix.", call. = FALSE)

  mode   <- match.arg(mode)
  layout <- match.arg(layout)

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
      layout       = layout,
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

  # --- Extract theme layer (applied before drawing) -------------------------
  th <- get_layer_params(obj, "theme")
  base_size <- if (!is.null(th$base_size)) th$base_size else 1

  # --- Scaling (exact match to plot_phylorug) --------------------------------
  per_tip <- min(0.15, max(0.03, 0.31 - 0.041 * log(ntip)))

  R_CEX     <- 4.0
  R_EDGE    <- 6.0
  R_SUPPORT <- 2.3
  R_LEG_CELL <- 0.011
  R_TH_SQ    <- 0.008
  R_LEG_TEXT <- 0.045
  R_TH_TEXT  <- 0.041

  dots <- obj$tree_args

  # Layout
  dots$type <- obj$layout
  if (obj$layout == "cladogram") dots$use.edge.length <- FALSE

  # User ... args > auto defaults (scaled by base_size)
  dots$cex            <- dots$cex            %||% (min(0.65, max(0.30, per_tip * R_CEX)) * base_size)
  dots$edge.width     <- dots$edge.width     %||% min(1.3, max(0.7, per_tip * R_EDGE))
  dots$label.offset   <- dots$label.offset   %||% 0.001
  dots$show.tip.label <- dots$show.tip.label %||% TRUE
  dots$font           <- dots$font           %||% 3
  dots$no.margin      <- dots$no.margin      %||% TRUE

  taxa_cex <- dots$cex

  # Theme settings
  dots$mar    <- dots$mar    %||% (if (!is.null(th$mar)) th$mar else c(0.5, 0.5, 0.5, 2.5))
  dots$family <- if (!is.null(th$family)) th$family else "Helvetica"

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

  # --- Dynamic y.lim (exact copy from plot_phylorug) -----------------------
  if (is.null(dots$y.lim) && has_layer(obj, "legend")) {
    din    <- graphics::par("din")
    mai    <- graphics::par("mai")
    plot_h <- din[2] - mai[1] - mai[3]

    w_ref        <- din[1] * 0.85 + din[2] * 0.15
    est_pos_cex  <- min(0.55, max(0.20, w_ref * R_LEG_TEXT))
    est_th_cex   <- min(0.50, max(0.18, w_ref * R_TH_TEXT))

    pos_line_h <- 0.2 * est_pos_cex
    th_line_h  <- 0.2 * est_th_cex
    pos_leg_h  <- n_tree * pos_line_h
    th_leg_h   <- if (obj$mode == "support") 6 * th_line_h else 0
    top_leg_h  <- max(pos_leg_h, th_leg_h)
    gap_inches <- 3 * max(pos_line_h, th_line_h)

    avail   <- max(1, plot_h - top_leg_h - gap_inches)
    top_pad <- ceiling((top_leg_h + gap_inches) * ntip / avail)

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


#' Extract params from a pre-draw layer (tree_style or theme)
#' @noRd
get_layer_params <- function(obj, type) {
  for (l in obj$layers) {
    if (identical(l$type, type)) return(l$params)
  }
  list()  # empty list if layer not added
}
