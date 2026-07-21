#' Draw a phylorug: a backbone tree with node rugs
#'
#' @param backbone The reference tree of class \code{"phylo"}.
#' @param npm The list returned by \code{\link{node_presence_matrix}}.
#' @param file Optional output file path (.pdf, .png, .jpg).
#' @param width,height Optional canvas dimensions in inches.
#' @param mode One of \code{"presence"} (default) or \code{"support"}.
#' @param support_col Integer 1-3. Which support matrix to use.
#' @param support_type Named character vector mapping trees to support measures.
#' @param thresholds Optional list overriding built-in bin thresholds.
#' @param n_rows,n_cols Optional grid shape at each node.
#' @param include_backbone Logical. Include backbone as cell 1. Default FALSE.
#' @param legend Logical. Draw the legend. Default TRUE.
#' @param show_support Logical. Show backbone support labels in red.
#' @param cell_scale Numeric multiplier on cell height. Default 0.3.
#' @param x_offset,y_offset Numeric grid shift. Default 0.
#' @param rug_position One of \code{"inside"} or \code{"outside"}.
#' @param dot_unanimous Logical. Dot for unanimous nodes. Default TRUE.
#' @param dot_col,dot_cex Dot appearance.
#' @param ... Passed to \code{ape::plot.phylo()}.
#' @return Invisibly the file path or NULL.
#' @export
plot_phylorug <- function(backbone, npm,
                          file             = NULL,
                          width            = NULL,
                          height           = NULL,
                          mode             = c("presence", "support"),
                          support_col      = 1,
                          support_type     = NULL,
                          thresholds       = NULL,
                          n_rows           = NULL,
                          n_cols           = NULL,
                          include_backbone = FALSE,
                          legend           = TRUE,
                          show_support     = NULL,
                          cell_scale       = 0.45,
                          x_offset         = 0,
                          y_offset         = 0,
                          rug_position     = c("inside", "outside"),
                          dot_unanimous    = TRUE,
                          dot_col          = "black",
                          dot_cex          = 0.45,
                          ...) {

  # --- 1. Validate -----------------------------------------------------------
  if (!inherits(backbone, "phylo"))
    stop("`backbone` must be a phylo object.", call. = FALSE)
  if (!is.list(npm) || !("presence" %in% names(npm)))
    stop("`npm` must be a valid node presence matrix.", call. = FALSE)

  mode         <- match.arg(mode)
  rug_position <- match.arg(rug_position)

  presence <- npm$presence
  if (ncol(presence) < 1L)
    stop("`npm` has no comparison trees to plot.", call. = FALSE)
  if (is.null(show_support))
    show_support <- !(mode == "support" && include_backbone)

  ntip <- ape::Ntip(backbone)

  # --- 2. Build matrices -----------------------------------------------------
  support <- NULL
  if (mode == "support") {
    support <- npm[[paste0("support_", support_col)]]
  }

  if (include_backbone) {
    bb_name <- "backbone"
    bb_p <- matrix(1, nrow = nrow(presence), ncol = 1,
                   dimnames = list(rownames(presence), bb_name))
    presence <- cbind(bb_p, presence)
    if (mode == "support") {
      bb_s <- node_support(backbone)[[paste0("support_", support_col)]]
      bb_col <- matrix(bb_s[as.integer(rownames(presence)) - ntip],
                       ncol = 1, dimnames = list(rownames(presence), bb_name))
      support <- cbind(bb_col, support)
    }
  }

  tree_names <- colnames(presence)
  n_tree     <- length(tree_names)

  if (is.null(n_cols)) n_cols <- choose_grid(n_tree)$n_cols
  if (is.null(n_rows)) n_rows <- ceiling(n_tree / n_cols)

  unanimous <- apply(presence, 1, function(row) all(row == 1, na.rm = FALSE))

  # ===================================================================
  # TWO-SOURCE SCALING
  #
  # 1. TREE elements scale with per_tip (height-based, log curve).
  #    Bigger trees → denser tips, thinner branches, smaller font.
  #
  # 2. LEGEND elements scale with canvas WIDTH (gentle variation).
  #    Width doesn't change as dramatically as height, so legends
  #    stay readable across all tree sizes.
  # ===================================================================
  dots <- list(...)

  # --- SOURCE 1: per_tip (master curve for tree elements) ---
  per_tip <- min(0.15, max(0.03, 0.31 - 0.041 * log(ntip)))

  # Tree ratios (per_tip × constant)
  R_CEX     <- 4.0                                               # ← RATIO 1
  R_EDGE    <- 6.0                                               # ← RATIO 2
  R_SUPPORT <- 2.3                                               # ← RATIO 3

  if (is.null(dots$cex))
    dots$cex <- min(0.65, max(0.30, per_tip * R_CEX))
  taxa_cex <- dots$cex

  if (is.null(dots$edge.width))
    dots$edge.width <- min(1.3, max(0.7, per_tip * R_EDGE))

  support_cex <- min(0.40, max(0.18, per_tip * R_SUPPORT))

  # --- Margins ---
  if (is.null(dots$mar)) dots$mar <- c(0.5, 0.5, 0.5, 2.5)
  dots$family <- "Helvetica"

  # --- SOURCE 2: canvas width (for legend elements, computed after device) ---
  # Legend ratios (canvas_width × constant). These are applied later
  # in section 10 once we know par("din")[1].
  R_LEG_CELL <- 0.011                                            # ← RATIO 4
  R_TH_SQ    <- 0.008                                            # ← RATIO 5
  R_LEG_TEXT <- 0.045                                             # ← RATIO 6
  R_TH_TEXT  <- 0.041                                             # ← RATIO 7

  # --- Margins ---
  if (is.null(dots$mar)) dots$mar <- c(0.5, 0.5, 0.5, 2.5)
  dots$family <- "Helvetica"

  # --- 3. Device routing -----------------------------------------------------
  if (!is.null(file)) {
    canvas <- auto_canvas(backbone = backbone, ntip = ntip,
                          n_tree = n_tree, mode = mode,
                          has_legend = legend, per_tip = per_tip,
                          taxa_cex = taxa_cex)
    w <- if (is.null(width))  canvas$width  else width
    h <- if (is.null(height)) canvas$height else height

    ext <- tolower(tools::file_ext(file))
    if (ext == "pdf") {
      grDevices::pdf(file, width = w, height = h)
    } else if (ext == "png") {
      grDevices::png(file, width = w, height = h, units = "in", res = 300)
    } else if (ext %in% c("jpg", "jpeg")) {
      grDevices::jpeg(file, width = w, height = h, units = "in", res = 300)
    } else {
      stop("Unsupported file extension. Use .pdf, .png, or .jpg", call. = FALSE)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  # --- 4. Dynamic y.lim ------------------------------------------------------
  if (is.null(dots$y.lim) && legend) {
    din    <- graphics::par("din")
    mai    <- graphics::par("mai")
    plot_h <- din[2] - mai[1] - mai[3]

    # Estimate legend text cex from canvas width
    w_ref        <- din[1] * 0.85 + din[2] * 0.15
    est_pos_cex  <- min(0.55, max(0.20, w_ref * R_LEG_TEXT))
    est_th_cex   <- min(0.50, max(0.18, w_ref * R_TH_TEXT))

    pos_line_h <- 0.2 * est_pos_cex
    th_line_h  <- 0.2 * est_th_cex
    pos_leg_h  <- n_tree * pos_line_h
    th_leg_h   <- if (mode == "support") 6 * th_line_h else 0
    top_leg_h  <- max(pos_leg_h, th_leg_h)
    gap_inches <- 3 * max(pos_line_h, th_line_h)

    avail   <- max(1, plot_h - top_leg_h - gap_inches)
    top_pad <- ceiling((top_leg_h + gap_inches) * ntip / avail)

    dots$y.lim <- c(-1.0, ntip + top_pad)
  }

  # --- 5. Draw tree ----------------------------------------------------------
  do.call(ape::plot.phylo, c(list(x = backbone), dots))
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  # --- 6. Cell geometry ------------------------------------------------------
  yy_tip <- sort(last_pp$yy[seq_len(ntip)])
  dy     <- stats::median(diff(yy_tip))
  pin    <- graphics::par("pin")
  cell_h <- dy * cell_scale
  cell_w <- cell_h * (diff(last_pp$x.lim) / pin[1]) /
    (diff(last_pp$y.lim) / pin[2])

  # --- 7. Node support labels ------------------------------------------------
  if (show_support && !is.null(backbone$node.label)) {
    labs <- backbone$node.label
    show <- which(!is.na(labs) & nzchar(labs))
    if (length(show) > 0) {
      node_ids <- show + ntip
      if (dot_unanimous) {
        keep <- !(node_ids %in% as.integer(rownames(presence)[unanimous]))
        show <- show[keep]
        node_ids <- node_ids[keep]
      }
      if (length(show) > 0) {
        num <- suppressWarnings(as.numeric(labs[show]))
        txt <- ifelse(is.na(num), labs[show],
                      format(round(num, 2), trim = TRUE))
        ape::nodelabels(text = txt, node = node_ids, frame = "none",
                        cex = support_cex, col = "red",
                        adj = c(1.1, 1.4))
      }
    }
  }

  # --- 8. Unanimous dots -----------------------------------------------------
  if (dot_unanimous && any(unanimous)) {
    ids <- as.integer(rownames(presence)[unanimous])
    graphics::points(last_pp$xx[ids], last_pp$yy[ids],
                     pch = 16, cex = dot_cex, col = dot_col)
  }

  # --- 9. Node rugs ----------------------------------------------------------
  variable <- if (dot_unanimous) !unanimous else rep(TRUE, nrow(presence))
  if (any(variable)) {
    plot_node_rug(
      presence     = presence[variable, , drop = FALSE],
      support      = if (is.null(support)) NULL else
        support[variable, , drop = FALSE],
      support_type = support_type,
      thresholds   = thresholds,
      cell_h       = cell_h,
      cell_w       = cell_w,
      n_cols       = n_cols,
      x_offset     = x_offset,
      y_offset     = y_offset,
      rug_position = rug_position,
      last_pp      = last_pp
    )
  }

  # --- 10. Legends -----------------------------------------------------------
  if (legend) {
    usr      <- graphics::par("usr")
    margin_x <- (usr[2] - usr[1]) * 0.015
    margin_y <- (usr[4] - usr[3]) * 0.015

    y0_top  <- usr[4] - margin_y
    x0_left <- usr[1] + margin_x

    # Legend sizes from canvas WIDTH (gentle scaling)
    # w_ref = 85% width + 15% height → width-dominant
    din   <- graphics::par("din")
    w_ref <- din[1] * 0.85 + din[2] * 0.15

    leg_cell_in  <- w_ref * R_LEG_CELL
    th_sq_in     <- w_ref * R_TH_SQ
    pos_text_cex <- min(0.55, max(0.20, w_ref * R_LEG_TEXT))
    th_text_cex  <- min(0.50, max(0.18, w_ref * R_TH_TEXT))

    leg_cell_h <- graphics::yinch(leg_cell_in)
    leg_cell_w <- graphics::xinch(leg_cell_in)
    th_sq_h    <- graphics::yinch(th_sq_in)
    th_sq_w    <- graphics::xinch(th_sq_in)

    # Truncate long names
    max_chars     <- 35L
    display_names <- ifelse(
      nchar(tree_names) > max_chars,
      paste0(substr(tree_names, 1, max_chars - 1), "\u2026"),
      tree_names
    )

    # --- Position legend (topleft) ---
    draw_position_legend(
      display_names, n_cols,
      cell_w   = leg_cell_w,
      cell_h   = leg_cell_h,
      x0       = x0_left,
      y0       = y0_top,
      text_cex = pos_text_cex
    )

    # --- Threshold legend (topright, support mode only) ---
    if (mode == "support") {
      th_width_in <- graphics::strwidth(
        "\u226598 (SH-aLRT/UFBoot2) or \u22650.99 (ASTRAL)",
        units = "inches", cex = th_text_cex
      ) + th_sq_in + 0.1

      x0_right <- usr[2] - margin_x - graphics::xinch(th_width_in)

      draw_threshold_legend(
        x0       = x0_right,
        y0       = y0_top,
        sq_h     = th_sq_h,
        sq_w     = th_sq_w,
        text_cex = th_text_cex
      )
    }
  }

  invisible(if (!is.null(file)) file else NULL)
}
