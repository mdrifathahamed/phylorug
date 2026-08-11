#' Draw a phylorug: a backbone tree with node rugs
#'
#' [plot_phylorug()] overlays clade stability grids (rugs) on a backbone
#' phylogeny, comparing how multiple analyses treat each internal node.
#' In presence mode, cells are black (recovered) or white (absent). In
#' support mode, cells are shaded by binned support strength. The function
#' handles canvas sizing, legend placement, and font scaling automatically.
#'
#' @param backbone A `phylo` object representing the backbone tree.
#'
#' @param npm A named list containing the `presence` and `support` matrices,
#'  exactly as returned by [node_presence_matrix()].
#'
#' @param file Optional character string specifying the output file name or full
#'  directory path. If left as `NULL` (the default), the plot renders in the
#'  active graphics device (e.g., RStudio) for quick drafting. To avoid aspect
#'  ratio distortion caused by GUI exports and to generate perfectly scaled,
#'  publication-ready figures, provide a file path here (must end in `.pdf`,
#'  `.png`, or `.jpg`) to utilize the package's internal scaling engine.
#'
#' @param width,height Numeric. Optional canvas dimensions in inches, applied
#'  only when exporting to a `file`.If left as `NULL` (the default), the
#'  package's internal engine dynamically calculates the optimal canvas
#'  dimensions based on the tree size and legend layout. Providing values here
#'  overrides the automatic scaling, which is useful for meeting strict journal
#'  dimension requirements.
#'
#' @param mode One of `"presence"` (default) or `"support"`.
#'
#' @param support_idx Integer (1, 2, or 3). Default is `1`. Specifies which
#'  single support matrix from the `npm` list to visualize when
#'  `mode = "support"`. This corresponds directly to your extraction order in
#'  [node_presence_matrix()]. For example, if you generated the data using
#'  `support_col = c(1, 2)`, passing `2` here tells the plotting engine to
#'  physically shade the grid cells using the second metric (stored in your
#'  list as `support_2`).
#'
#' @param support_type Named character vector mapping comparison trees to their
#'   support metrics (e.g., `"ufboot"`, `"sh_alrt"`, `"lpp"`). Required only
#'   when `mode = "support"`. Names must match your comparison tree names.
#'
#' @param thresholds Optional list overriding built-in bin thresholds. While
#'  default thresholds are provided based on common literature, it is highly
#'  recommended to define your own custom thresholds to suit your specific
#'  analytical framework.
#'
#' @param n_rows,n_cols Integer. Grid shape for the rug at each node. If
#'   `NULL` (the default), a roughly square grid is chosen automatically.
#'
#' @param include_backbone Logical. If `TRUE`, the backbone tree occupies
#'   cell 1 of every rug, showing its own presence (always `1`) or its own
#'   support value. Default `FALSE`, since the backbone defines the topology
#'   and trivially recovers every clade. Useful when the figure will be
#'   edited in a vector editor and every analysis must appear in the grid.
#'
#' @param legend Logical. Draw the legend. Default TRUE.
#'
#' @param show_support Logical. If `TRUE`, backbone node support labels are
#'   drawn beside each node. Default `FALSE`.
#'
#' @param cell_scale Numeric multiplier on cell height. Default 0.45.
#'
#' @param x_offset,y_offset Numeric grid shift. Default 0.
#'
#' @param rug_position One of `"inside"` (default) or `"outside"`. Controls
#'  where the node rug grid is placed relative to the backbone node: `"inside"`
#'  tucks the grid into the crook above-left (toward the root), while
#'  `"outside"` places the grid to the right of the node (toward the tips).
#'
#' @param dot_identical Logical. Default `TRUE`. If `TRUE`, draws a small dot on
#'  backbone nodes where every comparison tree is identical in recovering the
#'  clade.
#'
#' @param dot_col,dot_cex Colour and size of the identical-clade dot. Defaults
#'  are `"black"` and NULL (auto-scales).
#'
#' @param ... Additional arguments passed to [ape::plot.phylo()], such as
#'   `cex`, `edge.width`, `font`, or `label.offset`. These override the
#'   automatic scaling when provided.
#'
#' @param support_label_cex Numeric or `NULL`. Size of the backbone support
#'   labels. Default `NULL` auto-scales with tree size.
#'
#' @param support_label_col Colour of the backbone support labels. Default
#'   `"red"`.
#' @param rug_on_identical Logical. Default `FALSE`. When  it and
#'   `dot_identical`both are  `TRUE`, unanimous nodes receive both the dot and
#'   a support rug, making per-tree support strength visible even at universally
#'   recovered clades.
#' @param hide_unsupported Logical. Default `FALSE`. When `TRUE`, nodes where
#'   no comparison tree recovers the clade are left bare, no rug is drawn.
#'   The absence of both a dot and a rug signals that the clade is unique to
#'   the backbone topology.
#'
#' @returns Invisibly, the file path if a file was written, or `NULL` if plotted
#'   directly to the active graphics device (not recommended).
#'
#' @seealso [node_presence_matrix()] to build the input data,
#'   [check_taxa()] to verify taxon sets, and
#'   `plot_node_rug()` which handles the cell-level
#'   drawing(not used by the user).
#'
#' @export
#'
#' @examples
#' # Build the plotting input from real trees:
#' backbone <- sample_trees[["70p_uce"]]
#' others   <- sample_trees[names(sample_trees) != "70p_uce"]
#' npm <- node_presence_matrix(backbone, others, support_col = 1)
#'
#' # --- Presence mode --------------------------------------------------------
#' # Each cell shows whether an analysis recovered the backbone clade.
#' # Writing to a file uses the internal scaling engine for a clean figure.
#' tmp <- tempfile(fileext = ".pdf")
#' plot_phylorug(backbone, npm, file = tmp)
#' unlink(tmp)
#'
#' # --- Support mode ---------------------------------------------------------
#' # Cells are shaded by support strength. `support_type` tells plot_phylorug
#' # how to read each tree's values: ASTRAL trees carry local posterior
#' # probability ("lpp"), IQ-TREE trees carry UFBoot2 ("ufboot").
#' support_type <- c(
#'   "70p_ASTRAL_partition_entropy" = "lpp",
#'   "70p_ASTRAL_uce"               = "lpp",
#'   "70p_ghost"                    = "ufboot",
#'   "70p_partition_entropy"        = "ufboot"
#' )
#' tmp2 <- tempfile(fileext = ".pdf")
#' plot_phylorug(backbone, npm,
#'               file         = tmp2,
#'               mode         = "support",
#'               support_idx  = 1,
#'               support_type = support_type)
#' unlink(tmp2)
#'
#' # --- Some optional controls -----------------------------------------------
#' # include_backbone = TRUE adds the backbone as its own cell in every rug;
#' # rug_position = "outside" places grids toward the tips instead of the crook;
#' # hide_unsupported = TRUE leaves clades no analysis recovered bare.
#' tmp3 <- tempfile(fileext = ".pdf")
#' plot_phylorug(backbone, npm,
#'               file             = tmp3,
#'               include_backbone = TRUE,
#'               rug_position     = "outside")
#' unlink(tmp3)
plot_phylorug <- function(backbone, npm,
                          file             = NULL,
                          width            = NULL,
                          height           = NULL,
                          mode             = c("presence", "support"),
                          support_idx      = 1,
                          support_type     = NULL,
                          thresholds       = NULL,
                          n_rows           = NULL,
                          n_cols           = NULL,
                          include_backbone = FALSE,
                          legend           = TRUE,
                          show_support     = FALSE,
                          cell_scale       = 0.45,
                          x_offset         = 0,
                          y_offset         = 0,
                          rug_position     = c("inside", "outside"),
                          dot_identical    = TRUE,
                          dot_col          = "black",
                          dot_cex           = NULL,
                          support_label_cex = NULL,
                          support_label_col = "red",
                          rug_on_identical = FALSE,
                          hide_unsupported = FALSE,
                          ...) {

  # --- 1. Validate -----------------------------------------------------------
  if (!inherits(backbone, "phylo"))
    stop("`backbone` must be a phylo object.", call. = FALSE)
  if (!is.list(npm) || !("presence" %in% names(npm)))
    stop("`npm` must be a valid node presence matrix.", call. = FALSE)

  mode         <- match.arg(mode)
  rug_position <- match.arg(rug_position)

  # --- Support matrix validation ---
  if (mode == "support") {
    support_key <- paste0("support_", support_idx)
    if (!(support_key %in% names(npm))) {
      available_supports <- grep("^support_", names(npm), value = TRUE)

      msg <- sprintf(
        paste0(
          "Support matrix `%s` (`support_idx = %d`) not found in `npm`.\n",
          "Available support matrices: %s"
        ),
        support_key,
        support_idx,
        if (length(available_supports) > 0) {
          paste(available_supports, collapse = ", ")
        } else {
          "none"
        }
      )

      stop(msg, call. = FALSE)
    }
  }

  if (mode == "support" && is.null(support_type)) {
    stop(
      "`mode = \"support\"` requires `support_type` to interpret support ",
      "values. Without it, a value like 95 cannot be told apart between ",
      "UFBoot2, SH-aLRT, or a posterior probability. Supply a named vector ",
      "mapping each tree to its measure, e.g. c(iqtree = \"ufboot\").",
      call. = FALSE
    )
  }
  # ---------------------------------------------
  presence <- npm$presence
  if (ncol(presence) < 1L)
    stop("`npm` has no comparison trees to plot.", call. = FALSE)
  if (is.null(show_support))
    show_support <- !(mode == "support" && include_backbone)

  ntip <- ape::Ntip(backbone)

  # --- 2. Build matrices -----------------------------------------------------
  support <- NULL
  if (mode == "support") {
    support <- npm[[paste0("support_", support_idx)]]

    # ---- NEW: guard against empty support ----
    if (all(is.na(support))) {
      stop(
        "`mode = \"support\"` was requested, but the selected support matrix ",
        "(`support_", support_idx, "`) contains no support values (all NA). ",
        "This happens when the trees carry no node support, for example a ",
        "topology-only tree, or a BEAST/TreeAnnotator tree whose bracket ",
        "annotations ape could not import. Use `mode = \"presence\"` for these ",
        "trees, or supply trees whose node labels hold support values.",
        call. = FALSE
      )
    }

    empty_cols <- apply(support, 2, function(col) all(is.na(col)))
    if (any(empty_cols)) {
      warning(
        "These comparison tree(s) have no support values in `support_",
        support_idx, "` (their cells will show as not-recovered or ",
        "not-computed): ",
        paste(colnames(support)[empty_cols], collapse = ", "),
        ". Check that `support_type` maps them correctly, or that the trees ",
        "carry node support.",
        call. = FALSE
      )
    }
  }
  if (include_backbone) {
    bb_name <- "backbone"
    bb_p <- matrix(1, nrow = nrow(presence), ncol = 1,
                   dimnames = list(rownames(presence), bb_name))
    presence <- cbind(bb_p, presence)
    if (mode == "support") {
      if (!(bb_name %in% names(support_type))) {
        stop(
          "`include_backbone = TRUE` with `mode = \"support\"` requires a ",
          "\"backbone\" entry in `support_type`, e.g. ",
          "support_type = c(backbone = \"ufboot\", iqtree = \"sh_alrt\", ...)",
          call. = FALSE
        )
      }
      bb_s <- node_support(backbone)[[paste0("support_", support_idx)]]
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
  unsupported <- apply(presence, 1, function(row) all(row == 0, na.rm = TRUE))

  # ===================================================================
  # TWO-SOURCE SCALING
  #
  # 1. TREE elements scale with per_tip (height-based, log curve).
  #    Bigger trees means denser tips, thinner branches, smaller font.
  #
  # 2. LEGEND elements scale with canvas WIDTH (gentle variation).
  #    Width doesn't change as dramatically as height, so legends
  #    stay readable across all tree sizes.
  # ===================================================================
  dots <- list(...)

  # --- SOURCE 1: per_tip (master curve for tree elements) ---
  per_tip <- min(0.20, max(0.18, 0.35 - 0.041 * log(ntip)))

  # Tree ratios (per_tip * constant)
  R_CEX     <- 4.0
  R_EDGE    <- 10.0
  R_SUPPORT <- 3.0

  if (is.null(dots$cex))
    dots$cex <- min(0.80, max(0.55, per_tip * R_CEX))
  taxa_cex <- dots$cex

  if (is.null(dots$edge.width))
    dots$edge.width <- min(2.0, max(1.5, per_tip * R_EDGE))
 #--------------support_cex------------------------------------------------
  support_cex <- if (is.null(support_label_cex)) {
    min(0.60, max(0.40, per_tip * R_SUPPORT))
  } else {
    support_label_cex
  }
 #----------------------dot_cex-------------------------------------------
   dot_scale <- if (is.null(dot_cex)){
    min(1.2, max(0.80, per_tip * 4.5))
  } else {
    dot_cex
  }

  # --- SOURCE 2: canvas width (for legend elements, computed after device) ---
  # Legend ratios (canvas_width * constant). These are applied later
  # in section 10 once we know par("din")[1].
  R_LEG_CELL <- 0.018                                            #  RATIO 4
  R_TH_SQ    <- 0.008                                            #  RATIO 5
  R_LEG_TEXT <- 0.045                                             #  RATIO 6
  R_TH_TEXT  <- 0.041                                             #  RATIO 7

  # --- Margins ---
  if (is.null(dots$mar)) dots$mar <- c(0.5, 0.5, 0.5, 2.5)
  dots$family <- "sans"

  # --- 3. Device routing & Active Device Mirroring --------------------------
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
  } else {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    actual_w <- graphics::par("din")[1]
    shrink   <- min(1.0, actual_w / 12.0)

    dots$cex        <- taxa_cex * shrink
    taxa_cex        <- dots$cex
    support_cex     <- support_cex * shrink
    dot_scale       <- dot_scale * shrink
    dots$edge.width <- (dots$edge.width %||% 1.5) * max(0.5, shrink)

    if (is.null(dots$label.offset))
      dots$label.offset <- 0.001 * shrink

    max_lab  <- max(nchar(backbone$tip.label))
    r_mar    <- max(3.0, min(7.0, max_lab * 0.08))
    dots$mar <- c(0.5, 0.5, 0.5, r_mar)
  }
  # --- 4. Dynamic y.lim ------------------------------------------------------
  if (is.null(dots$y.lim) && legend) {
    din    <- graphics::par("din")
    mai    <- graphics::par("mai")
    plot_h <- din[2] - mai[1] - mai[3]

    # Estimate legend text cex from canvas width
    w_ref <- din[1]
    est_pos_cex  <- min(0.55, max(0.20, w_ref * R_LEG_TEXT))
    est_th_cex   <- min(0.50, max(0.18, w_ref * R_TH_TEXT))

    pos_line_h <- 0.2 * est_pos_cex
    th_line_h  <- 0.2 * est_th_cex
    pos_leg_h  <- n_tree * pos_line_h
    th_leg_h   <- th_leg_h   <- 6 * th_line_h
    top_leg_h  <- max(pos_leg_h, th_leg_h)
    gap_inches <- 4 * max(pos_line_h, th_line_h)

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
      if (dot_identical && !rug_on_identical) {
        keep <- !(node_ids %in% as.integer(rownames(presence)[unanimous]))
        show <- show[keep]
        node_ids <- node_ids[keep]
      }
      if (length(show) > 0) {
        num <- suppressWarnings(as.numeric(labs[show]))
        txt <- ifelse(is.na(num), labs[show],
                      format(round(num, 2), trim = TRUE))
        ape::nodelabels(text = txt, node = node_ids, frame = "none",
                        cex = support_cex, col = support_label_col,
                        adj = c(1.1, 1.4))
      }
    }
  }

  # --- 8. Unanimous dots (Robust scaling for large trees) -------------------
  if (dot_identical && any(unanimous)) {
    # Match the rows of presence that are unanimous
    uni_rows <- which(unanimous)

    # Map directly via the row index to the plotted internal node coordinates
    node_ids <- as.integer(rownames(presence)[uni_rows])

    # Ensure indices fall within bounds of the current plot layout
    valid_idx <- which(node_ids > ntip & node_ids <= length(last_pp$xx))

    if (length(valid_idx) > 0) {
      plot_nodes <- node_ids[valid_idx]
      graphics::points(
        last_pp$xx[plot_nodes],
        last_pp$yy[plot_nodes],
        pch = 16,
        cex = dot_scale,
        col = dot_col
      )
    }
  }
  # --- 9. Node rugs ----------------------------------------------------------
  variable <- if (dot_identical && !rug_on_identical) {
    !unanimous
  } else {
    rep(TRUE, nrow(presence))
  }

  if (hide_unsupported) {
    variable <- variable & !unsupported
  }
  if (any(variable)) {
    plot_node_rug(
      npm          = presence[variable, , drop = FALSE],
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

    # Legend sizes from canvas WIDTH only (independent of tree height)
    din   <- graphics::par("din")
    w_ref <- din[1]

    leg_cell_in  <- w_ref * R_LEG_CELL
    th_sq_in     <- w_ref * R_TH_SQ
    pos_text_cex <- min(0.55, max(0.20, w_ref * R_LEG_TEXT))
    th_text_cex  <- min(0.50, max(0.18, w_ref * R_TH_TEXT))

    # --- Tree-size adjustment ---
    # On large trees, shrink the legend squares and enlarge the legend text
    # so they stay proportionate. Anchored at ntip = 50 (small trees unchanged).
    size_factor <- min(1.0,  max(1.00, 1 - (ntip - 50) / 1000))
    text_factor <- min(1.15, max(0.85, 1 + (ntip - 50) / 1500))

    leg_cell_in  <- leg_cell_in  * size_factor
    th_sq_in     <- th_sq_in     * size_factor
    pos_text_cex <- pos_text_cex * text_factor
    th_text_cex  <- th_text_cex  * text_factor
    # --- end adjustment ---

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
      longest_label <- paste0(
        ">=80-97.9 (SH-aLRT) or >=95-97 (UFBoot2) ",
        "or >=0.95-0.98 (ASTRAL)"
      )
      th_width_in <- graphics::strwidth(
        longest_label, units = "inches", cex = th_text_cex
      ) + th_sq_in + 0.2

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
# auto_canvas -----------------------------------------------------------------
#' Compute canvas dimensions from actual tree properties
#'
#' Height is set by per_tip. Width is computed from three real measurements:
#'
#' 1. TREE DEPTH - if branch lengths exist, the max root-to-tip distance
#'    determines how much horizontal room the phylogram needs. Deeper trees
#'    get more space (log-scaled, so it doesn't explode). If no branch
#'    lengths, falls back to a cladogram estimate from tip count.
#'
#' 2. LABEL WIDTH - longest tip label * character width at taxa_cex.
#'
#' 3. LEGEND BAND - scales with n_tree and mode.
#'
#' This prevents unnecessary stretching when branches are short, and ensures
#' deep phylograms get enough room to show branch-length variation.
#'
#' @param backbone The backbone tree.
#' @param ntip Integer.
#' @param n_tree Integer.
#' @param mode Character.
#' @param has_legend Logical.
#' @param per_tip Numeric.
#' @param taxa_cex Numeric.
#' @noRd
auto_canvas <- function(backbone,
                        ntip,
                        n_tree     = 4L,
                        mode       = "presence",
                        has_legend = TRUE,
                        per_tip    = 0.15,
                        taxa_cex   = 0.6) {

  # --- Height ---
  height <- max(6, ntip * per_tip)

  # --- Width component 1: Tree depth (inches) ---
  has_bl <- !is.null(backbone$edge.length) &&
    length(backbone$edge.length) > 0 &&
    any(backbone$edge.length > 0, na.rm = TRUE)

  if (has_bl) {
    # Phylogram: use actual root-to-tip distance.
    # node.depth.edgelength() returns distances from root for each node/tip.
    # Max value = deepest root-to-tip path = total tree span on x-axis.
    depths    <- ape::node.depth.edgelength(backbone)
    max_depth <- max(depths)

    # Convert tree depth to inches via log scale.
    # Log prevents tiny-depth trees from being too narrow and deep trees

    tree_width <- max(2, min(8, 2 * log2(1 + max_depth * 50)))
  } else {
    # Cladogram: no branch lengths. Width from tip count (log-scaled).
    tree_width <- max(2, min(6, 1.5 * log2(ntip)))
  }

  # --- Width component 2: Label width (inches) ---
  max_nchar   <- max(nchar(backbone$tip.label))
  label_width <- max_nchar * 0.065 * taxa_cex

  # --- Width component 3: Rug space ---
  rug_width <- 0.5

  # --- Width component 4: Legend band ---
  leg_width <- 0
  if (has_legend) {
    leg_width <- max(2.5, n_tree * 0.25 + 1.5) + 3.5
  }

  # --- Total width ---
  width <- max(7, tree_width + label_width + rug_width + leg_width)

  list(width = width, height = height)
}

#' Choose a near-square grid shape for the per-node rug
#'
#' Decides how many rows and columns the rug grid at each node should have, so
#' that one cell per comparison tree forms a compact, near-square block. Given
#' `n_cells`, it takes `ceiling(sqrt(n_cells))` columns and then the rows needed
#' to hold them all: 4 trees give a 2x2 grid, 5 give 3 columns by 2 rows, 9 give
#' 3x3.
#'
#' @param n_cells Integer. Number of cells to lay out (one per comparison tree).
#'
#' @return A list with `n_rows` and `n_cols` (integers).
#'
#' @seealso [plot_phylorug()], which calls this to lay out each node's rug.
#'
#' @noRd
choose_grid <- function(n_cells) {
  n_cols <- ceiling(sqrt(n_cells))
  n_rows <- ceiling(n_cells / n_cols)
  list(n_rows = n_rows, n_cols = n_cols)
}
#' Draw the position legend (numbered analysis key)
#'
#' Draws the top-left legend of [plot_phylorug()]: a small numbered grid showing
#' which cell position maps to which analysis, followed by a key listing
#' "1 - analysis_name", "2 - analysis_name", and so on. Called once per plot
#' when `legend = TRUE`.
#'
#' @details
#' Cells are laid out left-to-right, top-to-bottom. For cell `k`, `row_idx` and
#' `col_idx` give its position, `xleft`/`xright` span one `cell_w` from `x0`, and
#' `ytop`/`ybottom` drop one `cell_h` per row downward from `y0` (y decreases
#' going down in user coordinates). Each cell is a white box with its number
#' centred inside. The text key is drawn to the right of the grid, starting half
#' a cell past the grid's right edge and top-aligned to `y0`.
#'
#' All coordinates are in user units; the caller ([plot_phylorug()]) converts
#' inches to user units before passing `cell_w`, `cell_h`, `x0`, and `y0`.
#'
#' @param analyses Character vector of comparison-tree names (already truncated
#'   to a maximum length by the caller).
#' @param n_cols Integer. Columns in the mini-grid, from `choose_grid()`.
#' @param cell_w,cell_h Numeric. Width and height of one legend cell, in user
#'   coordinates.
#' @param x0,y0 Numeric. Top-left anchor of the grid, in user coordinates.
#' @param text_cex Numeric. Font size for the cell numbers and the key text.
#'
#' @return Invisibly, the y-coordinate of the bottom of the grid, so the caller
#'   can stack content below it if needed.
#'
#' @seealso [plot_phylorug()] (Section 10) for the caller, and
#'   `draw_threshold_legend()` for the companion support-colour key.
#'
#' @noRd
draw_position_legend <- function(analyses, n_cols, cell_w, cell_h,
                                 x0, y0, text_cex = 0.5) {
  n_an   <- length(analyses)
  n_rows <- ceiling(n_an / n_cols)

  for (k in seq_len(n_an)) {
    row_idx <- ceiling(k / n_cols)
    col_idx <- ((k - 1) %% n_cols) + 1
    xleft   <- x0 + (col_idx - 1) * cell_w
    xright  <- xleft + cell_w
    ytop    <- y0 - (row_idx - 1) * cell_h
    ybottom <- ytop - cell_h
    graphics::rect(xleft, ybottom, xright, ytop,
                   col = "white", border = "black", lwd = 0.5)
    graphics::text((xleft + xright) / 2, (ytop + ybottom) / 2,
                   labels = k, cex = text_cex)
  }

  grid_right <- x0 + n_cols * cell_w
  key_x      <- grid_right + cell_w * 0.5
  key_lines  <- paste0(seq_len(n_an), " \u2013 ", analyses)
  graphics::text(key_x, y0,
                 labels = paste(key_lines, collapse = "\n"),
                 adj = c(0, 1), cex = text_cex, family = "sans")
  invisible(y0 - n_rows * cell_h)
}
#' Draw the support-threshold colour key
#'
#' Draws the top-right legend of [plot_phylorug()] in support mode: the colour
#' key explaining what each cell fill means. Black is very high support, greys
#' are high and moderate, yellow is low, white means the clade was not recovered
#' (monophyly not supported), and red means the clade was recovered but carries
#' no support value (not computed). Called only when `mode = "support"`.
#'
#' @details
#' The colour scheme here must stay in sync with `resolve_cell()` in
#' `plot_node_rug.R`. If a fill colour changes in one place it must change in
#' the other, or the legend will misdescribe the cells. The current scheme is
#' `#000000` very high, `#5F5E5A` high, `#B4B2A9` moderate, `#E8C547` low, white
#' not recovered, and `#D64545` not computed.
#'
#' Rows stack downward from `y0`: row `i` sits `(i - 1) * (sq_h + gap)` below the
#' top, where `gap` is 40 percent of a square's height. Each row is a filled
#' square on the left plus its label, vertically centred, `text_gap` (30 percent
#' of a square's width) to its right. All coordinates are in user units; the
#' caller converts inches to user units before passing `x0`, `y0`, `sq_h`, and
#' `sq_w`.
#'
#' @param x0,y0 Numeric. Top-left anchor of the legend, in user coordinates. The
#'   caller computes `x0` so the squares plus the longest label fit against the
#'   right edge of the plot.
#' @param sq_h,sq_w Numeric. Height and width of one colour square, in user
#'   coordinates.
#' @param text_cex Numeric. Font size for the labels.
#'
#' @return Invisibly, the y-coordinate below the last row.
#'
#' @seealso [plot_phylorug()] (Section 10) for the caller, and
#'   `draw_position_legend()` for the companion numbered analysis key.
#'
#' @noRd
draw_threshold_legend <- function(x0, y0, sq_h, sq_w, text_cex = 0.5) {
gap      <- sq_h * 0.4
text_gap <- sq_w * 0.3

rows <- list(
  list(fill = "#000000",
       label = ">=98 (SH-aLRT/UFBoot2) or >=0.99 (ASTRAL)"),
  list(fill = "#5F5E5A",
       label = ">=80-97.9 (SH-aLRT) or >=95-97 (UFBoot2) or >=0.95-0.98 (ASTRAL)"),
  list(fill = "#B4B2A9",
       label = ">=50-79.9 (SH-aLRT) or >=50-94 (UFBoot2) or >=0.5-0.95 (ASTRAL)"),
  list(fill = "#E8C547",
       label = "<50 (SH-aLRT/UFBoot2) or <0.5 (ASTRAL)"),
  list(fill = "white",
       label = "monophyly not supported"),
  list(fill = "#D64545",
       label = "not computed")
)

# Draw each row: a filled square on the left, its label to the right.
# Rows stack downward from (x0, y0), one square plus a gap per row.
for (i in seq_along(rows)) {
  r       <- rows[[i]]
  ytop    <- y0 - (i - 1) * (sq_h + gap)
  ybottom <- ytop - sq_h
  xright  <- x0 + sq_w

  graphics::rect(x0, ybottom, xright, ytop,
                 col = r$fill, border = "black", lwd = 0.5)
  graphics::text(xright + text_gap, (ytop + ybottom) / 2,
                 labels = r$label, adj = c(0, 0.5),
                 cex = text_cex, family = "sans")
}

invisible(y0 - length(rows) * (sq_h + gap))
}
