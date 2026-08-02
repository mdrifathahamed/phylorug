# gg_phylorug.R (v2)
# Fixes the cell/legend sizing bug: ggplot2 does not have a fixed
# pixel-per-data-unit relationship like base R's par("pin") does, so cell
# sizes must be computed relative to a *known* output size and locked in
# with coord_fixed(). This mirrors the base-R engine's pin/x.lim approach.

#' Draw a phylorug using ggplot2
#'
#' @param backbone A `phylo` object.
#' @param npm Named list from `node_presence_matrix()`.
#' @param mode `"presence"` or `"support"`.
#' @param support_idx Integer. Which support matrix to use. Default `1`.
#' @param support_type Named character vector mapping tree names to support
#'   measures. Required for `mode = "support"`.
#' @param thresholds Optional threshold overrides.
#' @param n_cols Integer or `NULL` for automatic grid shape.
#' @param rug_position `"outside"` (default) or `"inside"`.
#' @param dot_identical Logical. Draw a dot at unanimous nodes. Default `TRUE`.
#' @param tip_size Numeric text size for tip labels. Default `3`.
#' @param show_support Logical. Draw backbone support labels at non-unanimous
#'   nodes. Default `TRUE`.
#' @param support_label_size,support_label_col Styling for support labels.
#' @param legend Logical. Draw the position legend. Default `TRUE`.
#' @param cell_size_in Numeric. Cell edge length in inches. Default `0.09`.
#' @param file Optional output path (`.png`, `.pdf`, `.jpg`). If `NULL`
#'   (default), the `ggplot` object is returned without saving.
#' @param width,height Output dimensions in inches. `NULL` (default) uses
#'   `auto_canvas()` -- the same canvas-sizing engine `plot_phylorug()`
#'   uses -- to compute dimensions from the tree's actual properties
#'   (depth, tip count, label length, legend space). Supplying either
#'   overrides the automatic value for that dimension only.
#' @param dpi Resolution for raster output. Default `200`.
#'
#' @returns A `ggplot` object (invisibly, if `file` was given).
#' @export
gg_phylorug <- function(backbone, npm,
                        mode               = c("presence", "support"),
                        support_idx        = 1,
                        support_type       = NULL,
                        thresholds         = NULL,
                        n_cols             = NULL,
                        rug_position       = c("outside", "inside"),
                        dot_identical      = TRUE,
                        dot_size           = 3.2,
                        tip_size           = 3,
                        show_support       = TRUE,
                        support_label_size = 2.5,
                        support_label_col  = "red",
                        legend             = TRUE,
                        cell_size_in       = 0.09,
                        file               = NULL,
                        width              = NULL,
                        height             = NULL,
                        dpi                = 200) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for gg_phylorug().", call. = FALSE)

  mode         <- match.arg(mode)
  rug_position <- match.arg(rug_position)

  presence <- npm$presence
  if (ncol(presence) < 1L)
    stop("`npm` has no comparison trees to plot.", call. = FALSE)

  support <- NULL
  if (mode == "support") {
    support <- npm[[paste0("support_", support_idx)]]
    if (is.null(support))
      stop(sprintf("`support_%d` not found in `npm`.", support_idx), call. = FALSE)
  }

  tree_names <- colnames(presence)
  n_tree     <- length(tree_names)
  if (is.null(n_cols)) n_cols <- choose_grid(n_tree)$n_cols

  unanimous <- apply(presence, 1, function(row) all(row == 1, na.rm = FALSE))

  # --- Canvas sizing: reuse the SAME auto_canvas() engine plot_phylorug()
  #     uses, rather than a hardcoded default. -------------------------------
  if (is.null(width) || is.null(height)) {
    ntip_est     <- ape::Ntip(backbone)
    per_tip_est  <- min(0.20, max(0.18, 0.35 - 0.041 * log(ntip_est)))
    taxa_cex_est <- min(0.80, max(0.55, per_tip_est * 4.0))
    canvas <- auto_canvas(
      backbone = backbone, ntip = ntip_est, n_tree = n_tree,
      mode = mode, has_legend = legend,
      per_tip = per_tip_est, taxa_cex = taxa_cex_est
    )
    if (is.null(width))  width  <- canvas$width
    if (is.null(height)) height <- canvas$height
  }

  # --- Layout -----------------------------------------------------------------
  lo <- gg_tree_layout(backbone)
  ntip <- lo$ntip

  all_x <- c(lo$tips$x, lo$nodes$x)
  all_y <- c(lo$tips$y, lo$nodes$y)
  x_range <- diff(range(all_x))
  y_range <- diff(range(all_y))
  if (x_range <= 0) x_range <- 1
  if (y_range <= 0) y_range <- 1

  # --- Reserve real space for tip label text, measured the same way the
  #     base-R engine measures the right margin (graphics::strwidth() on a
  #     throwaway device), so labels don't get clipped at the canvas edge.
  #     ggplot2's `size` aesthetic for geom_text is in mm, converted to
  #     points via the constant .pt = 72.27/25.4.
  pt_size   <- tip_size * (72.27 / 25.4)
  tmp_lbl   <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp_lbl)
  graphics::par(family = "sans", ps = 12)
  label_w_in <- max(graphics::strwidth(
    lo$tips$label, units = "inches", cex = pt_size / 12
  ))
  grDevices::dev.off()
  unlink(tmp_lbl)
  label_w_in <- label_w_in * 1.15 + 0.05  # small safety margin

  # --- The core fix: derive data-units-per-inch from the KNOWN output size,
  #     then lock that scaling in with coord_fixed(). coord_fixed()'s ratio
  #     parameter is defined so that panel_height/panel_width = ratio *
  #     (y_range/x_range). To match our target output aspect (height/width),
  #     solve: ratio = (height/width) / (y_range/x_range)
  #                   = (height * x_range) / (width * y_range)
  #     Both the ratio AND the visible xlim use the label-padded width, so
  #     the reserved space is actually inside the panel, not clipped by it.
  units_per_in_x_raw <- x_range / width
  label_pad_units     <- label_w_in * units_per_in_x_raw
  x_max_padded        <- max(all_x) + label_pad_units
  x_range_padded       <- x_max_padded - min(all_x)

  fixed_ratio <- (height * x_range_padded) / (width * y_range)

  units_per_in_x <- x_range_padded / width
  units_per_in_y <- y_range / height

  cell_w <- cell_size_in * units_per_in_x
  cell_h <- cell_size_in * units_per_in_y

  # --- Rug cells at variable (non-unanimous) nodes ---------------------------
  variable_ids <- as.integer(rownames(presence))[!unanimous]
  var_presence <- presence[!unanimous, , drop = FALSE]
  var_support  <- if (is.null(support)) NULL else support[!unanimous, , drop = FALSE]

  cells <- NULL
  if (nrow(var_presence) > 0) {
    cells <- gg_rug_cells(
      presence = var_presence, support = var_support,
      support_type = support_type, thresholds = thresholds,
      node_xy = lo$nodes, cell_h = cell_h, cell_w = cell_w,
      n_cols = n_cols, rug_position = rug_position,
      max_x = max(all_x), max_y = max(all_y)
    )
  }

  # --- Unanimous dots ---------------------------------------------------------
  dot_df <- NULL
  if (dot_identical && any(unanimous)) {
    uni_ids <- as.integer(rownames(presence))[unanimous]
    dot_df  <- lo$nodes[lo$nodes$node %in% uni_ids, ]
  }

  # --- Support labels ----------------------------------------------------------
  label_df <- NULL
  if (show_support && !is.null(backbone$node.label)) {
    labs <- backbone$node.label
    show_idx <- which(!is.na(labs) & nzchar(labs))
    node_ids <- show_idx + ntip
    if (dot_identical) {
      uni_ids  <- as.integer(rownames(presence))[unanimous]
      keep     <- !(node_ids %in% uni_ids)
      show_idx <- show_idx[keep]
      node_ids <- node_ids[keep]
    }
    if (length(show_idx) > 0) {
      num <- suppressWarnings(as.numeric(labs[show_idx]))
      txt <- ifelse(is.na(num), labs[show_idx], format(round(num, 2), trim = TRUE))
      xy  <- lo$nodes[lo$nodes$node %in% node_ids, ]
      xy  <- xy[match(node_ids, xy$node), ]
      label_df <- data.frame(x = xy$x, y = xy$y, label = txt)
    }
  }

  # --- Build the plot ------------------------------------------------------
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = lo$edges_h,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend)
    ) +
    ggplot2::geom_segment(
      data = lo$edges_v,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend)
    ) +
    ggplot2::geom_text(
      data = lo$tips,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0, size = tip_size, fontface = "italic",
      nudge_x = 0.008 * x_range
    )

  if (!is.null(cells)) {
    p <- p +
      ggplot2::geom_rect(
        data = cells$rects,
        ggplot2::aes(
          xmin = .data$xmin, xmax = .data$xmax,
          ymin = .data$ymin, ymax = .data$ymax,
          fill = .data$fill
        ),
        colour = "grey40", linewidth = 0.15
      )

    if (!is.null(cells$cross)) {
      p <- p + ggplot2::geom_segment(
        data = cells$cross,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        colour = "grey30", linewidth = 0.3
      )
    }
    if (!is.null(cells$dots)) {
      p <- p + ggplot2::geom_point(
        data = cells$dots,
        ggplot2::aes(x = .data$x, y = .data$y),
        size = 0.4, colour = "black"
      )
    }
  }

  if (!is.null(dot_df)) {
    p <- p + ggplot2::geom_point(
      data = dot_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      size = dot_size, colour = "black"
    )
  }

  if (!is.null(label_df)) {
    p <- p + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = support_label_size, colour = support_label_col,
      hjust = 1.1, vjust = -0.6
    )
  }

  # --- Legend: cleaner styling -- background panel, bold header, and (for
  #     support mode) a threshold colour key, matching what the base-R
  #     engine's draw_threshold_legend() shows. ------------------------------
  if (legend) {
    leg_x0 <- min(all_x) - 0.02 * x_range
    leg_y0 <- max(all_y) + 0.07 * y_range
    leg_cell_w <- cell_w * 1.6
    leg_cell_h <- cell_h * 1.6
    header_gap <- leg_cell_h * 0.9

    leg_rects <- do.call(rbind, lapply(seq_len(n_tree), function(k) {
      row_idx <- ceiling(k / n_cols)
      col_idx <- ((k - 1) %% n_cols) + 1
      data.frame(
        xmin = leg_x0 + (col_idx - 1) * leg_cell_w,
        xmax = leg_x0 + col_idx * leg_cell_w,
        ymin = leg_y0 - header_gap - row_idx * leg_cell_h,
        ymax = leg_y0 - header_gap - (row_idx - 1) * leg_cell_h,
        label = as.character(k),
        lx = leg_x0 + (col_idx - 0.5) * leg_cell_w,
        ly = leg_y0 - header_gap - (row_idx - 0.5) * leg_cell_h
      )
    }))

    leg_text <- data.frame(
      x = leg_x0 + n_cols * leg_cell_w + 0.02 * x_range,
      y = seq(leg_y0 - header_gap - leg_cell_h / 2, by = -leg_cell_h, length.out = n_tree),
      label = tree_names
    )

    n_leg_rows <- ceiling(n_tree / n_cols)
    grid_w  <- n_cols * leg_cell_w

    tmp_leg <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp_leg)
    graphics::par(family = "sans")
    text_w      <- max(graphics::strwidth(tree_names, units = "inches", cex = 0.5))
    th_labels_w <- max(graphics::strwidth(
      c("Very high", "High", "Moderate", "Low", "Not evaluable"),
      units = "inches", cex = 0.5
    ))
    grDevices::dev.off()
    unlink(tmp_leg)

    text_w_units <- text_w * units_per_in_x
    panel_w <- grid_w + text_w_units + 0.06 * x_range
    panel_h <- header_gap + n_leg_rows * leg_cell_h + 0.03 * y_range

    panel_bg <- data.frame(
      xmin = leg_x0 - 0.02 * x_range,
      xmax = leg_x0 + panel_w,
      ymin = leg_y0 - panel_h,
      ymax = leg_y0 + 0.02 * y_range
    )

    header_df <- data.frame(
      x = leg_x0, y = leg_y0 - header_gap * 0.35, label = "Analyses"
    )

    p <- p +
      ggplot2::geom_rect(
        data = panel_bg,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ymin, ymax = .data$ymax),
        fill = "grey97", colour = "grey65", linewidth = 0.3
      ) +
      ggplot2::geom_text(
        data = header_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        hjust = 0, size = 3.0, fontface = "bold", colour = "grey20"
      ) +
      ggplot2::geom_rect(
        data = leg_rects,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ymin, ymax = .data$ymax),
        fill = "white", colour = "grey40", linewidth = 0.3
      ) +
      ggplot2::geom_text(
        data = leg_rects,
        ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$label),
        size = 2.3, colour = "grey20"
      ) +
      ggplot2::geom_text(
        data = leg_text,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = 2.3, hjust = 0, colour = "grey20"
      )

    # --- Threshold colour key (support mode only), top-right ----------------
    if (mode == "support") {
      th_states <- data.frame(
        tier  = c("Very high", "High", "Moderate", "Low", "Not evaluable"),
        fill  = c("#000000", "#5F5E5A", "#B4B2A9", "white", "white"),
        cross = c(FALSE, FALSE, FALSE, FALSE, TRUE)
      )
      th_w <- leg_cell_w
      th_h <- leg_cell_h
      th_x0 <- x_max_padded - 0.02 * x_range_padded
      th_y0 <- leg_y0

      th_rects <- data.frame(
        xmin = th_x0 - th_w, xmax = th_x0,
        ymin = th_y0 - header_gap - seq_len(5) * th_h,
        ymax = th_y0 - header_gap - (seq_len(5) - 1) * th_h,
        fill = th_states$fill
      )
      th_text <- data.frame(
        x = th_x0 - th_w - 0.015 * x_range,
        y = th_rects$ymin + th_h / 2,
        label = th_states$tier
      )
      th_bg <- data.frame(
        xmin = th_x0 - th_w - th_labels_w * units_per_in_x - 0.05 * x_range,
        xmax = th_x0 + 0.02 * x_range,
        ymin = min(th_rects$ymin) - 0.02 * y_range,
        ymax = th_y0 + 0.02 * y_range
      )

      p <- p +
        ggplot2::geom_rect(
          data = th_bg,
          ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                       ymin = .data$ymin, ymax = .data$ymax),
          fill = "grey97", colour = "grey65", linewidth = 0.3
        ) +
        ggplot2::geom_text(
          data = data.frame(x = th_bg$xmax - 0.01 * x_range, y = th_y0 - header_gap * 0.35,
                            label = "Support"),
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
          hjust = 1, size = 3.0, fontface = "bold", colour = "grey20"
        ) +
        ggplot2::geom_rect(
          data = th_rects,
          ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                       ymin = .data$ymin, ymax = .data$ymax, fill = .data$fill),
          colour = "grey40", linewidth = 0.3
        ) +
        ggplot2::geom_text(
          data = th_text,
          ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
          size = 2.3, hjust = 1, colour = "grey20"
        )
    }
  }

  p <- p +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_fixed(
      ratio = fixed_ratio,
      xlim  = c(min(all_x), x_max_padded),
      clip  = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))

  if (!is.null(file)) {
    ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = dpi)
    return(invisible(p))
  }

  p
}
