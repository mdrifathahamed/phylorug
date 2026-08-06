# utils-internal.R

choose_grid <- function(n_cells) {
  n_cols <- ceiling(sqrt(n_cells))
  n_rows <- ceiling(n_cells / n_cols)
  list(n_rows = n_rows, n_cols = n_cols)
}

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
  invisible(NULL)
}

