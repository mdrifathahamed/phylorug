# Tests for plot_phylorug() and its internal helpers

# ---- helpers ----------------------------------------------------------------
make_backbone <- function() {
  ape::read.tree(text = "(((A,B),C),(D,E));")
}

make_backbone_with_labels <- function() {
  ape::read.tree(text = "(((A,B)95,C)80,(D,E)100);")
}

make_npm <- function() {
  nodes <- as.character(6:9)
  cols  <- c("iqtree", "astral")
  presence <- matrix(
    c(1, 0, 1, 1,  1, 1, 0, 1),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_1 <- matrix(
    c(99, NA, 85, 70,  0.99, 0.60, NA, 0.97),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_2 <- matrix(
    c(95, NA, 80, 65,  NA, NA, NA, NA),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  out <- list(presence = presence, support_1 = support_1, support_2 = support_2)
  attr(out, "node_id")    <- 6:9
  attr(out, "pool_sizes") <- c(iqtree = 1L, astral = 1L)
  out
}

make_npm_all_na_support <- function() {
  nodes <- as.character(6:9)
  cols  <- c("iqtree", "astral")
  presence <- matrix(
    c(1, 0, 1, 1,  1, 1, 0, 1),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_1 <- matrix(
    NA_real_,
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  list(presence = presence, support_1 = support_1)
}

make_npm_partial_na_support <- function() {
  nodes <- as.character(6:9)
  cols  <- c("iqtree", "astral")
  presence <- matrix(
    c(1, 0, 1, 1,  1, 1, 0, 1),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_1 <- matrix(
    c(99, NA, 85, 70,  NA, NA, NA, NA),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  list(presence = presence, support_1 = support_1)
}

on_null_device <- function(code) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  force(code)
}


# ---- input validation -------------------------------------------------------
test_that("stops when backbone is not a phylo object", {
  expect_error(
    plot_phylorug("not_a_tree", make_npm()),
    "phylo"
  )
})

test_that("stops when npm has no presence element", {
  expect_error(
    plot_phylorug(make_backbone(), list(support_1 = matrix(1, 2, 2))),
    "presence"
  )
})

test_that("stops when npm is not a list", {
  expect_error(
    plot_phylorug(make_backbone(), matrix(1, 2, 2)),
    "presence"
  )
})

test_that("stops when mode is invalid", {
  expect_error(
    plot_phylorug(make_backbone(), make_npm(), mode = "invalid"),
    "should be one of"
  )
})

test_that("stops when rug_position is invalid", {
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), rug_position = "invalid")
    ),
    "should be one of"
  )
})

test_that("stops when support mode requested but support_type is NULL", {
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support", support_idx = 1)
    ),
    "requires `support_type`"
  )
})

test_that("stops when support_idx references a missing matrix", {
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support", support_idx = 3,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    ),
    "not found in `npm`"
  )
})

test_that("stops when all support values are NA", {
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm_all_na_support(),
                    mode = "support", support_idx = 1,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    ),
    "no support values"
  )
})

test_that("warns when some comparison trees have all-NA support", {
  expect_warning(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm_partial_na_support(),
                    mode = "support", support_idx = 1,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    ),
    "no support values"
  )
})


# ---- presence mode ----------------------------------------------------------
test_that("presence mode runs without error", {
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), make_npm()))
  )
})

test_that("presence mode returns invisible NULL on device", {
  result <- on_null_device(
    plot_phylorug(make_backbone(), make_npm(), mode = "presence")
  )
  expect_null(result)
})

test_that("presence mode ignores support_type silently", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "presence",
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})


# ---- support mode -----------------------------------------------------------
test_that("support mode runs with support_idx = 1", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support", support_idx = 1,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})

test_that("support mode runs with support_idx = 2", {
  expect_no_error(
    on_null_device(
      suppressWarnings(
        plot_phylorug(make_backbone(), make_npm(),
                      mode = "support", support_idx = 2,
                      support_type = c(iqtree = "sh_alrt", astral = "lpp"))
      )
    )
  )
})

test_that("support mode accepts custom thresholds", {
  custom <- list(
    sh_alrt = c(very_high = 98, high = 80, moderate = 50),
    lpp     = c(very_high = 0.99, high = 0.95, moderate = 0.5)
  )
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support",
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"),
                    thresholds = custom)
    )
  )
})


# ---- include_backbone -------------------------------------------------------
test_that("include_backbone = TRUE adds a backbone column in presence mode", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), include_backbone = TRUE)
    )
  )
})
test_that("include_backbone = TRUE works in support mode", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone_with_labels(), make_npm(),
                    mode = "support", support_idx = 1,
                    include_backbone = TRUE,
                    support_type = c(backbone = "ufboot",
                                     iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})


# ---- show_support -----------------------------------------------------------
test_that("show_support = TRUE draws node labels", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone_with_labels(), make_npm(),
                    show_support = TRUE)
    )
  )
})

test_that("show_support = FALSE hides node labels", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone_with_labels(), make_npm(),
                    show_support = FALSE)
    )
  )
})

test_that("show_support auto-resolves to FALSE when support + include_backbone", { # nolint: line_length_linter.
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone_with_labels(), make_npm(),
                    mode = "support", support_idx = 1,
                    include_backbone = TRUE,
                    support_type = c(backbone = "ufboot",
                                     iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})
test_that("include_backbone + support mode errors without backbone in support_type", { # nolint: line_length_linter.
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone_with_labels(), make_npm(),
                    mode = "support", support_idx = 1,
                    include_backbone = TRUE,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    ),
    "backbone.*entry in `support_type`"
  )
})
test_that("show_support with support_label_cex and support_label_col", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone_with_labels(), make_npm(),
                    show_support = TRUE,
                    support_label_cex = 0.5,
                    support_label_col = "blue")
    )
  )
})


# ---- legend -----------------------------------------------------------------
test_that("runs with legend = TRUE (default)", {
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), make_npm(), legend = TRUE))
  )
})

test_that("runs with legend = FALSE", {
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), make_npm(), legend = FALSE))
  )
})

test_that("support mode legend includes threshold key", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support", support_idx = 1,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"),
                    legend = TRUE)
    )
  )
})


# ---- dot_identical ----------------------------------------------------------
test_that("dot_identical = TRUE draws dots on unanimous nodes", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), dot_identical = TRUE)
    )
  )
})

test_that("dot_identical = FALSE shows rugs on all nodes", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), dot_identical = FALSE)
    )
  )
})

test_that("custom dot_col and dot_cex work", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    dot_col = "red", dot_cex = 1.5)
    )
  )
})


# ---- rug_on_identical -------------------------------------------------------
test_that("rug_on_identical = TRUE draws both dot and rug on unanimous nodes", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    dot_identical = TRUE, rug_on_identical = TRUE)
    )
  )
})


# ---- hide_unsupported -------------------------------------------------------
test_that("hide_unsupported = TRUE skips unsupported nodes", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), hide_unsupported = TRUE)
    )
  )
})

test_that("hide_unsupported + dot_identical combined", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    hide_unsupported = TRUE, dot_identical = TRUE)
    )
  )
})


# ---- grid shape -------------------------------------------------------------
test_that("runs with user-fixed grid shape", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), n_rows = 1, n_cols = 2)
    )
  )
})

test_that("auto grid shape works with default n_rows and n_cols", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), n_rows = NULL, n_cols = NULL)
    )
  )
})


# ---- rug_position -----------------------------------------------------------
test_that("rug_position = 'inside' works", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), rug_position = "inside")
    )
  )
})

test_that("rug_position = 'outside' works", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), rug_position = "outside")
    )
  )
})


# ---- cell_scale and offsets -------------------------------------------------
test_that("cell_scale adjusts rug size", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), cell_scale = 0.7)
    )
  )
})

test_that("x_offset and y_offset shift the rug", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    x_offset = 0.05, y_offset = 0.05)
    )
  )
})


# ---- passthrough to plot.phylo ----------------------------------------------
test_that("... passes through to plot.phylo", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    edge.width = 2, cex = 0.8, font = 3)
    )
  )
})


# ---- file output ------------------------------------------------------------
test_that("file = PDF creates a file and returns the path", {
  tmp <- tempfile(fileext = ".pdf")
  result <- plot_phylorug(make_backbone(), make_npm(), file = tmp)
  expect_true(file.exists(tmp))
  expect_equal(result, tmp)
  unlink(tmp)
})

test_that("file = PNG creates a file", {
  tmp <- tempfile(fileext = ".png")
  result <- plot_phylorug(make_backbone(), make_npm(), file = tmp)
  expect_true(file.exists(tmp))
  expect_equal(result, tmp)
  unlink(tmp)
})

test_that("file = JPG creates a file", {
  tmp <- tempfile(fileext = ".jpg")
  result <- plot_phylorug(make_backbone(), make_npm(), file = tmp)
  expect_true(file.exists(tmp))
  expect_equal(result, tmp)
  unlink(tmp)
})

test_that("unsupported file extension errors", {
  tmp <- tempfile(fileext = ".svg")
  expect_error(
    plot_phylorug(make_backbone(), make_npm(), file = tmp),
    "Unsupported file extension"
  )
})

test_that("custom width and height override auto_canvas", {
  tmp <- tempfile(fileext = ".pdf")
  expect_no_error(
    plot_phylorug(make_backbone(), make_npm(),
                  file = tmp, width = 10, height = 8)
  )
  expect_true(file.exists(tmp))
  unlink(tmp)
})


# ---- edge cases -------------------------------------------------------------
test_that("runs when every node is variable (no unanimous dots)", {
  npm <- make_npm()
  npm$presence[, 1] <- c(1, 0, 1, 0)
  npm$presence[, 2] <- c(0, 1, 0, 1)
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), npm))
  )
})

test_that("runs when every node is unanimous (dots only, no rug)", {
  npm <- make_npm()
  npm$presence[] <- 1
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), npm))
  )
})

test_that("runs with a single comparison tree", {
  nodes <- as.character(6:9)
  npm <- list(
    presence  = matrix(c(1, 0, 1, 1), nrow = 4, ncol = 1,
                       dimnames = list(nodes, "only")),
    support_1 = matrix(c(99, NA, 85, 70), nrow = 4, ncol = 1,
                       dimnames = list(nodes, "only"))
  )
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), npm))
  )
})

test_that("runs with a backbone that has no node labels", {
  bb <- make_backbone()
  bb$node.label <- NULL
  expect_no_error(
    on_null_device(plot_phylorug(bb, make_npm(), show_support = TRUE))
  )
})

test_that("runs with a backbone that has non-numeric node labels", {
  bb <- make_backbone()
  bb$node.label <- c("root", "clade_A", "clade_B", "clade_C")
  expect_no_error(
    on_null_device(plot_phylorug(bb, make_npm(), show_support = TRUE))
  )
})


# ---- internal: choose_grid --------------------------------------------------
test_that("choose_grid returns near-square for small counts", {
  expect_equal(choose_grid(1), list(n_rows = 1L, n_cols = 1L))
  expect_equal(choose_grid(4), list(n_rows = 2L, n_cols = 2L))
  expect_equal(choose_grid(9), list(n_rows = 3L, n_cols = 3L))
})

test_that("choose_grid handles non-square counts", {
  g5 <- choose_grid(5)
  expect_equal(g5$n_cols, 3L)
  expect_equal(g5$n_rows, 2L)
})

test_that("choose_grid handles 1 cell", {
  expect_equal(choose_grid(1), list(n_rows = 1L, n_cols = 1L))
})

test_that("choose_grid handles 2 cells", {
  g2 <- choose_grid(2)
  expect_equal(g2$n_cols, 2L)
  expect_equal(g2$n_rows, 1L)
})


# ---- internal: auto_canvas --------------------------------------------------
test_that("auto_canvas returns width and height", {
  bb     <- make_backbone()
  result <- auto_canvas(bb, ntip = 5L)
  expect_named(result, c("width", "height"))
  expect_true(result$width > 0)
  expect_true(result$height > 0)
})

test_that("auto_canvas height scales with tip count", {
  bb <- make_backbone()
  small <- auto_canvas(bb, ntip = 5L)
  large <- auto_canvas(bb, ntip = 50L)
  expect_true(large$height > small$height)
})

test_that("auto_canvas handles tree without branch lengths", {
  bb <- make_backbone()
  bb$edge.length <- NULL
  result <- auto_canvas(bb, ntip = 5L)
  expect_true(result$width > 0)
})

test_that("auto_canvas width increases with legend", {
  bb <- make_backbone()
  with_leg    <- auto_canvas(bb, ntip = 5L, has_legend = TRUE)
  without_leg <- auto_canvas(bb, ntip = 5L, has_legend = FALSE)
  expect_true(with_leg$width > without_leg$width)
})
test_that("auto_canvas width accounts for long tip labels", {
  bb1 <- make_backbone()
  long_tips <- paste0("Very_Long_Species_Name_", LETTERS[1:5])
  bb2 <- ape::read.tree(text = paste0(
    "(((", long_tips[1], ",", long_tips[2], "),",
    long_tips[3], "),(", long_tips[4], ",",
    long_tips[5], "));"
  ))
  short <- auto_canvas(bb1, ntip = 5L)
  long  <- auto_canvas(bb2, ntip = 5L)
  expect_true(long$width > short$width)
})


# ---- internal: draw_position_legend -----------------------------------------
test_that("draw_position_legend draws without error", {
  on_null_device({
    plot.new()
    plot.window(xlim = c(0, 10), ylim = c(0, 10))
    expect_no_error(
      draw_position_legend(
        analyses = c("iqtree", "astral"),
        n_cols = 2, cell_w = 0.5, cell_h = 0.5,
        x0 = 1, y0 = 9, text_cex = 0.5
      )
    )
  })
})

test_that("draw_position_legend returns invisible y bottom", {
  result <- on_null_device({
    plot.new()
    plot.window(xlim = c(0, 10), ylim = c(0, 10))
    draw_position_legend(
      analyses = c("a", "b", "c"),
      n_cols = 2, cell_w = 0.5, cell_h = 0.5,
      x0 = 1, y0 = 9, text_cex = 0.5
    )
  })
  expect_true(is.numeric(result))
  expect_true(result < 9)  # y moved down
})


# ---- internal: draw_threshold_legend ----------------------------------------
test_that("draw_threshold_legend draws without error", {
  on_null_device({
    plot.new()
    plot.window(xlim = c(0, 10), ylim = c(0, 10))
    expect_no_error(
      draw_threshold_legend(
        x0 = 1, y0 = 9, sq_h = 0.3, sq_w = 0.3, text_cex = 0.4
      )
    )
  })
})

test_that("draw_threshold_legend returns invisible y bottom", {
  result <- on_null_device({
    plot.new()
    plot.window(xlim = c(0, 10), ylim = c(0, 10))
    draw_threshold_legend(
      x0 = 1, y0 = 9, sq_h = 0.3, sq_w = 0.3, text_cex = 0.4
    )
  })
  expect_true(is.numeric(result))
  expect_true(result < 9)
})
