# Tests for utils-internal.R helpers

# ---- choose_grid ------------------------------------------------------------
test_that("choose_grid returns a roughly square grid", {
  g <- choose_grid(4)
  expect_equal(g$n_cols, 2)
  expect_equal(g$n_rows, 2)
})

test_that("choose_grid handles non-square counts", {
  g <- choose_grid(5)
  expect_equal(g$n_cols, 3)
  expect_equal(g$n_rows, 2)
  expect_true(g$n_cols * g$n_rows >= 5)
})

test_that("choose_grid handles 1 cell", {
  g <- choose_grid(1)
  expect_equal(g$n_cols, 1)
  expect_equal(g$n_rows, 1)
})

test_that("choose_grid handles 9 cells", {
  g <- choose_grid(9)
  expect_equal(g$n_cols, 3)
  expect_equal(g$n_rows, 3)
})

# ---- auto_canvas ------------------------------------------------------------
test_that("auto_canvas returns width and height", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(((A,B),C),(D,E));")
  result <- auto_canvas(tree, ntip = 5)
  expect_true(is.list(result))
  expect_true("width" %in% names(result))
  expect_true("height" %in% names(result))
  expect_true(result$width > 0)
  expect_true(result$height > 0)
})

test_that("auto_canvas height scales with ntip", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(((A,B),C),(D,E));")
  small <- auto_canvas(tree, ntip = 20, per_tip = 0.15)
  large <- auto_canvas(tree, ntip = 200, per_tip = 0.10)
  expect_true(large$height > small$height)
})

test_that("auto_canvas width is wider for support mode", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(((A,B),C),(D,E));")
  pres <- auto_canvas(tree, ntip = 50, mode = "presence")
  supp <- auto_canvas(tree, ntip = 50, mode = "support")
  expect_true(supp$width > pres$width)
})

test_that("auto_canvas respects has_legend = FALSE", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(((A,B),C),(D,E));")
  with_leg  <- auto_canvas(tree, ntip = 50, has_legend = TRUE)
  no_leg    <- auto_canvas(tree, ntip = 50, has_legend = FALSE)
  expect_true(with_leg$width > no_leg$width)
})

test_that("auto_canvas uses branch lengths when available", {
  skip_if_not_installed("ape")
  # Tree with branch lengths
  tree_bl <- ape::read.tree(text = "(((A:1,B:1):0.5,C:1.5):0.3,(D:0.8,E:0.8):1);")
  # Same topology, no branch lengths
  tree_no <- ape::read.tree(text = "(((A,B),C),(D,E));")

  c_bl <- auto_canvas(tree_bl, ntip = 5)
  c_no <- auto_canvas(tree_no, ntip = 5)

  # Both should produce valid dimensions (no errors)
  expect_true(c_bl$width > 0)
  expect_true(c_no$width > 0)
})

test_that("auto_canvas width scales with n_tree", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(((A,B),C),(D,E));")
  few  <- auto_canvas(tree, ntip = 50, n_tree = 3)
  many <- auto_canvas(tree, ntip = 50, n_tree = 12)
  expect_true(many$width >= few$width)
})

test_that("auto_canvas enforces minimum dimensions", {
  skip_if_not_installed("ape")
  tree <- ape::read.tree(text = "(A,B);")
  result <- auto_canvas(tree, ntip = 2, per_tip = 0.15)
  expect_true(result$height >= 6)
  expect_true(result$width >= 7)
})

# ---- draw_position_legend ---------------------------------------------------
test_that("draw_position_legend runs without error", {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  plot.new()
  expect_no_error(
    draw_position_legend(
      c("tree_A", "tree_B", "tree_C"), n_cols = 2,
      cell_w = 0.1, cell_h = 0.1, x0 = 0.1, y0 = 0.9, text_cex = 0.5
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

# ---- draw_threshold_legend --------------------------------------------------
test_that("draw_threshold_legend runs without error", {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  plot.new()
  expect_no_error(
    draw_threshold_legend(
      x0 = 0.5, y0 = 0.9,
      sq_h = 0.05, sq_w = 0.05, text_cex = 0.4
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})
