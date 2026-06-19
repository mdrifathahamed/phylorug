# Tests for choose_grid()

test_that("perfect squares give a square grid", {
  expect_equal(choose_grid(4),  list(n_rows = 2, n_cols = 2))
  expect_equal(choose_grid(9),  list(n_rows = 3, n_cols = 3))
  expect_equal(choose_grid(16), list(n_rows = 4, n_cols = 4))
})

test_that("non-squares round up to fit all cells", {
  # 5 cells: sqrt = 2.24 -> 3 cols, ceiling(5/3) = 2 rows
  expect_equal(choose_grid(5), list(n_rows = 2, n_cols = 3))
  # 8 cells: sqrt = 2.83 -> 3 cols, ceiling(8/3) = 3 rows
  expect_equal(choose_grid(8), list(n_rows = 3, n_cols = 3))
  # 12 cells: sqrt = 3.46 -> 4 cols, ceiling(12/4) = 3 rows
  expect_equal(choose_grid(12), list(n_rows = 3, n_cols = 4))
})

test_that("the grid always has at least as many cells as requested", {
  for (n in 1:30) {
    g <- choose_grid(n)
    expect_gte(g$n_rows * g$n_cols, n)
  }
})

test_that("the grid is never wastefully large", {
  # capacity should not exceed the requested count by a whole extra row
  for (n in 1:30) {
    g <- choose_grid(n)
    expect_lt(g$n_rows * g$n_cols, n + g$n_cols)
  }
})

test_that("single cell gives a 1 by 1 grid", {
  expect_equal(choose_grid(1), list(n_rows = 1, n_cols = 1))
})

test_that("returns a named list with n_rows and n_cols", {
  g <- choose_grid(6)
  expect_type(g, "list")
  expect_named(g, c("n_rows", "n_cols"))
})
