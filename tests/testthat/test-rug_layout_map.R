# Tests for rug_layout_map()

# ---- helper ----------------------------------------------------------------

make_rug_mt <- function(n_nodes = 3, tree_names = c("t1", "t2", "t3")) {
  mt <- matrix(
    c(seq_len(n_nodes), runif(n_nodes * length(tree_names))),
    nrow = n_nodes
  )
  colnames(mt) <- c("node_id", tree_names)
  mt
}

# ---- input validation ------------------------------------------------------

test_that("stops when rug_mt is not a matrix or data frame", {
  expect_error(
    rug_layout_map("not_a_matrix"),
    "must be a matrix or data frame"
  )
})

test_that("stops when rug_mt has fewer than two columns", {
  mt <- matrix(1:3, ncol = 1)
  colnames(mt) <- "node_id"
  expect_error(
    rug_layout_map(mt),
    "at least two columns"
  )
})

test_that("stops when the grid is too small for the trees", {
  mt <- make_rug_mt(2, c("t1", "t2", "t3", "t4"))
  expect_error(
    rug_layout_map(mt, n_rows = 1, n_cols = 2), # 2 cells, 4 trees
    "Increase"
  )
})

# ---- output shape ----------------------------------------------------------

test_that("returns one row per node-and-tree pair", {
  mt <- make_rug_mt(3, c("t1", "t2", "t3"))
  result <- rug_layout_map(mt)
  expect_equal(nrow(result), 3 * 3) # 3 nodes x 3 trees
})

test_that("returns the five expected columns", {
  result <- rug_layout_map(make_rug_mt())
  expect_named(result, c("node_id", "tree_name", "row", "col", "cell_index"))
})

test_that("tree_name is character, not a factor", {
  result <- rug_layout_map(make_rug_mt())
  expect_type(result$tree_name, "character")
})

# ---- grid position maths ---------------------------------------------------

test_that("cells fill left to right, top to bottom for a fixed 3-column grid", {
  mt <- make_rug_mt(1, c("t1", "t2", "t3", "t4", "t5"))
  result <- rug_layout_map(mt, n_rows = 2, n_cols = 3)
  # one node, so cell_index 1..5 map to:
  expect_equal(result$row, c(1, 1, 1, 2, 2))
  expect_equal(result$col, c(1, 2, 3, 1, 2))
  expect_equal(result$cell_index, 1:5)
})

test_that("node_id is carried through correctly", {
  mt <- make_rug_mt(2, c("t1", "t2"))
  result <- rug_layout_map(mt)
  # node ids are 1 and 2, each repeated once per tree (2 trees)
  expect_equal(result$node_id, c(1, 1, 2, 2))
})

test_that("tree_name matches the matrix columns in order", {
  mt <- make_rug_mt(1, c("alpha", "beta", "gamma"))
  result <- rug_layout_map(mt)
  expect_equal(result$tree_name, c("alpha", "beta", "gamma"))
})

# ---- automatic grid sizing -------------------------------------------------

test_that("auto-sizes the grid when no shape is given", {
  # 4 trees -> choose_grid gives 2x2
  mt <- make_rug_mt(1, c("t1", "t2", "t3", "t4"))
  result <- rug_layout_map(mt)
  expect_equal(max(result$col), 2)
  expect_equal(max(result$row), 2)
})

test_that("respects a user-fixed column count and derives the rows", {
  # 6 trees, user fixes 2 columns -> should derive 3 rows
  mt <- make_rug_mt(1, c("t1", "t2", "t3", "t4", "t5", "t6"))
  result <- rug_layout_map(mt, n_cols = 2)
  expect_equal(max(result$col), 2)
  expect_equal(max(result$row), 3)
})

test_that("accepts a data frame as well as a matrix", {
  df <- as.data.frame(make_rug_mt(2, c("t1", "t2")))
  result <- rug_layout_map(df)
  expect_equal(nrow(result), 2 * 2)
})
