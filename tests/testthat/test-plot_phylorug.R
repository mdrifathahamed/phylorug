# Tests for plot_phylorug()

# ---- helpers ----------------------------------------------------------------

make_test_tree <- function() {
  set.seed(1)
  ape::rtree(6)
}

# A support matrix with mixed scales and some NA, shaped like
# node_presence_matrix() output: node_id column then one column per analysis.
make_support_mt <- function(tree) {
  n_tip  <- ape::Ntip(tree)
  n_node <- ape::Nnode(tree)
  node_ids <- (n_tip + 1):(n_tip + n_node)
  a <- c(95, 80, NA, 100, 70)[seq_len(n_node)]
  b <- c(1.0, 0.6, 0.9, NA, 0.8)[seq_len(n_node)]
  cbind(node_id = node_ids, TreeA = a, TreeB = b)
}

# A presence matrix: every present cell is 1, absent is NA.
make_presence_mt <- function(tree) {
  n_tip  <- ape::Ntip(tree)
  n_node <- ape::Nnode(tree)
  node_ids <- (n_tip + 1):(n_tip + n_node)
  a <- c(1, 1, NA, 1, 1)[seq_len(n_node)]
  b <- c(1, NA, 1, 1, NA)[seq_len(n_node)]
  cbind(node_id = node_ids, TreeA = a, TreeB = b)
}

# Every analysis present at every node: all nodes are unanimous (dots only).
make_all_present_mt <- function(tree) {
  n_tip  <- ape::Ntip(tree)
  n_node <- ape::Nnode(tree)
  node_ids <- (n_tip + 1):(n_tip + n_node)
  cbind(
    node_id = node_ids,
    TreeA   = rep(1, n_node),
    TreeB   = rep(1, n_node)
  )
}

# At least one NA in every row: every node is variable (rug only, no dots).
make_all_variable_mt <- function(tree) {
  n_tip  <- ape::Ntip(tree)
  n_node <- ape::Nnode(tree)
  node_ids <- (n_tip + 1):(n_tip + n_node)
  cbind(
    node_id = node_ids,
    TreeA   = rep(NA_real_, n_node),
    TreeB   = rep(1, n_node)
  )
}

# Run code with a tree-capable null device active, so plotting draws nothing.
on_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(code)
}

# ---- input validation -------------------------------------------------------

test_that("stops when backbone is not a phylo object", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_error(
    on_null_device(plot_phylorug("not_a_tree", rug_mt)),
    "must be a phylogenetic tree"
  )
})

test_that("stops when rug_mt is not a matrix or data frame", {
  tree <- make_test_tree()
  expect_error(
    on_null_device(plot_phylorug(tree, "not_a_matrix")),
    "must be a matrix or data frame"
  )
})

test_that("stops when rug_mt has fewer than two columns", {
  tree <- make_test_tree()
  mt   <- matrix(1:3, ncol = 1, dimnames = list(NULL, "node_id"))
  expect_error(
    on_null_device(plot_phylorug(tree, mt)),
    "node_id column and at least one analysis"
  )
})

test_that("stops when hues count does not match the analyses", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_error(
    on_null_device(
      plot_phylorug(tree, rug_mt, colour = TRUE, hues = c("red", "blue", "green"))
    ),
    "one colour per analysis"
  )
})

test_that("stops when a user grid is too small for the analyses", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)   # 2 analyses
  expect_error(
    on_null_device(plot_phylorug(tree, rug_mt, n_rows = 1, n_cols = 1)),
    "cannot hold"
  )
})

# ---- return value -----------------------------------------------------------

test_that("returns invisible NULL", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  on_null_device({
    result <- plot_phylorug(tree, rug_mt)
    expect_null(result)
  })
})

# ---- runs cleanly across the branches ---------------------------------------

test_that("runs in support mode (default)", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(on_null_device(plot_phylorug(tree, rug_mt)))
})

test_that("runs in presence mode", {
  tree   <- make_test_tree()
  rug_mt <- make_presence_mt(tree)
  expect_no_error(on_null_device(plot_phylorug(tree, rug_mt)))
})

test_that("runs with show_values = TRUE", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(on_null_device(plot_phylorug(tree, rug_mt, show_values = TRUE)))
})

test_that("runs with legend = FALSE", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(on_null_device(plot_phylorug(tree, rug_mt, legend = FALSE)))
})

test_that("runs with gradient_legend = FALSE", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(
    on_null_device(plot_phylorug(tree, rug_mt, gradient_legend = FALSE))
  )
})

test_that("runs with dot_unanimous = FALSE", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(
    on_null_device(plot_phylorug(tree, rug_mt, dot_unanimous = FALSE))
  )
})

test_that("runs with custom hues", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(
    on_null_device(plot_phylorug(tree, rug_mt, hues = c("tomato", "steelblue")))
  )
})

test_that("runs with a user-fixed grid shape", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(
    on_null_device(plot_phylorug(tree, rug_mt, n_rows = 1, n_cols = 2))
  )
})

test_that("passes ... through to plot.phylo", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  expect_no_error(
    on_null_device(
      plot_phylorug(tree, rug_mt, cex = 0.5, edge.width = 0.8, no.margin = TRUE)
    )
  )
})

# ---- node-split edge cases --------------------------------------------------

test_that("runs when every node is unanimous (dots only, no rug)", {
  tree   <- make_test_tree()
  rug_mt <- make_all_present_mt(tree)
  expect_no_error(on_null_device(plot_phylorug(tree, rug_mt)))
})

test_that("runs when every node is variable (rug only, no dots)", {
  tree   <- make_test_tree()
  rug_mt <- make_all_variable_mt(tree)
  expect_no_error(on_null_device(plot_phylorug(tree, rug_mt)))
})

test_that("accepts a data frame rug_mt", {
  # The validation allows a data frame; this confirms the drawing path also
  # tolerates one. If this fails, it signals the function needs a matrix
  # internally (a real finding worth a coercion step), not a broken test.
  tree   <- make_test_tree()
  rug_mt <- as.data.frame(make_support_mt(tree))
  expect_no_error(on_null_device(plot_phylorug(tree, rug_mt)))
})
