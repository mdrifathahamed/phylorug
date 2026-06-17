# Tests for node_presence_matrix()

# ---- helpers ----------------------------------------------------------------

make_shared_tree <- function() {
  ape::read.tree(text = "(((A,B),C),(D,E));")
}

make_tree_list <- function() {
  list(
    tree1 = ape::read.tree(text = "(((A,B),C),(D,E));"),
    tree2 = ape::read.tree(text = "((A,(B,C)),(D,E));")
  )
}

# ---- input validation -------------------------------------------------------

test_that("stops when backbone is not a phylo object", {
  expect_error(
    node_presence_matrix("not_a_tree", make_tree_list()),
    "must be a phylogenetic tree"
  )
})

test_that("stops when trees is not a list", {
  expect_error(
    node_presence_matrix(make_shared_tree(), "not_a_list"),
    "must be a list"
  )
})

test_that("stops when support_col is not 1 or 2", {
  expect_error(
    node_presence_matrix(
      make_shared_tree(),
      make_tree_list(),
      support_col = 3
    ),
    "`support_col` must be 1 or 2"
  )
})

# ---- output structure -------------------------------------------------------

test_that("returns a matrix with the expected dimensions", {
  result <- node_presence_matrix(make_shared_tree(), make_tree_list())
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
  expect_equal(nrow(result), 4)
})

test_that("first column is node_id and columns are named after the trees", {
  result <- node_presence_matrix(make_shared_tree(), make_tree_list())
  expect_equal(colnames(result), c("node_id", "tree1", "tree2"))
})

test_that("generates generic tree names when the list is unnamed", {
  backbone  <- make_shared_tree()
  tree_list <- unname(make_tree_list())
  result    <- node_presence_matrix(backbone, tree_list)
  expect_equal(colnames(result), c("node_id", "tree_1", "tree_2"))
})

# ---- correct values ---------------------------------------------------------

test_that("cells are 1 for matching clades and NA for absent clades", {
  result <- node_presence_matrix(make_shared_tree(), make_tree_list())
  # tree1 is identical to the backbone, so every clade is present
  expect_true(all(result[, "tree1"] == 1))
  # tree2 differs, so at least one backbone clade is absent (NA)
  expect_true(any(is.na(result[, "tree2"])))
  # node_id values are positive integers
  expect_true(all(result[, "node_id"] > 0))
})

test_that("support mode returns numeric values or NA", {
  backbone <- ape::read.tree(text = "(((A,B)90,C)80,(D,E)100);")
  trees    <- list(
    t1 = ape::read.tree(text = "(((A,B)90,C)80,(D,E)100);")
  )
  result <- node_presence_matrix(backbone, trees, use_support = TRUE)
  vals   <- result[, "t1"]
  expect_true(all(is.na(vals) | is.numeric(vals)))
})
