# Tests for node_support()

# ── helpers ───────────────────────────────────────────────────────────────────

make_labeled_tree <- function(labels) {
  tree <- ape::rtree(length(labels) + 1)
  tree$node.label <- labels
  tree
}

# ── input validation ──────────────────────────────────────────────────────────

test_that("stops when tree is not a phylo object", {
  expect_error(
    node_support(list(node.label = c("95", "88"))),
    "must be a phylogenetic tree of class \"phylo\"."
  )
})

test_that("stops when round is a vector", {
  tree <- make_labeled_tree(c("95", "88"))
  expect_error(
    node_support(tree, round = c(1, 2)),
    "must be a single integer or NULL."
  )
})

test_that("stops when round is non-numeric", {
  tree <- make_labeled_tree(c("95", "88"))
  expect_error(
    node_support(tree, round = "two"),
    "must be a single integer or NULL."
  )
})

# ── single value labels ───────────────────────────────────────────────────────

test_that("single bootstrap values fill support_1 and leave support_2 NA", {
  tree   <- make_labeled_tree(c("95", "88"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(95, 88))
  expect_true(all(is.na(result$support_2)))
})

# ── compound labels ───────────────────────────────────────────────────────────

test_that("compound labels fill both columns", {
  tree   <- make_labeled_tree(c("0.95/0.88", "1.00/0.92"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(0.95, 1.00))
  expect_equal(result$support_2, c(0.88, 0.92))
})

# ── non-numeric and empty labels ──────────────────────────────────────────────

test_that("non-numeric labels become NA", {
  tree   <- make_labeled_tree(c("root", "88"))
  result <- node_support(tree)
  expect_true(is.na(result$support_1[1]))
  expect_equal(result$support_1[2], 88)
})

test_that("empty labels become NA", {
  tree   <- make_labeled_tree(c("", "88"))
  result <- node_support(tree)
  expect_true(is.na(result$support_1[1]))
  expect_equal(result$support_1[2], 88)
})

# ── NULL labels ───────────────────────────────────────────────────────────────

test_that("returns all-NA data frame when tree has no node labels", {
  tree            <- ape::rtree(5)
  tree$node.label <- NULL
  result          <- node_support(tree)
  expect_equal(nrow(result), ape::Nnode(tree))
  expect_true(all(is.na(result$support_1)))
  expect_true(all(is.na(result$support_2)))
})

# ── short label vector (padding) ──────────────────────────────────────────────

test_that("pads when label vector is shorter than node count", {
  tree            <- ape::rtree(5)
  n_node          <- ape::Nnode(tree)
  tree$node.label <- rep("90", n_node - 1)
  result          <- node_support(tree)
  expect_equal(nrow(result), n_node)
  expect_true(is.na(result$support_1[n_node]))
})

# ── rounding ──────────────────────────────────────────────────────────────────

test_that("rounds support values when round is supplied", {
  tree   <- make_labeled_tree(c("0.9567/0.8812", "1.0000/0.9234"))
  result <- node_support(tree, round = 2)
  expect_equal(result$support_1, c(0.96, 1.00))
  expect_equal(result$support_2, c(0.88, 0.92))
})

# ── return shape ──────────────────────────────────────────────────────────────

test_that("returns a data frame with the two expected columns", {
  tree   <- make_labeled_tree(c("95", "88"))
  result <- node_support(tree)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("support_1", "support_2"))
})
