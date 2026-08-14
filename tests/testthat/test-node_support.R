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
    "must be a phylogenetic tree"
  )
  expect_error(
    node_support("not a tree"),
    "must be a phylogenetic tree"
  )
  expect_error(
    node_support(data.frame(x = 1)),
    "must be a phylogenetic tree"
  )
})

test_that("stops when a multiPhylo is passed instead of a single tree", {
  pool <- structure(list(ape::rtree(5), ape::rtree(5)), class = "multiPhylo")
  expect_error(
    node_support(pool),
    "must be a phylogenetic tree"
  )
})

test_that("stops when digits is a vector", {
  tree <- make_labeled_tree(c("95", "88"))
  expect_error(
    node_support(tree, digits = c(1, 2)),
    "must be a single integer or NULL"
  )
})

test_that("stops when digits is non-numeric", {
  tree <- make_labeled_tree(c("95", "88"))
  expect_error(
    node_support(tree, digits = "two"),
    "must be a single integer or NULL"
  )
})


# ── single value labels ───────────────────────────────────────────────────────
test_that("single integer bootstrap fills support_1 only", {
  tree   <- make_labeled_tree(c("95", "88"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(95, 88))
  expect_true(all(is.na(result$support_2)))
  expect_true(all(is.na(result$support_3)))
})

test_that("single decimal posterior fills support_1 only", {
  tree   <- make_labeled_tree(c("0.95", "0.88"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(0.95, 0.88))
  expect_true(all(is.na(result$support_2)))
})


# ── compound labels ───────────────────────────────────────────────────────────
test_that("two-part compound labels fill support_1 and support_2", {
  tree   <- make_labeled_tree(c("0.95/0.88", "1.00/0.92"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(0.95, 1.00))
  expect_equal(result$support_2, c(0.88, 0.92))
  expect_true(all(is.na(result$support_3)))
})

test_that("three-part compound labels fill all three columns", {
  tree   <- make_labeled_tree(c("95/88/0.99", "100/92/0.87"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(95, 100))
  expect_equal(result$support_2, c(88, 92))
  expect_equal(result$support_3, c(0.99, 0.87))
})

test_that("more than three parts are truncated to three columns", {
  tree   <- make_labeled_tree(c("95/88/0.99/0.50", "100/92/0.87/0.60"))
  result <- node_support(tree)
  expect_equal(ncol(result), 3L)
  expect_equal(result$support_1, c(95, 100))
  expect_equal(result$support_2, c(88, 92))
  expect_equal(result$support_3, c(0.99, 0.87))
})

test_that("mixed numeric and non-numeric in compound label", {
  tree   <- make_labeled_tree(c("95/notanumber", "88/0.92"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(95, 88))
  expect_true(is.na(result$support_2[1]))
  expect_equal(result$support_2[2], 0.92)
})

test_that("trailing separator produces NA for missing part", {
  tree   <- make_labeled_tree(c("95/", "88/92"))
  result <- node_support(tree)
  expect_equal(result$support_1, c(95, 88))
  expect_true(is.na(result$support_2[1]))
  expect_equal(result$support_2[2], 92)
})

test_that("leading separator: empty piece is dropped, value shifts to support_1", { # nolint: line_length_linter.
  tree   <- make_labeled_tree(c("/95", "88/92"))
  result <- node_support(tree)
  expect_equal(result$support_1[1], 95)
  expect_true(is.na(result$support_2[1]))
  expect_equal(result$support_1[2], 88)
  expect_equal(result$support_2[2], 92)
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

test_that("NA in label vector becomes NA in output", {
  tree   <- make_labeled_tree(c(NA, "88"))
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
  expect_true(all(is.na(result$support_3)))
})


# ── short and long label vectors ──────────────────────────────────────────────
test_that("pads when label vector is shorter than node count", {
  tree            <- ape::rtree(5)
  n_node          <- ape::Nnode(tree)
  tree$node.label <- rep("90", n_node - 1)
  result          <- node_support(tree)
  expect_equal(nrow(result), n_node)
  expect_true(is.na(result$support_1[n_node]))
})

test_that("empty character(0) label vector is padded to all NA", {
  tree            <- ape::rtree(5)
  tree$node.label <- character(0)
  result          <- node_support(tree)
  expect_equal(nrow(result), ape::Nnode(tree))
  expect_true(all(is.na(result$support_1)))
})


# ── custom separator ─────────────────────────────────────────────────────────
test_that("custom separator works", {
  tree   <- make_labeled_tree(c("95_88", "100_92"))
  result <- node_support(tree, sep = "_")
  expect_equal(result$support_1, c(95, 100))
  expect_equal(result$support_2, c(88, 92))
})

test_that("separator not present in labels treats whole label as one value", {
  tree   <- make_labeled_tree(c("95", "88"))
  result <- node_support(tree, sep = "|")
  expect_equal(result$support_1, c(95, 88))
  expect_true(all(is.na(result$support_2)))
})


# ── rounding ──────────────────────────────────────────────────────────────────
test_that("rounds support values when digits is supplied", {
  tree   <- make_labeled_tree(c("0.9567/0.8812", "1.0000/0.9234"))
  result <- node_support(tree, digits = 2)
  expect_equal(result$support_1, c(0.96, 1.00))
  expect_equal(result$support_2, c(0.88, 0.92))
})

test_that("digits = 0 rounds to integers", {
  tree   <- make_labeled_tree(c("0.956", "0.881"))
  result <- node_support(tree, digits = 0)
  expect_equal(result$support_1, c(1, 1))
})

test_that("digits = NULL does not round", {
  tree   <- make_labeled_tree(c("0.9567"))
  result <- node_support(tree, digits = NULL)
  expect_equal(result$support_1, 0.9567)
})


# ── return structure ──────────────────────────────────────────────────────────
test_that("returns a data frame with three expected columns", {
  tree   <- make_labeled_tree(c("95", "88"))
  result <- node_support(tree)
  expect_s3_class(result, "data.frame")
  expect_named(result, c("support_1", "support_2", "support_3"))
})

test_that("row count equals number of internal nodes", {
  tree   <- make_labeled_tree(c("95", "88", "77", "66"))
  result <- node_support(tree)
  expect_equal(nrow(result), ape::Nnode(tree))
})

test_that("all columns are numeric type", {
  tree   <- make_labeled_tree(c("95/88/0.99", "100/92/0.87"))
  result <- node_support(tree)
  expect_type(result$support_1, "double")
  expect_type(result$support_2, "double")
  expect_type(result$support_3, "double")
})
