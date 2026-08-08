# Tests for sample_trees data
test_that("sample_trees exists and is a named list of phylo objects", {
  expect_type(sample_trees, "list")
  expect_true(length(sample_trees) > 0L)
  expect_true(all(vapply(sample_trees, inherits, logical(1), "phylo")))
  expect_false(is.null(names(sample_trees)))
})

test_that("all sample_trees share the same taxa", {
  taxa <- lapply(sample_trees, function(tr) sort(tr$tip.label))
  same <- vapply(taxa, identical, logical(1), taxa[[1]])
  expect_true(all(same))
})

test_that("all sample_trees have at least 3 tips", {
  tip_counts <- vapply(sample_trees, ape::Ntip, integer(1))
  expect_true(all(tip_counts >= 3L))
})

test_that("all sample_trees are rooted", {
  rooted <- vapply(sample_trees, ape::is.rooted, logical(1))
  expect_true(all(rooted))
})
