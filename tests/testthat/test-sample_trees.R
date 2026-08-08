# Tests for sample_trees data
test_that("sample_trees exists and is a named list of 5 phylo objects", {
  expect_type(sample_trees, "list")
  expect_length(sample_trees, 5L)
  expect_true(all(vapply(sample_trees, inherits, logical(1), "phylo")))
  expect_false(is.null(names(sample_trees)))
})

test_that("sample_trees has the expected analysis names", {
  expect_setequal(
    names(sample_trees),
    c("70p_uce", "70p_partition_entropy", "70p_ghost",
      "70p_ASTRAL_uce", "70p_ASTRAL_partition_entropy")
  )
})

test_that("all sample_trees have 20 tips", {
  tip_counts <- vapply(sample_trees, ape::Ntip, integer(1))
  expect_true(all(tip_counts == 20L))
})

test_that("all sample_trees share the same taxa", {
  taxa <- lapply(sample_trees, function(tr) sort(tr$tip.label))
  same <- vapply(taxa, identical, logical(1), taxa[[1]])
  expect_true(all(same))
})

test_that("IQ-TREE trees have node labels (support values)", {
  iqtree_names <- c("70p_uce", "70p_partition_entropy", "70p_ghost")
  for (nm in iqtree_names) {
    expect_false(is.null(sample_trees[[nm]]$node.label))
  }
})
