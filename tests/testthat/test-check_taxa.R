# Tests for check_taxa()

# ── input validation ──────────────────────────────────────────────────────────

test_that("stops when trees is not a list or multiPhylo", {
  mock_tree <- list(tip.label = c("A", "B", "C"))
  class(mock_tree) <- "phylo"
  expect_error(check_taxa(mock_tree), "must be a list")
})

test_that("stops when trees is an empty list", {
  expect_error(
    check_taxa(list()),
    "is empty"
  )
})

test_that("stops when trees contains NULL elements", {
  trees <- make_tree_list()
  bad_list <- list(
    tree1 = trees[[1]],
    tree2 = NULL,
    tree3 = trees[[3]]
  )
  expect_error(
    check_taxa(bad_list),
    "NULL elements"
  )
})

test_that("stops when trees is a compressed multiPhylo", {
  compressed <- make_compressed_multiphylo()
  expect_error(
    check_taxa(compressed),
    "compressed"
  )
})

# ── correct behaviour ─────────────────────────────────────────────────────────

test_that("returns TRUE when all trees have identical taxa", {
  expect_true(check_taxa(make_tree_list(), verbose = FALSE))
})

test_that("returns FALSE when one tree is missing a taxon", {
  trees <- make_tree_list()
  bad_tree <- trees[[1]]
  bad_tree$tip.label <- bad_tree$tip.label[-1]
  bad_list <- list(tree1 = trees[[1]], tree2 = bad_tree)
  expect_false(check_taxa(bad_list, verbose = FALSE))
})

test_that("returns FALSE when one tree has an extra taxon", {
  trees <- make_tree_list(n_tips = 5)
  bad_tree <- trees[[1]]
  bad_tree$tip.label[1] <- "Sp_EXTRA"
  bad_list <- list(tree1 = trees[[1]], tree2 = bad_tree)
  expect_false(check_taxa(bad_list, verbose = FALSE))
})

test_that("accepts an uncompressed multiPhylo object", {
  raw <- make_tree_list()
  class(raw) <- "multiPhylo"
  expect_true(check_taxa(raw, verbose = FALSE))
})

# ── verbose output ────────────────────────────────────────────────────────────

test_that("verbose = TRUE produces a message on success", {
  trees <- make_tree_list()
  expect_message(
    check_taxa(trees, verbose = TRUE)
  )
})

test_that("verbose = TRUE produces a message on mismatch", {
  trees <- make_tree_list()
  bad_tree <- trees[[1]]
  bad_tree$tip.label <- bad_tree$tip.label[-1]
  bad_list <- list(tree1 = trees[[1]], tree2 = bad_tree)
  expect_message(
    check_taxa(bad_list, verbose = TRUE),
    "mismatch"
  )
})

test_that("verbose = FALSE suppresses all messages", {
  trees <- make_tree_list()
  expect_no_message(
    check_taxa(trees, verbose = FALSE)
  )
})

# ── return type ───────────────────────────────────────────────────────────────

test_that("always returns a single logical value", {
  trees <- make_tree_list()
  result <- check_taxa(trees, verbose = FALSE)
  expect_type(result, "logical")
  expect_length(result, 1)
})
