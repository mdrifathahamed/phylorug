# Tests for translate_tips()

# ── helpers ───────────────────────────────────────────────────────────────────

make_translation_dict <- function() {
  data.frame(
    from = c("Sp_1", "Sp_2", "Sp_3"),
    to   = c("Species_one", "Species_two", "Species_three"),
    stringsAsFactors = FALSE
  )
}

# ── input validation ──────────────────────────────────────────────────────────

test_that("stops when trees is not a list or multiPhylo", {
  mock_tree <- list(tip.label = c("Sp_1", "Sp_2"))
  class(mock_tree) <- "phylo"
  expect_error(
    translate_tips(mock_tree, make_translation_dict()),
    "must be a list"
  )
})

test_that("stops when data is not a data frame", {
  trees <- make_tree_list()
  expect_error(
    translate_tips(trees, data = c("Sp_1", "Sp_2")),
    "must be a data frame"
  )
})

test_that("stops when from_col is not found in data", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_error(
    translate_tips(trees, dict, from_col = "wrong_col"),
    "not found"
  )
})

test_that("stops when to_col is not found in data", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_error(
    translate_tips(trees, dict, to_col = "wrong_col"),
    "not found"
  )
})

test_that("stops when from_col contains duplicate values", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = c("Sp_1", "Sp_1", "Sp_2"),
    to   = c("Species_one", "Species_one_dup", "Species_two"),
    stringsAsFactors = FALSE
  )
  expect_error(
    translate_tips(trees, dict),
    "unique"
  )
})

# ── correct behaviour ─────────────────────────────────────────────────────────

test_that("matched tips are translated correctly", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict, verbose = FALSE)
  translated <- result[[1]]$tip.label
  expect_true(all(c("Species_one", "Species_two", "Species_three") %in% translated))
})

test_that("unmatched tips are left unchanged", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = c("Sp_1", "Sp_2"),
    to   = c("Species_one", "Species_two"),
    stringsAsFactors = FALSE
  )
  result      <- translate_tips(trees, dict, verbose = FALSE)
  tip_labels  <- result[[1]]$tip.label
  expect_true("Sp_3" %in% tip_labels)
})

test_that("list names are preserved after translation", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict, verbose = FALSE)
  expect_equal(names(result), names(trees))
})

test_that("tree topology is unchanged after translation", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict, verbose = FALSE)
  expect_equal(result[[1]]$edge, trees[[1]]$edge)
})

test_that("branch lengths are unchanged after translation", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict, verbose = FALSE)
  expect_equal(result[[1]]$edge.length, trees[[1]]$edge.length)
})

# ── verbose output ────────────────────────────────────────────────────────────

test_that("verbose = TRUE produces a message", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_message(translate_tips(trees, dict, verbose = TRUE))
})

test_that("verbose = FALSE suppresses all messages", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_no_message(translate_tips(trees, dict, verbose = FALSE))
})

# ── return type ───────────────────────────────────────────────────────────────

test_that("returns a list of phylo objects", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict, verbose = FALSE)
  expect_type(result, "list")
  expect_true(all(vapply(result, inherits, logical(1), "phylo")))
})
