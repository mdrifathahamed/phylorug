# Tests for translate_tips()

# ── helpers ───────────────────────────────────────────────────────────────────
make_tree_list <- function() {
  tr <- ape::read.tree(text = "((Sp_1:0.1,Sp_2:0.2):0.3,Sp_3:0.4);")
  list(analysis_1 = tr)
}

make_pool_tree_list <- function() {
  tr1 <- ape::read.tree(text = "((Sp_1,Sp_2),Sp_3);")
  tr2 <- ape::read.tree(text = "((Sp_1,Sp_3),Sp_2);")
  pool <- structure(list(tr1, tr2), class = "multiPhylo")
  single <- ape::read.tree(text = "((Sp_1,Sp_2),Sp_3);")
  list(ml = single, parsimony = pool)
}

make_translation_dict <- function() {
  data.frame(
    from = c("Sp_1", "Sp_2", "Sp_3"),
    to   = c("Species_one", "Species_two", "Species_three"),
    stringsAsFactors = FALSE
  )
}
# ── input validation ──────────────────────────────────────────────────────────
test_that("stops when trees is not a list", {
  expect_error(
    translate_tips("not a list", make_translation_dict(),
                   from_col = "from", to_col = "to"),
    "must be a list"
  )
  expect_error(
    translate_tips(data.frame(x = 1), make_translation_dict(),
                   from_col = "from", to_col = "to"),
    "must be a list"
  )
})

test_that("stops when a single phylo is passed instead of a list", {
  tr <- ape::read.tree(text = "((Sp_1,Sp_2),Sp_3);")
  expect_error(
    translate_tips(tr, make_translation_dict(),
                   from_col = "from", to_col = "to"),
    "must be a list"
  )
})

test_that("stops on an empty list", {
  expect_error(
    translate_tips(list(), make_translation_dict(),
                   from_col = "from", to_col = "to"),
    "empty"
  )
})

test_that("stops when list contains a non-phylo element", {
  tr <- ape::read.tree(text = "((Sp_1,Sp_2),Sp_3);")
  bad_list <- list(good = tr, bad = "not a tree")
  expect_error(
    translate_tips(bad_list, make_translation_dict(),
                   from_col = "from", to_col = "to"),
    "Invalid elements at positions: 2"
  )
})

test_that("stops when data is not a data frame", {
  trees <- make_tree_list()
  expect_error(
    translate_tips(trees, data = c("Sp_1", "Sp_2"),
                   from_col = "from", to_col = "to"),
    "must be a data frame"
  )
})

test_that("stops when both from_col and to_col are missing", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_error(
    translate_tips(trees, dict),
    "required"
  )
})

test_that("stops when from_col is not found in data", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_error(
    translate_tips(trees, dict, from_col = "wrong_col", to_col = "to"),
    "not found"
  )
})

test_that("stops when to_col is not found in data", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_error(
    translate_tips(trees, dict, from_col = "from", to_col = "wrong_col"),
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
    translate_tips(trees, dict, from_col = "from", to_col = "to"),
    "unique"
  )
})

# ── correct behaviour ─────────────────────────────────────────────────────────
test_that("matched tips are translated correctly", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  translated <- result[[1]]$tip.label
  expect_true(all(c("Species_one", "Species_two", "Species_three") %in% translated)) # nolint: line_length_linter.
})

test_that("unmatched tips are left unchanged", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = c("Sp_1", "Sp_2"),
    to   = c("Species_one", "Species_two"),
    stringsAsFactors = FALSE
  )
  result     <- translate_tips(trees, dict,
                               from_col = "from", to_col = "to",
                               verbose = FALSE)
  tip_labels <- result[[1]]$tip.label
  expect_true("Sp_3" %in% tip_labels)
})

test_that("list names are preserved after translation", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_equal(names(result), names(trees))
})

test_that("tree topology is unchanged after translation", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_equal(result[[1]]$edge, trees[[1]]$edge)
})

test_that("branch lengths are unchanged after translation", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_equal(result[[1]]$edge.length, trees[[1]]$edge.length)
})

test_that("node labels are unchanged after translation", {
  tr <- ape::read.tree(text = "((Sp_1,Sp_2)95,Sp_3)100;")
  trees <- list(analysis_1 = tr)
  dict  <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_equal(result[[1]]$node.label, trees[[1]]$node.label)
})
# ── pool-specific behaviour ───────────────────────────────────────────────────
test_that("multiPhylo pool: all inner trees are translated", {
  trees  <- make_pool_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  pool <- result$parsimony
  expect_s3_class(pool, "multiPhylo")
  expect_true("Species_one" %in% pool[[1]]$tip.label)
  expect_true("Species_one" %in% pool[[2]]$tip.label)
  expect_false("Sp_1" %in% pool[[1]]$tip.label)
  expect_false("Sp_1" %in% pool[[2]]$tip.label)
})

test_that("multiPhylo pool preserves tree count", {
  trees  <- make_pool_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_length(result$parsimony, 2L)
})

test_that("mixed phylo + multiPhylo: both translated in one call", {
  trees  <- make_pool_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_true("Species_one" %in% result$ml$tip.label)
  expect_s3_class(result$ml, "phylo")
  expect_true("Species_one" %in% result$parsimony[[1]]$tip.label)
  expect_s3_class(result$parsimony, "multiPhylo")
})

test_that("pool_sizes attribute is preserved after translation", {
  trees <- make_pool_tree_list()
  attr(trees, "pool_sizes") <- c(ml = 1L, parsimony = 2L)
  dict  <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_equal(attr(result, "pool_sizes"), c(ml = 1L, parsimony = 2L))
})
# ── verbose output ────────────────────────────────────────────────────────────
test_that("verbose = TRUE produces a message", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_message(
    translate_tips(trees, dict,
                   from_col = "from", to_col = "to",
                   verbose = TRUE)
  )
})

test_that("verbose = FALSE suppresses all messages", {
  trees <- make_tree_list()
  dict  <- make_translation_dict()
  expect_no_message(
    translate_tips(trees, dict,
                   from_col = "from", to_col = "to",
                   verbose = FALSE)
  )
})
# ── return type ───────────────────────────────────────────────────────────────
test_that("returns a list of phylo objects for single-tree analyses", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_type(result, "list")
  expect_true(all(vapply(result, inherits, logical(1), "phylo")))
})

test_that("pool elements stay multiPhylo in output", {
  trees  <- make_pool_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_s3_class(result$parsimony, "multiPhylo")
  expect_s3_class(result$ml, "phylo")
})
# ── data edge cases ──────────────────────────────────────────────────────────
test_that("stops when from_col contains NA values", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = c("Sp_1", NA, "Sp_3"),
    to   = c("Species_one", "Species_two", "Species_three"),
    stringsAsFactors = FALSE
  )
  expect_error(
    translate_tips(trees, dict, from_col = "from", to_col = "to"),
    "NA values"
  )
})

test_that("stops when to_col contains NA values", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = c("Sp_1", "Sp_2", "Sp_3"),
    to   = c("Species_one", NA, "Species_three"),
    stringsAsFactors = FALSE
  )
  expect_error(
    translate_tips(trees, dict, from_col = "from", to_col = "to"),
    "NA values"
  )
})

test_that("zero-row data frame produces zero translations", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = character(0),
    to   = character(0),
    stringsAsFactors = FALSE
  )
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  # All tips unchanged
  expect_equal(result[[1]]$tip.label, c("Sp_1", "Sp_2", "Sp_3"))
})

test_that("factor columns from read.csv work correctly", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = factor(c("Sp_1", "Sp_2", "Sp_3")),
    to   = factor(c("Species_one", "Species_two", "Species_three"))
  )
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_true("Species_one" %in% result[[1]]$tip.label)
  expect_true("Species_two" %in% result[[1]]$tip.label)
})

# ── boundary conditions ──────────────────────────────────────────────────────
test_that("all tips match (100% translation)", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  tips <- result[[1]]$tip.label
  expect_false(any(tips %in% c("Sp_1", "Sp_2", "Sp_3")))
  expect_true(all(tips %in% c("Species_one", "Species_two", "Species_three")))
})

test_that("zero tips match (0% translation)", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = c("X_1", "X_2", "X_3"),
    to   = c("Wrong_one", "Wrong_two", "Wrong_three"),
    stringsAsFactors = FALSE
  )
  result <- translate_tips(trees, dict,
                           from_col = "from", to_col = "to",
                           verbose = FALSE)
  expect_equal(result[[1]]$tip.label, c("Sp_1", "Sp_2", "Sp_3"))
})

# ── message accuracy ─────────────────────────────────────────────────────────
test_that("verbose message shows correct counts for phylo", {
  trees <- make_tree_list()
  dict  <- data.frame(
    from = c("Sp_1", "Sp_2"),
    to   = c("Species_one", "Species_two"),
    stringsAsFactors = FALSE
  )
  expect_message(
    translate_tips(trees, dict,
                   from_col = "from", to_col = "to",
                   verbose = TRUE),
    "2 tips translated, 1 unchanged"
  )
})

test_that("verbose message shows correct counts for pool", {
  trees  <- make_pool_tree_list()
  dict   <- make_translation_dict()
  expect_message(
    translate_tips(trees, dict,
                   from_col = "from", to_col = "to",
                   verbose = TRUE),
    "6 tips translated, 0 unchanged \\(totals across pool\\)"
  )
})

test_that("verbose message includes analysis name", {
  trees  <- make_tree_list()
  dict   <- make_translation_dict()
  expect_message(
    translate_tips(trees, dict,
                   from_col = "from", to_col = "to",
                   verbose = TRUE),
    "analysis_1:"
  )
})
