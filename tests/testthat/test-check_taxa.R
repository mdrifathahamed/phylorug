# Tests for check_taxa()
# ---- helpers ----------------------------------------------------------------
# A tree with tips Sp_1 .. Sp_n
make_tree <- function(n_tips = 5L) {
  tr <- ape::rtree(n_tips)
  tr$tip.label <- paste0("Sp_", seq_len(n_tips))
  tr
}

# A named list of identical-taxa trees
make_tree_list <- function(n_trees = 3L, n_tips = 5L) {
  trees <- lapply(seq_len(n_trees), function(i) make_tree(n_tips))
  names(trees) <- paste0("tree", seq_len(n_trees))
  trees
}

# A tree with one taxon dropped
make_tree_missing <- function(n_tips = 5L, drop = "Sp_1") {
  tr <- make_tree(n_tips)
  ape::drop.tip(tr, drop)
}

# A tree with one extra taxon
make_tree_extra <- function(n_tips = 5L) {
  tr <- ape::rtree(n_tips + 1L)
  tr$tip.label <- c(paste0("Sp_", seq_len(n_tips)), "Sp_EXTRA")
  tr
}

# A pool of trees sharing the same taxa
make_pool <- function(n_trees = 3L, n_tips = 5L) {
  pool <- lapply(seq_len(n_trees), function(i) make_tree(n_tips))
  class(pool) <- "multiPhylo"
  pool
}

# A pool whose trees disagree about their own taxa -- a corrupt file
make_bad_pool <- function() {
  pool <- list(make_tree(5L), make_tree_missing(5L))
  class(pool) <- "multiPhylo"
  pool
}

make_compressed_multiphylo <- function() {
  trees <- make_tree_list()
  class(trees) <- "multiPhylo"
  ape::.compressTipLabel(trees)
}


# ---- input validation -------------------------------------------------------
test_that("stops when backbone is not a tree", {
  expect_error(
    check_taxa("not a tree", make_tree_list()),
    "`backbone` must be a `phylo` object"
  )
  expect_error(
    check_taxa(list(a = 1), make_tree_list()),
    "`backbone` must be a `phylo` object"
  )
})
test_that("stops when trees contains a compressed multiPhylo element", {
  bb <- make_tree()
  compressed <- make_compressed_multiphylo()
  expect_error(
    check_taxa(bb, list(pool = compressed)),
    "compressed multiPhylo"
  )
})
test_that("stops when trees is empty", {
  bb <- make_tree()
  expect_error(check_taxa(bb, list()), "is empty")
})

test_that("stops when trees contains NULL elements", {
  bb    <- make_tree()
  trees <- make_tree_list()
  # A list literal keeps the NULL; `trees$tree2 <- NULL` would delete it.
  bad <- list(tree1 = trees[[1]], tree2 = NULL, tree3 = trees[[3]])
  expect_error(check_taxa(bb, bad), "NULL elements")
})

test_that("stops when trees list contains a non-phylo element", {
  bb <- make_tree()
  bad <- list(good = make_tree(), bad = "not a tree")
  expect_error(check_taxa(bb, bad, verbose = FALSE))
})
test_that("stops when trees list contains a non-phylo element", {
  bb <- make_tree()
  bad <- list(good = make_tree(), bad = "not a tree")
  expect_error(
    check_taxa(bb, bad, verbose = FALSE),
    "Invalid elements at positions: 2"
  )
})

# ---- backbone edge cases -----------------------------------------------------
test_that("multiPhylo backbone with inconsistent taxa errors", {
  tr1 <- make_tree(5L)
  tr2 <- make_tree_missing(5L, "Sp_1")
  bad_bb <- structure(list(tr1, tr2), class = "multiPhylo")
  expect_error(
    check_taxa(bad_bb, make_tree_list(), verbose = FALSE),
    "do not share the same taxa"
  )
})

test_that("backbone with duplicate tip labels does not crash", {
  bb <- make_tree(5L)
  bb$tip.label[2] <- bb$tip.label[1]  # duplicate Sp_1
  trees <- make_tree_list()
  ok <- check_taxa(bb, trees, verbose = FALSE)
  expect_type(ok, "logical")
})

test_that("single-tip backbone works", {
  bb <- ape::read.tree(text = "(Sp_1:1);")
  tr <- ape::read.tree(text = "(Sp_1:2);")  # different edge length, not identical()
  expect_true(check_taxa(bb, list(single = tr), verbose = FALSE))
})


# ---- identical taxa ---------------------------------------------------------
test_that("returns TRUE when every analysis shares the backbone taxa", {
  bb    <- make_tree()
  trees <- make_tree_list()
  expect_true(check_taxa(bb, trees, verbose = FALSE))
})

test_that("backbone is silently removed when present in trees", {
  bb    <- make_tree()
  trees <- make_tree_list()
  trees[["backbone"]] <- bb          # add backbone into the list
  expect_true(check_taxa(bb, trees, verbose = FALSE))
  # diagnostics should not include a row for the backbone itself
  diag <- attr(check_taxa(bb, trees, verbose = FALSE), "diagnostics")
  expect_false("backbone" %in% diag$comparison)
})

test_that("stops when trees contains only the backbone", {
  bb <- make_tree()
  expect_error(
    check_taxa(bb, list(backbone = bb), verbose = FALSE),
    "contains only the backbone"
  )
})

test_that("taxon order does not matter", {
  # Comparison is on sorted sets, not on tip.label order.
  bb          <- make_tree()
  tr          <- make_tree()
  tr$tip.label <- rev(tr$tip.label)
  expect_true(check_taxa(bb, list(rev = tr), verbose = FALSE))
})

test_that("stops when trees is not a list", {
  bb <- make_tree()
  expect_error(check_taxa(bb, "not a list"), "must be a list")
})

test_that("a single comparison tree (length-1 list) works", {
  bb <- make_tree()
  tr <- make_tree()
  ok <- check_taxa(bb, list(only = tr), verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_true(ok)
  expect_equal(nrow(diag), 1L)
})


# ---- backbone removal edge cases --------------------------------------------
test_that("backbone appearing twice in trees is fully removed", {
  bb    <- make_tree()
  trees <- make_tree_list()
  trees[["bb_copy1"]] <- bb
  trees[["bb_copy2"]] <- bb
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_false("bb_copy1" %in% diag$comparison)
  expect_false("bb_copy2" %in% diag$comparison)
})

test_that("a tree with same topology but different object is NOT removed", {
  bb <- make_tree(5L)
  # Same taxa, but a distinct R object -- identical() should say FALSE
  # because tree structure (edge order, etc.) differs from rtree() calls.
  tr_same_taxa <- make_tree(5L)
  trees <- list(lookalike = tr_same_taxa)
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  # It should still appear as a comparison row, not be silently dropped
  expect_true("lookalike" %in% diag$comparison)
})


# ---- three-way classification -----------------------------------------------
test_that("an analysis missing a backbone taxon is classified as missing", {
  bb    <- make_tree(5L)
  trees <- list(good = make_tree(5L), gappy = make_tree_missing(5L, "Sp_1"))
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_false(ok)
  expect_equal(diag$status[diag$comparison == "gappy"], "missing")
  expect_equal(diag$missing[diag$comparison == "gappy"], "Sp_1")
  expect_equal(diag$extra[diag$comparison == "gappy"], "")
})

test_that("an analysis with extra taxa is classified as superset", {
  # Every backbone clade can still be evaluated, so this is not the same
  # problem as a missing taxon.
  bb    <- make_tree(5L)
  trees <- list(good = make_tree(5L), big = make_tree_extra(5L))
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_false(ok)
  expect_equal(diag$status[diag$comparison == "big"], "superset")
  expect_equal(diag$missing[diag$comparison == "big"], "")
  expect_equal(diag$extra[diag$comparison == "big"], "Sp_EXTRA")
})

test_that("an analysis both missing and gaining taxa is classified as missing", {
  # A missing backbone taxon is the serious case, and takes precedence.
  bb <- make_tree(5L)
  tr <- ape::rtree(5L)
  tr$tip.label <- c(paste0("Sp_", 2:5), "Sp_EXTRA")   # no Sp_1, plus an extra
  ok   <- check_taxa(bb, list(odd = tr), verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_equal(diag$status, "missing")
  expect_equal(diag$missing, "Sp_1")
  expect_equal(diag$extra, "Sp_EXTRA")
})

test_that("mixed statuses are all recorded", {
  bb <- make_tree(5L)
  trees <- list(
    same  = make_tree(5L),
    gappy = make_tree_missing(5L, "Sp_2"),
    big   = make_tree_extra(5L)
  )
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_false(ok)
  expect_equal(
    diag$status[match(c("same", "gappy", "big"), diag$comparison)],
    c("identical", "missing", "superset")
  )
})

test_that("all-superset trees still return FALSE (not identical)", {
  # check_taxa answers "are taxa sets identical", not "can I still score".
  bb    <- make_tree(5L)
  trees <- list(big1 = make_tree_extra(5L), big2 = make_tree_extra(5L))
  ok <- check_taxa(bb, trees, verbose = FALSE)
  expect_false(ok)
  diag <- attr(ok, "diagnostics")
  expect_true(all(diag$status == "superset"))
})

test_that("multiple missing taxa are all listed in the missing column", {
  bb <- make_tree(6L)
  tr <- ape::drop.tip(make_tree(6L), c("Sp_1", "Sp_3"))
  ok   <- check_taxa(bb, list(gappy = tr), verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_true(grepl("Sp_1", diag$missing))
  expect_true(grepl("Sp_3", diag$missing))
})

test_that("same taxon missing from multiple trees is independent per row", {
  bb <- make_tree(5L)
  trees <- list(
    gappy1 = make_tree_missing(5L, "Sp_1"),
    gappy2 = make_tree_missing(5L, "Sp_1")
  )
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_equal(diag$missing[diag$comparison == "gappy1"], "Sp_1")
  expect_equal(diag$missing[diag$comparison == "gappy2"], "Sp_1")
  expect_equal(nrow(diag), 2L)
})


# ---- diagnostics attribute --------------------------------------------------
test_that("diagnostics are attached even when everything matches", {
  bb    <- make_tree()
  trees <- make_tree_list()
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_s3_class(diag, "data.frame")
  expect_equal(nrow(diag), length(trees))
  expect_setequal(diag$comparison, names(trees))
  expect_true(all(diag$status == "identical"))
})

test_that("diagnostics are attached when verbose = FALSE", {
  # The old version threw the diagnostics away in silent mode.
  bb    <- make_tree(5L)
  trees <- list(gappy = make_tree_missing(5L, "Sp_1"))
  ok <- check_taxa(bb, trees, verbose = FALSE)
  expect_false(is.null(attr(ok, "diagnostics")))
})

test_that("diagnostics has the expected columns", {
  bb    <- make_tree()
  trees <- make_tree_list()
  diag <- attr(check_taxa(bb, trees, verbose = FALSE), "diagnostics")
  expect_named(
    diag,
    c("comparison", "status", "n_taxa", "missing", "extra")
  )
  expect_type(diag$n_taxa, "integer")
})

test_that("analyses with no names are given positional names", {
  bb    <- make_tree()
  trees <- unname(make_tree_list(2L))
  diag <- attr(check_taxa(bb, trees, verbose = FALSE), "diagnostics")
  expect_equal(diag$comparison, c("tree_1", "tree_2"))
})

test_that("n_taxa is correct for identical, superset, and missing trees", {
  bb <- make_tree(5L)
  trees <- list(
    same  = make_tree(5L),           # 5 taxa
    big   = make_tree_extra(5L),     # 6 taxa
    gappy = make_tree_missing(5L, "Sp_1")  # 4 taxa
  )
  diag <- attr(check_taxa(bb, trees, verbose = FALSE), "diagnostics")
  expect_equal(diag$n_taxa[diag$comparison == "same"], 5L)
  expect_equal(diag$n_taxa[diag$comparison == "big"], 6L)
  expect_equal(diag$n_taxa[diag$comparison == "gappy"], 4L)
})

test_that("missing and extra columns are empty strings when identical", {
  bb    <- make_tree(5L)
  trees <- list(same = make_tree(5L))
  diag  <- attr(check_taxa(bb, trees, verbose = FALSE), "diagnostics")
  expect_equal(diag$missing, "")
  expect_equal(diag$extra, "")
})


# ---- pools --------------------------------------------------------------
test_that("a pooled analysis is checked like any other", {
  bb    <- make_tree(5L)
  trees <- list(single = make_tree(5L), pool = make_pool(3L, 5L))
  expect_true(check_taxa(bb, trees, verbose = FALSE))
})

test_that("a pool missing a backbone taxon is caught", {
  bb <- make_tree(5L)
  pool <- lapply(1:3, function(i) make_tree_missing(5L, "Sp_1"))
  class(pool) <- "multiPhylo"
  ok   <- check_taxa(bb, list(gappy = pool), verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_false(ok)
  expect_equal(diag$status, "missing")
  expect_equal(diag$missing, "Sp_1")
})

test_that("a pool whose trees disagree about their own taxa is an error", {
  # Not a phylogenetic finding. A pool must be the equally optimal trees from
  # one search, run on one dataset.
  bb <- make_tree(5L)
  expect_error(
    check_taxa(bb, list(broken = make_bad_pool()), verbose = FALSE),
    "do not share the same taxa"
  )
})

test_that("a multiPhylo backbone is accepted", {
  bb    <- make_pool(2L, 5L)
  trees <- make_tree_list(2L, 5L)
  expect_true(check_taxa(bb, trees, verbose = FALSE))
})

test_that("a single-tree pool (multiPhylo of length 1) works", {
  tr <- make_tree(5L)
  pool <- structure(list(tr), class = "multiPhylo")
  bb <- make_tree(5L)
  ok <- check_taxa(bb, list(single_pool = pool), verbose = FALSE)
  expect_true(ok)
})

test_that("mixed phylo and multiPhylo elements both work in one call", {
  bb <- make_tree(5L)
  trees <- list(
    single = make_tree(5L),
    pool   = make_pool(3L, 5L)
  )
  ok   <- check_taxa(bb, trees, verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_true(ok)
  expect_equal(nrow(diag), 2L)
})


# ---- verbose output ---------------------------------------------------------
test_that("verbose reports success", {
  bb    <- make_tree()
  trees <- make_tree_list()
  expect_message(
    check_taxa(bb, trees, verbose = TRUE),
    "share the same"
  )
})

test_that("verbose names the analyses missing backbone taxa", {
  bb    <- make_tree(5L)
  trees <- list(good = make_tree(5L), gappy = make_tree_missing(5L, "Sp_1"))
  expect_message(
    check_taxa(bb, trees, verbose = TRUE),
    "MISSING backbone taxa: gappy"
  )
})

test_that("verbose says no action is needed for a superset", {
  bb    <- make_tree(5L)
  trees <- list(big = make_tree_extra(5L))
  expect_message(
    check_taxa(bb, trees, verbose = TRUE),
    "no action is\\s+needed"
  )
})

test_that("verbose = FALSE suppresses all messages", {
  bb    <- make_tree(5L)
  trees <- list(gappy = make_tree_missing(5L, "Sp_1"))
  expect_no_message(check_taxa(bb, trees, verbose = FALSE))
})

test_that("verbose fires both superset and missing messages when both occur", {
  bb <- make_tree(5L)
  trees <- list(
    big   = make_tree_extra(5L),
    gappy = make_tree_missing(5L, "Sp_1")
  )
  msgs <- testthat::capture_messages(
    check_taxa(bb, trees, verbose = TRUE)
  )
  combined <- paste(msgs, collapse = " ")
  expect_true(grepl("no action is\\s+needed", combined))
  expect_true(grepl("MISSING backbone taxa", combined))
})

test_that("success message includes the correct taxon count", {
  bb    <- make_tree(5L)
  trees <- make_tree_list(2L, 5L)
  expect_message(
    check_taxa(bb, trees, verbose = TRUE),
    "same 5 taxa"
  )
})


# ---- return value -----------------------------------------------------------
test_that("returns a single logical", {
  bb    <- make_tree()
  trees <- make_tree_list()
  result <- check_taxa(bb, trees, verbose = FALSE)
  expect_type(result, "logical")
  expect_length(result, 1L)
})


# ---- taxon matching semantics ------------------------------------------------
test_that("taxon matching is case-sensitive", {
  bb <- make_tree(3L)
  tr <- make_tree(3L)
  tr$tip.label <- tolower(tr$tip.label)  # "sp_1" instead of "Sp_1"
  ok   <- check_taxa(bb, list(lower = tr), verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_false(ok)
  expect_equal(diag$status, "missing")
})

test_that("leading/trailing whitespace in labels causes a mismatch", {
  bb <- make_tree(3L)
  tr <- make_tree(3L)
  tr$tip.label[1] <- paste0(" ", tr$tip.label[1])
  ok   <- check_taxa(bb, list(spaced = tr), verbose = FALSE)
  diag <- attr(ok, "diagnostics")
  expect_false(ok)
})
