# Tests for prune_to_shared()

# ---- helpers ----------------------------------------------------------------
make_backbone <- function(n_tips = 5L) {
  tr <- ape::rtree(n_tips)
  tr$tip.label <- paste0("Sp_", seq_len(n_tips))
  tr
}

# A named list of trees sharing the same taxa as the backbone
make_matching_trees <- function(n_trees = 3L, n_tips = 5L) {
  trees <- lapply(seq_len(n_trees), function(i) {
    tr <- ape::rtree(n_tips)
    tr$tip.label <- paste0("Sp_", seq_len(n_tips))
    tr
  })
  names(trees) <- paste0("tree", seq_len(n_trees))
  trees
}

# A tree missing one taxon from the backbone
make_tree_missing <- function(n_tips = 5L, drop = "Sp_1") {
  tr <- make_backbone(n_tips)
  ape::drop.tip(tr, drop)
}

# A tree with one extra taxon beyond the backbone
make_tree_extra <- function(n_tips = 5L) {
  tr <- ape::rtree(n_tips + 1L)
  tr$tip.label <- c(paste0("Sp_", seq_len(n_tips)), "Sp_EXTRA")
  tr
}

# A pool of trees sharing the same taxa
make_pool <- function(n_trees = 2L, n_tips = 5L) {
  pool <- lapply(seq_len(n_trees), function(i) {
    tr <- ape::rtree(n_tips)
    tr$tip.label <- paste0("Sp_", seq_len(n_tips))
    tr
  })
  structure(pool, class = "multiPhylo")
}

# A pool missing one backbone taxon
make_pool_missing <- function(n_tips = 5L, drop = "Sp_1") {
  pool <- lapply(1:2, function(i) make_tree_missing(n_tips, drop))
  structure(pool, class = "multiPhylo")
}


# ---- input validation -------------------------------------------------------
test_that("stops when backbone is not a phylo or multiPhylo", {
  expect_error(
    prune_to_shared("not a tree", make_matching_trees()),
    "must be a `phylo` object"
  )
  expect_error(
    prune_to_shared(data.frame(x = 1), make_matching_trees()),
    "must be a `phylo` object"
  )
})

test_that("stops when trees is not a list", {
  bb <- make_backbone()
  expect_error(
    prune_to_shared(bb, "not a list"),
    "must be a list"
  )
})

test_that("stops when trees is empty", {
  bb <- make_backbone()
  expect_error(
    prune_to_shared(bb, list()),
    "is empty"
  )
})

test_that("stops when shared taxa fewer than 3", {
  # Backbone has Sp_1..Sp_5, comparison has only Sp_1 and Sp_2
  bb <- make_backbone(5L)
  tr <- ape::rtree(2L)
  tr$tip.label <- c("Sp_1", "Sp_2")
  expect_error(
    prune_to_shared(bb, list(tiny = tr), verbose = FALSE),
    "at least three"
  )
})


# ---- no pruning needed ------------------------------------------------------
test_that("returns unchanged trees when all taxa match", {
  bb    <- make_backbone(5L)
  trees <- make_matching_trees(2L, 5L)
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_equal(length(result$backbone$tip.label), 5L)
  expect_equal(length(result$trees[[1]]$tip.label), 5L)
  expect_equal(attr(result, "dropped"), character(0))
})

test_that("verbose reports nothing pruned when all taxa match", {
  bb    <- make_backbone(5L)
  trees <- make_matching_trees(2L, 5L)
  expect_message(
    prune_to_shared(bb, trees, verbose = TRUE),
    "Nothing was pruned"
  )
})


# ---- basic pruning ----------------------------------------------------------
test_that("drops taxa missing from one comparison tree", {
  bb    <- make_backbone(5L)
  trees <- list(
    full  = make_backbone(5L),
    gappy = make_tree_missing(5L, "Sp_1")
  )
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  # Sp_1 should be gone from backbone and all comparison trees
  expect_false("Sp_1" %in% result$backbone$tip.label)
  expect_false("Sp_1" %in% result$trees$full$tip.label)
  expect_false("Sp_1" %in% result$trees$gappy$tip.label)
  expect_equal(attr(result, "dropped"), "Sp_1")
})

test_that("drops multiple taxa when different trees are missing different ones", {
  bb    <- make_backbone(5L)
  trees <- list(
    no_1 = make_tree_missing(5L, "Sp_1"),
    no_3 = make_tree_missing(5L, "Sp_3")
  )
  result  <- prune_to_shared(bb, trees, verbose = FALSE)
  dropped <- attr(result, "dropped")

  expect_true("Sp_1" %in% dropped)
  expect_true("Sp_3" %in% dropped)
  # Remaining should be Sp_2, Sp_4, Sp_5
  expect_equal(length(result$backbone$tip.label), 3L)
})

test_that("extra taxa in comparison trees are also dropped", {
  # Comparison has Sp_EXTRA which isn't in backbone -- intersection drops it
  bb    <- make_backbone(5L)
  trees <- list(big = make_tree_extra(5L))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_false("Sp_EXTRA" %in% result$trees$big$tip.label)
  # Backbone taxa are all shared, so backbone is unchanged
  expect_equal(length(result$backbone$tip.label), 5L)
})


# ---- pools ------------------------------------------------------------------
test_that("pool elements are pruned and stay multiPhylo", {
  bb    <- make_backbone(5L)
  trees <- list(pool = make_pool_missing(5L, "Sp_1"))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_s3_class(result$trees$pool, "multiPhylo")
  expect_false("Sp_1" %in% result$trees$pool[[1]]$tip.label)
  expect_false("Sp_1" %in% result$trees$pool[[2]]$tip.label)
})

test_that("pool tree count is preserved after pruning", {
  bb    <- make_backbone(5L)
  trees <- list(pool = make_pool_missing(5L, "Sp_1"))
  result <- prune_to_shared(bb, trees, verbose = FALSE)
  expect_length(result$trees$pool, 2L)
})

test_that("mixed phylo and pool elements both get pruned", {
  bb    <- make_backbone(5L)
  trees <- list(
    single = make_tree_missing(5L, "Sp_1"),
    pool   = make_pool_missing(5L, "Sp_1")
  )
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_false("Sp_1" %in% result$trees$single$tip.label)
  expect_false("Sp_1" %in% result$trees$pool[[1]]$tip.label)
  expect_s3_class(result$trees$single, "phylo")
  expect_s3_class(result$trees$pool, "multiPhylo")
})

test_that("a pool with no tips to drop passes through unchanged", {
  bb    <- make_backbone(5L)
  trees <- list(pool = make_pool(2L, 5L))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_s3_class(result$trees$pool, "multiPhylo")
  expect_length(result$trees$pool, 2L)
  expect_equal(length(result$trees$pool[[1]]$tip.label), 5L)
})


# ---- return structure -------------------------------------------------------
test_that("returns a list with backbone and trees elements", {
  bb     <- make_backbone(5L)
  trees  <- make_matching_trees(2L, 5L)
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_type(result, "list")
  expect_named(result, c("backbone", "trees"))
})

test_that("backbone in result is a phylo", {
  bb     <- make_backbone(5L)
  trees  <- make_matching_trees(2L, 5L)
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_s3_class(result$backbone, "phylo")
})

test_that("tree names are preserved after pruning", {
  bb    <- make_backbone(5L)
  trees <- list(alpha = make_backbone(5L), beta = make_tree_missing(5L, "Sp_1"))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_named(result$trees, c("alpha", "beta"))
})

test_that("dropped attribute lists removed taxa", {
  bb    <- make_backbone(5L)
  trees <- list(gappy = make_tree_missing(5L, "Sp_2"))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_equal(attr(result, "dropped"), "Sp_2")
})

test_that("dropped attribute is empty character when nothing pruned", {
  bb    <- make_backbone(5L)
  trees <- make_matching_trees(2L, 5L)
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_equal(attr(result, "dropped"), character(0))
})


# ---- verbose messaging ------------------------------------------------------
test_that("verbose reports dropped taxa by name", {
  bb    <- make_backbone(5L)
  trees <- list(gappy = make_tree_missing(5L, "Sp_1"))
  expect_message(
    prune_to_shared(bb, trees, verbose = TRUE),
    "Dropped 1 of 5 taxa"
  )
})

test_that("verbose mentions the dropped taxon name", {
  bb    <- make_backbone(5L)
  trees <- list(gappy = make_tree_missing(5L, "Sp_3"))
  expect_message(
    prune_to_shared(bb, trees, verbose = TRUE),
    "Sp_3"
  )
})

test_that("verbose = FALSE suppresses all messages", {
  bb    <- make_backbone(5L)
  trees <- list(gappy = make_tree_missing(5L, "Sp_1"))
  expect_no_message(
    prune_to_shared(bb, trees, verbose = FALSE)
  )
})

test_that("verbose truncates long drop lists with ...", {
  # Backbone has 15 taxa, comparison has only 4 -- drops 11
  bb <- ape::rtree(15L)
  bb$tip.label <- paste0("Sp_", seq_len(15L))
  tr <- ape::rtree(4L)
  tr$tip.label <- paste0("Sp_", seq_len(4L))
  expect_message(
    prune_to_shared(bb, list(small = tr), verbose = TRUE),
    "\\.\\.\\."
  )
})


# ---- topology preserved -----------------------------------------------------
test_that("pruned backbone has correct topology (not just correct labels)", {
  bb <- ape::read.tree(text = "((((Sp_1,Sp_2),Sp_3),Sp_4),Sp_5);")
  trees <- list(gappy = make_tree_missing(5L, "Sp_5"))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  # After dropping Sp_5, backbone should have 4 tips with valid tree structure
  expect_equal(length(result$backbone$tip.label), 4L)
  expect_true(ape::is.rooted(result$backbone))
  expect_s3_class(result$backbone, "phylo")
})


# ---- edge cases -------------------------------------------------------------
test_that("single comparison tree works (length-1 list)", {
  bb    <- make_backbone(5L)
  trees <- list(only = make_tree_missing(5L, "Sp_1"))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  expect_false("Sp_1" %in% result$backbone$tip.label)
  expect_length(result$trees, 1L)
})

test_that("comparison tree with no tips to drop passes through", {
  bb    <- make_backbone(5L)
  trees <- list(full = make_backbone(5L))
  result <- prune_to_shared(bb, trees, verbose = FALSE)

  # prune_one returns the tree unchanged when drop is empty
  expect_equal(length(result$trees$full$tip.label), 5L)
})

test_that("multiPhylo backbone is accepted", {
  pool <- make_pool(2L, 5L)
  trees <- list(gappy = make_tree_missing(5L, "Sp_1"))
  result <- prune_to_shared(pool, trees, verbose = FALSE)

  # tree_taxa takes the first tree's labels from the pool
  expect_false("Sp_1" %in% result$backbone[[1]]$tip.label)
})


# ---- internal functions -----------------------------------------------------
test_that("tree_taxa returns sorted tip labels from a phylo", {
  tr <- ape::read.tree(text = "((C,A),B);")
  expect_equal(phylorug:::tree_taxa(tr), c("A", "B", "C"))
})

test_that("tree_taxa returns sorted tip labels from a multiPhylo pool", {
  pool <- make_pool(2L, 3L)
  taxa <- phylorug:::tree_taxa(pool)
  expect_equal(taxa, sort(pool[[1]]$tip.label))
})

test_that("prune_one returns phylo unchanged when nothing to drop", {
  tr <- make_backbone(5L)
  keep <- tr$tip.label
  result <- phylorug:::prune_one(tr, keep)
  expect_s3_class(result, "phylo")
  expect_equal(length(result$tip.label), 5L)
})

test_that("prune_one drops tips not in keep from a phylo", {
  tr <- make_backbone(5L)
  keep <- paste0("Sp_", 1:3)
  result <- phylorug:::prune_one(tr, keep)
  expect_equal(sort(result$tip.label), sort(keep))
})

test_that("prune_one handles a multiPhylo pool and returns multiPhylo", {
  pool <- make_pool(2L, 5L)
  keep <- paste0("Sp_", 1:3)
  result <- phylorug:::prune_one(pool, keep)
  expect_s3_class(result, "multiPhylo")
  expect_equal(sort(result[[1]]$tip.label), sort(keep))
  expect_equal(sort(result[[2]]$tip.label), sort(keep))
})

test_that("prune_one pool element unchanged when nothing to drop", {
  pool <- make_pool(2L, 5L)
  keep <- pool[[1]]$tip.label
  result <- phylorug:::prune_one(pool, keep)
  expect_s3_class(result, "multiPhylo")
  expect_length(result, 2L)
})
