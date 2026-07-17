# Tests for node_presence_matrix()
#
# The function returns a list of two matrices: $presence and $support.

# ---- helpers ----------------------------------------------------------------

make_backbone <- function() {
  ape::read.tree(text = "(((A,B),C),(D,E));")
}

make_comparisons <- function() {
  list(
    tree1 = ape::read.tree(text = "(((A,B),C),(D,E));"),   # identical
    tree2 = ape::read.tree(text = "((A,(B,C)),(D,E));")     # different
  )
}

make_missing_taxon <- function() {
  list(gappy = ape::read.tree(text = "((A,B),(D,E));"))     # C is gone
}

make_supported <- function() {
  ape::read.tree(text = "(((A,B)90,C)80,(D,E)100);")
}

make_compound_support <- function() {
  ape::read.tree(text = "(((A,B)99/85,C)80/72,(D,E)100/98);")
}

make_triple_support <- function() {
  ape::read.tree(text = "(((A,B)99/85/0.999,C)80/72/0.95,(D,E)100/98/1);")
}

make_pool <- function() {
  pool <- list(
    ape::read.tree(text = "(((A,B),C),(D,E));"),
    ape::read.tree(text = "((A,(B,C)),(D,E));"),
    ape::read.tree(text = "(((A,B),C),(D,E));")
  )
  class(pool) <- "multiPhylo"
  pool
}

make_unrooted <- function() {
  ape::unroot(ape::read.tree(text = "(A,B,(C,(D,E)));"))
}


# ---- input validation -------------------------------------------------------

test_that("stops when backbone is not a phylo object", {
  expect_error(
    node_presence_matrix("not_a_tree", make_comparisons()),
    "must be a phylogenetic tree"
  )
})

test_that("stops when trees is not a list", {
  expect_error(
    node_presence_matrix(make_backbone(), "not_a_list"),
    "must be a list"
  )
})

test_that("stops when trees is empty", {
  expect_error(
    node_presence_matrix(make_backbone(), list()),
    "is empty"
  )
})

test_that("stops when support_col is out of range", {
  msg <- "`support_col` must be an integer or integer vector with values 1, 2, or 3."
  expect_error(
    node_presence_matrix(make_backbone(), make_comparisons(), support_col = 4),
    regexp = msg,
    fixed = TRUE
  )

  expect_error(
    node_presence_matrix(make_backbone(), make_comparisons(), support_col = 0),
    regexp = msg,
    fixed = TRUE
  )
})


# ---- rooted guard -----------------------------------------------------------

test_that("stops when backbone is unrooted", {
  expect_error(
    node_presence_matrix(make_unrooted(), make_comparisons()),
    "unrooted"
  )
})

test_that("stops when a comparison tree is unrooted", {
  expect_error(
    node_presence_matrix(make_backbone(), list(bad = make_unrooted())),
    "unrooted"
  )
})

test_that("names all unrooted comparison trees in the error", {
  trees <- list(
    good  = make_backbone(),
    bad_1 = make_unrooted(),
    bad_2 = make_unrooted()
  )
  expect_error(
    node_presence_matrix(make_backbone(), trees),
    "bad_1.*bad_2|bad_2.*bad_1"
  )
})


# ---- taxon gate -------------------------------------------------------------

test_that("stops when a comparison tree is missing a backbone taxon", {
  expect_error(
    node_presence_matrix(make_backbone(), make_missing_taxon()),
    "do not contain every backbone taxon"
  )
})

test_that("accepts a comparison tree with extra taxa", {
  trees <- list(big = ape::read.tree(text = "((((A,B),C),(D,E)),F);"))
  expect_no_error(node_presence_matrix(make_backbone(), trees))
})


# ---- output structure -------------------------------------------------------

test_that("returns a list with presence and support_1 by default", {
  result <- node_presence_matrix(make_backbone(), make_comparisons())
  expect_type(result, "list")
  expect_named(result, c("presence", "support_1"))
  expect_true(is.matrix(result$presence))
  expect_true(is.matrix(result$support_1))
})

test_that("both matrices have the same dimensions", {
  result <- node_presence_matrix(make_backbone(), make_comparisons())
  expect_equal(dim(result$presence), dim(result$support_1))
  expect_equal(ncol(result$presence), 2L)
  expect_equal(nrow(result$presence), make_backbone()$Nnode)
})

test_that("node IDs are in rownames, not a column", {
  result <- node_presence_matrix(make_backbone(), make_comparisons())

  expect_false("node_id" %in% colnames(result$presence))
  expect_equal(rownames(result$presence), as.character(6:9))
  expect_equal(rownames(result$support_1), as.character(6:9))
})

test_that("columns are named after the comparison trees", {
  result <- node_presence_matrix(make_backbone(), make_comparisons())
  expect_equal(colnames(result$presence), c("tree1", "tree2"))
  expect_equal(colnames(result$support_1), c("tree1", "tree2"))
})

test_that("generates positional names when the list is unnamed", {
  trees  <- unname(make_comparisons())
  result <- node_presence_matrix(make_backbone(), trees)
  expect_equal(colnames(result$presence), c("tree_1", "tree_2"))
})


# ---- presence matrix --------------------------------------------------------

test_that("presence: identical tree scores all 1", {
  result <- node_presence_matrix(make_backbone(), list(same = make_backbone()))
  expect_true(all(result$presence[, "same"] == 1))
})

test_that("presence: different tree has at least one 0", {
  trees  <- list(diff = ape::read.tree(text = "((A,(B,C)),(D,E));"))
  result <- node_presence_matrix(make_backbone(), trees)
  expect_true(any(result$presence[, "diff"] == 0))
})

test_that("presence: absence is 0, never NA", {
  trees  <- list(diff = ape::read.tree(text = "((A,(B,C)),(D,E));"))
  result <- node_presence_matrix(make_backbone(), trees)

  expect_false(any(is.na(result$presence)))
  expect_true(all(as.vector(result$presence) %in% c(0, 1)))
})

test_that("presence: values are between 0 and 1", {
  result <- node_presence_matrix(make_backbone(), make_comparisons())
  expect_true(all(result$presence >= 0 & result$presence <= 1))
})


# ---- pools ------------------------------------------------------------------

test_that("pool presence: identical trees give 1.0", {
  pool <- list(
    ape::read.tree(text = "(((A,B),C),(D,E));"),
    ape::read.tree(text = "(((A,B),C),(D,E));")
  )
  class(pool) <- "multiPhylo"

  result <- node_presence_matrix(make_backbone(), list(pool = pool))
  expect_true(all(result$presence[, "pool"] == 1))
})

test_that("pool presence: mixed recovery gives a fraction", {
  result <- node_presence_matrix(make_backbone(), list(pool = make_pool()))

  vals      <- result$presence[, "pool"]
  fractions <- vals[vals > 0 & vals < 1]
  expect_true(length(fractions) > 0)
})

test_that("pool presence: values reflect the pool proportion", {
  # make_pool() has 3 trees; 2 recover the (A,B) clade, 1 does not
  result <- node_presence_matrix(make_backbone(), list(pool = make_pool()))
  # every cell is a multiple of 1/3
  vals <- result$presence[, "pool"]
  expect_true(all(abs(vals * 3 - round(vals * 3)) < 1e-9))
})

test_that("pool and single tree can be mixed in one call", {
  trees <- list(
    single = ape::read.tree(text = "(((A,B),C),(D,E));"),
    pool   = make_pool()
  )
  result <- node_presence_matrix(make_backbone(), trees)

  expect_equal(ncol(result$presence), 2L)
  expect_true(all(result$presence[, "single"] %in% c(0, 1)))
})


# ---- support matrix ---------------------------------------------------------

test_that("support: present clade has a numeric value", {
  result <- node_presence_matrix(make_supported(), list(t1 = make_supported()))

  vals    <- result$support[, "t1"]
  present <- vals[!is.na(vals)]
  expect_true(length(present) > 0)
  expect_true(all(present > 0))
})

test_that("support: absent clade is NA", {
  trees  <- list(diff = ape::read.tree(text = "((A,(B,C)),(D,E));"))
  result <- node_presence_matrix(make_backbone(), trees)
  expect_true(any(is.na(result$support)))
})

test_that("support: a tree with no labels gives all NA", {
  # make_backbone() carries no node labels. Clades are present (topology
  # matches) but there is no value to read.
  result <- node_presence_matrix(make_backbone(), list(t1 = make_backbone()))
  expect_true(all(is.na(result$support)))
})

test_that("support: presence and support disagree on absent clades", {
  # Where a clade is absent, presence is 0 but support is NA. This is the whole
  # reason two matrices exist.
  trees  <- list(diff = ape::read.tree(text = "((A,(B,C)),(D,E));"))
  result <- node_presence_matrix(make_backbone(), trees)

  absent <- result$presence[, "diff"] == 0
  expect_true(all(is.na(result$support[absent, "diff"])))
})

test_that("support_col = 2 reads the second value", {
  bb <- make_compound_support()

  s1 <- node_presence_matrix(bb, list(t1 = make_compound_support()),
                             support_col = 1)$support_1
  s2 <- node_presence_matrix(bb, list(t1 = make_compound_support()),
                             support_col = 2)$support_1

  v1 <- s1[!is.na(s1)]
  v2 <- s2[!is.na(s2)]
  if (length(v1) > 0 && length(v2) > 0) {
    expect_false(identical(v1, v2))
  }
})

test_that("support_col = 3 reads the third value", {
  bb     <- make_triple_support()
  result <- node_presence_matrix(bb, list(t1 = make_triple_support()),
                                 support_col = 3)
  expect_false(all(is.na(result$support_1)))
  expect_false(all(is.na(result$support)))
})

test_that("support_col = 3 gives NA when the tree has fewer values", {
  # Single value per node, so column 3 does not exist. NA, not an error.
  bb     <- make_supported()
  result <- node_presence_matrix(bb, list(t1 = make_supported()),
                                 support_col = 3)
  expect_true(all(is.na(result$support_1)))
})


# ---- attributes -------------------------------------------------------------

test_that("node_id attribute matches backbone internal nodes", {
  bb     <- make_backbone()
  result <- node_presence_matrix(bb, make_comparisons())

  ntip     <- ape::Ntip(bb)
  expected <- (ntip + 1L):(ntip + bb$Nnode)
  expect_equal(attr(result, "node_id"), expected)
})

test_that("pool_sizes attribute records tree count per comparison", {
  trees  <- list(single = make_backbone(), pool = make_pool())
  result <- node_presence_matrix(make_backbone(), trees)
  sizes  <- attr(result, "pool_sizes")

  expect_equal(sizes[["single"]], 1L)
  expect_equal(sizes[["pool"]], 3L)
})

test_that("no phylorug_mode attribute (both matrices are always returned)", {
  result <- node_presence_matrix(make_backbone(), make_comparisons())
  expect_null(attr(result, "phylorug_mode"))
})
# ---- multi-column behaviour -------------------------------------------------
test_that("support_col = c(1,2) returns support_1 and support_2", {
  bb     <- make_compound_support()
  result <- node_presence_matrix(bb, list(t1 = make_compound_support()),
                                 support_col = c(1, 2))
  expect_named(result, c("presence", "support_1", "support_2"))
  expect_true(is.matrix(result$support_1))
  expect_true(is.matrix(result$support_2))
})

test_that("support_1 and support_2 differ when node labels have two fields", {
  bb     <- make_compound_support()
  result <- node_presence_matrix(bb, list(t1 = make_compound_support()),
                                 support_col = c(1, 2))
  v1 <- result$support_1[!is.na(result$support_1)]
  v2 <- result$support_2[!is.na(result$support_2)]
  expect_false(identical(v1, v2))
})
