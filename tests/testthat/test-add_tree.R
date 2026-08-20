# tests/testthat/test-add_tree.R
# Tests for add_tree()

# ---- Helpers ----------------------------------------------------------------

make_backbone <- function() {
  ape::read.tree(text = "(((A:1,B:1)90:1,C:2)80:1,(D:2,E:2)100:1);")
}

make_npm <- function() {
  bb <- make_backbone()
  t1 <- ape::read.tree(text = "(((A:1,B:1)85:1,C:2)75:1,(D:2,E:2)95:1);")
  t2 <- ape::read.tree(text = "((A:1,(B:1,C:1)60:1)70:1,(D:2,E:2)99:1);")
  node_presence_matrix(bb, list(t1 = t1, t2 = t2))
}

# =============================================================================
# 1. INPUT VALIDATION
# =============================================================================

test_that("add_tree errors on invalid npm", {
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree("not_a_list", bb, new, name = "x"), "npm")
  expect_error(add_tree(list(wrong = 1), bb, new, name = "x"), "npm")
})

test_that("add_tree errors on non-phylo backbone", {
  npm <- make_npm()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree(npm, "not_a_tree", new, name = "x"), "backbone")
})

test_that("add_tree errors on non-phylo new_tree", {
  npm <- make_npm()
  bb  <- make_backbone()
  expect_error(add_tree(npm, bb, "not_a_tree", name = "x"), "new_tree")
  expect_error(add_tree(npm, bb, 42, name = "x"), "new_tree")
})

test_that("add_tree errors on unrooted backbone", {
  npm <- make_npm()
  bb  <- ape::unroot(make_backbone())
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree(npm, bb, new, name = "x"), "unrooted")
})

test_that("add_tree errors on unrooted new_tree", {
  npm <- make_npm()
  bb  <- make_backbone()
  unr <- ape::unroot(ape::read.tree(text = "(A,B,(C,(D,E)));"))
  expect_error(add_tree(npm, bb, unr, name = "x"), "unrooted")
})

test_that("add_tree errors on invalid support_col", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree(npm, bb, new, name = "x", support_col = 4), "support_col")
  expect_error(add_tree(npm, bb, new, name = "x", support_col = 0), "support_col")
  expect_error(add_tree(npm, bb, new, name = "x", support_col = "a"), "support_col")
})

# =============================================================================
# 2. NAME VALIDATION
# =============================================================================

test_that("add_tree errors when name is missing", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree(npm, bb, new), "name")
})

test_that("add_tree errors on empty name", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree(npm, bb, new, name = ""), "non-empty")
})

test_that("add_tree errors on non-character name", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree(npm, bb, new, name = 42), "character")
})

test_that("add_tree errors on duplicate name", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_error(add_tree(npm, bb, new, name = "t1"), "already exists")
})

# =============================================================================
# 3. IDENTICAL TAXA — core functionality
# =============================================================================

test_that("add_tree appends one column to presence matrix", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  result <- add_tree(npm, bb, new, name = "new")
  expect_equal(ncol(result$presence), 3L)
  expect_true("new" %in% colnames(result$presence))
})

test_that("add_tree appends one column to support matrix", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B)88,C)72,(D,E)97);")
  result <- add_tree(npm, bb, new, name = "new")
  expect_equal(ncol(result$support_1), 3L)
  expect_true("new" %in% colnames(result$support_1))
})

test_that("add_tree does not modify existing columns", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  result <- add_tree(npm, bb, new, name = "new")
  # Existing columns unchanged
  expect_identical(result$presence[, "t1"], npm$presence[, "t1"])
  expect_identical(result$presence[, "t2"], npm$presence[, "t2"])
  expect_identical(result$support_1[, "t1"], npm$support_1[, "t1"])
  expect_identical(result$support_1[, "t2"], npm$support_1[, "t2"])
})

test_that("add_tree records correct presence for identical topology", {
  npm <- make_npm()
  bb  <- make_backbone()
  # Same topology as backbone — every clade present
  same <- ape::read.tree(text = "(((A,B),C),(D,E));")
  result <- add_tree(npm, bb, same, name = "same")
  new_col <- result$presence[, "same"]
  # Every clade should be 1 (present)
  expect_true(all(new_col == 1))
})

test_that("add_tree records correct presence for different topology", {
  npm <- make_npm()
  bb  <- make_backbone()
  # Different topology — (A,(B,C)) instead of ((A,B),C)
  diff <- ape::read.tree(text = "((A,(B,C)),(D,E));")
  result <- add_tree(npm, bb, diff, name = "diff")
  new_col <- result$presence[, "diff"]
  # Some clades should be 0 (absent)
  expect_true(any(new_col == 0))
  # Root clade and (D,E) should still be present
  expect_true(any(new_col == 1))
})

test_that("add_tree extracts support values correctly", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B)88,C)72,(D,E)97);")
  result <- add_tree(npm, bb, new, name = "new")
  present_rows <- which(result$presence[, "new"] == 1)
  if (length(present_rows) > 0) {
    support_vals <- result$support_1[present_rows, "new"]
    expect_true(any(!is.na(support_vals)))
  }
})

test_that("add_tree sets support to NA where clade is absent", {
  npm <- make_npm()
  bb  <- make_backbone()
  diff <- ape::read.tree(text = "((A,(B,C)70)85,(D,E)99);")
  result <- add_tree(npm, bb, diff, name = "diff")
  absent <- which(result$presence[, "diff"] == 0)
  if (length(absent) > 0) {
    expect_true(all(is.na(result$support_1[, "diff"][absent])))
  }
})

# =============================================================================
# 4. SUPERSET — extra taxa pruned internally
# =============================================================================
test_that("add_tree prunes extra taxa and matches clades correctly", {
  npm <- make_npm()
  bb  <- make_backbone()
  super <- ape::read.tree(text = "((((A,B),C),X),(D,E));")
  expect_message(
    add_tree(npm, bb, super, name = "super"),
    "pruned internally"
  )
  result <- suppressMessages(add_tree(npm, bb, super, name = "super2"))
  expect_equal(ncol(result$presence), 3L)
  new_col <- result$presence[, "super2"]
  expect_true(any(new_col == 1, na.rm = TRUE))
})
test_that("superset message lists extra taxa", {
  npm <- make_npm()
  bb  <- make_backbone()
  super <- ape::read.tree(text = "((((A,B),C),X),(D,E));")
  expect_message(
    add_tree(npm, bb, super, name = "super"),
    "X"
  )
})
test_that("superset does not produce false absences", {
  npm <- make_npm()
  bb  <- make_backbone()
  super <- ape::read.tree(text = "(((A,B),C),(D,E,X));")
  expect_message(
    add_tree(npm, bb, super, name = "super"),
    "pruned"
  )
  result <- suppressMessages(add_tree(npm, bb, super, name = "super2"))
  new_col <- result$presence[, "super2"]
  expect_true(sum(new_col == 1, na.rm = TRUE) >= 2)
})

# =============================================================================
# 5. MISSING TAXA — affected clades get NA
# =============================================================================

test_that("add_tree marks clades with missing taxa as NA", {
  npm <- make_npm()
  bb  <- make_backbone()
  gappy <- ape::read.tree(text = "((A,B),(D,E));")
  gappy <- ape::root(gappy, "A", resolve.root = TRUE)
  expect_message(
    add_tree(npm, bb, gappy, name = "gappy"),
    "missing"
  )
  result <- suppressMessages(add_tree(npm, bb, gappy, name = "gappy2"))
  new_col <- result$presence[, "gappy2"]
  expect_true(any(is.na(new_col)))
  expect_true(any(!is.na(new_col)))
})

test_that("missing taxa message lists the taxa", {
  npm <- make_npm()
  bb  <- make_backbone()
  gappy <- ape::read.tree(text = "((A,B),(D,E));")
  gappy <- ape::root(gappy, "A", resolve.root = TRUE)
  expect_message(
    add_tree(npm, bb, gappy, name = "gappy"),
    "C"
  )
})
test_that("missing taxa message suggests prune_to_shared", {
  npm <- make_npm()
  bb  <- make_backbone()
  gappy <- ape::read.tree(text = "((A,B),(D,E));")
  gappy <- ape::root(gappy, "A", resolve.root = TRUE)
  expect_message(
    add_tree(npm, bb, gappy, name = "gappy"),
    "prune_to_shared"
  )
})

test_that("support is NA where taxa are missing", {
  npm <- make_npm()
  bb  <- make_backbone()
  gappy <- ape::read.tree(text = "((A,B)90,(D,E)95);")
  gappy <- ape::root(gappy, "A", resolve.root = TRUE)
  expect_message(
    add_tree(npm, bb, gappy, name = "gappy"),
    "missing"
  )
  result <- suppressMessages(add_tree(npm, bb, gappy, name = "gappy2"))
  na_rows <- which(is.na(result$presence[, "gappy2"]))
  if (length(na_rows) > 0) {
    expect_true(all(is.na(result$support_1[, "gappy2"][na_rows])))
  }
})

# =============================================================================
# 6. MULTIPHY POOL
# =============================================================================

test_that("add_tree handles multiPhylo pool", {
  npm <- make_npm()
  bb  <- make_backbone()
  pool <- list(
    ape::read.tree(text = "(((A,B),C),(D,E));"),
    ape::read.tree(text = "((A,(B,C)),(D,E));")
  )
  class(pool) <- "multiPhylo"
  result <- add_tree(npm, bb, pool, name = "pool")
  expect_equal(ncol(result$presence), 3L)
  new_col <- result$presence[, "pool"]
  # Some values should be between 0 and 1 (partial recovery)
  has_partial <- any(new_col > 0 & new_col < 1, na.rm = TRUE)
  # (D,E) recovered by both — should be 1
  has_full <- any(new_col == 1, na.rm = TRUE)
  expect_true(has_partial || has_full)
})

# =============================================================================
# 7. ATTRIBUTES
# =============================================================================

test_that("add_tree updates pool_sizes attribute", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  result <- add_tree(npm, bb, new, name = "new")
  ps <- attr(result, "pool_sizes")
  expect_true("new" %in% names(ps))
  expect_equal(ps[["new"]], 1L)
})

test_that("add_tree stores support_type when supplied", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B)88,C)72,(D,E)97);")
  result <- add_tree(npm, bb, new, name = "new",
                     support_type = "ufboot")
  st <- attr(result, "support_type")
  expect_true("new" %in% names(st))
  expect_equal(st[["new"]], "ufboot")
})

test_that("add_tree does not add support_type when NULL", {
  npm <- make_npm()
  bb  <- make_backbone()
  new <- ape::read.tree(text = "(((A,B),C),(D,E));")
  result <- add_tree(npm, bb, new, name = "new")
  st <- attr(result, "support_type")
  expect_false("new" %in% names(st))
})

test_that("add_tree appends to existing support_type", {
  bb  <- make_backbone()
  t1  <- ape::read.tree(text = "(((A,B)85,C)75,(D,E)95);")
  npm <- node_presence_matrix(bb, list(t1 = t1),
                              support_type = c(t1 = "ufboot"))
  new <- ape::read.tree(text = "(((A,B)88,C)72,(D,E)97);")
  result <- add_tree(npm, bb, new, name = "new",
                     support_type = "lpp")
  st <- attr(result, "support_type")
  expect_equal(length(st), 2L)
  expect_equal(st[["t1"]], "ufboot")
  expect_equal(st[["new"]], "lpp")
})

# =============================================================================
# 8. MULTIPLE SUPPORT COLUMNS
# =============================================================================

test_that("add_tree extracts multiple support columns", {
  bb  <- make_backbone()
  t1  <- ape::read.tree(text = "(((A,B)85/90,C)75/80,(D,E)95/99);")
  npm <- node_presence_matrix(bb, list(t1 = t1), support_col = c(1, 2))
  new <- ape::read.tree(text = "(((A,B)88/92,C)72/78,(D,E)97/100);")
  result <- add_tree(npm, bb, new, name = "new", support_col = c(1, 2))
  expect_true("support_1" %in% names(result))
  expect_true("support_2" %in% names(result))
  expect_equal(ncol(result$support_1), 2L)
  expect_equal(ncol(result$support_2), 2L)
})

# =============================================================================
# 9. CHAINING — add_tree called multiple times
# =============================================================================

test_that("add_tree can be chained", {
  npm <- make_npm()
  bb  <- make_backbone()
  t3  <- ape::read.tree(text = "(((A,B),C),(D,E));")
  t4  <- ape::read.tree(text = "((A,(B,C)),(D,E));")
  npm <- add_tree(npm, bb, t3, name = "t3")
  npm <- add_tree(npm, bb, t4, name = "t4")
  expect_equal(ncol(npm$presence), 4L)
  expect_equal(colnames(npm$presence), c("t1", "t2", "t3", "t4"))
})

# =============================================================================
# 10. EDGE CASES
# =============================================================================

test_that("add_tree works with sample_trees data", {
  backbone <- sample_trees[["70p_uce"]]
  others   <- sample_trees[c("70p_partition_entropy", "70p_ghost")]
  npm      <- node_presence_matrix(backbone, others)
  result   <- add_tree(npm, backbone,
                       sample_trees[["70p_ASTRAL_uce"]],
                       name = "ASTRAL_uce")
  expect_equal(ncol(result$presence), 3L)
  expect_true("ASTRAL_uce" %in% colnames(result$presence))
})

test_that("add_tree handles tree with no support labels", {
  npm <- make_npm()
  bb  <- make_backbone()
  bare <- ape::read.tree(text = "(((A,B),C),(D,E));")
  result <- add_tree(npm, bb, bare, name = "bare")
  # Support should be all NA for bare tree
  expect_true(all(is.na(result$support_1[, "bare"]) |
                    result$presence[, "bare"] == 0))
})
