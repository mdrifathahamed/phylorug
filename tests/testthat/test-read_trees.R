# Tests for read_trees()

# ---- helpers ----------------------------------------------------------------

# One Newick file, one NEXUS file, one tree each.
make_valid_trees <- function() {
  tmp <- tempfile("phylorug_valid_")
  dir.create(tmp)

  writeLines(
    "(((A,B),C),(D,E));",
    file.path(tmp, "newick_tree.tre")
  )

  writeLines(
    c("#NEXUS", "BEGIN TREES;",
      "TREE 1 = (((A,B),C),(D,E));",
      "END;"),
    file.path(tmp, "nexus_tree.nex")
  )

  tmp
}

# A NEXUS file wearing a .tre extension, to prove detection is by content.
make_nexus_as_tre <- function() {
  tmp <- tempfile("phylorug_nexus_tre_")
  dir.create(tmp)
  writeLines(
    c("#NEXUS", "BEGIN TREES;",
      "TREE 1 = (((A,B),C),(D,E));",
      "END;"),
    file.path(tmp, "nexus_as_tre.tre")
  )
  tmp
}

# A .tre file that is not a tree at all.
make_trees_with_corrupt <- function() {
  tmp <- tempfile("phylorug_corrupt_")
  dir.create(tmp)

  writeLines("(((A,B),C),(D,E));", file.path(tmp, "good.tre"))
  writeLines("(((A,B),C),(D,E);",  file.path(tmp, "corrupt.tre"))

  tmp
}

# One single-tree file and one file holding several equally optimal trees, as
# POY or TNT would write them.
make_pool_trees <- function(n = 3L) {
  tmp <- tempfile("phylorug_pool_")
  dir.create(tmp)

  writeLines("(((A,B),C),(D,E));", file.path(tmp, "single.tre"))

  writeLines(
    c("(((A,B),C),(D,E));",
      "((A,B),(C,(D,E)));",
      "((A,(B,C)),(D,E));")[seq_len(n)],
    file.path(tmp, "pool.tre")
  )

  tmp
}

# A NEXUS character matrix: valid NEXUS, no TREES block, not a tree file.
make_dir_with_matrix <- function() {
  tmp <- tempfile("phylorug_matrix_")
  dir.create(tmp)

  writeLines("(((A,B),C),(D,E));", file.path(tmp, "good.tre"))

  writeLines(
    c("#NEXUS", "BEGIN DATA;",
      "DIMENSIONS NTAX=2 NCHAR=4;",
      "FORMAT DATATYPE=DNA;",
      "MATRIX",
      "A ACGT",
      "B ACGA",
      ";", "END;"),
    file.path(tmp, "alignment.nex")
  )

  tmp
}


# ---- input validation -------------------------------------------------------

test_that("stops when directory does not exist", {
  expect_error(
    read_trees("path/that/does/not/exist"),
    "does not exist"
  )
})

test_that("stops when ext is not a character vector", {
  expect_error(
    read_trees(tempdir(), ext = 5),
    "non-empty character vector"
  )
  expect_error(
    read_trees(tempdir(), ext = character(0)),
    "non-empty character vector"
  )
})

test_that("stops when no matching files found", {
  tmp <- tempfile("phylorug_empty_")
  dir.create(tmp)
  expect_error(
    read_trees(tmp, ext = "xyz"),
    "No files matching extensions"
  )
})

test_that("stops when format is invalid", {
  expect_error(
    read_trees(tempdir(), format = "invalid_format"),
    "should be one of"
  )
})


# ---- output structure -------------------------------------------------------

test_that("returns a named list of phylo objects", {
  tmp    <- make_valid_trees()
  result <- read_trees(tmp, ext = "tre", verbose = FALSE)

  expect_type(result, "list")
  expect_true(all(vapply(result, inherits, logical(1), "phylo")))
})

test_that("list names are filenames with the extension removed", {
  tmp    <- make_valid_trees()
  result <- read_trees(tmp, verbose = FALSE)

  # Names feed the rug legend, so they must be clean.
  expect_setequal(names(result), c("newick_tree", "nexus_tree"))
  expect_false(any(grepl("\\.", names(result))))
})

test_that("stops when names collide after removing extensions", {
  tmp <- tempfile("phylorug_dup_")
  dir.create(tmp)
  writeLines("(((A,B),C),(D,E));", file.path(tmp, "iqtree.tre"))
  writeLines("(((A,B),C),(D,E));", file.path(tmp, "iqtree.nwk"))

  expect_error(
    read_trees(tmp, verbose = FALSE),
    "not unique after removing file extensions"
  )
})


# ---- pools ------------------------------------------------------------------

test_that("a file with one tree yields a phylo", {
  tmp    <- make_pool_trees()
  result <- read_trees(tmp, verbose = FALSE)

  expect_s3_class(result$single, "phylo")
  expect_false(inherits(result$single, "multiPhylo"))
})

test_that("a file with several trees yields a multiPhylo", {
  tmp    <- make_pool_trees(n = 3L)
  result <- read_trees(tmp, verbose = FALSE)

  expect_s3_class(result$pool, "multiPhylo")
  expect_length(result$pool, 3L)
})

test_that("pool_sizes attribute counts trees per analysis", {
  tmp    <- make_pool_trees(n = 3L)
  result <- read_trees(tmp, verbose = FALSE)
  sizes  <- attr(result, "pool_sizes")

  expect_type(sizes, "integer")
  expect_equal(sizes[["single"]], 1L)
  expect_equal(sizes[["pool"]], 3L)
  expect_setequal(names(sizes), names(result))
})

test_that("single-tree and pooled analyses can be mixed in one directory", {
  # This is Caira et al. (2013): the parsimony arm gives several equally
  # optimal trees per condition, the ML arm gives one.
  tmp    <- make_pool_trees(n = 3L)
  result <- read_trees(tmp, verbose = FALSE)

  classes <- vapply(result, function(x) class(x)[1], character(1))
  expect_setequal(classes, c("phylo", "multiPhylo"))
})

test_that("verbose reports pooled analyses by name and size", {
  tmp <- make_pool_trees(n = 3L)
  expect_message(
    read_trees(tmp, verbose = TRUE),
    "pool \\(3 trees\\)"
  )
})


# ---- as_pool() and pool_size() ----------------------------------------------

test_that("as_pool wraps a phylo as a multiPhylo of length one", {
  tr   <- ape::read.tree(text = "(((A,B),C),(D,E));")
  pool <- phylorug:::as_pool(tr)

  expect_s3_class(pool, "multiPhylo")
  expect_length(pool, 1L)
  expect_s3_class(pool[[1]], "phylo")
})

test_that("as_pool passes a multiPhylo through unchanged", {
  trs  <- ape::read.tree(text = c("(((A,B),C),(D,E));", "((A,B),(C,(D,E)));"))
  pool <- phylorug:::as_pool(trs)

  expect_s3_class(pool, "multiPhylo")
  expect_length(pool, 2L)
})

test_that("as_pool rejects anything that is not a tree", {
  expect_error(phylorug:::as_pool(list(1, 2)), "Expected a `phylo`")
  expect_error(phylorug:::as_pool("not a tree"), "Expected a `phylo`")
})

test_that("pool_size returns 1 for a phylo, not its component count", {
  # A phylo is a list of 3-4 components, so length() would be wrong here.
  tr <- ape::read.tree(text = "(((A,B),C),(D,E));")
  expect_identical(phylorug:::pool_size(tr), 1L)
})

test_that("pool_size returns the tree count for a multiPhylo", {
  trs <- ape::read.tree(text = c("(((A,B),C),(D,E));", "((A,B),(C,(D,E)));"))
  expect_identical(phylorug:::pool_size(trs), 2L)
})


# ---- format detection -------------------------------------------------------

test_that("auto detects newick from file content", {
  tmp    <- make_valid_trees()
  result <- read_trees(tmp, ext = "tre", format = "auto", verbose = FALSE)

  expect_s3_class(result$newick_tree, "phylo")
})

test_that("auto detects nexus in a file with a .tre extension", {
  # Detection is by content, not by extension.
  tmp    <- make_nexus_as_tre()
  result <- read_trees(tmp, ext = "tre", format = "auto", verbose = FALSE)

  expect_s3_class(result$nexus_as_tre, "phylo")
})

test_that("auto detects nexus from #NEXUS content", {
  tmp    <- make_valid_trees()
  result <- read_trees(tmp, ext = "nex", format = "auto", verbose = FALSE)

  expect_s3_class(result$nexus_tree, "phylo")
})

test_that("forced newick format reads a newick file", {
  tmp    <- make_valid_trees()
  result <- read_trees(tmp, ext = "tre", format = "newick", verbose = FALSE)

  expect_s3_class(result$newick_tree, "phylo")
})

test_that("a single-tree NEXUS file yields a phylo, not a multiPhylo", {
  # ape::read.nexus() returns a multiPhylo even for one tree. Without the
  # unwrap, the return class would depend on the file format.
  tmp    <- make_valid_trees()
  result <- read_trees(tmp, ext = "nex", verbose = FALSE)

  expect_s3_class(result$nexus_tree, "phylo")
  expect_false(inherits(result$nexus_tree, "multiPhylo"))
})

test_that("detect_format stops on an empty file", {
  tmp <- tempfile()
  file.create(tmp)
  expect_error(
    phylorug:::detect_format(tmp),
    "appears to be empty"
  )
})

test_that("detect_format finds #NEXUS after a leading blank line or comment", {
  # MrBayes and PAUP* routinely emit these.
  tmp <- tempfile()
  writeLines(
    c("", "[ written by PAUP* ]", "#NEXUS", "BEGIN TREES;",
      "TREE 1 = (((A,B),C),(D,E));", "END;"),
    tmp
  )
  expect_equal(phylorug:::detect_format(tmp), "nexus")
})

test_that("detect_format returns none for a NEXUS file with no TREES block", {
  tmp <- tempfile()
  writeLines(
    c("#NEXUS", "BEGIN DATA;", "DIMENSIONS NTAX=2 NCHAR=4;",
      "MATRIX", "A ACGT", "B ACGA", ";", "END;"),
    tmp
  )
  expect_equal(phylorug:::detect_format(tmp), "none")
})

test_that("detect_format returns none for a file with no tree-like content", {
  tmp <- tempfile()
  writeLines(c("some", "plain", "text"), tmp)
  expect_equal(phylorug:::detect_format(tmp), "none")
})

test_that("detect_format defaults to newick", {
  tmp <- tempfile()
  writeLines("((A:0.1,B:0.2):0.3,C:0.4);", tmp)
  expect_equal(phylorug:::detect_format(tmp), "newick")
})


# ---- skipping non-tree files ------------------------------------------------

test_that("a NEXUS character matrix is skipped, not an error", {
  # An alignment sitting beside the trees is normal, and is not a broken
  # analysis. This is the concatenated_garli.nex case.
  tmp    <- make_dir_with_matrix()
  result <- read_trees(tmp, verbose = FALSE)

  expect_named(result, "good")
  expect_length(result, 1L)
})

test_that("skipped files are reported when verbose", {
  tmp <- make_dir_with_matrix()
  expect_message(
    read_trees(tmp, verbose = TRUE),
    "Skipped 1 file\\(s\\) containing no phylogenetic tree: alignment.nex"
  )
})

test_that("stops when files match ext but none contain a tree", {
  tmp <- tempfile("phylorug_nomatrix_")
  dir.create(tmp)
  writeLines(
    c("#NEXUS", "BEGIN DATA;", "MATRIX", "A ACGT", ";", "END;"),
    file.path(tmp, "alignment.nex")
  )

  expect_error(
    read_trees(tmp, verbose = FALSE),
    "none contained a tree"
  )
})


# ---- error handling ---------------------------------------------------------

test_that("errors on a tree file that cannot be parsed", {
  # A malformed .tre file is a broken analysis, and must be fatal -- unlike a
  # file that was never a tree at all.
  tmp <- make_trees_with_corrupt()
  expect_error(
    read_trees(tmp, ext = "tre", verbose = FALSE),
    "Could not parse tree file corrupt.tre"
  )
})


# ---- verbosity --------------------------------------------------------------

test_that("verbose = FALSE produces no messages", {
  tmp <- make_valid_trees()
  expect_silent(
    read_trees(tmp, ext = "nex", verbose = FALSE)
  )
})

test_that("verbose = TRUE reports analyses and trees separately", {
  tmp <- make_valid_trees()
  expect_message(
    read_trees(tmp, verbose = TRUE),
    "Read 2 analyses \\(2 trees\\)"
  )
})
