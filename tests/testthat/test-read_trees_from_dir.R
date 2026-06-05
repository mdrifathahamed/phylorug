# Tests for read_trees_from_dir()

# ---- helper -----------------------------------------------------------------
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

# Helper with corrupt file -- only for error handling tests

make_trees_with_corrupt <- function() {
  tmp <- tempfile("phylorug_corrupt_")
  dir.create(tmp)

  writeLines("(((A,B),C),(D,E));", file.path(tmp, "good.tre"))
  writeLines("THIS IS NOT A TREE", file.path(tmp, "corrupt.tre"))

  tmp
}

# ---- input validation -------------------------------------------------------

test_that("stops when directory does not exist", {
  expect_error(
    read_trees_from_dir("path/that/does/not/exist"),
    "does not exist"
  )
})

test_that("stops when no matching files found", {
  tmp <- tempdir()
  expect_error(
    read_trees_from_dir(tmp, ext = "xyz"),
    "No files matching extensions"
  )
})

test_that("stops when format is invalid", {
  expect_error(
    read_trees_from_dir(tempdir(), format = "invalid_format"),
    "should be one of"
  )
})

# ---- output structure -------------------------------------------------------

test_that("returns a named list of phylo objects", {
  tmp <- make_valid_trees()
  result <- read_trees_from_dir(tmp, ext = "tre", verbose = FALSE)

  expect_type(result, "list")
  expect_true(
    all(vapply(
      result[!vapply(result, is.null, logical(1))],
      inherits,
      logical(1),
      "phylo"
    ))
  )
})


test_that("list names matches filenames in the directory", {
  tmp <- make_valid_trees()
  result <- read_trees_from_dir(tmp, ext = "tre", verbose = FALSE)

  expect_equal(
    names(result),
    basename(list.files(tmp, pattern = "\\.tre$"))
  )
})

# ---- format detection -------------------------------------------------------

test_that("auto detects newick from file content", {
  tmp <- make_valid_trees()
  result <- read_trees_from_dir(
    tmp,
    ext     = "tre",
    format  = "auto",
    verbose = FALSE
  )

  expect_true(inherits(result$newick_tree.tre, "phylo"))
})

test_that("auto detects nexus from .tre file with #NEXUS content", {
  tmp <- make_nexus_as_tre() # helper does the setup
  result <- read_trees_from_dir(
    tmp,
    ext     = "tre",
    format  = "auto",
    verbose = FALSE
  )
  expect_true(inherits(result$nexus_as_tre.tre, "phylo"))
})

test_that("auto detects nexus from #NEXUS content", {
  tmp <- make_valid_trees()
  result <- read_trees_from_dir(
    tmp,
    ext     = "nex",
    format  = "auto",
    verbose = FALSE
  )

  expect_true(inherits(result$nexus_tree.nex, "phylo"))
})

test_that("forced newick format reads newick file correctly", {
  tmp <- make_valid_trees()
  result <- read_trees_from_dir(
    tmp,
    ext     = "tre",
    format  = "newick",
    verbose = FALSE
  )

  expect_true(inherits(result$newick_tree.tre, "phylo"))
})
#-----------------test for the detet_format internal helper -------------------
test_that("detect_format stops on empty file", {
  tmp <- tempfile()
  file.create(tmp)
  expect_error(
    phylorug:::detect_format(tmp),
    "appears to be empty"
  )
})

test_that("detect_format identifies NHX as newick", {
  tmp <- tempfile()
  writeLines(
    "((A:0.1,B:0.2)[&&NHX:bs=100]:0.3,C:0.4);",
    tmp
  )
  expect_equal(phylorug:::detect_format(tmp), "newick")
})
# ---- error handling ---------------------------------------------------------

test_that("errors when a file cannot be read", {
  tmp <- make_trees_with_corrupt()
  expect_error(
    read_trees_from_dir(
      tmp,
      ext     = "tre",
      verbose = FALSE
    ),
    "file\\(s\\) could not be read"
  )
})

test_that("verbose = FALSE produces no messages", {
  tmp <- make_valid_trees()
  expect_silent(
    read_trees_from_dir(
      tmp,
      ext     = "nex",
      verbose = FALSE
    )
  )
})
test_that("verbose = TRUE prints a message", {
  tmp <- make_valid_trees()
  expect_message(
    read_trees_from_dir(tmp, ext = "nex", verbose = TRUE),
    "Read"
  )
})
