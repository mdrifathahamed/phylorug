# Tests for plot_phylorug()
#
# Drawing needs a graphics device, so these open a throwaway PDF and check that
# each mode runs without error. Visual correctness is judged by eye on real
# figures, not here.

# ---- helpers ----------------------------------------------------------------

make_backbone <- function() {
  ape::read.tree(text = "(((A,B),C),(D,E));")
}

# A node_presence_matrix-style list: presence + support_1 + support_2,
# node ids as rownames. Matches the real builder's naming.
make_npm <- function() {
  nodes <- as.character(6:9)
  cols  <- c("iqtree", "astral")

  presence <- matrix(
    c(1, 0, 1, 1,
      1, 1, 0, 1),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_1 <- matrix(
    c(99, NA, 85, 70,
      0.99, 0.60, NA, 0.97),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_2 <- matrix(
    c(95, NA, 80, 65,
      NA, NA, NA, NA),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )

  out <- list(presence = presence, support_1 = support_1, support_2 = support_2)
  attr(out, "node_id")    <- 6:9
  attr(out, "pool_sizes") <- c(iqtree = 1L, astral = 1L)
  out
}

on_null_device <- function(code) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  })
  force(code)
}


# ---- input validation -------------------------------------------------------

test_that("stops when backbone is not a phylo object", {
  expect_error(
    plot_phylorug("not_a_tree", make_npm()),
    "must be a phylogenetic tree"
  )
})

test_that("stops when npm has no presence element", {
  expect_error(
    plot_phylorug(make_backbone(), list(support_1 = matrix(1, 2, 2))),
    "must be the list returned by"
  )
  expect_error(
    plot_phylorug(make_backbone(), matrix(1, 2, 2)),
    "must be the list returned by"
  )
})

test_that("support mode requires support_type", {
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), mode = "support")
    ),
    "requires `support_type`"
  )
})

test_that("support mode errors when the requested support column is absent", {
  npm <- make_npm()   # has support_1 and support_2 only
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone(), npm,
                    mode = "support",
                    support_col = 3,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    ),
    "has no `support_3`"
  )
})

test_that("stops when the grid cannot hold all trees", {
  expect_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), n_rows = 1, n_cols = 1)
    ),
    "cannot hold"
  )
})


# ---- presence mode ----------------------------------------------------------

test_that("presence mode returns invisible NULL", {
  result <- on_null_device(
    plot_phylorug(make_backbone(), make_npm(), mode = "presence")
  )
  expect_null(result)
})

test_that("presence mode runs by default", {
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), make_npm()))
  )
})

test_that("presence mode ignores support_type", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "presence",
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})


# ---- support mode -----------------------------------------------------------

test_that("support mode runs with support_col = 1", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support",
                    support_col = 1,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})

test_that("support mode runs with support_col = 2", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support",
                    support_col = 2,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})

test_that("support mode accepts custom thresholds", {
  custom <- list(
    sh_alrt = c(very_high = 98, high = 80, moderate = 50),
    lpp     = c(very_high = 0.99, high = 0.95, moderate = 0.5)
  )
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support",
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"),
                    thresholds = custom)
    )
  )
})


# ---- options ----------------------------------------------------------------

test_that("runs with legend = FALSE", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), legend = FALSE)
    )
  )
})

test_that("runs with dot_unanimous = FALSE", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), dot_unanimous = FALSE)
    )
  )
})

test_that("runs with a user-fixed grid shape", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), n_rows = 1, n_cols = 2)
    )
  )
})

test_that("runs with the inside rug position", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), rug_position = "inside")
    )
  )
})

test_that("passes ... through to plot.phylo", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), edge.width = 2)
    )
  )
})


# ---- edge cases -------------------------------------------------------------

test_that("runs when every node is variable (no unanimous dots)", {
  npm <- make_npm()
  npm$presence[, 1] <- c(1, 0, 1, 0)
  npm$presence[, 2] <- c(0, 1, 0, 1)

  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), npm))
  )
})

test_that("runs when every node is unanimous (dots only, no rug)", {
  npm <- make_npm()
  npm$presence[] <- 1

  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), npm))
  )
})
