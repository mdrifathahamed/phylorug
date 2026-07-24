# Tests for plot_phylorug()

# ---- helpers ----------------------------------------------------------------
make_backbone <- function() {
  ape::read.tree(text = "(((A,B),C),(D,E));")
}

make_npm <- function() {
  nodes <- as.character(6:9)
  cols  <- c("iqtree", "astral")
  presence <- matrix(
    c(1, 0, 1, 1,  1, 1, 0, 1),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_1 <- matrix(
    c(99, NA, 85, 70,  0.99, 0.60, NA, 0.97),
    nrow = 4, ncol = 2, dimnames = list(nodes, cols)
  )
  support_2 <- matrix(
    c(95, NA, 80, 65,  NA, NA, NA, NA),
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
  on.exit({ grDevices::dev.off(); unlink(tmp) })
  force(code)
}

# ---- input validation -------------------------------------------------------
test_that("stops when backbone is not a phylo object", {
  expect_error(
    plot_phylorug("not_a_tree", make_npm()),
    "phylo"
  )
})

test_that("stops when npm has no presence element", {
  expect_error(
    plot_phylorug(make_backbone(), list(support_1 = matrix(1, 2, 2))),
    "presence"
  )
  expect_error(
    plot_phylorug(make_backbone(), matrix(1, 2, 2)),
    "presence"
  )
})

# ---- presence mode ----------------------------------------------------------
test_that("presence mode runs without error", {
  expect_no_error(
    on_null_device(plot_phylorug(make_backbone(), make_npm()))
  )
})

test_that("presence mode returns invisible NULL", {
  result <- on_null_device(
    plot_phylorug(make_backbone(), make_npm(), mode = "presence")
  )
  expect_null(result)
})

test_that("presence mode ignores support_type silently", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "presence",
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})

# ---- support mode -----------------------------------------------------------
test_that("support mode runs with support_idx = 1", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support", support_idx = 1,
                    support_type = c(iqtree = "sh_alrt", astral = "lpp"))
    )
  )
})

test_that("support mode runs with support_idx = 2", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(),
                    mode = "support", support_idx = 2,
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

test_that("runs with dot_identical= FALSE", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), dot_identical= FALSE)
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

test_that("runs with inside rug position", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), rug_position = "inside")
    )
  )
})

test_that("passes ... through to plot.phylo", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), edge.width = 2, cex = 0.8)
    )
  )
})

# ---- file output ------------------------------------------------------------
test_that("file argument creates a PDF and returns the path", {
  tmp <- tempfile(fileext = ".pdf")
  result <- plot_phylorug(make_backbone(), make_npm(), file = tmp)
  expect_true(file.exists(tmp))
  expect_equal(result, tmp)
  unlink(tmp)
})

test_that("file argument creates a PNG", {
  tmp <- tempfile(fileext = ".png")
  result <- plot_phylorug(make_backbone(), make_npm(), file = tmp)
  expect_true(file.exists(tmp))
  unlink(tmp)
})

# ---- cell_scale and dot_cex ------------------------------------------------
test_that("cell_scale adjusts rug size", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), cell_scale = 0.7)
    )
  )
})

test_that("dot_cex adjusts dot size", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), dot_cex = 1.0)
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

test_that("runs with include_backbone = TRUE", {
  expect_no_error(
    on_null_device(
      plot_phylorug(make_backbone(), make_npm(), include_backbone = TRUE)
    )
  )
})
