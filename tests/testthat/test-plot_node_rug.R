# Tests for plot_node_rug()

# ---- helpers ----------------------------------------------------------------

make_test_tree <- function() {
  set.seed(1)
  ape::rtree(5)
}

# A support matrix shaped like node_presence_matrix() output:
# first column node_id, then one column per analysis.
make_support_mt <- function(tree, values = NULL) {
  n_tip  <- ape::Ntip(tree)
  n_node <- ape::Nnode(tree)
  node_ids <- (n_tip + 1):(n_tip + n_node)
  if (is.null(values)) {
    # two analyses with varying support, some NA
    a <- c(95, 80, NA, 100)[seq_len(n_node)]
    b <- c(1.0, 0.6, 0.9, NA)[seq_len(n_node)]
  } else {
    a <- values
    b <- values
  }
  mt <- cbind(node_id = node_ids, TreeA = a, TreeB = b)
  mt
}

# A presence matrix: every present cell is exactly 1, absent is NA.
make_presence_mt <- function(tree) {
  n_tip  <- ape::Ntip(tree)
  n_node <- ape::Nnode(tree)
  node_ids <- (n_tip + 1):(n_tip + n_node)
  a <- c(1, 1, NA, 1)[seq_len(n_node)]
  b <- c(1, NA, 1, 1)[seq_len(n_node)]
  cbind(node_id = node_ids, TreeA = a, TreeB = b)
}

# Draw a tree to an off-screen device and return its last_pp.
# Every test that draws cells needs a plotted tree first, because
# plot_node_rug reads node coordinates from the active device.
with_plotted_tree <- function(tree, code) {
  grDevices::pdf(NULL)            # null device: draws nothing to screen/disk
  on.exit(grDevices::dev.off())
  ape::plot.phylo(tree)
  force(code)
}

# ---- input validation -------------------------------------------------------

test_that("stops when rug_mt is not a matrix or data frame", {
  expect_error(
    plot_node_rug(
      rug_mt = "not_a_matrix",
      hues   = c("red", "blue"),
      cell_h = 0.1,
      cell_w = 0.1,
      n_cols = 2
    ),
    "must be a matrix or data frame"
  )
})

# ---- return value -----------------------------------------------------------

test_that("returns invisible NULL", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  with_plotted_tree(tree, {
    result <- plot_node_rug(
      rug_mt = rug_mt,
      hues   = c("#D85A30", "#185FA5"),
      cell_h = 0.1,
      cell_w = 0.1,
      n_cols = 2
    )
    expect_null(result)
  })
})

# ---- runs cleanly across the branches ---------------------------------------

test_that("runs without error in support mode", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  with_plotted_tree(tree, {
    expect_no_error(
      plot_node_rug(
        rug_mt = rug_mt,
        hues   = c("#D85A30", "#185FA5"),
        cell_h = 0.1,
        cell_w = 0.1,
        n_cols = 2
      )
    )
  })
})

test_that("runs without error with show_values = TRUE", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  with_plotted_tree(tree, {
    expect_no_error(
      plot_node_rug(
        rug_mt      = rug_mt,
        hues        = c("#D85A30", "#185FA5"),
        cell_h      = 0.1,
        cell_w      = 0.1,
        n_cols      = 2,
        show_values = TRUE
      )
    )
  })
})

test_that("runs without error in presence mode", {
  tree   <- make_test_tree()
  rug_mt <- make_presence_mt(tree)
  with_plotted_tree(tree, {
    expect_no_error(
      plot_node_rug(
        rug_mt = rug_mt,
        hues   = c("#D85A30", "#185FA5"),
        cell_h = 0.1,
        cell_w = 0.1,
        n_cols = 2
      )
    )
  })
})

test_that("show_values is ignored in presence mode (still runs cleanly)", {
  # presence cells are all 1, so the function forces show_values off itself
  tree   <- make_test_tree()
  rug_mt <- make_presence_mt(tree)
  with_plotted_tree(tree, {
    expect_no_error(
      plot_node_rug(
        rug_mt      = rug_mt,
        hues        = c("#D85A30", "#185FA5"),
        cell_h      = 0.1,
        cell_w      = 0.1,
        n_cols      = 2,
        show_values = TRUE   # asked for, but should be suppressed internally
      )
    )
  })
})

test_that("handles an all-NA column without error", {
  tree   <- make_test_tree()
  n_node <- ape::Nnode(tree)
  node_ids <- (ape::Ntip(tree) + 1):(ape::Ntip(tree) + n_node)
  rug_mt <- cbind(
    node_id = node_ids,
    TreeA   = rep(NA_real_, n_node),
    TreeB   = c(0.9, 0.5, 1.0, 0.7)[seq_len(n_node)]
  )
  with_plotted_tree(tree, {
    expect_no_error(
      plot_node_rug(
        rug_mt = rug_mt,
        hues   = c("#D85A30", "#185FA5"),
        cell_h = 0.1,
        cell_w = 0.1,
        n_cols = 2
      )
    )
  })
})

test_that("handles a single analysis (one-column rug)", {
  tree   <- make_test_tree()
  n_node <- ape::Nnode(tree)
  node_ids <- (ape::Ntip(tree) + 1):(ape::Ntip(tree) + n_node)
  rug_mt <- cbind(
    node_id = node_ids,
    TreeA   = c(0.9, 0.5, 1.0, 0.7)[seq_len(n_node)]
  )
  with_plotted_tree(tree, {
    expect_no_error(
      plot_node_rug(
        rug_mt = rug_mt,
        hues   = "#D85A30",
        cell_h = 0.1,
        cell_w = 0.1,
        n_cols = 1
      )
    )
  })
})

# ---- last_pp argument -------------------------------------------------------

test_that("uses a supplied last_pp instead of fetching from the device", {
  tree   <- make_test_tree()
  rug_mt <- make_support_mt(tree)
  with_plotted_tree(tree, {
    last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
    expect_no_error(
      plot_node_rug(
        rug_mt  = rug_mt,
        hues    = c("#D85A30", "#185FA5"),
        cell_h  = 0.1,
        cell_w  = 0.1,
        n_cols  = 2,
        last_pp = last_pp
      )
    )
  })
})
