# Tests for plot_node_rug() and its helpers
#
# resolve_cell returns list(fill, pattern, border) where pattern is
# "none", "cross", or "dots". No $hatch field.
# bin_fill returns list(fill, pattern), not a bare string.

# ---- resolve_tier -----------------------------------------------------------
test_that("tier 1 when support is NULL", {
  expect_equal(resolve_tier(NULL, NULL), 1L)
  expect_equal(resolve_tier(NULL, c(t1 = "ufboot")), 1L)
})

test_that("tier 2 when support given but support_type is NULL", {
  m <- matrix(1, 2, 2)
  expect_equal(resolve_tier(m, NULL), 2L)
})

test_that("tier 3 when both support and support_type are given", {
  m <- matrix(1, 2, 2)
  expect_equal(resolve_tier(m, c(t1 = "ufboot")), 3L)
})

# ---- resolve_cell: not evaluable --------------------------------------------
test_that("NA presence gets cross pattern in every tier", {
  for (tier in 1:3) {
    cell <- resolve_cell(NA_real_, NA_real_, NA_character_, NULL, tier)
    expect_equal(cell$pattern, "cross")
    expect_equal(cell$fill, "white")
  }
})

# ---- resolve_cell: absent ---------------------------------------------------
test_that("absent clade is white with no pattern in tiers 1 and 2", {
  c1 <- resolve_cell(0, NA_real_, NA_character_, NULL, 1L)
  c2 <- resolve_cell(0, NA_real_, NA_character_, NULL, 2L)
  expect_equal(c1$fill, "white")
  expect_equal(c2$fill, "white")
  expect_equal(c1$pattern, "none")
  expect_equal(c2$pattern, "none")
})

test_that("absent clade is white with no pattern in tier 3", {
  c3 <- resolve_cell(0, NA_real_, "ufboot", NULL, 3L)
  expect_equal(c3$fill, "white")
  expect_equal(c3$pattern, "none")
})

# ---- resolve_cell: present, tiers 1 and 2 -----------------------------------
test_that("present clade is solid black in tiers 1 and 2", {
  c1 <- resolve_cell(1, NA_real_, NA_character_, NULL, 1L)
  expect_equal(c1$fill, grDevices::rgb(0, 0, 0))
  expect_equal(c1$pattern, "none")
})

test_that("pool proportion is grey by proportion in tiers 1 and 2", {
  c1 <- resolve_cell(0.5, NA_real_, NA_character_, NULL, 1L)
  expect_equal(c1$fill, grDevices::rgb(0.5, 0.5, 0.5))
})

# ---- resolve_cell: present, tier 3 ------------------------------------------
test_that("tier 3 present-but-no-support is solid black", {
  c3 <- resolve_cell(1, NA_real_, "ufboot", NULL, 3L)
  expect_equal(c3$fill, "black")
  expect_equal(c3$pattern, "none")
})

test_that("tier 3 shades a present clade by its binned support", {
  # UFBoot 99 -> very high -> bin 4
  c_high <- resolve_cell(1, 99, "ufboot", NULL, 3L)
  expect_equal(c_high$fill, bin_fill(4L)$fill)
  expect_equal(c_high$pattern, bin_fill(4L)$pattern)

  # UFBoot 40 -> low -> bin 1
  c_low <- resolve_cell(1, 40, "ufboot", NULL, 3L)
  expect_equal(c_low$fill, bin_fill(1L)$fill)
  expect_equal(c_low$pattern, bin_fill(1L)$pattern)
})

# ---- bin_support: the thresholds --------------------------------------------
test_that("bin_support returns NA for NA input", {
  expect_true(is.na(bin_support(NA_real_, "ufboot", NULL)))
})

test_that("UFBoot binning follows the default thresholds", {
  expect_equal(bin_support(96, "ufboot", NULL), 4L)
  expect_equal(bin_support(95, "ufboot", NULL), 4L)
  expect_equal(bin_support(70, "ufboot", NULL), 2L)
  expect_equal(bin_support(40, "ufboot", NULL), 1L)
})

test_that("SH-aLRT binning follows the default thresholds", {
  expect_equal(bin_support(99, "sh_alrt", NULL), 4L)
  expect_equal(bin_support(85, "sh_alrt", NULL), 3L)
  expect_equal(bin_support(60, "sh_alrt", NULL), 2L)
  expect_equal(bin_support(30, "sh_alrt", NULL), 1L)
})

test_that("ASTRAL LPP binning follows the default thresholds", {
  expect_equal(bin_support(0.99, "lpp", NULL), 4L)
  expect_equal(bin_support(0.96, "lpp", NULL), 3L)
  expect_equal(bin_support(0.70, "lpp", NULL), 2L)
  expect_equal(bin_support(0.30, "lpp", NULL), 1L)
})

test_that("the same numeric value bins differently by support type", {
  expect_equal(bin_support(85, "sh_alrt", NULL), 3L)
  expect_equal(bin_support(85, "ufboot",  NULL), 2L)
})

test_that("user thresholds override the defaults", {
  custom <- list(ufboot = c(very_high = 99, high = 90, moderate = 70))
  expect_equal(bin_support(95, "ufboot", custom), 3L)
  expect_equal(bin_support(80, "ufboot", custom), 2L)
  expect_equal(bin_support(60, "ufboot", custom), 1L)
})

test_that("unknown support type falls back to a percentage scale", {
  expect_equal(bin_support(96, "mystery", NULL), 4L)
  expect_equal(bin_support(85, "mystery", NULL), 3L)
})

# ---- bin_fill ---------------------------------------------------------------
test_that("bin_fill returns a list with fill and pattern", {
  result <- bin_fill(4L)
  expect_true(is.list(result))
  expect_true("fill" %in% names(result))
  expect_true("pattern" %in% names(result))
})

test_that("bin_fill maps each bin to a distinct fill colour", {
  fills <- vapply(1:4, function(b) bin_fill(b)$fill, character(1))
  expect_length(unique(fills), 4L)
})

test_that("bin_fill returns black fill for NA (present but unquantified)", {
  result <- bin_fill(NA_integer_)
  expect_equal(result$fill, "black")
  expect_equal(result$pattern, "none")
})

test_that("bin 1 (low support) uses dots pattern", {
  result <- bin_fill(1L)
  expect_equal(result$pattern, "dots")
})

test_that("bins 2-4 use no pattern", {
  for (b in 2:4) {
    expect_equal(bin_fill(b)$pattern, "none")
  }
})

# ---- plot_node_rug: input validation ----------------------------------------
test_that("stops when presence is not a matrix", {
  expect_error(
    plot_node_rug(presence = "not a matrix", cell_h = 1, cell_w = 1, n_cols = 1),
    "must be a matrix"
  )
})

test_that("stops when support dims differ from presence", {
  p <- matrix(1, 3, 2, dimnames = list(c("4", "5", "6"), c("a", "b")))
  s <- matrix(1, 2, 2)
  expect_error(
    plot_node_rug(p, support = s, cell_h = 1, cell_w = 1, n_cols = 2),
    "same dimensions"
  )
})

# ---- plot_node_rug: draws without error on a real device --------------------
test_that("tier 1 draws without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  presence <- matrix(
    c(1, 0, 1, 1,
      1, 1, 0, 1),
    nrow = 4, ncol = 2,
    dimnames = list(as.character(6:9), c("iqtree", "astral"))
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  expect_no_error(
    plot_node_rug(
      presence = presence, cell_h = 0.2, cell_w = 0.2,
      n_cols = 2, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("tier 3 draws without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  presence <- matrix(
    c(1, 0, 1, 1,  1, 1, 0, 1),
    nrow = 4, ncol = 2,
    dimnames = list(as.character(6:9), c("iqtree", "astral"))
  )
  support <- matrix(
    c(99, NA, 85, 70,  0.99, 0.60, NA, 0.97),
    nrow = 4, ncol = 2,
    dimnames = list(as.character(6:9), c("iqtree", "astral"))
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  expect_no_error(
    plot_node_rug(
      presence = presence, support = support,
      support_type = c(iqtree = "sh_alrt", astral = "lpp"),
      cell_h = 0.2, cell_w = 0.2, n_cols = 2, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})
