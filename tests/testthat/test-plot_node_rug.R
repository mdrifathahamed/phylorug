# Tests for plot_node_rug() and its helpers
#
# resolve_tier: 1L = presence (support NULL), 2L = support (support given)
# resolve_cell: returns list(fill, pattern, border)
# bin_fill: returns list(fill, pattern)

# ---- resolve_tier -----------------------------------------------------------
test_that("tier 1 when support is NULL", {
  expect_equal(resolve_tier(NULL, NULL), 1L)
  expect_equal(resolve_tier(NULL, c(t1 = "ufboot")), 1L)
})

test_that("tier 2 when support is given (regardless of support_type)", {
  m <- matrix(1, 2, 2)
  expect_equal(resolve_tier(m, NULL), 2L)
  expect_equal(resolve_tier(m, c(t1 = "ufboot")), 2L)
})


# ---- resolve_cell: NA presence (not evaluable) ------------------------------
test_that("NA presence returns red fill in tier 1", {
  cell <- resolve_cell(NA_real_, NA_real_, NA_character_, NULL, 1L)
  expect_equal(cell$fill, "#D64545")
  expect_equal(cell$pattern, "none")
  expect_equal(cell$border, "grey40")
})

test_that("NA presence returns red fill in tier 2", {
  cell <- resolve_cell(NA_real_, NA_real_, "ufboot", NULL, 2L)
  expect_equal(cell$fill, "#D64545")
  expect_equal(cell$pattern, "none")
})


# ---- resolve_cell: absent (presence = 0) ------------------------------------
test_that("absent clade is white in tier 1", {
  cell <- resolve_cell(0, NA_real_, NA_character_, NULL, 1L)
  expect_equal(cell$fill, "white")
  expect_equal(cell$pattern, "none")
  expect_equal(cell$border, "grey40")
})

test_that("absent clade is white in tier 2", {
  cell <- resolve_cell(0, NA_real_, "ufboot", NULL, 2L)
  expect_equal(cell$fill, "white")
  expect_equal(cell$pattern, "none")
})


# ---- resolve_cell: present, tier 1 (presence mode) -------------------------
test_that("tier 1: fully present clade (p=1) is solid black", {
  cell <- resolve_cell(1, NA_real_, NA_character_, NULL, 1L)
  expect_equal(cell$fill, grDevices::rgb(0, 0, 0))
  expect_equal(cell$pattern, "none")
})

test_that("tier 1: pool proportion is grey by proportion", {
  cell <- resolve_cell(0.5, NA_real_, NA_character_, NULL, 1L)
  expect_equal(cell$fill, grDevices::rgb(0.5, 0.5, 0.5))
})

test_that("tier 1: p=0.25 gives lighter grey than p=0.75", {
  light <- resolve_cell(0.25, NA_real_, NA_character_, NULL, 1L)
  dark  <- resolve_cell(0.75, NA_real_, NA_character_, NULL, 1L)
  # Lower p = lighter (more white), higher p = darker (more black)
  expect_true(light$fill > dark$fill)  # rgb string comparison
})


# ---- resolve_cell: present, tier 2 (support mode) --------------------------
test_that("tier 2: present but no support value is red (not-computed)", {
  cell <- resolve_cell(1, NA_real_, "ufboot", NULL, 2L)
  expect_equal(cell$fill, "#D64545")
  expect_equal(cell$pattern, "none")
})

test_that("tier 2: present with support is shaded by bin", {
  # UFBoot 99 >= 98 → very high → bin 4
  cell <- resolve_cell(1, 99, "ufboot", NULL, 2L)
  expect_equal(cell$fill, bin_fill(4L)$fill)
  expect_equal(cell$pattern, bin_fill(4L)$pattern)
})

test_that("tier 2: low support uses yellow fill", {
  # UFBoot 40 < 50 → low → bin 1
  cell <- resolve_cell(1, 40, "ufboot", NULL, 2L)
  expect_equal(cell$fill, "#E8C547")
})

test_that("tier 2: pool proportion with support still uses binned fill", {
  # p=0.67 (partial pool recovery) with support value
  cell <- resolve_cell(0.67, 96, "ufboot", NULL, 2L)
  expect_equal(cell$fill, bin_fill(3L)$fill)  # 96 >= 95 → high → bin 3
})


# ---- bin_support: NA passthrough --------------------------------------------
test_that("bin_support returns NA for NA input", {
  expect_true(is.na(bin_support(NA_real_, "ufboot", NULL)))
})


# ---- bin_support: UFBoot thresholds (very_high=98, high=95, moderate=50) ----
test_that("UFBoot binning follows the default thresholds", {
  expect_equal(bin_support(99, "ufboot", NULL), 4L)   # >= 98 very high
  expect_equal(bin_support(98, "ufboot", NULL), 4L)   # boundary: exactly 98
  expect_equal(bin_support(97, "ufboot", NULL), 3L)   # < 98, >= 95 high
  expect_equal(bin_support(95, "ufboot", NULL), 3L)   # boundary: exactly 95
  expect_equal(bin_support(70, "ufboot", NULL), 2L)   # < 95, >= 50 moderate
  expect_equal(bin_support(50, "ufboot", NULL), 2L)   # boundary: exactly 50
  expect_equal(bin_support(49, "ufboot", NULL), 1L)   # < 50 low
  expect_equal(bin_support(0,  "ufboot", NULL), 1L)   # zero
})


# ---- bin_support: SH-aLRT thresholds (very_high=98, high=80, moderate=50) ---
test_that("SH-aLRT binning follows the default thresholds", {
  expect_equal(bin_support(99, "sh_alrt", NULL), 4L)   # >= 98
  expect_equal(bin_support(85, "sh_alrt", NULL), 3L)   # >= 80
  expect_equal(bin_support(60, "sh_alrt", NULL), 2L)   # >= 50
  expect_equal(bin_support(30, "sh_alrt", NULL), 1L)   # < 50
})


# ---- bin_support: ASTRAL LPP thresholds (very_high=0.99, high=0.95, mod=0.5)
test_that("ASTRAL LPP binning follows the default thresholds", {
  expect_equal(bin_support(0.99, "lpp", NULL), 4L)
  expect_equal(bin_support(0.96, "lpp", NULL), 3L)
  expect_equal(bin_support(0.70, "lpp", NULL), 2L)
  expect_equal(bin_support(0.30, "lpp", NULL), 1L)
})


# ---- bin_support: posterior (same scale as LPP) -----------------------------
test_that("posterior binning follows the default thresholds", {
  expect_equal(bin_support(0.99, "posterior", NULL), 4L)
  expect_equal(bin_support(0.96, "posterior", NULL), 3L)
  expect_equal(bin_support(0.70, "posterior", NULL), 2L)
  expect_equal(bin_support(0.30, "posterior", NULL), 1L)
})


# ---- bin_support: cross-type comparison -------------------------------------
test_that("same numeric value bins differently by support type", {
  expect_equal(bin_support(85, "sh_alrt", NULL), 3L)   # 85 >= 80 → high
  expect_equal(bin_support(85, "ufboot",  NULL), 2L)   # 85 < 95 → moderate
})


# ---- bin_support: custom thresholds -----------------------------------------
test_that("user thresholds override the defaults", {
  custom <- list(ufboot = c(very_high = 99, high = 90, moderate = 70))
  expect_equal(bin_support(99, "ufboot", custom), 4L)
  expect_equal(bin_support(95, "ufboot", custom), 3L)
  expect_equal(bin_support(80, "ufboot", custom), 2L)
  expect_equal(bin_support(60, "ufboot", custom), 1L)
})

test_that("custom thresholds for one type do not affect others", {
  custom <- list(ufboot = c(very_high = 99, high = 90, moderate = 70))
  # sh_alrt should still use defaults since custom doesn't include it
  expect_equal(bin_support(85, "sh_alrt", custom), 3L)
})


# ---- bin_support: unknown type fallback -------------------------------------
test_that("unknown support type falls back to a percentage scale", {
  expect_equal(bin_support(96, "mystery", NULL), 4L)   # >= 95
  expect_equal(bin_support(85, "mystery", NULL), 3L)   # >= 80
  expect_equal(bin_support(60, "mystery", NULL), 2L)   # >= 50
  expect_equal(bin_support(30, "mystery", NULL), 1L)   # < 50
})


# ---- bin_support: boundary values -------------------------------------------
test_that("bin_support handles exact boundary values correctly", {
  # Exactly at each threshold for ufboot
  expect_equal(bin_support(98, "ufboot", NULL), 4L)    # >= 98 → very high
  expect_equal(bin_support(95, "ufboot", NULL), 3L)    # >= 95, < 98 → high
  expect_equal(bin_support(50, "ufboot", NULL), 2L)    # >= 50, < 95 → moderate
})

test_that("bin_support handles zero and negative values", {
  expect_equal(bin_support(0, "ufboot", NULL), 1L)
})


# ---- bin_fill ---------------------------------------------------------------
test_that("bin_fill returns a list with fill and pattern", {
  result <- bin_fill(4L)
  expect_true(is.list(result))
  expect_named(result, c("fill", "pattern"))
})

test_that("bin_fill maps each bin to a distinct fill colour", {
  fills <- vapply(1:4, function(b) bin_fill(b)$fill, character(1))
  expect_length(unique(fills), 4L)
})

test_that("bin_fill returns black for NA (present but unquantified)", {
  result <- bin_fill(NA_integer_)
  expect_equal(result$fill, "black")
  expect_equal(result$pattern, "none")
})

test_that("bin 4 (very high) is black", {
  expect_equal(bin_fill(4L)$fill, "#000000")
})

test_that("bin 3 (high) is dark grey", {
  expect_equal(bin_fill(3L)$fill, "#5F5E5A")
})

test_that("bin 2 (moderate) is light grey", {
  expect_equal(bin_fill(2L)$fill, "#B4B2A9")
})

test_that("bin 1 (low) is yellow", {
  expect_equal(bin_fill(1L)$fill, "#E8C547")
})

test_that("all bins use no pattern", {
  for (b in 1:4) {
    expect_equal(bin_fill(b)$pattern, "none")
  }
})


# ---- default_thresholds -----------------------------------------------------
test_that("default_thresholds returns named vector for known types", {
  for (type in c("ufboot", "sh_alrt", "lpp", "posterior")) {
    th <- default_thresholds(type)
    expect_named(th, c("very_high", "high", "moderate"))
    expect_true(th[["very_high"]] > th[["high"]])
    expect_true(th[["high"]] > th[["moderate"]])
  }
})

test_that("default_thresholds falls back for unknown type", {
  th <- default_thresholds("mystery")
  expect_named(th, c("very_high", "high", "moderate"))
})

test_that("default_thresholds handles NULL support_type via %||%", {
  th <- default_thresholds(NULL)
  # Falls back to ufboot via %||%
  expect_equal(th, default_thresholds("ufboot"))
})


# ---- plot_node_rug: input validation ----------------------------------------
test_that("stops when npm is not a matrix", {
  expect_error(
    plot_node_rug(npm = "not a matrix", cell_h = 1, cell_w = 1, n_cols = 1),
    "must be a matrix"
  )
})

test_that("stops when npm is a data frame", {
  expect_error(
    plot_node_rug(npm = data.frame(a = 1), cell_h = 1, cell_w = 1, n_cols = 1),
    "must be a matrix"
  )
})

test_that("stops when support dims differ from npm", {
  p <- matrix(1, 3, 2, dimnames = list(c("4", "5", "6"), c("a", "b")))
  s <- matrix(1, 2, 2)
  expect_error(
    plot_node_rug(npm = p, support = s, cell_h = 1, cell_w = 1, n_cols = 2),
    "same dimensions"
  )
})

test_that("rug_position rejects invalid values", {
  p <- matrix(1, 2, 2, dimnames = list(c("4", "5"), c("a", "b")))
  expect_error(
    plot_node_rug(npm = p, cell_h = 1, cell_w = 1, n_cols = 2,
                  rug_position = "invalid"),
    "should be one of"
  )
})


# ---- plot_node_rug: draws without error on a real device --------------------
test_that("tier 1 (presence mode) draws without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
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
      npm = npm, cell_h = 0.2, cell_w = 0.2,
      n_cols = 2, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("tier 2 (support mode) draws without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
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
      npm = npm, support = support,
      support_type = c(iqtree = "sh_alrt", astral = "lpp"),
      cell_h = 0.2, cell_w = 0.2, n_cols = 2, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("inside and outside rug_position both draw without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
    c(1, 0, 1, 1),
    nrow = 4, ncol = 1,
    dimnames = list(as.character(6:9), "t1")
  )
  for (pos in c("inside", "outside")) {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(tmp)
    ape::plot.phylo(backbone)
    last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
    expect_no_error(
      plot_node_rug(
        npm = npm, cell_h = 0.2, cell_w = 0.2,
        n_cols = 1, rug_position = pos, last_pp = last_pp
      )
    )
    grDevices::dev.off()
    unlink(tmp)
  }
})

test_that("x_offset and y_offset shift the grid without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
    c(1, 0, 1, 1),
    nrow = 4, ncol = 1,
    dimnames = list(as.character(6:9), "t1")
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  expect_no_error(
    plot_node_rug(
      npm = npm, cell_h = 0.2, cell_w = 0.2,
      n_cols = 1, x_offset = 0.05, y_offset = 0.05, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("pool proportions (fractional presence) draw without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
    c(1, 0.67, 0.33, 1),
    nrow = 4, ncol = 1,
    dimnames = list(as.character(6:9), "pool")
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  expect_no_error(
    plot_node_rug(
      npm = npm, cell_h = 0.2, cell_w = 0.2,
      n_cols = 1, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("NA presence values draw without error", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
    c(1, NA, 0, 1),
    nrow = 4, ncol = 1,
    dimnames = list(as.character(6:9), "t1")
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  expect_no_error(
    plot_node_rug(
      npm = npm, cell_h = 0.2, cell_w = 0.2,
      n_cols = 1, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("multi-column grid (n_cols > 1) with many trees draws without error", { # nolint: line_length_linter.
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
    rep(c(1, 0, 1, 1), 4),
    nrow = 4, ncol = 4,
    dimnames = list(as.character(6:9), paste0("t", 1:4))
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  expect_no_error(
    plot_node_rug(
      npm = npm, cell_h = 0.15, cell_w = 0.15,
      n_cols = 2, last_pp = last_pp
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("last_pp = NULL fetches from ape environment", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
    c(1, 0, 1, 1),
    nrow = 4, ncol = 1,
    dimnames = list(as.character(6:9), "t1")
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  # Don't pass last_pp — let it fetch from ape environment
  expect_no_error(
    plot_node_rug(
      npm = npm, cell_h = 0.2, cell_w = 0.2,
      n_cols = 1, last_pp = NULL
    )
  )
  grDevices::dev.off()
  unlink(tmp)
})

test_that("returns invisible NULL", {
  skip_if_not_installed("ape")
  backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
  npm <- matrix(
    c(1, 0, 1, 1),
    nrow = 4, ncol = 1,
    dimnames = list(as.character(6:9), "t1")
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  ape::plot.phylo(backbone)
  last_pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  result <- plot_node_rug(
    npm = npm, cell_h = 0.2, cell_w = 0.2,
    n_cols = 1, last_pp = last_pp
  )
  expect_null(result)
  grDevices::dev.off()
  unlink(tmp)
})


# ---- draw_cell --------------------------------------------------------------
test_that("draw_cell draws without error on a device", {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  plot.new()
  expect_no_error(
    draw_cell(0, 0, 1, 1, list(fill = "black", pattern = "none", border = "grey40")) # nolint: line_length_linter.
  )
  grDevices::dev.off()
  unlink(tmp)
})


# ---- %||% operator ---------------------------------------------------------
test_that("%||% returns x when x is not NULL", {
  expect_equal(5 %||% 10, 5)
  expect_equal("a" %||% "b", "a")
})

test_that("%||% returns y when x is NULL", {
  expect_equal(NULL %||% 10, 10)
  expect_equal(NULL %||% "fallback", "fallback")
})
