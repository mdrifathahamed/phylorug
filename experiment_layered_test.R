# experiment_layered_real_data.R
# Test the layered plotting system on real datasets

library(devtools)
load_all()

# =============================================================================
# TEST 1: Culicomorpha — presence mode
# =============================================================================
cat("=== TEST 1: Culicomorpha ===\n")

sample_trees <- read_trees("data/Culicomorpha")
sample_trees <- lapply(sample_trees, function(tr) {
  ape::root(tr, outgroup = "Pseudosmittia_sp", resolve.root = TRUE)
})

backbone <- sample_trees[["Matrix1-smart_partitioning"]]
others   <- sample_trees[names(sample_trees) != "Matrix1-smart_partitioning"]
others   <- lapply(others, function(tr) ape::keep.tip(tr, backbone$tip.label))

check_taxa(backbone, others)
rugmt <- node_presence_matrix(backbone, others, support_col = c(1, 2))

support_type <- c(
  "Matrix1-kpi_ASTRAL"            = "lpp",
  "Matrix1-smart_ASTRAL"          = "lpp",
  "Matrix1-kpi_partitioning"      = "sh_alrt",
  "Matrix1-kpi_PMSF(H1.guide)"   = "sh_alrt",
  "Matrix1-kpi_PMSF(H2.guide)"   = "sh_alrt",
  "Matrix1_LG_C20_F_R"           = "sh_alrt",
  "Matrix1-smart_PMSF(H1.guide)" = "sh_alrt",
  "Matrix1-smart_PMSF(H2.guide)" = "sh_alrt",
  "Matrix2_partitioning"          = "sh_alrt"
)

# --- OLD way ---
plot_phylorug(backbone, rugmt,
              file = "culico_OLD_presence.pdf",
              cell_scale = 0.45, dot_cex = 0.65, font = 3)

# --- NEW layered way — presence ---
p1 <- phylorug_plot(backbone, rugmt,
                    file = "culico_NEW_presence.pdf",
                    font = 3) +
  rug_layer(cell_scale = 0.45) +
  dot_layer(col = "black", cex = 0.65) +
  support_labels(col = "red", cex = 0.3, font = 3) +
  legend_layer(position_pos = "topleft", text_cex = 0.35)
print(p1)

# --- OLD support ---
plot_phylorug(backbone, rugmt,
              file = "culico_OLD_support.pdf",
              width = 10, height = 10,
              mode = "support", support_type = support_type,
              show_support = TRUE, support_label_cex = 0.25,
              support_label_col = "Red", cell_scale = 0.25,
              font = 3, rug_position = "inside", edge.width = 1.5)

# --- NEW support ---
p2 <- phylorug_plot(backbone, rugmt,
                    file = "culico_NEW_support.pdf",
                    width = 10, height = 10,
                    mode = "support", support_type = support_type,
                    font = 3, edge.width = 1.5) +
  rug_layer(cell_scale = 0.25, rug_position = "inside") +
  dot_layer(col = "black", cex = 0.65) +
  support_labels(col = "Red", cex = 0.25, font = 3) +
  legend_layer(text_cex = 0.3)
print(p2)

cat("Culicomorpha done.\n\n")


# =============================================================================
# TEST 2: Sample figure (70p, ~54 tips)
# =============================================================================
cat("=== TEST 2: Sample figure ===\n")

trees_70p <- read_trees("data/sample_figure")
trees_70p <- lapply(trees_70p, function(tr) {
  ape::drop.tip(tr, c("Lucanus capreolus COL3897", "Xylonichus sp COL066",
                      "Anoplostethus sp COL892", "Trichaulax sp COL1050"))
})

backbone <- trees_70p[["70p_uce"]]
others   <- trees_70p[names(trees_70p) != "70p_uce"]
others   <- lapply(others, function(tr) ape::keep.tip(tr, backbone$tip.label))

check_taxa(backbone, others)
rugmt <- node_presence_matrix(backbone, others, support_col = c(1, 2))

support_type <- c(
  "70p_ASTRAL_partition_entropy" = "lpp",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ghost"                    = "ufboot",
  "70p_partition_entropy"        = "ufboot"
)

# --- OLD ---
plot_phylorug(backbone, rugmt,
              file = "sample_OLD_presence.pdf",
              cell_scale = 0.45, dot_cex = 0.65, font = 3)

# --- NEW ---
p3 <- phylorug_plot(backbone, rugmt,
                    file = "sample_NEW_presence.pdf",
                    font = 3) +
  rug_layer(cell_scale = 0.45) +
  dot_layer(cex = 0.65) +
  support_labels(font = 3) +
  legend_layer()
print(p3)

# --- OLD support ---
plot_phylorug(backbone, rugmt,
              file = "sample_OLD_support.pdf",
              mode = "support", support_type = support_type,
              show_support = TRUE, support_label_cex = 0.39,
              support_label_col = "black", cell_scale = 0.75,
              dot_cex = 0.65, font = 3)

# --- NEW support ---
p4 <- phylorug_plot(backbone, rugmt,
                    file = "sample_NEW_support.pdf",
                    mode = "support", support_type = support_type,
                    font = 3) +
  rug_layer(cell_scale = 0.75) +
  dot_layer(cex = 0.65) +
  support_labels(col = "black", cex = 0.39, font = 3) +
  legend_layer(position_pos = "topleft", threshold_pos = "topright")
print(p4)

cat("Sample figure done.\n\n")


# =============================================================================
# TEST 3: Beetles (70p, ~315 tips)
# =============================================================================
cat("=== TEST 3: Beetles ===\n")

trees_b <- read_trees("data/bettles/70p")
trees_b <- lapply(trees_b, function(tr) {
  tr <- ape::root(tr, outgroup = c("NicorbUCE", "NicvesUCE"), resolve.root = TRUE)
  ape::drop.tip(tr, c("NicorbUCE", "NicvesUCE"))
})

biogeo  <- readxl::read_excel("data/biogeo.xlsx", sheet = "BioGeo")
trees_b <- translate_tips(trees_b, biogeo, from = "from", to = "to")

backbone <- trees_b[["70p_uce"]]
others   <- trees_b[names(trees_b) != "70p_uce"]
others   <- lapply(others, function(tr) ape::keep.tip(tr, backbone$tip.label))

check_taxa(backbone, others)
rugmt <- node_presence_matrix(backbone, others, support_col = c(1, 2))

support_type <- c(
  "70p_ASTRAL_partition_entropy" = "lpp",
  "70p_ASTRAL_uce"               = "lpp",
  "70p_ghost"                    = "ufboot",
  "70p_partition_entropy"        = "ufboot"
)

# --- OLD ---
plot_phylorug(backbone, rugmt,
              file = "beetles_OLD_presence.pdf",
              cell_scale = 0.45, dot_cex = 0.65, font = 3)

# --- NEW ---
p5 <- phylorug_plot(backbone, rugmt,
                    file = "beetles_NEW_presence.pdf",
                    font = 3) +
  rug_layer(cell_scale = 0.45) +
  dot_layer(cex = 0.65) +
  support_labels(font = 3) +
  legend_layer()
print(p5)

# --- OLD support ---
plot_phylorug(backbone, rugmt,
              file = "beetles_OLD_support.pdf",
              mode = "support", support_type = support_type,
              show_support = TRUE, support_label_cex = 0.39,
              support_label_col = "black", cell_scale = 0.45,
              dot_cex = 0.65, font = 3)

# --- NEW support ---
p6 <- phylorug_plot(backbone, rugmt,
                    file = "beetles_NEW_support.pdf",
                    mode = "support", support_type = support_type,
                    font = 3) +
  rug_layer(cell_scale = 0.45) +
  dot_layer(cex = 0.65) +
  support_labels(col = "black", cex = 0.39, font = 3) +
  legend_layer(text_cex = 0.35)
print(p6)

cat("Beetles done.\n\n")

cat("=============================================\n")
cat("  ALL DONE — compare OLD vs NEW PDFs\n")
cat("=============================================\n")
