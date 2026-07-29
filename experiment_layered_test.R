# experiment_layered_test.R
# Test script for the layered phylorug plotting system
# Run this from the phylorug project root after devtools::load_all()
#
# This script compares the old monolithic plot_phylorug() with the new
# layered system side by side.

library(devtools)
load_all()  # loads phylorug from source

# ============================================================================
# 1. TINY TREE — smoke test
# ============================================================================

backbone <- ape::read.tree(text = "((((A:1,B:1):0.5,C:1.5):0.3,(D:1,E:1):0.8):0.2,F:2);")
backbone$node.label <- c("", "95", "80", "100", "70")

# Two fake comparison trees with different topologies
tree_1 <- ape::read.tree(text = "((((A:1,B:1):0.5,C:1.5):0.3,(D:1,E:1):0.8):0.2,F:2);")
tree_1$node.label <- c("", "90", "75", "100", "65")

tree_2 <- ape::read.tree(text = "(((A:1,(B:1,C:1.5):0.5):0.3,(D:1,E:1):0.8):0.2,F:2);")
tree_2$node.label <- c("", "88", "60", "99", "55")

trees <- list(tree_1 = tree_1, tree_2 = tree_2)

npm <- node_presence_matrix(backbone, trees)

cat("=== npm structure ===\n")
str(npm, max.level = 1)
cat("\n=== presence matrix ===\n")
print(npm$presence)

# --- Old way (monolithic) ---
cat("\n--- Old monolithic plot_phylorug() ---\n")
tmp_old <- tempfile(fileext = ".pdf")
plot_phylorug(backbone, npm, file = tmp_old)
cat("Old plot saved to:", tmp_old, "\n")

# --- New way (layered) ---
cat("\n--- New layered system ---\n")
tmp_new <- tempfile(fileext = ".pdf")

p <- phylorug_plot(backbone, npm, file = tmp_new) +
  rug_layer() +
  dot_layer() +
  support_labels() +
  legend_layer()

print(p)  # this triggers rendering
cat("Layered plot saved to:", tmp_new, "\n")


# ============================================================================
# 2. LAYER COMBINATIONS — test composability
# ============================================================================

cat("\n--- Testing layer combinations ---\n")

# Just tree + dots (no rugs)
tmp_dots <- tempfile(fileext = ".pdf")
p_dots <- phylorug_plot(backbone, npm, file = tmp_dots) +
  dot_layer(col = "blue", cex = 1.0)
print(p_dots)
cat("Dots-only plot:", tmp_dots, "\n")

# Tree + rugs only (no dots, no labels)
tmp_rugs <- tempfile(fileext = ".pdf")
p_rugs <- phylorug_plot(backbone, npm, file = tmp_rugs) +
  rug_layer(rug_position = "outside")
print(p_rugs)
cat("Rugs-only plot:", tmp_rugs, "\n")

# Tree + labels only
tmp_labs <- tempfile(fileext = ".pdf")
p_labs <- phylorug_plot(backbone, npm, file = tmp_labs) +
  support_labels(col = "blue", cex = 0.8)
print(p_labs)
cat("Labels-only plot:", tmp_labs, "\n")

# Full stack with custom settings
tmp_full <- tempfile(fileext = ".pdf")
p_full <- phylorug_plot(backbone, npm, file = tmp_full) +
  rug_layer(cell_scale = 0.6, rug_position = "outside") +
  dot_layer(col = "darkgreen") +
  support_labels(col = "purple") +
  legend_layer()
print(p_full)
cat("Full custom plot:", tmp_full, "\n")


# ============================================================================
# 3. ON-SCREEN RENDERING — no file, just plot
# ============================================================================

cat("\n--- On-screen rendering (if interactive) ---\n")
if (interactive()) {
  p_screen <- phylorug_plot(backbone, npm) +
    rug_layer() +
    dot_layer() +
    support_labels()

  print(p_screen)
  cat("Check your plot window!\n")
} else {
  cat("Not interactive, skipping screen rendering.\n")
}


# ============================================================================
# 4. THE KEY TEST — modifying after creation
# ============================================================================

cat("\n--- Modify-after-creation test ---\n")
tmp_mod <- tempfile(fileext = ".pdf")

# Create base
p_base <- phylorug_plot(backbone, npm, file = tmp_mod)

# Add layers one at a time
p_base <- p_base + rug_layer()
p_base <- p_base + dot_layer()

# Can inspect layers before rendering
cat("Number of layers:", length(p_base$layers), "\n")
cat("Layer types:",
    paste(vapply(p_base$layers, function(l) l$type, character(1)),
          collapse = ", "), "\n")

print(p_base)
cat("Modified plot saved to:", tmp_mod, "\n")


# ============================================================================
# 5. COMPARISON CHECKLIST
# ============================================================================

cat("\n")
cat("=============================================\n")
cat("  VISUAL COMPARISON CHECKLIST\n")
cat("=============================================\n")
cat("Open the PDFs and check:\n")
cat("  [ ] Rug cells render at correct nodes\n")
cat("  [ ] Dots appear at unanimous nodes\n")
cat("  [ ] Support labels appear at variable nodes\n")
cat("  [ ] Legend renders in top-left\n")
cat("  [ ] Rugs-only: no dots, no labels\n")
cat("  [ ] Dots-only: no rugs, no labels\n")
cat("  [ ] Labels-only: no rugs, no dots\n")
cat("  [ ] Old vs new: identical visual output\n")
cat("=============================================\n")
