# gg_test.R -- disposable smoke test for gg_phylorug()
library(devtools)
load_all()
library(ggplot2)

# Toy example, matches the roxygen @examples in plot_phylorug
backbone <- ape::read.tree(text = "(((A,B),C),(D,E));")
backbone$node.label <- c("100", "88", "72", "95")

npm <- list(
  presence = matrix(
    c(1, 0, 1, 1, 1, 1, 0, 1),
    nrow = 4, ncol = 2,
    dimnames = list(as.character(6:9), c("tree_1", "tree_2"))
  )
)

# --- Presence mode -----------------------------------------------------------
p1 <- gg_phylorug(backbone, npm, mode = "presence")
ggsave("gg_test_presence.png", p1, width = 6, height = 4, dpi = 150)

# --- Support mode --------------------------------------------------------------
npm$support_1 <- matrix(
  c(99, NA, 85, 70, 0.99, 0.6, NA, 0.97),
  nrow = 4, ncol = 2,
  dimnames = list(as.character(6:9), c("tree_1", "tree_2"))
)

p2 <- gg_phylorug(
  backbone, npm, mode = "support",
  support_type = c(tree_1 = "sh_alrt", tree_2 = "lpp")
)
ggsave("gg_test_support.png", p2, width = 6, height = 4, dpi = 150)

cat("Done. Check gg_test_presence.png and gg_test_support.png\n")

# --- Real data (uncomment once toy example looks right) ---------------------
# sample_trees <- read_trees("data/sample_figure")
# backbone <- sample_trees[["70p_uce"]]
# others   <- sample_trees[names(sample_trees) != "70p_uce"]
# rugmt    <- node_presence_matrix(backbone, others, support_col = c(1, 2))
# p3 <- gg_phylorug(backbone, rugmt, mode = "presence")
# ggsave("gg_test_real_presence.png", p3, width = 10, height = 10, dpi = 200)
