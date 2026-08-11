
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phylorug <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/mdrifathahamed/phylorug/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdrifathahamed/phylorug/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/mdrifathahamed/phylorug/branch/master/graphs/badge.svg)](https://codecov.io/github/mdrifathahamed/phylorug)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status:
WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

<!-- badges: end -->

------------------------------------------------------------------------

**phylorug** is an R package for comparing and visualizing clade
recovery and support values across multiple phylogenetic trees on a
single reference topology using rug plots. Modern phylogenomic studies
routinely produce trees from different inference pipelines such as
concatenation versus coalescent methods, site-homogeneous versus
site-heterogeneous substitution models, different data filtering
strategies — each with support values in incompatible formats (ultrafast
bootstrap, SH-aLRT, posterior probability, local posterior probability).
**phylorug** maps all of these onto one tree, letting the researcher see
at a glance which nodes are robust across methods and which are not.

### Presence mode

<figure>
<img src="man/figures/README-presence.png" alt="Presence mode" />
<figcaption aria-hidden="true">Presence mode</figcaption>
</figure>

### Support mode

<figure>
<img src="man/figures/README-support.png" alt="Support mode" />
<figcaption aria-hidden="true">Support mode</figcaption>
</figure>

### Installation

**phylorug** is not yet on CRAN. To install the development version from
GitHub:

``` r
# install.packages("devtools")
devtools::install_github("mdrifathahamed/phylorug")
library(phylorug)
```

### Quick example

``` r
library(phylorug)

# Load bundled example trees
data(sample_trees)

# Build the node presence matrix
npm <- node_presence_matrix(sample_trees)

# Presence mode
plot_phylorug(npm, mode = "presence")

# Support mode
plot_phylorug(npm, mode = "support",
              support_type <- c(
                "70p_ASTRAL_partition_entropy" = "lpp",
                "70p_ASTRAL_uce"               = "lpp",
                "70p_ghost"                    = "ufboot",
                "70p_partition_entropy"        = "ufboot"
              ),)
```

### Two modes

**Presence mode** shows whether each analysis recovers a given clade
(black = recovered, white = not recovered). This is the direct
descendant of the “Navajo rug” plots introduced by Wheeler (1995) and
formalized by Giribet (2003).

**Support mode** shows how strongly each analysis supports a given
clade, with greyscale shading from low (light) to high (dark) support.
Support values from different formats are automatically normalized to a
common scale.

### Features

- Reads Newick and Nexus tree files via `read_trees()`
- Handles compound IQ-TREE labels (e.g. `SH-aLRT/UFBoot2`)
- Translates tip labels between naming conventions with
  `translate_tips()`
- Prunes to shared taxa across analyses with `prune_to_shared()`
- Validates taxon overlap with `check_taxa()`
- Treats a bare `multiPhylo` as one analysis (a pool of tied-optimal
  trees)
- Explicit errors over silent failures throughout

### Dependencies

**phylorug** depends on [ape](https://cran.r-project.org/package=ape)
and [phangorn](https://cran.r-project.org/package=phangorn). No other
external dependencies are required.

### Help

An overview of the package with links to all function documentation:

``` r
?phylorug
```

### Getting help

If you find a bug or have a feature request, please open an issue on
[GitHub](https://github.com/mdrifathahamed/phylorug/issues).

### Citation

If you use **phylorug** in a publication, please cite:

> Ahamed, M.R., Tarasov, S., and Arias, J.S. (2026). phylorug: Visualize
> Clade Recovery and Support Across Phylogenetic Trees. R package
> version 0.1.0. <https://github.com/mdrifathahamed/phylorug>

### Acknowledgements

This package was developed at the [Finnish Museum of Natural History
(LUOMUS)](https://www.luomus.fi/en), University of Helsinki, as part of
an MSc thesis in Ecology and Evolutionary Biology under the supervision
of Dr. Sergei Tarasov.

The rug-plot concept for phylogenetic sensitivity analysis traces to
Wheeler (1995), was formalized by Giribet (2003), and automated by
Sanders (2010, Cladescan) and Machado (2015, YBYRÁ).
