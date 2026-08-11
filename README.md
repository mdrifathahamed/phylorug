
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phylorug <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/mdrifathahamed/phylorug/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mdrifathahamed/phylorug/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Project Status:
WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/mdrifathahamed/phylorug/graph/badge.svg)](https://app.codecov.io/gh/mdrifathahamed/phylorug)

<!-- badges: end -->

------------------------------------------------------------------------

\##Overview

Phylogenetic trees now routinely include hundreds or thousands of taxa,
and as sequencing technologies become cheaper and datasets grow, trees
will only get larger. Additionally, modern phylogenomic studies
routinely produce multiple species trees for the same underlying
hypothesis. Researchers combine different data types (morphology, UCEs,
transcriptomes, whole genomes), apply different analytical strategies
(concatenation with varying partitioning schemes, coalescent methods,
site-heterogeneous models), use different alignment approaches, and
employ different inference software (IQ-TREE, ASTRAL, MrBayes, RAxML).
Each combination yields a tree with its own support values in its own
format, such as ultrafast bootstrap, SH-aLRT, posterior probability, and
local posterior probability, and the central question becomes: which
nodes hold up across methods, and how strongly?

The rug plot concept for comparing clade recovery across analytical
conditions was introduced by Wheeler (1995) and popularized as “Navajo
rugs” by Giribet (2003). Software to automate these plots, Cladescan
(Sanders 2010) and YBYRÁ (Machado 2015), both of which are no longer
maintained, were designed for parsimony parameter sensitivity rather
than modern multi-method phylogenomic comparison. No existing tool maps
heterogeneous support values from multiple inference pipelines onto a
single tree within R.

phylorug maps clade recovery and support values from any number of
independently inferred phylogenetic trees onto a single reference
topology as a rug plot. It operates in two modes:

### Two modes

### Presence mode

Black/white cells showing whether each tree recovers a given clade (the
direct descendant of Wheeler’s space plots, but phylorug is not limited
to parameters) ![Presence mode](man/figures/README-presence.png)

### Support mode

![Support mode](man/figures/README-support.png) Black, grey, yellow,
white, and red cells showing how strongly each analysis supports a given
clade, with automatic normalization across heterogeneous support
formats.

The package reads tree files directly, handles compound IQ-TREE labels,
translates tip names between naming conventions, prunes to shared taxa,
and produces publication-ready figures from R  no manual placement, no
external software, no disposable scripts.

### Features

- Reads Newick and Nexus tree files via `read_trees()`
- Handles compound IQ-TREE labels (e.g. `SH-aLRT/UFBoot2`)
- Translates tip labels between naming conventions with
  `translate_tips()`
- Prunes to shared taxa across analyses with `prune_to_shared()`
- Validates taxon overlap with `check_taxa()`
- Treats a bare `multiPhylo` as one analysis (a pool of tied-optimal
  trees)
- Explicit errors over silent failures throughout \### Installation

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
