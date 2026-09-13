# phylorug: Visualize Clade Recovery and Support Across Phylogenetic Trees

Generates rug plot visualizations directly on tree nodes to illustrate
clade recovery and support from multiple phylogenetic and phylogenomic
analyses. See
[`vignette("phylorug")`](https://mdrifathahamed.github.io/phylorug/articles/phylorug.md)
for a detailed guide.

## Core functions

- [`read_trees`](https://mdrifathahamed.github.io/phylorug/reference/read_trees.md):

  Read tree files from disk

- [`translate_tips`](https://mdrifathahamed.github.io/phylorug/reference/translate_tips.md):

  Rename specimen codes to species names

- [`check_taxa`](https://mdrifathahamed.github.io/phylorug/reference/check_taxa.md):

  Diagnose taxon consistency across trees

- [`prune_to_shared`](https://mdrifathahamed.github.io/phylorug/reference/prune_to_shared.md):

  Prune trees to shared taxa

- [`node_presence_matrix`](https://mdrifathahamed.github.io/phylorug/reference/node_presence_matrix.md):

  Build the presence/support matrix

- [`plot_phylorug`](https://mdrifathahamed.github.io/phylorug/reference/plot_phylorug.md):

  Draw the rug plot on a reference tree

- [`add_tree`](https://mdrifathahamed.github.io/phylorug/reference/add_tree.md):

  Add a tree to an existing matrix

## References

Wheeler, W.C. (1995). Sequence alignment, parameter sensitivity, and the
phylogenetic analysis of molecular data. *Systematic Biology*, 44(3),
321-331.
[doi:10.1093/sysbio/44.3.321](https://doi.org/10.1093/sysbio/44.3.321)

Giribet, G. (2003). Stability in phylogenetic formulations and its
relationship to nodal support. *Systematic Biology*, 52(4), 554-564.
[doi:10.1080/10635150390223730](https://doi.org/10.1080/10635150390223730)

Sanders, J.G. (2010). Program note: Cladescan, a program for automated
phylogenetic sensitivity analysis. *Cladistics*, 26(1), 114-116.
[doi:10.1111/j.1096-0031.2009.00280.x](https://doi.org/10.1111/j.1096-0031.2009.00280.x)

Machado, D.J. (2015). YBYRA facilitates comparison of large phylogenetic
trees. *BMC Bioinformatics*, 16, 204.
[doi:10.1186/s12859-015-0642-9](https://doi.org/10.1186/s12859-015-0642-9)

Wolfe, J.M. et al. (2019). A phylogenomic framework, evolutionary
timeline and genomic resources for comparative studies of decapod
crustaceans. *Proceedings of the Royal Society B*, 286(1901), 20190079.
[doi:10.1098/rspb.2019.0079](https://doi.org/10.1098/rspb.2019.0079)

Fu, Y. et al. (2025). Phylogenomic insights into the higher level
relationships within Culicomorpha (Diptera). *Insect Systematics and
Diversity*, 9(6), ixaf056.
[doi:10.1093/isd/ixaf056](https://doi.org/10.1093/isd/ixaf056)

Montanaro, G., Lopes, F. et al. (2026). Phylogenomics resolves a
200-year-old puzzle: a revised tribal classification of Afro-Eurasian
dung beetles (Coleoptera: Scarabaeinae). *bioRxiv*.
[doi:10.64898/2026.07.22.740134](https://doi.org/10.64898/2026.07.22.740134)

## See also

Useful links:

- <https://github.com/mdrifathahamed/phylorug>

- <https://mdrifathahamed.github.io/phylorug/>

- Report bugs at <https://github.com/mdrifathahamed/phylorug/issues>

## Author

**Maintainer**: Md Rifath Ahamed <rifath.ahamed@helsinki.fi>
([ORCID](https://orcid.org/0009-0005-6158-8658)) \[copyright holder\]

Authors:

- Md Rifath Ahamed <rifath.ahamed@helsinki.fi>
  ([ORCID](https://orcid.org/0009-0005-6158-8658)) \[copyright holder\]

- Sergei Tarasov <sergei.tarasov@helsinki.fi>
  ([ORCID](https://orcid.org/0000-0002-1737-9403))

Other contributors:

- J. Salvador Arias <js.arias@conicet.gov.ar>
  ([ORCID](https://orcid.org/0000-0002-3717-435X)) \[contributor\]
