#' Import phylogenetic trees from a directory
#'
#' Scans a directory for phylogenetic tree files, parses each one, and returns
#' them as a named list. Format is detected from file content, so a directory
#' may mix Newick and NEXUS files.
#'
#' @details
#' Each file is one analysis. A file holding a single tree is returned as a
#' `"phylo"` ; a file holding several equally optimal trees from one search
#' (as from POY, TNT, or PAUP*) is returned as a `"multiPhylo"` and is
#' scored as a pool. How a pool's clade recovery is summarised (as a
#' continuous proportion, or binarised at a threshold) is decided later, by
#' [node_presence_matrix()], not here.
#'
#' Do not supply posterior samples, bootstrap replicates, or sets of gene trees.
#' These are distributions rather than analyses and must be summarised before
#' use. A file holding more than 100 trees is an error.
#'
#' Support values written as Newick node labels are read normally. BEAST-style
#' bracket annotations (`[&posterior=0.98]`) are discarded by \pkg{ape};
#' the topology is unaffected, and a message is emitted.
#'
#' Support values are imported exactly as written in the tree file;
#' [read_trees()] does not recompute or verify them. When a file holds
#' several tied-optimal trees, be aware that upstream programs differ in how
#' they summarise support across such trees: TNT and POY4 default to the more
#' conservative strict-consensus approach, whereas PAUP* and PHYLIP default to
#' the frequency-within-replicates approach, which Simmons and Freudenstein
#' (2011) showed can inflate apparent support for unsupported clades. This
#' choice is made by the upstream software before the file reaches
#' [read_trees()] and cannot be corrected here.
#'
#' Files matching `ext` that contain no tree, such as a NEXUS character
#' matrix, are skipped with a message. A tree file that fails to parse is an
#' error.
#'
#' @param dir Path to the directory containing the tree files. Defaults to the
#'   working directory.
#'
#' @param ext Character vector of file extensions to search for, without the
#'   leading dot. Matching is case-insensitive. You may supply your own, for
#'   example `ext = "treefile"` to read only IQ-TREE output. If you do not,
#'   the default filters the directory to the common tree file extensions, so
#'   that alignments, log files and configuration files sitting beside the trees
#'   are not read. Files with no extension, such as classic RAxML output, cannot
#'   be matched.
#'
#' @param format The parsing strategy to use. Defaults to `"auto"`, which
#'   detects the format of each file individually, so a directory may mix
#'   formats. Use `"newick"` or `"nexus"` to force one parser for all
#'   files.
#'
#' @param verbose Logical. If `TRUE` (default), reports trees read, files
#'   skipped, and pooled analyses. Errors are always reported.
#'
#' @return A named list, one element per tree file, named after the filenames
#'   with extensions removed. Single-tree files yield `"phylo"` objects;
#'   multi-tree files yield `"multiPhylo"`. Carries a `"pool_sizes"`
#'   attribute giving the number of trees per analysis.
#'
#' @references
#' Simmons, M.P. & Freudenstein, J.V. (2011). Spurious 99% bootstrap and
#' jackknife support for unsupported clades. \emph{Molecular Phylogenetics and
#' Evolution}, 61(1), 177-191. \doi{10.1016/j.ympev.2011.06.003}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # One tree per analysis. Each element is a `phylo`.
#' trees <- read_trees("path/to/your/trees")
#' ape::plot.phylo(trees[[1]])
#'
#' # Trees per analysis
#' attr(trees, "pool_sizes")
#'
#' # Restrict to one extension
#' trees <- read_trees("path/to/your/trees", ext = "treefile")
#'
#' # Force Newick parsing
#' trees <- read_trees("path/to/your/trees", format = "newick")
#' }
read_trees <- function(dir     = ".",
                       ext     = c("tre", "tree", "treefile", "nwk",
                                   "newick", "nex", "nexus", "contree"),
                       format  = c("auto", "newick", "nexus"),
                       verbose = TRUE) {
  format <- match.arg(format)

  if (!dir.exists(dir)) {
    stop("`dir` does not exist: ", dir, call. = FALSE)
  }
  if (!is.character(ext) || length(ext) == 0L) {
    stop(
      "`ext` must be a non-empty character vector of file extensions, ",
      "without the leading dot, e.g. ext = c(\"tre\", \"nex\").",
      call. = FALSE
    )
  }

  ext         <- ext[order(nchar(ext), decreasing = TRUE)]
  ext_pattern <- paste0("\\.(", paste(ext, collapse = "|"), ")$")

  files <- list.files(
    dir,
    pattern     = ext_pattern,
    full.names  = TRUE,
    ignore.case = TRUE
  )
  if (length(files) == 0L) {
    stop(
      "No files matching extensions (", paste(ext, collapse = ", "),
      ") found in: ", dir,
      call. = FALSE
    )
  }

  parsed <- lapply(files, function(f) {
    read_one_analysis(f, format = format)
  })

  is_skipped <- vapply(parsed, is.null, logical(1))
  if (any(is_skipped)) {
    if (verbose) {
      message(
        "Skipped ", sum(is_skipped),
        " file(s) containing no phylogenetic tree: ",
        paste(basename(files[is_skipped]), collapse = ", ")
      )
    }
    parsed <- parsed[!is_skipped]
    files  <- files[!is_skipped]
  }

  if (length(parsed) == 0L) {
    stop(
      "No phylogenetic trees found in: ", dir,
      ". Files matched `ext` but none contained a tree.",
      call. = FALSE
    )
  }

  nm <- tools::file_path_sans_ext(basename(files))
  if (anyDuplicated(nm)) {
    dup <- unique(nm[duplicated(nm)])
    stop(
      "Tree names are not unique after removing file extensions: ",
      paste(dup, collapse = ", "),
      ". Rename the files, or narrow `ext` so that only one is read.",
      call. = FALSE
    )
  }
  names(parsed) <- nm

  pool_sizes        <- vapply(parsed, pool_size, integer(1))
  names(pool_sizes) <- nm

  if (verbose) {
    is_pool <- pool_sizes > 1L
    message(
      "Read ", length(parsed), " analyses (", sum(pool_sizes),
      " trees) from: ", dir
    )
    if (any(is_pool)) {
      message(
        "Scored as pools: ",
        paste0(nm[is_pool], " (", pool_sizes[is_pool], " trees)",
               collapse = ", ")
      )
    }
  }

  attr(parsed, "pool_sizes") <- pool_sizes
  parsed
}


#' Coerce one analysis to a pool of trees
#'
#' Internal. Every function that iterates over the trees of an analysis calls
#' this first, so that single-tree and multi-tree analyses share one code path.
#' A \code{phylo} becomes a \code{multiPhylo} of length one; a \code{multiPhylo}
#' passes through unchanged.
#'
#' @noRd
as_pool <- function(x) {
  if (inherits(x, "multiPhylo")) {
    return(x)
  }
  if (inherits(x, "phylo")) {
    return(structure(list(x), class = "multiPhylo"))
  }
  stop(
    "Expected a `phylo` or `multiPhylo` object, got class ",
    paste(class(x), collapse = "/"), ".",
    call. = FALSE
  )
}


#' Number of trees in one analysis
#'
#' @noRd
pool_size <- function(x) {
  if (inherits(x, "phylo")) 1L else length(x)
}


#' Read a single file into one analysis
#'
#' Returns a \code{phylo} (single tree), a \code{multiPhylo} (pool), or
#' \code{NULL} when the file contains no tree and should be skipped.
#'
#' @noRd
read_one_analysis <- function(f, format) {

  # Maximum trees in one file. A file holding more is a posterior sample, a
  # bootstrap set, or trees pooled across several analytical conditions --
  # none of which is a single analysis.
  pool_max <- 100L

  fmt <- if (format == "auto") detect_format(f) else format
  # detect_format() returns "none" for a file that matched `ext` but holds no
  # tree, e.g. a PAUP*-format NEXUS character matrix.
  if (identical(fmt, "none")) {
    return(NULL)
  }

  reader <- switch(fmt,
                   nexus  = ape::read.nexus,
                   newick = ape::read.tree,
                   stop("Unrecognised format \"", fmt, "\" for file: ", basename(f),
                        call. = FALSE)
  )

  tr <- tryCatch(
    suppressWarnings(reader(f)),
    error = function(e) {
      stop(
        "Could not parse tree file ", basename(f), ": ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (is.null(tr)) {
    stop("No tree could be read from: ", basename(f), call. = FALSE)
  }
  if (!inherits(tr, "phylo") && !inherits(tr, "multiPhylo")) {
    stop(
      basename(f), " did not yield a phylogenetic tree (got class ",
      paste(class(tr), collapse = "/"), ").",
      call. = FALSE
    )
  }

  # ape returns a multiPhylo even for a NEXUS file holding exactly one tree.
  # Unwrap it, so that a single tree is a phylo whatever the file format.
  if (inherits(tr, "multiPhylo") && length(tr) == 1L) {
    tr <- tr[[1L]]
  }

  n <- pool_size(tr)
  if (n == 0L) {
    stop("No tree could be read from: ", basename(f), call. = FALSE)
  }

  # A file holding many trees is one analysis with a tied optimum. A file
  # holding very many is a posterior sample, a bootstrap set, or trees pooled
  # across several analytical conditions, none of which is a single analysis.
  if (n > pool_max) {
    stop(
      basename(f), " contains ", n, " trees, exceeding `pool_max` (", pool_max,
      "). This is usually a posterior sample, a bootstrap set, or trees pooled ",
      "across several analytical conditions. phylorug treats each file as one ",
      "analysis. If these are genuinely the equally optimal trees from one ",
      "search, raise `pool_max`.",
      call. = FALSE
    )
  }

  # ape discards BEAST-style bracket annotations, so support is silently lost
  # unless we say so. Only report when there is genuinely nothing to read.
  no_labels <- all(vapply(
    as_pool(tr),
    function(x) is.null(x$node.label),
    logical(1)
  ))
  if (no_labels && has_beast_annotations(f)) {
    message(
      basename(f), " contains BEAST-style node annotations. ape reads ",
      "topology only, so posterior probabilities are not imported. The tree ",
      "is fully usable for presence/absence rugs; see ?read_trees to import ",
      "support values."
    )
  }

  tr
}


#' Detect the format of a tree file
#'
#' Returns \code{"nexus"}, \code{"newick"}, or \code{"none"} when the file
#' contains no tree and should be skipped.
#'
#' @noRd
detect_format <- function(path) {
  txt <- readLines(path, warn = FALSE)

  if (length(txt) == 0L || all(!nzchar(trimws(txt)))) {
    stop("File appears to be empty: ", basename(path), call. = FALSE)
  }

  # #NEXUS need not be on the first physical line: files routinely open with a
  # blank line or a bracketed comment block. Check the head, not only line 1.
  head_lines <- utils::head(txt, 50L)
  if (any(grepl("^\\s*#NEXUS", head_lines, ignore.case = TRUE))) {
    # NEXUS is a container format. A character matrix is a valid NEXUS file
    # with no TREES block, and is not a tree file.
    if (!any(grepl("BEGIN\\s+TREES", txt, ignore.case = TRUE))) {
      return("none")
    }
    return("nexus")
  }

  # Otherwise assume Newick, but only if something tree-like is present. A
  # Newick tree must contain a parenthesis and terminate with a semicolon.
  has_paren <- any(grepl("(", txt, fixed = TRUE))
  has_semi  <- any(grepl(";", txt, fixed = TRUE))
  if (!has_paren || !has_semi) {
    return("none")
  }

  "newick"
}


#' Test whether a file carries BEAST-style node annotations
#'
#' BEAST and TreeAnnotator write node metadata as bracketed comments, e.g.
#' \code{[&posterior=0.98,rate=1.01]}. The Newick grammar treats brackets as
#' comments, so ape strips them and \code{node.label} comes back \code{NULL}.
#' The annotations must therefore be detected in the raw file text, before
#' parsing.
#'
#' @noRd
has_beast_annotations <- function(path) {
  txt <- readLines(path, warn = FALSE)
  any(grepl("[&", txt, fixed = TRUE))
}
