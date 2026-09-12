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
#' scored as a pool. How many pool trees must recover a clade before it
#' counts as present is controlled by the `pool_threshold` argument in
#' [node_presence_matrix()].
#'
#' Do not supply posterior samples, bootstrap replicates, or sets of gene trees.
#' These are distributions rather than analyses and must be summarised in the
#' upstream analysis before use (e.g. into a consensus tree). Files holding
#' more than 100 trees are rejected with an error.
#'
#' Support values stored as numeric node labels in the tree file (e.g.
#' `((A,B)95,(C,D)100);`) are read normally. BEAST-style bracket
#' annotations (`[&posterior=0.98]`) are discarded by \pkg{ape}; the
#' topology is unaffected, and a message is emitted. Support values are imported
#' exactly as written in the tree file; [read_trees()] does not recompute or
#' verify them.
#'
#' Files matching `ext` that contain no tree, such as a NEXUS character
#' matrix, are skipped with a message. A tree file that fails to parse is an
#' error.
#'
#' @param dir Path to the directory containing the tree files. Defaults to the
#'   working directory.
#'
#' @param ext Character vector of file extensions to look for, without the
#'   leading dot (e.g. `"tre"`, not `".tre"`). Case-insensitive. The default
#'   covers common tree extensions (`tre`, `tree`, `treefile`, `nwk`, `nex`,
#'   `nexus`, `contree`). Pass your own to narrow or widen the search.
#'
#' @param format The parsing strategy to use. Defaults to `"auto"`, which
#'   detects the format of each file individually, so a directory may mix
#'   formats. Use `"newick"` or `"nexus"` to force one parser for all
#'   files.
#'
#' @param verbose Logical. If `TRUE` (default), reports trees read, files
#'   skipped, and pooled analyses. Errors are always reported.
#'
#' @return A named list, one element per tree file, named after the file names
#'   with extensions removed. Single-tree files yield `"phylo"` objects;
#'   multi-tree files yield `"multiPhylo"`. Carries a `"pool_sizes"`
#'   attribute giving the number of trees per analysis.
#'
#' @export
#'
#' @examples
#' # phylorug ships example  raw tree files in its `extdata` directory.
#' # Point read_trees() at one of those folders:
#' dir <- system.file("extdata", "beetles_70p", package = "phylorug")
#'
#' # With only the directory, read_trees() detects each file's format from its
#' # contents (NEXUS or Newick, the two formats this version supports) and
#' # matches the default extension filter (tre, tree, treefile, nwk, ...).
#'
#' trees <- read_trees(dir)
#'
#' # The result is a named list, one element per analysis, ready to pass to
#' # [check_taxa()] or [node_presence_matrix()]:
#'
#' names(trees)
#'
#' # If a folder holds many files and you want only some, narrow `ext`
#' # to one extension, or to a set of them:
#' trees <- read_trees(dir, ext = "tre")
#' trees <- read_trees(dir, ext = c("tre", "tree", "treefile"))
#'
#' # For full control you can also force a single parser. Note this applies
#' # one format to every matched file, so a file of a different format would
#' # not be read.
#' trees <- read_trees(dir, format = "newick")
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
  has_dot <- grepl("^\\.", ext)
  if (any(has_dot)) {
    stop(
      "`ext` must not include a leading dot. ",
      "Use ext = \"tre\", not ext = \".tre\".",
      call. = FALSE
    )
  }
  has_blank <- !nzchar(trimws(ext))
  if (any(has_blank)) {
    stop(
      "`ext` must not contain empty or whitespace-only strings.",
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
      ". Rename the files, so that only one is read.",
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

#' Read a single file into one analysis
#'
#' Returns a `"phylo"` (single tree), a `"multiPhylo"` (pool), or `NULL` when
#' the file contains no tree and should be skipped.
#' @noRd
read_one_analysis <- function(f, format) {
  # A file holding more than 100 tres is a posterior sample, a bootstrap set, or
  # trees pooled across several analytical conditions none of which is a single
  # analysis.
  pool_max <- 100L
  fmt <- if (format == "auto") detect_format(f) else format
  if (identical(fmt, "none")) {
    return(NULL)
  }
  reader <- switch(
    fmt,
    nexus  = ape::read.nexus,
    newick = ape::read.tree,
    stop(
      "Unrecognised format \"", fmt,
      "\" for file: ", basename(f),
      call. = FALSE
    )
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
    stop("No tree could be read from: ", basename(f), call. = FALSE) # nocov
  }
  if (!inherits(tr, "phylo") && !inherits(tr, "multiPhylo")) {
    stop(                                                          # nocov start
      basename(f), " did not yield a phylogenetic tree (got class ",
      paste(class(tr), collapse = "/"), ").",
      call. = FALSE
    )                                                               # nocov end
  }
  if (inherits(tr, "multiPhylo") && length(tr) == 1L) {
    tr <- tr[[1L]]
  }
  n <- pool_size(tr)
  if (n == 0L) {
    stop("No tree could be read from: ", basename(f), call. = FALSE) # nocov
  }
  if (n > pool_max) {
    stop(
      basename(f), " contains ", n, " trees, which exceeds the limit of 100. ",
      "This is usually a posterior sample, a bootstrap set, or trees pooled ",
      "across several analytical conditions. phylorug treats each file as one ",
      "analysis. Summarise these trees (e.g. into a consensus) before use.",
      call. = FALSE
    )
  }
  no_labels <- all(vapply(
    as_pool(tr),
    function(x) is.null(x$node.label),
    logical(1)
  ))
  if (no_labels && has_beast_annotations(f)) {
    message(
      basename(f), " contains BEAST-style node annotations. ape reads ",
      "topology only, so posterior probabilities are not imported. The tree ",
      "is fully usable for presence/absence rugs. To import support values, ",
      "use treeio::read.beast() and convert with as.phylo()."
    )
  }
  tr
}
#' Coerce one analysis to a pool of trees
#'
#' Internal. Every function that iterates over the trees of an analysis calls
#' this first, so that single-tree and multi-tree analyses share one code path.
#' A `"phylo"` becomes a `"multiPhylo"` of length one; a `"multiPhylo"`
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
#' Detect the format of a tree file
#' Returns `"nexus"`, `"newick"`, or `"none"` when the file
#' contains no tree and should be skipped.
#' @noRd
detect_format <- function(path) {
  txt <- readLines(path, warn = FALSE)
  if (length(txt) == 0L || all(!nzchar(trimws(txt)))) {
    stop("File appears to be empty: ", basename(path), call. = FALSE)
  }
  head_lines <- utils::head(txt, 50L)
  if (any(grepl("^\\s*#NEXUS", head_lines, ignore.case = TRUE))) {
    if (!any(grepl("BEGIN\\s+TREES", txt, ignore.case = TRUE))) {
      return("none")
    }
    return("nexus")
  }
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
#' `[&posterior=0.98,rate=1.01]`. The Newick grammar treats brackets as
#' comments, so ape strips them during parsing and `node.label` comes back
#' `NULL`. The annotations can only be found by searching the raw file text,
#' not the parsed tree.
#'
#' @noRd
has_beast_annotations <- function(path) {
  txt <- readLines(path, warn = FALSE)
  any(grepl("[&", txt, fixed = TRUE))
}
