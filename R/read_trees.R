#' Read phylogenetic trees from a directory
#'
#' Looks for phylogenetic tree files with a matching formats and extension
#' inside the assigned directory. Once found, returns them as a named list of
#' `"phylo"` objects. Supports both Newick and NEXUS formats.
#'
#' @param dir A character string specifying the target directory containing the
#'  phylogenetic tree files.
#'
#' @param ext A character string specifying file extension to look for. Defaults
#'  to `"tre"`.
#'
#' @param format A character string specifying file format of the phylogenetic
#'  trees. Use `"auto"` (default) to detect format from file content, `"newick"`
#'  to force Newick parsing, or `"nexus"` to force NEXUS parsing.
#'
#' @param verbose Logical. If `TRUE` (default), prints a message reporting
#'  the number of trees that were read successfully. Set to `FALSE` for silent
#'  operation.
#'
#' @return A named list of `"phylo"` objects. List elements are named after the
#'  filenames of the original phylogenetic tree files. Elements for files that
#'  could not be read are set to `NULL`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simplest usage -- auto-detects format from file content
#' trees <- read_trees("path/to/your/trees")
#'
#' # Specify extension explicitly
#' trees <- read_trees(
#'   dir = "path/to/your/trees",
#'   ext = "nex"
#' )
#'
#' # Force Newick parsing
#' trees <- read_trees(
#'   dir    = "path/to/your/trees",
#'   format = "newick"
#' )
#' }
read_trees <- function(dir,
                       ext     = "tre",
                       format  = c("auto", "newick", "nexus"),
                       verbose = TRUE) {
  format <- match.arg(format)
  # Validate directory
  if (!dir.exists(dir)) {
    stop("`dir` does not exist: ", dir, call. = FALSE)
  }

  # List matching files
  ext_pattern <- paste0(
    "\\.(", paste(ext, collapse = "|"), ")$"
  )

  files <- list.files(
    dir,
    pattern     = ext_pattern,
    full.names  = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {
    stop(
      "No files matching extensions (",
      paste(ext, collapse = ", "),
      ") found in: ", dir,
      call. = FALSE
    )
  }

  # Read each file -- detect format per file when auto
  trees <- lapply(files, function(f) {
    tryCatch(
      {
        actual_reader <- if (format == "auto") {
          fmt <- detect_format(f)
          switch(fmt,
            nexus  = ape::read.nexus,
            newick = ape::read.tree
          )
        } else {
          switch(format,
            newick = ape::read.tree,
            nexus  = ape::read.nexus
          )
        }
        suppressWarnings(actual_reader(f))
      },
      error = function(e) {
        message("Could not read: ", basename(f), " - ", e$message)
        NULL
      }
    )
  })

  # Name list elements by filename
  names(trees) <- basename(files)

  # Report how many trees were read successfully
  if (verbose) {
    n_read <- sum(!vapply(trees, is.null, logical(1)))
    message(
      "Read ", n_read, " of ", length(files),
      " trees from: ", dir
    )
  }

  # Reports which files failed to read
  failed      <- sum(vapply(trees, is.null, logical(1)))
  failed_names <- names(trees)[vapply(trees, is.null, logical(1))]

  if (failed > 0) {
    stop(
      failed, " file(s) could not be read: ",
      paste(failed_names, collapse = ", "),
      call. = FALSE
    )
  }

  trees
}
