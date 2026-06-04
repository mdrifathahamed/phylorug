#' Read multiple phylogenetic trees from a directory
#'
#' Searches a specified directory for phylogenetic tree files with a matching
#' file extension and imports them into R as a named list of tree objects.
#' Supports both Newick and NEXUS formats.
#'
#' @param dir Character string. Path to the directory containing the tree files.
#'
#' @param ext Character string. File extension to search for, without the
#'  leading dot. Defaults to `"tre"`.
#'
#' @param format Character string. File format of the trees. Use `"auto"`
#'   (default) to detect format from file extension, `"newick"` to force Newick
#'   parsing, or `"nexus"` to force
#'   NEXUS parsing.
#' @param verbose Logical. If `TRUE` (default), prints a message reporting how
#'   many trees were read successfully. Set to `FALSE` for silent operation
#'   inside pipelines or loops.
#'
#' @return A named list of objects of class `"phylo"`. List names are the
#'  filenames the trees were read from. Elements for files that could not be
#'  read are set to `NULL`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' trees <- read_trees_from_dir(
#'   dir    = "path/to/your/trees",
#'   ext    = "tre",
#'   format = "newick"
#' )
#' }
read_trees_from_dir <- function(dir,
                                ext = "tre",
                                format = c("auto", "newick", "nexus"),
                                verbose = TRUE) {
  format <- match.arg(format)
  # Validate directory
  if (!dir.exists(dir)) {
    stop("`dir` does not exist: ", dir, call. = FALSE)
  }

  # List matching files
  files <- list.files(
    dir,
    pattern    = paste0("\\.", ext, "$"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No .", ext, " files found in: ", dir, call. = FALSE)
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

  # Report how many files failed to read
  failed <- sum(vapply(trees, is.null, logical(1)))
  if (failed > 0) {
    warning(
      failed, " file(s) could not be read and were set to NULL. ",
      "Run check_same_taxa() before proceeding.",
      call. = FALSE
    )
  }

  trees
}
