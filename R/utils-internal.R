#detect Format automatically
detect_format <- function(path) {
  lines <- readLines(path, n = 5, warn = FALSE)

  # Empty file check
  if (length(lines) == 0 || all(!nzchar(trimws(lines)))) {
    stop("File appears to be empty: ", basename(path),
         call. = FALSE)
  }

  # NEXUS -- always starts with #NEXUS
  if (grepl("^\\s*#NEXUS", lines[1], ignore.case = TRUE)) {
    return("nexus")
  }

  # NHX -- contains [&&NHX: marker anywhere in first 5 lines
  # Still read with ape::read.tree() since NHX is Newick superset
  if (any(grepl("[&&NHX:", lines, fixed = TRUE))) {
    return("newick")
  }

  # Default -- standard Newick
  "newick"
}
