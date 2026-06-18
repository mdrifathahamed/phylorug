# detect Format automatically _used in read_trees_from_dir()
detect_format <- function(path) {
  lines <- readLines(path, n = 5, warn = FALSE)

  # Empty file check
  if (length(lines) == 0 || all(!nzchar(trimws(lines)))) {
    stop("File appears to be empty: ", basename(path),
      call. = FALSE
    )
  }

  # NEXUS -- always starts with #NEXUS
  if (grepl("^\\s*#NEXUS", lines[1], ignore.case = TRUE)) {
    return("nexus")
  }

  if (any(grepl("\\[&&NHX:", lines, ignore.case = TRUE))) {
    return("newick")
  }

  # Default - Newick
  "newick"
}
# normalize_support()
# Rescale a vector of support values to the 0-1 range.
# Bootstrap-style values (0-100) are divided by 100; values already on a
# 0-1 scale are returned unchanged. Detection is per vector: if any value
# exceeds 1, the whole vector is treated as 0-100. All-NA vectors pass
# through untouched. Called at colour time, once per tree column.
normalize_support <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  if (any(x > 1, na.rm = TRUE)) {
    x / 100
  } else {
    x
  }
}
