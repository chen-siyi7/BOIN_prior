script_root <- function() {
  script <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script) != 1L) stop("Run this entry point with Rscript.")
  dirname(dirname(normalizePath(sub("^--file=", "", script))))
}

arguments <- function(allowed = character()) {
  args <- commandArgs(TRUE)
  if (length(args) %% 2L) stop("Use named arguments followed by their values.")
  if (!length(args)) return(list())
  keys <- args[seq(1, length(args), by = 2L)]
  if (any(!startsWith(keys, "--")) || anyDuplicated(keys)) stop("Invalid arguments.")
  names <- substring(keys, 3L)
  if (length(allowed) && any(!names %in% allowed)) stop("Unknown argument: ",
    paste(setdiff(names, allowed), collapse = ", "))
  setNames(as.list(args[seq(2, length(args), by = 2L)]), names)
}

check_inventory <- function(root) {
  manifest <- read.csv(file.path(root, "MANIFEST.csv"), stringsAsFactors = FALSE)
  stopifnot(all(file.exists(file.path(root, manifest$path))))
  size <- file.info(file.path(root, manifest$path))$size
  actual <- vapply(file.path(root, manifest$path), function(path) {
    digest::digest(file = path, algo = "sha256")
  }, character(1))
  if (!all(size == manifest$bytes) || !identical(unname(actual), manifest$sha256))
    stop("A distributed file differs from the release inventory.")
  paths <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  paths <- paths[!grepl("^\\.git/|^outputs/", paths)]
  paths <- paths[!basename(paths) %in% c("MANIFEST.csv", ".DS_Store", ".Rhistory", ".RData")]
  paths <- paths[!grepl("\\.(o|so|dll|dylib)$", paths)]
  paths <- paths[!file.info(file.path(root, paths))$isdir]
  if (!setequal(paths, manifest$path)) stop("Distributed files and inventory entries differ.")
  nrow(manifest)
}
