
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- script_root()
paths <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
paths <- paths[!grepl("^\\.git/|^outputs/", paths)]
paths <- paths[!basename(paths) %in% c("MANIFEST.csv", ".DS_Store", ".Rhistory", ".RData")]
paths <- paths[!grepl("\\.(o|so|dll|dylib)$", paths)]
paths <- sort(paths[!file.info(file.path(root, paths))$isdir])
inventory <- data.frame(path = paths, bytes = file.info(file.path(root, paths))$size,
  sha256 = vapply(file.path(root, paths), function(path) digest::digest(file = path, algo = "sha256"), ""))
write.csv(inventory, file.path(root, "MANIFEST.csv"), row.names = FALSE)
stopifnot(check_inventory(root) == length(paths))
cat("Updated inventory:", length(paths), "files.\n")
