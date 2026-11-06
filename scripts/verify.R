script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- script_root()
source(file.path(root, "R/project.R"))
cat("Checking", check_inventory(root), "distributed files.\n")
project <- boin_project(root)
examples <- read.csv(file.path(root, "tests/lookup_examples.csv"))
for (i in seq_len(nrow(examples))) {
  row <- examples[i, ]
  counts <- as.numeric(strsplit(row$counts, ",", fixed = TRUE)[[1]])
  dlts <- as.numeric(strsplit(row$dlts, ",", fixed = TRUE)[[1]])
  result <- project$select_dose(counts, dlts, row$table)
  stopifnot(result$BOIN == row$BOIN, result$recommendation == row$final)
}
invalid <- try(project$select_dose(rep(0, 6), rep(0, 6)), silent = TRUE)
stopifnot(inherits(invalid, "try-error"))
files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE)
files <- files[!grepl("/\\.git/", files)]
stopifnot(!any(grepl("\\.(py|pyc|pyo|ipynb)$", files)))
other_languages <- c("py", "pyc", "pyo", "ipynb", "jl", "js", "ts", "sh", "bash", "m", "c", "java")
stopifnot(!any(tolower(tools::file_ext(files)) %in% other_languages))
for (path in files[grepl("\\.R$", files)]) {
  parse(path)
}
cat("Lookup examples, unreachable-data rejection and R/Rcpp source checks pass.\n")
