packages <- c("Rcpp", "gmp", "jsonlite", "digest", "Matrix", "highs")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
for (package in packages) cat(package, as.character(packageVersion(package)), "\n")
