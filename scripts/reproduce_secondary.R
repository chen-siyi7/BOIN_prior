script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- script_root()
source(file.path(root, "R/project.R"))
args <- arguments(c("output", "mode"))
if (is.null(args$output)) stop("Supply --output with a new directory outside the code folder.")
mode <- if (is.null(args$mode)) "all" else args$mode
stopifnot(mode %in% c("all", "evaluate", "gains", "logit", "cse", "benchmark"))
output <- boin_output(args$output, root)
started <- proc.time()[["elapsed"]]
timing <- list()
phase <- function(name, calculation) {
  before <- proc.time()[["elapsed"]]
  result <- force(calculation)
  timing[[length(timing) + 1L]] <<- data.frame(phase = name, seconds = proc.time()[["elapsed"]] - before)
  write.csv(do.call(rbind, timing), file.path(output, "timing.csv"), row.names = FALSE)
  invisible(result)
}
project <- phase("Input loading and native compilation", boin_project(root, native = TRUE))
checks <- list()
if (mode %in% c("all", "evaluate")) {
  checks$secondary <- phase("Secondary risks and outcomes", project$verify_secondary(output))
  checks$bound_comparisons <- nrow(phase("Risk-bound comparisons", project$bound_comparisons(output)))
  checks$intervals <- nrow(phase("Generator intervals", project$generator_intervals(output)))
}
if (mode %in% c("all", "gains")) {
  checks$family_gains <- nrow(phase("Conditional-family gain bounds", project$verify_family_gains(output)))
  checks$robust_ceilings <- nrow(phase("Robust gain ceilings", project$verify_robust_ceilings(output)))
}
if (mode %in% c("all", "logit"))
  checks$logistic <- nrow(phase("Logistic posterior and calibration", project$verify_logit_comparators(output)))
if (mode %in% c("all", "cse", "benchmark"))
  checks$CSE <- nrow(phase("Adaptive comparison", project$run_cse(output, replay = mode != "benchmark")))
jsonlite::write_json(list(mode = mode, checks = checks, seconds = proc.time()[["elapsed"]] - started,
  R = as.character(getRversion()), platform = R.version$platform),
  file.path(output, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Secondary calculations completed.\n")
