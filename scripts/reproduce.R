script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- script_root()
source(file.path(root, "R/project.R"))
args <- arguments(c("output", "mode"))
if (is.null(args$output)) stop("Supply --output with a new directory outside this code folder.")
mode <- if (is.null(args$mode)) "all" else args$mode
if (!mode %in% c("all", "verify", "evaluate", "clinical", "gains"))
  stop("Mode must be all, verify, evaluate, clinical or gains.")
output <- boin_output(args$output, root)
started <- proc.time()[["elapsed"]]
timing <- list()
phase <- function(name, calculation) {
  before <- proc.time()[["elapsed"]]
  result <- force(calculation)
  timing[[length(timing) + 1L]] <<- data.frame(phase = name,
    seconds = proc.time()[["elapsed"]] - before)
  write.csv(do.call(rbind, timing), file.path(output, "timing.csv"), row.names = FALSE)
  invisible(result)
}
project <- phase("Input loading and native compilation", boin_project(root, native = TRUE))
checks <- list()
if (mode %in% c("all", "verify")) {
  cat("Rebuilding catalogs, curve generators and partition coverage.\n")
  checks$catalogs <- nrow(phase("Catalogs and input checks", project$verify_inputs(output)))
  checks$point_enclosures <- nrow(phase("Point likelihood enclosure checks", project$verify_point_enclosures(output)))
  cat("Checking core and clinical risks with native integer arithmetic.\n")
  checks$core_inequalities <- nrow(phase("Core risk and endpoint verification", project$verify_core_risks(output)))
  checks$clinical_inequalities <- nrow(phase("Clinical risk verification", project$verify_clinical_risks(output)))
  checks$continuation <- phase("Exact expansion check", project$verify_continuation(output))
}
if (mode %in% c("all", "gains")) {
  cat("Recomputing finite-prior gains with exact rational arithmetic.\n")
  gains <- phase("Exact core gains and table sizes", project$verify_exact_gains(output))
  checks$gain_bounds <- nrow(gains$gains)
  checks$minimum_tables <- length(gains$minimum)
  checks$budget_retention <- length(gains$retention)
  cat("Recomputing clinical score gains with exact rational arithmetic.\n")
  checks$clinical_objectives <- nrow(phase("Exact clinical score gains", project$verify_clinical_gains(output)))
  cat("Checking nominal gain ceilings with independent likelihood enclosures.\n")
  checks$gain_ceilings <- nrow(phase("Independent nominal gain ceilings", project$verify_gain_ceilings(output)))
  checks$family_gain_rows <- nrow(phase("Conditional-family gain bounds", project$verify_family_gains(output)))
  checks$robust_gain_ceilings <- nrow(phase("Independent robust gain ceilings", project$verify_robust_ceilings(output)))
}
if (mode %in% c("all", "evaluate")) {
  cat("Recalculating the independent evaluation.\n")
  tables <- phase("Independent evaluation curves", project$evaluate_rules(output))
  for (name in names(tables)) {
    keys <- if (name == "summary") c("category", "prior", "rule") else
      if (name == "positions") c("rule", "label") else c("prior", "rule")
    project$compare_numeric_table(tables[[name]], project$read_reference_csv(paste0(name, ".csv")), keys)
  }
  checks$evaluation_curves <- nrow(readRDS(file.path(root, "data/evaluation.rds"))$curves)
  checks$generator_intervals <- nrow(phase("Paired generator intervals", project$generator_intervals(output)))
  checks$secondary <- phase("Matched methods and supporting designs", project$verify_secondary(output))
  checks$bound_comparisons <- nrow(phase("Risk-bound comparisons", project$bound_comparisons(output)))
}
if (mode %in% c("all", "clinical")) {
  cat("Checking continuous historical prior sensitivity in R.\n")
  phase("Continuous historical prior sensitivity", project$continuous_prior_check(output))
  checks$continuous_prior <- TRUE
  cat("Recalculating protocol-linked clinical consequences.\n")
  phase("Clinical outcomes and expansion consequences", project$run_clinical(output))
  checks$clinical_consequences <- TRUE
}
if (mode == "all") {
  phase("Practical benefit translations", project$practical_impact(output))
  checks$practical_impact <- TRUE
  checks$logistic_priors <- nrow(phase("Logistic posterior and calibration checks", project$verify_logit_comparators(output)))
  checks$adaptive_comparisons <- nrow(phase("Adaptive ordered-likelihood and CSE replay", project$run_cse(output)))
}
jsonlite::write_json(list(mode = mode, seconds = proc.time()[["elapsed"]] - started,
  checks = checks, R = as.character(getRversion()), platform = R.version$platform,
  package_versions = setNames(lapply(c("Rcpp", "gmp", "jsonlite", "digest"),
    function(package) as.character(packageVersion(package))), c("Rcpp", "gmp", "jsonlite", "digest"))),
  file.path(output, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 17)
cat("Completed", mode, "in", round(proc.time()[["elapsed"]] - started, 2), "seconds.\n")
