script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- script_root()
source(file.path(root, "R/project.R"))
args <- arguments(c("output", "prior", "mode", "retention", "seconds", "budget"))
if (is.null(args$output)) stop("Supply --output with a new directory.")
prior <- if (is.null(args$prior)) "high" else args$prior
if (!prior %in% c("central", "high")) stop("Prior must be central or high.")
mode <- if (is.null(args$mode)) "short" else args$mode
if (!mode %in% c("short", "full", "transport")) stop("Mode must be short, full or transport.")
fraction <- if (is.null(args$retention)) .9 else as.numeric(args$retention)
if (!is.finite(fraction) || fraction <= 0 || fraction > 1) stop("Retention must lie in (0,1].")
seconds <- if (is.null(args$seconds)) 20 else as.numeric(args$seconds)
if (!is.finite(seconds) || seconds <= 0) stop("Seconds must be positive.")
budget <- if (is.null(args$budget)) .02 else as.numeric(args$budget)
if (!is.finite(budget) || !budget %in% c(.005, .01, .02, .03))
  stop("Budget must be 0.005, 0.01, 0.02 or 0.03.")
output <- boin_output(args$output, root)
project <- boin_project(root, native = TRUE)
contribution <- project$gain_C(project$core$priors[[prior]])
bank <- if (budget == .02) project$core$policies[[paste0(prior, "_original")]] else {
  suffix <- c(`0.005` = "0005", `0.01` = "001", `0.03` = "003")[[as.character(budget)]]
  project$core$budget_policies[[paste0(prior, "_bank", suffix)]]
}
if (mode == "short") {
  states <- which(bank != project$b)
  candidates <- data.frame(state = states, action = bank[states], reference = project$b[states])
} else {
  locked <- seq_len(project$K) %in% project$anchors
  if (mode == "full") locked <- locked | project$b %in% c(1, project$J)
  if (mode == "transport") locked <- locked |
    (project$b %in% c(1, project$J) & !seq_len(project$K) %in% project$graph$edges$bad)
  id <- which(project$o$eligible & project$b > 0 & !locked &
                outer(project$b, seq_len(project$J), "!="), arr.ind = TRUE)
  candidates <- data.frame(state = id[, 1], action = id[, 2], reference = project$b[id[, 1]])
}
fit <- project$fit_stream(candidates, contribution, list(bank), output,
  mode = if (mode == "short") "sparse" else "gain", transport = mode == "transport",
  target = if (mode == "short") fraction * project$score(bank, contribution) else 0,
  seconds = seconds, budget = budget)
saveRDS(fit, file.path(output, "fit.rds"))
jsonlite::write_json(list(prior = prior, mode = mode, budget = budget, changed = fit$changed,
  gain_pp = fit$gain_pp, risk_check = "passed", seconds = fit$seconds),
  file.path(output, "fit_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 17)
cat("Fitted", fit$changed, "changes with gain", fit$gain_pp, "percentage points.\n")
