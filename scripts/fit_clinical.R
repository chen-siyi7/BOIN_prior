script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- script_root()
source(file.path(root, "R/project.R"))
args <- arguments(c("output", "rule", "seconds"))
if (is.null(args$output)) stop("Supply --output with a new directory.")
rule <- if (is.null(args$rule)) "primary_MN" else args$rule
seconds <- if (is.null(args$seconds)) 10 else as.numeric(args$seconds)
if (!is.finite(seconds) || seconds <= 0) stop("Seconds must be positive.")
project <- boin_project(root, native = TRUE)
clinical <- project$clinical
if (!rule %in% names(clinical$settings)) stop("Unknown clinical rule.")
output <- boin_output(args$output, root)
model <- clinical$model
partition <- clinical$partition
settings <- clinical$settings[[rule]]
probabilities <- model$planning
mass <- project$exact_mass(model$catalog, probabilities)
components <- list(A = abs(probabilities - .35) <= .05,
                   R = probabilities > .40, B = probabilities < .30)
preferences <- if (settings$robust) clinical$preferences else
  data.frame(lambda = 1, underdose = .25, none = .5)
objectives <- list()
for (prior in settings$priors) {
  weighted <- sweep(mass, 2, clinical$priors[[prior]], "*")
  terms <- lapply(components, function(component) weighted %*% component)
  for (i in seq_len(nrow(preferences))) {
    p <- preferences[i, ]
    objectives[[paste(prior, i, sep = "/")]] <- cbind(
      terms$A - p$lambda * terms$R - p$underdose * terms$B, -p$none * rowSums(weighted))
  }
}
j <- model$design$J
bounds <- project$box_dp_extrema(as.matrix(model$catalog[paste0("n", 1:j)]),
  as.matrix(model$catalog[paste0("x", 1:j)]), model$catalog$log_path_coefficient,
  partition$low, partition$up, as.integer(seq_len(nrow(model$catalog))))
locked <- if (settings$locked) model$reference %in% c(1L, j) else rep(FALSE, length(model$reference))
fit <- project$clinical_cost_fit(model, partition, bounds$U, bounds$L, objectives,
  clinical$specs[[rule]], locked, settings$none, clinical$policies, output, seconds = seconds)
audit <- project$audit_clinical(model, partition, fit$action, fit$spec)
write.csv(audit, file.path(output, "risk_audit.csv"), row.names = FALSE)
write.csv(data.frame(state = seq_along(fit$action), action = fit$action),
          file.path(output, "saved_actions.csv"), row.names = FALSE)
saveRDS(fit, file.path(output, "fit.rds"))
jsonlite::write_json(list(rule = rule, changed = fit$changed, gain_pp = fit$gain,
  risk_check = "passed", seconds = fit$seconds), file.path(output, "fit_summary.json"),
  pretty = TRUE, auto_unbox = TRUE, digits = 17)
cat("Fitted", rule, "with", fit$changed, "changes and an independent risk check.\n")
