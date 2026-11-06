run_clinical <- function(output) {
  env <- new.env(parent = environment())
  env$output <- output
  sys.source(file.path(root, "R/clinical_consequences.R"), env)
  result <- read.csv(file.path(output, "protocol_consequences.csv"))
  compare_numeric_table(result, read_reference_csv("protocol_consequences.csv"),
                         c("rule", "prior"), 1e-8)
  thresholds <- read.csv(file.path(output, "benefit_thresholds.csv"))
  compare_numeric_table(thresholds, read_reference_csv("benefit_thresholds.csv"),
                         "specification", 1e-8)
  invisible(result)
}

practical_impact <- function(output) {
  paired <- read.csv(file.path(output, "paired_changes.csv"))
  selected <- paired[paired$rule %in% c("central_sparse90", "high_sparse90"), ]
  stopifnot(max(abs(selected$corrected - selected$lost - selected$net_PCS_gain)) < 1e-10)
  selection <- selected[c("prior", "rule")]
  for (name in c("changed", "corrected", "lost", "other", "net_PCS_gain"))
    selection[[paste0(name, "_per_1000_trials")]] <- 10 * selected[[name]]
  write.csv(selection, file.path(output, "practical_selection.csv"), row.names = FALSE)
  compare_numeric_table(selection, read_reference_csv("practical_selection.csv"),
                         c("prior", "rule"), 1e-8)
  saved <- if (file.exists(file.path(output, "clinical_exact_gains.json")))
    jsonlite::read_json(file.path(output, "clinical_exact_gains.json"), simplifyVector = TRUE)$results else
    read_reference_json("cost_audited_gains.json")$results
  result <- read.csv(file.path(output, "protocol_consequences.csv"))
  settings <- list(
    MN_fixed = list("primary_MN", "MN_0.5", TRUE),
    MN_flexible = list("expanded_MN", "MN_0.5", TRUE),
    DD_fixed = list("primary_DD", "DD_0.5", TRUE),
    DD_flexible = list("expanded_DD", "DD_0.5", TRUE),
    MN_preferences = list("uncertain_MN", c("MN_0.25", "MN_0.5", "MN_1"), FALSE),
    sources_preferences = list("uncertain_sources", names(clinical$priors), FALSE))
  burden <- expansion <- list()
  for (key in names(settings)) {
    spec <- settings[[key]]
    z <- saved[saved$rule == spec[[1]] & saved$prior %in% spec[[2]], ]
    if (spec[[3]]) z <- z[as.numeric(z$lambda) == 1 &
      as.numeric(z$underdose) == .25 & as.numeric(z$none) == .5, ]
    at <- match(paste(z$rule, z$prior), paste(result$rule, result$prior))
    probability <- result$change_percent[at] / 100
    lower <- as.numeric(z$lower_pp)
    limit <- lower / probability
    burden[[key]] <- data.frame(
      specification = key, rule = spec[[1]],
      review_cost_limit_score_points = min(limit),
      binding_prior = z$prior[which.min(limit)],
      available_fixed_burden_at_cost_5 = min(lower - 5 * probability))
    if (key %in% c("MN_flexible", "DD_flexible", "sources_preferences")) {
      rows <- result[result$rule == spec[[1]] & result$prior %in% spec[[2]], ]
      expansion[[key]] <- data.frame(
        specification = key, rule = rows$rule, prior = rows$prior,
        changes_per_1000_trials = 10 * rows$change_percent,
        patient_change_per_1000_trials = 1000 * rows$patient_change,
        DLT_change_per_1000_trials = 1000 * rows$DLT_change)
    }
  }
  burden <- do.call(rbind, burden)
  expansion <- do.call(rbind, expansion)
  rownames(burden) <- rownames(expansion) <- NULL
  write.csv(burden, file.path(output, "practical_review_cost.csv"), row.names = FALSE)
  write.csv(expansion, file.path(output, "practical_expansion.csv"), row.names = FALSE)
  compare_numeric_table(burden, read_reference_csv("practical_review_cost.csv"),
                         "specification", 1e-8)
  compare_numeric_table(expansion, read_reference_csv("practical_expansion.csv"),
                         c("specification", "rule", "prior"), 1e-8)
  invisible(list(selection = selection, burden = burden, expansion = expansion))
}
