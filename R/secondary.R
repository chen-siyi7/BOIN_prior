
secondary_metrics <- function(model, curves, policies) {
  j <- model$design$J
  probabilities <- as.matrix(curves[paste0("p", seq_len(j))])
  result <- lapply(policies, function(a) matrix(0, nrow(curves), 7L,
    dimnames = list(NULL, c("PCS", "A", "R", "B", "N", "DLT", "changed"))))
  for (first in seq(1, nrow(curves), by = 128L)) {
    ids <- first:min(first + 127L, nrow(curves))
    p <- probabilities[ids, , drop = FALSE]
    mass <- exact_mass(model$catalog, p)
    labels <- curves$mtd[ids]
    for (key in names(policies)) {
      a <- policies[[key]]
      selected <- vapply(seq_len(j), function(dose)
        colSums(mass[a == dose, , drop = FALSE]), numeric(length(ids)))
      result[[key]][ids, ] <- cbind(
        PCS = rowSums(selected * (col(selected) == labels)),
        A = rowSums(selected * (abs(p - model$design$phi) <= .05)),
        R = rowSums(selected * (p > model$design$phi + .05)),
        B = rowSums(selected * (p < model$design$phi - .05)),
        N = colSums(mass[a < 0L, , drop = FALSE]),
        DLT = rowSums(selected * p),
        changed = colSums(mass[a != model$reference, , drop = FALSE]))
    }
  }
  result
}

secondary_risk <- function(model, partition, policies) {
  j <- model$design$J
  n <- as.matrix(model$catalog[paste0("n", seq_len(j))])
  x <- as.matrix(model$catalog[paste0("x", seq_len(j))])
  do.call(rbind, lapply(names(policies), function(key) {
    a <- policies[[key]]
    check_actions(model, a)
    states <- which(a != model$reference)
    h <- policy_contrast(a[states], model$reference[states], j)
    z <- fixedpoint_interval_audit(n, x, as.double(model$catalog$path_coefficient),
      partition$low, partition$up, partition$feasible, as.integer(states), h,
      as.integer(round(200 * model$design$phi)), 100L, 200L, 10000L)
    stopifnot(all(z$pass))
    cbind(rule = key, z)
  }))
}

verify_secondary <- function(output) {
  verify_secondary_generators()
  protected <- secondary$policies[!startsWith(names(secondary$policies), "uncapped")]
  risks <- secondary_risk(o, part, protected)
  write.csv(risks, file.path(output, "secondary_core_risk_checks.csv"), row.names = FALSE)
  comparisons <- secondary_metrics(o, secondary$comparison_curves,
    c(list(BOIN = b), secondary$policies))
  curves <- secondary$comparison_curves
  rows <- list()
  for (library in unique(curves$library)) for (prior in names(core$priors)) {
    ids <- which(curves$library == library)
    w <- core$priors[[prior]][curves$mtd[ids]] / 96
    for (key in names(comparisons)) {
      value <- comparisons[[key]][ids, , drop = FALSE]
      base <- comparisons$BOIN[ids, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(library = library, policy = key, prior = prior,
        baseline_PCS = 100 * sum(w * base[, "PCS"]), PCS = 100 * sum(w * value[, "PCS"]),
        gain_pp = 100 * sum(w * (value[, "PCS"] - base[, "PCS"])),
        worst_PCS_pp = 100 * min(value[, "PCS"] - base[, "PCS"]))
    }
  }
  comparison <- do.call(rbind, rows)
  expected <- secondary$comparison_reference
  compare_numeric_table(comparison[comparison$policy %in% unique(expected$policy), ],
    expected, c("library", "policy", "prior"))
  write.csv(comparison, file.path(output, "matched_comparisons_by_library.csv"), row.names = FALSE)
  pooled <- do.call(rbind, lapply(split(comparison, interaction(comparison$policy, comparison$prior)), function(d)
    data.frame(rule = d$policy[1], prior = d$prior[1], BOIN_PCS = mean(d$baseline_PCS),
      PCS = mean(d$PCS), gain_pp = mean(d$gain_pp), worst_pp = min(d$worst_PCS_pp))))
  rownames(pooled) <- NULL
  write.csv(pooled, file.path(output, "matched_comparisons.csv"), row.names = FALSE)
  saveRDS(list(curves = curves, metrics = comparisons), file.path(output, "matched_evaluation.rds"))

  support_rows <- support_risks <- support_planning <- list()
  for (kind in names(secondary$supporting)) {
    item <- secondary$supporting[[kind]]
    model <- item$model
    j <- model$design$J
    verify_partition(item$partition, j, model$design$phi, if (j == 8L) 5L else 8L)
    support_risks[[kind]] <- cbind(design = kind,
      secondary_risk(model, item$partition, item$policies))
    d <- item$curves
    names(d)[4:(3 + j)] <- paste0("p", seq_len(j))
    z <- secondary_metrics(model, d, c(list(BOIN = model$reference), item$policies))
    planning_mass <- exact_mass(model$catalog, model$planning)
    for (prior in names(model$weights)) for (key in names(item$policies)) {
      action <- item$policies[[key]]
      difference <- vapply(seq_len(nrow(model$planning)), function(c)
        sum(planning_mass[, c] * (as.integer(action == model$labels[c]) -
          as.integer(model$reference == model$labels[c]))), 0.)
      gain <- 100 * sum(model$weights[[prior]] * difference)
      fit <- item$fits[[key]]
      if (!is.null(fit) && endsWith(key, prior)) stopifnot(abs(gain - fit$gain_pp) < 1e-8)
      support_planning[[length(support_planning) + 1L]] <- data.frame(design = kind,
        rule = key, prior = prior, planning_gain_pp = gain)
    }
    for (prior in names(item$q)) for (key in names(item$policies)) {
      ids <- which(d$mtd > 0)
      q <- item$q[[prior]][-1]
      w <- q[d$mtd[ids]] / (96 * sum(q))
      value <- z[[key]][ids, "PCS"]
      base <- z$BOIN[ids, "PCS"]
      support_rows[[length(support_rows) + 1L]] <- data.frame(setting = kind, policy = key,
        evaluation_prior = prior, baseline_PCS = 100 * sum(w * base), PCS = 100 * sum(w * value),
        gain_pp = 100 * sum(w * (value - base)), worst_PCS_pp = 100 * min(value - base))
    }
    saveRDS(list(curves = d, metrics = z), file.path(output, paste0("supporting_", kind, ".rds")))
    lookup <- model$catalog[c(paste0("n", seq_len(j)), paste0("x", seq_len(j)))]
    lookup$BOIN <- model$reference
    for (key in names(item$policies)) lookup[[key]] <- item$policies[[key]]
    write.csv(lookup, file.path(output, paste0("supporting_", kind, "_decisions.csv")), row.names = FALSE)
    cat("Supporting design checked:", kind, "\n")
    flush.console()
  }
  support_table <- do.call(rbind, support_rows)
  compare_numeric_table(support_table, secondary$supporting_reference,
    c("setting", "policy", "evaluation_prior"))
  support_risk_table <- do.call(rbind, support_risks)
  write.csv(support_table, file.path(output, "supporting_designs.csv"), row.names = FALSE)
  write.csv(support_risk_table, file.path(output, "supporting_risk_checks.csv"), row.names = FALSE)
  write.csv(do.call(rbind, support_planning), file.path(output, "supporting_planning_gains.csv"), row.names = FALSE)

  d_budget <- readRDS(file.path(root, "data/evaluation.rds"))$curves
  z_budget <- secondary_metrics(o, d_budget, c(list(BOIN = b), core$budget_policies))
  budget_rows <- lapply(names(core$budget_policies), function(key) {
    prior <- strsplit(key, "_", fixed = TRUE)[[1]][1]
    ids <- which(d_budget$category == "bracketing")
    w <- core$priors[[prior]][d_budget$mtd[ids]] / 160
    value <- z_budget[[key]][ids, , drop = FALSE]
    base <- z_budget$BOIN[ids, , drop = FALSE]
    data.frame(rule = key, PCS = 100 * sum(w * value[, "PCS"]),
      gain = 100 * sum(w * (value[, "PCS"] - base[, "PCS"])),
      worst = 100 * min(value[, "PCS"] - base[, "PCS"]),
      changed = 100 * sum(w * value[, "changed"]))
  })
  budget_table <- do.call(rbind, budget_rows)
  compare_numeric_table(budget_table, read_reference_csv("budget_evaluation.csv"), "rule")
  write.csv(budget_table, file.path(output, "budget_evaluation.csv"), row.names = FALSE)

  robust_actions <- c(list(BOIN = b, central = core$policies$central_original,
    high = core$policies$high_original), secondary$policies[c("high_mode", "high_mean", "high_family", "joint_family")])
  d <- secondary$robust$evaluation$curves
  z <- secondary_metrics(o, d, robust_actions)
  for (key in names(z)) stopifnot(max(abs(z[[key]] - secondary$robust$evaluation$metrics[[key]])) < 1e-10)
  saveRDS(list(curves = d, metrics = z), file.path(output, "family_evaluation.rds"))
  clinical_evaluation <- readRDS(file.path(root, "data/secondary/clinical_evaluation.rds"))
  clinical_metrics <- secondary_metrics(clinical$model, clinical_evaluation$curves,
    c(list(BOIN = clinical$model$reference), clinical$policies))
  for (key in names(clinical_metrics))
    stopifnot(max(abs(clinical_metrics[[key]] - clinical_evaluation$metrics[[key]])) < 1e-10)
  saveRDS(list(curves = clinical_evaluation$curves, metrics = clinical_metrics),
    file.path(output, "clinical_independent_evaluation.rds"))
  list(secondary_risk_inequalities = nrow(risks), supporting_risk_inequalities = nrow(support_risk_table),
    matched_rows = nrow(pooled), supporting_rows = nrow(support_table),
    supporting_planning_rows = length(support_planning), budget_rows = nrow(budget_table),
    family_evaluation_curves = nrow(d), clinical_evaluation_curves = nrow(clinical_evaluation$curves))
}
