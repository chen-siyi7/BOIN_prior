
verify_clinical_gains <- function(output) {
  model <- clinical$model
  rules <- clinical$policies
  baseline <- model$reference
  probabilities <- lapply(1:3, function(j) gmp::as.bigq(model$planning[, j]))
  weights <- lapply(clinical$priors, function(w) {
    value <- gmp::as.bigq(w)
    value / sum(value)
  })
  category <- lapply(probabilities, function(p) {
    ifelse(p >= gmp::as.bigq("3/10") & p <= gmp::as.bigq("2/5"), 1L,
           ifelse(p > gmp::as.bigq("2/5"), 2L, 3L))
  })
  category[[4]] <- rep(4L, nrow(model$planning))
  column <- function(action) if (action > 0) action else 4L
  components <- lapply(rules, function(a) {
    lapply(weights, function(w) gmp::as.bigq(rep(0, 4L)))
  })
  states <- sort(unique(unlist(lapply(rules, function(a) which(a != baseline)))))
  for (state in states) {
    if (match(state, states) %% 8L == 1L) {
      cat("Clinical gain states:", match(state, states), "of", length(states), "\n")
      flush.console()
    }
    mass <- gmp::as.bigq(rep(model$catalog$path_coefficient[state], nrow(model$planning)))
    for (j in 1:3) {
      n <- model$catalog[[paste0("n", j)]][state]
      x <- model$catalog[[paste0("x", j)]][state]
      mass <- mass * probabilities[[j]]^x * (1 - probabilities[[j]])^(n - x)
    }
    before <- category[[column(baseline[state])]]
    for (prior in names(weights)) {
      weighted <- mass * weights[[prior]]
      for (rule in names(rules)) {
        if (rules[[rule]][state] == baseline[state]) next
        after <- category[[column(rules[[rule]][state])]]
        changed <- after != before
        for (component in 1:4) {
          positive <- which(changed & after == component)
          negative <- which(changed & before == component)
          if (length(positive)) components[[rule]][[prior]][component] <-
            components[[rule]][[prior]][component] + sum(weighted[positive])
          if (length(negative)) components[[rule]][[prior]][component] <-
            components[[rule]][[prior]][component] - sum(weighted[negative])
        }
      }
    }
  }
  preferences <- rbind(clinical$preferences, data.frame(lambda = 1, underdose = .25, none = .5))
  results <- objectives <- list()
  for (rule in names(rules)) {
    relevant <- list()
    settings <- clinical$settings[[rule]]
    for (prior in names(weights)) for (i in seq_len(nrow(preferences))) {
      p <- preferences[i, ]
      score <- c(gmp::as.bigq(1), -decimal_rational(c(p$lambda, p$underdose, p$none)))
      gain <- 100 * sum(components[[rule]][[prior]] * score)
      interval <- enclose_rational(gain)
      results[[length(results) + 1L]] <- data.frame(rule = rule, prior = prior,
        lambda = p$lambda, underdose = p$underdose, none = p$none,
        lower_pp = sprintf("%.17g", interval$lower), upper_pp = sprintf("%.17g", interval$upper))
      if (prior %in% settings$priors &&
          (if (settings$robust) i < nrow(preferences) else i == nrow(preferences))) {
        relevant[[length(relevant) + 1L]] <- gain
      }
    }
    worst <- min(do.call(c, relevant))
    interval <- enclose_rational(worst)
    objectives[[rule]] <- data.frame(rule = rule,
      objective_lower_pp = sprintf("%.17g", interval$lower),
      objective_upper_pp = sprintf("%.17g", interval$upper))
  }
  results <- do.call(rbind, results)
  objectives <- do.call(rbind, objectives)
  reference <- read_reference_json("cost_audited_gains.json")
  for (name in c("lambda", "underdose", "none", "lower_pp", "upper_pp"))
    reference$results[[name]] <- as.numeric(reference$results[[name]])
  compare_numeric_table(transform(results, lower_pp = as.numeric(lower_pp), upper_pp = as.numeric(upper_pp)),
    transform(reference$results, lower_pp = as.numeric(lower_pp), upper_pp = as.numeric(upper_pp)),
    c("rule", "prior", "lambda", "underdose", "none"), 1e-10)
  jsonlite::write_json(list(objectives = objectives, results = results,
    scope = "Exact normalized saved binary64 finite priors and rational preference vertices."),
    file.path(output, "clinical_exact_gains.json"), pretty = TRUE, auto_unbox = TRUE, digits = 17)
  write.csv(results, file.path(output, "clinical_exact_gains.csv"), row.names = FALSE)
  objectives
}
