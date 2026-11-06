
verify_gain_ceilings <- function(output) {
  inputs <- readRDS(file.path(root, "data/gain_ceilings.rds"))
  rows <- list()
  add_bound <- function(name, value, saved) {
    interval <- enclose_rational(value)
    stopifnot(value >= 0, abs(interval$upper - as.numeric(saved)) < 1e-7)
    rows[[name]] <<- data.frame(rule = name, lower_pp = sprintf("%.17g", interval$lower),
      upper_pp = sprintf("%.17g", interval$upper), saved_upper_pp = as.numeric(saved), pass = TRUE)
  }
  sum_by_state <- function(value, state) {
    sum(do.call(c, lapply(split(seq_along(state), state), function(ids)
      max(c(gmp::as.bigq(0), value[ids])))))
  }
  validate_point <- function(point, metric, phi) {
    p <- gmp::as.bigq(point)
    j <- length(p)
    stopifnot(all(p >= 0 & p <= 1), all(p[-1] >= p[-j]))
    if (metric$kind %in% c(1L, 3L)) {
      dose <- metric$index
      if (dose > 1L) stopifnot(p[dose - 1L] + p[dose] <= 2 * phi)
      if (dose < j) stopifnot(p[dose] + p[dose + 1L] >= 2 * phi)
      if (metric$kind == 3L) stopifnot(p[1] <= phi, p[j] >= phi)
    }
    if (metric$kind == 2L) {
      l <- metric$left
      u <- metric$right
      low <- phi - gmp::as.bigq("1/20")
      high <- phi + gmp::as.bigq("1/20")
      if (l > 0) stopifnot(p[l] <= low)
      if (l < j) stopifnot(p[l + 1L] >= low)
      if (u > 0) stopifnot(p[u] <= high)
      if (u < j) stopifnot(p[u + 1L] >= high)
    }
    invisible(TRUE)
  }
  candidate_key <- function(state, action) paste(state, action, sep = "/")
  expected <- which(o$eligible & b > 0 & !b %in% c(1L, J) &
    !seq_along(b) %in% anchors & outer(b, seq_len(J), "!="), arr.ind = TRUE)
  candidates <- inputs$core$candidates
  stopifnot(!anyDuplicated(candidate_key(candidates$state, candidates$action)),
    setequal(candidate_key(expected[, 1], expected[, 2]),
             candidate_key(candidates$state, candidates$action)),
    all(candidates$reference == b[candidates$state]))
  states <- sort(unique(candidates$state))
  cat("Auditing the nominal planning objective on", length(states), "candidate states.\n")
  flush.console()
  conditional <- exact_conditional_masses(states)
  indices <- sort(unique(unlist(lapply(inputs$core$prices, function(p) p$point + 1L))))
  mass <- point_mass_bounds(nn, xx, as.double(o$catalog$path_coefficient),
                            inputs$core$grid[indices, , drop = FALSE])
  mass_low <- gmp::as.bigq(as.vector(mass$low))
  mass_up <- gmp::as.bigq(as.vector(mass$up))
  for (prior in names(inputs$core$prices)) {
    q <- exact_location_prior(prior)
    reduced <- gmp::as.bigq(rep(0, nrow(candidates)))
    at <- match(candidates$state, states)
    for (j in seq_len(J)) reduced <- reduced + 100 * q[j] * conditional[[j]][at] *
      (as.integer(candidates$action == j) - as.integer(candidates$reference == j))
    prices <- inputs$core$prices[[prior]]
    constant <- gmp::as.bigq(0)
    for (i in seq_len(nrow(prices))) {
      price <- gmp::as.bigq(prices$price[i])
      stopifnot(price >= 0)
      point <- prices$point[i] + 1L
      tag <- prices$metric[i]
      j <- as.integer(sub("^(C|Tp|Tm)", "", tag))
      metric <- list(kind = if (startsWith(tag, "C")) 1L else 0L, index = j)
      validate_point(inputs$core$grid[point, ], metric, gmp::as.bigq("1/4"))
      h <- if (startsWith(tag, "C"))
        as.integer(candidates$reference == j) - as.integer(candidates$action == j) else
        (as.integer(candidates$action > j) - as.integer(candidates$reference > j)) *
        if (startsWith(tag, "Tp")) 1L else -1L
      index <- candidates$state + K * (match(point, indices) - 1L)
      bound <- mass_low[index]
      bound[h < 0] <- mass_up[index[h < 0]]
      reduced <- reduced - 50 * price * h * bound
      constant <- constant + price
    }
    upper <- constant + sum_by_state(reduced, candidates$state)
    reference <- inputs$core$reference$results
    add_bound(paste0("core_", prior), upper, reference$decimal_upper_pp[reference$prior == prior])
    cat("Nominal gain ceiling checked for", prior, "prior.\n")
    flush.console()
  }
  model <- clinical$model
  j <- model$design$J
  k <- length(model$reference)
  source_names <- c("MN_0.5", "DD_0.5")
  weights <- lapply(clinical$priors[source_names], function(w) {
    exact <- gmp::as.bigq(w)
    enclose_rational(exact / sum(exact))
  })
  utilities <- vapply(1:j, function(dose) {
    p <- gmp::as.bigq(model$planning[, dose])
    as.integer(ifelse(p >= gmp::as.bigq("3/10") & p <= gmp::as.bigq("2/5"), 4L,
                      ifelse(p > gmp::as.bigq("2/5"), -4L, -1L)))
  }, integer(nrow(model$planning)))
  score <- clinical_score_bounds(as.matrix(model$catalog[paste0("n", 1:j)]),
    as.matrix(model$catalog[paste0("x", 1:j)]), as.double(model$catalog$path_coefficient),
    model$planning, do.call(cbind, lapply(weights, `[[`, "lower")),
    do.call(cbind, lapply(weights, `[[`, "upper")), cbind(utilities, -2L))
  score_low <- gmp::as.bigq(as.vector(score$low))
  score_up <- gmp::as.bigq(as.vector(score$up))
  indices <- sort(unique(unlist(lapply(inputs$clinical$rules, function(r) r$prices$point + 1L))))
  mass <- point_mass_bounds(as.matrix(model$catalog[paste0("n", 1:j)]),
    as.matrix(model$catalog[paste0("x", 1:j)]), as.double(model$catalog$path_coefficient),
    inputs$clinical$points[indices, , drop = FALSE])
  mass_low <- gmp::as.bigq(as.vector(mass$low))
  mass_up <- gmp::as.bigq(as.vector(mass$up))
  reference <- read_reference_json("clinical_gain_bounds.json")$results
  for (rule in names(inputs$clinical$rules)) {
    data <- inputs$clinical$rules[[rule]]
    candidates <- data$candidates
    settings <- clinical$settings[[rule]]
    locked <- settings$locked & model$reference %in% c(1L, j)
    expected <- which(model$eligible & model$reference > 0 & !locked &
      outer(model$reference, 1:j, "!="), arr.ind = TRUE)
    keys <- candidate_key(expected[, 1], expected[, 2])
    if (settings$none) keys <- c(keys, candidate_key(which(model$reference > 0), -1L))
    stopifnot(setequal(keys, candidate_key(candidates$state, candidates$action)),
      !anyDuplicated(candidate_key(candidates$state, candidates$action)),
      all(candidates$reference == model$reference[candidates$state]))
    source <- match(if (endsWith(rule, "MN")) "MN_0.5" else "DD_0.5", source_names)
    column <- ifelse(candidates$action > 0, candidates$action, j + 1L)
    after <- candidates$state + k * (column - 1L) + k * (j + 1L) * (source - 1L)
    before <- candidates$state + k * (candidates$reference - 1L) + k * (j + 1L) * (source - 1L)
    reduced <- 100 * (score_up[after] - score_low[before])
    h_matrix <- contrasts(candidates$action, candidates$reference, data$metrics)
    constant <- gmp::as.bigq(0)
    for (i in seq_len(nrow(data$prices))) {
      m <- data$prices$metric[i] + 1L
      point <- data$prices$point[i] + 1L
      price <- gmp::as.bigq(data$prices$price[i])
      stopifnot(price >= 0)
      validate_point(inputs$clinical$points[point, ], data$metrics[m, ], gmp::as.bigq("7/20"))
      h <- h_matrix[, m]
      index <- candidates$state + k * (match(point, indices) - 1L)
      bound <- mass_low[index]
      bound[h < 0] <- mass_up[index[h < 0]]
      reduced <- reduced - 100 * price * h * bound
      constant <- constant + 100 * decimal_rational(data$metrics$cap[m]) * price
    }
    upper <- constant + sum_by_state(reduced, candidates$state)
    add_bound(rule, upper, reference$directed_upper_pp[reference$rule == rule])
    cat("Clinical gain ceiling checked for", rule, ".\n")
    flush.console()
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  write.csv(result, file.path(output, "gain_ceiling_checks.csv"), row.names = FALSE)
  result
}
