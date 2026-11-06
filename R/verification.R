check_actions <- function(model, action, locked = TRUE, anchors = TRUE, none = FALSE) {
  reference <- model$reference
  j <- model$design$J
  if (length(action) != length(reference) || any(!action %in% c(-1L, seq_len(j))))
    stop("A rule has invalid actions.")
  stopifnot(all(action[reference < 0] == -1L))
  if (!none) stopifnot(all(action[reference > 0] > 0))
  ids <- which(action > 0)
  stopifnot(all(model$eligible[cbind(ids, action[ids])]))
  if (locked) stopifnot(all(action[reference %in% c(1L, j)] == reference[reference %in% c(1L, j)]))
  if (anchors) {
    fixed <- vapply(seq_len(j), function(dose) which(reference == dose)[1], integer(1))
    stopifnot(all(action[fixed] == reference[fixed]))
  }
  invisible(TRUE)
}

verify_point_enclosures <- function(output) {
  rows <- list()
  for (name in c("core", "clinical")) {
    model <- if (name == "core") core$model else clinical$model
    j <- model$design$J
    states <- unique(c(1L, nrow(model$catalog), which.max(model$catalog$path_coefficient),
                      as.integer(seq(1, nrow(model$catalog), length.out = 12L))))
    ids <- as.integer(seq(1, nrow(model$planning), length.out = 8L))
    points <- rbind(model$planning[ids, , drop = FALSE], rep(0, j), rep(1, j),
                    rep(1e-10, j), seq(.01, .99, length.out = j))
    enrollment <- as.matrix(model$catalog[paste0("n", 1:j)])[states, , drop = FALSE]
    events <- as.matrix(model$catalog[paste0("x", 1:j)])[states, , drop = FALSE]
    coefficient <- model$catalog$path_coefficient[states]
    interval <- point_mass_bounds(enrollment, events, as.double(coefficient), points)
    for (i in seq_len(nrow(points))) {
      mass <- gmp::as.bigq(coefficient)
      for (dose in 1:j) {
        p <- gmp::as.bigq(points[i, dose])
        mass <- mass * p^events[, dose] * (1 - p)^(enrollment[, dose] - events[, dose])
      }
      stopifnot(all(gmp::as.bigq(interval$low[, i]) <= mass),
                all(gmp::as.bigq(interval$up[, i]) >= mass))
    }
    rows[[name]] <- data.frame(design = name, exact_comparisons = length(states) * nrow(points),
                               pass = TRUE)
  }
  result <- do.call(rbind, rows)
  write.csv(result, file.path(output, "point_enclosure_checks.csv"), row.names = FALSE)
  result
}

verify_core_risks <- function(output) {
  rows <- flows <- list()
  for (key in names(core$policies)) {
    action <- core$policies[[key]]
    transport <- grepl("transport", key) || key == "endpoint_pair"
    check_actions(o, action, locked = !transport)
    states <- which(action != b)
    h <- policy_contrast(action[states], b[states], J)
    if (transport) {
      h[, c(1L, J)] <- 0L
      flows[[key]] <- cbind(rule = key, check_endpoint_flow(o, action, graph))
    }
    z <- fixedpoint_interval_audit(nn, xx, as.double(o$catalog$path_coefficient),
      part$low, part$up, part$feasible, as.integer(states), h, 50L, 100L, 200L, 10000L)
    if (transport) z <- z[!z$metric %in% c("C1", paste0("C", J)), ]
    stopifnot(all(z$pass))
    rows[[key]] <- cbind(rule = key, z)
  }
  for (key in names(core$budget_policies)) {
    action <- core$budget_policies[[key]]
    check_actions(o, action)
    stem <- sub(".*_(bank|short)", "", key)
    budget <- as.numeric(stem) / 10^(nchar(stem) - 1L)
    states <- which(action != b)
    h <- policy_contrast(action[states], b[states], J)
    z <- fixedpoint_interval_audit(nn, xx, as.double(o$catalog$path_coefficient),
      part$low, part$up, part$feasible, as.integer(states), h,
      50L, 100L, as.integer(round(10000 * budget)), 10000L)
    stopifnot(all(z$pass))
    rows[[key]] <- cbind(rule = key, z)
  }
  risks <- do.call(rbind, rows)
  write.csv(risks, file.path(output, "core_risk_checks.csv"), row.names = FALSE)
  write.csv(do.call(rbind, flows), file.path(output, "endpoint_flow_checks.csv"), row.names = FALSE)
  risks
}

verify_clinical_risks <- function(output) {
  rows <- list()
  for (key in names(clinical$policies)) {
    action <- clinical$policies[[key]]
    settings <- clinical$settings[[key]]
    check_actions(clinical$model, action, locked = settings$locked, anchors = FALSE,
                  none = settings$none)
    z <- audit_clinical(clinical$model, clinical$partition, action, clinical$specs[[key]])
    rows[[key]] <- cbind(rule = key, z)
  }
  result <- do.call(rbind, rows)
  write.csv(result, file.path(output, "clinical_risk_checks.csv"), row.names = FALSE)
  result
}

verify_continuation <- function(output) {
  reference <- read_reference_csv("continuation_cases.csv")
  eliminated <- function(n, x) {
    if (n < 3L) return(FALSE)
    p <- gmp::as.bigq("7/20")
    tail <- sum(gmp::as.bigq(choose(n + 1L, 0:x)) * p^(0:x) * (1 - p)^(n + 1L - 0:x))
    tail > gmp::as.bigq("19/20")
  }
  forward <- function(p, n, x) {
    state <- list()
    state[[as.character(x)]] <- gmp::as.bigq(1)
    patients <- events <- gmp::as.bigq(0)
    if (n == 18L) return(c(0, 0))
    for (current in seq(n, 15L, by = 3L)) {
      next_state <- list()
      for (key in names(state)) {
        dlt <- as.integer(key)
        if (eliminated(current, dlt)) next
        mass <- state[[key]]
        patients <- patients + 3 * mass
        events <- events + 3 * p * mass
        for (y in 0:3) {
          target <- as.character(dlt + y)
          value <- mass * choose(3, y) * p^y * (1 - p)^(3 - y)
          if (is.null(next_state[[target]])) next_state[[target]] <- value else
            next_state[[target]] <- next_state[[target]] + value
        }
      }
      state <- next_state
    }
    as.numeric(c(patients, events))
  }
  difference <- numeric(nrow(reference))
  for (i in seq_len(nrow(reference))) {
    row <- reference[i, ]
    p <- decimal_rational(row$p)
    value <- forward(p, row$n, row$x)
    difference[i] <- max(abs(value - c(row$patients, row$DLT)))
  }
  stopifnot(max(difference) < 1e-11)
  result <- data.frame(cases = nrow(reference), maximum_difference = max(difference), pass = TRUE)
  write.csv(result, file.path(output, "continuation_checks.csv"), row.names = FALSE)
  result
}
