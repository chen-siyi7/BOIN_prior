
decimal_rational <- function(value) {
  text <- as.character(value)
  if (any(!grepl("^[+-]?[0-9]+([.][0-9]*)?([eE][+-]?[0-9]+)?$", text)))
    stop("An exact decimal value is malformed.")
  result <- lapply(text, function(item) {
    parts <- strsplit(tolower(item), "e", fixed = TRUE)[[1]]
    exponent <- if (length(parts) == 2L) as.integer(parts[2]) else 0L
    mantissa <- parts[1]
    decimal <- regexpr(".", mantissa, fixed = TRUE)[1]
    places <- if (decimal < 0) 0L else nchar(mantissa) - decimal
    digits <- gsub(".", "", mantissa, fixed = TRUE)
    negative <- startsWith(digits, "-")
    digits <- sub("^[+-]", "", digits)
    digits <- sub("^0+", "", digits)
    if (!nzchar(digits)) digits <- "0"
    numerator <- gmp::as.bigz(digits)
    if (negative) numerator <- -numerator
    power <- places - exponent
    if (power >= 0) gmp::as.bigq(numerator, gmp::as.bigz(10)^power) else
      gmp::as.bigq(numerator * gmp::as.bigz(10)^(-power))
  })
  do.call(c, result)
}

exact_location_prior <- function(name) {
  values <- switch(name,
    central = c("2/100", "4/100", "42/100", "42/100", "8/100", "2/100"),
    high = c("1/100", "2/100", "4/100", "8/100", "20/100", "65/100"),
    equal = rep("1/6", 6L), stop("Unknown location prior."))
  gmp::as.bigq(values)
}

enclose_rational <- function(value) {
  lower <- upper <- as.numeric(value)
  repeat {
    at <- which(gmp::as.bigq(lower) > value)
    if (!length(at)) break
    lower[at] <- adjacent_double(lower[at], FALSE)
  }
  repeat {
    at <- which(gmp::as.bigq(upper) < value)
    if (!length(at)) break
    upper[at] <- adjacent_double(upper[at], TRUE)
  }
  list(lower = lower, upper = upper)
}

exact_conditional_masses <- function(states) {
  p <- o$planning
  weights <- gmp::as.bigq(o$weights$central)
  totals <- lapply(seq_len(J), function(j) sum(weights[o$labels == j]))
  n <- nn[states, , drop = FALSE]
  x <- xx[states, , drop = FALSE]
  coefficient <- gmp::as.bigq(o$catalog$path_coefficient[states])
  value <- lapply(seq_len(J), function(j) gmp::as.bigq(rep(0, length(states))))
  for (c in which(o$weights$central > 0)) {
    mass <- coefficient
    for (j in seq_len(J)) {
      prob <- gmp::as.bigq(p[c, j])
      mass <- mass * prob^x[, j] * (1 - prob)^(n[, j] - x[, j])
    }
    label <- o$labels[c]
    value[[label]] <- value[[label]] + mass * weights[c] / totals[[label]]
  }
  value
}

verify_exact_gains <- function(output) {
  policies <- c(core$policies, core$budget_policies)
  states <- sort(unique(unlist(lapply(policies, function(a) which(a != b)))))
  conditional <- exact_conditional_masses(states)
  contributions <- function(a, prior) {
    q <- exact_location_prior(prior)
    gain <- gmp::as.bigq(rep(0, length(states)))
    for (j in seq_len(J)) {
      sign <- as.integer(a[states] == j) - as.integer(b[states] == j)
      gain <- gain + q[j] * conditional[[j]] * sign
    }
    100 * gain
  }
  gains <- exact <- list()
  for (key in names(policies)) for (prior in c("central", "high", "equal")) {
    individual <- contributions(policies[[key]], prior)
    total <- sum(individual)
    interval <- enclose_rational(total)
    exact[[paste(key, prior, sep = "/")]] <- total
    gains[[length(gains) + 1L]] <- data.frame(
      rule = key, prior = prior, lower_pp = sprintf("%.17g", interval$lower),
      upper_pp = sprintf("%.17g", interval$upper))
  }
  minimum <- retention <- list()
  for (prior in c("central", "high")) for (fraction in c(.8, .9)) {
    bank <- paste0(prior, "_original")
    key <- paste0(prior, "_sparse", round(100 * fraction))
    stopifnot(all(policies[[key]] == b | policies[[key]] == policies[[bank]]))
    target <- exact[[paste(bank, prior, sep = "/")]] * decimal_rational(fraction)
    individual <- contributions(policies[[bank]], prior)
    bounds <- enclose_rational(individual)$upper
    bounds <- sort(pmax(bounds, 0), decreasing = TRUE)
    count <- sum(policies[[key]] != b)
    previous <- sum(gmp::as.bigq(head(bounds, count - 1L)))
    attained <- exact[[paste(key, prior, sep = "/")]]
    stopifnot(attained >= target, previous < target)
    minimum[[length(minimum) + 1L]] <- data.frame(
      rule = key, changes = count, one_fewer_gain_upper_pp = enclose_rational(previous)$upper,
      required_gain_lower_pp = enclose_rational(target)$lower, minimum_count_proved = TRUE)
  }
  for (key in names(core$budget_policies)[grepl("_short", names(core$budget_policies))]) {
    prior <- strsplit(key, "_", fixed = TRUE)[[1]][1]
    bank <- sub("_short", "_bank", key, fixed = TRUE)
    stopifnot(all(policies[[key]] == b | policies[[key]] == policies[[bank]]))
    stopifnot(exact[[paste(key, prior, sep = "/")]] >=
                gmp::as.bigq("9/10") * exact[[paste(bank, prior, sep = "/")]])
    retention[[length(retention) + 1L]] <- data.frame(rule = key, retained_90_percent = TRUE)
  }
  table <- do.call(rbind, gains)
  write.csv(table, file.path(output, "exact_gain_bounds.csv"), row.names = FALSE)
  write.csv(do.call(rbind, minimum), file.path(output, "minimum_table_checks.csv"), row.names = FALSE)
  write.csv(do.call(rbind, retention), file.path(output, "budget_retention_checks.csv"), row.names = FALSE)
  list(gains = table, minimum = minimum, retention = retention)
}
