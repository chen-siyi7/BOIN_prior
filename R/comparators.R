
normal_quadrature <- function(n) {
  matrix <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) matrix[i, i + 1L] <- matrix[i + 1L, i] <- sqrt(i)
  eigen <- eigen(matrix, symmetric = TRUE)
  weights <- eigen$vectors[1, ]^2
  keep <- weights > 0
  list(node = eigen$values[keep], weight = weights[keep] / sum(weights[keep]))
}

calibrate_logit_prior <- function() {
  dose <- c(7, 14, 29, 51, 103, 188)
  z <- log(dose / 51)
  phi <- .25
  probabilities <- c(.025, .5, .975)
  b_low <- log(.05) / log1p(-phi)
  a_high <- log(.21) / log(phi)
  med_low <- 1 - .5^(1 / b_low)
  med_high <- .5^(1 / a_high)
  slope <- (qlogis(med_high) - qlogis(med_low)) / (z[6] - z[1])
  intercept <- qlogis(med_low) - slope * z[1]
  med <- plogis(intercept + slope * z)
  a <- ifelse(med >= .5, log(.5) / log(med), 1)
  beta <- ifelse(med < .5, log(.5) / log1p(-med), 1)
  target <- t(vapply(seq_along(dose), function(j) qbeta(probabilities, a[j], beta[j]), numeric(3)))
  rule <- normal_quadrature(121L)
  implied <- function(par) {
    value <- exp(par[2] + exp(par[4]) * rule$node)
    t(vapply(z, function(zz) vapply(probabilities, function(probability) {
      if (zz == 0) return(plogis(qnorm(probability, par[1], exp(par[3]))))
      f <- function(t) sum(rule$weight * pnorm((t - zz * value - par[1]) / exp(par[3]))) - probability
      plogis(uniroot(f, c(-1e6, 1e6), tol = 1e-9)$root)
    }, 0.), numeric(3)))
  }
  loss <- function(par) sum((implied(par) - target)^2)
  fits <- lapply(list(c(intercept, log(slope), log(1), log(.7)),
    c(-1, 0, log(1.5), log(1)), c(-1.5, .3, log(.7), log(.5))), function(p)
    optim(p, loss, method = "L-BFGS-B", lower = c(-5, -4, log(.01), log(.01)),
      upper = c(3, 3, log(4), log(3)), control = list(maxit = 500, factr = 1e7)))
  fit <- fits[[which.min(vapply(fits, `[[`, 0., "value"))]]
  broad <- optim(fit$par[1:3], function(p) loss(c(p, log(.75))), method = "L-BFGS-B",
    lower = c(-5, -4, 0), upper = c(3, 3, log(4)), control = list(maxit = 500, factr = 1e7))
  stopifnot(fit$convergence == 0L, broad$convergence == 0L)
  setting <- function(p) list(d = dose, d_ref = 51, mu0 = p[1], mu1 = p[2], sd0 = exp(p[3]), sd1 = exp(p[4]))
  list(settings = list(quantile = setting(fit$par), weak = setting(c(broad$par, log(.75)))),
    target = target, implied = implied(fit$par), objectives = c(quantile = fit$value, weak = broad$value))
}

logit_posterior <- function(model, setting, nodes) {
  rule <- normal_quadrature(nodes)
  grid <- expand.grid(i = seq_along(rule$node), j = seq_along(rule$node))
  intercept <- setting$mu0 + setting$sd0 * rule$node[grid$i]
  slope <- exp(setting$mu1 + setting$sd1 * rule$node[grid$j])
  weights <- rule$weight[grid$i] * rule$weight[grid$j]
  eta <- outer(slope, log(setting$d / setting$d_ref)) + intercept
  softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  logp <- -softplus(-eta)
  logq <- -softplus(eta)
  p <- plogis(eta)
  j <- model$design$J
  n <- as.matrix(model$catalog[paste0("n", seq_len(j))])
  x <- as.matrix(model$catalog[paste0("x", seq_len(j))])
  out <- matrix(0, nrow(n), j)
  for (first in seq(1, nrow(n), by = 128L)) {
    ids <- first:min(first + 127L, nrow(n))
    value <- x[ids, , drop = FALSE] %*% t(logp) +
      (n[ids, , drop = FALSE] - x[ids, , drop = FALSE]) %*% t(logq)
    value <- sweep(value, 2, log(weights), "+")
    value <- exp(value - apply(value, 1, max))
    out[ids, ] <- (value / rowSums(value)) %*% p
  }
  stopifnot(all(is.finite(out)), all(out >= 0 & out <= 1))
  out
}

logit_actions <- function(posterior, locked = FALSE) {
  loss <- abs(posterior - o$design$phi)
  loss[!o$eligible] <- Inf
  action <- max.col(-loss, ties.method = "first")
  action[b < 0] <- -1L
  if (locked) {
    fixed <- b %in% c(1L, J) | seq_along(b) %in% anchors
    action[fixed] <- b[fixed]
  }
  as.integer(action)
}

verify_logit_comparators <- function(output) {
  calibration <- calibrate_logit_prior()
  for (name in names(calibration$settings)) {
    saved <- secondary$comparators$settings[[name]]
    stopifnot(max(abs(unlist(calibration$settings[[name]]) - unlist(saved))) < 1e-7)
  }
  saveRDS(calibration, file.path(output, "logistic_prior_calibration.rds"))
  rows <- list()
  for (name in c("quantile", "weak")) {
    setting <- secondary$comparators$settings[[name]]
    low <- logit_posterior(o, setting, 181L)
    high <- logit_posterior(o, setting, 241L)
    saved <- secondary$comparators$posterior_means[[paste0(name, "_unlocked")]]
    disagreements <- sum(logit_actions(low) != logit_actions(high))
    stopifnot(disagreements == 0L, max(abs(low - high)) < 1e-5,
      max(abs(high - saved)) < 1e-10,
      identical(logit_actions(high), secondary$comparators$proposals[[paste0(name, "_unlocked")]]))
    proposal <- logit_actions(high, locked = TRUE)
    for (prior in c("central", "high")) {
      key <- paste(prior, "logit", name, sep = "_")
      action <- secondary$policies[[key]]
      stopifnot(all(action == b | action == proposal))
    }
    rows[[name]] <- data.frame(prior = name, low_nodes = 181L, high_nodes = 241L,
      action_disagreements = disagreements, maximum_mean_difference = max(abs(low - high)),
      saved_mean_difference = max(abs(high - saved)))
    cat("Logistic posterior and filtered bank checked:", name, "\n")
    flush.console()
  }
  
  setting <- secondary$comparators$settings$quantile
  calibration <- secondary$comparators$calibration
  probabilities <- c(.025, .5, .975)
  quantiles <- as.matrix(calibration[c("implied_025", "implied_50", "implied_975")])
  difference <- numeric()
  for (j in seq_len(J)) for (k in 1:3) {
    integrand <- function(t) dnorm(t) * pnorm((qlogis(quantiles[j, k]) -
      exp(setting$mu1 + setting$sd1 * t) * log(setting$d[j] / setting$d_ref) - setting$mu0) / setting$sd0)
    value <- integrate(integrand, -12, 12, rel.tol = 1e-10, subdivisions = 1000L)$value
    difference <- c(difference, value - probabilities[k])
  }
  stopifnot(max(abs(difference)) < 1e-5)
  result <- do.call(rbind, rows)
  write.csv(result, file.path(output, "logistic_quadrature_checks.csv"), row.names = FALSE)
  write.csv(data.frame(maximum_CDF_difference = max(abs(difference))),
    file.path(output, "logistic_calibration_check.csv"), row.names = FALSE)
  result
}
