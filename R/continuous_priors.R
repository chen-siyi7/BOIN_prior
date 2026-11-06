
continuous_prior_check <- function(output) {
  model <- clinical$model
  catalog <- model$catalog
  enrollment <- as.matrix(catalog[paste0("n", 1:3)])
  events <- as.matrix(catalog[paste0("x", 1:3)])
  rules <- clinical$policies
  preferences <- rbind(clinical$preferences, data.frame(lambda = 1, underdose = .25, none = .5))
  rows <- normalization <- list()
  integrate_ordered <- function(a, b, nodes, weights) {
    k <- nrow(a)
    density <- function(j) matrix(dbeta(rep(nodes, each = k),
      shape1 = rep(a[, j], length(nodes)), shape2 = rep(b[, j], length(nodes))), nrow = k)
    cdf <- function(j, value, lower = TRUE) matrix(pbeta(rep(value, each = k),
      shape1 = rep(a[, j], length(value)), shape2 = rep(b[, j], length(value)),
      lower.tail = lower), nrow = k)
    middle <- density(2)
    first <- cdf(1, nodes)
    last <- cdf(3, nodes, FALSE)
    total <- as.vector((middle * first * last) %*% weights)
    components <- array(0, c(k, 3L, 3L))
    ranges <- rbind(c(.3, .4), c(.4, 1), c(0, .3))
    for (component in 1:3) {
      low <- ranges[component, 1]
      up <- ranges[component, 2]
      first_segment <- pmax(cdf(1, pmin(nodes, up)) - cdf(1, pmin(nodes, low)), 0)
      last_segment <- pmax(cdf(3, rep(up, length(nodes))) - cdf(3, pmax(nodes, low)), 0)
      components[, 1, component] <- as.vector((middle * first_segment * last) %*% weights)
      components[, 2, component] <- as.vector((middle * first * last) %*%
        (weights * (nodes >= low & nodes <= up)))
      components[, 3, component] <- as.vector((middle * first * last_segment) %*% weights)
    }
    list(total = total, components = components)
  }
  for (order in c(64L, 128L, 256L)) {
    index <- seq_len(order - 1L)
    off_diagonal <- index / sqrt(4 * index^2 - 1)
    jacobi <- matrix(0, order, order)
    jacobi[cbind(index, index + 1L)] <- off_diagonal
    jacobi[cbind(index + 1L, index)] <- off_diagonal
    decomposition <- eigen(jacobi, symmetric = TRUE)
    sorted <- order(decomposition$values)
    nodes <- decomposition$values[sorted]
    weights <- 2 * decomposition$vectors[1, sorted]^2
    segments <- rbind(c(0, .3), c(.3, .4), c(.4, 1))
    shifted_nodes <- unlist(lapply(1:3, function(i)
      mean(segments[i, ]) + diff(segments[i, ]) / 2 * nodes))
    shifted_weights <- unlist(lapply(1:3, function(i) diff(segments[i, ]) / 2 * weights))
    for (prior in names(clinical$priors)) {
      power <- if (prior == "none") 0 else as.numeric(sub(".*_", "", prior))
      x3 <- if (startsWith(prior, "MN")) 5L else if (prior == "none") 0L else 3L
      historic_events <- c(0, 1, x3)
      historic_enrollment <- c(0, 6, 12)
      prior_a <- 1 + power * historic_events
      prior_b <- 1 + power * (historic_enrollment - historic_events)
      constant <- integrate_ordered(matrix(prior_a, 1), matrix(prior_b, 1),
                                    shifted_nodes, shifted_weights)$total
      a <- sweep(1 + events, 2, power * historic_events, "+")
      b <- sweep(1 + enrollment - events, 2, power * (historic_enrollment - historic_events), "+")
      integrated <- integrate_ordered(a, b, shifted_nodes, shifted_weights)
      scale <- exp(catalog$log_path_coefficient + rowSums(lbeta(a, b)) -
                     sum(lbeta(prior_a, prior_b))) / constant
      probability <- integrated$total * scale
      components <- integrated$components * scale
      stopifnot(abs(sum(probability) - 1) < 1e-6)
      outcome <- function(action) {
        value <- matrix(0, length(action), 4L)
        ids <- which(action > 0)
        for (component in 1:3) value[ids, component] <-
          components[cbind(ids, action[ids], component)]
        value[action < 0, 4] <- probability[action < 0]
        value
      }
      baseline <- outcome(model$reference)
      for (rule in names(rules)) {
        delta <- colSums(outcome(rules[[rule]]) - baseline)
        gain <- 100 * (delta[1] - preferences$lambda * delta[2] -
                        preferences$underdose * delta[3] - preferences$none * delta[4])
        rows[[length(rows) + 1L]] <- data.frame(order = order, prior = prior, rule = rule,
          preferences, gain_pp = gain)
      }
      normalization[[length(normalization) + 1L]] <- data.frame(
        order = order, prior = prior, total_mass = sum(probability))
    }
  }
  result <- do.call(rbind, rows)
  compare_numeric_table(result, read_reference_csv("continuous_prior_gains.csv"),
    c("order", "prior", "rule", "lambda", "underdose", "none"), 1e-7)
  write.csv(result, file.path(output, "continuous_prior_gains.csv"), row.names = FALSE)
  write.csv(do.call(rbind, normalization), file.path(output, "continuous_normalization.csv"), row.names = FALSE)
  invisible(result)
}
