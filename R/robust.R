
family_masses <- function(states, directed = FALSE) {
  info <- secondary$robust$families
  groups <- match(paste(o$labels, secondary$robust$support_family, sep = "/"),
                  paste(info$label, info$family, sep = "/"))
  weights <- gmp::as.bigq(o$weights$central)
  totals <- lapply(seq_len(nrow(info)), function(f) sum(weights[groups == f & !is.na(groups)]))
  stopifnot(all(vapply(totals, function(x) x > 0, FALSE)))
  value_low <- value_up <- lapply(seq_len(nrow(info)), function(f) gmp::as.bigq(rep(0, length(states))))
  bounds <- point_mass_bounds(nn[states, , drop = FALSE], xx[states, , drop = FALSE],
    as.double(o$catalog$path_coefficient[states]), o$planning)
  for (c in which(o$weights$central > 0)) {
    f <- groups[c]
    w <- weights[c] / totals[[f]]
    value_low[[f]] <- value_low[[f]] + w * gmp::as.bigq(bounds$low[, c])
    value_up[[f]] <- value_up[[f]] + w * gmp::as.bigq(bounds$up[, c])
  }
  list(low = value_low, up = value_up, info = info, groups = groups, weights = weights, totals = totals)
}

verify_family_gains <- function(output) {
  rules <- c(core$policies[c("central_original", "high_original", "central_sparse80",
    "central_sparse90", "high_sparse80", "high_sparse90")], secondary$policies[c("high_family", "joint_family")])
  states <- sort(unique(unlist(lapply(rules, function(a) which(a != b)))))
  masses <- family_masses(states, directed = TRUE)
  rows <- family_rows <- list()
  for (key in names(rules)) {
    action <- rules[[key]][states]
    low <- up <- gmp::as.bigq(rep(0, nrow(masses$info)))
    for (f in seq_len(nrow(masses$info))) {
      j <- masses$info$label[f]
      h <- as.integer(action == j) - as.integer(b[states] == j)
      low[f] <- 100 * sum(h * ifelse_bigq(h >= 0, masses$low[[f]], masses$up[[f]]))
      up[f] <- 100 * sum(h * ifelse_bigq(h >= 0, masses$up[[f]], masses$low[[f]]))
    }
    for (prior in names(core$priors)) {
      q <- exact_location_prior(prior)
      worst_low <- worst_up <- nominal_low <- gmp::as.bigq(0)
      nominal_up <- gmp::as.bigq(0)
      for (j in seq_len(J)) {
        ids <- which(masses$info$label == j)
        worst_low <- worst_low + q[j] * min(low[ids])
        worst_up <- worst_up + q[j] * min(up[ids])
        den <- sum(do.call(c, masses$totals[ids]))
        share <- do.call(c, masses$totals[ids]) / den
        nominal_low <- nominal_low + q[j] * sum(share * low[ids])
        nominal_up <- nominal_up + q[j] * sum(share * up[ids])
      }
      rows[[length(rows) + 1L]] <- data.frame(rule = key, prior = prior,
        nominal_lower_pp = enclose_rational(nominal_low)$lower,
        nominal_upper_pp = enclose_rational(nominal_up)$upper,
        worst_lower_pp = enclose_rational(worst_low)$lower,
        worst_upper_pp = enclose_rational(worst_up)$upper)
    }
    family_rows[[key]] <- data.frame(rule = key, masses$info,
      lower_pp = enclose_rational(low)$lower, upper_pp = enclose_rational(up)$upper)
  }
  result <- do.call(rbind, rows)
  for (key in c("high_family", "joint_family")) {
    target <- secondary$robust$construction$worst_family_pp[
      secondary$robust$construction$objective == if (key == "high_family") "high" else "joint"]
    r <- result[result$rule == key, ]
    lower <- if (key == "high_family") r$worst_lower_pp[r$prior == "high"] else min(r$worst_lower_pp)
    upper <- if (key == "high_family") r$worst_upper_pp[r$prior == "high"] else min(r$worst_upper_pp)
    stopifnot(lower <= target + 1e-10, upper >= target - 1e-10, upper - lower < 1e-6)
  }
  write.csv(result, file.path(output, "family_gain_bounds.csv"), row.names = FALSE)
  write.csv(do.call(rbind, family_rows), file.path(output, "conditional_family_gain_bounds.csv"), row.names = FALSE)
  result
}

ifelse_bigq <- function(condition, yes, no) {
  value <- no
  value[condition] <- yes[condition]
  value
}

verify_robust_ceilings <- function(output) {
  inputs <- readRDS(file.path(root, "data/gain_ceilings.rds"))$core
  cand <- inputs$candidates
  fixed <- b %in% c(1L, J) | seq_len(K) %in% anchors
  expected <- which(o$eligible & b > 0 & !fixed & outer(b, seq_len(J), "!="), arr.ind = TRUE)
  key <- function(state, action) paste(state, action, sep = "/")
  stopifnot(!anyDuplicated(key(cand$state, cand$action)),
    setequal(key(cand$state, cand$action), key(expected[, 1], expected[, 2])),
    all(cand$reference == b[cand$state]))
  states <- sort(unique(cand$state))
  info <- secondary$robust$families
  groups <- match(paste(o$labels, secondary$robust$support_family, sep = "/"),
    paste(info$label, info$family, sep = "/"))
  base_weights <- gmp::as.bigq(o$weights$central)
  totals <- lapply(seq_len(nrow(info)), function(f) sum(base_weights[which(groups == f)]))
  point_bounds <- point_mass_bounds(nn[states, , drop = FALSE], xx[states, , drop = FALSE],
    as.double(o$catalog$path_coefficient[states]), o$planning)
  risk_indices <- sort(unique(unlist(lapply(secondary$robust$prices, function(p) p$point + 1L))))
  risk_bounds <- point_mass_bounds(nn[states, , drop = FALSE], xx[states, , drop = FALSE],
    as.double(o$catalog$path_coefficient[states]), inputs$grid[risk_indices, , drop = FALSE])
  rows <- membership <- list()
  for (name in c("high", "joint")) {
    prior <- secondary$robust$priors[[name]]
    pi <- gmp::as.bigq(prior$location_prices)
    stopifnot(all(pi >= 0), sum(pi) > 0)
    pi <- pi / sum(pi)
    q <- gmp::as.bigq(rep(0, J))
    for (k in seq_along(pi)) q <- q + pi[k] * exact_location_prior(prior$locations[k])
    shares <- gmp::as.bigq(prior$family_prices)
    stopifnot(all(shares >= 0))
    for (j in seq_len(J)) {
      ids <- which(info$label == j)
      stopifnot(sum(shares[ids]) > 0)
      shares[ids] <- shares[ids] / sum(shares[ids])
    }
    weights <- gmp::as.bigq(rep(0, nrow(o$planning)))
    for (c in which(o$weights$central > 0)) {
      f <- groups[c]
      weights[c] <- q[o$labels[c]] * shares[f] * base_weights[c] / totals[[f]]
    }
    stopifnot(sum(weights) == 1, sum(q) == 1, all(weights >= 0))
    lo <- hi <- lapply(seq_len(J), function(j) gmp::as.bigq(rep(0, length(states))))
    for (c in which(as.numeric(weights) > 0)) {
      j <- o$labels[c]
      lo[[j]] <- lo[[j]] + weights[c] * gmp::as.bigq(point_bounds$low[, c])
      hi[[j]] <- hi[[j]] + weights[c] * gmp::as.bigq(point_bounds$up[, c])
    }
    at <- match(cand$state, states)
    reduced <- gmp::as.bigq(rep(0, nrow(cand)))
    for (j in seq_len(J)) reduced <- reduced + 100 * (
      as.integer(cand$action == j) * hi[[j]][at] - as.integer(cand$reference == j) * lo[[j]][at])
    prices <- secondary$robust$prices[[name]]
    constant <- gmp::as.bigq(0)
    for (k in seq_len(nrow(prices))) {
      price <- gmp::as.bigq(prices$price[k])
      stopifnot(price >= 0)
      tag <- prices$metric[k]
      j <- as.integer(sub("^(C|Tp|Tm)", "", tag))
      p <- gmp::as.bigq(inputs$grid[prices$point[k] + 1L, ])
      stopifnot(all(p >= 0 & p <= 1), all(p[-1] >= p[-J]))
      if (startsWith(tag, "C")) {
        stopifnot(j > 1L, j < J, p[j - 1L] + p[j] <= gmp::as.bigq("1/2"),
          p[j] + p[j + 1L] >= gmp::as.bigq("1/2"))
        h <- as.integer(cand$reference == j) - as.integer(cand$action == j)
      } else h <- (as.integer(cand$action > j) - as.integer(cand$reference > j)) *
        if (startsWith(tag, "Tp")) 1L else -1L
      column <- match(prices$point[k] + 1L, risk_indices)
      bound <- gmp::as.bigq(risk_bounds$low[at, column])
      bound[h < 0] <- gmp::as.bigq(risk_bounds$up[at[h < 0], column])
      reduced <- reduced - 50 * price * h * bound
      constant <- constant + price
    }
    upper <- constant + sum(do.call("c", lapply(split(seq_len(nrow(cand)), cand$state),
      function(ids) max(c(gmp::as.bigq(0), reduced[ids])))))
    value <- enclose_rational(upper)$upper
    saved <- as.numeric(secondary$robust$upper_reference$directed_upper_pp[
      secondary$robust$upper_reference$objective == name])
    stopifnot(abs(value - saved) < 1e-5)
    rows[[name]] <- data.frame(objective = name, upper_pp = value,
      original_upper_pp = saved, prices = nrow(prices), prior_membership_pass = TRUE)
    membership[[name]] <- list(locations = prior$locations, location_mixture = as.character(pi),
      q = as.character(q), family_shares = as.character(shares), total_weight = as.character(sum(weights)))
    cat("Robust gain ceiling checked:", name, "\n")
    flush.console()
  }
  result <- do.call(rbind, rows)
  write.csv(result, file.path(output, "robust_gain_ceilings.csv"), row.names = FALSE)
  jsonlite::write_json(membership, file.path(output, "robust_prior_membership.json"), pretty = TRUE, auto_unbox = TRUE)
  result
}
