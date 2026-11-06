
fit_family_rule <- function(prior, output, seconds = 60, rounds = 10L, use_saved = TRUE) {
  stopifnot(prior %in% c("high", "joint"), seconds > 0, rounds > 0)
  if (!requireNamespace("highs", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE))
    stop("Fitting requires highs and Matrix.")
  info <- secondary$robust$families
  mass <- exact_mass(o$catalog, o$planning)
  contributions <- lapply(seq_len(nrow(info)), function(f) {
    ids <- which(o$labels == info$label[f] & secondary$robust$support_family == info$family[f] & o$weights$central > 0)
    weights <- o$weights$central[ids]
    value <- matrix(0, K, J)
    value[, info$label[f]] <- as.vector(mass[, ids, drop = FALSE] %*% (weights / sum(weights)))
    value
  })
  rm(mass)
  q <- if (prior == "high") matrix(core$priors$high, 1L) else do.call(rbind, core$priors)
  family_values <- function(a) vapply(contributions, function(c) score(a, c), 0.)
  value <- function(a) {
    gains <- family_values(a)
    minimum <- vapply(seq_len(J), function(j) min(gains[info$label == j]), 0.)
    min(q %*% minimum)
  }
  initial <- if (use_saved) c(list(b), core$policies[c("central_original", "high_original")],
    secondary$policies[c("high_family", "joint_family")]) else list(b)
  best <- initial[[which.max(vapply(initial, value, 0.))]]
  locked <- b %in% c(1L, J) | seq_len(K) %in% anchors
  id <- which(o$eligible & b > 0 & !locked & outer(b, seq_len(J), "!="), arr.ind = TRUE)
  cand <- data.frame(state = id[, 1], action = id[, 2], reference = b[id[, 1]])
  gains <- vapply(contributions, function(c) 100 * (c[cbind(cand$state, cand$action)] -
    c[cbind(cand$state, cand$reference)]), numeric(nrow(cand)))
  metrics <- c(seq_len(J), (3L * J + 1L):P)
  caps <- rep(.02 - 1e-6, P)
  caps[c(1L, J)] <- 0
  used <- stream_bounds(best)
  active <- do.call(rbind, lapply(metrics, function(m) {
    ids <- relevant(m)
    data.frame(metric = m, cell = head(ids[order(used[ids, m], decreasing = TRUE)], 4L))
  }))
  variables <- nrow(cand)
  states <- match(cand$state, unique(cand$state))
  choice <- Matrix::sparseMatrix(i = states, j = seq_len(variables), x = 1)
  auxiliary <- Matrix::sparseMatrix(i = seq_len(nrow(info)), j = info$label, x = 1,
    dims = c(nrow(info), J))
  history <- list()
  upper <- Inf
  for (round in seq_len(rounds)) {
    costs <- t(active_costs(cand, active))
    matrix <- rbind(
      cbind(choice, Matrix::Matrix(0, nrow(choice), J + 1L, sparse = TRUE)),
      cbind(Matrix::Matrix(costs, sparse = TRUE), Matrix::Matrix(0, nrow(costs), J + 1L, sparse = TRUE)),
      cbind(Matrix::Matrix(-t(gains), sparse = TRUE), auxiliary, 0),
      cbind(Matrix::Matrix(0, nrow(q), variables, sparse = TRUE), Matrix::Matrix(-q, sparse = TRUE), 1))
    rhs <- c(rep(1, nrow(choice)), caps[active$metric] / .02,
             rep(0, nrow(info) + nrow(q)))
    started <- proc.time()[["elapsed"]]
    fit <- highs::highs_solve(L = c(rep(0, variables + J), -1),
      lower = c(rep(0, variables), rep(-100, J), 0),
      upper = c(rep(1, variables), rep(100, J + 1L)), A = matrix,
      lhs = rep(-Inf, nrow(matrix)), rhs = rhs,
      types = c(rep("I", variables), rep("C", J + 1L)),
      control = highs::highs_control(time_limit = seconds, mip_rel_gap = .001,
        presolve = "on", log_to_console = FALSE))
    x <- fit$primal_solution
    valid <- length(x) == variables + J + 1L && all(is.finite(x)) &&
      max(abs(x[seq_len(variables)] - round(x[seq_len(variables)]))) < 1e-5 &&
      all(as.vector(matrix %*% x) <= rhs + 1e-5)
    proposed <- b
    if (valid) {
      chosen <- which(round(x[seq_len(variables)]) == 1L)
      stopifnot(!anyDuplicated(cand$state[chosen]))
      proposed[cand$state[chosen]] <- cand$action[chosen]
    }
    bound <- fit$info$mip_dual_bound
    if (length(bound) && is.finite(bound)) upper <- min(upper, -bound)
    used <- stream_bounds(proposed)
    new <- do.call(rbind, lapply(metrics, function(m) {
      ids <- relevant(m)
      bad <- ids[used[ids, m] > caps[m] + 1e-12]
      data.frame(metric = rep(m, min(4L, length(bad))),
        cell = head(bad[order(used[bad, m], decreasing = TRUE)], 4L))
    }))
    if (!nrow(new) && value(proposed) > value(best)) best <- proposed
    history[[round]] <- data.frame(round = round, active = nrow(active), proposed_pp = value(proposed),
      best_pp = value(best), numerical_upper_pp = upper, risk_pass = !nrow(new),
      solver_seconds = proc.time()[["elapsed"]] - started, status = fit$status)
    write.csv(do.call(rbind, history), file.path(output, "history.csv"), row.names = FALSE)
    if (!nrow(new)) break
    active <- unique(rbind(active, new))
  }
  check <- secondary_risk(o, part, list(fitted = best))
  write.csv(check, file.path(output, "risk_checks.csv"), row.names = FALSE)
  write.csv(cand, file.path(output, "candidates.csv"), row.names = FALSE)
  write.csv(active, file.path(output, "active_rows.csv"), row.names = FALSE)
  write.csv(data.frame(state = seq_len(K), action = best), file.path(output, "saved_actions.csv"), row.names = FALSE)
  list(action = best, gain_pp = value(best), numerical_upper_pp = upper,
    history = do.call(rbind, history), prior = prior, initial = if (use_saved) "Saved feasible rules" else "BOIN")
}
