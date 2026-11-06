
run_cse <- function(output, replay = TRUE, seconds = 60, cell_limit = 200000L) {
  for (name in c("cse", "cse_resource"))
    Rcpp::sourceCpp(file.path(root, "src", paste0(name, ".cpp")), env = environment(), showOutput = FALSE)
  policies <- core$policies[c("central_original", "high_original")]
  changed <- sort(unique(unlist(lapply(policies, function(a) which(a != b)))))
  h <- do.call(cbind, lapply(policies, function(a) policy_contrast(a[changed], b[changed], J)))
  keep <- c(2:(J - 1L), (3L * J + 1L):P)
  h <- h[, c(keep, keep + P), drop = FALSE]
  orders <- c(1.01, 1.025, 1.05, 1.1, 1.2, 1.4, 1.7, 2, 2.5, 3, 4, 6, 8,
              12, 16, 24, 32, 48, 64, 96, 128, 256, 512, 1024)
  masks <- function(low, up) {
    possible <- rowSums(t(apply(low, 1, cummax)) > up) == 0
    mask <- cbind(label_mask(low, up, .25)[, 2:(J - 1L), drop = FALSE],
                  matrix(TRUE, nrow(low), 2L * (J - 1L)))
    mask[!possible, ] <- FALSE
    cbind(mask, mask)
  }
  score <- function(bounds, mask) {
    bounds[!mask] <- -Inf
    apply(bounds, 1, max)
  }
  calculate <- function(low, up, method) {
    if (method == "Ordered likelihood") {
      z <- box_dp_extrema(nn, xx, lp, low, up, as.integer(changed))
      return(t(crossprod(pmax(h, 0), z$U) - crossprod(pmax(-h, 0), z$L)))
    }
    ordered_low <- t(apply(low, 1, cummax))
    ordered_up <- t(apply(up, 1, function(x) rev(cummin(rev(x)))))
    center <- pmin(pmax((ordered_low + ordered_up) / 2, 1e-12), 1 - 1e-12)
    mass <- exp(nn[changed, , drop = FALSE] %*% t(log1p(-center)) +
      xx[changed, , drop = FALSE] %*% t(log(center) - log1p(-center)) + lp[changed])
    plus <- t(crossprod(pmax(h, 0), mass))
    minus <- t(crossprod(pmax(-h, 0), mass))
    if (method == "CSE independent") cse_bounds(low, up, center, apply(nn, 2, max), orders, plus, minus)$upper else
      cse_resource_bounds(low, up, center, unique(nn), orders, plus, minus)$upper
  }
  original <- read.csv(file.path(root, "data/secondary/original_cse_summary.csv"))
  initial <- make_partition(J, .25)
  rows <- list()
  for (method in c("Ordered likelihood", "CSE independent", "CSE stopping")) {
    key <- switch(method, "Ordered likelihood" = "signed", "CSE independent" = "cse", "CSE stopping" = "sequential")
    target <- original$final_cells[original$method == method]
    stopifnot(length(target) == 1L)
    low <- initial$low
    up <- initial$up
    started <- proc.time()[["elapsed"]]
    bound <- calculate(low, up, method)
    history <- list()
    repeat {
      mask <- masks(low, up)
      values <- score(bound, mask)
      bad <- which(values > .02)
      history[[length(history) + 1L]] <- data.frame(cells = nrow(low), unresolved = length(bad),
        upper_pp = 100 * max(values), seconds = proc.time()[["elapsed"]] - started)
      if (!length(bad) || (replay && nrow(low) >= target) ||
          (!replay && proc.time()[["elapsed"]] - started >= seconds) || nrow(low) >= cell_limit) break
      room <- if (replay) target - nrow(low) else cell_limit - nrow(low)
      ids <- head(bad[order(values[bad], decreasing = TRUE)], min(1000L, room))
      lo <- rbind(low[ids, , drop = FALSE], low[ids, , drop = FALSE])
      hi <- rbind(up[ids, , drop = FALSE], up[ids, , drop = FALSE])
      for (k in seq_along(ids)) {
        axis <- which.max(up[ids[k], ] - low[ids[k], ])
        middle <- (up[ids[k], axis] + low[ids[k], axis]) / 2
        hi[k, axis] <- middle
        lo[k + length(ids), axis] <- middle
      }
      child <- pmin(calculate(lo, hi, method), bound[rep(ids, 2L), , drop = FALSE])
      stay <- setdiff(seq_len(nrow(low)), ids)
      low <- rbind(low[stay, , drop = FALSE], lo)
      up <- rbind(up[stay, , drop = FALSE], hi)
      bound <- rbind(bound[stay, , drop = FALSE], child)
    }
    independent_seconds <- NA_real_
    if (method == "Ordered likelihood" && !length(bad)) {
      started_audit <- proc.time()[["elapsed"]]
      partition <- list(low = low, up = up, feasible = label_mask(low, up, .25))
      verify_partition(partition, J, .25)
      audit <- secondary_risk(o, partition, policies)
      independent_seconds <- proc.time()[["elapsed"]] - started_audit
      write.csv(audit, file.path(output, "adaptive_independent_risk_checks.csv"), row.names = FALSE)
    }
    row <- data.frame(method = method, initial_cells = nrow(initial$low), final_cells = nrow(low),
      evaluated_cells = nrow(initial$low) + 2 * (nrow(low) - nrow(initial$low)),
      upper_pp = 100 * max(values), seconds = proc.time()[["elapsed"]] - started,
      independent_seconds = independent_seconds, target_reached = !length(bad),
      scope = if (replay) "Replay to original refinement count" else "Fresh resource-limited benchmark")
    if (replay) {
      ref <- original[original$method == method, ]
      stopifnot(row$final_cells == ref$final_cells, abs(row$upper_pp - ref$upper_pp) < 1e-7)
    }
    rows[[method]] <- row
    write.csv(do.call(rbind, history), file.path(output, paste0("adaptive_", key, "_history.csv")), row.names = FALSE)
    write.csv(do.call(rbind, rows), file.path(output, "adaptive_comparison.csv"), row.names = FALSE)
    cat("Adaptive comparison complete:", method, "\n")
    flush.console()
  }
  
  check <- list()
  p0 <- c(.06, .12, .20, .32, .45, .60)
  mass0 <- exact_mass(o$catalog, matrix(p0, 1L))[, 1]
  for (k in 1:4) {
    p <- c(.04, .10, .22, .35, .50, .70) * (1 + (k - 1) / 100)
    mass <- exact_mass(o$catalog, matrix(p, 1L))[, 1]
    moment <- c(sum(mass^2 / mass0), sum(mass0^2 / mass))
    upper <- exp(cse_resource_moments(matrix(p, 1L), matrix(p, 1L), matrix(p0, 1L), unique(nn), 2))
    stopifnot(all(moment <= upper * (1 + 1e-12)))
    check[[k]] <- data.frame(case = k, forward = moment[1], forward_upper = upper[1],
      reverse = moment[2], reverse_upper = upper[2], pass = TRUE)
  }
  write.csv(do.call(rbind, check), file.path(output, "cse_moment_checks.csv"), row.names = FALSE)
  do.call(rbind, rows)
}
