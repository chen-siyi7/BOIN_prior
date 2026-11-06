
evaluate_rules <- function(output, rules = core$policies) {
  frozen <- readRDS(file.path(root, "data/evaluation.rds"))
  curves <- frozen$curves
  probabilities <- as.matrix(curves[paste0("p", seq_len(J))])
  actions <- c(list(BOIN = b), rules)
  metric_names <- c("PCS", "A", "R", "B", "N", "DLT", "changed", "below", "above",
                    "corrected", "lost", "other")
  result <- lapply(actions, function(a) {
    matrix(0, nrow(curves), length(metric_names), dimnames = list(NULL, metric_names))
  })
  for (first in seq(1, nrow(curves), by = 256L)) {
    ids <- first:min(first + 255L, nrow(curves))
    p <- probabilities[ids, , drop = FALSE]
    mass <- exact_mass(o$catalog, p)
    labels <- curves$mtd[ids]
    for (key in names(actions)) {
      a <- actions[[key]]
      selected <- vapply(seq_len(J), function(j) colSums(mass[a == j, , drop = FALSE]),
                         numeric(length(ids)))
      value <- cbind(
        PCS = rowSums(selected * (col(selected) == labels)),
        A = rowSums(selected * (abs(p - o$design$phi) <= .05)),
        R = rowSums(selected * (p > o$design$phi + .05)),
        B = rowSums(selected * (p < o$design$phi - .05)),
        N = colSums(mass[a < 0, , drop = FALSE]),
        DLT = rowSums(selected * p),
        changed = colSums(mass[a != b, , drop = FALSE]),
        below = rowSums(selected * (col(selected) < labels)),
        above = rowSums(selected * (col(selected) > labels)))
      corrected <- lost <- other <- numeric(length(ids))
      for (j in seq_len(J)) {
        at <- which(labels == j)
        if (!length(at)) next
        corrected[at] <- colSums(mass[a == j & b != j, at, drop = FALSE])
        lost[at] <- colSums(mass[a != j & b == j, at, drop = FALSE])
        other[at] <- colSums(mass[a != b & a != j & b != j, at, drop = FALSE])
      }
      result[[key]][ids, ] <- cbind(value, corrected, lost, other)
    }
  }
  summary <- paired <- absolute <- positions <- list()
  for (category in unique(curves$category)) {
    priors <- if (category == "bracketing") names(core$priors) else "within_class"
    ids <- which(curves$category == category)
    for (prior in priors) {
      w <- if (prior == "within_class") rep(1 / length(ids), length(ids)) else {
        q <- core$priors[[prior]]
        q[curves$mtd[ids]] / as.vector(table(curves$mtd[ids])[as.character(curves$mtd[ids])])
      }
      for (key in names(actions)) {
        value <- result[[key]][ids, , drop = FALSE]
        baseline <- result$BOIN[ids, , drop = FALSE]
        average <- colSums(value * w)
        base_average <- colSums(baseline * w)
        bracket <- category == "bracketing"
        summary[[length(summary) + 1L]] <- data.frame(
          category = category, prior = prior, rule = key,
          BOIN_PCS = if (bracket) 100 * base_average["PCS"] else NA_real_,
          PCS = if (bracket) 100 * average["PCS"] else NA_real_,
          PCS_gain = if (bracket) 100 * (average["PCS"] - base_average["PCS"]) else NA_real_,
          within = 100 * average["A"], toxic = 100 * average["R"],
          under = 100 * average["B"], none = 100 * average["N"],
          override = 100 * average["changed"],
          mean_DLT = 100 * average["DLT"] / (1 - average["N"]),
          worst_PCS = if (bracket) 100 * min(value[, "PCS"] - baseline[, "PCS"]) else NA_real_)
        if (bracket && prior %in% c("central", "high") &&
            key %in% c("BOIN", paste0(prior, c("_original", "_sparse80", "_sparse90")))) {
          absolute[[length(absolute) + 1L]] <- data.frame(
            prior = prior, rule = key, below = 100 * average["below"],
            MTD = 100 * average["PCS"], above = 100 * average["above"],
            none = 100 * average["N"], within = 100 * average["A"],
            toxic = 100 * average["R"], mean_DLT = 100 * average["DLT"] / (1 - average["N"]))
          if (key != "BOIN") {
            paired[[length(paired) + 1L]] <- data.frame(
              prior = prior, rule = key, changed = 100 * average["changed"],
              corrected = 100 * average["corrected"], lost = 100 * average["lost"],
              other = 100 * average["other"],
              net_PCS_gain = 100 * (average["corrected"] - average["lost"]),
              net_gain_per_changed_trial = 100 * (average["corrected"] - average["lost"]) /
                average["changed"])
          }
        }
      }
    }
  }
  for (key in names(actions)) for (j in seq_len(J)) {
    ids <- which(curves$category == "bracketing" & curves$mtd == j)
    positions[[length(positions) + 1L]] <- data.frame(
      rule = key, label = j, BOIN_PCS = 100 * mean(result$BOIN[ids, "PCS"]),
      PCS = 100 * mean(result[[key]][ids, "PCS"]),
      gain = 100 * mean(result[[key]][ids, "PCS"] - result$BOIN[ids, "PCS"]))
  }
  tables <- list(summary = do.call(rbind, summary), paired_changes = do.call(rbind, paired),
                 absolute_outcomes = do.call(rbind, absolute), positions = do.call(rbind, positions))
  for (name in names(tables)) {
    rownames(tables[[name]]) <- NULL
    write.csv(tables[[name]], file.path(output, paste0(name, ".csv")), row.names = FALSE)
  }
  saveRDS(list(curves = curves, metrics = result), file.path(output, "evaluation.rds"))
  tables
}

compare_numeric_table <- function(actual, expected, keys, tolerance = 1e-8) {
  key <- function(x) do.call(paste, c(x[keys], sep = "/"))
  at <- match(key(actual), key(expected))
  if (anyNA(at) || nrow(actual) != nrow(expected)) stop("Table keys differ from the reference.")
  for (column in intersect(names(actual), names(expected))) {
    if (is.numeric(actual[[column]])) {
      a <- actual[[column]]
      e <- expected[[column]][at]
      if (!identical(is.na(a), is.na(e)) ||
          any(abs(a[is.finite(a)] - e[is.finite(a)]) > tolerance))
        stop("Reference mismatch in ", column)
    }
  }
  invisible(TRUE)
}
