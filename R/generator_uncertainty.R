
generator_intervals <- function(output) {
  evaluation <- readRDS(file.path(root, "data/evaluation.rds"))
  curves <- evaluation$curves
  ids <- which(curves$category == "bracketing")
  strata <- split(seq_along(ids), interaction(curves$mtd[ids], curves$family[ids], drop = TRUE))
  stopifnot(length(strata) == 24L, all(lengths(strata) == 40L))
  keys <- c("central_original", "central_sparse80", "central_sparse90",
            "high_original", "high_sparse80", "high_sparse90")
  result <- do.call(rbind, lapply(keys, function(key) {
    prior <- if (startsWith(key, "central")) "central" else "high"
    difference <- evaluation$metrics[[key]][ids, "PCS"] - evaluation$metrics$BOIN[ids, "PCS"]
    means <- vapply(strata, function(i) mean(difference[i]), 0.)
    variances <- vapply(strata, function(i) var(difference[i]) / length(i), 0.)
    weights <- vapply(strata, function(i) core$priors[[prior]][curves$mtd[ids[i[1]]]] / 4, 0.)
    stopifnot(abs(sum(weights) - 1) < 1e-12)
    value <- 100 * sum(weights * means)
    se <- 100 * sqrt(sum(weights^2 * variances))
    data.frame(rule = key, mean_pp = value, SE_pp = se,
               lower_pp = value - 1.96 * se, upper_pp = value + 1.96 * se)
  }))
  write.csv(result, file.path(output, "primary_generator_intervals.csv"), row.names = FALSE)
  result
}
