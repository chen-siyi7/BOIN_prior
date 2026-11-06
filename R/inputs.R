verify_partition <- function(partition, j, phi, bins = 8L) {
  initial <- make_partition(j, phi, bins)
  edges <- sort(unique(c(initial$low, initial$up)))
  key <- function(values) apply(values, 1, paste, collapse = "/")
  root_keys <- key(matrix(match(initial$low, edges), ncol = j))
  midpoints <- (partition$low + partition$up) / 2
  leaf_keys <- key(matrix(findInterval(midpoints, edges), ncol = j))
  groups <- split(seq_len(nrow(partition$low)), leaf_keys)
  stopifnot(setequal(names(groups), root_keys), !anyDuplicated(root_keys))
  cover <- function(low, up, ids) {
    if (!length(ids)) {
      stopifnot(!all(cummax(low) <= up))
      return(invisible(TRUE))
    }
    if (length(ids) == 1L && all(partition$low[ids, ] == low) &&
        all(partition$up[ids, ] == up)) return(invisible(TRUE))
    axis <- which.max(up - low)
    middle <- (low[axis] + up[axis]) / 2
    left <- ids[partition$up[ids, axis] <= middle]
    right <- ids[partition$low[ids, axis] >= middle]
    stopifnot(middle > low[axis], middle < up[axis],
              length(left) + length(right) == length(ids), !length(intersect(left, right)))
    split_up <- up
    split_up[axis] <- middle
    split_low <- low
    split_low[axis] <- middle
    cover(low, split_up, left)
    cover(split_low, up, right)
  }
  for (i in seq_along(root_keys)) {
    cover(initial$low[i, ], initial$up[i, ], groups[[root_keys[i]]])
  }
  invisible(TRUE)
}

verify_inputs <- function(output) {
  supporting <- readRDS(file.path(root, "data/supporting.rds"))
  models <- c(list(core = core$model, clinical = clinical$model), supporting)
  rows <- list()
  for (name in names(models)) {
    model <- models[[name]]
    d <- model$design
    started <- proc.time()[["elapsed"]]
    rebuilt <- enumerate_catalog(d$N, d$J, d$phi, d$cohort, d$start, d$stop_at)
    columns <- c(paste0("n", seq_len(d$J)), paste0("x", seq_len(d$J)),
                 "path_coefficient", "boin_selection", "eligible_mask")
    stopifnot(identical(rebuilt$catalog[columns], model$catalog[columns]))
    rows[[name]] <- data.frame(design = name, states = nrow(rebuilt$catalog),
      seconds = proc.time()[["elapsed"]] - started, pass = TRUE)
  }
  for (name in c("core", "clinical")) {
    object <- get(name)
    verify_partition(object$partition, object$model$design$J, object$model$design$phi)
  }
  curves <- generate_evaluation_curves(6L, .25)
  frozen <- readRDS(file.path(root, "data/evaluation.rds"))$curves
  stopifnot(identical(curves, frozen))
  set.seed(301026)
  planning <- t(apply(matrix(runif(10000L * 3L), ncol = 3L), 1, sort))
  stopifnot(identical(unname(planning), unname(clinical$model$planning)))
  priors <- list(none = rep(1 / nrow(planning), nrow(planning)))
  for (source in c("MN", "DD")) for (power in c(.25, .5, 1)) {
    x3 <- if (source == "MN") 5L else 3L
    log_likelihood <- power * (log(planning[, 2]) + 5 * log1p(-planning[, 2]) +
                                x3 * log(planning[, 3]) + (12 - x3) * log1p(-planning[, 3]))
    weights <- exp(log_likelihood - max(log_likelihood))
    priors[[paste(source, power, sep = "_")]] <- weights / sum(weights)
  }
  for (name in names(priors)) stopifnot(max(abs(priors[[name]] - clinical$priors[[name]])) < 1e-14)
  table <- do.call(rbind, rows)
  write.csv(table, file.path(output, "catalog_checks.csv"), row.names = FALSE)
  write.csv(data.frame(core_cells = nrow(core$partition$low),
    clinical_cells = nrow(clinical$partition$low), coverage = TRUE,
    evaluation_curves = nrow(curves), clinical_planning_curves = nrow(planning),
    clinical_priors = length(priors)), file.path(output, "input_checks.csv"), row.names = FALSE)
  table
}
