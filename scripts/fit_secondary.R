script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script)), "common.R"))
root <- script_root()
source(file.path(root, "R/project.R"))
args <- arguments(c("output", "mode", "prior", "design", "proposal", "seconds", "rounds", "start"))
if (is.null(args$output)) stop("Supply --output with a new directory outside the code folder.")
mode <- if (is.null(args$mode)) "robust" else args$mode
prior <- if (is.null(args$prior)) "high" else args$prior
seconds <- if (is.null(args$seconds)) 60 else as.numeric(args$seconds)
rounds <- if (is.null(args$rounds)) 10L else as.integer(args$rounds)
start <- if (is.null(args$start)) "saved" else args$start
stopifnot(mode %in% c("robust", "matched", "supporting"), start %in% c("saved", "boin"),
  is.finite(seconds), seconds > 0, is.finite(rounds), rounds > 0)
output <- boin_output(args$output, root)
project <- boin_project(root, native = TRUE)
if (mode == "robust") {
  stopifnot(prior %in% c("high", "joint"))
  fit <- project$fit_family_rule(prior, output, seconds, rounds, use_saved = start == "saved")
} else {
  stopifnot(prior %in% c("central", "high"))
  if (mode == "supporting") {
    kind <- if (is.null(args$design)) "N36" else args$design
    stopifnot(kind %in% names(project$secondary$supporting))
    item <- project$secondary$supporting[[kind]]
    project$o <- item$model
    project$part <- item$partition
    project$b <- item$model$reference
    project$J <- item$model$design$J
    project$K <- length(project$b)
    project$P <- 5L * project$J - 2L
    project$nn <- as.matrix(project$o$catalog[paste0("n", seq_len(project$J))])
    project$xx <- as.matrix(project$o$catalog[paste0("x", seq_len(project$J))])
    project$lp <- project$o$catalog$log_path_coefficient
    project$anchors <- vapply(seq_len(project$J), function(j) which(project$b == j)[1], 0L)
    key <- paste0("shared_", if (prior == "high") "high_boundary" else "central")
    q <- item$q[[if (prior == "high") "high_boundary" else "central"]][-1]
    bank <- item$policies[[key]]
    fixed <- project$b %in% c(1L, project$J) | seq_len(project$K) %in% project$anchors
    id <- which(project$o$eligible & project$b > 0 & !fixed &
      outer(project$b, seq_len(project$J), "!="), arr.ind = TRUE)
    cand <- data.frame(state = id[, 1], action = id[, 2], reference = project$b[id[, 1]])
  } else {
    proposal_name <- if (is.null(args$proposal)) "quantile" else args$proposal
    stopifnot(proposal_name %in% c("quantile", "weak", "mode", "mean"))
    if (proposal_name %in% c("quantile", "weak")) {
      proposal <- project$logit_actions(project$secondary$comparators$posterior_means[[paste0(proposal_name, "_unlocked")]], locked = TRUE)
      key <- paste(prior, "logit", proposal_name, sep = "_")
    } else {
      key <- paste(proposal_name, prior, sep = "_")
      
      m <- project$exact_mass(project$o$catalog, project$o$planning)
      c <- project$gain_C(project$core$priors[[prior]])
      if (proposal_name == "mode") {
        c[!project$o$eligible] <- -Inf
        proposal <- max.col(c, ties.method = "first")
      } else {
        w <- project$o$weights$central
        for (j in seq_len(project$J)) {
          ids <- which(project$o$labels == j)
          w[ids] <- project$core$priors[[prior]][j] * w[ids] / sum(w[ids])
        }
        wm <- sweep(m, 2, w, "*")
        posterior <- (wm %*% project$o$planning) / pmax(rowSums(wm), .Machine$double.xmin)
        proposal <- project$logit_actions(posterior, locked = TRUE)
      }
      fixed <- project$b %in% c(1L, project$J) | seq_len(project$K) %in% project$anchors | project$b < 0
      proposal[fixed] <- project$b[fixed]
    }
    states <- which(proposal != project$b)
    cand <- data.frame(state = states, action = proposal[states], reference = project$b[states])
    bank <- project$secondary$policies[[key]]
    q <- project$core$priors[[prior]]
    stopifnot(all(bank == project$b | bank == proposal))
  }
  contribution <- project$gain_C(q)
  initial <- if (start == "saved") list(bank) else list(project$b)
  fit <- project$fit_stream(cand, contribution, initial, output, seconds = seconds, rounds = rounds)
}
saveRDS(fit, file.path(output, "fit.rds"))
cat("Verified fresh", mode, "fit with gain", fit$gain_pp, "percentage points.\n")
