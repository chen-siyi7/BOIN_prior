supporting_planning <- function(j, phi) {
  values <- list()
  labels <- integer()
  for (dose in seq_len(j)) for (slope in c(.3, .5, .8, 1.2, 1.8, 2.6)) {
    offsets <- if (dose == 1L) -c(.04, .08, .15, .22) else
      if (dose == j) c(.04, .08, .15, .22) else c(-.22, -.08, .08, .22)
    for (offset in offsets) {
      values[[length(values) + 1L]] <- plogis(qlogis(phi) + slope * (seq_len(j) - dose + offset))
      labels <- c(labels, dose)
    }
  }
  list(p = do.call(rbind, values), label = labels)
}

generate_supporting_curves <- function(j, phi) {
  probabilities <- metadata <- list()
  for (dose in seq_len(j)) for (family in c("step", "plateau", "irregular", "probit")) for (k in 1:24) {
    spacing <- switch(family, step = runif(1, .06, .25), plateau = runif(1, .006, .035),
                      irregular = runif(1, .025, .18), probit = runif(1, .04, .35))
    p <- numeric(j)
    p[dose] <- phi
    if (family == "probit") p <- pnorm(qnorm(phi) + spacing * (seq_len(j) - dose)) else {
      if (dose > 1L) for (i in (dose - 1L):1L)
        p[i] <- p[i + 1L] * exp(-spacing * if (family == "irregular") runif(1, .3, 2.5) else 1 / phi)
      if (dose < j) for (i in (dose + 1L):j)
        p[i] <- 1 - (1 - p[i - 1L]) * exp(-spacing * if (family == "irregular") runif(1, .3, 2.5) else 1 / (1 - phi))
      if (family == "step") {
        if (dose > 1L) p[1:(dose - 1L)] <- p[1:(dose - 1L)] * .65
        if (dose < j) p[(dose + 1L):j] <- 1 - (1 - p[(dose + 1L):j]) * .65
      }
    }
    distance <- min(abs(p[-dose] - phi))
    shift <- sample(c(-1, 1), 1) * runif(1, 0, .4) * distance
    if (dose == 1L) shift <- -abs(shift)
    if (dose == j) shift <- abs(shift)
    p <- pmin(pmax(p + shift, 1e-6), 1 - 1e-6)
    probabilities[[length(probabilities) + 1L]] <- p
    metadata[[length(metadata) + 1L]] <- data.frame(mtd = dose, family = family, profile = k)
  }
  for (k in 1:4) {
    p <- if (k <= 2) seq(phi + .03, if (k == 1) .75 else .95, length.out = j) else
      seq(.005, phi - if (k == 3) .03 else .12, length.out = j)
    probabilities[[length(probabilities) + 1L]] <- p
    metadata[[length(metadata) + 1L]] <- data.frame(mtd = 0L, family = "nonbracketing", profile = k)
  }
  p <- do.call(rbind, probabilities)
  colnames(p) <- paste0("p", seq_len(j))
  cbind(do.call(rbind, metadata), p)
}

verify_secondary_generators <- function() {
  for (i in 1:3) {
    rebuilt <- generate_matched_curves(301030 + i)
    saved <- subset(secondary$comparison_curves, library == paste0("rep", i))
    stopifnot(max(abs(rebuilt$P - as.matrix(saved[paste0("p", 1:6)]))) < 1e-14,
      identical(rebuilt$meta$family, saved$family), all(rebuilt$meta$mtd == saved$mtd))
  }
  set.seed(2026110215)
  rebuilt <- generate_family_evaluation_curves(6L, .25)
  saved <- secondary$robust$evaluation$curves
  stopifnot(identical(rebuilt, saved))
  rebuilt <- generate_family_evaluation_curves(3L, .35)
  saved <- readRDS(file.path(root, "data/secondary/clinical_evaluation.rds"))$curves
  stopifnot(identical(rebuilt, saved))
  set.seed(301027)
  for (kind in c("primary", "clinical", "N36", "J8", "cohort2")) {
    j <- switch(kind, clinical = 3L, J8 = 8L, cohort2 = 4L, 6L)
    phi <- if (kind == "clinical") .35 else .25
    rebuilt <- generate_supporting_curves(j, phi)
    if (kind %in% names(secondary$supporting)) {
      item <- secondary$supporting[[kind]]
      saved <- item$curves
      stopifnot(max(abs(as.matrix(rebuilt[paste0("p", seq_len(j))]) -
        as.matrix(saved[, 4:(3 + j)]))) < 1e-14, identical(rebuilt$family, saved$family),
        all(rebuilt$mtd == saved$mtd))
      planning <- supporting_planning(j, phi)
      stopifnot(identical(unname(planning$p), unname(item$model$planning)),
        all(planning$label == item$model$labels))
    }
  }
  invisible(TRUE)
}
