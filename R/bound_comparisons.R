
bound_comparisons <- function (output) 
{
    ref <- o
    policies <- core$policies[c("central_original", "high_original")]
    states <- sort(unique(unlist(lapply(policies, function(a) which(a != b)))))
    nn <- as.matrix(ref$catalog[paste0("n", 1:J)])
    xx <- as.matrix(ref$catalog[paste0("x", 1:J)])
    t0 <- proc.time()[["elapsed"]]
    global <- get_global(ref$catalog[states, ], J)
    tg <- proc.time()[["elapsed"]] - t0
    lo <- t(apply(part$low, 1, cummax))
    up <- t(apply(part$up, 1, function(x) rev(cummin(rev(x)))))
    stopifnot(all(lo <= up))
    t0 <- proc.time()[["elapsed"]]
    CU <- CL <- matrix(1, length(states), nrow(lo))
    for (k in seq_along(states)) for (j in 1:J) {
        n <- nn[states[k], j]
        x <- xx[states[k], j]
        if (!n) 
            next
        factor <- function(p) {
            v <- numeric(length(p))
            if (x > 0) 
                v <- v + x * log(p)
            if (n > x) 
                v <- v + (n - x) * log1p(-p)
            v
        }
        l <- factor(lo[, j])
        u <- factor(up[, j])
        at <- pmin(pmax(x/n, lo[, j]), up[, j])
        CU[k, ] <- CU[k, ] * exp(factor(at))
        CL[k, ] <- CL[k, ] * exp(pmin(l, u))
    }
    CU <- CU * ref$catalog$path_coefficient[states]
    CL <- CL * ref$catalog$path_coefficient[states]
    tc <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    z <- box_dp_extrema(nn, xx, ref$catalog$log_path_coefficient, part$low, part$up, as.integer(states))
    to <- proc.time()[["elapsed"]] - t0
    metrics <- c(2:(J - 1), (3 * J + 1):(5 * J - 2))
    rows <- list()
    stopifnot(all(CU + 1e-12 >= z$U), all(CL <= z$L + 1e-12))
    for (key in names(policies)) {
        H <- policy_contrast(policies[[key]][states], b[states], J)
        bounds <- list(global_adverse = matrix(colSums(pmax(H, 0) * global), 1), ordered_adverse = t(crossprod(pmax(H, 
            0), z$U)), coordinate_signed = t(crossprod(pmax(H, 0), CU) - crossprod(pmax(-H, 0), CL)), 
            ordered_signed = t(crossprod(pmax(H, 0), z$U) - crossprod(pmax(-H, 0), z$L)))
        for (method in names(bounds)) {
            v <- bounds[[method]]
            largest <- max(vapply(metrics, function(m) {
                ids <- if (nrow(v) == 1) 1 else if (m <= J) which(part$feasible[, m]) else 1:nrow(v)
                max(v[ids, m])
            }, 0))
            rows[[length(rows) + 1]] <- data.frame(policy = key, method = method, largest_upper_pp = 100 * 
                largest, shared_bound_seconds = switch(method, global_adverse = tg, ordered_adverse = to, 
                coordinate_signed = tc, ordered_signed = to), states = length(states), cells = nrow(part$low))
        }
    }
    result <- do.call(rbind, rows)
    names(result)[names(result) == "policy"] <- "rule"
    stopifnot(all(result$largest_upper_pp[result$method == "ordered_signed"] <= 2 + 1e-10), all(result$largest_upper_pp[result$method != 
        "ordered_signed"] > 2))
    write.csv(result, file.path(output, "bound_comparisons.csv"), row.names = FALSE)
    result
}
