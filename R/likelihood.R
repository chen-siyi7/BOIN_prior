key_of <- function(n, x) paste(c(n, x), collapse = ",")

add_to <- function(environment, key, value) {
    if (exists(key, environment, inherits = FALSE))
        environment[[key]] <- environment[[key]] + value
    else environment[[key]] <- value
}

enumerate_catalog <- function(Nmax, J = 6L, phi = 0.25, cohort = 3L, start = 1L, stop_at = Inf) {
    tab <- elim_table(phi, Nmax)
    bnd <- boin_boundaries(phi)
    terminal <- new.env(hash = TRUE, parent = emptyenv())
    frontier <- new.env(hash = TRUE, parent = emptyenv())
    frontier[[paste(c(start, integer(2L * J)), collapse = ",")]] <- 1
    for (step in seq_len(Nmax%/%cohort)) {
        next_frontier <- new.env(hash = TRUE, parent = emptyenv())
        for (key in ls(frontier, all.names = TRUE)) {
            v <- as.integer(strsplit(key, ",", fixed = TRUE)[[1L]])
            current <- v[1L]
            n <- v[1L + seq_len(J)]
            x <- v[1L + J + seq_len(J)]
            for (dlt in 0:cohort) {
                nn <- n
                xx <- x
                nn[current] <- nn[current] + cohort
                xx[current] <- xx[current] + dlt
                a <- frontier[[key]] * choose(cohort, dlt)
                el <- elim_flags(nn, xx, phi, tab = tab)
                if (el[1L] || step == Nmax%/%cohort || max(nn) >= stop_at) {
                  add_to(terminal, key_of(nn, xx), a)
                }
                else {
                  rate <- xx[current]/nn[current]
                  next_dose <- if (rate <= bnd["le"]) 
                    current + 1L
                  else if (rate > bnd["ld"]) 
                    current - 1L
                  else current
                  next_dose <- max(1L, min(J, next_dose))
                  while (next_dose > 1L && el[next_dose]) next_dose <- next_dose - 1L
                  add_to(next_frontier, paste(c(next_dose, nn, xx), collapse = ","), a)
                }
            }
        }
        frontier <- next_frontier
        counts <- c(unlist(as.list(frontier), use.names = FALSE), unlist(as.list(terminal), use.names = FALSE))
        stopifnot(all(counts == floor(counts)), max(counts) < 2^53)
    }
    keys <- ls(terminal, all.names = TRUE)
    mat <- t(vapply(strsplit(keys, ",", fixed = TRUE), as.integer, integer(2L * J)))
    catlg <- as.data.frame(mat)
    names(catlg) <- c(paste0("n", seq_len(J)), paste0("x", seq_len(J)))
    catlg$path_coefficient <- unlist(mget(keys, terminal), use.names = FALSE)
    stopifnot(max(catlg$path_coefficient) < 2^53)
    catlg$log_path_coefficient <- log(catlg$path_coefficient)
    catlg$state_id <- seq_len(nrow(catlg))
    catlg$boin_selection <- vapply(seq_len(nrow(catlg)), function(s) select_boin(mat[s, seq_len(J)], 
        mat[s, J + seq_len(J)], phi, tab = tab), 0L)
    catlg$eligible_mask <- vapply(seq_len(nrow(catlg)), function(s) {
        n <- mat[s, seq_len(J)]
        x <- mat[s, J + seq_len(J)]
        ids <- which(n > 0L & !elim_flags(n, x, phi, tab = tab))
        as.integer(sum(2^(ids - 1L)))
    }, 0L)
    list(catalog = catlg, design = list(J = J, N = Nmax, phi = phi, cohort = cohort, start = start, stop_at = stop_at))
}

exact_mass <- function(catlg, P) {
    J <- ncol(P)
    nn <- as.matrix(catlg[paste0("n", 1:J)])
    xx <- as.matrix(catlg[paste0("x", 1:J)])
    v <- matrix(catlg$log_path_coefficient, nrow(catlg), nrow(P))
    for (j in 1:J) for (s in seq_len(nrow(catlg))) {
        if (xx[s, j] > 0L) 
            v[s, ] <- v[s, ] + xx[s, j] * log(P[, j])
        if (nn[s, j] > xx[s, j]) 
            v[s, ] <- v[s, ] + (nn[s, j] - xx[s, j]) * log1p(-P[, j])
    }
    M <- exp(v)
    stopifnot(max(abs(colSums(M) - 1)) < 1e-09)
    M
}

get_global <- function(catlg, J) {
    nn <- as.matrix(catlg[paste0("n", 1:J)])
    xx <- as.matrix(catlg[paste0("x", 1:J)])
    vapply(seq_len(nrow(catlg)), function(s) {
        ids <- which(nn[s, ] > 0)
        p <- wpava(xx[s, ids]/nn[s, ids], nn[s, ids])
        v <- 0
        for (k in seq_along(ids)) {
            i <- ids[k]
            if (xx[s, i] > 0) 
                v <- v + xx[s, i] * log(p[k])
            if (xx[s, i] < nn[s, i]) 
                v <- v + (nn[s, i] - xx[s, i]) * log1p(-p[k])
        }
        exp(catlg$log_path_coefficient[s] + v)
    }, 0)
}

label_mask <- function(low, up, phi) {
    J <- ncol(low)
    sapply(1:J, function(j) {
        l <- low[, j]
        u <- up[, j]
        if (j > 1L) 
            u <- pmin(u, 2 * phi - low[, j - 1L])
        if (j < J) 
            l <- pmax(l, 2 * phi - up[, j + 1L])
        l <= u + 1e-14
    })
}

make_partition <- function(J, phi, bins_count = 8L) {
    edges <- if (bins_count == 8L) 
        c(0, phi * (1:5)/5, 1.28 * phi, 0.5, 1)
    else c(0, phi/2, phi, 0.5, 0.75, 1)
    recurse <- function(prefix = integer(), smallest = 1L) {
        if (length(prefix) == J) 
            return(matrix(prefix, nrow = 1L))
        do.call(rbind, lapply(smallest:(length(edges) - 1L), function(i) recurse(c(prefix, i), i)))
    }
    bins <- recurse()
    nr <- nrow(bins)
    low <- matrix(edges[bins], nrow = nr)
    up <- matrix(edges[bins + 1L], nrow = nr)
    list(low = low, up = up, feasible = label_mask(low, up, phi), bins = bins)
}
