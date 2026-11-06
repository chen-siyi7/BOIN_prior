contrasts <- function(a, b, meta) {
    h <- matrix(0L, length(a), nrow(meta))
    for (m in 1:nrow(meta)) {
        z <- meta[m, ]
        l <- z$left
        v <- z$right
        j <- z$index
        aw <- as.integer(a > l & a <= v)
        bw <- as.integer(b > l & b <= v)
        h[, m] <- if (z$tag == "N") 
            as.integer(a < 0) - as.integer(b < 0)
        else if (startsWith(z$tag, "C")) 
            as.integer(b == j) - as.integer(a == j)
        else if (startsWith(z$tag, "Tp")) 
            as.integer(a > j) - as.integer(b > j)
        else if (startsWith(z$tag, "Tm")) 
            as.integer(b > j) - as.integer(a > j)
        else if (startsWith(z$tag, "A")) 
            bw - aw
        else if (startsWith(z$tag, "R")) 
            as.integer(a > v) - as.integer(b > v)
        else bw - aw + as.integer(a > v) - as.integer(b > v)
    }
    h
}

all_bounds <- function(a, b, U, L, spec) {
    ids <- which(a != b)
    if (!length(ids)) 
        return(matrix(0, ncol(U), nrow(spec$meta)))
    h <- contrasts(a[ids], b[ids], spec$meta)
    t(crossprod(pmax(h, 0), U[ids, , drop = FALSE]) - crossprod(pmax(-h, 0), L[ids, , drop = FALSE]))
}

audit_clinical <- function(o, part, a, spec) {
    ids <- which(a != o$reference)
    h <- contrasts(a[ids], o$reference[ids], spec$meta)
    z <- bracket_interval_audit(as.matrix(o$catalog[paste0("n", 1:o$design$J)]), as.matrix(o$catalog[paste0("x", 
        1:o$design$J)]), as.double(o$catalog$path_coefficient), part$low, part$up, spec$mask, as.integer(ids), 
        h, as.integer(spec$meta$kind), as.integer(spec$meta$index), as.integer(spec$meta$left), as.integer(spec$meta$right), 
        as.integer(round(200 * o$design$phi)), as.integer(round(100 * (o$design$phi - 0.05))), as.integer(round(100 * 
            (o$design$phi + 0.05))), 100L, as.integer(round(10000 * spec$meta$cap)), as.character(spec$meta$tag), 
        10000L)
    stopifnot(all(z$pass))
    z
}
