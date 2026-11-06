policy_contrast <- function(a, b, J) {
    h <- matrix(0L, length(a), 5L * J - 2L)
    for (j in 1:J) {
        h[, j] <- (b == j) - (a == j)
        h[, J + j] <- (a > j) - (b > j)
        h[, 2 * J + j] <- (a > 0L & a < j) - (b > 0L & b < j)
    }
    for (j in 1:(J - 1L)) {
        h[, 3 * J + j] <- (a > j) - (b > j)
        h[, 4 * J - 1L + j] <- -h[, 3 * J + j]
    }
    h
}

policy_bounds <- function(a, b, U, L, feasible) {
    J <- ncol(feasible)
    changed <- which(a != b)
    h <- policy_contrast(a[changed], b[changed], J)
    if (!length(changed)) 
        return(matrix(0, ncol(U), 5L * J - 2L))
    t(crossprod(pmax(h, 0), U[changed, , drop = FALSE]) - crossprod(pmax(-h, 0), L[changed, , drop = FALSE]))
}
