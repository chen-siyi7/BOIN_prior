
boin_boundaries <- function(phi, psaf = 0.6, ptox = 1.4) {
    ps <- psaf * phi
    pt <- ptox * phi
    le <- log((1 - ps)/(1 - phi))/log((phi * (1 - ps))/(ps * (1 - phi)))
    ld <- log((1 - phi)/(1 - pt))/log((pt * (1 - phi))/(phi * (1 - pt)))
    c(le = le, ld = ld)
}

elim_flags <- function(n, x, phi, cutoff = 0.95, tab = NULL) {
    if (is.null(tab)) {
        own <- (n >= 3) & ((1 - pbeta(phi, 1 + x, 1 + n - x)) > cutoff)
    }
    else {
        own <- tab[cbind(n + 1, x + 1)]
    }
    as.logical(cumsum(own) > 0)
}

elim_table <- function(phi, Nmax, cutoff = 0.95) {
    tab <- matrix(FALSE, Nmax + 1, Nmax + 1)
    for (nn in 3:Nmax) for (xx in 0:nn) tab[nn + 1, xx + 1] <- (1 - pbeta(phi, 1 + xx, 1 + nn - xx)) > 
        cutoff
    tab
}

simulate_one_trial <- function(p_true, phi, N, cohort = 3, bnd = NULL, tab = NULL) {
    J <- length(p_true)
    if (is.null(bnd)) 
        bnd <- boin_boundaries(phi)
    n <- integer(J)
    x <- integer(J)
    cur <- 1
    for (c_i in seq_len(N%/%cohort)) {
        d <- rbinom(1, cohort, p_true[cur])
        n[cur] <- n[cur] + cohort
        x[cur] <- x[cur] + d
        el <- elim_flags(n, x, phi, tab = tab)
        if (el[1]) 
            break
        rate <- x[cur]/n[cur]
        nxt <- if (rate <= bnd["le"]) 
            cur + 1
        else if (rate > bnd["ld"]) 
            cur - 1
        else cur
        if (nxt > J) 
            nxt <- J
        if (nxt < 1) 
            nxt <- 1
        while (nxt > 1 && el[nxt]) nxt <- nxt - 1
        cur <- nxt
    }
    list(n = n, x = x)
}

wpava <- function(y, w) {
    lv <- numeric(0)
    lw <- numeric(0)
    ln <- integer(0)
    for (i in seq_along(y)) {
        cv <- y[i]
        cw <- w[i]
        cn <- 1L
        while (length(lv) > 0 && lv[length(lv)] > cv) {
            k <- length(lv)
            cw2 <- lw[k] + cw
            cv <- (lv[k] * lw[k] + cv * cw)/cw2
            cn <- ln[k] + cn
            cw <- cw2
            lv <- lv[-k]
            lw <- lw[-k]
            ln <- ln[-k]
        }
        lv <- c(lv, cv)
        lw <- c(lw, cw)
        ln <- c(ln, cn)
    }
    rep(lv, ln)
}

iso_estimate <- function(n, x, phi, tab = NULL) {
    J <- length(n)
    iso <- rep(NA_real_, J)
    tr <- which(n > 0 & !elim_flags(n, x, phi, tab = tab))
    if (length(tr) == 0) 
        return(iso)
    ph <- (x + 0.05)/(n + 0.1)
    v <- (x + 0.05) * (n - x + 0.05)/((n + 0.1)^2 * (n + 0.1 + 1))
    iso[tr] <- wpava(ph[tr], 1/v[tr])
    iso
}

.pick <- function(idx, score, iso, phi) {
    smax <- max(score)
    tied <- idx[abs(score - smax) < 1e-12]
    if (length(tied) == 1L) 
        return(tied)
    offset <- match(tied, idx) * 1e-10
    tied[which.min(abs(iso[tied] + offset - phi))]
}

select_boin <- function(n, x, phi, iso = NULL, tab = NULL) {
    J <- length(n)
    if (is.null(iso)) 
        iso <- iso_estimate(n, x, phi, tab = tab)
    idx <- which(n > 0 & !elim_flags(n, x, phi, tab = tab))
    if (length(idx) == 0) 
        return(-1L)
    .pick(idx, -abs(iso[idx] - phi), iso, phi)
}
