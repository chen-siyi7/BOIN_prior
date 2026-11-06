
generate_matched_curves <- function(seed) {
    J <- 6L
    phi <- 0.25
    set.seed(seed)
    pl <- meta <- list()
    for (j in 1:J) for (family in c("step", "plateau", "irregular", "probit")) for (k in 1:24) {
        spacing <- switch(family, step = runif(1, 0.06, 0.25), plateau = runif(1, 0.006, 0.035), irregular = runif(1, 
            0.025, 0.18), probit = runif(1, 0.04, 0.35))
        p <- numeric(J)
        p[j] <- phi
        if (family == "probit") 
            p <- pnorm(qnorm(phi) + spacing * ((1:J) - j))
        else {
            if (j > 1) 
                for (i in (j - 1):1) p[i] <- p[i + 1] * exp(-spacing * if (family == "irregular") runif(1, 
                  0.3, 2.5) else 1/phi)
            if (j < J) 
                for (i in (j + 1):J) p[i] <- 1 - (1 - p[i - 1]) * exp(-spacing * if (family == "irregular") runif(1, 
                  0.3, 2.5) else 1/(1 - phi))
            if (family == "step") {
                if (j > 1) 
                  p[1:(j - 1)] <- p[1:(j - 1)] * 0.65
                if (j < J) 
                  p[(j + 1):J] <- 1 - (1 - p[(j + 1):J]) * 0.65
            }
        }
        distance <- min(abs(p[-j] - phi))
        shift <- sample(c(-1, 1), 1) * runif(1, 0, 0.4) * distance
        if (j == 1) 
            shift <- -abs(shift)
        if (j == J) 
            shift <- abs(shift)
        p <- pmin(pmax(p + shift, 1e-06), 1 - 1e-06)
        stopifnot(which.min(abs(p - phi)) == j, p[1] <= phi, p[J] >= phi, all(diff(p) >= 0))
        pl[[length(pl) + 1]] <- p
        meta[[length(meta) + 1]] <- data.frame(mtd = j, family = family, profile = k)
    }
    list(P = do.call(rbind, pl), meta = do.call(rbind, meta))
}

generate_family_evaluation_curves <- function(J, phi) {
    families <- c("ordered_beta", "gamma_gaps", "logit_mixture", "window_plateau")
    ps <- info <- list()
    for (f in families) for (j in 1:J) for (k in 1:40) {
        if (f == "ordered_beta") {
            repeat {
                p <- sort(rbeta(J, 0.5, 1.5))
                if (p[1] <= phi && p[J] >= phi && which.min(abs(p - phi)) == j) 
                  break
            }
        }
        else {
            if (f == "gamma_gaps") {
                gap <- cumsum(c(0, rgamma(J - 1, 0.5, 1)))
                p <- phi + gap - gap[j]
                if (j > 1) 
                  p[1:(j - 1)] <- phi - (phi - p[1:(j - 1)])/(phi - p[1]) * runif(1, 0.1, phi - 0.001)
                if (j < J) 
                  p[(j + 1):J] <- phi + (p[(j + 1):J] - phi)/(p[J] - phi) * runif(1, 0.1, 1 - phi - 0.001)
            }
            if (f == "logit_mixture") {
                loc <- cumsum(c(0, rgamma(J - 1, 0.7, 1)))
                loc <- loc - loc[j]
                slopes <- exp(runif(2, log(0.1), log(3)))
                mix <- runif(1, 0.1, 0.9)
                shift <- runif(1, -2, 2)
                intercept <- uniroot(function(z) mix * plogis(z + shift) + (1 - mix) * plogis(z - shift) - 
                  phi, c(-20, 20), tol = 1e-12)$root
                p <- mix * plogis(intercept + shift + slopes[1] * loc) + (1 - mix) * plogis(intercept - 
                  shift + slopes[2] * loc)
            }
            if (f == "window_plateau") {
                gap <- cumsum(c(0, 10^runif(J - 1, -4, -0.5)))
                p <- phi + gap - gap[j]
                if (j > 1 && p[1] < 0.001) 
                  p[1:(j - 1)] <- phi - (phi - p[1:(j - 1)])/(phi - p[1]) * (phi - 0.001)
                if (j < J && p[J] > 0.999) 
                  p[(j + 1):J] <- phi + (p[(j + 1):J] - phi)/(p[J] - phi) * (0.999 - phi)
            }
            distance <- min(abs(p[-j] - phi))
            d <- runif(1, -0.35, 0.35) * distance
            if (j == 1) 
                d <- -abs(d)
            if (j == J) 
                d <- abs(d)
            p <- pmin(pmax(p + d, 1e-10), 1 - 1e-10)
        }
        stopifnot(all(diff(p) >= 0), p[1] <= phi, p[J] >= phi, which.min(abs(p - phi)) == j)
        ps[[length(ps) + 1]] <- p
        info[[length(info) + 1]] <- data.frame(category = "bracketing", family = f, mtd = j, profile = k)
    }
    for (category in c("all_high", "all_low")) for (f in families) for (k in 1:(40 * J)) {
        z <- switch(f, ordered_beta = sort(rbeta(J, 0.5, 1.5)), gamma_gaps = {
            g <- cumsum(rgamma(J + 1, 0.5, 1))
            g[1:J]/g[J + 1]
        }, logit_mixture = {
            x <- sort(rnorm(J))
            plogis(runif(1, -2, 2) + exp(runif(1, -2, 1)) * x)
        }, window_plateau = {
            x <- cumsum(c(0, 10^runif(J - 1, -4, -1)))
            pmin(1, runif(1, 0.05, 0.75) + x)
        })
        p <- if (category == "all_high") 
            phi + 0.06 + (1 - phi - 0.06) * z
        else (phi - 0.06) * z
        p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
        stopifnot(all(diff(p) >= 0), if (category == "all_high") 
            all(p > phi + 0.05)
        else all(p < phi - 0.05))
        ps[[length(ps) + 1]] <- p
        info[[length(info) + 1]] <- data.frame(category = category, family = f, mtd = 0, profile = k)
    }
    p <- do.call(rbind, ps)
    colnames(p) <- paste0("p", 1:J)
    cbind(do.call(rbind, info), p)
}
