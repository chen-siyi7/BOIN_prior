
clinical_cost_fit <- function(o, part, U, L, Cs, spec, locked, none, initial, out, seconds = 10, rounds = 8) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    b <- o$reference
    K <- length(b)
    J <- ncol(o$eligible)
    cols <- function(a) ifelse(a > 0, a, J + 1L)
    E <- cbind(o$eligible, rep(none, K))
    aa <- c(1:J, -1L)
    id <- which(E & b > 0 & !locked & outer(b, aa, "!="), arr.ind = TRUE)
    cand <- data.frame(state = id[, 1], action = aa[id[, 2]], reference = b[id[, 1]], column = id[, 2])
    stopifnot(nrow(cand) > 0)
    stopifnot(all(vapply(Cs, ncol, 0L) == J + 1L))
    scorevec <- function(a) vapply(Cs, function(C) 100 * sum(C[cbind(1:K, cols(a))] - C[cbind(1:K, cols(b))]), 
        0)
    score <- function(a) min(scorevec(a))
    meta <- spec$meta
    P <- nrow(meta)
    cap <- meta$cap - ifelse(meta$cap > 0, 1e-06, 0)
    feasible <- function(a) {
        if (any(a[locked] != b[locked]) || any(a[b < 0] != -1) || (!none && any(a[b > 0] < 0))) 
            return(FALSE)
        good <- which(a > 0)
        if (!all(o$eligible[cbind(good, a[good])])) 
            return(FALSE)
        u <- all_bounds(a, b, U, L, spec)
        all(u[spec$mask] <= matrix(cap, nrow(u), P, byrow = TRUE)[spec$mask] + 1e-12)
    }
    initial <- Filter(feasible, c(list(BOIN = b), initial))
    best <- initial[[which.max(vapply(initial, score, 0))]]
    bestscore <- score(best)
    upper <- Inf
    history <- list()
    H <- contrasts(cand$action, cand$reference, meta)
    G <- vapply(Cs, function(C) 100 * (C[cbind(cand$state, cand$column)] - C[cbind(cand$state, cols(cand$reference))]), 
        numeric(nrow(cand)))
    used <- all_bounds(best, b, U, L, spec)
    active <- data.frame(metric = integer(), cell = integer())
    for (m in 1:P) {
        ids <- which(spec$mask[, m])
        active <- rbind(active, data.frame(metric = m, cell = head(ids[order(used[ids, m], decreasing = TRUE)], 
            4)))
    }
    write.csv(cand, file.path(out, "candidates.csv"), row.names = FALSE)
    started <- proc.time()[["elapsed"]]
    for (round in 1:rounds) {
        costs <- vapply(1:nrow(active), function(i) {
            m <- active$metric[i]
            r <- active$cell[i]
            h <- H[, m]
            ifelse(h > 0, h * U[cand$state, r], ifelse(h < 0, h * L[cand$state, r], 0))/0.02
        }, numeric(nrow(cand)))
        jsonlite::write_json(list(rows = nrow(active), objectives = ncol(G), upper = as.numeric(cap[active$metric]/0.02)), 
            file.path(out, "shape.json"), auto_unbox = TRUE)
        for (key in c("active_costs", "gains")) {
            con <- file(file.path(out, paste0(key, ".bin")), "wb")
            writeBin(as.double(if (key == "gains") G else costs), con, size = 8, endian = "little")
            close(con)
        }
        txt <- solve_clinical_files(out, seconds)
        if (!is.null(attr(txt, "status"))) 
            stop(paste(txt, collapse = "\n"))
        si <- jsonlite::read_json(file.path(out, "solver.json"), simplifyVector = TRUE)
        if (!is.null(si$numerical_upper_pp)) 
            upper <- min(upper, si$numerical_upper_pp)
        chosen <- scan(file.path(out, "solution.txt"), quiet = TRUE) == 1L
        stopifnot(!anyDuplicated(cand$state[chosen]))
        a <- b
        a[cand$state[chosen]] <- cand$action[chosen]
        u <- all_bounds(a, b, U, L, spec)
        new <- data.frame(metric = integer(), cell = integer())
        for (m in 1:P) {
            ids <- which(spec$mask[, m] & u[, m] > cap[m] + 1e-12)
            if (length(ids)) 
                new <- rbind(new, data.frame(metric = m, cell = head(ids[order(u[ids, m], decreasing = TRUE)], 
                  4)))
        }
        repaired <- a
        for (step in 0:K) {
            excess <- sweep(u, 2, cap, "-")
            excess[!spec$mask] <- -Inf
            if (max(excess) <= 1e-12) 
                break
            at <- which.max(excess)
            r <- (at - 1)%%nrow(u) + 1
            m <- (at - 1)%/%nrow(u) + 1
            states <- which(repaired != b)
            hh <- contrasts(repaired[states], b[states], meta)
            pos <- which(hh[, m] > 0 & U[states, r] > 0)
            stopifnot(length(pos) > 0)
            ss <- states[pos]
            now <- scorevec(repaired)
            loss <- vapply(ss, function(s) min(now) - min(now - vapply(Cs, function(C) 100 * (C[s, cols(repaired[s])] - 
                C[s, cols(b[s])]), 0)), 0)
            k <- pos[which.min(loss/pmin(hh[pos, m] * U[ss, r], max(excess)))]
            state <- states[k]
            h <- hh[k, ]
            for (mm in which(h != 0)) u[, mm] <- u[, mm] - if (h[mm] > 0) 
                h[mm] * U[state, ]
            else h[mm] * L[state, ]
            repaired[state] <- b[state]
        }
        stopifnot(feasible(repaired))
        if (score(repaired) > bestscore + 1e-12) {
            best <- repaired
            bestscore <- score(best)
        }
        history[[round]] <- data.frame(round = round, score = bestscore, upper = upper, violations = nrow(new), 
            seconds = si$seconds)
        write.csv(do.call(rbind, history), file.path(out, "history.csv"), row.names = FALSE)
        cat(out, "round", round, "gain", bestscore, "upper", upper, "\n")
        flush.console()
        if (upper - bestscore < 0.005 || (!nrow(new) && round >= 2)) 
            break
        before <- nrow(active)
        active <- unique(rbind(active, new))
        if (before == nrow(active) && round >= 3) 
            break
    }
    stopifnot(feasible(best))
    unlink(file.path(out, c("active_costs.bin", "gains.bin")))
    list(action = best, gain = bestscore, gains = scorevec(best), upper = upper, changed = sum(best != 
        b), seconds = proc.time()[["elapsed"]] - started, spec = spec, history = do.call(rbind, history))
}
