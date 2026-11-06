gain_C <- function(q) {
    M <- exact_mass(o$catalog, o$planning)
    C <- matrix(0, K, J)
    for (j in 1:J) {
        ids <- which(o$labels == j)
        w <- o$weights$central[ids]
        w <- w/sum(w)
        C[, j] <- q[j] * as.vector(M[, ids, drop = FALSE] %*% w)
    }
    C
}

score <- function(a, C) 100 * sum(C[cbind(1:K, pmax(a, 1))] - C[cbind(1:K, pmax(b, 1))])

stream_bounds <- function(a, block = 256L) {
    states <- which(a != b)
    ans <- matrix(0, nrow(part$low), P)
    if (!length(states)) 
        return(ans)
    H <- policy_contrast(a[states], b[states], J)
    for (first in seq(1, nrow(part$low), by = block)) {
        ids <- first:min(first + block - 1L, nrow(part$low))
        z <- box_dp_extrema(nn, xx, lp, part$low[ids, , drop = FALSE], part$up[ids, , drop = FALSE], 
            as.integer(states))
        ans[ids, ] <- t(crossprod(pmax(H, 0), z$U) - crossprod(pmax(-H, 0), z$L))
    }
    ans
}

active_costs <- function(cand, active) {
    states <- sort(unique(cand$state))
    cells <- sort(unique(active$cell))
    z <- box_dp_extrema(nn, xx, lp, part$low[cells, , drop = FALSE], part$up[cells, , drop = FALSE], 
        as.integer(states))
    H <- policy_contrast(cand$action, cand$reference, J)
    si <- match(cand$state, states)
    vapply(1:nrow(active), function(i) {
        m <- active$metric[i]
        ci <- match(active$cell[i], cells)
        h <- H[, m]
        ifelse(h > 0, z$U[si, ci], ifelse(h < 0, -z$L[si, ci], 0))/0.02
    }, numeric(nrow(cand)))
}

relevant <- function(m) if (m <= J) which(part$feasible[, m]) else seq_len(nrow(part$low))

fit_stream <- function(cand, C, initial, out, mode = "gain", transport = FALSE, target = 0, seconds = 30,
    rounds = 6, budget = .02) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    started <- proc.time()[["elapsed"]]
    metrics <- c(2:(J - 1), (3L * J + 1L):P)
    caps <- rep(budget - 1e-06, P)
    if (!transport) 
        metrics <- c(1:J, (3L * J + 1L):P)
    caps[c(1, J)] <- 0
    initial <- c(list(b), initial)
    if (mode == "sparse") 
        initial <- Filter(function(a) score(a, C) >= target - 1e-10, initial)
    stopifnot(length(initial) > 0)
    best <- if (mode == "sparse") 
        initial[[which.min(vapply(initial, function(a) sum(a != b), 0L))]]
    else initial[[which.max(vapply(initial, score, 0, C = C))]]
    upper <- Inf
    lower_changes <- 0
    history <- list()
    active <- data.frame(metric = integer(), cell = integer())
    used <- stream_bounds(best)
    for (m in metrics) {
        ids <- relevant(m)
        active <- rbind(active, data.frame(metric = m, cell = head(ids[order(used[ids, m], decreasing = TRUE)], 
            4)))
    }
    G <- 100 * (C[cbind(cand$state, cand$action)] - C[cbind(cand$state, cand$reference)])
    write.csv(cand, file.path(out, "candidates.csv"), row.names = FALSE)
    for (round in seq_len(rounds)) {
        costs <- active_costs(cand, active)
        for (key in c("active_costs", "gains")) {
            con <- file(file.path(out, paste0(key, ".bin")), "wb")
            writeBin(as.double(if (key == "gains") G else costs), con, size = 8, endian = "little")
            close(con)
        }
        write.csv(active, file.path(out, "active_rows.csv"), row.names = FALSE)
        jsonlite::write_json(list(rows = nrow(active), upper = unname(caps[active$metric]/0.02), mode = mode, 
            transport = transport, target = target), file.path(out, "shape.json"), auto_unbox = TRUE)
        txt <- solve_stream_files(out, seconds, graph)
        if (!is.null(attr(txt, "status"))) 
            stop(paste(txt, collapse = "\n"))
        info <- jsonlite::read_json(file.path(out, "solver.json"), simplifyVector = TRUE)
        if (!is.null(info$gain_upper_pp)) 
            upper <- min(upper, info$gain_upper_pp)
        if (!is.null(info$minimum_changes_bound)) 
            lower_changes <- max(lower_changes, info$minimum_changes_bound)
        chosen <- scan(file.path(out, "solution.txt"), quiet = TRUE) == 1L
        stopifnot(!anyDuplicated(cand$state[chosen]))
        a <- b
        a[cand$state[chosen]] <- cand$action[chosen]
        if (transport) {
            write.csv(data.frame(state = 1:K, action = a), file.path(out, "proposed_actions.csv"), row.names = FALSE)
            txt <- audit_transport_file(out, "proposed_actions.csv", o, graph)
            if (!is.null(attr(txt, "status"))) 
                stop(paste(txt, collapse = "\n"))
        }
        used <- stream_bounds(a)
        new <- data.frame(metric = integer(), cell = integer())
        for (m in metrics) {
            ids <- relevant(m)
            bad <- ids[used[ids, m] > caps[m] + 1e-12]
            if (length(bad)) 
                new <- rbind(new, data.frame(metric = m, cell = head(bad[order(used[bad, m], decreasing = TRUE)], 
                  4)))
        }
        eligible <- !nrow(new) && (mode != "sparse" || score(a, C) >= target - 1e-08)
        better <- if (mode == "sparse") 
            sum(a != b) < sum(best != b)
        else score(a, C) > score(best, C) + 1e-12
        if (eligible && better) 
            best <- a
        history[[round]] <- data.frame(round = round, active = nrow(active), proposed_gain = score(a, 
            C), best_gain = score(best, C), changed = sum(best != b), gain_upper = upper, minimum_changes_bound = lower_changes, 
            risk_pass = !nrow(new), solver_seconds = info$seconds)
        write.csv(do.call(rbind, history), file.path(out, "history.csv"), row.names = FALSE)
        cat(out, "round", round, "gain", score(best, C), "changed", sum(best != b), "upper", upper, "\n")
        flush.console()
        if (!nrow(new)) 
            break
        active <- unique(rbind(active, new))
    }
    write.csv(data.frame(state = 1:K, action = best), file.path(out, "saved_actions.csv"), row.names = FALSE)
    states <- which(best != b)
    H <- policy_contrast(best[states], b[states], J)
    if (transport) 
        H[, c(1, J)] <- 0L
    audit <- fixedpoint_interval_audit(nn, xx, as.double(o$catalog$path_coefficient), part$low, part$up, 
        part$feasible, as.integer(states), H, 50L, 100L,
        as.integer(round(10000 * budget)), 10000L, 0L, 0L)
    if (transport) 
        audit <- audit[!audit$metric %in% c("C1", paste0("C", J)), ]
    stopifnot(all(audit$pass))
    write.csv(audit, file.path(out, "risk_audit.csv"), row.names = FALSE)
    if (transport) {
        txt <- audit_transport_file(out, "saved_actions.csv", o, graph)
        if (!is.null(attr(txt, "status"))) 
            stop(paste(txt, collapse = "\n"))
    }
    unlink(file.path(out, c("active_costs.bin", "gains.bin")))
    list(action = best, gain_pp = score(best, C), changed = sum(best != b), gain_upper_pp = upper, minimum_changes_bound = lower_changes, 
        seconds = proc.time()[["elapsed"]] - started, history = do.call(rbind, history))
}
