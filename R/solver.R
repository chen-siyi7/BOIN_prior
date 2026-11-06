
solve_stream_files <- function(path, seconds, graph) {
  if (!requireNamespace("highs", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE))
    stop("Rule fitting requires the R packages highs and Matrix.")
  meta <- jsonlite::read_json(file.path(path, "shape.json"), simplifyVector = TRUE)
  candidates <- read.csv(file.path(path, "candidates.csv"))
  variables <- nrow(candidates)
  read_doubles <- function(name, count) {
    con <- file(file.path(path, name), "rb")
    on.exit(close(con))
    out <- readBin(con, "double", n = count, size = 8L, endian = "little")
    if (length(out) != count) stop("Incomplete optimization input.")
    out
  }
  costs <- t(matrix(read_doubles("active_costs.bin", variables * meta$rows),
                    nrow = variables))
  gains <- read_doubles("gains.bin", variables)
  states <- match(candidates$state, unique(candidates$state))
  choice <- Matrix::sparseMatrix(i = states, j = seq_len(variables), x = 1)
  a <- rbind(choice, Matrix::Matrix(costs, sparse = TRUE))
  upper <- c(rep(1, nrow(choice)), meta$upper)
  flows <- 0L
  flow_upper <- numeric()
  if (isTRUE(meta$transport)) {
    edges <- graph$edges
    flows <- nrow(edges)
    key <- function(endpoint, state, kind) paste(endpoint, state, kind, sep = "/")
    nodes <- unique(c(key(edges$endpoint, edges$bad, "bad"),
                      key(edges$endpoint, edges$good, "good")))
    ri <- ci <- integer()
    values <- numeric()
    add <- function(row, column, value) {
      ri <<- c(ri, row)
      ci <<- c(ci, column)
      values <<- c(values, value)
    }
    capacities <- unlist(graph$capacities)
    for (v in seq_len(variables)) {
      s <- candidates$state[v]
      for (endpoint in c(1L, J)) {
        bad <- match(key(endpoint, s, "bad"), nodes)
        good <- match(key(endpoint, s, "good"), nodes)
        if (!is.na(bad) && candidates$reference[v] == endpoint)
          add(bad, v, capacities[as.character(s)])
        if (!is.na(good) && candidates$action[v] == endpoint)
          add(good, v, -capacities[as.character(s)])
      }
    }
    for (e in seq_len(flows)) {
      z <- edges[e, ]
      add(match(key(z$endpoint, z$bad, "bad"), nodes), variables + e, -1)
      add(match(key(z$endpoint, z$good, "good"), nodes), variables + e, 1)
      flow_upper[e] <- min(capacities[as.character(c(z$bad, z$good))])
    }
    constraints <- Matrix::sparseMatrix(i = ri, j = ci, x = values,
                                       dims = c(length(nodes), variables + flows))
    a <- rbind(cbind(a, Matrix::Matrix(0, nrow(a), flows, sparse = TRUE)), constraints)
    upper <- c(upper, rep(0, length(nodes)))
  }
  if (meta$mode == "sparse") {
    a <- rbind(a, Matrix::Matrix(matrix(c(-gains, rep(0, flows)), nrow = 1), sparse = TRUE))
    upper <- c(upper, -meta$target)
    objective <- c(rep(1, variables), rep(0, flows))
  } else {
    objective <- c(-gains, rep(0, flows))
  }
  started <- proc.time()[["elapsed"]]
  fit <- highs::highs_solve(L = objective, lower = rep(0, variables + flows),
    upper = c(rep(1, variables), flow_upper), A = a,
    lhs = rep(-Inf, nrow(a)), rhs = upper,
    types = c(rep("I", variables), rep("C", flows)),
    control = highs::highs_control(time_limit = seconds, mip_rel_gap = .001,
                                   presolve = "on", log_to_console = FALSE))
  x <- fit$primal_solution
  has_solution <- length(x) == variables + flows && all(is.finite(x)) &&
                  max(abs(x[seq_len(variables)] - round(x[seq_len(variables)]))) < 1e-5 &&
                  all(as.vector(a %*% x) <= upper + 1e-5)
  chosen <- if (has_solution) round(x[seq_len(variables)]) else rep(0L, variables)
  bound <- fit$info$mip_dual_bound
  if (is.null(bound) || !is.finite(bound)) bound <- NULL
  info <- list(status = fit$status, message = fit$status_message,
    has_solution = has_solution, seconds = proc.time()[["elapsed"]] - started,
    choices = variables, flows = flows,
    gain_upper_pp = if (meta$mode == "gain" && !is.null(bound)) -bound else NULL,
    minimum_changes_bound = if (meta$mode == "sparse") bound else NULL)
  writeLines(as.character(chosen), file.path(path, "solution.txt"))
  jsonlite::write_json(info, file.path(path, "solver.json"), auto_unbox = TRUE,
                       pretty = TRUE, digits = 17, null = "null")
  jsonlite::toJSON(info, auto_unbox = TRUE, null = "null")
}

solve_clinical_files <- function(path, seconds) {
  if (!requireNamespace("highs", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE))
    stop("Clinical rule fitting requires highs and Matrix.")
  meta <- jsonlite::read_json(file.path(path, "shape.json"), simplifyVector = TRUE)
  candidates <- read.csv(file.path(path, "candidates.csv"))
  variables <- nrow(candidates)
  read_doubles <- function(name, count) {
    con <- file(file.path(path, name), "rb")
    on.exit(close(con))
    value <- readBin(con, "double", n = count, size = 8L, endian = "little")
    if (length(value) != count) stop("Incomplete optimization input.")
    value
  }
  costs <- t(matrix(read_doubles("active_costs.bin", variables * meta$rows), nrow = variables))
  gains <- t(matrix(read_doubles("gains.bin", variables * meta$objectives), nrow = variables))
  states <- match(candidates$state, unique(candidates$state))
  choice <- Matrix::sparseMatrix(i = states, j = seq_len(variables), x = 1)
  constraints <- rbind(cbind(choice, 0), cbind(Matrix::Matrix(costs, sparse = TRUE), 0),
                      cbind(Matrix::Matrix(-gains, sparse = TRUE), 1))
  rhs <- c(rep(1, nrow(choice)), meta$upper, rep(0, nrow(gains)))
  started <- proc.time()[["elapsed"]]
  fit <- highs::highs_solve(L = c(rep(0, variables), -1),
    lower = c(rep(0, variables), -Inf), upper = c(rep(1, variables), Inf),
    A = constraints, lhs = rep(-Inf, nrow(constraints)), rhs = rhs,
    types = c(rep("I", variables), "C"),
    control = highs::highs_control(time_limit = seconds, mip_rel_gap = .001,
                                   presolve = "on", log_to_console = FALSE))
  x <- fit$primal_solution
  valid <- length(x) == variables + 1L && all(is.finite(x)) &&
    max(abs(x[seq_len(variables)] - round(x[seq_len(variables)]))) < 1e-5 &&
    all(as.vector(constraints %*% x) <= rhs + 1e-5)
  chosen <- if (valid) round(x[seq_len(variables)]) else rep(0L, variables)
  bound <- fit$info$mip_dual_bound
  if (is.null(bound) || !is.finite(bound)) bound <- NULL
  info <- list(status = fit$status, message = fit$status_message, has_solution = valid,
    seconds = proc.time()[["elapsed"]] - started,
    numerical_upper_pp = if (!is.null(bound)) -bound else NULL)
  writeLines(as.character(chosen), file.path(path, "solution.txt"))
  jsonlite::write_json(info, file.path(path, "solver.json"), pretty = TRUE,
                       auto_unbox = TRUE, null = "null", digits = 17)
  jsonlite::toJSON(info, auto_unbox = TRUE, null = "null")
}
