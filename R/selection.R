
select_dose <- function(counts, dlts, table = "central_sparse90") {
  choices <- c(names(core$policies), names(clinical$policies))
  if (!table %in% choices) stop("Unknown decision table.")
  model <- if (table %in% names(clinical$policies)) clinical$model else core$model
  action <- if (table %in% names(clinical$policies)) clinical$policies[[table]] else core$policies[[table]]
  j <- model$design$J
  enrollment <- as.matrix(model$catalog[paste0("n", seq_len(j))])
  events <- as.matrix(model$catalog[paste0("x", seq_len(j))])
  if (length(counts) != j || length(dlts) != j ||
      any(!is.finite(c(counts, dlts))) || any(counts != floor(counts)) ||
      any(dlts != floor(dlts)) || any(counts < 0) || any(dlts < 0) ||
      any(dlts > counts)) stop("Supply ", j, " valid integer enrollment and DLT counts.")
  matches <- which(rowSums(abs(sweep(enrollment, 2, counts))) == 0 &
                   rowSums(abs(sweep(events, 2, dlts))) == 0)
  if (length(matches) != 1L) stop("The count pattern is outside this protocol's terminal catalog.")
  state <- matches[1]
  data.frame(state = state, BOIN = model$reference[state], recommendation = action[state],
             changed = model$reference[state] != action[state])
}

audit_transport_file <- function(path, filename, model, graph) {
  action <- read.csv(file.path(path, filename))$action
  z <- check_endpoint_flow(model, action, graph)
  jsonlite::write_json(list(results = z), file.path(path,
    paste0(tools::file_path_sans_ext(filename), "_endpoint_audit.json")),
    pretty = TRUE, auto_unbox = TRUE, digits = 17)
  jsonlite::toJSON(z, dataframe = "rows")
}

check_endpoint_flow <- function(model, action, graph) {
  k <- length(model$reference)
  capacity <- integer(k)
  capacity[as.integer(names(graph$capacities))] <- as.integer(unlist(graph$capacities))
  edges <- as.matrix(graph$edges[c("endpoint", "bad", "good")])
  storage.mode(edges) <- "integer"
  endpoint_flow_audit(as.matrix(model$catalog[paste0("n", 1:model$design$J)]),
    as.matrix(model$catalog[paste0("x", 1:model$design$J)]),
    as.double(model$catalog$path_coefficient), as.integer(model$reference),
    as.integer(action), capacity, edges)
}
