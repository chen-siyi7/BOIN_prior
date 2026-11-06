
boin_project <- function(root = ".", native = FALSE) {
  root <- normalizePath(root, mustWork = TRUE)
  required <- c("jsonlite", "digest")
  if (native) required <- c(required, "Rcpp", "gmp")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the required R packages: ", paste(missing, collapse = ", "))
  env <- new.env(parent = globalenv())
  env$root <- root
  env$core <- readRDS(file.path(root, "data/core.rds"))
  env$clinical <- readRDS(file.path(root, "data/clinical.rds"))
  env$secondary <- readRDS(file.path(root, "data/secondary.rds"))
  env$read_reference_csv <- function(name) read.csv(file.path(root, "data/reference", name))
  env$read_reference_json <- function(name) {
    jsonlite::read_json(file.path(root, "data/reference", name), simplifyVector = TRUE)
  }
  for (name in c("boin", "likelihood", "curves", "inputs", "risk_helpers", "solver", "fit_helpers",
                 "clinical_risk", "clinical_fit", "selection", "evaluation", "gain_checks", "clinical_gains", "gain_ceilings",
                 "verification", "continuous_priors", "practical_impact",
                 "secondary", "secondary_curves", "secondary_generator_models", "robust", "comparators", "cse", "generator_uncertainty", "secondary_fit", "bound_comparisons")) {
    sys.source(file.path(root, "R", paste0(name, ".R")), env)
  }
  env$o <- env$core$model
  env$part <- env$core$partition
  env$old <- list(qs = env$core$priors, policies = setNames(env$core$policies[c("central_original", "high_original")], c("central_locked", "high_locked")))
  env$b <- env$o$reference
  env$J <- env$o$design$J
  env$K <- length(env$b)
  env$P <- 5L * env$J - 2L
  env$nn <- as.matrix(env$o$catalog[paste0("n", seq_len(env$J))])
  env$xx <- as.matrix(env$o$catalog[paste0("x", seq_len(env$J))])
  env$lp <- env$o$catalog$log_path_coefficient
  env$anchors <- vapply(seq_len(env$J), function(j) which(env$b == j)[1], integer(1))
  env$graph <- env$core$graph
  if (native) {
    for (name in c("ordered_likelihood", "risk_audit", "clinical_audit",
                   "probabilities", "endpoint_flow", "point_bounds")) {
      Rcpp::sourceCpp(file.path(root, "src", paste0(name, ".cpp")), env = env,
                      showOutput = FALSE, verbose = FALSE)
    }
    env$exact_mass <- function(catalog, probabilities) {
      j <- ncol(probabilities)
      mass <- env$terminal_mass_cpp(
        as.matrix(catalog[paste0("n", seq_len(j))]),
        as.matrix(catalog[paste0("x", seq_len(j))]),
        as.double(catalog$path_coefficient), as.matrix(probabilities))
      if (max(abs(colSums(mass) - 1)) > 1e-9) stop("Terminal probabilities do not sum to one.")
      mass
    }
  }
  env
}

boin_output <- function(path, root) {
  if (dir.exists(path) || file.exists(path)) stop("Use a new output directory.")
  parent <- normalizePath(dirname(path), mustWork = TRUE)
  output <- file.path(parent, basename(path))
  root <- normalizePath(root)
  if (identical(output, root) || startsWith(output, paste0(root, "/")))
    stop("Use an output directory outside the supplied code folder.")
  dir.create(output)
  output
}
