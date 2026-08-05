#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(Matrix)
})

root <- normalizePath(getwd(), mustWork = TRUE)
quadra_root <- Sys.getenv("QUADRA_ROOT", file.path(root, "quadra"))
source_file <- file.path(root, "comparison", "bayes", "quadra_ar_parameterizations.cpp")
binary <- tempfile("quadra-ar-parameterizations-")
includes <- c(
  paste0("-I", quadra_root),
  paste0("-I", file.path(quadra_root, "external", "eigen")),
  paste0("-I", file.path(quadra_root, "external", "had")),
  paste0("-I", file.path(quadra_root, "external", "LBFGSpp", "include"))
)
status <- system2(Sys.getenv("CXX", "c++"),
                  c("-std=c++17", "-O3", includes, "-o", binary, source_file))
if (status != 0L) stop("Failed to compile Quadra AR parameterization parity")
on.exit(unlink(binary), add = TRUE)
diagnostics_dir <- file.path(root, "comparison", "results", "reports")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
quadra_output <- system2(binary, diagnostics_dir, stdout = TRUE)
quadra_csv <- quadra_output[grepl("^(engine,|quadra,)", quadra_output)]
quadra <- read.csv(text = paste(quadra_csv, collapse = "\n"))

y <- c(-0.1, 0.6, 1.0, 0.5)
anchor <- -0.2

run_rtmb <- function(parameterization) {
  data <- list(y = y, anchor = anchor, process_sd = 0.8, observation_sd = 0.5)
  parameters <- list(unconstrained_phi = 0.4, latent = rep(0, 4))
  nll <- function(parameters) {
    getAll(data, parameters)
    phi <- 2 * plogis(unconstrained_phi) - 1
    scale <- process_sd * sqrt(1 - phi^2)
    states <- latent * 0
    value <- latent[1] * 0
    if (parameterization %in% c("centered", "centered_manual")) {
      states <- latent
      previous <- anchor
      for (i in seq_along(states)) {
        value <- value - dnorm(states[i], phi * previous, scale, log = TRUE)
        previous <- states[i]
      }
    } else {
      previous <- anchor
      for (i in seq_along(states)) {
        if (parameterization == "standard_innovations") {
          value <- value - dnorm(latent[i], 0, 1, log = TRUE)
          states[i] <- phi * previous + scale * latent[i]
        } else {
          value <- value - dnorm(latent[i], 0, scale, log = TRUE)
          states[i] <- phi * previous + latent[i]
        }
        previous <- states[i]
      }
    }
    value - sum(dnorm(y, states, observation_sd, log = TRUE))
  }

  object <- MakeADFun(nll, parameters, random = "latent", silent = TRUE)
  objective <- object$fn(object$par)
  gradient <- unname(object$gr(object$par)[1])
  full_mode <- object$env$last.par.best
  mode <- unname(full_mode[object$env$random])
  phi <- 2 * plogis(object$par[1]) - 1
  scale <- 0.8 * sqrt(1 - phi^2)
  if (parameterization %in% c("centered", "centered_manual")) {
    states <- mode
  } else {
    states <- numeric(4)
    previous <- anchor
    for (i in seq_along(states)) {
      innovation <- if (parameterization == "standard_innovations") scale * mode[i] else mode[i]
      states[i] <- phi * previous + innovation
      previous <- states[i]
    }
  }
  hessian <- object$env$spHess(random = TRUE)
  entries <- summary(hessian)
  bandwidth <- max(abs(entries$i - entries$j))
  structure <- if (bandwidth == 1) "tridiagonal" else "dense"
  data.frame(
    engine = "rtmb", parameterization = parameterization,
    objective = objective, gradient = gradient, structure = structure,
    backend = "rtmb_sparse_ad", state0 = states[1], state1 = states[2],
    state2 = states[3], state3 = states[4]
  )
}

parameterizations <- c("centered", "centered_manual", "scaled_innovations",
                       "standard_innovations")
rtmb <- do.call(rbind, lapply(parameterizations, run_rtmb))
results <- rbind(quadra, rtmb)
paired <- merge(subset(results, engine == "quadra"),
                subset(results, engine == "rtmb"), by = "parameterization",
                suffixes = c("_quadra", "_rtmb"))
state_columns <- paste0("state", 0:3)
state_differences <- vapply(state_columns, function(column) max(abs(
  paired[[paste0(column, "_quadra")]] - paired[[paste0(column, "_rtmb")]]
)), numeric(1))

stopifnot(
  max(abs(paired$objective_quadra - paired$objective_rtmb)) < 1e-8,
  max(abs(paired$gradient_quadra - paired$gradient_rtmb)) < 1e-8,
  max(state_differences) < 1e-8,
  max(abs(quadra$objective - quadra$objective[1])) < 1e-8,
  max(abs(quadra$gradient - quadra$gradient[1])) < 1e-8,
  all(quadra$structure[quadra$parameterization %in% c("centered", "centered_manual")] == "tridiagonal"),
  all(quadra$structure[quadra$parameterization %in% c("scaled_innovations", "standard_innovations")] == "dense")
)

dir.create(file.path("comparison", "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("comparison", "results", "raw"), recursive = TRUE,
           showWarnings = FALSE)
write.csv(results, file.path("comparison", "results", "raw", "bayes_ar_parameterization_parity.csv"),
          row.names = FALSE)
print(results, row.names = FALSE, digits = 12)
cat("PASS: AR Bayesian parameterizations agree across Quadra and RTMB.\n")
