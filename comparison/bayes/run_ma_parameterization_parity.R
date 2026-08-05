#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(RTMB); library(Matrix)})
root <- normalizePath(getwd(), mustWork = TRUE)
quadra_root <- Sys.getenv("QUADRA_ROOT", file.path(root, "quadra"))
source_file <- file.path(root, "comparison", "bayes", "quadra_ma_parameterizations.cpp")
binary <- tempfile("quadra-ma-parameterizations-")
includes <- c(paste0("-I", quadra_root),
              paste0("-I", file.path(quadra_root, "external", "eigen")),
              paste0("-I", file.path(quadra_root, "external", "had")),
              paste0("-I", file.path(quadra_root, "external", "LBFGSpp", "include")))
status <- system2(Sys.getenv("CXX", "c++"),
                  c("-std=c++17", "-O3", includes, "-o", binary, source_file))
if (status != 0L) stop("Failed to compile Quadra MA parameterization parity")
on.exit(unlink(binary), add = TRUE)
diagnostics_dir <- file.path(root, "comparison", "results", "reports")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
quadra_output <- system2(binary, diagnostics_dir, stdout = TRUE)
quadra_csv <- quadra_output[grepl("^(engine,|quadra,)", quadra_output)]
quadra <- read.csv(text = paste(quadra_csv, collapse = "\n"))

y <- c(0.7, 1.0, 0.1, 0.9)
run_rtmb <- function(parameterization) {
  data <- list(y = y, innovation_sd = 0.8, observation_sd = 0.5)
  parameters <- list(unconstrained_theta = -0.6, latent = rep(0, 4))
  nll <- function(parameters) {
    getAll(data, parameters)
    theta <- 2 * plogis(unconstrained_theta) - 1
    innovations <- latent * 0
    states <- latent * 0
    value <- latent[1] * 0
    if (parameterization == "innovation_good") {
      innovations <- latent
      for (i in seq_along(states)) {
        value <- value - dnorm(innovations[i], 0, innovation_sd, log = TRUE)
        states[i] <- innovations[i] + if (i == 1) 0 else theta * innovations[i - 1]
      }
    } else {
      states <- latent
      innovations[1] <- states[1]
      value <- value - dnorm(innovations[1], 0, innovation_sd, log = TRUE)
      for (i in 2:length(states)) {
        innovations[i] <- states[i] - theta * innovations[i - 1]
        value <- value - dnorm(innovations[i], 0, innovation_sd, log = TRUE)
      }
    }
    value - sum(dnorm(y, states, observation_sd, log = TRUE))
  }
  object <- MakeADFun(nll, parameters, random = "latent", silent = TRUE)
  objective <- object$fn(object$par)
  gradient <- unname(object$gr(object$par)[1])
  full_mode <- object$env$last.par.best
  mode <- unname(full_mode[object$env$random])
  theta <- 2 * plogis(object$par[1]) - 1
  states <- if (parameterization == "inverse_filter_bad") mode else
    vapply(seq_along(mode), function(i) mode[i] + if (i == 1) 0 else theta * mode[i - 1], numeric(1))
  hessian <- object$env$spHess(random = TRUE)
  entries <- summary(hessian)
  bandwidth <- max(abs(entries$i - entries$j))
  data.frame(engine = "rtmb", parameterization = parameterization,
             objective = objective, gradient = gradient,
             structure = if (bandwidth == 1) "tridiagonal" else "dense",
             backend = "rtmb_sparse_ad", state0 = states[1], state1 = states[2],
             state2 = states[3], state3 = states[4])
}

parameterizations <- c("innovation_good", "inverse_filter_bad")
rtmb <- do.call(rbind, lapply(parameterizations, run_rtmb))
results <- rbind(quadra, rtmb)
paired <- merge(subset(results, engine == "quadra"), subset(results, engine == "rtmb"),
                by = "parameterization", suffixes = c("_quadra", "_rtmb"))
state_columns <- paste0("state", 0:3)
state_difference <- max(vapply(state_columns, function(column) max(abs(
  paired[[paste0(column, "_quadra")]] - paired[[paste0(column, "_rtmb")]]
)), numeric(1)))
stopifnot(max(abs(paired$objective_quadra - paired$objective_rtmb)) < 1e-8,
          max(abs(paired$gradient_quadra - paired$gradient_rtmb)) < 1e-8,
          state_difference < 1e-8,
          abs(diff(quadra$objective)) < 1e-8,
          abs(diff(quadra$gradient)) < 1e-8,
          quadra$structure[quadra$parameterization == "innovation_good"] == "tridiagonal",
          quadra$structure[quadra$parameterization == "inverse_filter_bad"] == "dense")
dir.create(file.path("comparison", "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("comparison", "results", "raw"), recursive = TRUE,
           showWarnings = FALSE)
write.csv(results, file.path("comparison", "results", "raw", "bayes_ma_parameterization_parity.csv"), row.names = FALSE)
print(results, row.names = FALSE, digits = 12)
cat("PASS: MA Bayesian parameterizations agree across Quadra and RTMB.\n")
