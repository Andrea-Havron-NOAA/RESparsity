#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: rtmb_ar1_scaling_worker.R N")
n <- as.integer(args[[1]])

time <- seq_len(n)
dat <- list(
  y = 0.5 * sin(0.13 * time) + 0.25 * cos(0.07 * time),
  process_sd = 0.8,
  observation_sd = 0.5
)
parameters <- list(unconstrained_phi = 0.4, x = rep(0, n))
nll <- function(parameters) {
  getAll(dat, parameters)
  phi <- 2 * plogis(unconstrained_phi) - 1
  value <- -dnorm(x[1], 0, process_sd / sqrt(1 - phi^2), log = TRUE)
  for (i in 2:length(x))
    value <- value - dnorm(x[i], phi * x[i - 1], process_sd, log = TRUE)
  value - sum(dnorm(y, x, observation_sd, log = TRUE))
}

time_repeated <- function(fn, repetitions) {
  start <- proc.time()[["elapsed"]]
  for (i in seq_len(repetitions)) value <- fn()
  list(value = value,
       mean_ms = 1000 * (proc.time()[["elapsed"]] - start) / repetitions)
}

setup <- time_repeated(function() object <<- MakeADFun(
  nll, parameters, random = "x", silent = TRUE
), 1L)
first <- time_repeated(function() object$fn(object$par), 1L)
repetitions <- if (n <= 30) 1000L else if (n <= 100) 300L else if (n <= 300) 100L else 30L
warm <- time_repeated(function() object$fn(object$par), repetitions)
gradient_repetitions <- if (n <= 30) 100L else if (n <= 100) 50L else if (n <= 300) 20L else 10L
gradient_benchmark <- time_repeated(function() object$gr(object$par), gradient_repetitions)

optimization_benchmark <- time_repeated(function() nlminb(
  object$par, object$fn, object$gr,
  control = list(iter.max = 30, eval.max = 100, rel.tol = 1e-8)
), 1L)
optimization <- optimization_benchmark$value
optimization_ms <- optimization_benchmark$mean_ms
estimate <- unname(optimization$par[1])
optimization_converged <- as.integer(optimization$convergence == 0L)

hessian <- object$env$spHess(random = TRUE)
entries <- summary(hessian)
row <- data.frame(
  engine = "rtmb", n = n, setup_ms = setup$mean_ms,
  first_laplace_ms = first$mean_ms, warm_laplace_mean_ms = warm$mean_ms,
  gradient_setup_ms = NA_real_,
  gradient_ms = gradient_benchmark$mean_ms,
  gradient = unname(gradient_benchmark$value[1]),
  gradient_method = "automatic_differentiation",
  hessian_nnz = length(hessian@x),
  hessian_density = length(hessian@x) / n^2,
  hessian_unique_nnz = length(hessian@x),
  hessian_bandwidth = max(abs(entries$i - entries$j)),
  objective = first$value,
  objective_phase_ms = NA_real_, factorization_phase_ms = NA_real_,
  sensitivity_phase_ms = NA_real_, hdot_phase_ms = NA_real_,
  hdot_validation_ms = NA_real_, hdot_direction_setup_ms = NA_real_,
  hdot_reverse_ms = NA_real_, hdot_contraction_ms = NA_real_,
  trace_phase_ms = NA_real_, gradient_internal_total_ms = NA_real_,
  optimization_ms = optimization_ms, estimate = estimate,
  optimization_converged = optimization_converged,
  repetitions = repetitions, gradient_repetitions = gradient_repetitions,
  success = as.integer(all(is.finite(c(first$value, gradient_benchmark$value)))),
  detected_structure = "tridiagonal", selected_backend = "rtmb_sparse_ad",
  symbolic_reuse = NA_integer_, active_directions = NA_integer_,
  hdot_workers = NA_integer_, objective_tape_rebuilds = NA_integer_,
  hdot_tape_rebuilds = NA_integer_
)
write.csv(row, stdout(), row.names = FALSE)
