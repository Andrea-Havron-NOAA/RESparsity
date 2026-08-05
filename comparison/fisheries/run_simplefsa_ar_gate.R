#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(RTMB); library(Matrix)})
root <- normalizePath(getwd(), mustWork = TRUE)
gate <- Sys.getenv("SIMPLEFSA_GATE", "core")
full_fixed <- identical(gate, "full")
ma_good <- identical(gate, "ma_good")
ma_bad <- identical(gate, "ma_bad")
ma_model <- ma_good || ma_bad
full_fixed <- full_fixed || ma_model
expanded_nuisance <- identical(gate, "nuisance") || full_fixed
load(file.path(root, "simplefsa", "fsa.RData"))
observations_file <- tempfile("simplefsa-observations-", fileext = ".csv")
write.csv(data.frame(year = dat$year, fleet = dat$fleet, age = dat$age,
                     value = dat$obs), observations_file, row.names = FALSE,
          quote = FALSE)
on.exit(unlink(observations_file), add = TRUE)

quadra_root <- Sys.getenv("QUADRA_ROOT", file.path(root, "quadra"))
source_file <- file.path(root, "comparison", "fisheries", "quadra_simplefsa_ar_gate.cpp")
binary <- tempfile("quadra-simplefsa-ar-")
includes <- c(paste0("-I", quadra_root),
              paste0("-I", file.path(quadra_root, "external", "eigen")),
              paste0("-I", file.path(quadra_root, "external", "had")),
              paste0("-I", file.path(quadra_root, "external", "LBFGSpp", "include")))
status <- system2(Sys.getenv("CXX", "c++"),
                  c("-std=c++17", "-O3", includes, "-o", binary, source_file))
if (status != 0L) stop("Failed to compile Quadra simple FSA AR gate")
on.exit(unlink(binary), add = TRUE)
diagnostics_dir <- file.path(root, "comparison", "results", "reports")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
raw_dir <- file.path(root, "comparison", "results", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
diagnostics_base <- file.path(
  diagnostics_dir,
  if (ma_good) "simplefsa_ma_good_full" else if (ma_bad)
    "simplefsa_ma_bad_full" else if (full_fixed) "simplefsa_ar_full_fixed" else if (expanded_nuisance)
    "simplefsa_ar_nuisance" else "simplefsa_ar_scale"
)
binary_args <- c(observations_file, diagnostics_base)
if (expanded_nuisance) binary_args <- c(binary_args, gate)
quadra_output <- system2(binary, binary_args,
                         stdout = TRUE)
header <- grep("^engine,", quadra_output)
record <- grep("^quadra,", quadra_output)
if (length(header) != 1L || length(record) != 1L)
  stop("Quadra simple FSA output did not contain one CSV record")
quadra <- read.csv(text = paste(quadra_output[c(header, record)], collapse = "\n"))

data <- list(year = dat$year, fleet = dat$fleet, age = dat$age, obs = dat$obs,
             M = dat$M, surveyTime = dat$surveyTime)
parameters <- if (full_fixed) {
  fixed_parameters <- list(logSdR = 0, logMuR = 0, logSdCatch = 0,
       logSdSurvey = 0, logQ = rep(0, 5), logN1Y = rep(0, 7),
       logFY = rep(0, 45), logFA = rep(0, 4))
  fixed_parameters <- c(list(if (ma_model) 0 else 0), fixed_parameters)
  names(fixed_parameters)[1] <- if (ma_model) "tthetaR" else "tphiR"
  c(fixed_parameters,
    if (ma_good) list(eps = rep(0, ncol(dat$M))) else if (ma_bad)
      list(logN1A = rep(0, ncol(dat$M))) else
      list(logN1A = rep(0, ncol(dat$M) - 1)))
} else if (expanded_nuisance) {
  list(tphiR = 0, logSdR = 0, logMuR = 0, logSdCatch = 0,
       logSdSurvey = 0, logQ = rep(0, 5),
       logN1A = rep(0, ncol(dat$M) - 1))
} else {
  list(tphiR = 0, logSdR = 0,
       logN1A = rep(0, ncol(dat$M) - 1))
}
nll <- function(parameters) {
  getAll(data, parameters)
  phi <- 2 * plogis(if (ma_model) tthetaR else tphiR) - 1
  recruitment_sd <- exp(logSdR)
  recruitment_mean <- if (expanded_nuisance) logMuR else 0
  catch_sd <- if (expanded_nuisance) exp(logSdCatch) else 1
  survey_sd <- if (expanded_nuisance) exp(logSdSurvey) else 1
  conditional_sd <- recruitment_sd * sqrt(1 - phi^2)
  if (!ma_model) {
    recruits <- logN1A
    value <- -dnorm(0, recruitment_mean, recruitment_sd, log = TRUE)
    previous <- 0
    for (i in seq_along(recruits)) {
      conditional_mean <- recruitment_mean + phi * (previous - recruitment_mean)
      value <- value - dnorm(recruits[i], conditional_mean, conditional_sd,
                             log = TRUE)
      previous <- recruits[i]
    }
  } else if (ma_good) {
    value <- -sum(dnorm(eps, 0, recruitment_sd, log = TRUE))
    recruits <- recruitment_mean + eps[2:45] + phi * eps[1:44]
  } else {
    innovations <- logN1A * 0
    innovations[1] <- logN1A[1]
    for (i in 2:45) innovations[i] <- logN1A[i] - phi * innovations[i - 1]
    value <- -sum(dnorm(innovations, 0, recruitment_sd, log = TRUE))
    recruits <- logN1A[2:45] + recruitment_mean
  }
  logN <- matrix(recruits[1] * 0, nrow = 7, ncol = 45)
  fishing <- if (full_fixed) exp(outer(c(logFA, rep(0, 3)), logFY, "+")) else
    matrix(1, nrow = 7, ncol = 45)
  if (full_fixed) logN[, 1] <- logN1Y
  for (y in 2:45) {
    logN[1, y] <- recruits[y - 1]
    for (a in 2:7)
      logN[a, y] <- logN[a - 1, y - 1] - fishing[a - 1, y - 1] -
        M[a - 1, y - 1]
  }
  logPred <- obs * 0
  for (i in seq_along(obs)) {
    a <- age[i]
    y <- year[i] - 1963 + 1
    if (fleet[i] == 1) {
      total <- fishing[a, y] + M[a, y]
      logPred[i] <- log(fishing[a, y]) - log(total) +
        log(1 - exp(-total)) + logN[a, y]
    } else {
      q <- if (expanded_nuisance && a <= 5) logQ[a] else 0
      logPred[i] <- q - (fishing[a, y] + M[a, y]) * surveyTime + logN[a, y]
    }
  }
  sdvec <- ifelse(fleet == 1, catch_sd, survey_sd)
  value - sum(dnorm(log(obs), logPred, sdvec, log = TRUE))
}
random_name <- if (ma_good) "eps" else "logN1A"
object <- MakeADFun(nll, parameters, random = random_name, silent = TRUE)
objective <- object$fn(object$par)
gradient <- unname(object$gr(object$par))
fit_control <- if (full_fixed) {
  list(iter.max = 2000, eval.max = 5000, rel.tol = 1e-10)
} else {
  list(iter.max = 100, eval.max = 500, rel.tol = 1e-10)
}
fit_time <- system.time(fit <- nlminb(object$par, object$fn, object$gr,
                                      control = fit_control))[["elapsed"]]
fit_gradient <- unname(object$gr(fit$par))
if (full_fixed && ma_model) {
  seeded_start <- tempfile("simplefsa-fixed-start-", fileext = ".txt")
  writeLines(format(unname(fit$par), digits = 17, scientific = TRUE),
             seeded_start)
  on.exit(unlink(seeded_start), add = TRUE)
  quadra_output <- system2(binary, c(binary_args, seeded_start), stdout = TRUE)
  header <- grep("^engine,", quadra_output)
  record <- grep("^quadra,", quadra_output)
  if (length(header) != 1L || length(record) != 1L)
    stop("Seeded Quadra fisheries output did not contain one CSV record")
  quadra <- read.csv(text = paste(quadra_output[c(header, record)],
                                  collapse = "\n"))
}
if (full_fixed) {
  parameter_header <- grep("^parameter,name,", quadra_output)
  parameter_rows <- grep("^parameter,", quadra_output)
  parameter_rows <- setdiff(parameter_rows, parameter_header)
  quadra_parameters <- read.csv(text = paste(
    quadra_output[c(parameter_header, parameter_rows)], collapse = "\n"))
  fixed_names <- c(if (ma_model) "tthetaR" else "tphiR", "logSdR", "logMuR", "logSdCatch", "logSdSurvey",
                   paste0("logQ_", 1:5), paste0("logN1Y_", 1:7),
                   paste0("logFY_", 1:45), paste0("logFA_", 1:4))
  rtmb_parameters <- data.frame(
    parameter = "parameter", name = fixed_names,
    initial_gradient = gradient, estimate = unname(fit$par))
  parameter_results <- merge(quadra_parameters, rtmb_parameters, by = "name",
                             suffixes = c("_quadra", "_rtmb"), sort = FALSE)
  parameter_results$gradient_difference <- abs(
    parameter_results$initial_gradient_quadra -
      parameter_results$initial_gradient_rtmb)
  parameter_results$estimate_difference <- abs(
      parameter_results$estimate_quadra - parameter_results$estimate_rtmb)
  if (ma_model) {
    probe_fields <- strsplit(quadra_output[grep("^probe,", quadra_output)],
                             ",", fixed = TRUE)[[1]]
    probe_par <- object$par
    probe_par[1] <- 0.75
    probe_differences <- c(
      objective = abs(as.numeric(probe_fields[3]) - object$fn(probe_par)),
      gradient = abs(as.numeric(probe_fields[4]) - object$gr(probe_par)[1]))
    print(probe_differences, digits = 6)
    stopifnot(max(probe_differences) < 1e-6)
  }
  rtmb <- data.frame(
    engine = "rtmb", objective = objective,
    gradient_norm = sqrt(sum(gradient^2)), fit_objective = fit$objective,
    fit_gradient_norm = sqrt(sum(fit_gradient^2)),
    optimization_converged = as.integer(fit$convergence == 0L),
    structure = if (ma_bad) "dense" else "tridiagonal",
    backend = "rtmb_sparse_ad",
    fixed_effects = length(fit$par), random_effects = if (ma_model) 45 else 44,
    active_directions = NA_integer_, hdot_workers = NA_integer_)
  results <- rbind(quadra, rtmb)
  print(results, row.names = FALSE, digits = 14)
  cat(sprintf("max initial-gradient difference: %.4e\n",
              max(parameter_results$gradient_difference)))
  cat(sprintf("max estimate difference: %.4e\n",
              max(parameter_results$estimate_difference)))
  if (ma_model && max(parameter_results$estimate_difference) >= 1e-3)
    print(head(parameter_results[order(-parameter_results$estimate_difference), ],
               12), row.names = FALSE, digits = 8)
  stopifnot(all(results$optimization_converged == 1L),
            abs(diff(results$objective)) < 1e-7,
            max(parameter_results$gradient_difference) < 1e-6,
            max(parameter_results$estimate_difference) < 1e-3,
            abs(diff(results$fit_objective)) < 1e-6,
            quadra$structure == if (ma_bad) "dense" else "tridiagonal")
  result_stem <- if (ma_good) "simplefsa_ma_good_full" else if (ma_bad)
    "simplefsa_ma_bad_full" else "simplefsa_ar_full_fixed"
  write.csv(results, file.path(raw_dir, paste0(result_stem, "_gate.csv")),
            row.names = FALSE)
  write.csv(parameter_results,
            file.path(raw_dir, paste0(result_stem, "_parameters.csv")),
            row.names = FALSE)
  cat(sprintf("RTMB optimization elapsed: %.3f s\n", fit_time))
  cat(sprintf("PASS: complete simple FSA %s fixed blocks agree across Quadra and RTMB.\n",
              if (ma_good) "MA-good" else if (ma_bad) "MA-bad" else "AR"))
  quit(save = "no", status = 0L)
}
if (expanded_nuisance) {
  rtmb <- data.frame(
    engine = "rtmb", objective = objective,
    gradient_phi = gradient[1], gradient_log_sd = gradient[2],
    gradient_mu = gradient[3], gradient_catch_sd = gradient[4],
    gradient_survey_sd = gradient[5], gradient_q1 = gradient[6],
    gradient_q2 = gradient[7], gradient_q3 = gradient[8],
    gradient_q4 = gradient[9], gradient_q5 = gradient[10],
    estimate_phi = unname(fit$par[1]),
    estimate_log_sd = unname(fit$par[2]), estimate_mu = unname(fit$par[3]),
    estimate_catch_sd = unname(fit$par[4]),
    estimate_survey_sd = unname(fit$par[5]),
    estimate_q1 = unname(fit$par[6]), estimate_q2 = unname(fit$par[7]),
    estimate_q3 = unname(fit$par[8]), estimate_q4 = unname(fit$par[9]),
    estimate_q5 = unname(fit$par[10]), fit_objective = fit$objective,
    fit_gradient_norm = sqrt(sum(fit_gradient^2)),
    optimization_converged = as.integer(fit$convergence == 0L),
    structure = "tridiagonal", backend = "rtmb_sparse_ad",
    random_effects = 44, active_directions = NA_integer_,
    hdot_workers = NA_integer_)
  results <- rbind(quadra, rtmb)
  compared <- c("objective", "gradient_phi", "gradient_log_sd", "gradient_mu",
                "gradient_catch_sd", "gradient_survey_sd", "estimate_phi",
                "estimate_log_sd", "estimate_mu", "estimate_catch_sd",
                "estimate_survey_sd", paste0("gradient_q", 1:5),
                paste0("estimate_q", 1:5), "fit_objective")
  differences <- setNames(vapply(compared, function(name)
    abs(diff(results[[name]])), numeric(1)), compared)
  print(results, row.names = FALSE, digits = 14)
  print(differences, digits = 5)
  stopifnot(all(results$optimization_converged == 1L),
            max(differences[c("objective", "gradient_phi", "gradient_log_sd",
                              "gradient_mu", "gradient_catch_sd",
                              "gradient_survey_sd", paste0("gradient_q", 1:5))]) < 1e-6,
            max(differences[c("estimate_phi", "estimate_log_sd", "estimate_mu",
                              "estimate_catch_sd", "estimate_survey_sd",
                              paste0("estimate_q", 1:5))]) < 1e-4,
            differences[["fit_objective"]] < 1e-7,
            quadra$structure == "tridiagonal")
  dir.create(file.path("comparison", "results"), recursive = TRUE,
             showWarnings = FALSE)
  write.csv(results, file.path(raw_dir, "simplefsa_ar_nuisance_gate.csv"),
            row.names = FALSE)
  cat(sprintf("RTMB optimization elapsed: %.3f s\n", fit_time))
  cat("PASS: simple FSA AR nuisance-block gate agrees across Quadra and RTMB.\n")
  quit(save = "no", status = 0L)
}

rtmb <- data.frame(engine = "rtmb", objective = objective,
                   gradient_phi = gradient[1], gradient_log_sd = gradient[2],
                   estimate_phi = unname(fit$par[1]),
                   estimate_log_sd = unname(fit$par[2]),
                   fit_objective = fit$objective,
                   fit_gradient_phi = fit_gradient[1],
                   fit_gradient_log_sd = fit_gradient[2],
                   optimization_converged = as.integer(fit$convergence == 0L),
                   structure = "tridiagonal", backend = "rtmb_sparse_ad",
                   random_effects = 44, active_directions = NA_integer_,
                   hdot_workers = NA_integer_)
results <- rbind(quadra, rtmb)
differences <- c(objective = abs(diff(results$objective)),
                 gradient_phi = abs(diff(results$gradient_phi)),
                 gradient_log_sd = abs(diff(results$gradient_log_sd)),
                 estimate_phi = abs(diff(results$estimate_phi)),
                 estimate_log_sd = abs(diff(results$estimate_log_sd)),
                 fit_objective = abs(diff(results$fit_objective)),
                 fit_gradient_phi = abs(diff(results$fit_gradient_phi)),
                 fit_gradient_log_sd = abs(diff(results$fit_gradient_log_sd)))
print(results, row.names = FALSE, digits = 14)
print(differences, digits = 5)
stopifnot(all(results$optimization_converged == 1L),
          differences[["objective"]] < 1e-7,
          differences[["gradient_phi"]] < 1e-6,
          differences[["gradient_log_sd"]] < 1e-6,
          differences[["estimate_phi"]] < 1e-5,
          differences[["estimate_log_sd"]] < 1e-5,
          differences[["fit_objective"]] < 1e-7,
          quadra$structure == "tridiagonal")
dir.create(file.path("comparison", "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(raw_dir, "simplefsa_ar_gate.csv"),
          row.names = FALSE)
cat(sprintf("RTMB optimization elapsed: %.3f s\n", fit_time))
cat("PASS: simple FSA AR correlation/scale Laplace gate agrees across Quadra and RTMB.\n")
