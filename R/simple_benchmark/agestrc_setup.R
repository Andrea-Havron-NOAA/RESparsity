library(RTMB)
library(Matrix)
library(bench)

# model and simulator
source("R/simple_benchmark/agestrc_model.R")
source("R/simple_benchmark/agestrc_fsa_model.R")

## ---- constants -------------------------------------------------------------

# biology and survey design are from the simplefsa data set
fsa <- local({
  load("simplefsa/fsa.RData")
  dat
})
n_a <- nrow(fsa$M) # oldest age is the plus group, only n_t varies
M_a <- fsa$M[, 1] # natural mortality at age, constant across years in the data
stopifnot(all(fsa$M == M_a))
srv_ages <- sort(unique(fsa$age[fsa$fleet == 2])) - min(fsa$age) + 1
survey_time <- fsa$surveyTime

n_t_grid <- c(20, 40, 60, 80, 100, 150, 200)
n_t_max <- max(n_t_grid) # sim once this long, subset to n_t
n_t_warm <- 20
n_iter <- 5 # bench::mark iterations
random_seed <- 1

# control for nlminb
ctrl <- list(iter.max = 1000, eval.max = 1000)

# shared F age effect for ages 5+, as in the simulation component
map <- list(logF_a = factor(c(1:4, rep(NA, n_a - 4))))

## ---- treatments ------------------------------------------------------------
treatments <- list(
  "none" = list(
    plus = "none",
    surv = "deterministic",
    sparse_lse = TRUE,
    sim = "np"
  ),
  "naive" = list(
    plus = "naive",
    surv = "deterministic",
    sparse_lse = TRUE,
    sim = "pg"
  ),
  "logspace (LSE0)" = list(
    plus = "logspace",
    surv = "deterministic",
    sparse_lse = TRUE,
    sim = "pg"
  ),
  "aux_plus" = list(
    plus = "naive",
    surv = "aux_plus",
    sparse_lse = TRUE,
    sim = "pg"
  ),
  "aux_all" = list(
    plus = "naive",
    surv = "aux_all",
    sparse_lse = TRUE,
    sim = "pg"
  )
)

# states of deviations
rec_res <- c("state", "innovation")

## ---- model objects ---------------------------------------------------------

#' Data, pars and random names for one model
#'
#' Only the pars a config uses, so nothing is left unidentified.
#'
#' @param ln_C Log catch at age, n_a x n_t
#' @param ln_I Log survey index for the survey ages, x n_t (NA = no survey)
#' @param n_t Number of years
#' @param plus,surv,sparse_lse Plus group and survival options of the model
#' @param rec_re "state" or "innovation"
#' @param rec_proc "ar1" or "ma1"
#' @param start Starting values (logR0, ln_Ninit, logQ_a, logF_t, logsigma,
#'   trans). NULL starts everything at zero, as the simulation component does.
agestrc_build <- function(ln_C, ln_I, n_t, plus, surv, sparse_lse, rec_re,
                          rec_proc, start = NULL) {

  data <- list(n_a = n_a, n_t = n_t, ln_C = ln_C, ln_I = ln_I, M_a = M_a,
               srv_ages = srv_ages, survey_time = survey_time,
               plus = plus, rec_proc = rec_proc, rec_re = rec_re, surv = surv,
               tiny = 1e-4, sparse_lse = sparse_lse, sd_prior = TRUE)

  if (is.null(start))
    start <- list(logR0 = 0, ln_Ninit = rep(0, n_a),
                  logQ_a = rep(0, length(srv_ages)), logF_t = 0, logsigma = 0,
                  trans = 0)

  pars <- list(logR0 = start$logR0, logsigmaR = start$logsigma,
               logQ_a = start$logQ_a, logsigmaC = start$logsigma,
               logsigmaI = start$logsigma, logF_a = rep(0, n_a),
               logF_t = rep(start$logF_t, n_t), ln_Ninit = start$ln_Ninit)

  # process coefficient
  if (rec_proc == "ar1") pars$trans_phiR <- start$trans
  if (rec_proc == "ma1") pars$trans_thetaR <- start$trans

  # recruitment latent: years 2..n_t for AR(1); MA(1) carries one extra element
  if (rec_proc == "ar1" && rec_re == "state") pars$ln_R <- rep(start$logR0, n_t - 1)
  if (rec_proc == "ar1" && rec_re == "innovation") pars$epsR <- rep(0, n_t - 1)
  if (rec_proc == "ma1" && rec_re == "innovation") pars$epsR <- rep(0, n_t)
  if (rec_proc == "ma1" && rec_re == "state") pars$ln_R <- rep(0, n_t)
  random <- if (rec_re == "state") "ln_R" else "epsR"

  # numbers-at-age latents start at the numbers implied by the starting values
  if (surv != "deterministic") {
    Z0 <- exp(outer(pars$logF_a, pars$logF_t, "+")) + M_a
    lnN0 <- det_lnN(pars$ln_Ninit, rep(pars$logR0, n_t - 1), Z0,
                    plus_group = plus != "none")
  }

  # plus group only
  if (surv == "aux_plus") {
    pars$aux_t <- lnN0[n_a, 2:n_t]
    random <- c(random, "aux_t")
  }

  # ages 2+ latent: pinned sd (aux_all) or estimated sd (state)
  if (surv == "aux_all" || surv == "state") {
    pars$ln_N2 <- lnN0[2:n_a, 2:n_t]
    random <- c(random, "ln_N2")
  }
  if (surv == "state") pars$logsigmaS <- log(0.1)

  list(data = data, pars = pars, random = random)
}

#' MakeADFun() for one config
#'
#' @param cfg From agestrc_build() or agestrc_config()
agestrc_obj <- function(cfg) {
  MakeADFun(cmb(jnll_agestrc_fsa, cfg$data), cfg$pars, map = map,
            random = cfg$random, silent = TRUE)
}

## ---- simulation ------------------------------------------------------------

#' simplefsa observations as matrices (ln_I is NA where there is no survey)
agestrc_fsa_data <- function() {
  n_t <- ncol(fsa$M)
  a <- fsa$age - min(fsa$age) + 1
  y <- fsa$year - min(fsa$year) + 1
  ln_C <- matrix(NA_real_, n_a, n_t)
  ln_I <- matrix(NA_real_, length(srv_ages), n_t)
  k <- fsa$fleet == 1
  ln_C[cbind(a[k], y[k])] <- log(fsa$obs[k])
  k <- fsa$fleet == 2
  ln_I[cbind(match(a[k], srv_ages), y[k])] <- log(fsa$obs[k])
  list(ln_C = ln_C, ln_I = ln_I, n_t = n_t)
}

#' Values the benchmark data are simulated from
#'
#' The aligned model with a deterministic plus group is fitted to the simplefsa
#' data and its estimates are used as the truth, so the simulated data resemble
#' the data the simulation component is built on.
#'
#' @param rec_proc "ar1" or "ma1"
agestrc_fsa_truth <- function(rec_proc) {
  d <- agestrc_fsa_data()
  rec_re <- if (rec_proc == "ar1") "state" else "innovation"
  cfg <- agestrc_build(d$ln_C, d$ln_I, d$n_t, plus = "naive",
                       surv = "deterministic", sparse_lse = TRUE,
                       rec_re = rec_re, rec_proc = rec_proc)
  obj <- agestrc_obj(cfg)
  fit <- nlminb(obj$par, obj$fn, obj$gr, control = ctrl)
  est <- obj$env$parList(fit$par)
  list(M_a = M_a, srv_ages = srv_ages, survey_time = survey_time,
       logF_a = est$logF_a, logF_t = est$logF_t, logQ_a = est$logQ_a,
       sigmaC = exp(est$logsigmaC), sigmaI = exp(est$logsigmaI),
       logR0 = est$logR0, sigmaR = exp(est$logsigmaR),
       phi = if (rec_proc == "ar1") 2 * plogis(est$trans_phiR) - 1 else NA,
       theta = if (rec_proc == "ma1") 2 * plogis(est$trans_thetaR) - 1 else NA,
       ln_Ninit = est$ln_Ninit, nll = fit$objective, conv = fit$convergence)
}

#' Simulate the two datasets, one long series each
#'
#' @param rec_proc "ar1" or "ma1"
agestrc_sims <- function(rec_proc) {
  truth <- agestrc_fsa_truth(rec_proc)
  sim_one <- function(plus_group) {
    sim_agestrc_fsa(truth, n_t = n_t_max, rec_proc = rec_proc,
                    plus_group = plus_group, seed = random_seed)
  }
  list(np = sim_one(FALSE),
       pg = sim_one(TRUE),
       truth = truth)
}

## ---- one configuration -----------------------------------------------------

#' Data, pars and random names for one benchmark config
#'
#' @param sims From agestrc_sims()
#' @param treatment Name in treatments
#' @param rec_re "state" or "innovation"
#' @param rec_proc "ar1" or "ma1"
#' @param n_t Number of years
agestrc_config <- function(sims, treatment, rec_re, rec_proc, n_t) {

  trt <- treatments[[treatment]]
  sim <- sims[[trt$sim]]

  # scale parameters start at the simulated values, the rest at generic values
  start <- list(logR0 = sims$truth$logR0, ln_Ninit = sims$truth$ln_Ninit,
                logQ_a = sims$truth$logQ_a, logF_t = log(0.4),
                logsigma = log(0.3), trans = qlogis(0.75))

  agestrc_build(sim$lC_at[, seq_len(n_t)], sim$lI_at[, seq_len(n_t)], n_t,
                plus = trt$plus, surv = trt$surv, sparse_lse = trt$sparse_lse,
                rec_re = rec_re, rec_proc = rec_proc, start = start)
}

## ---- warm-up ---------------------------------------------------------------

#' Warm-up, to remove the effects of C++ linking and memory allocation
#'
#' @param rec_proc "ar1" or "ma1"
#' @param reps Number of passes
agestrc_warmup <- function(rec_proc, reps = 3) {
  sims <- agestrc_sims(rec_proc)
  for (ii in seq_len(reps)) {
    for (treatment in names(treatments)) {
      for (rec_re in rec_res) {
        cfg <- agestrc_config(sims, treatment, rec_re, rec_proc, n_t_warm)
        obj <- agestrc_obj(cfg)
        nlminb(obj$par, obj$fn, obj$gr, control = ctrl)
      }
    }
  }
  invisible(NULL)
}

## ---- diagnostics -----------------------------------------------------------

#' Fit every config once, to see what is being timed
#'
#' @param rec_proc "ar1" or "ma1"
agestrc_diagnostics <- function(rec_proc) {

  sims <- agestrc_sims(rec_proc)
  acc <- list()

  for (n_t in n_t_grid) {
    for (treatment in names(treatments)) {
      for (rec_re in rec_res) {

        cfg <- agestrc_config(sims, treatment, rec_re, rec_proc, n_t) # get config
        obj <- agestrc_obj(cfg) # get obj
        fit <- nlminb(obj$par, obj$fn, obj$gr, control = ctrl) # optimize
        H <- obj$env$spHess(obj$env$last.par.best, random = TRUE) # get sparsity pattern
        pd <- sdreport(obj)$pdHess # Hessian of the fixed effects positive definite at the optimum

        # results
        acc[[length(acc) + 1]] <-
          data.frame(rec_proc = rec_proc, n_t = n_t, treatment = treatment,
                     rec_re = rec_re, sim = treatments[[treatment]]$sim,
                     n_re = nrow(H), nnzero = nnzero(H),
                     pct_dense = 100 * nnzero(H) / nrow(H)^2,
                     nll = fit$objective,
                     max_gr = max(abs(obj$gr(fit$par))),
                     conv = fit$convergence, msg = fit$message, pd_hess = pd)
      }
    }
  }

  do.call(rbind, acc)
}

## ---- timing ----------------------------------------------------------------

#' Time MakeADFun() and the fit against n_t
#'
#' nlminb() is not timed alone - refitting one object warm-starts the inner
#' Newton solve, which flatters whoever goes second. Fresh object inside the
#' timed expression keeps every fit a cold start, nlminb comes back as
#' total - obj.
#'
#' @param rec_proc "ar1" or "ma1"
agestrc_timings <- function(rec_proc) {

  sims <- agestrc_sims(rec_proc)

  bench::press(
    n_t = n_t_grid,
    treatment = names(treatments),
    rec_re = rec_res,
    {
      gc(reset = TRUE)
      cfg <- agestrc_config(sims, treatment, rec_re, rec_proc, n_t)

      bench::mark(
        obj = {
          # AD object
          obj <- agestrc_obj(cfg)
        },

        total = {
          # AD object and fit, cold start
          obj <- agestrc_obj(cfg)
          fit <- nlminb(obj$par, obj$fn, obj$gr, control = ctrl)
          fit$convergence
        },

        iterations = n_iter, time_unit = "s", filter_gc = FALSE, check = FALSE
      )
    }
  )
}

## ---- reshaping -------------------------------------------------------------

#' Flatten bench::press() output for plotting purposes
#'
#' @param results From agestrc_timings()
#' @param rec_proc For labelling
agestrc_tidy <- function(results, rec_proc) {

  df <- data.frame(
    rec_proc = rec_proc,
    n_t = results$n_t,
    treatment = results$treatment,
    rec_re = results$rec_re,
    func = results$expression |> attr(which = "description"),
    median_time = as.numeric(results$median),
    memory = as.numeric(results$mem_alloc)
  )

  # one row per config, both expressions side by side
  key <- c("rec_proc", "n_t", "treatment", "rec_re")
  wide <- merge(df[df$func == "obj", c(key, "median_time", "memory")],
                df[df$func == "total", c(key, "median_time", "memory")],
                by = key, suffixes = c("_obj", "_total"))

  # MakeADFun, nlminb as the difference, measured total
  panel <- function(func, median_time, memory) {
    data.frame(wide[key], func = func, median_time = median_time,
               memory = memory, row.names = NULL)
  }
  df_long <- rbind(
    panel("obj", wide$median_time_obj, wide$memory_obj),
    panel("opt", wide$median_time_total - wide$median_time_obj,
          wide$memory_total - wide$memory_obj),
    panel("Total", wide$median_time_total, wide$memory_total)
  )

  df_long$func <- factor(df_long$func, levels = c("obj", "opt", "Total"))
  df_long$treatment <- factor(df_long$treatment, levels = names(treatments))
  df_long$rec_re <- factor(df_long$rec_re, levels = rec_res)
  df_long
}

# panel labels
facet_names <- c("obj" = "MakeADFun", "opt" = "nlminb", "Total" = "Total")
rec_re_names <- c("state" = "Recruitment: states",
                  "innovation" = "Recruitment: deviations")
