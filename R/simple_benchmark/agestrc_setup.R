library(RTMB)
library(Matrix)
library(bench)

# model and simulator
source("R/simple_benchmark/agestrc_model.R")

## ---- constants -------------------------------------------------------------

n_a <- 6 # oldest age is the plus group, only n_t varies
n_t_grid <- c(20, 40, 60, 80, 100, 150, 200)
n_t_max <- max(n_t_grid) # sim once this long, subset to n_t
n_t_warm <- 20
n_iter <- 5 # bench::mark iterations
random_seed <- 1

# control for nlminb
ctrl <- list(iter.max = 1000, eval.max = 1000)

# fix M, the selectivity intercept and the catch SD
map <- list(logF_a = factor(c(NA, seq_len(n_a - 1))), logM = factor(NA),
            logsigmaC = factor(NA))

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

## ---- simulation ------------------------------------------------------------

#' Simulate the two datasets, one long series each
#'
#' @param rec_proc "ar1" or "ma1"
agestrc_sims <- function(rec_proc) {
  sim_one <- function(plus_group, sigmaS) {
    sim_age_strc_data(n_a = n_a, n_t = n_t_max, M = 0.2, sigmaS = sigmaS,
                      plus_group = plus_group, rec_proc = rec_proc,
                      phiR = 0.6, thetaR = 0.6, seed = random_seed)
  }
  list(np = sim_one(FALSE, 0),
       pg = sim_one(TRUE,  0))
}

## ---- one configuration -----------------------------------------------------

#' Data, pars and random names for one config
#'
#' Only the pars a config uses, so nothing is left unidentified.
#'
#' @param sims From agestrc_sims()
#' @param treatment Name in treatments
#' @param rec_re "state" or "innovation"
#' @param rec_proc "ar1" or "ma1"
#' @param n_t Number of years
agestrc_config <- function(sims, treatment, rec_re, rec_proc, n_t) {

  trt <- treatments[[treatment]]
  sim <- sims[[trt$sim]]

  data <- list(n_a = n_a, n_t = n_t,
               ln_C = sim$lC_at[, seq_len(n_t)],
               ln_I = sim$lI_at[, seq_len(n_t)],
               plus = trt$plus, rec_proc = rec_proc, rec_re = rec_re,
               obs = "iid", F_form = "separable", surv = trt$surv,
               tiny = 1e-4, sparse_lse = trt$sparse_lse)

  pars <- list(logM = log(0.2), logR0 = 0, logsigmaR = log(0.5),
               logQ_a = rep(0, n_a), logsigmaC = log(0.1),
               logsigmaI = log(0.3), logF_a = rep(0, n_a),
               logF_t = rep(log(0.2), n_t),
               ln_Ninit = -0.2 * (1:(n_a - 1)))

  # process coefficient
  if (rec_proc == "ar1") pars$trans_phiR <- qlogis(0.75)
  if (rec_proc == "ma1") pars$trans_thetaR <- qlogis(0.75)

  # recruitment latent
  rec_par <- if (rec_re == "state") "ln_R" else "epsR"
  pars[[rec_par]] <- rep(0, n_t)
  random <- rec_par

  # plus group only
  if (trt$surv == "aux_plus") {
    pars$aux_t <- rep(-0.2 * (n_a - 1), n_t - 1)
    random <- c(random, "aux_t")
  }

  # ages 2+ latent, pinned
  if (trt$surv == "aux_all") {
    pars$ln_N2 <- matrix(-0.2 * (1:(n_a - 1)), n_a - 1, n_t - 1)
    random <- c(random, "ln_N2")
  }

  list(data = data, pars = pars, random = random)
}

#' MakeADFun() for one config
#'
#' @param cfg From agestrc_config()
agestrc_obj <- function(cfg) {
  MakeADFun(cmb(jnll_agestrc, cfg$data), cfg$pars, map = map,
            random = cfg$random, silent = TRUE)
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
#' aux pins sigmaS at tiny and routinely returns conv = 1 at the right nll, so
#' read nll and max_gr next to conv. nll only comparable within a sim.
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

        # results
        acc[[length(acc) + 1]] <-
          data.frame(rec_proc = rec_proc, n_t = n_t, treatment = treatment,
                     rec_re = rec_re, sim = treatments[[treatment]]$sim,
                     n_re = nrow(H), nnzero = nnzero(H),
                     pct_dense = 100 * nnzero(H) / nrow(H)^2,
                     nll = fit$objective,
                     max_gr = max(abs(obj$gr(fit$par))),
                     conv = fit$convergence)
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
