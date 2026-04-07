## ============================================================
## MA bad: logN1A sampled directly, eps derived via inverse MA filter
## ============================================================

library(RTMB)
library(tmbstan)
library(rstan)
library(tictoc)
library(ggplot2)
library(gridExtra)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

load("simplefsa/fsa.RData")

## ---- Parameters ------------------------------------------
par <- list(
  logN1Y      = rep(0, nrow(dat$M)),
  logN1A      = rep(0, ncol(dat$M)),
  logFY       = rep(0, ncol(dat$M)),
  logFA       = rep(0, nrow(dat$M)),
  logSdCatch  = 0,
  logQ        = rep(0, length(unique(dat$age[dat$fleet == 2]))),
  logSdSurvey = 0,
  logSdR      = 0,
  tthetaR     = 0,
  logMuR      = 0
)

## ---- RTMB NLL --------------------------------------------
nll <- function(par) {
  getAll(par, dat)

  na <- max(age) - min(age) + 1
  ny <- max(year) - min(year) + 1

  F <- exp(outer(logFA, logFY, "+"))

  ## Derive eps from logN1A via inverse MA filter
  theta <- 2 * plogis(tthetaR) - 1
  eps <- numeric(length(logN1A)) # length ny
  eps[1] <- logN1A[1]
  for (i in 2:length(logN1A)) {
    eps[i] <- logN1A[i] - theta * eps[i - 1]
  }

  ## Prior on derived eps
  ans <- -sum(dnorm(eps, 0, exp(logSdR), log = TRUE))
  ## Priors on log-SD pars -- must match Stan
  ans <- ans - dnorm(logSdCatch, 0, 1, log = TRUE)
  ans <- ans - dnorm(logSdSurvey, 0, 1, log = TRUE)
  ans <- ans - dnorm(logSdR, 0, 1, log = TRUE)

  logN <- matrix(0, nrow = na, ncol = ny)
  logN[, 1] <- logN1Y
  for (y in 2:ny) {
    logN[1, y] <- logN1A[y] + logMuR
    for (a in 2:na) {
      logN[a, y] <- logN[a - 1, y - 1] - F[a - 1, y - 1] - M[a - 1, y - 1]
    }
  }

  obs <- OBS(obs)
  logObs <- log(obs)
  logPred <- numeric(length(logObs))
  sdvec <- numeric(length(logObs))
  for (i in 1:length(logObs)) {
    a <- age[i] - min(age) + 1
    y <- year[i] - min(year) + 1
    if (fleet[i] == 1) {
      logPred[i] <- log(F[a, y]) - log(F[a, y] + M[a, y]) +
        log(1 - exp(-F[a, y] - M[a, y])) + logN[a, y]
      sdvec[i] <- exp(logSdCatch)
    } else {
      logPred[i] <- logQ[a] - (F[a, y] + M[a, y]) * surveyTime + logN[a, y]
      sdvec[i] <- exp(logSdSurvey)
    }
  }
  ans <- ans - sum(dnorm(logObs, logPred, sdvec, log = TRUE))

  logssb <- log(apply(exp(logN) * stockMeanWeight * propMature, 2, sum))
  ADREPORT(logssb)
  ans
}

## ---- Build TMB object ------------------------------------
obj <- MakeADFun(nll, par,
  map    = list(logFA = factor(c(1:4, NA, NA, NA))),
  silent = TRUE,
  random = "logN1A"
)
opt <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(iter.max = 1000, eval.max = 1000)
)
sdrep <- sdreport(obj)

## ---- Stan model ------------------------------------------
## logN1A is sampled directly; eps derived inside transformed parameters
## via the same inverse MA filter as RTMB.
na_dat <- max(dat$age) - min(dat$age) + 1
ny_dat <- max(dat$year) - min(dat$year) + 1
n_logFA_free <- 4
n_survey_ages <- length(unique(dat$age[dat$fleet == 2]))
survey_age_positions <- sort(unique(dat$age[dat$fleet == 2])) - min(dat$age) + 1
stopifnot(all(survey_age_positions == seq_len(n_survey_ages)))

stan_model_code <- "
data {
  int<lower=1> n_obs;
  vector[n_obs] obs;
  array[n_obs] int age;
  array[n_obs] int year;
  array[n_obs] int fleet;
  int<lower=1> na;
  int<lower=1> ny;
  int<lower=1> n_logFA_free;
  int<lower=1> n_survey_ages;
  int<lower=1> min_age;
  int<lower=1> min_year;
  matrix[na, ny] M;
  matrix[na, ny] stockMeanWeight;
  matrix[na, ny] propMature;
  real surveyTime;
}
parameters {
  vector[na] logN1Y;
  vector[ny] logN1A;                // direct recruitment log-numbers
  vector[ny] logFY;
  vector[n_logFA_free] logFA_free;
  real logSdCatch;
  vector[n_survey_ages] logQ;
  real logSdSurvey;
  real logSdR;
  real tthetaR;
  real logMuR;
}
transformed parameters {
  matrix[na, ny] F;
  matrix[na, ny] logN;
  vector[na] logFA;
  real theta = 2 * inv_logit(tthetaR) - 1;
  real sdR   = exp(logSdR);
  // Derive eps via inverse MA filter (length ny)
  vector[ny] eps;
  eps[1] = logN1A[1];
  for (i in 2:ny)
    eps[i] = logN1A[i] - theta * eps[i - 1];

  for (a in 1:n_logFA_free)      logFA[a] = logFA_free[a];
  for (a in (n_logFA_free+1):na) logFA[a] = 0.0;

  for (y in 1:ny)
    for (a in 1:na)
      F[a, y] = exp(logFA[a] + logFY[y]);

  for (a in 1:na) logN[a, 1] = logN1Y[a];
  for (y in 2:ny) {
    logN[1, y] = logN1A[y] + logMuR;
    for (a in 2:na)
      logN[a, y] = logN[a-1,y-1] - F[a-1,y-1] - M[a-1,y-1];
  }
}
model {
  eps         ~ normal(0, sdR);
  logSdCatch  ~ normal(0, 1);
  logSdSurvey ~ normal(0, 1);
  logSdR      ~ normal(0, 1);
  for (i in 1:n_obs) {
    int a_idx = age[i]  - min_age  + 1;
    int y_idx = year[i] - min_year + 1;
    real logPred;
    real sd;
    if (fleet[i] == 1) {
      logPred = log(F[a_idx,y_idx]) - log(F[a_idx,y_idx] + M[a_idx,y_idx])
              + log1m_exp(-F[a_idx,y_idx] - M[a_idx,y_idx]) + logN[a_idx,y_idx];
      sd = exp(logSdCatch);
    } else {
      logPred = logQ[a_idx] - (F[a_idx,y_idx] + M[a_idx,y_idx]) * surveyTime
              + logN[a_idx,y_idx];
      sd = exp(logSdSurvey);
    }
    log(obs[i]) ~ normal(logPred, sd);
  }
}

"
stan_mod <- stan_model(model_code = stan_model_code)

stan_data <- list(
  n_obs = length(dat$obs),
  obs = dat$obs,
  age = as.integer(dat$age),
  year = as.integer(dat$year),
  fleet = as.integer(dat$fleet),
  na = na_dat,
  ny = ny_dat,
  n_logFA_free = n_logFA_free,
  n_survey_ages = n_survey_ages,
  min_age = min(dat$age),
  min_year = min(dat$year),
  M = dat$M,
  stockMeanWeight = dat$stockMeanWeight,
  propMature = dat$propMature,
  surveyTime = dat$surveyTime
)

## ---- NLL equivalence check -------------------------------
# check_nll <- function(fit_stan, obj_tmb, i = 10, j = 20, verbose = TRUE) {
#   lp            <- rstan::extract(fit_stan, "lp__", permuted = TRUE)$lp__
#   delta_lp_stan <- lp[i] - lp[j]
#
#   ap        <- rstan::extract(fit_stan, permuted = TRUE)
#   ref       <- obj_tmb$env$last.par.best
#   par_names <- names(ref)
#
#   build <- function(row) {
#     pvec <- ref
#     pvec["logSdCatch"]  <- ap$logSdCatch[row]
#     pvec["logSdSurvey"] <- ap$logSdSurvey[row]
#     pvec["logSdR"]      <- ap$logSdR[row]
#     pvec["tthetaR"]     <- ap$tthetaR[row]
#     pvec["logMuR"]      <- ap$logMuR[row]
#     pvec[par_names == "logN1Y"]  <- ap$logN1Y[row, ]
#     pvec[par_names == "logFY"]   <- ap$logFY[row, ]
#     pvec[par_names == "logFA"]   <- ap$logFA_free[row, ]
#     pvec[par_names == "logQ"]    <- ap$logQ[row, ]
#     pvec[par_names == "logN1A"]  <- ap$logN1A[row, ]  # ny elements
#     pvec
#   }
#
#   nll_i         <- obj_tmb$env$f(build(i))
#   nll_j         <- obj_tmb$env$f(build(j))
#   delta_nll_tmb <- -(nll_i - nll_j)
#   diff          <- delta_lp_stan - delta_nll_tmb
#
#   if (verbose) {
#     cat(sprintf("TMB  joint log-lik[%d] - log-lik[%d] = %.6f\n", i, j, delta_nll_tmb))
#     cat(sprintf("Difference (Stan - TMB)              = %.6f\n", diff))
#   }
#   invisible(list(delta_lp_stan = delta_lp_stan, delta_nll_tmb = delta_nll_tmb, diff = diff))
# }
# fit_check <- sampling(stan_mod, data = stan_data, chains = 2, iter = 500, refresh = 0)
# check_nll(fit_check, obj, i = 10,  j = 20)
# check_nll(fit_check, obj, i = 37,  j = 100)
# check_nll(fit_check, obj, i = 50,  j = 200)

## ---- Comparison loop -------------------------------------
n_sim_iter <- 100
iter_hmc <- 4000
adapt_delta <- 0.99
n_chains <- 4
set.seed(1234)

load("simplefsa/simMA.RData")

results <- data.frame()

for (i in 1:n_sim_iter) {
  dat <- simMA[[i]]

  stan_data <- list(
    n_obs = length(dat$obs),
    obs = dat$obs,
    age = as.integer(dat$age),
    year = as.integer(dat$year),
    fleet = as.integer(dat$fleet),
    na = na_dat,
    ny = ny_dat,
    n_logFA_free = n_logFA_free,
    n_survey_ages = n_survey_ages,
    min_age = min(dat$age),
    min_year = min(dat$year),
    M = dat$M,
    stockMeanWeight = dat$stockMeanWeight,
    propMature = dat$propMature,
    surveyTime = dat$surveyTime
  )

  obj <- MakeADFun(nll, par,
    map = list(logFA = factor(c(1:4, NA, NA, NA))),
    silent = TRUE, random = "logN1A"
  )

  tic()
  opt <- nlminb(obj$par, obj$fn, obj$gr,
    control = list(iter.max = 1000, eval.max = 1000)
  )
  mle_tmb <- toc(quiet = TRUE)

  tic()
  fit_tmb <- tmbstan(obj,
    chains = n_chains, cores = 1, iter = iter_hmc,
    control = list(adapt_delta = adapt_delta),
    init = "last.par.best"
  )
  tt_tmb <- toc(quiet = TRUE)
  p_tmb <- summary(fit_tmb)$summary

  tic()
  fit_stan <- sampling(stan_mod,
    data = stan_data, chains = n_chains,
    cores = n_chains, iter = iter_hmc,
    control = list(adapt_delta = adapt_delta), refresh = 0
  )
  tt_stan <- toc(quiet = TRUE)
  p_stan <- summary(fit_stan)$summary

  row_tmbmle <- data.frame(
    i = i, method = "tmb",
    time_s = as.numeric(mle_tmb$toc - mle_tmb$tic),
    max_rhat = NA,
    pars_hi_rhat = NA,
    mean_neff = NA,
    min_neff = NA
  )

  row_tmb <- data.frame(
    i = i, method = "tmbstan",
    time_s = as.numeric(tt_tmb$toc - tt_tmb$tic),
    max_rhat = max(p_tmb[, "Rhat"], na.rm = TRUE),
    pars_hi_rhat = sum(p_tmb[, "Rhat"] > 1.15, na.rm = TRUE),
    mean_neff = mean(p_tmb[, "n_eff"], na.rm = TRUE),
    min_neff = min(p_tmb[, "n_eff"], na.rm = TRUE)
  )
  row_stan <- data.frame(
    i = i, method = "stan",
    time_s = as.numeric(tt_stan$toc - tt_stan$tic),
    max_rhat = max(p_stan[, "Rhat"], na.rm = TRUE),
    pars_hi_rhat = sum(p_stan[, "Rhat"] > 1.15, na.rm = TRUE),
    mean_neff = mean(p_stan[, "n_eff"], na.rm = TRUE),
    min_neff = min(p_stan[, "n_eff"], na.rm = TRUE)
  )

  results <- rbind(results, row_tmbmle, row_tmb, row_stan)
}

saveRDS(results, "bayes/ma_bad_results.rds")
