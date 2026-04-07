## ============================================================
## AR bad: x ~ N(0,1) non-centred, manual loop with sqrt(1-phi^2) scaling
## Includes stationary prior on logN1Y[1]
## random = "x" -- dense Hessian
## Saves: ar_bad_results.rds
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

## ---- 0. Parameters ------------------------------------------
par <- list(
  logN1Y      = rep(0, nrow(dat$M)),
  x           = rep(0, ncol(dat$M) - 1),
  logFY       = rep(0, ncol(dat$M)),
  logFA       = rep(0, nrow(dat$M)),
  logSdCatch  = 0,
  logQ        = rep(0, length(unique(dat$age[dat$fleet == 2]))),
  logSdSurvey = 0,
  logSdR      = 0,
  tphiR       = 0,
  logMuR      = 0
)

## ---- RTMB NLL --------------------------------------------
nll <- function(par) {
  getAll(par, dat)

  na <- max(age) - min(age) + 1
  ny <- max(year) - min(year) + 1

  F <- exp(outer(logFA, logFY, "+"))
  sdAR <- exp(logSdR)
  phiAR <- 2 * plogis(tphiR) - 1

  ## x ~ N(0,1) non-centred prior
  ans <- -sum(dnorm(x, 0, 1, log = TRUE))
  ## Stationary marginal prior on initial log-recruitment
  ans <- ans - dnorm(logN1Y[1], logMuR, sdAR, log = TRUE)
  ## Priors on log-SD pars -- must match Stan
  ans <- ans - dnorm(logSdCatch, 0, 1, log = TRUE)
  ans <- ans - dnorm(logSdSurvey, 0, 1, log = TRUE)
  ans <- ans - dnorm(logSdR, 0, 1, log = TRUE)

  ## AR(1) with stationary variance scaling: innovation sd = sdAR*sqrt(1-phi^2)
  logN1A <- numeric(ny - 1)
  logN1A[1] <- phiAR * logN1Y[1] + logMuR * (1 - phiAR) +
    sdAR * sqrt(1 - phiAR^2) * x[1]
  for (i in 2:(ny - 1)) {
    logN1A[i] <- phiAR * logN1A[i - 1] + logMuR * (1 - phiAR) +
      sdAR * sqrt(1 - phiAR^2) * x[i]
  }

  logN <- matrix(0, nrow = na, ncol = ny)
  logN[, 1] <- logN1Y
  for (y in 2:ny) {
    logN[1, y] <- logN1A[y - 1]
    for (a in 2:na) {
      logN[a, y] <- logN[a - 1, y - 1] - F[a - 1, y - 1] - M[a - 1, y - 1]
    }
  }

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
  random = "x"
)
opt <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(iter.max = 1000, eval.max = 1000)
)
sdrep <- sdreport(obj)

## ---- Stan model ------------------------------------------
## x ~ std_normal() non-centred; logN1A computed inside transformed parameters.
## logN1Y[1] has stationary marginal prior N(logMuR, sdAR).
na_dat <- max(dat$age) - min(dat$age) + 1
ny_dat <- max(dat$year) - min(dat$year) + 1
n_logFA_free <- 4
n_survey_ages <- length(unique(dat$age[dat$fleet == 2]))
survey_age_positions <- sort(unique(dat$age[dat$fleet == 2])) - min(dat$age) + 1

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
  vector[ny - 1] x;
  vector[ny] logFY;
  vector[n_logFA_free] logFA_free;
  real logSdCatch;
  vector[n_survey_ages] logQ;
  real logSdSurvey;
  real logSdR;
  real tphiR;
  real logMuR;
}
transformed parameters {
  matrix[na, ny] F;
  matrix[na, ny] logN;
  vector[na] logFA;
  vector[ny - 1] logN1A;
  real phi  = 2 * inv_logit(tphiR) - 1;
  real sdAR = exp(logSdR);

  for (a in 1:n_logFA_free)      logFA[a] = logFA_free[a];
  for (a in (n_logFA_free+1):na) logFA[a] = 0.0;

  for (y in 1:ny)
    for (a in 1:na)
      F[a, y] = exp(logFA[a] + logFY[y]);

  // AR(1) with stationary variance scaling
  logN1A[1] = phi * logN1Y[1] + logMuR * (1 - phi) +
              sdAR * sqrt(1 - phi^2) * x[1];
  for (i in 2:(ny - 1))
    logN1A[i] = phi * logN1A[i-1] + logMuR * (1 - phi) +
                sdAR * sqrt(1 - phi^2) * x[i];

  for (a in 1:na) logN[a, 1] = logN1Y[a];
  for (y in 2:ny) {
    logN[1, y] = logN1A[y - 1];
    for (a in 2:na)
      logN[a, y] = logN[a-1,y-1] - F[a-1,y-1] - M[a-1,y-1];
  }
}
model {
  x           ~ std_normal();
  logN1Y[1]   ~ normal(logMuR, sdAR);
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
#     pvec["tphiR"]       <- ap$tphiR[row]
#     pvec["logMuR"]      <- ap$logMuR[row]
#     pvec[par_names == "logN1Y"] <- ap$logN1Y[row, ]
#     pvec[par_names == "logFY"]  <- ap$logFY[row, ]
#     pvec[par_names == "logFA"]  <- ap$logFA_free[row, ]
#     pvec[par_names == "logQ"]   <- ap$logQ[row, ]
#     pvec[par_names == "x"]      <- ap$x[row, ]
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
#
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

load("simplefsa/simAR.RData")
results <- data.frame()

for (i in 1:n_sim_iter) {
  dat <- simAR[[i]]

  obj <- MakeADFun(nll, par,
    map = list(logFA = factor(c(1:4, NA, NA, NA))),
    silent = TRUE, random = "x"
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

saveRDS(results, "bayes/ar_bad_results.rds")
