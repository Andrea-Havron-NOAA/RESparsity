library(RTMB)
library(tmbstan)
library(rstan)
library(tictoc)

rstan_options(auto_write = TRUE)
options(mc.cores = 1)

load("simplefsa/fsa.RData")

## ---- Parameter list --------------------------------------
## eps: centred innovation, eps ~ N(0, exp(logSdR))
par <- list(
  logN1Y = rep(0, nrow(dat$M)),
  eps = rep(0, ncol(dat$M)),
  logFY = rep(0, ncol(dat$M)),
  logFA = rep(0, nrow(dat$M)),
  logSdCatch = 0,
  logQ = rep(0, length(unique(dat$age[dat$fleet == 2]))),
  logSdSurvey = 0,
  logSdR = 0,
  tthetaR = 0,
  logMuR = 0
)

## ---- RTMB NLL --------------------------------------------
nll <- function(par) {
  getAll(par, dat)

  na <- max(age) - min(age) + 1
  ny <- max(year) - min(year) + 1

  ## --- Fishing mortality
  F <- exp(outer(logFA, logFY, "+"))

  ## --- Priors
  ## Centered prior on innovations
  ans <- -sum(dnorm(eps, 0, sd = exp(logSdR), log = TRUE))
  ## Priors on log-SD parameters
  ans <- ans - dnorm(logSdCatch, 0, 1, log = TRUE)
  ans <- ans - dnorm(logSdSurvey, 0, 1, log = TRUE)
  ans <- ans - dnorm(logSdR, 0, 1, log = TRUE)

  ## --- MA(1) recruitment
  ## Stan: logN[1,y] = logMuR + eps[y] + theta*eps[y-1]  for y=2..ny
  ## RTMB: logN1A[i] = logMuR + eps[i+1] + theta*eps[i]  for i=1..(ny-1)
  theta <- 2 * plogis(tthetaR) - 1
  logN1A <- numeric(ny - 1)
  for (i in 1:(ny - 1)) {
    logN1A[i] <- logMuR + eps[i + 1] + theta * eps[i]
  }

  ## --- Population dynamics
  logN <- matrix(0, nrow = na, ncol = ny)
  logN[, 1] <- logN1Y
  for (y in 2:ny) {
    logN[1, y] <- logN1A[y - 1]
    for (a in 2:na) {
      logN[a, y] <- logN[a - 1, y - 1] - F[a - 1, y - 1] - M[a - 1, y - 1]
    }
  }

  ## --- Observation likelihood
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

  ## --- Derived quantity
  logssb <- log(apply(exp(logN) * stockMeanWeight * propMature, 2, sum))
  ADREPORT(logssb)

  ans
}

## ---- Build TMB object ------------------------------------
obj <- MakeADFun(nll, par,
  map    = list(logFA = factor(c(1:4, NA, NA, NA))),
  silent = TRUE,
  random = "eps"
)
opt <- nlminb(obj$par, obj$fn, obj$gr,
  control = list(iter.max = 1000, eval.max = 1000)
)
sdrep <- sdreport(obj)

## ---- Stan model ------------------------------------------
na_dat <- max(dat$age) - min(dat$age) + 1
ny_dat <- max(dat$year) - min(dat$year) + 1

stan_model_code <- "
data {
  int<lower=1> n_obs;
  vector[n_obs] obs;
  array[n_obs] int age;
  array[n_obs] int year;
  array[n_obs] int fleet;
  int<lower=1> na;
  int<lower=1> ny;
  int<lower=1> n_logFA_free;   // free logFA pars (= 4, rest fixed at 0 via map)
  int<lower=1> n_survey_ages;  // number of ages in survey fleet (= 5)
  int<lower=1> min_age;
  int<lower=1> min_year;
  matrix[na, ny] M;
  matrix[na, ny] stockMeanWeight;
  matrix[na, ny] propMature;
  real surveyTime;
}
parameters {
  vector[na] logN1Y;
  vector[ny] eps;                  // centered: eps ~ N(0, sdR)
  vector[ny] logFY;
  vector[n_logFA_free] logFA_free; // only free logFA (last na-n_logFA_free fixed=0)
  real logSdCatch;
  vector[n_survey_ages] logQ;      // only survey ages, matching TMB
  real logSdSurvey;
  real logSdR;
  real tthetaR;
  real logMuR;
}
transformed parameters {
  matrix[na, ny] F;
  matrix[na, ny] logN;
  vector[na] logFA;              // full logFA: free entries + zeros
  real theta = 2 * inv_logit(tthetaR) - 1;
  real sdR   = exp(logSdR);

  // Construct full logFA: first n_logFA_free are estimated, rest fixed at 0
  for (a in 1:n_logFA_free) logFA[a] = logFA_free[a];
  for (a in (n_logFA_free+1):na) logFA[a] = 0.0;

  // Fishing mortality
  for (y in 1:ny)
    for (a in 1:na)
      F[a, y] = exp(logFA[a] + logFY[y]);

  // Population dynamics
  for (a in 1:na) logN[a, 1] = logN1Y[a];
  for (y in 2:ny) {
    logN[1, y] = logMuR + eps[y] + theta * eps[y - 1];
    for (a in 2:na)
      logN[a, y] = logN[a - 1, y - 1] - F[a - 1, y - 1] - M[a - 1, y - 1];
  }
}
model {
  // Priors -- must match TMB exactly
  eps ~ normal(0, sdR);
  logSdCatch  ~ normal(0, 1);
  logSdSurvey ~ normal(0, 1);
  logSdR      ~ normal(0, 1);

  // Likelihood
  for (i in 1:n_obs) {
    int a_idx = age[i]  - min_age  + 1;
    int y_idx = year[i] - min_year + 1;
    real logPred;
    real sd;
    if (fleet[i] == 1) {
      logPred = log(F[a_idx, y_idx])
              - log(F[a_idx, y_idx] + M[a_idx, y_idx])
              + log1m_exp(-F[a_idx, y_idx] - M[a_idx, y_idx])
              + logN[a_idx, y_idx];
      sd = exp(logSdCatch);
    } else {
      logPred = logQ[a_idx]
              - (F[a_idx, y_idx] + M[a_idx, y_idx]) * surveyTime
              + logN[a_idx, y_idx];
      sd = exp(logSdSurvey);
    }
    log(obs[i]) ~ normal(logPred, sd);
  }
}
//generated quantities {
//  vector[ny] ssb;
//  for (y in 1:ny)
//    ssb[y] = sum(rows_dot_product(exp(logN[, y]),
//                 stockMeanWeight[, y] .* propMature[, y]));
//}
"
stan_mod <- stan_model(model_code = stan_model_code)

n_logFA_free <- sum(!is.na(levels(factor(c(1:4, NA, NA, NA))))) # = 4
n_survey_ages <- length(unique(dat$age[dat$fleet == 2])) # = 5

stan_data <- list(
  n_obs           = length(dat$obs),
  obs             = dat$obs,
  age             = as.integer(dat$age),
  year            = as.integer(dat$year),
  fleet           = as.integer(dat$fleet),
  na              = na_dat,
  ny              = ny_dat,
  n_logFA_free    = n_logFA_free,
  n_survey_ages   = n_survey_ages,
  min_age         = min(dat$age),
  min_year        = min(dat$year),
  M               = dat$M,
  stockMeanWeight = dat$stockMeanWeight,
  propMature      = dat$propMature,
  surveyTime      = dat$surveyTime
)

# ## ---- NLL equivalence check -------------------------------
# check_nll_equivalence_ma <- function(fit_stan, obj_tmb, i = 10, j = 20,
#                                      verbose = TRUE) {
#   ## Stan lp__ delta
#   lp <- rstan::extract(fit_stan, "lp__", permuted = TRUE)$lp__
#   delta_lp_stan <- lp[i] - lp[j]
#   if (verbose) {
#     cat(sprintf("Stan   lp__[%d] - lp__[%d] = %.6f\n", i, j, delta_lp_stan))
#   }
#
#   ## Extract posterior samples
#   ap <- rstan::extract(fit_stan, permuted = TRUE)
#
#   ## TMB full par vector reference (names and ordering)
#   ref <- obj_tmb$env$last.par.best
#   par_names <- names(ref)
#
#   build_tmb_parvec <- function(row) {
#     pvec <- ref
#     ## Scalars
#     pvec["logSdCatch"] <- ap$logSdCatch[row]
#     pvec["logSdSurvey"] <- ap$logSdSurvey[row]
#     pvec["logSdR"] <- ap$logSdR[row]
#     pvec["tthetaR"] <- ap$tthetaR[row]
#     pvec["logMuR"] <- ap$logMuR[row]
#     ## Vectors
#     pvec[par_names == "logN1Y"] <- ap$logN1Y[row, ]
#     pvec[par_names == "logFY"] <- ap$logFY[row, ]
#     ## logFA: TMB stores only the 4 free entries (map removes last 3)
#     pvec[par_names == "logFA"] <- ap$logFA_free[row, ]
#     ## logQ: TMB has 5 entries (survey ages only), Stan now matches
#     pvec[par_names == "logQ"] <- ap$logQ[row, ]
#     pvec[par_names == "eps"] <- ap$eps[row, ]
#     pvec
#   }
#
#   tmb_pv_i <- build_tmb_parvec(i)
#   tmb_pv_j <- build_tmb_parvec(j)
#
#   ## Joint NLL via obj$env$f() (no Laplace marginalisation)
#   nll_i <- obj_tmb$env$f(tmb_pv_i)
#   nll_j <- obj_tmb$env$f(tmb_pv_j)
#   delta_nll_tmb <- -(nll_i - nll_j)
#
#   diff <- delta_lp_stan - delta_nll_tmb
#
#   invisible(list(
#     delta_lp_stan = delta_lp_stan,
#     delta_nll_tmb = delta_nll_tmb,
#     difference = diff
#   ))
# }
#
# ## Run the check -- do this before the full loop
# fit_check <- sampling(stan_mod,
#   data = stan_data,
#   chains = 2, iter = 500, refresh = 0
# )
# check_nll_equivalence_ma(fit_check, obj, i = 10, j = 20)
# check_nll_equivalence_ma(fit_check, obj, i = 37, j = 100)
# check_nll_equivalence_ma(fit_check, obj, i = 50, j = 200)

## ---- Full tmbstan vs Stan comparison ---------------------
n_sim_iter <- 10
iter_hmc <- 4000
adapt_delta <- 0.99
n_chains <- 4
set.seed(1234)

orig_dat <- dat

sim_df <- data.frame()
sim_df_stan <- data.frame()

for (i in 1:n_sim_iter) {
  dat <- orig_dat
  dat$obs <- dat$obs * exp(rnorm(length(dat$obs), 0, 0.01))
  stan_data$obs <- dat$obs

  ## Rebuild TMB obj with jittered data
  obj <- MakeADFun(nll, par,
    map    = list(logFA = factor(c(1:4, NA, NA, NA))),
    silent = TRUE,
    random = "eps"
  )
  opt <- nlminb(obj$par, obj$fn, obj$gr,
    control = list(iter.max = 1000, eval.max = 1000)
  )

  ## tmbstan
  tic()
  fit_tmb <- tmbstan(obj,
    chains  = n_chains,
    cores   = 1,
    iter    = iter_hmc,
    control = list(adapt_delta = adapt_delta),
    init    = "last.par.best"
  )
  tt_tmb <- toc(quiet = TRUE)

  p_tmb <- summary(fit_tmb)$summary
  row_tmb <- data.frame(
    i            = i,
    method       = "tmbstan",
    time_s       = as.numeric(tt_tmb$toc - tt_tmb$tic),
    max_rhat     = max(p_tmb[, "Rhat"], na.rm = TRUE),
    pars_hi_rhat = sum(p_tmb[, "Rhat"] > 1.05, na.rm = TRUE),
    mean_neff    = mean(p_tmb[, "n_eff"], na.rm = TRUE),
    min_neff     = min(p_tmb[, "n_eff"], na.rm = TRUE)
  )

  ## Stan
  tic()
  fit_stan <- sampling(stan_mod,
    data    = stan_data,
    chains  = n_chains,
    iter    = iter_hmc,
    control = list(adapt_delta = adapt_delta),
    refresh = 0
  )
  tt_stan <- toc(quiet = TRUE)

  p_stan <- summary(fit_stan)$summary
  row_stan <- data.frame(
    i            = i,
    method       = "stan",
    time_s       = as.numeric(tt_stan$toc - tt_stan$tic),
    max_rhat     = max(p_stan[, "Rhat"], na.rm = TRUE),
    pars_hi_rhat = sum(p_stan[, "Rhat"] > 1.05, na.rm = TRUE),
    mean_neff    = mean(p_stan[, "n_eff"], na.rm = TRUE),
    min_neff     = min(p_stan[, "n_eff"], na.rm = TRUE)
  )

  sim_df <- rbind(sim_df, row_tmb)
  sim_df_stan <- rbind(sim_df_stan, row_stan)
}

results <- rbind(sim_df, sim_df_stan)
saveRDS(results, "bayes/ma_centered_results.rds")
