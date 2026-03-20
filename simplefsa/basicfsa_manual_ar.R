load("fsa.RData") # gets "dat"
library(RTMB)
library(tmbstan)
library(rstan)
library(broom.mixed)
library(ggplot2)
library(gridExtra)
library(tictoc)
library(rstan)

par <- list(
  logN1Y=rep(0,nrow(dat$M)),
  logN1A=rep(0,ncol(dat$M)-1),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logSdCatch=0,
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logSdSurvey=0,
  logSdR=0,
  tphiR=0,
  logMuR=0
)

nll <- function(par) {
  getAll(par, dat)

  na <- max(age) - min(age) + 1
  ny <- max(year) - min(year) + 1

  ## setup F
  F <- exp(outer(logFA, logFY, "+"))

  ## setup N (Non-centered AR1)
  phi <- 2 * plogis(tphiR) - 1
  sdR <- exp(logSdR)

  # 1. Penalize innovations (z ~ N(0,1))
  ans <- -sum(dnorm(z, 0, 1, log = TRUE))

  # 2. Build the recruitment deviations manually
  # Stationary variance for the first recruitment year
  rec_dev <- numeric(ny - 1)
  rec_dev[1] <- z[1] * (sdR / sqrt(1 - phi^2))

  for(i in 2:(ny - 1)) {
    rec_dev[i] <- phi * rec_dev[i-1] + z[i] * sdR
  }

  # 3. Map to logN matrix
  logN <- matrix(0, nrow = na, ncol = ny)
  logN[, 1] <- logN1Y
  for(y in 2:ny) {
    logN[1, y] <- logMuR + rec_dev[y - 1] # Derived recruitment
    for(a in 2:na) {
      logN[a, y] <- logN[a - 1, y - 1] - F[a - 1, y - 1] - M[a - 1, y - 1]
    }
  }

  # Match to observations
  logObs <- log(obs)
  logPred <- numeric(length(logObs))
  sdvec <- numeric(length(logObs))
  for(i in 1:length(logObs)) {
    a <- age[i] - min(age) + 1
    y <- year[i] - min(year) + 1
    if(fleet[i] == 1) {
      logPred[i] <- log(F[a, y]) - log(F[a, y] + M[a, y]) + log(1 - exp(-F[a, y] - M[a, y])) + logN[a, y]
      sdvec[i] <- exp(logSdCatch)
    } else {
      logPred[i] <- logQ[a] - (F[a, y] + M[a, y]) * surveyTime + logN[a, y]
      sdvec[i] <- exp(logSdSurvey)
    }
  }

  ans <- ans - sum(dnorm(logObs, logPred, sdvec, TRUE))

  logssb <- log(apply(exp(logN) * stockMeanWeight * propMature, 2, sum))
  REPORT(logssb)
  ADREPORT(logssb)
  return(ans)
}



# This is for creating the list for Stan
na <- max(dat$age) - min(dat$age) + 1
ny <- max(dat$year) - min(dat$year) + 1

stan_data <- list(
  n_obs           = length(dat$obs),
  obs             = dat$obs,
  age             = as.integer(dat$age),
  year            = as.integer(dat$year),
  fleet           = as.integer(dat$fleet),
  na              = na,
  ny              = ny,
  min_age         = min(dat$age),
  min_year        = min(dat$year),
  M               = dat$M,
  stockMeanWeight = dat$stockMeanWeight,
  propMature      = dat$propMature,
  surveyTime      = dat$surveyTime
)

n_sim_iter <- 10
set.seed(42)
# preserve original data
orig_dat <- dat

sim_df <- data.frame(i = 1:n_sim_iter,
                     max_cor = NA,
                     mean_cor = NA,
                     sampling_time = NA,
                     pars_hi_rhat = NA,
                     max_r_hat = NA,
                     mean_n_eff = NA)
sim_df_stan <- sim_df
for(i in 1:n_sim_iter) {
  dat <- orig_dat
  # jitter data a little
  dat$obs <- dat$obs * exp( rnorm(length(dat$obs), 0, 0.01) )
  # also update stan data
  stan_data$obs <- dat$obs

  par$z <- rep(0, length(par$logN1A))
  map_nc <- list(
    logFA = factor(c(1:4, NA, NA, NA)),
    logN1A = factor(rep(NA, length(par$logN1A))) # Map off the centered parameter
  )

  # Create the objective function with 'z' as the random effect
  obj <- MakeADFun(nll,
                   par,
                   map = map_nc,
                   random = "z",
                   silent = TRUE)

  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))

  sdrep <- sdreport(obj)

  # fit the model
  tic()
  stan_fit <- tmbstan(obj, chains=4,
                      iter=3000, # 2500 burn in
                      control=list(adapt_delta=0.99))
  tictoc <- toc()

  # diagnostics look pretty good
  pars <- summary(stan_fit)$summary
  pars_hi_rhat <- which(pars[,"Rhat"] > 1.15)

  sim_df$pars_hi_rhat[i] <- length(pars_hi_rhat)
  sim_df$max_r_hat[i] <- max(pars[,"Rhat"])
  sim_df$sampling_time[i] <- as.numeric(tictoc$toc - tictoc$tic)
  sim_df$mean_n_eff[i] <- mean(pars[,"n_eff"])

  pars <- rstan::extract(stan_fit, permuted = FALSE)
  pars <- rbind(pars[,1,], pars[,2,], pars[,3,], pars[,4,])
  cormat <- cor(pars)
  sim_df$mean_cor[i] <- mean(cormat[upper.tri(cormat)])
  sim_df$max_cor[i] <- max(cormat[upper.tri(cormat)])

  # also fit the fully bayesian model with stan
  tic()
  fit <- stan(
    file = "stan_manual_ar.stan",
    data = stan_data,
    chains = 4,
    iter = 3000,
    control = list(adapt_delta = 0.99)
  )
  tictoc <- toc()

  pars <- summary(fit)$summary
  pars_hi_rhat <- which(pars[,"Rhat"] > 1.15)

  sim_df_stan$pars_hi_rhat[i] <- length(pars_hi_rhat)
  sim_df_stan$max_r_hat[i] <- max(pars[,"Rhat"])
  sim_df_stan$sampling_time[i] <- as.numeric(tictoc$toc - tictoc$tic)
  sim_df_stan$mean_n_eff[i] <- mean(pars[,"n_eff"])

  pars <- rstan::extract(fit, permuted = FALSE)
  pars <- rbind(pars[,1,], pars[,2,], pars[,3,], pars[,4,])
  cormat <- cor(pars)
  sim_df_stan$mean_cor[i] <- mean(cormat[upper.tri(cormat)])
  sim_df_stan$max_cor[i] <- max(cormat[upper.tri(cormat)])

}

sim_df$sampling <- "tmbstan"
sim_df_stan$sampling <- "stan"

saveRDS(rbind(sim_df$sampling, sim_df_stan$sampling), "ar_results_manual.rds")


