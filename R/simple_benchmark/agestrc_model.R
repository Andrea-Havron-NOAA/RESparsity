## Age-structured model and simulator for the treatment benchmarks. Copied from
## the sparsity project (R/functions.R, R/estimation.R), minus the commented-out
## spatial block. Refresh by hand if either changes over there.

#' Log sum exp trick to recover sparsity
#'
#' @param x Vector
#' @param sparse Boolean on whether to use sparse mode
logspace_sum <- function(x, sparse = TRUE) {
  if(inherits(x, "advector")) {
    if(sparse) RTMB:::LSE0(x) else RTMB:::LSE(x)
  } else {
    M <- max(x)
    log(sum(exp(x - M))) + M
  }
}

#' Logistic curve at age (selectivity / maturity)
#'
#' @param a Vector of ages
#' @param a50 Age at 50%
#' @param k Slope
logistic_at <- function(a, a50, k) 1 / (1 + exp(-k * (a - a50)))

#' Stationary AR(1) correlation matrix
#'
#' @param n Dimension
#' @param rho Correlation between adjacent elements
ar1_cor <- function(n, rho) {
  idx <- seq_len(n)
  rho^abs(outer(idx, idx, "-"))
}

#' Precision matrix of a stationary AR(1) correlation matrix
#'
#' Exact tridiagonal inverse of [ar1_cor()], returned sparse, so the separable
#' precision below never has to be formed by dense inversion.
#'
#' @param n Dimension
#' @param rho Correlation between adjacent elements
ar1_prec <- function(n, rho) {
  if(n == 1) return(Matrix::Matrix(1, 1, 1, sparse = TRUE))
  d <- c(1, rep(1 + rho^2, n - 2), 1)
  Q <- Matrix::bandSparse(n, k = c(0, 1),
                          diagonals = list(d, rep(-rho, n - 1)),
                          symmetric = TRUE)
  Q / (1 - rho^2)
}

#' Draw a separable AR(1) x AR(1) field over age and year
#'
#' Draws `eps` with `cov(vec(eps)) = Sigma_t %x% Sigma_a`, where
#' `Sigma_a = sd^2 * ar1_cor(n_a, rho_a)` and `Sigma_t = ar1_cor(n_t, rho_t)`
#' (`vec` stacks ages within a year). All of the marginal variance is carried by
#' the age factor, so `sd(eps[a, t]) = sd` in every cell.
#'
#' Age margin by Cholesky, year margin by its stationary AR(1) recursion, which
#' avoids a dense n_t x n_t factor. `rho_a = rho_t = 0` gives iid N(0, sd^2).
#'
#' @param n_a Number of ages
#' @param n_t Number of years
#' @param sd Marginal standard deviation of each cell
#' @param rho_a AR(1) correlation between adjacent ages
#' @param rho_t AR(1) correlation between adjacent years
#' @return An `n_a` x `n_t` matrix of mean-zero deviates
rmvn_sep_ar1 <- function(n_a, n_t, sd, rho_a = 0, rho_t = 0) {
  if(sd == 0) return(matrix(0, n_a, n_t))       # degenerate: switches the process off
  L_a <- t(chol(sd^2 * ar1_cor(n_a, rho_a)))
  z   <- L_a %*% matrix(rnorm(n_a * n_t), n_a, n_t)  # AR(1) in age, iid in year
  eps <- matrix(0, n_a, n_t)
  eps[, 1] <- z[, 1]                                 # stationary start
  s <- sqrt(1 - rho_t^2)
  for(t in seq_len(n_t)[-1]) eps[, t] <- rho_t * eps[, t - 1] + s * z[, t]
  eps
}

#' Separable margins and precisions of an AR(1) x AR(1) field
#'
#' @inheritParams rmvn_sep_ar1
#' @return List with the dense age covariance margin and the two sparse
#'   tridiagonal precision margins. The joint covariance is
#'   `ar1_cor(n_t, rho_t) %x% Sigma_a` and the joint precision is `Q_t %x% Q_a`.
#'   The year covariance margin is deliberately not materialised: it is dense
#'   n_t x n_t, while `Q_t` is tridiagonal at any length.
sep_ar1_struct <- function(n_a, n_t, sd, rho_a = 0, rho_t = 0) {
  list(Sigma_a = sd^2 * ar1_cor(n_a, rho_a),
       Q_a     = ar1_prec(n_a, rho_a) / sd^2,
       Q_t     = ar1_prec(n_t, rho_t),
       sd = sd, rho_a = rho_a, rho_t = rho_t, n_a = n_a, n_t = n_t)
}

#' Simulate data from an age structured model
#'
#' Ages are 1:n_a with `n_a` a plus group. Fishing mortality, survival process
#' error, and both observation errors are separable AR(1) x AR(1) fields over
#' age and year: `cov(vec(eps)) = Sigma_t %x% Sigma_a`, precision `Q_t %x% Q_a`.
#' The survival field defaults to iid.
#'
#' @param n_a Number of ages (age `n_a` is a plus group)
#' @param n_t Number of years
#' @param M Natural mortality
#' @param plus_group Treat the oldest age as a plus group; FALSE makes it a true age
#' @param logR0 Mean log recruitment (also sets the year-1 equilibrium)
#' @param sigmaR Recruitment deviation SD (iid, not part of any Kronecker field)
#' @param sigmaI Marginal SD of the index-at-age observation error
#' @param sigmaC Marginal SD of the catch-at-age observation error
#' @param sigmaF Marginal SD of the log-F-at-age deviations
#' @param sigmaS SD of the survival process error, ages 2+ from year 2 on.
#'   Set to 0 for a deterministic survival recursion.
#' @param rhoI_a,rhoI_t Index error AR(1) correlation across ages / across years
#' @param rhoC_a,rhoC_t Catch error AR(1) correlation across ages / across years
#' @param rhoF_a,rhoF_t log-F AR(1) correlation across ages / across years
#' @param rec_proc Recruitment process: "iid", "ar1", or "ma1". `sigmaR` is the
#'   innovation SD, so the AR(1) marginal SD is `sigmaR / sqrt(1 - phiR^2)`.
#' @param phiR AR(1) coefficient for recruitment (`rec_proc = "ar1"`)
#' @param thetaR MA(1) coefficient for recruitment (`rec_proc = "ma1"`)
#' @param rhoS_a,rhoS_t Survival error AR(1) correlation across ages / years
#' @param q Survey catchability scalar
#' @param a50_f,k_f Fishery selectivity: age at 50% and slope
#' @param a50_s,k_s Survey selectivity: age at 50% and slope
#' @param a50_m,k_m Maturity: age at 50% and slope (declared but unused --
#'   nothing in the model is in weight, so there is no SSB)
#' @param seed Simulation seed
#'
#' @return A list holding every simulated quantity. Matrices are `n_a` x `n_t`;
#'   `*_t` vectors are length `n_t`; `*_a` vectors are length `n_a`.
#'   Observations are returned on both the natural (`C_at`, `I_at`) and log
#'   (`lC_at`, `lI_at`) scale; the error-free values carry a `_true` suffix.
#'   `cov` holds the separable margins and precisions for the three fields.
sim_age_strc_data <- function(n_a = 7,
                              n_t = 30,
                              M = 0.2,
                              plus_group = TRUE,
                              logR0   = 0,
                              sigmaR  = 0.6,
                              sigmaI  = 0.3,
                              sigmaC  = 0.1,
                              sigmaF  = 0.2,
                              sigmaS = 0.1,
                              rec_proc = c("iid", "ar1", "ma1"),
                              phiR = 0.6,
                              thetaR = 0.6,
                              rhoI_a  = 0.4,
                              rhoI_t = 0.2,
                              rhoC_a  = 0.3,
                              rhoC_t = 0.2,
                              rhoF_a  = 0.6,
                              rhoF_t = 0.8,
                              rhoS_a = 0,
                              rhoS_t = 0,
                              q = 1,
                              a50_f = 3,
                              k_f = 1,
                              a50_s = 2,
                              k_s = 1,
                              a50_m = 4,
                              k_m = 1,
                              seed = 1
                              ) {

  # simulate
  set.seed(seed)

  ages  <- seq_len(n_a)
  years <- seq_len(n_t)

  ## ---- biology and vulnerability at age -------------------------------------
  sel_a <- logistic_at(ages, a50_f, k_f)  # fishery selectivity
  srv_a <- logistic_at(ages, a50_s, k_s)  # survey selectivity
  Qa <- q * srv_a # catchability at age

  ## ---- apical F -------------------------------
  F_t <- 0.05 + 0.35 * sin(pi * (years - 1) / (n_t - 1))
  F_t <- rep_len(F_t, n_t)

  ## ---- deviations -----------------------------------------------------------
  epsR <- rnorm(n_t) # iid
  epsF <- rmvn_sep_ar1(n_a, n_t, sigmaF, rhoF_a, rhoF_t) # separabel 2dar1 by age and year
  epsC <- rmvn_sep_ar1(n_a, n_t, sigmaC, rhoC_a, rhoC_t) # separabel 2dar1 by age and year
  epsI <- rmvn_sep_ar1(n_a, n_t, sigmaI, rhoI_a, rhoI_t) # separabel 2dar1 by age and year
  epsS <- rmvn_sep_ar1(n_a, n_t, sigmaS, rhoS_a, rhoS_t) # survival, iid by default

  ## ---- Construct F and Z ----------------------------------------------------
  lF_at <- log(outer(sel_a, F_t)) + epsF
  F_at  <- exp(lF_at)
  Z_at  <- F_at + M

  # y = 1 (numbers equations)
  lN <- matrix(0, n_a, n_t)

  # recruitment process
  rec_proc <- match.arg(rec_proc)
  if (rec_proc == "iid") {
    lN[1, ] <- logR0 + sigmaR * epsR
  } else if (rec_proc == "ar1") {
    lN[1, 1] <- logR0 + sigmaR / sqrt(1 - phiR^2) * epsR[1] # stationary start
    for (t in 2:n_t)
      lN[1, t] <- phiR * lN[1, t - 1] + logR0 * (1 - phiR) + sigmaR * epsR[t]
  } else {
    epsMA <- sigmaR * epsR
    lN[1, 1] <- logR0 + epsMA[1]
    for (t in 2:n_t)
      lN[1, t] <- logR0 + epsMA[t] + thetaR * epsMA[t - 1]
  }
  lN[2:n_a, 1] <- logR0 - M * 1:(n_a - 1)
  if (plus_group)
    lN[n_a, 1] <- lN[n_a, 1] - log(1 - exp(-M)) # geometric series

  # y > 2 (numbers equations), survival process error epsS on ages 2+
  for (t in 2:n_t) {
    for (a in 2:(n_a - 1))
      lN[a, t] <- lN[a - 1, t - 1] - Z_at[a - 1, t - 1] + epsS[a, t] # non plus group
    if (plus_group) {
      lN[n_a, t] <- log(exp(lN[n_a - 1, t - 1] - Z_at[n_a - 1, t - 1]) + exp(lN[n_a  , t -
                                                                                  1] - Z_at[n_a  , t - 1])) + epsS[n_a, t]
    } else {
      lN[n_a, t] <- lN[n_a - 1, t - 1] - Z_at[n_a - 1, t - 1] + epsS[n_a, t]
    }
  } # end t loop

  N_at <- exp(lN)

  # Catch equation (Baranov) w/ observation error
  lC_at_true <- log(F_at / Z_at * (1 - exp(-Z_at))) + lN
  lC_at      <- lC_at_true + epsC

  # Index-at-age w/ observation error
  lI_at_true <- log(Qa) + lN
  lI_at      <- lI_at_true + epsI

  ## ---- derived quantities ---------------------------------------------------
  R_t    <- N_at[1, ] # recruitment
  Fbar_t <- apply(F_at, 2, max) # apical F

  # return list
  list(
    # dimensions
    n_a = n_a,
    n_t = n_t,
    ages = ages,
    years = years,

    # observations (natural and log scale)
    I_at = exp(lI_at),
    lI_at = lI_at,
    C_at = exp(lC_at),
    lC_at = lC_at,

    # error-free counterparts
    I_at_true = exp(lI_at_true),
    lI_at_true = lI_at_true,
    C_at_true = exp(lC_at_true),
    lC_at_true = lC_at_true,

    # states
    lN = lN,
    N_at = N_at,
    F_at = F_at,
    lF_at = lF_at,
    Z_at = Z_at,

    # deviations actually drawn
    epsR = epsR,
    epsF = epsF,
    epsC = epsC,
    epsI = epsI,
    epsS = epsS,

    # covariances
    cov = list(
      F = sep_ar1_struct(n_a, n_t, sigmaF, rhoF_a, rhoF_t),
      C = sep_ar1_struct(n_a, n_t, sigmaC, rhoC_a, rhoC_t),
      I = sep_ar1_struct(n_a, n_t, sigmaI, rhoI_a, rhoI_t),
      S = sep_ar1_struct(n_a, n_t, sigmaS, rhoS_a, rhoS_t)
    ),

    # at-age inputs
    sel_a = sel_a,
    srv_a = srv_a,
    Qa = Qa,

    # derived time series
    R_t = R_t,
    Fbar_t = Fbar_t,
    F_t = F_t,

    # SDs kept at the top level for convenience
    sigmaR = sigmaR,
    sigmaI = sigmaI,
    sigmaC = sigmaC,
    sigmaF = sigmaF,
    sigmaS = sigmaS,
    rec_proc = rec_proc,
    phiR = phiR,
    thetaR = thetaR,
    plus_group = plus_group,
    seed = seed,

    # true parameter vector (estimation scale)
    true = c(
      logM = log(M),
      logR0 = logR0,
      logsigmaR = log(sigmaR),
      logsigmaI = log(sigmaI),
      logsigmaC = log(sigmaC),
      logsigmaF = log(sigmaF),
      logsigmaS = log(sigmaS),
      logq = log(q),
      trans_rhoI_a = atanh(rhoI_a),
      trans_rhoI_t = atanh(rhoI_t),
      trans_rhoC_a = atanh(rhoC_a),
      trans_rhoC_t = atanh(rhoC_t),
      trans_rhoF_a = atanh(rhoF_a),
      trans_rhoF_t = atanh(rhoF_t),
      trans_rhoS_a = atanh(rhoS_a),
      trans_rhoS_t = atanh(rhoS_t),
      trans_phiR   = qlogis((phiR + 1) / 2),
      trans_thetaR = qlogis((thetaR + 1) / 2)
    )
  )

}

#' Bind data to a two-argument objective, giving the single-argument function
#' that RTMB::MakeADFun() expects
cmb <- function(f, d) function(p) f(d, p)

#' Joint negative log likelihood of an age-structured model
#'
#' @param data List: n_a, n_t, ln_C, ln_I, plus, rec_proc (iid|ma1|ar1), rec_re (state|innovation), obs, F_form,
#'   surv (deterministic|aux_plus|aux_all|state), tiny, sparse_lse
#' @param pars Parameter list; only what the configuration uses need be present
jnll_agestrc <- function(data, pars) {

  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")

  RTMB::getAll(data, pars)
  stopifnot(rec_proc %in% c("iid", "ma1", "ar1"),
            rec_re %in% c("state", "innovation"),
            plus %in% c("none", "naive", "logspace"),
            obs %in% c("iid", "sep_ar1"), F_form %in% c("separable", "rw"),
            surv %in% c("deterministic", "aux_plus", "aux_all", "state"))

  # initialize
  ages <- seq_len(n_a)
  years <- seq_len(n_t)
  jnll_agestrc <- 0

  # parameter transforms
  M <- exp(logM)
  Q_a <- exp(logQ_a)
  sigmaR <- exp(logsigmaR)

  # transform more stuff
  sigmaC <- exp(logsigmaC)
  sigmaI <- exp(logsigmaI)

  # survival as a state parameter or using auxillary variable
  if (surv == "state") sigmaS <- exp(logsigmaS)
  if (surv == "aux_all" || surv == "aux_plus") sigmaS <- tiny # pinned

  # fishing mortality form
  if (F_form == "rw") {
    sigmaF <- exp(logsigmaF)
    rho_aF <- tanh(trans_rhoF_a)
  }


  # how to do obsevations
  if (obs == "sep_ar1") {
    rho_aC <- tanh(trans_rhoC_a)
    rho_yC <- tanh(trans_rhoC_t)
    rho_aI <- tanh(trans_rhoI_a)
    rho_yI <- tanh(trans_rhoI_t)
  }

  # separable by age and time
  if (F_form == "separable") {
    ln_Fat <- logF_a + rep(logF_t, each = n_a)
    dim(ln_Fat) <- c(n_a, n_t)
  }

  F_at <- exp(ln_Fat)
  Z_at <- M + F_at

  # recruitment: one explicit block per process x parameterisation
  lnN1 <- numeric(n_t)

  # iid, innovation: standard-normal innovations, lnN1 derived
  if (rec_proc == "iid" && rec_re == "innovation") {
    lnN1 <- logR0 + sigmaR * epsR
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(epsR, 0, 1, TRUE))
  }

  # iid, state: log recruitment is the latent
  if (rec_proc == "iid" && rec_re == "state") {
    lnN1 <- ln_R
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(ln_R, logR0, sigmaR, TRUE))
  }

  # MA(1), innovation: innovations are the latent, lnN1 derived, eps_0 = 0
  if (rec_proc == "ma1" && rec_re == "innovation") {
    theta <- 2 * plogis(trans_thetaR) - 1 # bound -1 and 1
    lnN1[1] <- logR0 + epsR[1]
    for (t in 2:n_t) lnN1[t] <- logR0 + epsR[t] + theta * epsR[t - 1]
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(epsR, 0, sigmaR, TRUE))
  }

  # MA(1), state: log recruitment is the latent; inverting the filter makes eps_t read every earlier one
  if (rec_proc == "ma1" && rec_re == "state") {
    theta <- 2 * plogis(trans_thetaR) - 1 # bound -1 and 1
    lnN1 <- ln_R
    e <- numeric(n_t)
    e[1] <- ln_R[1] - logR0
    for (t in 2:n_t) e[t] <- ln_R[t] - logR0 - theta * e[t - 1]
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(e, 0, sigmaR, TRUE))
  }

  # AR(1), innovation: white noise drives the recursion, lnN1 derived
  if (rec_proc == "ar1" && rec_re == "innovation") {
    phi <- 2 * plogis(trans_phiR) - 1 # bound -1 and 1
    lnN1[1] <- logR0 + sigmaR / sqrt(1 - phi^2) * epsR[1]
    for (t in 2:n_t) lnN1[t] <- phi * lnN1[t - 1] + logR0 * (1 - phi) + sigmaR * epsR[t]
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(epsR, 0, 1, TRUE))
  }

  # AR(1), state: log recruitment is the latent; sigmaR is the innovation sd, dautoreg wants the marginal
  if (rec_proc == "ar1" && rec_re == "state") {
    phi <- 2 * plogis(trans_phiR) - 1 # bound -1 and 1
    lnN1 <- ln_R
    jnll_agestrc <- jnll_agestrc - RTMB::dautoreg(ln_R - logR0, mu = 0, phi = phi, TRUE, scale = sigmaR / sqrt(1 - phi^2))
  }

  # setup numbers at age process
  ln_N <- array(0, dim = c(n_a, n_t))
  ln_N[1, ] <- lnN1 # recruitment
  ln_N[2:n_a, 1] <- ln_Ninit # initialize model
  if (surv == "state" || surv == "aux_all") ln_N[2:n_a, 2:n_t] <- ln_N2 # all ages 2+ as latents
  if (surv == "aux_plus") ln_N[n_a, 2:n_t] <- aux_t # plus group only as a latent

  ln_Npred <- array(0, dim = c(n_a, n_t))
  for (t in 2:n_t) {

    # not the plus group
    for (a in 2:(n_a - 1)) ln_Npred[a, t] <- ln_N[a - 1, t - 1] - Z_at[a - 1, t - 1]

    # plus group: none (oldest age is a true age), naive log-sum-exp, or the log-space trick
    ln_plus <- c(ln_N[n_a - 1, t - 1] - Z_at[n_a - 1, t - 1], ln_N[n_a, t - 1] - Z_at[n_a, t - 1])
    if (plus == "none") ln_Npred[n_a, t] <- ln_plus[1] # just delete the plus group / ignore
    if (plus == "naive") ln_Npred[n_a, t] <- log(sum(exp(ln_plus))) # dense version
    if (plus == "logspace") ln_Npred[n_a, t] <- logspace_sum(ln_plus, sparse_lse) # lse mode

    # deterministic: the state is the prediction, needed by the next year
    if(surv == "deterministic") ln_N[2:n_a, t] <- ln_Npred[2:n_a, t]
    if(surv == "aux_plus") ln_N[2:(n_a - 1), t] <- ln_Npred[2:(n_a - 1), t] # ages below the plus group stay deterministic
  }

  # predictions
  ln_Cpred <- log(F_at / Z_at * (1 - exp(-Z_at))) + ln_N
  ln_Ipred <- logQ_a + ln_N

  ## processes ----------------------------------------------------------------
  # all ages 2+ as latents: sd estimated (state) or pinned (aux_all)
  if (surv == "state" || surv == "aux_all") {
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(ln_N[2:n_a, 2:n_t], ln_Npred[2:n_a, 2:n_t], sigmaS, TRUE))
  }

  # plus group only as a latent, pinned
  if (surv == "aux_plus") {
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(ln_N[n_a, 2:n_t], ln_Npred[n_a, 2:n_t], sigmaS, TRUE))
  }

  # log F random walk, AR1 across ages
  if (F_form == "rw") {
    dF <- ln_Fat[, 2:n_t] - ln_Fat[, 1:(n_t - 1)]
    dim(dF) <- c(n_a, n_t - 1)
    f1 <- function(x) RTMB::dautoreg(x, mu = 0, phi = rho_aF, TRUE)
    f2 <- function(x) RTMB::dautoreg(x, mu = 0, phi = 0, TRUE)
    jnll_agestrc <- jnll_agestrc - RTMB::dseparable(f1, f2)(dF, scale = sigmaF)
  }

  ## observations -------------------------------------------------------------

  # iid
  if (obs == "iid") {
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(ln_C - ln_Cpred, 0, sigmaC, TRUE))
    jnll_agestrc <- jnll_agestrc - sum(RTMB::dnorm(ln_I - ln_Ipred, 0, sigmaI, TRUE))
  }

  # separable AR(1) age x AR(1) year
  if (obs == "sep_ar1") {
    f1 <- function(x) RTMB::dautoreg(x, mu = 0, phi = rho_aC, TRUE)
    f2 <- function(x) RTMB::dautoreg(x, mu = 0, phi = rho_yC, TRUE)
    jnll_agestrc <- jnll_agestrc - RTMB::dseparable(f1, f2)(ln_C - ln_Cpred, scale = sigmaC)

    f1 <- function(x) RTMB::dautoreg(x, mu = 0, phi = rho_aI, TRUE)
    f2 <- function(x) RTMB::dautoreg(x, mu = 0, phi = rho_yI, TRUE)
    jnll_agestrc <- jnll_agestrc - RTMB::dseparable(f1, f2)(ln_I - ln_Ipred, scale = sigmaI)
  }

  # reporting
  RTMB::REPORT(ln_N)
  RTMB::REPORT(ln_Npred)
  RTMB::REPORT(ln_Cpred)
  RTMB::REPORT(ln_Ipred)
  RTMB::REPORT(F_at)

  return(jnll_agestrc)
}
