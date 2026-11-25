source("R/simple_benchmark/agestrc_setup.R")

args <- commandArgs(trailingOnly = TRUE)
rec_proc <- args[1]
stopifnot(rec_proc %in% c("ar1", "ma1"))

trts <- c("none", "naive", "aux_plus", "aux_all")
out <- sprintf("R/simple_benchmark/agestrc_replicates_%s.rds", rec_proc)

sims <- local({
  nm <- if (rec_proc == "ar1") "simAR" else "simMA"
  load(sprintf("simplefsa/%s.RData", nm))
  get(nm)
})
n_rep <- if (length(args) > 1) as.integer(args[2]) else length(sims)

#' One replicate as matrices (ln_I is NA where there is no survey)
replicate_mats <- function(d) {
  n_t <- ncol(d$M)
  a <- d$age - min(d$age) + 1
  y <- d$year - min(d$year) + 1
  ln_C <- matrix(NA_real_, n_a, n_t)
  ln_I <- matrix(NA_real_, length(srv_ages), n_t)
  k <- d$fleet == 1
  ln_C[cbind(a[k], y[k])] <- log(d$obs[k])
  k <- d$fleet == 2
  ln_I[cbind(match(a[k], srv_ages), y[k])] <- log(d$obs[k])
  list(ln_C = ln_C, ln_I = ln_I, n_t = n_t)
}

# starting values as the series-length analysis (agestrc_config())
truth <- agestrc_fsa_truth(rec_proc)
start <- list(logR0 = truth$logR0, ln_Ninit = truth$ln_Ninit,
              logQ_a = truth$logQ_a, logF_t = log(0.4),
              logsigma = log(0.3), trans = qlogis(0.75))

replicate_config <- function(m, treatment, rec_re) {
  trt <- treatments[[treatment]]
  agestrc_build(m$ln_C, m$ln_I, m$n_t, plus = trt$plus, surv = trt$surv,
                sparse_lse = trt$sparse_lse, rec_re = rec_re,
                rec_proc = rec_proc, start = start)
}

#' Fit one configuration from a cold start; timings exclude the diagnostics
replicate_fit <- function(cfg) {
  gc()
  t0 <- bench::hires_time()
  obj <- agestrc_obj(cfg)
  t1 <- bench::hires_time()
  fit <- nlminb(obj$par, obj$fn, obj$gr, control = ctrl)
  t2 <- bench::hires_time()

  H <- obj$env$spHess(obj$env$last.par.best, random = TRUE)
  pd <- tryCatch(sdreport(obj)$pdHess, error = function(e) NA)
  data.frame(obj_time = t1 - t0, opt_time = t2 - t1, total_time = t2 - t0,
             n_re = nrow(H), pct_dense = 100 * nnzero(H) / nrow(H)^2,
             nll = fit$objective, max_gr = max(abs(obj$gr(fit$par))),
             conv = fit$convergence, msg = fit$message, pd_hess = pd)
}

grid <- expand.grid(treatment = trts, rec_re = rec_res, stringsAsFactors = FALSE)

# warm-up: every configuration once
m1 <- replicate_mats(sims[[1]])
for (k in seq_len(nrow(grid)))
  invisible(replicate_fit(replicate_config(m1, grid$treatment[k], grid$rec_re[k])))

set.seed(random_seed) # order of the configurations within a data set
acc <- list()
for (i in seq_len(n_rep)) {
  m <- replicate_mats(sims[[i]])
  for (k in sample(nrow(grid))) {
    res <- replicate_fit(replicate_config(m, grid$treatment[k], grid$rec_re[k]))
    acc[[length(acc) + 1]] <- cbind(rec_proc = rec_proc, i = i, grid[k, ], res,
                                    row.names = NULL)
  }
  saveRDS(do.call(rbind, acc), out)
  message(sprintf("%s data set %d of %d", rec_proc, i, n_rep))
}

sessionInfo()
