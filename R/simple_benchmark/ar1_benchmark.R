library(RTMB)
library(bench)

# Setup models for sparse and manual parameterizations

ar1_sparse_dautoreg <- function(par) {
  getAll(dat, par)
  phi <- 2 * plogis(tPhi) - 1
  time_steps <- length(x)
  sd <- exp(logSigma)
  sdo <- exp(logSigmaObs)
  nll <- 0

  nll <- nll - dautoreg(x, phi = phi, scale = sqrt(sd * sd / (1 - phi * phi)),
                        log = TRUE)
  nll <- nll - sum(dnorm(y, x, sdo, TRUE))

  nll
}

ar1_manual <- function(par) {
  getAll(dat, par)
  phi <- 2 * plogis(tPhi) - 1
  time_steps <- length(x)
  sd <- exp(logSigma)
  sdo <- exp(logSigmaObs)
  nll <- 0

  nll <- nll - dnorm(x[1], 0, sqrt(sd * sd / (1 - phi * phi)), log = TRUE)
  nll <- nll - sum(dnorm(x[-1], 0, sd, log = TRUE))
  lam <- numeric(length(x))
  lam[1] <- 0 + x[1]
  for(i in 2:time_steps){
    lam[i] <- phi * lam[i - 1] + x[i]
  }
  nll <- nll - sum(dnorm(y, lam, sdo, TRUE))

  nll
}

# Warm-up: run models a number of times to remove the effects of C++ linking
# and memory allocation
for (ii in 1:10) {

  #simulate data
  n <- 30
  set.seed(ii)
  sd <- 1
  sdo <- 1
  phi <- 0.7

  # init condition
  x <- rep(0, n)
  x[1] <- rnorm(1, sd = sqrt(sd * sd / (1 - phi * phi)))

  # simulate ar1 process
  for (i in 2:n) {
    x[i] <- phi * x[i - 1] + rnorm(1, sd = sd)
  }

  # simulate observations
  y <- x + rnorm(length(x), sd = sdo)
  dat <- data.frame(y = y)

  # dat needs to be globally accessible
  dat <<- dat
  par <- list(logSigma = 0, tPhi = 0,
              logSigmaObs = 0, x = rep(0, length(x)))

  # run sparse model
  obj_sparse_dautoreg <- MakeADFun(ar1_sparse_dautoreg, par, random = "x",
                                   silent = TRUE)
  opt_sparse_dautoreg <- nlminb(obj_sparse_dautoreg$par,
                                obj_sparse_dautoreg$fn,
                                obj_sparse_dautoreg$gr,
                                control = list(iter.max = 1000,
                                               eval.max = 1000))

  # run manual model
  obj_manual <- MakeADFun(ar1_manual, par, random = "x", silent = TRUE)
  opt_manual <- nlminb(obj_manual$par, obj_manual$fn, obj_manual$gr,
                       control = list(iter.max = 1000, eval.max = 1000))
}


results <- bench::press(
  n = c(30, 50, 100, 200, 300, 400, 500),
  {
    set.seed(n)
    sd <- 1
    sdo <- 1
    phi <- 0.7

    # init condition
    x <- rep(0, n)
    x[1] <- rnorm(1, sd = sqrt(sd * sd / (1 - phi * phi)))

    # simulate ar1 process
    for (i in 2:n) {
      x[i] <- phi * x[i - 1] + rnorm(1, sd = sd)
    }

    # simulate observations
    y <- x + rnorm(length(x), sd = sdo)
    dat <- data.frame(y = y)

    # dat needs to be globally accessible
    dat <<- dat
    par <- list(logSigma = 0, tPhi = 0,
                logSigmaObs = 0, x = rep(0, length(x)))

    results <- bench::mark(
      sparse = {
        # run sparse model
        obj_sparse_dautoreg <- MakeADFun(ar1_sparse_dautoreg, par,
                                         random = "x", silent = TRUE)
        opt_sparse_dautoreg <- nlminb(obj_sparse_dautoreg$par,
                                      obj_sparse_dautoreg$fn,
                                      obj_sparse_dautoreg$gr,
                                      control=list(iter.max = 1000,
                                                   eval.max = 1000))
      },
      manual = {
        # run manual model
        obj_manual <- MakeADFun(ar1_manual, par, random = "x", silent = TRUE)
        opt_manual <- nlminb(obj_manual$par, obj_manual$fn, obj_manual$gr,
                             control = list(iter.max = 1000, eval.max = 1000))
      },
      iterations = 100, time_unit = "s", filter_gc = FALSE
    )
  }
)

# set-up dataframe
df <- data.frame(
  n = results$n,
  expression = rep(c("process", "deviations"), 7),
  median_time = results$median,
  memory = results$mem_alloc
)

# plot results
library(ggplot2)
png("figures/ar1_benchmark.png", width = 6, height = 4, units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = median_time |> log(), color = expression)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of time steps (N)", y = "Log of Median time, 
       log(seconds)", color = "Method") +
  theme_minimal()
dev.off()

png("figures/ar1_memory.png", width = 6, height = 4, units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = memory, color = expression)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of time steps (N)", y = "Memory allocation", 
       color = "Method") +
  theme_minimal()
dev.off()
