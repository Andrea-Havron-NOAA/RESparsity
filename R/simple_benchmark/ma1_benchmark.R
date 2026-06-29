library(RTMB)
library(bench)

# Randomly draw seed using Sys.time()
# random_seed <- as.integer(Sys.time()) %% 100000
# print(random_seed)
random_seed <- 38139

# Setup models for process and deviation parameterizations

ma1_process <- function(par) {
  getAll(dat, par)
  phi <- 2 * plogis(tPhi) - 1
  time_steps <- length(x)
  sd <- exp(logSigma)
  sdo <- exp(logSigmaObs)
  nll <- 0

  eps <- numeric(time_steps)
  eps[1] <- eps0
  for (i in 1:time_steps) {
    eps[i + 1] <- x[i] - phi * eps[i]
  }

  nll <- nll - sum(dnorm(eps, 0, sd, log = TRUE))
  nll <- nll - sum(dnorm(y, x, sdo, TRUE))

  nll
}

ma1_deviation <- function(par) {
  getAll(dat, par)
  phi <- 2 * plogis(tPhi) - 1
  time_steps <- length(x) - 1
  sd <- exp(logSigma)
  sdo <- exp(logSigmaObs)
  nll <- 0

  nll <- nll - sum(dnorm(x, 0, sd, log = TRUE))
  lam <- numeric(length(x) - 1)
  for(i in 1:time_steps){
    lam[i] <- x[i + 1] + phi * x[i]
  }
  nll <- nll - sum(dnorm(y, lam, sdo, TRUE))

  nll
}

#True parameters and initial condition
sd <- 1
sdo <- 1
phi <- 0.7
set.seed(random_seed)

# Warm-up: run models a number of times to remove the effects of C++ linking
# and memory allocation
for (ii in 1:10) {

  #simulate data
  n <- 30
  set.seed(ii)

  # simulate ma1 process
  eps <- rnorm(n + 1, sd = sd)
  x <- rep(0, n)
  for (i in 1:n) {
    x[i] <- eps[i + 1] + phi * eps[i]
  }

  # simulate observations
  y <- x + rnorm(length(x), sd = sdo)
  dat <- data.frame(y = y)

  # dat needs to be globally accessible
  dat <<- dat
  par_process <- list(logSigma = 0, tPhi = 0, eps0 = 0,
                      logSigmaObs = 0, x = rep(0, length(x)))
  par_deviation <- list(logSigma = 0, tPhi = 0,
                        logSigmaObs = 0, x = rep(0, length(x) + 1))

  # run process model
  obj_process <- MakeADFun(ma1_process, par_process,
                           random = c("x", "eps0"), silent = TRUE)
  opt_process <- nlminb(obj_process$par, obj_process$fn, obj_process$gr,
                        control = list(iter.max = 1000, eval.max = 1000))

  # run deviation model
  obj_deviation <- MakeADFun(ma1_deviation, par_deviation,
                             random = "x", silent = TRUE)
  opt_deviation <- nlminb(obj_deviation$par, obj_deviation$fn,
                          obj_deviation$gr,
                          control = list(iter.max = 1000, eval.max = 1000))
}


set.seed(random_seed)
# state vector
n_sim <- 500
x_sim <- rep(0, n_sim)


# simulate ma1 process
eps <- rnorm(n_sim + 1, sd = sd)
x <- rep(0, n_sim)
for (i in 1:n_sim) {
  x[i] <- eps[i + 1] + phi * eps[i]
}

# simulate observations
y <- x + rnorm(length(x), sd = sdo)

results <- bench::press(
  n = c(30, 50, 100, 200, 300, 400, 500),
  {
    gc(reset = TRUE)
    dat <- data.frame(y = y[1:n])

    # dat needs to be globally accessible
    dat <<- dat
    par_process <- list(logSigma = 0, tPhi = 0, eps0 = 0,
                        logSigmaObs = 0, x = rep(0, length(x)))
    par_deviation <- list(logSigma = 0, tPhi = 0,
                          logSigmaObs = 0, x = rep(0, length(x) + 1))

    results <- bench::mark(

      process_obj = {
        # run process model
        obj_process <- MakeADFun(ma1_process, par_process,
                                 random = "x", silent = TRUE)
      },

      process_opt = {
        opt_process <- nlminb(obj_process$par, obj_process$fn,
                              obj_process$gr,
                              control = list(iter.max = 1000,
                                             eval.max = 1000))
      },

      deviations_obj = {
        # run deviations model
        obj_deviations <- MakeADFun(ma1_deviation, par_deviation,
                                    random = "x", silent = TRUE)
      },

      deviations_opt = {
        opt_deviations <- nlminb(obj_deviations$par,
                                 obj_deviations$fn,
                                 obj_deviations$gr,
                                 control = list(iter.max = 1000,
                                                eval.max = 1000))
      },
      iterations = 100, time_unit = "s", filter_gc = FALSE, check = FALSE
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

df$method <- rep(c("process", "process", "deviations", "deviations"), 7)
df$func <- rep(c("obj", "opt"), 14)

df_total <- df |>
  dplyr::group_by(n, method) |>
  dplyr::summarize(total_time = sum(median_time),
                   total_memory = sum(memory))
df_total$func = "Total"
facet_names = c("obj" = "MakeADFun", "opt" = "nlminb", "Total" = "Total")

# plot results
library(ggplot2)
png("figures/ma1_benchmark.png", width = 6, height = 4, units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = median_time |> log(), color = method)) +
  geom_line() +
  geom_point() +
  geom_line(data = df_total, mapping = aes(x = n, y = total_time |> log(), color = method)) +
  geom_point(data = df_total, mapping = aes(x = n, y = total_time |> log(), color = method)) +
  labs(x = "Number of time steps (N)", y = "Log of Median time, 
       log(seconds)", color = "Method") +
  theme_minimal() + facet_wrap(~func, labeller = as_labeller(facet_names)) +
  theme(strip.text = element_text(face = "bold"))
dev.off()

png("figures/ma1_memory.png", width = 6, height = 4, units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = memory, color = method)) +
  geom_line() +
  geom_point() +
  geom_line(data = df_total, mapping = aes(x = n, y = total_memory, color = method)) +
  geom_point(data = df_total, mapping = aes(x = n, y = total_memory, color = method)) +
  labs(x = "Number of time steps (N)", y = "Memory allocation", 
       color = "Method") +
  theme_minimal() + facet_wrap(~func, labeller = as_labeller(facet_names)) +
  theme(strip.text = element_text(face = "bold"))
dev.off()

