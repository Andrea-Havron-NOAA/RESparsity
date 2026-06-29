library(RTMB)
library(bench)

# Randomly draw seed using Sys.time()
# random_seed <- as.integer(Sys.time()) %% 100000
# print(random_seed)
random_seed <- 38139

# Setup models for process and deviations parameterizations

ar1_process_dautoreg <- function(par) {
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

ar1_deviations <- function(par) {
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

#True parameters and initial condition
sd <- 1
sdo <- 1
phi <- 0.7
set.seed(random_seed)
x_init <- rnorm(1, sd = sqrt(sd * sd / (1 - phi * phi)))

# Warm-up: run models a number of times to remove the effects of C++ linking
# and memory allocation
for (ii in 1:10) {

  #simulate data
  n <- 30
  set.seed(ii)

  # state vector
  x <- rep(0, n)
  x[1] <- x_init

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

  # run process model
  obj_process_dautoreg <- MakeADFun(ar1_process_dautoreg, par, random = "x",
                                   silent = TRUE)
  opt_process_dautoreg <- nlminb(obj_process_dautoreg$par,
                                obj_process_dautoreg$fn,
                                obj_process_dautoreg$gr,
                                control = list(iter.max = 1000,
                                               eval.max = 1000))

  # run deviations model
  obj_deviations <- MakeADFun(ar1_deviations, par, random = "x", silent = TRUE)
  opt_deviations <- nlminb(obj_deviations$par, obj_deviations$fn,
                           obj_deviations$gr, 
                           control = list(iter.max = 1000, eval.max = 1000))
}

set.seed(random_seed)
# state vector
n_sim <- 500
x_sim <- rep(0, n_sim)
x_sim[1] <- x_init

# simulate ar1 process
for (i in 2:n_sim) {
  x_sim[i] <- phi * x_sim[i - 1] + rnorm(1, sd = sd)
}

# simulate observations
y <- x_sim + rnorm(n_sim, sd = sdo)

results <- bench::press(
  n = c(30, 50, 100, 200, 300, 400, 500),
  {
    gc(reset = TRUE)
    dat <- data.frame(y = y[1:n])
    # dat needs to be globally accessible
    dat <<- dat
    par <- list(logSigma = 0, tPhi = 0,
                logSigmaObs = 0, x = rep(0, n))


    results <- bench::mark(
      process_obj = {
        # run process model
        obj_process_dautoreg <- MakeADFun(ar1_process_dautoreg, par,
                                          random = "x", silent = TRUE)
      },

      process_opt = {
        opt_process_dautoreg <- nlminb(obj_process_dautoreg$par,
                                       obj_process_dautoreg$fn,
                                       obj_process_dautoreg$gr,
                                       control = list(iter.max = 1000,
                                                      eval.max = 1000))
        opt_process_dautoreg$convergence
      },

      deviations_obj = {
        # run deviations model
        obj_deviations <- MakeADFun(ar1_deviations, par,
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
expr <-  results$expression |> attr(which = "description")
# set-up dataframe
df <- data.frame(
  n = results$n,
  expression = expr,
  median_time = results$median,
  memory = results$mem_alloc
)

df$method <- rep(c("process", "process", "deviations", "deviations"), 7)
df$func <- rep(c("obj", "opt"), 14)

df_total <- df |>
  dplyr::group_by(n, method) |>
  dplyr::summarize(total_time = sum(median_time),
                   total_memory = sum(memory))
df_total$func <- "Total"
facet_names <- c("obj" = "MakeADFun", "opt" = "nlminb", "Total" = "Total")

# plot results
library(ggplot2)
png("figures/ar1_benchmark.png", width = 6, height = 4, units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = median_time, color = method)) +
  geom_line() +
  geom_point() +
  geom_line(data = df_total,
            mapping = aes(x = n, y = total_time, color = method)) +
  geom_point(data = df_total,
             mapping = aes(x = n, y = total_time, color = method)) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Number of time steps (N) - Log Scale", y = "Median time, 
       (seconds) - Log Scale", color = "Method") +
  theme_minimal() + facet_wrap(~func, labeller = as_labeller(facet_names)) +
  theme(strip.text = element_text(face = "bold"))
dev.off()

png("figures/ar1_memory.png", width = 6, height = 4, units = "in", res = 300)
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
