library(RTMB)
library(bench)

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

#True parameters
sd <- 1
sdo <- 1
phi <- 0.7


# Simulate one large dataset
n_sim = 500
set.seed(random_seed)

# simulate ma1 process
eps <- rnorm(n_sim + 1, sd = sd)
x <- rep(0, n_sim)
for (i in 1:n_sim) {
  x[i] <- eps[i + 1] + phi * eps[i]
}

# simulate observations
y <- x + rnorm(length(x), sd = sdo)

# Warm-up: run models a number of times to remove the effects of C++ linking
# and memory allocation
n <- 30
for (ii in 1:10) {

  dat <- data.frame(y = y[1:n])

  # dat needs to be globally accessible
  dat <<- dat
  par_process <- list(logSigma = 0, tPhi = 0, eps0 = 0,
                      logSigmaObs = 0, x = rep(0, n))
  par_deviation <- list(logSigma = 0, tPhi = 0,
                        logSigmaObs = 0, x = rep(0, n + 1))

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


full_results <- bench::press(
  n = c(30, 50, 100, 200, 300, 400, 500),
  {
    gc(reset = TRUE)
    dat <- data.frame(y = y[1:n])

    # dat needs to be globally accessible
    dat <<- dat
    par_process <- list(logSigma = 0, tPhi = 0, eps0 = 0,
                        logSigmaObs = 0, x = rep(0, n))
    par_deviation <- list(logSigma = 0, tPhi = 0,
                          logSigmaObs = 0, x = rep(0, n + 1))

    results <- bench::mark(

      process_obj = {
        # run process model
        obj_process <- MakeADFun(ma1_process, par_process,
                                 random = c("x", "eps0"), silent = TRUE)
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
      iterations = 100, time_unit = "s", check = FALSE
    )
  }
)

# set-up dataframe
df <- full_results |>
  dplyr::mutate(expression = as.character(expression)) |> 
  # 2. Split expression into 'method' and 'func' columns by the underscore
  tidyr::separate_wider_delim(
    cols = expression, 
    delim = "_", 
    names = c("method", "func")
  ) |>
  dplyr::select(!where(is.list))

df_total <- df |>
  dplyr::group_by(n, method) |>
  dplyr::summarize(total_time = sum(median),
                   total_memory = sum(mem_alloc))
df_total$func <- "Total"
facet_names <- c("obj" = "MakeADFun", "opt" = "nlminb", "Total" = "Total")


# plot results
library(ggplot2)
png("figures/ma1_benchmark.png", width = 6, height = 4,
    units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = median |> log(), color = method)) +
  geom_line() +
  geom_point() +
  geom_line(data = df_total,
            mapping = aes(x = n, y = total_time |> log(), color = method)) +
  geom_point(data = df_total,
             mapping = aes(x = n, y = total_time |> log(), color = method)) +
  labs(x = "Number of time steps (N)", y = "Log of Median time, 
       log(seconds)", color = "Method") +
  theme_minimal() + facet_wrap(~func, labeller = as_labeller(facet_names)) +
  theme(strip.text = element_text(face = "bold"))
dev.off()

png("figures/ma1_memory.png", width = 6, height = 4,
    units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = mem_alloc, color = method)) +
  geom_line() +
  geom_point() +
  geom_line(data = df_total,
            mapping = aes(x = n, y = total_memory, color = method)) +
  geom_point(data = df_total,
             mapping = aes(x = n, y = total_memory, color = method)) +
  labs(x = "Number of time steps (N)", y = "Memory allocation",
       color = "Method") +
  theme_minimal() + facet_wrap(~func, labeller = as_labeller(facet_names)) +
  theme(strip.text = element_text(face = "bold"))
dev.off()
