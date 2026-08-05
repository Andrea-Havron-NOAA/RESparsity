suppressPackageStartupMessages(library(RTMB))

root <- normalizePath(file.path(getwd()), mustWork = TRUE)
quadra_root <- Sys.getenv("QUADRA_ROOT", file.path(root, "quadra"))
source_file <- file.path(root, "comparison", "ar1", "quadra_ar1_parity.cpp")
binary <- tempfile("quadra-ar1-parity-")

include_flags <- c(
  paste0("-I", quadra_root),
  paste0("-I", file.path(quadra_root, "external", "eigen")),
  paste0("-I", file.path(quadra_root, "external", "had")),
  paste0("-I", file.path(quadra_root, "external", "LBFGSpp", "include"))
)

compile_status <- system2(
  Sys.getenv("CXX", "c++"),
  c("-std=c++17", "-O2", include_flags, "-o", binary, source_file)
)
if (compile_status != 0L) {
  stop("Failed to compile the Quadra AR(1) parity executable.")
}
on.exit(unlink(binary), add = TRUE)

diagnostics_dir <- file.path(root, "comparison", "results", "reports")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
quadra_output <- system2(binary,
                         file.path(diagnostics_dir, "observed_ar1_fixed.md"),
                         stdout = TRUE)
quadra <- read.csv(text = paste(quadra_output, collapse = "\n"))

dat <- list(
  x = c(-0.2, 0.4, 1.1, 0.7),
  mean = 0.3,
  innovation_sd = 0.8
)
parameters <- list(unconstrained_phi = 0.4)

nll <- function(parameters) {
  getAll(dat, parameters)
  phi <- 2 * plogis(unconstrained_phi) - 1
  value <- -dnorm(
    x[1], mean,
    innovation_sd / sqrt(1 - phi^2),
    log = TRUE
  )
  for (i in 2:length(x)) {
    conditional_mean <- mean + phi * (x[i - 1] - mean)
    value <- value - dnorm(x[i], conditional_mean, innovation_sd, log = TRUE)
  }
  value
}

object <- MakeADFun(nll, parameters, silent = TRUE)
rtmb <- data.frame(
  engine = "rtmb",
  objective = object$fn(object$par),
  gradient = unname(object$gr(object$par)[1]),
  phi = 2 * plogis(parameters$unconstrained_phi) - 1
)

results <- rbind(quadra, rtmb)
objective_difference <- abs(diff(results$objective))
gradient_difference <- abs(diff(results$gradient))
phi_difference <- abs(diff(results$phi))

print(results, row.names = FALSE, digits = 16)
cat(sprintf("absolute objective difference: %.3e\n", objective_difference))
cat(sprintf("absolute gradient difference:  %.3e\n", gradient_difference))
cat(sprintf("absolute phi difference:       %.3e\n", phi_difference))

stopifnot(
  objective_difference < 1e-10,
  gradient_difference < 1e-8,
  phi_difference < 1e-14
)

cat("PASS: Quadra and RTMB AR(1) objective and gradient agree.\n")
