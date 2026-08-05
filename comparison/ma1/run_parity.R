suppressPackageStartupMessages(library(RTMB))

root <- normalizePath(getwd(), mustWork = TRUE)
quadra_root <- Sys.getenv("QUADRA_ROOT", file.path(root, "quadra"))
source_file <- file.path(root, "comparison", "ma1", "quadra_ma1_parity.cpp")
binary <- tempfile("quadra-ma1-parity-")

include_flags <- c(
  paste0("-I", quadra_root),
  paste0("-I", file.path(quadra_root, "external", "eigen")),
  paste0("-I", file.path(quadra_root, "external", "had")),
  paste0("-I", file.path(quadra_root, "external", "LBFGSpp", "include"))
)
status <- system2(
  Sys.getenv("CXX", "c++"),
  c("-std=c++17", "-O2", include_flags, "-o", binary, source_file)
)
if (status != 0L) stop("Failed to compile the Quadra MA(1) parity executable.")
on.exit(unlink(binary), add = TRUE)

diagnostics_dir <- file.path(root, "comparison", "results", "reports")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
quadra <- read.csv(text = paste(system2(
  binary, file.path(diagnostics_dir, "observed_ma1_fixed.md"), stdout = TRUE
), collapse = "\n"))

dat <- list(
  innovations = c(-0.3, 0.2, 0.8, -0.1, 0.4),
  y = c(0.7, 1.0, 0.1, 0.9),
  mean = 0.3,
  innovation_sd = 0.8,
  observation_sd = 0.5
)

nll <- function(parameters) {
  getAll(dat, parameters)
  theta <- 2 * plogis(unconstrained_theta) - 1
  process <- numeric(length(innovations) - 1)
  for (i in seq_along(process)) {
    process[i] <- mean + innovations[i + 1] + theta * innovations[i]
  }
  -sum(dnorm(innovations, 0, innovation_sd, log = TRUE)) -
    sum(dnorm(y, process, observation_sd, log = TRUE))
}

parameters <- list(unconstrained_theta = -0.6)
object <- MakeADFun(nll, parameters, silent = TRUE)
rtmb <- data.frame(
  engine = "rtmb",
  objective = object$fn(object$par),
  gradient = unname(object$gr(object$par)[1]),
  theta = 2 * plogis(parameters$unconstrained_theta) - 1
)

results <- rbind(quadra, rtmb)
differences <- c(
  objective = abs(diff(results$objective)),
  gradient = abs(diff(results$gradient)),
  theta = abs(diff(results$theta))
)

print(results, row.names = FALSE, digits = 16)
print(differences, digits = 4)
stopifnot(
  differences[["objective"]] < 1e-10,
  differences[["gradient"]] < 1e-8,
  differences[["theta"]] < 1e-14
)
cat("PASS: Quadra and RTMB MA(1) objective and gradient agree.\n")
