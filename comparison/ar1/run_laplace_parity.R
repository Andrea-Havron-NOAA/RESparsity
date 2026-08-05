suppressPackageStartupMessages(library(RTMB))

root <- normalizePath(getwd(), mustWork = TRUE)
quadra_root <- Sys.getenv("QUADRA_ROOT", file.path(root, "quadra"))
source_file <- file.path(root, "comparison", "ar1", "quadra_ar1_laplace.cpp")
binary <- tempfile("quadra-ar1-laplace-")
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
if (status != 0L) stop("Failed to compile the Quadra latent AR(1) executable.")
on.exit(unlink(binary), add = TRUE)
diagnostics_dir <- file.path(root, "comparison", "results", "reports")
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
quadra_output <- system2(binary,
                         file.path(diagnostics_dir, "latent_ar1_laplace"),
                         stdout = TRUE)
header <- grep("^engine,", quadra_output)
if (length(header) != 1L || header == length(quadra_output)) {
  stop("Quadra latent AR(1) output did not contain the expected CSV record.")
}
quadra <- read.csv(
  text = paste(quadra_output[c(header, header + 1L)], collapse = "\n")
)

dat <- list(y = c(-0.1, 0.6, 1.0, 0.5), process_sd = 0.8, observation_sd = 0.5)
parameters <- list(unconstrained_phi = 0.4, x = rep(0, length(dat$y)))
nll <- function(parameters) {
  getAll(dat, parameters)
  phi <- 2 * plogis(unconstrained_phi) - 1
  value <- -dnorm(x[1], 0, process_sd / sqrt(1 - phi^2), log = TRUE)
  for (i in 2:length(x)) {
    value <- value - dnorm(x[i], phi * x[i - 1], process_sd, log = TRUE)
  }
  value - sum(dnorm(y, x, observation_sd, log = TRUE))
}

object <- MakeADFun(nll, parameters, random = "x", silent = TRUE)
objective <- object$fn(object$par)
gradient <- unname(object$gr(object$par)[1])
full_mode <- object$env$last.par.best
random_index <- object$env$random
u_hat <- unname(full_mode[random_index])
joint <- object$env$f(full_mode)
hessian <- object$env$spHess(random = TRUE)
logdet <- as.numeric(determinant(as.matrix(hessian), logarithm = TRUE)$modulus)

rtmb <- data.frame(
  engine = "rtmb", objective = objective, joint = joint, logdet = logdet,
  gradient = gradient,
  u0 = u_hat[1], u1 = u_hat[2], u2 = u_hat[3], u3 = u_hat[4]
)
results <- rbind(quadra, rtmb)
metrics <- setdiff(names(results), "engine")
differences <- vapply(metrics, function(name) abs(diff(results[[name]])), numeric(1))

print(results, row.names = FALSE, digits = 16)
print(differences, digits = 4)
stopifnot(
  differences[["objective"]] < 1e-8,
  differences[["joint"]] < 1e-8,
  differences[["logdet"]] < 1e-8,
  differences[["gradient"]] < 1e-8,
  max(differences[c("u0", "u1", "u2", "u3")]) < 1e-8
)
cat("PASS: Quadra and RTMB latent AR(1) Laplace calculations agree.\n")
