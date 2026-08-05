suppressPackageStartupMessages({
  library(RTMB)
  library(Matrix)
})

root <- normalizePath(getwd(), mustWork = TRUE)
quadra_root <- Sys.getenv("QUADRA_ROOT", file.path(root, "quadra"))
source_file <- file.path(root, "comparison", "ar1", "quadra_ar1_scaling.cpp")
binary <- tempfile("quadra-ar1-scaling-")
results_dir <- file.path(root, "comparison", "results", "runtime-comparisons")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

include_flags <- c(
  paste0("-I", quadra_root),
  paste0("-I", file.path(quadra_root, "external", "eigen")),
  paste0("-I", file.path(quadra_root, "external", "had")),
  paste0("-I", file.path(quadra_root, "external", "LBFGSpp", "include"))
)
status <- system2(
  Sys.getenv("CXX", "c++"),
  c("-std=c++17", "-O3", include_flags, "-o", binary, source_file)
)
if (status != 0L) stop("Failed to compile the Quadra scaling benchmark.")
on.exit(unlink(binary), add = TRUE)

run_measured <- function(command, args) {
  stdout_file <- tempfile("comparison-stdout-")
  stderr_file <- tempfile("comparison-time-")
  on.exit(unlink(c(stdout_file, stderr_file)), add = TRUE)
  status <- system2(
    "/usr/bin/time", c("-l", command, args),
    stdout = stdout_file, stderr = stderr_file
  )
  if (status != 0L) {
    stop("Measured command failed: ", command, "\n",
         paste(readLines(stderr_file, warn = FALSE), collapse = "\n"))
  }
  timing <- readLines(stderr_file, warn = FALSE)
  rss_line <- grep("maximum resident set size", timing, value = TRUE)
  if (length(rss_line) != 1L) stop("Peak RSS was not reported for ", command)
  rss_bytes <- as.numeric(strsplit(trimws(rss_line), "[[:space:]]+")[[1]][1])
  list(
    output = readLines(stdout_file, warn = FALSE),
    peak_rss_mib = rss_bytes / 1024^2
  )
}

sizes <- c(30L, 100L, 300L, 1000L)
quadra <- do.call(rbind, lapply(sizes, function(n) {
  measured <- run_measured(binary, as.character(n))
  header <- grep("^engine,n,", measured$output)
  record <- grep("^quadra,", measured$output)
  if (length(header) != 1L || length(record) != 1L)
    stop("Could not locate Quadra benchmark CSV output for n=", n)
  row <- read.csv(text = paste(measured$output[c(header, record)], collapse = "\n"))
  row$peak_rss_mib <- measured$peak_rss_mib
  row
}))

rtmb_worker <- file.path(root, "comparison", "ar1", "rtmb_ar1_scaling_worker.R")
rtmb <- do.call(rbind, lapply(sizes, function(n) {
  measured <- run_measured(
    file.path(R.home("bin"), "Rscript"), c(rtmb_worker, as.character(n))
  )
  row <- read.csv(text = paste(measured$output, collapse = "\n"))
  row$peak_rss_mib <- measured$peak_rss_mib
  row
}))
results <- rbind(quadra, rtmb)
results <- results[order(results$n, results$engine), ]

paired <- merge(
  subset(results, engine == "quadra"),
  subset(results, engine == "rtmb"),
  by = "n", suffixes = c("_quadra", "_rtmb")
)
stopifnot(
  all(results$success == 1L),
  all(paired$hessian_unique_nnz_quadra == paired$hessian_unique_nnz_rtmb),
  all(paired$hessian_bandwidth_quadra == 1L),
  all(paired$hessian_bandwidth_rtmb == 1L),
  max(abs(paired$objective_quadra - paired$objective_rtmb)) < 1e-8,
  max(abs(paired$gradient_quadra - paired$gradient_rtmb)) < 1e-8,
  all(paired$optimization_converged_quadra == 1L),
  all(paired$optimization_converged_rtmb == 1L),
  max(abs(paired$estimate_quadra - paired$estimate_rtmb)) < 1e-5
)

ratios <- data.frame(
  n = paired$n,
  warm_laplace_quadra_over_rtmb =
    paired$warm_laplace_mean_ms_quadra / paired$warm_laplace_mean_ms_rtmb,
  gradient_quadra_over_rtmb =
    paired$gradient_ms_quadra / paired$gradient_ms_rtmb,
  peak_rss_quadra_over_rtmb =
    paired$peak_rss_mib_quadra / paired$peak_rss_mib_rtmb,
  optimization_quadra_over_rtmb =
    paired$optimization_ms_quadra / paired$optimization_ms_rtmb
)

output <- file.path(results_dir, "ar1_scaling.csv")
write.csv(results, output, row.names = FALSE)
write.csv(
  ratios, file.path(results_dir, "ar1_scaling_ratios.csv"), row.names = FALSE
)
print(results, row.names = FALSE, digits = 5)
cat("\nQuadra / RTMB timing ratios\n")
print(ratios, row.names = FALSE, digits = 4)
cat("Wrote ", output, "\n", sep = "")

plot_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  file.path("comparison", "ar1", "plot_scaling.R")
)
if (plot_status != 0L) stop("Failed to plot the AR(1) scaling benchmark.")
