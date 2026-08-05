runtime_dir <- file.path("comparison", "results", "runtime-comparisons")
results_file <- file.path(runtime_dir, "ar1_scaling.csv")
results <- read.csv(results_file, stringsAsFactors = FALSE)
results$engine <- factor(results$engine, levels = c("rtmb", "quadra"))
colors <- c(rtmb = "#0072B2", quadra = "#D55E00")

png(
  file.path(runtime_dir, "ar1_scaling.png"),
  width = 1800, height = 1200, res = 180
)
par(mfrow = c(2, 3), mar = c(4.2, 4.6, 2.2, 1.2))

draw_panel <- function(column, title, y_label) {
  valid <- is.finite(results[[column]]) & results[[column]] > 0
  plot(
    NA, xlim = range(results$n[valid]), ylim = range(results[[column]][valid]),
    log = "xy", xlab = "Latent states", ylab = y_label, main = title
  )
  grid()
  for (engine in levels(results$engine)) {
    selected <- results$engine == engine & is.finite(results[[column]]) &
      results[[column]] > 0
    lines(
      results$n[selected], results[[column]][selected], type = "b", pch = 19,
      col = colors[[engine]], lwd = 2
    )
  }
  legend("topleft", levels(results$engine), col = colors, lty = 1, pch = 19, bty = "n")
}

draw_panel("warm_laplace_mean_ms", "Warm Laplace evaluation", "Mean time (ms)")
draw_panel("gradient_ms", "Exact marginal gradient", "Mean time (ms)")
draw_panel("optimization_ms", "Fixed-effect optimization", "Time (ms)")
draw_panel("peak_rss_mib", "Isolated process peak RSS", "Peak RSS (MiB)")
dev.off()

cat("Wrote comparison/results/runtime-comparisons/ar1_scaling.png\n")
