## Manuscript figures for the age-structured benchmark (Experiment 2).
##
## Re-plots R/simple_benchmark/{ar1,ma1}_agestrc_benchmark.rds with AR(1) and
## MA(1) on one panel grid and without the logspace (LSE0) treatment, which is
## not part of the manuscript. Nothing is re-run here; the rds files are read
## as saved by the benchmark scripts.

library(dplyr)
library(ggplot2)

files <- c("AR(1)" = "R/simple_benchmark/ar1_agestrc_benchmark.rds",
           "MA(1)" = "R/simple_benchmark/ma1_agestrc_benchmark.rds")

trt_labels <- c("none" = "No plus group",
                "naive" = "Plus group",
                "aux_plus" = "Aux: plus group",
                "aux_all" = "Aux: ages 2+")
re_labels <- c("state" = "Random effect: states",
               "innovation" = "Random effect: innovations")

read_part <- function(part) {
  bind_rows(lapply(names(files), function(p) {
    readRDS(files[[p]])[[part]] |> mutate(process = p)
  })) |>
    filter(treatment != "logspace (LSE0)") |>
    mutate(treatment = factor(trt_labels[as.character(treatment)], levels = trt_labels),
           rec_re = factor(re_labels[as.character(rec_re)], levels = re_labels),
           process = factor(process, levels = names(files)))
}

timings <- read_part("timings")
diagnostics <- read_part("diagnostics")

## ---- total time ------------------------------------------------------------
p_time <- timings |>
  filter(func == "Total") |>
  ggplot(aes(x = n_t, y = median_time, colour = treatment)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(process ~ rec_re) +
  labs(x = "Number of years - log scale",
       y = "Median total time (seconds) - log scale",
       colour = "Plus group treatment") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

png("figures/agestrc_time.png", width = 8, height = 6.5, units = "in", res = 300)
print(p_time)
dev.off()

## ---- sparsity --------------------------------------------------------------
p_sparse <- diagnostics |>
  ggplot(aes(x = n_t, y = pct_dense, colour = treatment)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  facet_grid(process ~ rec_re) +
  labs(x = "Number of years - log scale",
       y = "Inner Hessian, % nonzero",
       colour = "Plus group treatment") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")

png("figures/agestrc_sparsity.png", width = 8, height = 6.5, units = "in", res = 300)
print(p_sparse)
dev.off()

## ---- numbers quoted in the text -------------------------------------------
cat("\nTotal time at 200 years (s):\n")
timings |>
  filter(func == "Total", n_t == 200) |>
  select(process, rec_re, treatment, median_time) |>
  arrange(process, rec_re, treatment) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat("\nInner Hessian % nonzero at 200 years:\n")
diagnostics |>
  filter(n_t == 200) |>
  select(process, rec_re, treatment, n_re, pct_dense, conv) |>
  arrange(process, rec_re, treatment) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)
