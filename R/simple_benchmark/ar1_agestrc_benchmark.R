## AR(1) recruitment, age-structured model timing tests

source("R/simple_benchmark/agestrc_setup.R")

rec_proc <- "ar1"

# warm-up
agestrc_warmup(rec_proc)

## ---- fit stuff once --------------------------------------------------------------

diag <- agestrc_diagnostics(rec_proc)
print(diag, row.names = FALSE, digits = 4)

## ---- time ------------------------------------------------------------------

full_results <- agestrc_timings(rec_proc)
df <- agestrc_tidy(full_results, rec_proc)

saveRDS(list(timings = df, diagnostics = diag),
        "R/simple_benchmark/ar1_agestrc_benchmark.rds")

## ---- plot ------------------------------------------------------------------

library(ggplot2)
png("figures/ar1_agestrc_benchmark.png", width = 9, height = 6,
    units = "in", res = 300)
print(
  df |>
    ggplot(mapping = aes(x = n_t, y = median_time, color = treatment)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Number of time steps (N) - Log Scale",
         y = "Median time, (seconds) - Log Scale", color = "Treatment") +
    theme_minimal() +
    facet_grid(rec_re ~ func,
               labeller = labeller(func = facet_names, rec_re = rec_re_names)) +
    theme(strip.text = element_text(face = "bold"))
)
dev.off()

png("figures/ar1_agestrc_memory.png", width = 9, height = 6,
    units = "in", res = 300)
print(
  df |>
    ggplot(mapping = aes(x = n_t, y = memory, color = treatment)) +
    geom_line() +
    geom_point() +
    labs(x = "Number of time steps (N)", y = "Memory allocation",
         color = "Treatment") +
    theme_minimal() +
    facet_grid(rec_re ~ func, scales = "free_y",
               labeller = labeller(func = facet_names, rec_re = rec_re_names)) +
    theme(strip.text = element_text(face = "bold"))
)
dev.off()

# the sparsity behind the timings
png("figures/ar1_agestrc_sparsity.png", width = 9, height = 4,
    units = "in", res = 300)
print(
  diag |>
    transform(treatment = factor(treatment, levels = names(treatments)),
              rec_re = factor(rec_re, levels = rec_res)) |>
    ggplot(mapping = aes(x = n_t, y = pct_dense, color = treatment)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(x = "Number of time steps (N) - Log Scale",
         y = "Inner Hessian, % dense", color = "Treatment") +
    theme_minimal() +
    facet_wrap(~rec_re, labeller = as_labeller(rec_re_names)) +
    theme(strip.text = element_text(face = "bold"))
)
dev.off()
