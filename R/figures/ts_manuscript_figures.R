library(dplyr)
library(ggplot2)

files <- c("AR(1)" = "R/simple_benchmark/ar1_benchmark.rds",
           "MA(1)" = "R/simple_benchmark/ma1_benchmark.rds")

facet_names <- c("obj" = "MakeADFun", "opt" = "nlminb", "Total" = "Total")
re_labels <- c("process" = "states", "deviations" = "innovations")

# set-up dataframe
df <- bind_rows(lapply(names(files), function(p) {
  results <- readRDS(files[[p]])
  bind_rows(
    data.frame(n = results$timings$n,
               method = results$timings$method,
               func = results$timings$func,
               median_time = results$timings$median_time),
    data.frame(n = results$totals$n,
               method = results$totals$method,
               func = "Total",
               median_time = results$totals$total_time)
  ) |>
    mutate(process = p)
}))

df$func <- factor(facet_names[df$func], levels = facet_names)
df$method <- factor(re_labels[df$method], levels = re_labels)
df$process <- factor(df$process, levels = names(files))

# plot results
png("figures/ts_time.png", width = 8, height = 6, units = "in", res = 300)
df |>
  ggplot(mapping = aes(x = n, y = median_time, color = method)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  facet_grid(process ~ func) +
  labs(x = "Number of time steps (n) - log scale",
       y = "Median time (seconds) - log scale",
       color = "Random effect") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
dev.off()
