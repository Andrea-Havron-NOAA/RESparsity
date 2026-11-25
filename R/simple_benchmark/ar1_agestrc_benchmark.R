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
