#!/usr/bin/env Rscript

runtime_dir <- file.path("comparison", "results", "runtime-comparisons")
dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)

studies <- data.frame(
  study = c("ar_bad", "ar_good", "ar_good_eps", "ar_good_nodautoreg",
            "ma_bad", "ma_good"),
  process = c("AR1", "AR1", "AR1", "AR1", "MA1", "MA1"),
  latent_parameterization = c(
    "noncentered_standard_innovations", "centered_latent_states_dautoreg",
    "centered_scaled_innovations", "centered_latent_states_manual",
    "centered_latent_states_inverse_filter", "centered_innovations"
  ),
  stringsAsFactors = FALSE
)

summaries <- lapply(seq_len(nrow(studies)), function(index) {
  metadata <- studies[index, ]
  path <- file.path("bayes", paste0(metadata$study, "_results.rds"))
  archived <- readRDS(path)
  do.call(rbind, lapply(split(archived, archived$method), function(rows) {
    data.frame(
      study = metadata$study,
      process = metadata$process,
      latent_parameterization = metadata$latent_parameterization,
      method = rows$method[[1]],
      replicates = nrow(rows),
      timed_replicates = sum(is.finite(rows$time_s)),
      median_time_s = if (any(is.finite(rows$time_s)))
        median(rows$time_s[is.finite(rows$time_s)]) else NA_real_,
      median_max_rhat = if (any(is.finite(rows$max_rhat)))
        median(rows$max_rhat[is.finite(rows$max_rhat)]) else NA_real_,
      median_mean_neff = if (any(is.finite(rows$mean_neff)))
        median(rows$mean_neff[is.finite(rows$mean_neff)]) else NA_real_,
      median_min_neff = if (any(is.finite(rows$min_neff)))
        median(rows$min_neff[is.finite(rows$min_neff)]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
})

result <- do.call(rbind, summaries)
result <- result[order(result$process, result$study, result$method), ]
output <- file.path(runtime_dir, "bayes_archived_summary.csv")
write.csv(result, output, row.names = FALSE)
print(result, row.names = FALSE, digits = 5)
cat("Wrote ", output, "\n", sep = "")
