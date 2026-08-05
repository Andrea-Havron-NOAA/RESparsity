scripts <- c(
  file.path("comparison", "ar1", "run_parity.R"),
  file.path("comparison", "ma1", "run_parity.R"),
  file.path("comparison", "ar1", "run_laplace_parity.R"),
  file.path("comparison", "bayes", "run_ar_parameterization_parity.R"),
  file.path("comparison", "bayes", "run_ma_parameterization_parity.R"),
  file.path("comparison", "fisheries", "run_simplefsa_ar_gate.R"),
  file.path("comparison", "fisheries", "run_simplefsa_ar_nuisance_gate.R"),
  file.path("comparison", "fisheries", "run_simplefsa_ar_full_fixed_gate.R"),
  file.path("comparison", "fisheries", "run_simplefsa_ma_good_full_gate.R"),
  file.path("comparison", "fisheries", "run_simplefsa_ma_bad_full_gate.R"),
  file.path("comparison", "fisheries", "audit_naa_sources.R")
)

for (script in scripts) {
  cat("\n== ", script, " ==\n", sep = "")
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (status != 0L) {
    stop("Comparison failed: ", script)
  }
}

parity_reports <- c(
  "latent_ar1_laplace",
  "bayes_ar_centered",
  "bayes_ar_centered_manual",
  "bayes_ar_scaled_innovations",
  "bayes_ar_standard_innovations",
  "bayes_ma_innovation_good",
  "bayes_ma_inverse_filter_bad"
)
for (report in parity_reports) {
  base <- file.path("comparison", "results", "reports", report)
  text <- paste(readLines(paste0(base, ".txt"), warn = FALSE), collapse = "\n")
  csv <- paste(readLines(paste0(base, ".csv"), warn = FALSE), collapse = "\n")
  if (!grepl("marginal_fixed_gradient_norm:", text, fixed = TRUE) ||
      !grepl("latent_mode_converged:", text, fixed = TRUE) ||
      grepl("\nconverged:", text, fixed = TRUE) ||
      !grepl("parity_point,marginal_fixed_gradient_norm", csv, fixed = TRUE) ||
      !grepl("latent_mode,converged", csv, fixed = TRUE) ||
      grepl("optimization,converged", csv, fixed = TRUE)) {
    stop("Ambiguous parity-report convergence semantics: ", report)
  }
}

source(file.path("comparison", "write_diagnostics_index.R"))

cat("\nPASS: all Quadra–RTMB comparison gates passed.\n")
