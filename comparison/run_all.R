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

source(file.path("comparison", "write_diagnostics_index.R"))

cat("\nPASS: all Quadra–RTMB comparison gates passed.\n")
