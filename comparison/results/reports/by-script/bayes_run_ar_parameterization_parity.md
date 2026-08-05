# Diagnostics for `comparison/bayes/run_ar_parameterization_parity.R`

This page maps the executable comparison script to its Quadra diagnostic reports and the legacy RESparsity scripts it represents.

| Case | Quadra diagnostics | Legacy/source script | Status |
|---|---|---|---|
| centered | [bayes_ar_centered.md](../bayes_ar_centered.md) | `bayes/ar_good_comparison.R` | `passing` |
| centered manual | [bayes_ar_centered_manual.md](../bayes_ar_centered_manual.md) | `bayes/ar_good_nodautoreg_comparison.R` | `passing` |
| scaled innovations | [bayes_ar_scaled_innovations.md](../bayes_ar_scaled_innovations.md) | `bayes/ar_good_eps_comparison.R` | `passing` |
| standard innovations | [bayes_ar_standard_innovations.md](../bayes_ar_standard_innovations.md) | `bayes/ar_bad_comparison.R` | `passing` |

## What the reports contain

Each linked report records Quadra's objective/gradient state, curvature health, complete signed eigenvalue spectrum, Hessian inertia, effective structure, and backend-factorization selection. Latent models additionally include uncertainty and latent-state diagnostics.
