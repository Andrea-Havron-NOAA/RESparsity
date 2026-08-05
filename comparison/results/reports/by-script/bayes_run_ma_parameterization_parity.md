# Diagnostics for `comparison/bayes/run_ma_parameterization_parity.R`

This page maps the executable comparison script to its Quadra diagnostic reports and the legacy RESparsity scripts it represents.

| Case | Quadra diagnostics | Legacy/source script | Status |
|---|---|---|---|
| innovation good | [bayes_ma_innovation_good.md](../bayes_ma_innovation_good.md) | `bayes/ma_good_comparison.R` | `passing` |
| inverse-filter bad | [bayes_ma_inverse_filter_bad.md](../bayes_ma_inverse_filter_bad.md) | `bayes/ma_bad_comparison.R` | `passing` |

## What the reports contain

Each linked report records Quadra's objective/gradient state, curvature health, complete signed eigenvalue spectrum, Hessian inertia, effective structure, and backend-factorization selection. Latent models additionally include uncertainty and latent-state diagnostics.
