# Diagnostics for `comparison/ar1/run_laplace_parity.R`

This page maps the executable comparison script to its Quadra diagnostic reports and the legacy RESparsity scripts it represents.

| Case | Quadra diagnostics | Legacy/source script | Status |
|---|---|---|---|
| latent centered AR(1) | [latent_ar1_laplace.md](../latent_ar1_laplace.md) | `R/ar1_simtest.R` | `passing` |

## What the reports contain

Each linked report records Quadra's objective/gradient state, curvature health, complete signed eigenvalue spectrum, Hessian inertia, effective structure, and backend-factorization selection. Latent models additionally include uncertainty and latent-state diagnostics.
