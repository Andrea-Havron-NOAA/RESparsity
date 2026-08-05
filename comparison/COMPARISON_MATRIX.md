# RESparsity comparison matrix

Comparisons advance only after objective/gradient parity passes for identical
data, parameter maps, transforms, constants, and latent parameterization.
Script-to-report provenance is indexed in
`comparison/results/reports/README.md` and the machine-readable
`comparison/results/reports/script_diagnostics_manifest.csv`.

| Gate | Models | Engines | Required outputs | Status |
|---|---|---|---|---|
| Distribution parity | AR(1), MA(1) | Quadra, RTMB | objective, gradient, transform | passing |
| Latent Laplace parity | centered AR(1) | Quadra, RTMB | mode, joint, logdet, marginal objective/gradient | passing |
| Isolated scaling | centered AR(1) | Quadra, RTMB | estimates, timings, RSS, Quadra diagnostics | passing through 1,000 states |
| Parameterization parity | AR good/bad/eps/manual | Quadra, RTMB | marginal objective/gradient, transformed latent mode, structure | passing |
| Parameterization parity | MA good/bad | Quadra, RTMB | marginal objective/gradient, transformed latent mode, structure | passing |
| Fisheries MAP/Laplace | simple FSA AR correlation + recruitment scale gate | Quadra, RTMB | all 440 observations; objective/gradient, optimized estimates/objective, Hessian structure | passing |
| Fisheries MAP/Laplace | AR observation/recruitment nuisance block | Quadra, RTMB | 10 fixed effects, 44 latent recruits, full gradients/estimates, diagnostics | passing |
| Fisheries MAP/Laplace | complete `basicfsa_ar_good.R` mapped fixed blocks | Quadra, RTMB | 66 fixed effects, 44 latent recruits, long-form parameter parity, full diagnostics | passing |
| Fisheries MAP/Laplace | complete MA good/bad mapped fixed blocks | Quadra, RTMB | 66 fixed effects, 45 latent effects, off-zero kernel probe, long-form parity, diagnostics | passing |
| Fisheries MAP/Laplace | NAA variants | Quadra, RTMB | source audit before numerical parity | source-blocked; audit report generated |
| Bayesian smoke test | six archived Bayes cases | Quadra kernel, tmbstan, Stan | log density differences, short-chain diagnostics | Quadra–RTMB kernels passing; Stan/tmbstan smoke execution pending |
| Bayesian full study | six cases, 100 replicates | tmbstan, Stan; Quadra diagnostics alongside | runtime, RSS, R-hat, ESS | deferred until smoke gates pass |

Quadra does not currently expose a native HMC sampler. The Bayesian comparison
therefore separates two questions:

1. Does Quadra reproduce the same posterior kernel, derivatives, transforms,
   MAP, and Laplace diagnostics?
2. How do tmbstan and Stan sampling efficiency change across the same good/bad
   parameterizations?

This avoids presenting Laplace optimization and posterior sampling as though
they were the same inference method.
