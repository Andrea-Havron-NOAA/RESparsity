# Quadra–RESparsity study summary

## Executive conclusion

Quadra reproduced every valid implemented RTMB target in this study. Objectives
and derivatives agree at numerical precision for the distribution, latent
Laplace, Bayesian-parameterization, and fisheries gates. Optimized fisheries
solutions also agree within the declared tolerances. No implemented gate found
evidence of a Quadra correctness defect.

For the matched centered AR(1) scaling case, Quadra's exact marginal gradient
is faster than RTMB at 30 and 100 latent states and within 9% at 300 states. At
1,000 states the measured times were 4.113 ms for Quadra and 1.700 ms for RTMB,
while peak RSS was 23.50 MiB and 307.88 MiB, respectively. The result is a large
memory advantage with a high-dimension gradient-time gap still to optimize.
These measurements describe this model, toolchain, and machine—not all possible
model structures.

## What was compared

| Study block | Coverage | Result |
|---|---|---|
| Observed AR(1) and MA(1) | objective, transform, AD gradient | passing |
| Latent centered AR(1) | mode, joint objective, log determinant, marginal objective and exact gradient | passing |
| AR Bayesian parameterizations | centered, manual centered, scaled innovations, standard innovations | passing |
| MA Bayesian parameterizations | innovation and inverse-filter forms | passing |
| Fisheries AR | scale gate, nuisance block, all 66 mapped fixed effects and 44 latent recruits | passing |
| Fisheries MA | complete good and bad forms, 66 fixed effects and 45 latent effects | passing |
| NAA variants | static source audit | source-blocked |
| Stan/tmbstan sampling | archived results only | not a native Quadra sampler comparison |

The detailed coverage contract is in [`COMPARISON_MATRIX.md`](COMPARISON_MATRIX.md).

## Numerical findings

- Observed AR(1) and MA(1) objective and gradient differences are at machine
  precision.
- The latent AR(1) conditional mode, joint objective, Hessian log determinant,
  Laplace objective, and exact marginal gradient agree with RTMB.
- All implemented AR and MA parameterizations reproduce the same marginal
  target after their intended transformations.
- Full fisheries AR and MA fits converge in both engines. Initial objective and
  gradient probes agree, and fitted parameter differences remain within each
  gate's declared tolerance.
- Quadra automatically selects tridiagonal factorization for the structured AR
  and innovation-form MA models and dense LDLT for the intentionally dense
  inverse-filter MA model.

## Final scaling result

| Latent states | Quadra exact gradient ms | RTMB gradient ms | Quadra/RTMB | Quadra RSS MiB | RTMB RSS MiB |
|---:|---:|---:|---:|---:|---:|
| 30 | 0.040 | 0.210 | 0.193 | 1.86 | 279.89 |
| 100 | 0.223 | 0.360 | 0.619 | 4.08 | 276.80 |
| 300 | 0.816 | 0.750 | 1.088 | 6.61 | 282.56 |
| 1,000 | 4.113 | 1.700 | 2.419 | 23.50 | 307.88 |

The benchmark uses exact marginal derivatives only. Finite differences are not
timed and do not participate in parity acceptance or optimization.

The main performance changes developed during the study were persistent tape
reuse, automatic Hessian-structure selection, structured log-determinants,
selected-inverse trace contraction, active-direction discovery, parallel Hdot
workers, flat Hessian/Hdot storage, one-time topology probing, and precomputed
reverse-propagation destinations. In the final run the 1,000-state directional
reverse phase was 0.436 ms; the remaining scaling gap lies outside that isolated
reverse phase and is visible in the full exact-gradient timing.

## Diagnostics

Every implemented case has a Markdown report containing convergence state,
curvature health, Hessian inertia, conditioning, sparsity and bandwidth,
backend selection, latent-state summaries, and the complete ordered Hessian
spectrum. Start with the generated [`results/reports/README.md`](results/reports/README.md)
for the script-to-report map.

Generated artifacts are separated into:

- [`results/reports/`](results/reports/): human-readable reports and diagnostic
  companion tables.
- [`results/raw/`](results/raw/): parity and long-form parameter tables.
- [`results/runtime-comparisons/`](results/runtime-comparisons/): timing, RSS,
  ratios, archived runtime summaries, and the scaling plot.

## Boundaries and unresolved work

The checked-in NAA scripts are internally incomplete: parameters are missing,
the declared shape of `logNAA` conflicts with its indexing, and the bad model's
random-effect declaration references an absent parameter. The
[`NAA source audit`](results/reports/naa_source_audit.md) records the required
source decisions. No replacement assumptions were invented.

Quadra does not expose a native HMC sampler. Consequently, this study validates
Quadra posterior kernels, derivatives, MAP/Laplace calculations, and
diagnostics, but it does not claim a Quadra-versus-Stan sampling comparison.
The archived Stan/tmbstan table is retained only as historical context.

Runtime results are isolated-process measurements from one machine and should
be rerun on target hardware before making portable performance claims. The
single-fit optimization timings are noisier than the repeated gradient timings.

## Reproduction

Run from the RESparsity repository root:

```sh
# All numerical and source-audit gates
Rscript comparison/run_all.R

# Isolated runtime and peak-RSS benchmark; also regenerates the PNG
Rscript comparison/ar1/run_scaling_benchmark.R
```

`run_all.R` regenerates the diagnostics index after every successful suite.
The scaling benchmark is separate because it launches isolated measured
processes and is slower than the correctness suite.

## Study disposition

The implemented study is complete. Quadra is numerically validated against
RTMB for every executable target currently in scope, has a substantial memory
advantage, and remains competitive through 300 states on the matched structured
case. Future work should begin as a
new study phase: repair and specify the NAA sources, add larger and block-banded
scaling models, and define a sampling comparison only if Quadra gains or is
paired with an explicitly selected posterior sampler.
