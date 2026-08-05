# Observed MA(1) Fixed-Effect Diagnostics

Quadra fixed-effect curvature diagnostics at the parity evaluation point. This model has no latent Hessian; the spectrum below is the complete fixed-effect objective Hessian spectrum.

## Evaluation

- Objective: `5.23915657861654`
- Gradient norm: `0.0864325587395479`
- Parameter count: `1`

## Curvature

- Positive definite: `yes`
- Minimum eigenvalue: `0.6280521194557`
- Maximum eigenvalue: `0.6280521194557`
- Condition number: `1`
- Hessian inertia (positive / near-zero / negative): `1 / 0 / 0`
- Normalized spectral entropy: `0`
- Participation ratio: `1`
- Stable rank: `1`

## Quadra Backend Selection

| Quantity | Selection |
|---|---|
| Detected Hessian structure | `diagonal` |
| Factorization backend | `diagonal` |
| Solver recommendation | `Newton` |
| Bandwidth | `0` |
| Expected complexity | `O(n)` |
| Symbolic reuse supported | `no` |
| Selection reason | zero off-diagonal bandwidth |

## Full Eigenvalue Spectrum

| Rank | Eigenvalue | Cumulative positive-curvature share |
|---:|---:|---:|
| 1 | `0.6280521194557` | `1.000000` |
