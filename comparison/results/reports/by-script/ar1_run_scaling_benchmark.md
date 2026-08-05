# Diagnostics for `comparison/ar1/run_scaling_benchmark.R`

This benchmark is intentionally separate from `comparison/run_all.R` because it launches isolated timed processes at multiple state dimensions.

## Artifacts

- [Scaling results CSV](../../runtime-comparisons/ar1_scaling.csv)
- [Scaling ratio CSV](../../runtime-comparisons/ar1_scaling_ratios.csv)
- [Scaling plot](../../runtime-comparisons/ar1_scaling.png)

## Legacy/source scripts

- `R/simple_benchmark/ar1_benchmark.R`

Quadra backend selection, structure, tape reuse, active directions, Hdot workers, timings, and peak RSS are recorded in the scaling CSV. Full eigenspectra are intentionally emitted by the smaller parity/diagnostic gates rather than during isolated timing runs.
