# Quadra–RTMB comparisons

The comparison suite first establishes numerical parity, then measures
performance. Run commands from the RESparsity repository root.

For final conclusions, scope boundaries, and the reproducibility handoff, see
[`STUDY_SUMMARY.md`](STUDY_SUMMARY.md).

```sh
Rscript comparison/run_all.R

# Or run the gates separately:
Rscript comparison/ar1/run_parity.R
Rscript comparison/ar1/run_laplace_parity.R
Rscript comparison/ma1/run_parity.R
Rscript comparison/bayes/run_ar_parameterization_parity.R
Rscript comparison/bayes/run_ma_parameterization_parity.R
Rscript comparison/ar1/run_scaling_benchmark.R
Rscript comparison/bayes/summarize_archived.R
```

See `comparison/COMPARISON_MATRIX.md` for the staged fisheries and Bayesian
comparison gates.

For script-to-report provenance, start at
`comparison/results/reports/README.md`. Each executable comparison script
has a landing page mapping its cases to the generated diagnostics and the
legacy RESparsity source script it represents. The machine-readable mapping is
`comparison/results/reports/script_diagnostics_manifest.csv`.

The AR(1) parity test compiles a temporary Quadra executable and compares its
objective, automatic-differentiation gradient, and transformed correlation
against RTMB for identical data and parameters. Set `QUADRA_ROOT` to use a
different Quadra checkout; otherwise the clean `quadra` worktree is used.

The latent AR(1) test additionally integrates out the state vector and checks
the conditional mode, joint objective, Hessian log-determinant, Laplace
objective, and marginal gradient.

The MA(1) parity test checks the innovation parameterization used by
RESparsity, including the innovation density, deterministic MA(1) transform,
observation likelihood, and derivative with respect to `theta`.

The AR(1) scaling benchmark runs every engine and state dimension in an
isolated process and writes
`comparison/results/runtime-comparisons/ar1_scaling.csv`. It
separates setup, first and warm Laplace evaluations, gradient evaluation, and
Hessian structure, and records peak RSS in MiB. Quadra uses the full exact
marginal gradient exposed by `quadra::stats::ExactLaplaceEvaluator`; RTMB uses
automatic differentiation.
The exact Quadra timing includes the derivative of the Hessian log-determinant
and is reported separately from its one-time structural setup cost.

For bandwidth-one Hessians, Quadra now contracts the log-determinant derivative
against an O(n) tridiagonal selected inverse instead of materializing inverse
columns. Derivative replay remains an important exact-gradient phase at higher
state dimensions. The final matched benchmark is faster than RTMB at 30 and
100 states, near parity at 300, and 2.42x slower at 1,000. Fixed-effect optimization uses
`quadra::stats::optimize_laplace`, an exact L-BFGS path that evaluates the
marginal objective and gradient together and reuses the evaluator's warm
latent mode and persistent tapes. The benchmark exercises that path at every
state dimension. No finite-difference path is timed or used for parity
acceptance.

The Quadra Laplace timings use the stateful public
`quadra::stats::LaplaceEvaluator`. It records the model tape once, warm-starts
the latent mode, detects diagonal/tridiagonal structure automatically, and
reuses the selected structured log-determinant backend. A direct-value guard
plus a nearby deterministic probe rebuilds the tape and reselects the backend
if parameter-dependent control flow changes the model topology.

## Current matched AR(1) result

The isolated-process run uses the same data, parameterization, latent states,
and marginal likelihood in Quadra and RTMB. Objective differences are below
`1e-8` and exact-gradient differences are below `1e-8` at every size.

| states | Quadra exact gradient ms | RTMB gradient ms | Quadra RSS MiB | RTMB RSS MiB |
|---:|---:|---:|---:|---:|
| 30 | 0.040 | 0.210 | 1.86 | 279.89 |
| 100 | 0.223 | 0.360 | 4.08 | 276.80 |
| 300 | 0.816 | 0.750 | 6.61 | 282.56 |
| 1,000 | 4.113 | 1.700 | 23.50 | 307.88 |

Quadra's diagnostics independently identify every random-effect Hessian as
tridiagonal, select the tridiagonal backend, retain one active fixed-effect
direction and one Hdot worker, and report zero objective/Hdot tape rebuilds.
At 1,000 states Quadra uses about 7.6% of RTMB's peak memory, while its exact
gradient takes 2.42x RTMB's time in this run. Both optimizers run and converge at
every size with estimate agreement below `1e-5`.

The exact Hdot timer is split into validation/replay, direction setup,
directional reverse propagation, and contraction. At 1,000 states the reverse
sweep accounts for nearly all of the separately timed Hdot phase; the full
gradient timing also includes objective, factorization, sensitivity, and tape
management work. Hoisting shared topology, registered-edge,
operation, and accumulator access out of the hot traversal reduced the
1,000-state reverse phase from roughly 0.805 ms to 0.712 ms. Precomputing the
immutable first- and second-parent destination slots then removed per-edge hash
lookups and reduced that phase further to roughly 0.436 ms in the final run.

After the discovery sweep, HAD also converts intermediate second-order edges
to stable flat slots. Later Hessian sweeps update contiguous scalar storage
instead of performing BTree queries and insertions. This reduced the 1,000-state
warm AR(1) Laplace evaluation from roughly 23 ms to below 1 ms on the benchmark
machine. Exact Hdot propagation now uses paired flat Hessian/Hdot slot arrays,
and the exact evaluator reuses the persistent sweep's joint fixed gradient and
mixed Hessian. Together these changes reduced the 1,000-state exact-gradient
evaluation from roughly 168 ms to about 42 ms. Optional explicit vectorization
of contiguous directional workspace clearing is guarded by
`QUADRA_ENABLE_SIMD`; the scalar implementation remains the default reference
path.

`ExactLaplaceEvaluator` also owns a persistent total-Hdot tape. Its graph,
intermediate slots, and propagation plan are built during structural setup and
replayed for later parameter values. The same direct-value and nearby-probe
guard rebuilds this tape when control flow changes. During development this
reduced the 1,000-state exact-gradient evaluation from about 42 ms to about
12 ms.

The exact evaluator now profiles its objective/mode, factorization,
mode-sensitivity solve, Hdot propagation, and trace-contraction phases. Reusing
the constructor's discovery evaluator as the persistent objective evaluator
removed a duplicate first-use tape recording. Single-fixed-effect models also
skip numeric active-direction discovery because there is no direction to
prune. On the same 1,000-state AR(1) case, these changes reduced structural
setup from about 122 ms to about 52 ms and the reported exact-gradient time to
below 3 ms in that development run. After rebasing onto Quadra's promoted
low-memory engine, the final matched 1,000-state measurement is 4.113 ms with
peak RSS reduced to 23.50 MiB. Multi-fixed-effect models retain automatic
direction discovery.

When discovery finds multiple active fixed-effect directions, the evaluator
now partitions them across persistent Hdot workers. Each worker owns its own AD
graph, so directional propagation is thread-safe; implicit mode sensitivities
are solved once before dispatch. `ExactLaplaceGradientEngineOptions::hdot_workers`
controls the count (`0` selects up to four workers automatically, `1` forces
the serial reference path). Single-direction models never create worker
threads. The automatic cap prevents every active direction from duplicating a
full persistent graph on machines with many cores.

Worker initialization now records one canonical Hdot tape. Other workers share
its immutable flat-edge registry and per-vertex propagation-slot plan while
keeping mutable sweep state independent. Active-direction discovery is executed on the canonical
persistent tape as well; the older temporary discovery providers have been
removed from this public path. In the Quadra-only four-fixed-effect state-space
study, this reduced 1,000-state exact-evaluator setup from roughly 355 ms to
67 ms and kept the three-worker setup premium near 3%. HAD replay operation
metadata—opcode, operands, and constants—is now separated from mutable vertex
state and shared too. The 1,000-state three-worker benchmark avoids about
1.15 MB of duplicated operation storage.

Generated per-case Markdown diagnostics are written to
`comparison/results/reports/`. Reports include curvature health,
conditioning, uncertainty structure, effective sparsity and bandwidth, latent
state summaries, and the complete ordered Hessian eigenvalue spectrum with
cumulative curvature shares. `README.md` in that directory indexes the
available reports.
