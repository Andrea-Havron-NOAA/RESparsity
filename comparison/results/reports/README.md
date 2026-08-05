# Quadra Comparison Diagnostics

Use the script mapping below when reviewing or sharing results. Every executable parity script has a landing page that identifies its case-level diagnostics and the original RESparsity script being represented.

## Comparison script mapping

| Comparison script | Script diagnostics page | Case reports |
|---|---|---:|
| `comparison/ar1/run_parity.R` | [ar1_run_parity.md](by-script/ar1_run_parity.md) | 1 |
| `comparison/ma1/run_parity.R` | [ma1_run_parity.md](by-script/ma1_run_parity.md) | 1 |
| `comparison/ar1/run_laplace_parity.R` | [ar1_run_laplace_parity.md](by-script/ar1_run_laplace_parity.md) | 1 |
| `comparison/bayes/run_ar_parameterization_parity.R` | [bayes_run_ar_parameterization_parity.md](by-script/bayes_run_ar_parameterization_parity.md) | 4 |
| `comparison/bayes/run_ma_parameterization_parity.R` | [bayes_run_ma_parameterization_parity.md](by-script/bayes_run_ma_parameterization_parity.md) | 2 |
| `comparison/fisheries/run_simplefsa_ar_gate.R` | [fisheries_run_simplefsa_ar_gate.md](by-script/fisheries_run_simplefsa_ar_gate.md) | 1 |
| `comparison/fisheries/run_simplefsa_ar_nuisance_gate.R` | [fisheries_run_simplefsa_ar_nuisance_gate.md](by-script/fisheries_run_simplefsa_ar_nuisance_gate.md) | 1 |
| `comparison/fisheries/audit_naa_sources.R` | [fisheries_audit_naa_sources.md](by-script/fisheries_audit_naa_sources.md) | 2 |
| `comparison/fisheries/run_simplefsa_ma_good_full_gate.R` | [fisheries_run_simplefsa_ma_good_full_gate.md](by-script/fisheries_run_simplefsa_ma_good_full_gate.md) | 1 |
| `comparison/fisheries/run_simplefsa_ma_bad_full_gate.R` | [fisheries_run_simplefsa_ma_bad_full_gate.md](by-script/fisheries_run_simplefsa_ma_bad_full_gate.md) | 1 |
| `comparison/fisheries/run_simplefsa_ar_full_fixed_gate.R` | [fisheries_run_simplefsa_ar_full_fixed_gate.md](by-script/fisheries_run_simplefsa_ar_full_fixed_gate.md) | 1 |

## Coverage outside `run_all.R`

| Script | Script page | Coverage status |
|---|---|---|
| `comparison/ar1/run_scaling_benchmark.R` | [ar1_run_scaling_benchmark.md](by-script/ar1_run_scaling_benchmark.md) | Passing through 1,000 states; separate benchmark because it is substantially slower |
| `comparison/bayes/summarize_archived.R` | [bayes_summarize_archived.md](by-script/bayes_summarize_archived.md) | Historical Stan/tmbstan results only; not a native Quadra sampler comparison |

## Remaining legacy scripts

The following scripts are not claimed as covered by the current gates:

- `simplefsa/basicfsa.R`: initial abundance and fishing-mortality blocks.
- `simplefsa/basicfsa_ar_bad.R`: full fisheries noncentered AR model.
- `simplefsa/basicNAA_ar_good.R` and `simplefsa/basicNAA_ar_bad.R`: source-blocked age-structured variants; see [`naa_source_audit.md`](naa_source_audit.md).
- `R/combined_baby*.R`: combined historical study variants.

These remain staged work; they are listed here to prevent a passing smaller gate from being mistaken for full-script coverage.
