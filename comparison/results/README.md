# Comparison outputs

Generated comparison artifacts are grouped by purpose:

The study-level conclusions and limitations are in
[`../STUDY_SUMMARY.md`](../STUDY_SUMMARY.md).

- [`reports/`](reports/README.md): human-readable model diagnostics, full-spectrum analyses, factorization selections, and their compact diagnostic tables.
- [`raw/`](raw/): machine-readable parity, parameter, and source-audit tables.
- [`runtime-comparisons/`](runtime-comparisons/): timing, peak-RSS, scaling-ratio, and plot artifacts.

Run `Rscript comparison/run_all.R` from the repository root to regenerate the suite. The report index maps every comparison script to its corresponding Markdown report.
