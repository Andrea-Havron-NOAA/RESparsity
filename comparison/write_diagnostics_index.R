#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
diagnostics_dir <- file.path(root, "comparison", "results", "reports")
script_dir <- file.path(diagnostics_dir, "by-script")
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(
  comparison_script = c(
    "comparison/ar1/run_parity.R",
    "comparison/ma1/run_parity.R",
    "comparison/ar1/run_laplace_parity.R",
    rep("comparison/bayes/run_ar_parameterization_parity.R", 4),
    rep("comparison/bayes/run_ma_parameterization_parity.R", 2),
    "comparison/fisheries/run_simplefsa_ar_gate.R"
  ),
  case = c(
    "observed AR(1)", "observed MA(1)", "latent centered AR(1)",
    "centered", "centered manual", "scaled innovations",
    "standard innovations", "innovation good", "inverse-filter bad",
    "simple FSA AR correlation and recruitment scale"
  ),
  diagnostic_report = c(
    "observed_ar1_fixed.md", "observed_ma1_fixed.md",
    "latent_ar1_laplace.md", "bayes_ar_centered.md",
    "bayes_ar_centered_manual.md", "bayes_ar_scaled_innovations.md",
    "bayes_ar_standard_innovations.md", "bayes_ma_innovation_good.md",
    "bayes_ma_inverse_filter_bad.md", "simplefsa_ar_scale.md"
  ),
  legacy_source = c(
    "R/ar1.R", "R/ma1.R", "R/ar1_simtest.R",
    "bayes/ar_good_comparison.R", "bayes/ar_good_nodautoreg_comparison.R",
    "bayes/ar_good_eps_comparison.R", "bayes/ar_bad_comparison.R",
    "bayes/ma_good_comparison.R", "bayes/ma_bad_comparison.R",
    "simplefsa/basicfsa_ar_good.R"
  ),
  status = "passing",
  stringsAsFactors = FALSE
)
manifest <- rbind(manifest, data.frame(
  comparison_script =
    "comparison/fisheries/run_simplefsa_ar_nuisance_gate.R",
  case = "simple FSA AR observation/recruitment nuisance block",
  diagnostic_report = "simplefsa_ar_nuisance.md",
  legacy_source = "simplefsa/basicfsa_ar_good.R",
  status = "passing",
  stringsAsFactors = FALSE
))
manifest <- rbind(manifest, data.frame(
  comparison_script = rep("comparison/fisheries/audit_naa_sources.R", 2),
  case = c("NAA AR good legacy source audit",
           "NAA AR bad legacy source audit"),
  diagnostic_report = rep("naa_source_audit.md", 2),
  legacy_source = c("simplefsa/basicNAA_ar_good.R",
                    "simplefsa/basicNAA_ar_bad.R"),
  status = "source_blocked",
  stringsAsFactors = FALSE
))
manifest <- rbind(manifest, data.frame(
  comparison_script = c(
    "comparison/fisheries/run_simplefsa_ma_good_full_gate.R",
    "comparison/fisheries/run_simplefsa_ma_bad_full_gate.R"),
  case = c("complete simple FSA MA innovation fixed blocks",
           "complete simple FSA MA inverse-filter fixed blocks"),
  diagnostic_report = c("simplefsa_ma_good_full.md",
                         "simplefsa_ma_bad_full.md"),
  legacy_source = c("simplefsa/basicfsa_ma_good_orig.R",
                    "simplefsa/basicfsa_ma_bad_orig.R"),
  status = "passing",
  stringsAsFactors = FALSE
))
manifest <- rbind(manifest, data.frame(
  comparison_script =
    "comparison/fisheries/run_simplefsa_ar_full_fixed_gate.R",
  case = "complete simple FSA AR mapped fixed blocks",
  diagnostic_report = "simplefsa_ar_full_fixed.md",
  legacy_source = "simplefsa/basicfsa_ar_good.R",
  status = "passing",
  stringsAsFactors = FALSE
))

manifest_path <- file.path(diagnostics_dir, "script_diagnostics_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)

slug <- function(path) {
  value <- sub("^comparison/", "", path)
  value <- sub("\\.R$", "", value)
  gsub("[^A-Za-z0-9]+", "_", value)
}

script_pages <- character()
for (script in unique(manifest$comparison_script)) {
  rows <- manifest[manifest$comparison_script == script, , drop = FALSE]
  page <- paste0(slug(script), ".md")
  script_pages <- c(script_pages, page)
  lines <- c(
    paste0("# Diagnostics for `", script, "`"),
    "",
    paste0("This page maps the executable comparison script to its Quadra ",
           "diagnostic reports and the legacy RESparsity scripts it represents."),
    "",
    "| Case | Quadra diagnostics | Legacy/source script | Status |",
    "|---|---|---|---|"
  )
  for (i in seq_len(nrow(rows))) {
    lines <- c(lines, sprintf(
      "| %s | [%s](../%s) | `%s` | `%s` |",
      rows$case[i], rows$diagnostic_report[i], rows$diagnostic_report[i],
      rows$legacy_source[i], rows$status[i]
    ))
  }
  lines <- c(
    lines, "", "## What the reports contain", "",
    "Each linked report records Quadra's objective/gradient state, curvature health, complete signed eigenvalue spectrum, Hessian inertia, effective structure, and backend-factorization selection. Latent models additionally include uncertainty and latent-state diagnostics.",
    ""
  )
  writeLines(lines, file.path(script_dir, page))
}

writeLines(c(
  "# Diagnostics for `comparison/ar1/run_scaling_benchmark.R`", "",
  "This benchmark is intentionally separate from `comparison/run_all.R` because it launches isolated timed processes at multiple state dimensions.", "",
  "## Artifacts", "",
  "- [Scaling results CSV](../../runtime-comparisons/ar1_scaling.csv)",
  "- [Scaling ratio CSV](../../runtime-comparisons/ar1_scaling_ratios.csv)",
  "- [Scaling plot](../../runtime-comparisons/ar1_scaling.png)", "",
  "## Legacy/source scripts", "",
  "- `R/simple_benchmark/ar1_benchmark.R`", "",
  "Quadra backend selection, structure, tape reuse, active directions, Hdot workers, timings, and peak RSS are recorded in the scaling CSV. Full eigenspectra are intentionally emitted by the smaller parity/diagnostic gates rather than during isolated timing runs."
), file.path(script_dir, "ar1_run_scaling_benchmark.md"))

writeLines(c(
  "# Diagnostics for `comparison/bayes/summarize_archived.R`", "",
  "This script summarizes historical Stan/tmbstan outputs. It does not execute a Quadra sampler and therefore does not produce a Quadra Hessian diagnostics report.", "",
  "## Artifact", "",
  "- [Archived Bayesian summary CSV](../../runtime-comparisons/bayes_archived_summary.csv)", "",
  "## Legacy/source scripts", "",
  "- `bayes/ar_bad_comparison.R`", "- `bayes/ar_good_comparison.R`",
  "- `bayes/ar_good_eps_comparison.R`",
  "- `bayes/ar_good_nodautoreg_comparison.R`",
  "- `bayes/ma_bad_comparison.R`", "- `bayes/ma_good_comparison.R`", "",
  "Use the AR/MA parameterization parity pages for the matching Quadra kernel and curvature diagnostics."
), file.path(script_dir, "bayes_summarize_archived.md"))

lines <- c(
  "# Quadra Comparison Diagnostics",
  "",
  "Use the script mapping below when reviewing or sharing results. Every executable parity script has a landing page that identifies its case-level diagnostics and the original RESparsity script being represented.",
  "",
  "## Comparison script mapping",
  "",
  "| Comparison script | Script diagnostics page | Case reports |",
  "|---|---|---:|"
)
scripts <- unique(manifest$comparison_script)
for (i in seq_along(scripts)) {
  count <- sum(manifest$comparison_script == scripts[i])
  lines <- c(lines, sprintf("| `%s` | [%s](by-script/%s) | %d |",
                            scripts[i], script_pages[i], script_pages[i], count))
}

lines <- c(
  lines, "", "## Coverage outside `run_all.R`", "",
  "| Script | Script page | Coverage status |",
  "|---|---|---|",
  "| `comparison/ar1/run_scaling_benchmark.R` | [ar1_run_scaling_benchmark.md](by-script/ar1_run_scaling_benchmark.md) | Passing through 1,000 states; separate benchmark because it is substantially slower |",
  "| `comparison/bayes/summarize_archived.R` | [bayes_summarize_archived.md](by-script/bayes_summarize_archived.md) | Historical Stan/tmbstan results only; not a native Quadra sampler comparison |",
  "", "## Remaining legacy scripts", "",
  "The following scripts are not claimed as covered by the current gates:",
  "",
  "- `simplefsa/basicfsa.R`: initial abundance and fishing-mortality blocks.",
  "- `simplefsa/basicfsa_ar_bad.R`: full fisheries noncentered AR model.",
  "- `simplefsa/basicNAA_ar_good.R` and `simplefsa/basicNAA_ar_bad.R`: source-blocked age-structured variants; see [`naa_source_audit.md`](naa_source_audit.md).",
  "- `R/combined_baby*.R`: combined historical study variants.",
  "",
  "These remain staged work; they are listed here to prevent a passing smaller gate from being mistaken for full-script coverage."
)
writeLines(lines, file.path(diagnostics_dir, "README.md"))
