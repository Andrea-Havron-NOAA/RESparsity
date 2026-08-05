#!/usr/bin/env Rscript

status <- system2(file.path(R.home("bin"), "Rscript"),
                  file.path("comparison", "fisheries", "run_simplefsa_ar_gate.R"),
                  env = "SIMPLEFSA_GATE=ma_bad")
if (status != 0L) stop("Complete simple FSA MA-bad gate failed")
