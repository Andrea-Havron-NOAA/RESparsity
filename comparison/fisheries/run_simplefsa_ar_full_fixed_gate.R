#!/usr/bin/env Rscript

status <- system2(file.path(R.home("bin"), "Rscript"),
                  file.path("comparison", "fisheries",
                            "run_simplefsa_ar_gate.R"),
                  env = "SIMPLEFSA_GATE=full")
if (status != 0L) stop("Complete simple FSA AR fixed-block gate failed")
