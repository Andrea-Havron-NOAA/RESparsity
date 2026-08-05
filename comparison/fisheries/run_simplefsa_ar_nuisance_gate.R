#!/usr/bin/env Rscript

status <- system2(file.path(R.home("bin"), "Rscript"),
                  file.path("comparison", "fisheries",
                            "run_simplefsa_ar_gate.R"),
                  env = "SIMPLEFSA_GATE=nuisance")
if (status != 0L) stop("Simple FSA AR nuisance-block gate failed")
