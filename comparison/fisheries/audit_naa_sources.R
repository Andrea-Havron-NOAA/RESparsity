#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
good_path <- file.path(root, "simplefsa", "basicNAA_ar_good.R")
bad_path <- file.path(root, "simplefsa", "basicNAA_ar_bad.R")
good <- readLines(good_path, warn = FALSE)
bad <- readLines(bad_path, warn = FALSE)

checks <- data.frame(
  script = c("simplefsa/basicNAA_ar_good.R",
             rep("simplefsa/basicNAA_ar_bad.R", 4)),
  issue = c(
    "logNAA is declared as a flat vector but indexed with [age, year]",
    "sdAR is referenced but never defined",
    "phiAR is referenced but never defined",
    "logNAA is referenced but absent from the parameter list",
    "random='logN1A' names a parameter absent from the parameter list"
  ),
  detected = c(
    any(grepl("logNAA\\s*=\\s*rep", good)) &&
      any(grepl("logNAA\\[a,y\\]", good)),
    any(grepl("sdAR", bad)) && !any(grepl("sdAR\\s*<-", bad)),
    any(grepl("phiAR", bad)) && !any(grepl("phiAR\\s*<-", bad)),
    any(grepl("logNAA", bad)) && !any(grepl("logNAA\\s*=", bad)),
    any(grepl('random="logN1A"', bad, fixed = TRUE)) &&
      !any(grepl("logN1A\\s*=", bad))
  ),
  stringsAsFactors = FALSE
)
stopifnot(all(checks$detected))

reports_dir <- file.path(root, "comparison", "results", "reports")
raw_dir <- file.path(root, "comparison", "results", "raw")
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(checks, file.path(raw_dir, "naa_source_audit.csv"), row.names = FALSE)

lines <- c(
  "# NAA Legacy Source Audit", "",
  "The NAA pair cannot yet be used as a faithful Quadra–RTMB parity target because the checked-in legacy scripts are internally incomplete. No replacement model assumptions were invented.",
  "", "## Findings", "",
  "| Legacy script | Blocking issue | Detected |", "|---|---|---|",
  sprintf("| `%s` | %s | `%s` |", checks$script, checks$issue,
          ifelse(checks$detected, "yes", "no")),
  "", "## Required source decisions", "",
  "1. Give `logNAA` an explicit `6 x 44` matrix shape, or replace two-dimensional indexing with the intended flat indexing.",
  "2. Define the noncentered AR quantities `sdAR` and `phiAR` in the bad model.",
  "3. Decide whether `x`, `logN1A`, or another block is the intended random effect in the bad model.",
  "4. Add the missing `logNAA` parameter/state construction to the bad model.",
  "", "Once those choices are resolved, this audit should be replaced by objective, gradient, mode, optimization, spectrum, and backend-selection parity gates."
)
writeLines(lines, file.path(reports_dir, "naa_source_audit.md"))
cat("PASS: NAA source audit reproduced all known blocking inconsistencies.\n")
