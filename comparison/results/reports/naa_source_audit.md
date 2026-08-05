# NAA Legacy Source Audit

The NAA pair cannot yet be used as a faithful Quadra–RTMB parity target because the checked-in legacy scripts are internally incomplete. No replacement model assumptions were invented.

## Findings

| Legacy script | Blocking issue | Detected |
|---|---|---|
| `simplefsa/basicNAA_ar_good.R` | logNAA is declared as a flat vector but indexed with [age, year] | `yes` |
| `simplefsa/basicNAA_ar_bad.R` | sdAR is referenced but never defined | `yes` |
| `simplefsa/basicNAA_ar_bad.R` | phiAR is referenced but never defined | `yes` |
| `simplefsa/basicNAA_ar_bad.R` | logNAA is referenced but absent from the parameter list | `yes` |
| `simplefsa/basicNAA_ar_bad.R` | random='logN1A' names a parameter absent from the parameter list | `yes` |

## Required source decisions

1. Give `logNAA` an explicit `6 x 44` matrix shape, or replace two-dimensional indexing with the intended flat indexing.
2. Define the noncentered AR quantities `sdAR` and `phiAR` in the bad model.
3. Decide whether `x`, `logN1A`, or another block is the intended random effect in the bad model.
4. Add the missing `logNAA` parameter/state construction to the bad model.

Once those choices are resolved, this audit should be replaced by objective, gradient, mode, optimization, spectrum, and backend-selection parity gates.
