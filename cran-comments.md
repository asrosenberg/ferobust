<!-- Before submitting: run devtools::check_win_devel() and rhub::rhub_check(),
     then update the test-environment list and note results below. -->

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.3 — `R CMD check --as-cran`
* win-builder: R-devel (2026-06-08 r90120 ucrt) — 1 NOTE (new submission)
* R-hub: Linux / Windows / macOS (pending)

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
  "New submission." This is the first release of ferobust.
  The same NOTE flags possibly misspelled words in DESCRIPTION: Imbens, Manski,
  Mismeasured, mismeasured. These are false positives -- Imbens and Manski are
  author surnames (the Imbens-Manski confidence interval), and "mismeasured" is
  standard measurement-error terminology.
* checking for future file timestamps ... NOTE
  "unable to verify current time." This is a local issue reaching the time
  server, not a package problem.

`fixest` is used only conditionally — in examples, tests, and the vignette, each
guarded with `requireNamespace()` — so it is in Suggests rather than Imports.

## Downstream dependencies

None. This is a new package with no reverse dependencies.
