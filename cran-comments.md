<!-- Before submitting: run devtools::check_win_devel() and rhub::rhub_check(),
     then update the test-environment list and note results below. -->

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.3 — `R CMD check --as-cran`
* win-builder: R-devel (pending)
* R-hub: Linux / Windows / macOS (pending)

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
  "New submission." This is the first release of ferobust.
* checking for future file timestamps ... NOTE
  "unable to verify current time." This is a local issue reaching the time
  server, not a package problem.

`fixest` is used only conditionally — in examples, tests, and the vignette, each
guarded with `requireNamespace()` — so it is in Suggests rather than Imports.

## Downstream dependencies

None. This is a new package with no reverse dependencies.
