## Resubmission

This is a resubmission. The previous incoming pretest flagged the `.github`
directory (a GitHub Actions workflow used for R-hub checks) as "included in
error." It is now listed in `.Rbuildignore` and is no longer part of the tarball.

## Test environments

* local: macOS (aarch64-apple-darwin20), R 4.5.3 — `R CMD check --as-cran`
* win-builder: R-devel and Debian (CRAN incoming pretest)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission. The NOTE also flags possibly misspelled words in DESCRIPTION:
  Imbens, Manski, Mismeasured, mismeasured. These are false positives -- Imbens
  and Manski are author surnames (the Imbens-Manski confidence interval), and
  "mismeasured" is standard measurement-error terminology.

`fixest` is used only conditionally -- in examples, tests, and the vignette, each
guarded with `requireNamespace()` -- so it is in Suggests rather than Imports.

## Downstream dependencies

None. This is a new package with no reverse dependencies.
