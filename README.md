# ferobust

Measurement-error diagnostics for fixed-effects panel regression.

For a continuous, slow-moving regressor measured with error, fixed effects can
*amplify* attenuation: they discard the between-unit variation that carries most
of the signal and identify the coefficient from low-reliability within-unit
movement. A coefficient that shrinks under fixed effects is therefore ambiguous.
It may be confounding removed, or it may be attenuation. `ferobust` implements the
diagnostic and reporting workflow from the working paper that tells these
apart, using quantities the analyst already has.

It computes the corrected within-reliability, partial-identification bounds with
Imbens-Manski confidence intervals, the breakdown reliability, and an
autocorrelation frontier with a certification rule, and routes each application to
a verdict: rescue, certified, complement, sign flip, or not identified.

## Installation

```r
# install.packages("remotes")
remotes::install_github("asrosenberg/ferobust")

# to include the vignette, build it on install (needs pandoc, which
# RStudio and Quarto provide):
remotes::install_github("asrosenberg/ferobust", build_vignettes = TRUE)
```

## Quick start

```r
library(ferobust)
data(panel_demo)

# reliability = "auto" reads the measurement-model posterior (here, x_sd)
# and computes the reliability itself, so you supply no reliability number
audit(y ~ x | unit + year, data = panel_demo,
      key_var = "x", reliability = "auto")
```

This reports the empirical ICC of the regressor, pooled OLS / fixed effects /
random effects, the corrected within-reliability, the identified set, the
Imbens-Manski interval, and a verdict. On `panel_demo` the bounds exclude zero, so
the case is a *rescue*: the fixed-effects estimate is attenuated, and the
corrected bounds recover an effect the bare coefficient understated.

When no reliability is available, the autocorrelation frontier bounds it from the
regressor's own persistence:

```r
frontier_audit(y ~ x | unit + year, data = panel_demo,
               key_var = "x", psi_max = c(0, 0.5, 0.7))
```

### V-Dem indices

For V-Dem variables, `reliability = "auto"` works directly: the package picks up
the variable's `_sd` column, or its `_codelow`/`_codehigh` credible interval, and
de-attenuates accordingly. No reliability number to look up.

## Learn more

```r
vignette("within-reliability-diagnostic", package = "ferobust")
```

## Reference

Rosenberg, A. "Reliable Panel Regression: A Default Workflow for
Slow-Moving, Mismeasured Variables." Unpublished manuscript.

## License

MIT (c) Andrew S. Rosenberg
