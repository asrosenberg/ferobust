# ferobust 0.1.0

First release. `ferobust` implements the within-reliability diagnostic for
fixed-effects panel regression with a slow-moving, mismeasured regressor: it asks
how much of the shrinkage from pooled OLS to fixed effects is measurement-error
attenuation rather than confounding removed.

* `audit()` runs the full workflow and returns the empirical ICC, the pooled /
  fixed effects / random effects estimates, the corrected within-reliability, the
  partial-identification bounds with an Imbens-Manski confidence interval, the
  breakdown reliability, and a routed verdict (rescue, complement, sign flip, or
  not identified).
* `reliability = "auto"` reads a measurement-model posterior (an `_sd` column, or
  a `_codelow` / `_codehigh` credible interval as V-Dem ships) and computes the
  reliability itself, so the analyst supplies no reliability number.
* `frontier_audit()` bounds the within-reliability from the regressor's own
  persistence when no measurement model is available.
* `breakdown_gamma()` reports the differential-measurement-error sensitivity of a
  rescue: the correlation between the regressor's measurement error and the
  outcome that would push the lower bound back to zero.
* Ships a synthetic panel (`panel_demo`) and the vignette
  `within-reliability-diagnostic`.
