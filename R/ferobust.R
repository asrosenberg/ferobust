#' Measurement Error Diagnostic for Fixed Effects Regression
#'
#' Computes the ICC, within-reliability, partial identification bounds,
#' Imbens-Manski confidence interval, breakdown reliability, and FE-FD
#' comparison for a fixed effects specification with a potentially
#' mismeasured independent variable.
#'
#' @param beta_pooled Pooled OLS coefficient on the key IV.
#' @param se_pooled Standard error of the pooled coefficient.
#' @param beta_fe Fixed effects coefficient on the key IV.
#' @param se_fe Standard error of the FE coefficient.
#' @param icc Intraclass correlation of the key IV. If NULL, computed from
#'   \code{x} and \code{unit}.
#' @param icc_time ICC of the time component (for two-way FE). Default 0.
#' @param lambda_range Numeric vector of length 2 giving the assumed range
#'   of overall reliability. Default \code{c(0.70, 0.95)}.
#' @param beta_fd First-differenced coefficient (optional, for FE-FD comparison).
#' @param x Numeric vector of the key IV (used to compute ICC if \code{icc} is NULL).
#' @param unit Unit identifier (used with \code{x} to compute ICC).
#' @param alpha Significance level for the confidence interval. Default 0.05.
#'
#' @return A list of class \code{ferobust} containing:
#'   \item{sign_test}{Logical: do pooled and FE have the same sign?}
#'   \item{shrinkage}{Fraction of pooled coefficient lost under FE.}
#'   \item{icc}{Intraclass correlation (one-way).}
#'   \item{icc_time}{Time ICC (if provided).}
#'   \item{lambda_w}{Within-reliability at each endpoint of lambda_range.}
#'   \item{bounds}{Identified set at each endpoint of lambda_range.}
#'   \item{im_ci}{Imbens-Manski confidence interval.}
#'   \item{breakdown}{Breakdown reliability lambda*.}
#'   \item{fe_fd}{List with FE, FD, and whether |FE| > |FD|.}
#'
#' @examples
#' # Resource curse (Haber & Menaldo 2011)
#' ferobust(beta_pooled = -0.055, se_pooled = 0.009,
#'          beta_fe = -0.018, se_fe = 0.007,
#'          icc = 0.83, lambda_range = c(0.70, 0.95))
#'
#' @export
ferobust <- function(beta_pooled, se_pooled, beta_fe, se_fe,
                     icc = NULL, icc_time = 0,
                     lambda_range = c(0.70, 0.95),
                     beta_fd = NULL,
                     x = NULL, unit = NULL,
                     alpha = 0.05) {

  # Compute ICC if not provided
  if (is.null(icc)) {
    if (is.null(x) || is.null(unit)) {
      stop("Either provide icc or both x and unit.")
    }
    icc <- compute_icc(x, unit)
  }

  # Sign test
  same_sign <- sign(beta_pooled) == sign(beta_fe)

  # Shrinkage
  shrinkage <- 1 - abs(beta_fe) / abs(beta_pooled)

  # Within-reliability at lambda_range endpoints
  icc_total <- icc + icc_time
  lw <- within_reliability(lambda_range, icc_total)
  names(lw) <- paste0("lambda_", lambda_range)

  # Bounds at each lambda
  bounds <- lapply(lambda_range, function(lam) {
    lw_val <- within_reliability(lam, icc_total)
    corrected <- beta_fe / lw_val
    lo <- min(corrected, beta_pooled)
    hi <- max(corrected, beta_pooled)
    c(lower = lo, upper = hi)
  })
  names(bounds) <- paste0("lambda_", lambda_range)

  # Imbens-Manski CI at lambda_min
  lw_min <- within_reliability(lambda_range[1], icc_total)
  beta_L <- beta_fe / lw_min
  beta_U <- beta_pooled
  se_L <- se_fe / lw_min
  se_U <- se_pooled
  lo <- min(beta_L, beta_U)
  hi <- max(beta_L, beta_U)
  se_lo <- ifelse(lo == beta_L, se_L, se_U)
  se_hi <- ifelse(hi == beta_U, se_U, se_L)
  cv <- im_critical_value(hi - lo, max(se_lo, se_hi), alpha)
  im_ci <- c(lower = unname(lo - cv * se_lo), upper = unname(hi + cv * se_hi))

  # Breakdown reliability
  bd <- breakdown_reliability(beta_pooled, beta_fe, icc_total)

  # FE-FD comparison
  fe_fd <- NULL
  if (!is.null(beta_fd)) {
    fe_fd <- list(
      beta_fe = beta_fe,
      beta_fd = beta_fd,
      fe_larger = abs(beta_fe) > abs(beta_fd)
    )
  }

  out <- list(
    sign_test = same_sign,
    shrinkage = shrinkage,
    icc = icc,
    icc_time = icc_time,
    lambda_range = lambda_range,
    lambda_w = lw,
    bounds = bounds,
    im_ci = im_ci,
    breakdown = bd,
    fe_fd = fe_fd,
    beta_pooled = beta_pooled,
    beta_fe = beta_fe,
    alpha = alpha
  )
  class(out) <- "ferobust"
  out
}


#' @export
print.ferobust <- function(x, ...) {
  cat("ferobust diagnostic\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("Pooled:     %8.4f    FE:     %8.4f\n", x$beta_pooled, x$beta_fe))
  cat(sprintf("Shrinkage:  %5.1f%%\n", x$shrinkage * 100))
  cat(sprintf("Same sign:  %s\n", ifelse(x$sign_test, "yes", "NO (sign flip)")))
  cat(sprintf("ICC (unit): %5.3f", x$icc))
  if (!is.null(x$icc_2way) && !is.na(x$icc_2way)) {
    cat(sprintf("    ICC (unit+year absorbed): %5.3f  [bounds use this two-way ICC]",
                x$icc_2way))
  } else if (!is.null(x$icc_time) && !is.na(x$icc_time) && x$icc_time > 0) {
    cat(sprintf("    ICC (time): %5.3f", x$icc_time))
  }
  cat("\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")

  if (!x$sign_test) {
    cat("Sign flip detected. Bounds cannot be applied.\n")
    cat("Interpretation: confounding, not ME.\n")
    return(invisible(x))
  }

  cat("Within-reliability:\n")
  for (i in seq_along(x$lambda_range)) {
    cat(sprintf("  at lambda = %.2f:  lambda_w = %.3f\n",
                x$lambda_range[i], x$lambda_w[i]))
  }
  cat(sprintf("Breakdown reliability (lambda*): %.3f\n", x$breakdown))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat("Identified set:\n")
  for (i in seq_along(x$lambda_range)) {
    b <- x$bounds[[i]]
    cat(sprintf("  at lambda = %.2f:  [%.4f, %.4f]\n",
                x$lambda_range[i], b["lower"], b["upper"]))
  }
  cat(sprintf("Imbens-Manski %d%% CI: [%.4f, %.4f]\n",
              round((1 - x$alpha) * 100), x$im_ci["lower"], x$im_ci["upper"]))
  cat(sprintf("Excludes zero:       %s\n",
              ifelse((x$im_ci["lower"] > 0) | (x$im_ci["upper"] < 0), "yes", "no")))
  if (!is.null(x$breakdown_gamma) && !all(is.na(x$breakdown_gamma))) {
    g <- range(x$breakdown_gamma, na.rm = TRUE)
    if (g[1] > 1) {
      cat("Differential-ME breakdown |gamma|: > 1 (robust to any admissible error-outcome correlation)\n")
    } else {
      cat(sprintf("Differential-ME breakdown |gamma|: [%.3f, %.3f]\n", g[1], g[2]))
      cat("  (error-outcome correlation that pushes the lower bound to zero)\n")
    }
  }

  if (!is.null(x$fe_fd)) {
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat(sprintf("FE-FD comparison: |FE| %s |FD|  (%s classical ME)\n",
                ifelse(x$fe_fd$fe_larger, ">", "<"),
                ifelse(x$fe_fd$fe_larger, "consistent with", "inconsistent with")))
  }

  # (A1) consistency checks via breakdown reliability and Oster delta
  if (!is.null(x$A1_check) || !is.null(x$oster)) {
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat("(A1) consistency checks:\n")
    if (!is.null(x$A1_check)) {
      cat(sprintf("  breakdown lambda* vs lambda_max:  %s\n",
                  toupper(x$A1_check$verdict)))
      cat(sprintf("    %s\n", x$A1_check$note))
    }
    if (!is.null(x$oster) && !is.na(x$oster$delta)) {
      cat(sprintf("  Oster delta (R^2_max = %.3f):     |delta| = %.2f  (%s)\n",
                  x$oster$R2_max, abs(x$oster$delta),
                  ifelse(x$oster$robust_threshold, "ROBUST: |delta| > 1",
                         "FRAGILE: |delta| <= 1")))
      cat(sprintf("    %s\n", x$oster$note))
    } else if (!is.null(x$oster)) {
      cat(sprintf("  Oster delta:                       %s\n", x$oster$note))
    }
  }

  if (!is.null(x$scope_warnings) && length(x$scope_warnings) > 0) {
    cat(paste(rep("-", 50), collapse = ""), "\n")
    cat("Scope warnings:\n")
    for (w in x$scope_warnings) cat("  - ", w, "\n", sep = "")
  }

  invisible(x)
}


#' Compute within-reliability (corrected formula)
#'
#' Under classical ME, the empirical ICC of the observed regressor satisfies
#' \code{icc_hat = ICC_true * lambda}, so substituting into the Griliches-Hausman
#' expression yields \code{lambda_w = (lambda - icc_hat) / (1 - icc_hat)}.
#' This is the form the working paper derives and recommends; it uses
#' only the assumed reliability and the empirical ICC of the observed data.
#'
#' @param lambda Assumed overall reliability (scalar or vector).
#' @param icc_hat Empirical ICC of the observed regressor.
#' @return Within-reliability, floored at 1e-3 when icc_hat exceeds lambda
#'   (data-inconsistency flag).
#' @export
within_reliability <- function(lambda, icc_hat) {
  pmax((lambda - icc_hat) / (1 - icc_hat), 1e-3)
}


#' Audit a panel regression for measurement-error amplification
#'
#' Formula-interface wrapper for \code{\link{ferobust}}: takes a formula and
#' data frame, fits pooled OLS, fixed effects, and random effects, and returns
#' a ferobust diagnostic with the corrected within-reliability, partial
#' identification bounds, Imbens-Manski CI, sign-test verdict, and two checks
#' of assumption (A1) of Proposition 2: a breakdown-reliability check (does
#' lambda* fall inside the defensible reliability range?) and Oster (2019)
#' delta with the Altonji-Elder-Taber sign assumption (|delta| > 1 is
#' "robust").
#'
#' @param formula A formula of the form \code{y ~ x + controls | unit + year},
#'   following the \code{fixest} convention: the part after \code{|} lists
#'   fixed-effects variables.
#' @param data A data frame.
#' @param reliability Either a numeric vector of length 2 giving the assumed
#'   reliability range, or a named numeric vector with elements \code{low},
#'   \code{mid}, \code{high}.
#' @param key_var Character name of the key regressor (the one the diagnostic
#'   reports for). Defaults to the first regressor.
#' @param re Logical; if \code{TRUE} (the default), also fit a random-effects
#'   model with year dummies as a single-number point summary. Requires
#'   \pkg{lme4}; skipped with a warning if unavailable.
#' @param posterior_sd Optional posterior-SD column name used when
#'   \code{reliability = "auto"}. If omitted, \code{audit()} looks for a column
#'   named \code{key_var_sd}, or for the 68\% credible interval
#'   \code{key_var_codelow}/\code{key_var_codehigh} (e.g.\ V-Dem), from which the
#'   posterior SD is taken as half the interval width.
#' @param replicate_var Optional repeated/proxy-measure column name used when
#'   \code{reliability = "auto"} and no posterior SD is available.
#' @param replicate_lambda Optional known reliability of \code{replicate_var}.
#'   If omitted, the repeated-measures route assumes equal reliability.
#' @param reliability_lookup Optional lookup table passed to
#'   \code{\link{estimate_lambda}} when \code{reliability = "auto"}.
#' @param reliability_default Default reliability vector used when
#'   \code{reliability = "auto"} finds no posterior SD, repeated measure, or
#'   lookup match.
#' @return An object of class \code{ferobust}.
#'
#' \code{ferobust_lm} is an alias for \code{audit}: same arguments, same
#' return value. Use whichever name reads better in the caller's pipeline.
#'
#' @examples
#' data(panel_demo)
#' if (requireNamespace("fixest", quietly = TRUE)) {
#'   # reliability = "auto" reads the x_sd posterior column automatically
#'   audit(y ~ x | unit + year, data = panel_demo,
#'         key_var = "x", reliability = "auto", re = FALSE)
#' }
#' @export
audit <- function(formula, data,
                  reliability = c(low = 0.85, mid = 0.90, high = 0.95),
                  key_var = NULL,
                  re = TRUE,
                  posterior_sd = NULL,
                  replicate_var = NULL,
                  replicate_lambda = NULL,
                  reliability_lookup = NULL,
                  reliability_default = c(low = 0.85, mid = 0.90, high = 0.95)) {
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required. Install with install.packages('fixest').",
         call. = FALSE)
  }
  if (isTRUE(re) && !requireNamespace("lme4", quietly = TRUE)) {
    warning("Package 'lme4' not available; skipping the random-effects summary.",
            call. = FALSE)
    re <- FALSE
  }

  # Parse formula: y ~ regressors | fixed_effects
  fparts <- strsplit(deparse1(formula), "\\|")[[1]]
  if (length(fparts) != 2) {
    stop("Formula must have form 'y ~ regressors | fixed_effects'.", call. = FALSE)
  }
  rhs_main <- trimws(fparts[1])
  fe_vars  <- trimws(strsplit(fparts[2], "\\+")[[1]])

  # Use R's formula machinery so transformed regressors (log(x), I(x^2)),
  # interactions, and multi-token terms are parsed correctly.
  main_form <- stats::as.formula(rhs_main)
  rhs_regs  <- attr(stats::terms(main_form), "term.labels")
  lhs       <- all.vars(main_form)[1]
  if (is.null(key_var)) key_var <- rhs_regs[1]

  # Unit FE is the first listed; time FE (if any) is the second.
  unit_var <- fe_vars[1]
  time_var <- if (length(fe_vars) >= 2) fe_vars[2] else NA_character_

  # Cluster variable for SEs
  data <- as.data.frame(data)
  extra_vars <- c(
    if (is.character(posterior_sd) && length(posterior_sd) == 1) posterior_sd else character(0),
    if (is.character(replicate_var) && length(replicate_var) == 1) replicate_var else character(0)
  )
  needed_vars <- unique(c(all.vars(main_form), fe_vars, extra_vars))
  needed_vars <- intersect(needed_vars, names(data))
  data <- data[stats::complete.cases(data[, needed_vars, drop = FALSE]), ]

  # Pooled OLS with year FE if a second FE variable is given
  if (length(fe_vars) >= 2) {
    pool_form <- stats::as.formula(paste0(rhs_main, " + factor(", fe_vars[2], ")"))
  } else {
    pool_form <- stats::as.formula(rhs_main)
  }
  m_pool <- fixest::feols(pool_form, data = data,
                          vcov = stats::as.formula(paste0("~", unit_var)))
  # FE
  m_fe <- fixest::feols(formula, data = data,
                        vcov = stats::as.formula(paste0("~", unit_var)))
  # Pull coefficients on the key variable
  b_P  <- stats::coef(m_pool)[key_var]
  se_P <- sqrt(stats::vcov(m_pool)[key_var, key_var])
  b_FE <- stats::coef(m_fe)[key_var]
  se_FE <- sqrt(stats::vcov(m_fe)[key_var, key_var])

  # RE summary (unit random + year fixed). Optional and failure-tolerant: a
  # convergence problem on a large panel degrades the RE point summary to NA
  # rather than aborting the whole diagnostic.
  b_RE <- NA_real_; se_RE <- NA_real_
  if (isTRUE(re)) {
    re_form <- if (length(fe_vars) >= 2) {
      stats::as.formula(paste0(rhs_main, " + factor(", fe_vars[2], ") + (1 | ",
                               unit_var, ")"))
    } else {
      stats::as.formula(paste0(rhs_main, " + (1 | ", unit_var, ")"))
    }
    m_re <- tryCatch(
      lme4::lmer(re_form, data = data, REML = FALSE,
                 control = lme4::lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 50000))),
      error = function(e) { warning("RE fit failed: ", conditionMessage(e),
                                    call. = FALSE); NULL })
    if (!is.null(m_re)) {
      b_RE  <- lme4::fixef(m_re)[key_var]
      se_RE <- sqrt(stats::vcov(m_re)[key_var, key_var])
    }
  }

  # ICC of the key regressor (evaluate the term so log(x) etc. resolve)
  key_vals <- tryCatch(eval(parse(text = key_var), envir = data),
                       error = function(e) data[[key_var]])
  icc_hat <- compute_icc(key_vals, data[[unit_var]])
  # Two-way absorbed ICC (the attenuation-relevant share for a two-way FE
  # estimate); reported alongside the one-way ICC when a time FE is present.
  icc_2way <- if (!is.na(time_var) && time_var %in% names(data)) {
    tryCatch(compute_icc_2way(key_vals, data[[unit_var]], data[[time_var]]),
             error = function(e) NA_real_)
  } else NA_real_

  # Resolve reliability when requested. "auto" is intentionally conservative:
  # it uses posterior SDs or an explicitly supplied replicate/lookup when
  # available, otherwise it falls back to the paper's convention.
  lambda_source <- NULL
  if (inherits(reliability, "ferobust_lambda")) {
    lambda_source <- reliability
    reliability <- reliability$reliability
  } else if (is.character(reliability) && length(reliability) == 1 &&
             identical(tolower(reliability), "auto")) {
    lambda_source <- estimate_lambda(
      data = data,
      variable = key_var,
      posterior_sd = posterior_sd,
      replicate = replicate_var,
      replicate_lambda = replicate_lambda,
      lookup = reliability_lookup,
      default = reliability_default,
      quiet = TRUE
    )
    reliability <- lambda_source$reliability
  }

  # Build result. For a two-way FE spec, the attenuation-relevant ICC is the
  # share absorbed jointly by unit and time, so the bounds use the two-way ICC
  # when it is available (passed as icc + icc_time); the one-way unit ICC is
  # still reported separately as icc_oneway.
  lambda_range <- range(reliability)
  icc_time_val <- if (!is.na(icc_2way)) max(0, icc_2way - icc_hat) else 0
  diag <- ferobust(beta_pooled = b_P, se_pooled = se_P,
                   beta_fe = b_FE, se_fe = se_FE,
                   icc = icc_hat, icc_time = icc_time_val,
                   lambda_range = lambda_range)
  diag$icc_oneway <- icc_hat

  diag$beta_re <- b_RE
  diag$se_re   <- se_RE
  diag$key_var <- key_var
  diag$reliability_vec <- reliability
  diag$lambda_source <- lambda_source
  diag$icc_2way <- icc_2way
  diag$icc_time <- if (is.na(icc_2way)) NA_real_ else max(0, icc_2way - icc_hat)

  # Panel-shape diagnostics for small-T / small-N scope warnings
  unit_obs <- table(data[[unit_var]])
  T_med    <- as.numeric(stats::median(unit_obs))
  n_units  <- length(unit_obs)
  diag$n_units <- n_units
  diag$T_median <- T_med
  diag$scope_warnings <- character(0)
  if (T_med < 20) {
    diag$scope_warnings <- c(diag$scope_warnings,
      sprintf("small-T (median T = %g < 20): bounds may under-cover by 5-7pp; report Imbens-Manski CI rather than bare bounds.",
              T_med))
  }
  if (n_units < 30 && T_med < 15) {
    diag$scope_warnings <- c(diag$scope_warnings,
      sprintf("small-N short-T (N = %d, median T = %g): bounds coverage falls to ~85%%; rely on the Imbens-Manski CI.",
              n_units, T_med))
  }

  # ---- (A1) consistency check via breakdown reliability vs. lambda_max ----
  # (A1) says pooled OLS is biased away from zero, so beta lies between FE/lambda_w
  # and pooled OLS. If the breakdown reliability lambda* exceeds the upper end of
  # the defensible reliability range, the bounds invert at every defensible lambda,
  # which is a sign that (A1) and (A2) are jointly inconsistent in this data.
  lambda_max <- max(lambda_range)
  lambda_min <- min(lambda_range)
  if (!is.na(diag$breakdown)) {
    if (diag$breakdown <= lambda_max) {
      diag$A1_check <- list(
        verdict = "consistent",
        breakdown = diag$breakdown,
        lambda_max = lambda_max,
        note = sprintf("breakdown lambda* = %.3f <= lambda_max = %.2f; bounds non-degenerate at upper end.",
                       diag$breakdown, lambda_max)
      )
    } else {
      diag$A1_check <- list(
        verdict = "inconsistent",
        breakdown = diag$breakdown,
        lambda_max = lambda_max,
        note = sprintf("breakdown lambda* = %.3f > lambda_max = %.2f; (A1) and (A2) are jointly inconsistent at every defensible reliability.",
                       diag$breakdown, lambda_max)
      )
    }
  } else if (!diag$sign_test) {
    diag$A1_check <- list(verdict = "not applicable",
                          note = "sign flip; (A1)-based bounds do not apply.")
  } else {
    diag$A1_check <- list(verdict = "violated",
                          note = "same-sign but |FE| > |pooled OLS|; (A1) is rejected before ME correction. Suggests non-classical ME or that the cross-sectional pattern understates the within-unit effect.")
  }

  # ---- Oster delta as a partial (A1) check via Altonji-Elder-Taber ----
  # Oster (2019) delta is the proportional ratio of selection on unobservables
  # vs. observables that would wipe out the long-model effect. Under the
  # Altonji-Elder-Taber sign assumption (unobservables push in the same direction
  # as observables), |delta| > 1 is the conventional robustness threshold.
  # Short regression: y ~ x (treatment alone, no FE, no controls).
  # Long regression:  y ~ x + controls + factor(unit) + factor(year) (FE expanded).
  oster_delta_check <- tryCatch({
    other_regs <- setdiff(rhs_regs, key_var)
    short_form <- stats::as.formula(paste0("`", lhs, "` ~ `", key_var, "`"))
    m_short <- stats::lm(short_form, data = data)
    b_short <- stats::coef(m_short)[paste0("`", key_var, "`")]
    if (is.na(b_short)) b_short <- stats::coef(m_short)[key_var]
    R2_short <- summary(m_short)$r.squared

    long_rhs <- c(paste0("`", key_var, "`"),
                  if (length(other_regs) > 0) paste0("`", other_regs, "`") else character(0),
                  paste0("factor(`", fe_vars, "`)"))
    long_form <- stats::as.formula(paste0("`", lhs, "` ~ ",
                                           paste(long_rhs, collapse = " + ")))
    m_long  <- stats::lm(long_form, data = data)
    b_long  <- stats::coef(m_long)[paste0("`", key_var, "`")]
    if (is.na(b_long)) b_long <- stats::coef(m_long)[key_var]
    R2_long <- summary(m_long)$r.squared

    R2_max <- min(1.3 * R2_long, 1)
    num <- as.numeric(b_long * (R2_max - R2_long))
    den <- as.numeric((b_short - b_long) * (R2_long - R2_short))
    delta <- if (abs(den) < 1e-12) NA_real_ else num / den

    list(delta = as.numeric(delta),
         b_short = as.numeric(b_short), b_long = as.numeric(b_long),
         R2_short = R2_short, R2_long = R2_long, R2_max = R2_max,
         robust_threshold = !is.na(delta) && abs(delta) > 1,
         note = if (is.na(delta)) {
           "delta undefined (denominator near zero)."
         } else if (abs(delta) > 1) {
           sprintf("|delta| = %.2f > 1: under Altonji-Elder-Taber, unobservables would need to be MORE important than observables to overturn the long-model effect; partial evidence consistent with (A1).",
                   abs(delta))
         } else {
           sprintf("|delta| = %.2f <= 1: unobservables LESS important than observables can overturn the effect; (A1) is fragile by Oster's criterion.",
                   abs(delta))
         })
  }, error = function(e) {
    list(delta = NA_real_, note = sprintf("Oster delta computation failed: %s", conditionMessage(e)))
  })
  diag$oster <- oster_delta_check

  # Diagnostic verdict, in the paper's taxonomy:
  # sign flip / complement / rescue / not identified. ("certified" is a frontier
  # verdict; see frontier_audit().)
  im_excl0 <- !is.na(diag$im_ci["lower"]) &&
    ((diag$im_ci["lower"] > 0) || (diag$im_ci["upper"] < 0))
  diag$verdict <- if (!diag$sign_test) {
    "sign flip"
  } else if (abs(b_FE) > abs(b_P)) {
    "complement"
  } else if (im_excl0) {
    "rescue"
  } else {
    "not identified"
  }

  # ---- Differential-ME breakdown gamma (confounding margin) -----------------
  # For a rescue, zero-exclusion rests on the lower (pooled) bound. A single
  # differential-ME correlation gamma between the regressor's measurement error
  # and the outcome biases the pooled estimate by gamma * sqrt(1 - lambda) *
  # sd_e / sd_x, with sd_x the SD of the key regressor at the pooled (time-
  # absorbed) level and sd_e the pooled residual SD. The lower bound reaches zero
  # at |gamma| = |b_P| * sd_x / (sqrt(1 - lambda) * sd_e), reported across the
  # reliability range. This is the confounding-margin sensitivity, defined only
  # for the rescue verdict; values above 1 mean no admissible correlation can
  # overturn it.
  diag$breakdown_gamma <- NA_real_
  if (identical(diag$verdict, "rescue")) {
    px_form <- if (!is.na(time_var) && time_var %in% names(data)) {
      stats::as.formula(paste0(key_var, " ~ 1 | ", time_var))
    } else {
      stats::as.formula(paste0(key_var, " ~ 1"))
    }
    sd_x_P <- tryCatch(stats::sd(stats::resid(fixest::feols(px_form, data = data))),
                       error = function(e) NA_real_)
    sd_e_P <- stats::sd(stats::resid(m_pool))
    if (!is.na(sd_x_P) && sd_x_P > 0 && !is.na(sd_e_P) && sd_e_P > 0) {
      diag$breakdown_gamma <- abs(b_P) * sd_x_P / (sqrt(1 - reliability) * sd_e_P)
    }
  }

  if (length(diag$scope_warnings) > 0) {
    for (w in diag$scope_warnings) warning(w, call. = FALSE)
  }

  diag
}


#' Differential-measurement-error breakdown gamma for a rescue
#'
#' For a \dQuote{rescue} verdict, zero-exclusion rests on the lower (pooled)
#' bound of the identified set. This returns the breakdown \eqn{|\gamma|}, the
#' correlation between the regressor's measurement error and the outcome that
#' would push that bound to zero, evaluated across the reliability range. It is
#' the confounding-margin sensitivity reported in the paper: a thin wrapper
#' around \code{\link{audit}} (random effects skipped for speed) that returns
#' \code{NA} for any verdict other than a rescue.
#'
#' @inheritParams audit
#' @return A named numeric vector of \eqn{|\gamma|} thresholds, one per point in
#'   \code{reliability}, or \code{NA} when the verdict is not a rescue.
#' @examples
#' if (requireNamespace("fixest", quietly = TRUE)) {
#'   breakdown_gamma(y ~ x | unit + year, data = panel_demo,
#'                   key_var = "x", reliability = "auto")
#' }
#' @export
breakdown_gamma <- function(formula, data,
                            reliability = c(low = 0.85, mid = 0.90, high = 0.95),
                            key_var = NULL) {
  audit(formula, data, reliability = reliability, key_var = key_var,
        re = FALSE)$breakdown_gamma
}

#' @rdname audit
#' @export
ferobust_lm <- function(formula, data,
                        reliability = c(low = 0.85, mid = 0.90, high = 0.95),
                        key_var = NULL,
                        re = TRUE,
                        posterior_sd = NULL,
                        replicate_var = NULL,
                        replicate_lambda = NULL,
                        reliability_lookup = NULL,
                        reliability_default = c(low = 0.85, mid = 0.90, high = 0.95)) {
  audit(formula = formula, data = data,
        reliability = reliability, key_var = key_var, re = re,
        posterior_sd = posterior_sd,
        replicate_var = replicate_var,
        replicate_lambda = replicate_lambda,
        reliability_lookup = reliability_lookup,
        reliability_default = reliability_default)
}


#' Plug-and-play frontier audit for a panel regression
#'
#' Runs \code{\link{audit}} and adds the autocorrelation frontier for
#' within-reliability. This is the single-series diagnostic for point-estimate
#' variables when no posterior SD, repeated measure, or literature reliability is
#' available. The output reports a lower bound on \code{lambda_w} across
#' assumptions about the maximum persistence of measurement error.
#'
#' @param formula A formula of the form \code{y ~ x + controls | unit + year}.
#'   The first fixed effect is treated as the unit and the second as time.
#' @param data A data frame.
#' @param key_var Character name of the key regressor. Defaults to the first
#'   regressor.
#' @param psi_max Numeric vector of maximum allowed one-period
#'   measurement-error persistence values.
#' @param max_lag Maximum lag used in the autocorrelation frontier.
#' @param reliability Optional conventional reliability range still passed to
#'   \code{\link{audit}} for the standard diagnostic table.
#' @return An object of class \code{ferobust_frontier}.
#' @export
frontier_audit <- function(formula, data,
                           key_var = NULL,
                           psi_max = c(0, 0.3, 0.5, 0.7, 0.9, 1),
                           max_lag = 5,
                           reliability = c(low = 0.85, mid = 0.90, high = 0.95)) {
  fparts <- strsplit(deparse1(formula), "\\|")[[1]]
  if (length(fparts) != 2) {
    stop("Formula must have form 'y ~ regressors | unit + time'.",
         call. = FALSE)
  }
  rhs_main <- trimws(fparts[1])
  fe_vars <- trimws(strsplit(fparts[2], "\\+")[[1]])
  if (length(fe_vars) < 2) {
    stop("frontier_audit() needs both unit and time fixed effects, e.g. y ~ x | country + year.",
         call. = FALSE)
  }

  lhs <- trimws(strsplit(rhs_main, "~")[[1]][1])
  rhs_regs <- trimws(strsplit(trimws(strsplit(rhs_main, "~")[[1]][2]), "\\+")[[1]])
  if (is.null(key_var)) key_var <- rhs_regs[1]

  unit_var <- fe_vars[1]
  time_var <- fe_vars[2]
  data <- as.data.frame(data)
  needed <- unique(c(lhs, rhs_regs, fe_vars))
  d <- data[stats::complete.cases(data[, needed, drop = FALSE]), ]

  diag <- audit(formula = formula, data = d,
                reliability = reliability, key_var = key_var)
  frontier <- within_reliability_frontier(
    x = d[[key_var]],
    unit = d[[unit_var]],
    time = d[[time_var]],
    psi_max = psi_max,
    max_lag = max_lag
  )

  frontier$corrected_fe_frontier <- ifelse(
    frontier$lambda_w_lower > 0,
    diag$beta_fe / frontier$lambda_w_lower,
    ifelse(diag$beta_fe >= 0, Inf, -Inf)
  )

  bounds <- t(vapply(seq_len(nrow(frontier)), function(i) {
    endpoints <- c(diag$beta_fe, diag$beta_pooled,
                   frontier$corrected_fe_frontier[i])
    c(lower = min(endpoints, na.rm = TRUE),
      upper = max(endpoints, na.rm = TRUE))
  }, numeric(2)))
  frontier$frontier_lower <- bounds[, "lower"]
  frontier$frontier_upper <- bounds[, "upper"]

  # Certification corollary (Cor. 2): with matching signs and |FE| < |pooled|,
  # the frontier rules out the measurement-error reading iff the frontier floor
  # exceeds the shrinkage ratio r = |FE / pooled|.
  r_shrink <- abs(diag$beta_fe) / abs(diag$beta_pooled)
  same_sign_shrink <- diag$sign_test && abs(diag$beta_fe) < abs(diag$beta_pooled)
  frontier$shrinkage_ratio <- r_shrink
  frontier$verdict <- ifelse(
    !same_sign_shrink, "n/a (not a same-sign shrinkage case)",
    ifelse(frontier$lambda_w_lower > r_shrink,
           "certify (ME rescue ruled out)", "rescue possible"))

  out <- list(
    diagnostic = diag,
    frontier = frontier,
    shrinkage_ratio = r_shrink,
    same_sign_shrink = same_sign_shrink,
    key_var = key_var,
    unit_var = unit_var,
    time_var = time_var,
    max_lag = max_lag
  )
  class(out) <- "ferobust_frontier"
  out
}


#' Print a frontier audit
#' @param x A \code{ferobust_frontier} object.
#' @param ... Ignored.
#' @export
print.ferobust_frontier <- function(x, ...) {
  diag <- x$diagnostic
  cat("ferobust frontier audit\n")
  cat(paste(rep("-", 58), collapse = ""), "\n")
  cat(sprintf("Key variable: %s\n", x$key_var))
  cat(sprintf("Panel:        %s x %s\n", x$unit_var, x$time_var))
  cat(sprintf("Pooled:       % .4f\n", diag$beta_pooled))
  cat(sprintf("FE:           % .4f\n", diag$beta_fe))
  cat(sprintf("Same sign:    %s\n", ifelse(diag$sign_test, "yes", "NO")))
  cat(sprintf("ICC:          %.3f\n", diag$icc))
  cat(sprintf("Shrinkage r = |FE/pooled|: %.3f   (Cor. 2: certify when floor > r)\n",
              x$shrinkage_ratio))
  cat(paste(rep("-", 58), collapse = ""), "\n")
  cat("Within-reliability frontier (lambda_w lower bound):\n")
  for (i in seq_len(nrow(x$frontier))) {
    row <- x$frontier[i, ]
    corr <- if (is.finite(row$corrected_fe_frontier)) {
      sprintf("% .4f", row$corrected_fe_frontier)
    } else {
      ifelse(row$corrected_fe_frontier > 0, "+Inf", "-Inf")
    }
    lower <- if (is.finite(row$frontier_lower)) {
      sprintf("% .4f", row$frontier_lower)
    } else {
      ifelse(row$frontier_lower > 0, "+Inf", "-Inf")
    }
    upper <- if (is.finite(row$frontier_upper)) {
      sprintf("% .4f", row$frontier_upper)
    } else {
      ifelse(row$frontier_upper > 0, "+Inf", "-Inf")
    }
    cat(sprintf("  psi_max = %.1f: lambda_w >= %.3f, set = [%s, %s]  ->  %s\n",
                row$psi_max, row$lambda_w_lower, lower, upper, row$verdict))
  }
  if (!diag$sign_test) {
    cat(paste(rep("-", 58), collapse = ""), "\n")
    cat("Sign flip detected: frontier bounds are diagnostic only; report FE.\n")
  }
  invisible(x)
}


#' Compute the empirical ICC from data
#'
#' Returns the empirical intraclass correlation \code{ICC_hat = Var(xbar_i) /
#' Var(x)}, the share of the regressor's variance that lies between units. This
#' is the obs-weighted between/total decomposition (each observation contributes
#' its unit mean), which equals the standard one-way ANOVA ICC up to the
#' finite-sample noise correction and is the quantity the within-reliability
#' formula requires. (An earlier version used the unweighted variance of the
#' unit means, which differs on unbalanced panels.)
#'
#' @param x Numeric vector.
#' @param unit Unit identifier.
#' @return Scalar ICC in \code{[0, 1]}.
#' @export
compute_icc <- function(x, unit) {
  ok <- !is.na(x) & !is.na(unit)
  x <- suppressWarnings(as.numeric(x[ok]))
  unit <- unit[ok]
  if (length(x) < 2 || stats::var(x) <= 0) return(NA_real_)
  means <- base::tapply(x, unit, mean, na.rm = TRUE)
  # Between/total: variance of the unit-mean series across observations
  # (each obs weighted by its unit), i.e. Var(xbar_i) in the decomposition.
  stats::var(means[as.character(unit)], na.rm = TRUE) / stats::var(x, na.rm = TRUE)
}

#' Compute the jointly-absorbed (two-way) ICC for unit + time fixed effects
#'
#' Returns the share of the regressor's variance absorbed jointly by unit and
#' time fixed effects, \code{1 - Var(within)/Var(x)}, where the within residual
#' is the two-way demeaned regressor. This is the attenuation-relevant ICC for a
#' two-way fixed-effects estimate (Equation for lambda_w (two-way) in the paper).
#'
#' @param x Numeric vector.
#' @param unit Unit identifier.
#' @param time Time identifier.
#' @return Scalar two-way absorbed ICC in \code{[0, 1]}.
#' @export
compute_icc_2way <- function(x, unit, time) {
  ok <- !is.na(x) & !is.na(unit) & !is.na(time)
  x <- suppressWarnings(as.numeric(x[ok])); unit <- unit[ok]; time <- time[ok]
  if (length(x) < 3 || stats::var(x) <= 0) return(NA_real_)
  df <- data.frame(x = x, u = factor(unit), t = factor(time))
  if (requireNamespace("fixest", quietly = TRUE)) {
    r <- stats::resid(fixest::feols(x ~ 1 | u + t, data = df))
  } else {
    r <- stats::resid(stats::lm(x ~ u + t, data = df))
  }
  1 - stats::var(r) / stats::var(x)
}


#' Bound within-reliability from autocorrelation
#'
#' Computes an assumption-indexed lower bound on the FE attenuation factor
#' \code{lambda_w = Var(x*_within) / Var(x_within)}. If measurement-error
#' autocorrelation is bounded above by \code{psi_max^k} at lag \code{k}, then
#' the observed within autocorrelation \code{rho_k} implies
#' \code{lambda_w >= (rho_k - psi_max^k) / (1 - psi_max^k)}. The function takes
#' the maximum bound over lags.
#'
#' This is a frontier, not a point estimate. At \code{psi_max = 1}, persistent
#' signal and persistent measurement error are observationally equivalent and
#' the bound collapses to zero.
#'
#' @param x Numeric vector of the key regressor.
#' @param unit Unit identifier.
#' @param time Time identifier used to order observations within unit.
#' @param psi_max Numeric vector of maximum allowed one-period measurement-error
#'   persistence values.
#' @param max_lag Maximum lag to use in the autocorrelation frontier.
#' @return A data frame with \code{psi_max}, the lower bound
#'   \code{lambda_w_lower}, and the lag that binds.
#' @export
within_reliability_frontier <- function(x, unit, time,
                                        psi_max = c(0, 0.3, 0.5, 0.7, 0.9),
                                        max_lag = 5) {
  ok <- !is.na(x) & !is.na(unit) & !is.na(time)
  x <- suppressWarnings(as.numeric(x[ok]))
  unit <- unit[ok]
  time <- time[ok]
  if (length(x) < 3 || stats::var(x) <= 0) {
    stop("Need at least three non-missing observations with positive variance.",
         call. = FALSE)
  }
  if (any(!is.finite(psi_max)) || any(psi_max < 0) || any(psi_max > 1)) {
    stop("psi_max must contain values in [0, 1].", call. = FALSE)
  }
  if (!is.finite(max_lag) || max_lag < 1) {
    stop("max_lag must be a positive integer.", call. = FALSE)
  }
  max_lag <- as.integer(max_lag)

  x_dm <- x - stats::ave(x, unit, FUN = function(z) mean(z, na.rm = TRUE))
  v0 <- sum(x_dm^2, na.rm = TRUE)
  if (!is.finite(v0) || v0 <= 0) {
    stop("Within-unit variance is zero; frontier is undefined.", call. = FALSE)
  }

  # Gap-aware within-unit autocorrelation: rho(k) = sum_t z_t z_{t-k} / sum z^2,
  # summing only over pairs whose time index differs by exactly k, so calendar
  # gaps are respected (matches the paper's within_acf).
  rhos <- within_acf(x_dm, unit, time, max_lag, v0)

  rows <- lapply(psi_max, function(psi) {
    if (psi >= 1) {
      return(data.frame(psi_max = psi, lambda_w_lower = 0,
                        binding_lag = NA_integer_))
    }
    bounds <- vapply(seq_len(max_lag), function(k) {
      denom <- 1 - psi^k
      if (!is.finite(rhos[k]) || denom <= 1e-8) return(NA_real_)
      max(0, (rhos[k] - psi^k) / denom)
    }, numeric(1))
    if (all(is.na(bounds))) {
      data.frame(psi_max = psi, lambda_w_lower = NA_real_,
                 binding_lag = NA_integer_)
    } else {
      j <- which.max(bounds)
      data.frame(psi_max = psi,
                 lambda_w_lower = min(bounds[j], 0.999),
                 binding_lag = j)
    }
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}


within_acf <- function(x_dm, unit, time, max_lag, v0) {
  byu <- split(data.frame(z = x_dm, t = as.numeric(time)), unit)
  vapply(seq_len(max_lag), function(k) {
    num <- sum(vapply(byu, function(d) {
      if (nrow(d) <= k) return(0)
      m <- match(d$t - k, d$t)          # for each obs, the obs exactly k periods earlier
      ok <- !is.na(m)
      if (!any(ok)) return(0)
      sum(d$z[ok] * d$z[m[ok]], na.rm = TRUE)
    }, numeric(1)))
    num / v0
  }, numeric(1))
}


#' Compute breakdown reliability
#'
#' The breakdown reliability lambda* is the assumed reliability at which the
#' corrected debiased FE estimate equals pooled OLS — i.e., the reliability
#' threshold above which the bounds collapse onto a point. Under the corrected
#' formula \code{lambda_w = (lambda - icc_hat)/(1 - icc_hat)}, this is
#' \code{lambda* = (icc_hat + |beta_fe| * (1 - icc_hat) / |beta_pooled|)}.
#'
#' @param beta_pooled Pooled OLS coefficient.
#' @param beta_fe FE coefficient.
#' @param icc_hat Empirical ICC.
#' @return Breakdown reliability lambda*.
#' @export
breakdown_reliability <- function(beta_pooled, beta_fe, icc_hat) {
  if (abs(beta_pooled) <= abs(beta_fe) || sign(beta_pooled) != sign(beta_fe)) {
    return(NA_real_)
  }
  ratio <- abs(beta_fe) / abs(beta_pooled)
  icc_hat + ratio * (1 - icc_hat)
}


#' Imbens-Manski critical value
#' @param delta Width of the identified set.
#' @param se_max Maximum SE of the bounds.
#' @param alpha Significance level.
#' @return Critical value c_alpha.
im_critical_value <- function(delta, se_max, alpha = 0.05) {
  C <- delta / se_max
  f <- function(cv) stats::pnorm(C + cv) - stats::pnorm(-cv) - (1 - alpha)
  if (f(0) >= 0) return(0)
  stats::uniroot(f, c(0, 5))$root
}


#' Estimate or resolve the overall reliability lambda
#'
#' \code{estimate_lambda()} is a practical reliability resolver. It estimates
#' lambda from measurement-model posterior standard deviations when available,
#' from an explicitly supplied repeated/proxy measure when requested, from a
#' user-supplied lookup table when available, and otherwise returns the
#' convention default \code{c(low = 0.85, mid = 0.90, high = 0.95)}.
#'
#' It deliberately does not estimate lambda from a single point-estimate panel
#' series. With persistent measurement error, a lone observed series cannot
#' distinguish persistent signal from persistent noise without additional
#' structure.
#'
#' @param data Optional data frame.
#' @param variable Character variable name in \code{data}. If \code{data} is
#'   omitted, provide \code{x}.
#' @param x Numeric point-estimate vector.
#' @param posterior_sd Numeric posterior SD vector, or a column name in
#'   \code{data}. If omitted and \code{variable_sd} exists in \code{data}, that
#'   column is used automatically.
#' @param replicate Numeric repeated/proxy measure, or a column name in
#'   \code{data}. Errors in the two measures must be independent conditional on
#'   the latent construct.
#' @param replicate_lambda Known reliability of \code{replicate}. If omitted,
#'   the repeated-measures estimate assumes equal reliability across the two
#'   measures.
#' @param lookup Optional lookup table or named vector/list with known
#'   reliabilities. Data frames should have a \code{variable} column plus either
#'   \code{low/mid/high}, \code{lambda_min/lambda/lambda_max}, or \code{lambda}.
#' @param default Named numeric vector returned when no estimating route is
#'   available. Default is \code{c(low = 0.85, mid = 0.90, high = 0.95)}.
#' @param margin Half-width used to turn a point reliability into a range.
#' @param quiet If \code{FALSE}, print the chosen source.
#' @return An object of class \code{ferobust_lambda}.
#' @export
estimate_lambda <- function(data = NULL, variable = NULL, x = NULL,
                            posterior_sd = NULL,
                            replicate = NULL,
                            replicate_lambda = NULL,
                            lookup = NULL,
                            default = c(low = 0.85, mid = 0.90, high = 0.95),
                            margin = 0.03,
                            quiet = FALSE) {
  if (!is.null(data)) {
    data <- as.data.frame(data)
  }

  if (!is.null(data) && !is.null(variable)) {
    if (!variable %in% names(data)) {
      stop(sprintf("Variable '%s' is not in data.", variable), call. = FALSE)
    }
    x <- data[[variable]]

    if (is.null(posterior_sd)) {
      sd_name <- paste0(variable, "_sd")
      lo_name <- paste0(variable, "_codelow")
      hi_name <- paste0(variable, "_codehigh")
      if (sd_name %in% names(data)) {
        posterior_sd <- sd_name
      } else if (all(c(lo_name, hi_name) %in% names(data))) {
        # V-Dem and similar measurement models often ship a 68% credible
        # interval (codelow, codehigh) instead of a posterior SD column. The
        # half-width of a one-standard-deviation (68%) interval is the
        # posterior SD, so we derive it directly.
        posterior_sd <- (suppressWarnings(as.numeric(data[[hi_name]])) -
                         suppressWarnings(as.numeric(data[[lo_name]]))) / 2
      }
    }
  }

  if (!is.null(data) && is.character(posterior_sd) &&
      length(posterior_sd) == 1 && posterior_sd %in% names(data)) {
    posterior_sd <- data[[posterior_sd]]
  }

  if (!is.null(data) && is.character(replicate) &&
      length(replicate) == 1 && replicate %in% names(data)) {
    replicate <- data[[replicate]]
  }

  if (!is.null(x) && !is.null(posterior_sd)) {
    out <- lambda_from_posterior_sd(x, posterior_sd, margin = margin)
    if (!quiet) print(out)
    return(out)
  }

  if (!is.null(x) && !is.null(replicate)) {
    out <- lambda_from_replicate(x, replicate,
                                 replicate_lambda = replicate_lambda,
                                 margin = margin)
    if (!quiet) print(out)
    return(out)
  }

  if (!is.null(lookup) && !is.null(variable)) {
    out <- lambda_from_lookup(lookup, variable, margin = margin)
    if (!is.null(out)) {
      if (!quiet) print(out)
      return(out)
    }
  }

  out <- make_lambda_result(
    reliability = normalize_reliability(default),
    method = "default",
    source = "convention",
    note = "No posterior SD, repeated measure, or lookup match was available; using the convention default."
  )
  if (!quiet) print(out)
  out
}


#' @rdname estimate_lambda
#' @export
resolve_lambda <- estimate_lambda


#' Estimate lambda from posterior standard deviations
#'
#' @param x Numeric point-estimate vector.
#' @param posterior_sd Numeric posterior SD vector.
#' @param margin Half-width used to form a reliability range around the point
#'   estimate.
#' @return An object of class \code{ferobust_lambda}.
#' @export
lambda_from_posterior_sd <- function(x, posterior_sd, margin = 0.03) {
  x <- suppressWarnings(as.numeric(x))
  posterior_sd <- suppressWarnings(as.numeric(posterior_sd))
  ok <- is.finite(x) & is.finite(posterior_sd)
  if (sum(ok) < 2) {
    stop("Need at least two complete observations to estimate lambda.",
         call. = FALSE)
  }

  total_var <- stats::var(x[ok])
  if (!is.finite(total_var) || total_var <= 0) {
    stop("The point-estimate variance is zero; lambda is undefined.",
         call. = FALSE)
  }

  error_var <- mean(posterior_sd[ok]^2)
  lambda <- clamp_lambda(1 - error_var / total_var)
  make_lambda_result(
    reliability = reliability_range(lambda, margin),
    lambda = lambda,
    method = "posterior_sd",
    source = "measurement_model",
    n_obs = sum(ok),
    note = "Computed as 1 - mean(posterior variance) / var(point estimate). This captures model-reported measurement uncertainty, not conceptual validity error."
  )
}


#' Estimate lambda from a repeated or proxy measure
#'
#' @param x Numeric point-estimate vector.
#' @param replicate Numeric repeated/proxy measure of the same latent construct.
#' @param replicate_lambda Known reliability of the replicate. If omitted, equal
#'   reliability of the two measures is assumed.
#' @param margin Half-width used to form a reliability range around the point
#'   estimate.
#' @return An object of class \code{ferobust_lambda}.
#' @export
lambda_from_replicate <- function(x, replicate, replicate_lambda = NULL,
                                  margin = 0.05) {
  x <- suppressWarnings(as.numeric(x))
  replicate <- suppressWarnings(as.numeric(replicate))
  ok <- is.finite(x) & is.finite(replicate)
  if (sum(ok) < 3) {
    stop("Need at least three complete paired observations to estimate lambda.",
         call. = FALSE)
  }

  r <- stats::cor(x[ok], replicate[ok])
  if (!is.finite(r)) {
    stop("The repeated-measure correlation is undefined.", call. = FALSE)
  }

  if (is.null(replicate_lambda)) {
    lambda <- abs(r)
    method <- "replicate_equal_reliability"
    note <- "Computed as abs(correlation) under independent errors and equal reliability of the two measures."
  } else {
    if (!is.finite(replicate_lambda) || replicate_lambda <= 0 ||
        replicate_lambda > 1) {
      stop("replicate_lambda must be in (0, 1].", call. = FALSE)
    }
    lambda <- r^2 / replicate_lambda
    method <- "replicate_known_reliability"
    note <- "Computed as cor(x, z)^2 / lambda_z under independent errors and a known reliability for the repeated measure."
  }

  lambda <- clamp_lambda(lambda)
  make_lambda_result(
    reliability = reliability_range(lambda, margin),
    lambda = lambda,
    method = method,
    source = "repeated_measure",
    n_obs = sum(ok),
    correlation = r,
    note = note
  )
}


#' Print a resolved lambda
#'
#' @param x A \code{ferobust_lambda} object.
#' @param ... Ignored.
#' @export
print.ferobust_lambda <- function(x, ...) {
  r <- x$reliability
  cat("ferobust lambda\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  cat(sprintf("method: %s\n", x$method))
  cat(sprintf("source: %s\n", x$source))
  cat(sprintf("lambda: [%.3f, %.3f, %.3f]\n",
              unname(r[1]), unname(r[2]), unname(r[3])))
  if (!is.null(x$n_obs)) cat(sprintf("n:      %d\n", x$n_obs))
  if (!is.null(x$correlation)) {
    cat(sprintf("cor:    %.3f\n", x$correlation))
  }
  if (!is.null(x$note)) cat(sprintf("note:   %s\n", x$note))
  invisible(x)
}


lambda_from_lookup <- function(lookup, variable, margin = 0.03) {
  if (is.data.frame(lookup)) {
    if (!"variable" %in% names(lookup)) {
      stop("lookup data frames must include a 'variable' column.",
           call. = FALSE)
    }
    row <- lookup[lookup$variable == variable, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    row <- row[1, , drop = FALSE]

    low_col <- intersect(c("low", "lambda_low", "lambda_min", "min"),
                         names(row))
    mid_col <- intersect(c("mid", "lambda_mid", "lambda", "estimate"),
                         names(row))
    high_col <- intersect(c("high", "lambda_high", "lambda_max", "max"),
                          names(row))

    if (length(low_col) && length(mid_col) && length(high_col)) {
      reliability <- c(low = row[[low_col[1]]],
                       mid = row[[mid_col[1]]],
                       high = row[[high_col[1]]])
    } else if (length(mid_col)) {
      reliability <- reliability_range(row[[mid_col[1]]], margin)
    } else {
      stop("lookup match needs low/mid/high columns or a lambda column.",
           call. = FALSE)
    }
  } else if (is.list(lookup) && !is.null(lookup[[variable]])) {
    reliability <- normalize_reliability(lookup[[variable]], margin = margin)
  } else if (is.numeric(lookup) && !is.null(names(lookup)) &&
             variable %in% names(lookup)) {
    reliability <- reliability_range(lookup[[variable]], margin)
  } else {
    return(NULL)
  }

  make_lambda_result(
    reliability = normalize_reliability(reliability, margin = margin),
    method = "lookup",
    source = "measurement_literature",
    note = sprintf("Matched '%s' in the user-supplied reliability lookup.",
                   variable)
  )
}


make_lambda_result <- function(reliability, method, source, note = NULL,
                               lambda = NULL, n_obs = NULL,
                               correlation = NULL) {
  structure(
    list(reliability = normalize_reliability(reliability),
         lambda = lambda,
         method = method,
         source = source,
         note = note,
         n_obs = n_obs,
         correlation = correlation),
    class = "ferobust_lambda"
  )
}


normalize_reliability <- function(x, margin = 0.03) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 1) {
    return(reliability_range(x, margin))
  }
  if (length(x) == 2) {
    x <- c(low = min(x), mid = mean(x), high = max(x))
  } else if (length(x) >= 3) {
    x <- x[seq_len(3)]
    x <- c(low = min(x), mid = stats::median(x), high = max(x))
  } else {
    stop("Reliability must have length 1, 2, or 3.", call. = FALSE)
  }
  clamp_lambda(x)
}


reliability_range <- function(lambda, margin = 0.03) {
  lambda <- clamp_lambda(lambda)
  c(low = clamp_lambda(lambda - margin),
    mid = lambda,
    high = clamp_lambda(lambda + margin))
}


clamp_lambda <- function(lambda) {
  pmin(pmax(lambda, 0), 0.999)
}


#' Reliability of a first-differenced regressor
#'
#' Differencing a persistent series with classical noise lowers reliability.
#' For an AR(1) signal with persistence \code{phi} and classical noise, the
#' reliability of the first difference is
#' \code{lambda(1 - phi) / (lambda(1 - phi) + (1 - lambda))}. Use this when the
#' regressor is a change/difference but only a \emph{levels} reliability
#' \code{lambda} is available (the Andersen-Doucette case): apply it, then pass
#' the result and the change regressor's own ICC to \code{\link{within_reliability}}.
#'
#' @param lambda Levels reliability (scalar or vector).
#' @param phi AR(1) persistence of the level signal.
#' @return Reliability of the first difference.
#' @examples
#' lambda_delta(lambda = c(0.85, 0.95), phi = 0.92)
#' @export
lambda_delta <- function(lambda, phi) {
  (lambda * (1 - phi)) / (lambda * (1 - phi) + (1 - lambda))
}


#' Within-reliability under serially-correlated measurement error
#'
#' The i.i.d. within-reliability of \code{\link{within_reliability}} puts all
#' noise in the within dimension and is a worst case. If a share \code{pi} of
#' the measurement-error variance is persistent (between-unit), the
#' within-reliability rises to
#' \code{(lambda - icc_hat + pi (1 - lambda)) / (1 - icc_hat)}, recovering the
#' i.i.d. formula at \code{pi = 0} and approaching 1 as \code{pi -> 1}. Sweep
#' \code{pi} to check whether a verdict survives persistent coder error.
#'
#' @param lambda Overall reliability.
#' @param icc_hat Empirical ICC of the observed regressor.
#' @param pi Persistent share of the measurement-error variance (the error's
#'   own ICC), in \code{[0, 1]}.
#' @return Within-reliability accounting for persistent error.
#' @export
lambda_w_serial <- function(lambda, icc_hat, pi = 0) {
  pmax((lambda - icc_hat + pi * (1 - lambda)) / (1 - icc_hat), 1e-3)
}
