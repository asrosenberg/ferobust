#' Synthetic panel for demonstrating the within-reliability diagnostic
#'
#' A simulated unit-year panel with a slow-moving, mismeasured regressor and its
#' measurement-model posterior uncertainty, built so that pooled OLS and fixed
#' effects shrink in the way the diagnostic targets. Running the workflow on it
#' returns a \dQuote{rescue}: the corrected bounds exclude zero even though the
#' fixed-effects estimate is attenuated.
#'
#' The data were generated from \eqn{y_{it} = \alpha_i + \beta x^*_{it} +
#' \varepsilon_{it}} with true \eqn{\beta = 0.5}. The latent regressor \eqn{x^*}
#' is persistent (high intraclass correlation), the observed \code{x} adds
#' classical measurement error with per-observation posterior standard deviation
#' \code{x_sd}, and the unit effects are mildly correlated with \eqn{\bar x^*_i}
#' so pooled OLS is inflated. See \code{data-raw/panel_demo.R} for the generator.
#'
#' @format A data frame with 3000 rows (120 units by 25 years) and 7 variables:
#' \describe{
#'   \item{unit}{Unit identifier, \code{"U001"} to \code{"U120"}.}
#'   \item{year}{Calendar year, 2000 to 2024.}
#'   \item{y}{Outcome.}
#'   \item{x}{Observed (mismeasured) slow-moving regressor.}
#'   \item{x_sd}{Posterior standard deviation of \code{x} from the measurement model.}
#'   \item{x_codelow}{Lower bound of the 68 percent credible interval, \code{x - x_sd}.}
#'   \item{x_codehigh}{Upper bound of the 68 percent credible interval, \code{x + x_sd}.}
#' }
#' @examples
#' if (requireNamespace("fixest", quietly = TRUE)) {
#'   # reliability = "auto" reads the x_sd posterior column automatically
#'   audit(y ~ x | unit + year, data = panel_demo,
#'         key_var = "x", reliability = "auto", re = FALSE)
#' }
"panel_demo"
