## Generator for the bundled `panel_demo` dataset.
## Run from the package root:  source("data-raw/panel_demo.R")
set.seed(20260601)
N <- 120L; T <- 25L; beta <- 0.5
mu  <- rnorm(N, 0, 1)                       # unit means -> high ICC (slow moving)
phi <- 0.75; sw <- sqrt(0.20 * (1 - phi^2)) # persistent within-unit AR(1)
sigma_u <- 0.46                             # classical ME -> reliability ~ 0.85
rho <- 0.12                                 # mild positive confounding in pooled
unit <- rep(sprintf("U%03d", 1:N), each = T)
year <- rep(2000:(2000 + T - 1L), times = N)
xstar <- numeric(N * T)
for (i in 1:N) {
  w <- numeric(T); w[1] <- rnorm(1, 0, sw / sqrt(1 - phi^2))
  for (t in 2:T) w[t] <- phi * w[t - 1] + rnorm(1, 0, sw)
  xstar[((i - 1) * T + 1):(i * T)] <- mu[i] + w
}
xsd <- pmax(0.08, rnorm(N * T, sigma_u, 0.10 * sigma_u))   # per-obs posterior SD
x   <- xstar + rnorm(N * T, 0, xsd)                        # observed, mismeasured
y   <- rho * mu[match(unit, sprintf("U%03d", 1:N))] + beta * xstar + rnorm(N * T, 0, 1)
panel_demo <- data.frame(unit, year, y, x, x_sd = xsd,
                         x_codelow = x - xsd, x_codehigh = x + xsd,
                         stringsAsFactors = FALSE)
save(panel_demo, file = "data/panel_demo.rda", compress = "xz")
