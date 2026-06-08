test_that("compute_icc is between/total and exact on a balanced panel", {
  set.seed(1)
  N <- 50; T <- 10
  mu <- rnorm(N, 0, 2); unit <- rep(1:N, each = T)
  x <- rep(mu, each = T) + rnorm(N*T, 0, 1)          # var(between)=4, var(within)=1 => ICC ~ 0.8
  expect_equal(compute_icc(x, unit), 0.8, tolerance = 0.12)
  # equals the obs-weighted between/total decomposition by construction
  m <- tapply(x, unit, mean)
  expect_equal(unname(compute_icc(x, unit)),
               unname(var(m[as.character(unit)]) / var(x)), tolerance = 1e-10)
})

test_that("within_reliability matches the corrected formula and floors at data-inconsistency", {
  expect_equal(within_reliability(0.90, 0.75), (0.90-0.75)/(1-0.75), tolerance = 1e-9)
  expect_equal(within_reliability(0.50, 0.80), 1e-3)   # lambda < icc -> floored
})

test_that("breakdown_reliability solves FE/lambda_w = pooled and is NA for complements", {
  bd <- breakdown_reliability(beta_pooled = 0.184, beta_fe = 0.153, icc_hat = 0.736)
  lw <- within_reliability(bd, 0.736)
  expect_equal(0.153/lw, 0.184, tolerance = 1e-6)
  expect_true(is.na(breakdown_reliability(0.05, 0.18, 0.5)))  # |FE| > |pooled|
})

test_that("lambda_delta and lambda_w_serial recover their boundary cases", {
  expect_equal(lambda_delta(0.9, 0), 0.9, tolerance = 1e-9)   # phi=0 -> levels reliability
  expect_lt(lambda_delta(0.9, 0.95), 0.5)                     # high persistence -> low
  expect_equal(lambda_w_serial(0.9, 0.75, 0), within_reliability(0.9, 0.75), tolerance = 1e-9)
  expect_gt(lambda_w_serial(0.9, 0.75, 0.5), within_reliability(0.9, 0.75))
})

test_that("within_reliability_frontier floor is in [0,1) and falls with psi_max", {
  set.seed(2); N<-40; T<-15; unit<-rep(1:N,each=T); time<-rep(1:T,N)
  s <- as.vector(replicate(N, as.numeric(arima.sim(list(ar=0.8), T))))
  x <- s + rnorm(N*T, 0, 0.5)
  fr <- within_reliability_frontier(x, unit, time, psi_max=c(0,0.5,0.9), max_lag=4)
  expect_true(all(fr$lambda_w_lower >= 0 & fr$lambda_w_lower < 1))
  expect_gte(fr$lambda_w_lower[1], fr$lambda_w_lower[3])      # weaker assumption -> tighter bound
})
