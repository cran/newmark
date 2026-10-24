# Cases of test-kmax-mixture.R and test-kmax-quantiles.R of newmark 1.7.1
# (lib/tests/testthat at commit f7d11bb) rewritten for summarizeKmax(): every
# realization of a model has Dn = root / ky, so its kmax at Da is root / Da.
LEVELS <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)

mixtureRoots <- function(roots, samples, Da, ky = c(0.05, 0.10, 0.20, 0.30, 0.40)) {
  rbindlist(lapply(names(roots), function(model) {
    solveKmax(ky = ky, Dn = matrix(roots[[model]] / ky, nrow = samples[[model]], ncol = length(ky), byrow = TRUE),
              PGA = NULL, model = model, Da = Da)
  }))
}

test_that("the pool saturates the realizations whose ceiling the target exceeds", {
  Roots <- mixtureRoots(c(AM88 = 0.10, BT07 = 0.30), c(AM88 = 20L, BT07 = 20L), Da = c(0.5, 1, 4))
  OUT <- summarizeKmax(Roots, weights = c(AM88 = 0.5, BT07 = 0.5), levels = LEVELS)
  # Da = 4 exceeds the AM88 ceiling (0.10 / 0.05 = 2): that half of the pool
  # is null demand.
  expect_equal(OUT[level == "mean", kmax], c(0.4, 0.2, 0.0375))
  expect_equal(OUT[level == "0.95", kmax], 0.30 / c(0.5, 1, 4))
  expect_identical(unique(OUT$level), c("mean", "0.05", "0.10", "0.16", "0.50", "0.84", "0.90", "0.95"))
  expect_identical(unique(OUT$n), 40L)
})

test_that("one model keeps its own roots", {
  OUT <- summarizeKmax(mixtureRoots(c(AM88 = 0.10), c(AM88 = 23L), Da = c(0.5, 1, 4)), weights = c(AM88 = 1),
                       levels = LEVELS)
  expect_equal(OUT[level == "mean", kmax], c(0.2, 0.1, 0))
  expect_equal(OUT[level == "0.95", kmax], c(0.2, 0.1, 0))
})

test_that("the pool mixes the model distributions with their normalized weights", {
  Roots <- mixtureRoots(c(AM88 = 0.10, BT07 = 0.30), c(AM88 = 37L, BT07 = 37L), Da = 1)
  OUT <- summarizeKmax(Roots, weights = c(AM88 = 0.8, BT07 = 0.2), levels = LEVELS)
  expect_equal(OUT[level == "mean", kmax], 0.14)
  expect_equal(OUT[level == "0.50", kmax], 0.10)
  expect_equal(OUT[level == "0.95", kmax], 0.30)
})

test_that("the summary ignores realization order and weight scale", {
  ky <- c(0.05, 0.10, 0.20, 0.40, 0.80, 1.00)
  curves <- function(roots) matrix(roots, nrow = length(roots), ncol = length(ky)) / matrix(ky, length(roots), length(ky), byrow = TRUE)
  Baseline <- summarizeKmax(rbind(solveKmax(ky = ky, Dn = curves(c(0.10, 0.20, 0.30, 0.40)), PGA = NULL, model = "AM88", Da = 1),
                                  solveKmax(ky = ky, Dn = curves(c(0.50, 0.60, 0.70, 0.80)), PGA = NULL, model = "BT07", Da = 1)),
                            weights = c(AM88 = 0.5, BT07 = 0.5), levels = c(0.25, 0.50, 0.75))
  Candidate <- summarizeKmax(rbind(solveKmax(ky = ky, Dn = curves(c(0.80, 0.70, 0.60, 0.50)), PGA = NULL, model = "BT07", Da = 1),
                                   solveKmax(ky = ky, Dn = curves(c(0.10, 0.20, 0.30, 0.40)), PGA = NULL, model = "AM88", Da = 1)),
                             weights = c(BT07 = 5, AM88 = 5), levels = c(0.25, 0.50, 0.75))
  expect_equal(Candidate, Baseline)
  expect_equal(Baseline[level == "mean", kmax], 0.45)
  expect_equal(Baseline[level != "mean", kmax], c(0.35, 0.50, 0.65))
})

test_that("one model of unit weight gives the type-7 sample quantiles", {
  ky <- c(0.05, 0.10, 0.20, 0.40, 0.80, 1.00)
  Root <- c(0.11, 0.27, 0.19, 0.43, 0.31, 0.24, 0.38)
  Dn <- matrix(Root, length(Root), length(ky)) / matrix(ky, length(Root), length(ky), byrow = TRUE)
  Roots <- solveKmax(ky = ky, Dn = Dn, PGA = NULL, model = "BT07", Da = 1)
  OUT <- summarizeKmax(Roots, weights = c(BT07 = 1), levels = LEVELS)
  expect_equal(OUT[level != "mean", kmax], unname(stats::quantile(Roots$kmax, LEVELS, type = 7)))
  expect_equal(OUT[level == "mean", kmax], mean(Roots$kmax))
})

test_that("pools with missing models or unshared realizations fail", {
  Roots <- mixtureRoots(c(AM88 = 0.10), c(AM88 = 3L), Da = 1)
  expect_error(summarizeKmax(Roots, weights = c(AM88 = 0.5, BT07 = 0.5), levels = LEVELS), "missing positive-weight models")
  Roots <- mixtureRoots(c(AM88 = 0.10, BT07 = 0.30), c(AM88 = 3L, BT07 = 2L), Da = 1)
  expect_error(summarizeKmax(Roots, weights = c(AM88 = 0.5, BT07 = 0.5), levels = LEVELS), "shared realization coverage")
  expect_error(summarizeKmax(Roots, weights = c(AM88 = 1), levels = c(0.1, 0.101)), "collide")
  # Coverage must hold within every Da, not only over the whole table: here
  # both models share realizations 1 to 3 but at different Da.
  Roots <- rbind(mixtureRoots(c(AM88 = 0.10), c(AM88 = 3L), Da = 1), mixtureRoots(c(BT07 = 0.30), c(BT07 = 3L), Da = 2))
  expect_error(summarizeKmax(Roots, weights = c(AM88 = 0.5, BT07 = 0.5), levels = LEVELS), "within every Da")
  Roots <- mixtureRoots(c(AM88 = 0.10), c(AM88 = 3L), Da = 1)
  expect_error(summarizeKmax(rbind(Roots, Roots[1L]), weights = c(AM88 = 1), levels = LEVELS), "one row per model")
})

test_that("a sampled PGA leaves invertible curves unchanged and no random number is drawn", {
  ky <- c(0.05, 0.10, 0.20, 0.30, 0.40)
  Dn <- matrix(0.10 / ky, nrow = 20L, ncol = 5L, byrow = TRUE)
  set.seed(730)
  Seed <- .Random.seed
  Without <- solveKmax(ky = ky, Dn = Dn, PGA = NULL, model = "AM88", Da = c(0.5, 1, 4))
  With <- solveKmax(ky = ky, Dn = Dn, PGA = rep(0.75, 20L), model = "AM88", Da = c(0.5, 1, 4))
  summarizeKmax(With, weights = c(AM88 = 1), levels = LEVELS)
  expect_identical(With, Without)
  expect_identical(.Random.seed, Seed)
})

test_that("tied pools keep full precision in the general summary and text precision in the compatibility route", {
  Roots <- data.table(IDn = "BT07", realization = 1:3, Da = 1, kmax = c(0, 0, 0.1 + 0.2), beyond = FALSE)
  General <- summarizeKmax(Roots, weights = c(BT07 = 1), levels = 0.95)[level == "0.95", kmax]
  Legacy <- .textWeightedQuantile(Roots$kmax, rep(1, 3L), 0.95)
  expect_false(identical(General, Legacy))
  expect_equal(General, Legacy)
  expect_identical(Legacy, .weightedQuantile(c(0, 0, 0.3), rep(1, 3L), 0.95))
})

test_that("the beyond share of a fully beyond family pool stays at one", {
  # The family weights of buildKmaxTable() normalized again by the pool: added
  # in extended precision their share exceeds 1 and fails the beyond bound.
  Roots <- CJ(IDn = c("AM88", "JB07", "SR08", "BT07"), realization = 1:20, Da = 1)[, `:=`(kmax = seq_len(.N) / 100, beyond = TRUE)]
  OUT <- summarizeKmax(Roots, weights = c(AM88 = 1, JB07 = 1, SR08 = 3, BT07 = 1) / 6, levels = 0.5)
  expect_identical(unique(OUT$beyond), 1)
})
