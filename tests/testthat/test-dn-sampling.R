# Oracle: fixtures/dn-reference.rds, written by newmark 1.7.1 in its own
# process from the inputs of its lib/tests/testthat/test-sampleSaCorr.R
# (SHA-256 5e3e37dbe624f973...), copied below, and from the roxygen
# examples of its rhoBJ.R and buildQSpline.R.
reference <- function() readRDS(test_path("fixtures", "dn-reference.rds"))

samplerQuantiles <- function(Tn, mu) {
  p <- c("0.05", "0.16", "0.50", "0.84", "0.95", "mean")
  f <- c(0.40, 0.62, 0.86, 1.22, 1.70, 1.00)
  data.table::rbindlist(Map(function(Tj, muj) data.table::data.table(Tn = Tj, p = p, Sa = muj * f), Tn, mu))
}

samplerMixture <- function() {
  OUT <- data.table::data.table(Tn = c(0, 0.4, 1), k = "oq", sigmaLn = c(0.20, 0.45, 0.75), w = 1)
  OUT[, muLn := log(c(0.80, 0.55, 0.42)) - sigmaLn^2 / 2][]
}

# The reference-period correlation of newmark 1.7.1 sampleSaCorr() as a matrix.
referenceCorrelation <- function(rho) {
  C <- diag(length(rho) + 1L)
  C[1, -1] <- C[-1, 1] <- rho
  C
}

asMatrix <- function(sample) matrix(sample$Sa, ncol = data.table::uniqueN(sample$Tn))

test_that("correlateResiduals reproduces the 1.7.1 coefficients and keeps its conventions", {
  expectReproduces(c(correlateResiduals(Tn = 0.3, T0 = 0.01), correlateResiduals(Tn = 0.5, T0 = 0.5),
                     correlateResiduals(Tn = 0.3, T0 = 0.5), correlateResiduals(Tn = 1, T0 = 0.01)),
                   reference()$rho)
  expectReproduces(correlateResiduals(Tn = c(0.3, 1), T0 = 0.01), reference()$rho[c(1L, 4L)])
  expect_identical(correlateResiduals(Tn = 0.5, T0 = 0.3), correlateResiduals(Tn = 0.3, T0 = 0.5))
  expect_identical(correlateResiduals(Tn = 0, T0 = 1), correlateResiduals(Tn = 0.01, T0 = 1))
  expect_warning(correlateResiduals(Tn = 20, T0 = 1), "0.01 to 10 s")
  expect_error(correlateResiduals(Tn = c(0.1, 0.2), T0 = c(0.1, 0.2, 0.3)), "equal length")
})

test_that("sampleSpectra reproduces the 1.7.1 quantile draws in any row order", {
  R <- reference()
  Q <- samplerQuantiles(c(0, 0.6, 0.9), c(0.80, 0.55, 0.42))
  set.seed(42)
  OUT <- sampleSpectra(Q, Tn = c(0, 0.6, 0.9), correlation = referenceCorrelation(c(0.62, 0.38)), NS = 47L, zTrunc = NA)
  expectReproduces(asMatrix(OUT), R$q42)
  expect_identical(OUT$realization, rep(1:47, 3L))
  expect_equal(OUT[, mean(Sa), by = Tn]$V1, c(0.80, 0.55, 0.42), tolerance = 1e-13)
  set.seed(42)
  expect_identical(sampleSpectra(Q[rev(seq_len(nrow(Q)))], Tn = c(0, 0.6, 0.9),
                                 correlation = referenceCorrelation(c(0.62, 0.38)), NS = 47L, zTrunc = NA), OUT)
  set.seed(43)
  expectReproduces(asMatrix(sampleSpectra(Q, Tn = c(0.9, 0, 0.6), correlation = referenceCorrelation(c(0.38, 0.62)),
                                          NS = 47L, zTrunc = NA)), R$q43)
  C <- matrix(c(1.00, 0.62, 0.38, 0.62, 1.00, 0.54, 0.38, 0.54, 1.00), nrow = 3L, byrow = TRUE)
  set.seed(44)
  expectReproduces(asMatrix(sampleSpectra(samplerQuantiles(c(0, 1.3 * 0.4, 1.5 * 0.4), c(0.80, 0.55, 0.42)),
                                          Tn = c(0, 1.3 * 0.4, 1.5 * 0.4), correlation = C, NS = 47L, zTrunc = NA)), R$full44)
  set.seed(51)
  expectReproduces(asMatrix(sampleSpectra(samplerQuantiles(c(0, 1), c(0.80, 0.42)), Tn = c(0, 1),
                                          correlation = referenceCorrelation(0.5), NS = 7L, zTrunc = NA)), R$q51)
})

test_that("analytic mixtures, atoms and truncation reproduce the 1.7.1 draws", {
  R <- reference()
  set.seed(45)
  expectReproduces(asMatrix(sampleSpectra(samplerMixture(), Tn = c(0, 0.52, 0.60), correlation = referenceCorrelation(c(0.62, 0.38)),
                                          NS = 47L, zTrunc = NA)), R$a45)
  Atoms <- data.table::data.table(Tn = rep(c(0, 1), each = 2L), k = rep(c("a", "b"), 2L),
                                  muLn = rep(log(c(0.2, 0.8)), 2L), sigmaLn = 0, w = 0.5)
  set.seed(47)
  expectReproduces(asMatrix(sampleSpectra(Atoms, Tn = c(0, 1), correlation = referenceCorrelation(0.5), NS = 1L, zTrunc = NA)),
                   R$atom47)
  Mixture <- data.table::data.table(Tn = c(0, 0, 1, 1), k = c("a", "b", "a", "b"), muLn = c(0.20, 0.35, -0.90, -0.70),
                                    sigmaLn = c(0.80, 0.60, 0.75, 0.55), w = 0.5)
  set.seed(20260806)
  OUT <- sampleSpectra(Mixture, Tn = c(0, 1), correlation = referenceCorrelation(correlateResiduals(Tn = 1, T0 = 0.01)),
                       NS = 4000, zTrunc = 2)
  expectReproduces(asMatrix(OUT), R$trunc)
  Truncated <- Mixture[, .(mean = sum(w * exp(muLn + sigmaLn^2 / 2) * (pnorm(2 - sigmaLn) - pnorm(-2 - sigmaLn)) /
                                        (pnorm(2) - pnorm(-2)))), by = Tn]$mean
  expect_equal(OUT[, mean(Sa), by = Tn]$V1, Truncated, tolerance = 1e-12)
})

test_that("the inverse of ten equal atoms at p = 0.9 is the ninth, as added in double", {
  # Ten weights of 0.1 add to 0.99999999999999989 in double, and nine of them
  # normalized by that sum reach 0.9; normalized by the long-double sum, 1,
  # they add to 0.89999999999999991 and would select the tenth atom.
  expect_identical(.invertMixture(p = 0.9, muLn = (1:10) / 10, sigmaLn = rep(0, 10L), w = rep(0.1, 10L)), 9 / 10)
})

test_that("the recentring mean stays finite when the double sum overflows", {
  x <- c(1e308, 1e308)
  expect_true(is.finite(.portableMean(x)))
  expect_identical(.portableMean(x), mean(x))
})

test_that("sampleIntensities reproduces the paired intensities of 1.7.1", {
  R <- reference()
  set.seed(44)
  OUT <- sampleIntensities(samplerQuantiles(c(0, 1.3 * 0.4, 1.5 * 0.4), c(0.80, 0.55, 0.42)), Ts = 0.4, NS = 47L, zTrunc = NA)
  expectReproduces(as.list(OUT[, .(PGA, Sa13, Sa15, PGV, AI)]), R$scenario44)
  set.seed(46)
  OUT <- sampleIntensities(samplerMixture(), Ts = 0.4, NS = 47L, zTrunc = NA)
  expectReproduces(as.list(OUT[, .(PGA, Sa13, Sa15, PGV, AI)]), R$scenario46)
  expect_identical(OUT$PGV, estimatePGV(OUT$PGA))
  expect_identical(OUT$AI, estimateAI(OUT$PGA))
})

test_that("invalid requests fail before consuming the random stream", {
  set.seed(49)
  Seed <- .Random.seed
  Incomplete <- samplerMixture()[, w := NULL]
  expect_error(sampleSpectra(Incomplete, Tn = c(0, 1), correlation = referenceCorrelation(0.5), NS = 5L, zTrunc = NA),
               "incomplete schema")
  expect_error(sampleSpectra(samplerQuantiles(c(0, 1), c(0.80, 0.42))[p %in% c("0.16", "0.84", "mean")], Tn = c(0, 1),
                             correlation = referenceCorrelation(0.5), NS = 5L, zTrunc = NA), "3 or more quantile rows")
  expect_error(sampleSpectra(samplerQuantiles(c(0, 1), c(0.80, 0.42)), Tn = c(0, 1), correlation = referenceCorrelation(0.5),
                             NS = 5L, zTrunc = 2), "analytic mixture")
  expect_error(sampleSpectra(samplerMixture(), Tn = 0, correlation = diag(1), NS = 10, zTrunc = -1), "zTrunc")
  expect_error(sampleSpectra(samplerMixture(), Tn = c(0, 1), correlation = matrix(c(1, 0.5, 0.4, 1), 2), NS = 5L, zTrunc = NA),
               "symmetric")
  expect_identical(.Random.seed, Seed)
})

test_that("invalid quantiles and means fail before consuming the random stream", {
  Q <- samplerQuantiles(c(0, 1), c(0.80, 0.42))
  draw <- function(spectrum) {
    sampleSpectra(spectrum, Tn = c(0, 1), correlation = referenceCorrelation(0.5), NS = 5L, zTrunc = NA)
  }
  set.seed(49)
  Seed <- .Random.seed
  expect_error(draw(data.table::copy(Q)[p == "mean", Sa := -Sa]), "mean row")
  expect_error(draw(data.table::copy(Q)[p == "mean" & Tn == 1, Sa := NA_real_]), "mean row")
  expect_error(draw(data.table::copy(Q)[p == "0.84", p := "1.84"]), "strictly between 0 and 1")
  expect_error(draw(data.table::copy(Q)[2L, p := NA_character_]), "strictly between 0 and 1")
  expect_error(draw(data.table::copy(Q)[1L, Sa := NA_real_]), "finite positive accelerations")
  expect_identical(.Random.seed, Seed)
  # A period without a mean row is recentered on its sample mean, which
  # leaves its draws unscaled.
  set.seed(51)
  With <- asMatrix(draw(Q))
  set.seed(51)
  Without <- asMatrix(draw(Q[p != "mean"]))
  expect_identical(With, sweep(Without, 2L, Q[p == "mean", Sa] / apply(Without, 2L, mean), "*"))
})

test_that("buildQuantileSpline rejects probabilities outside (0, 1) and non-positive accelerations", {
  buildSpline <- function(p, Sa) buildQuantileSpline(data.frame(p = p, Sa = Sa))
  expect_error(buildSpline(c(0, 0.50, 0.84), c(0.3, 0.5, 0.8)), "strictly between 0 and 1")
  expect_error(buildSpline(c(0.16, 0.50, 1), c(0.3, 0.5, 0.8)), "strictly between 0 and 1")
  expect_error(buildSpline(c("0.16", "high", "0.84"), c(0.3, 0.5, 0.8)), "strictly between 0 and 1")
  expect_error(buildSpline(c(0.16, NA, 0.50, 0.84), c(0.3, 0.4, 0.5, 0.8)), "strictly between 0 and 1")
  expect_error(buildSpline(c(0.16, 0.50, 0.84), c(0.3, 0, 0.8)), "finite positive accelerations")
  expect_error(buildSpline(c(0.16, 0.50, 0.84), c(0.3, Inf, 0.8)), "finite positive accelerations")
})

test_that("buildQuantileSpline and interpolateSpectrum reproduce 1.7.1", {
  R <- reference()
  Q <- buildQuantileSpline(data.frame(p = c(0.16, 0.50, 0.84), Sa = c(0.3, 0.5, 0.8)))
  expectReproduces(Q(c(0, 0.1, 0.16, 0.3, 0.5, 0.84, 0.9, 1)), R$spline)
  expect_warning(buildQuantileSpline(data.frame(p = c(0.16, 0.5, 0.5, 0.84), Sa = c(0.3, 0.5, 0.5, 0.8))), "duplicate p rows")
  expect_error(buildQuantileSpline(data.frame(p = c(0.16, 0.5, 0.84), Sa = c(0.3, 0.5, 0.8), Tn = c(0, 0, 1))), "several values")
  Spectrum <- data.frame(Tn = rep(c(0.2, 0.8), each = 3L), p = rep(c("0.16", "0.84", "mean"), 2L), Sa = c(1, 2, 1.5, 4, 8, 6))
  expectReproduces(as.list(interpolateSpectrum(Spectrum, Tn = 0.4)), R$interpolate)
  # Beyond the tabulated periods the end ordinates hold, up to exp(log(x)).
  expect_equal(interpolateSpectrum(Spectrum, Tn = 2)$Sa, c(4, 8, 6), tolerance = 1e-14)
})

test_that("buildKyGrid spans 0.01 g to the mean PGA or 0.80 g", {
  expectReproduces(buildKyGrid(PGA = samplerMixture()[Tn == 0, sum(w * exp(muLn + sigmaLn^2 / 2))], n = 30L),
                   reference()$kyAnalytic)
  Grid <- buildKyGrid(PGA = 1.2, n = 30L)
  expect_equal(range(Grid), c(0.01, 1.2), tolerance = 1e-15)
  expect_identical(buildKyGrid(PGA = 0.3, n = 30L), buildKyGrid(PGA = 0.8, n = 30L))
  # One point is the lower bound alone, whatever the PGA.
  expect_identical(c(buildKyGrid(PGA = 0.3, n = 1L), buildKyGrid(PGA = 1.2, n = 1L)), rep(exp(log(0.01)), 2L))
  expect_error(buildKyGrid(PGA = 0, n = 30L), "positive")
})
