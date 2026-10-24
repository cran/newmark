# Conformity with test-fitSaF.R of newmark 1.7.1: same fixtures, seeds and
# arguments; the 1.7.1 outputs come from fixtures/siteResponse-reference.rds,
# written by the installed newmark 1.7.1 (build d57b64f) in its own R process,
# never together with newmark 5.
siteReference <- function() readRDS(test_path("fixtures", "siteResponse-reference.rds"))

siteColumns <- function(x) as.list(x[, .(Tn, p, Sa, SaF)])

test_that("integrateSiteSpectrum reproduces the fitSaF quadrature of newmark 1.7.1 exactly", {
  Reference <- siteReference()
  set.seed(920)
  Seed <- .Random.seed
  for (Order in c(32, 64)) {
    OUT <- integrateSiteSpectrum(Reference$uhs, VS30 = 270, Vref = 760, model = "ST20", order = Order)
    expectReproduces(siteColumns(OUT), siteColumns(Reference[[paste0("quadrature", Order)]]))
  }
  expect_identical(.Random.seed, Seed)
  expect_identical(names(OUT), c("Tn", "p", "Sa", "SaF"))
  Monotone <- OUT[p != "mean"][order(as.numeric(p)), .(!is.unsorted(SaF)), by = Tn]
  expect_true(all(Monotone$V1))
})

test_that("integrateSiteSpectrum resolves the steep S2L1W upper tail as newmark 1.7.1", {
  Reference <- siteReference()
  Spectra <- lapply(c(64, 128, 256), function(Order) {
    integrateSiteSpectrum(Reference$steep, VS30 = 180, Vref = 760, model = "ST20", order = Order)
  })
  expectReproduces(siteColumns(Spectra[[1L]]), siteColumns(Reference$steep64))
  expectReproduces(siteColumns(Spectra[[2L]]), siteColumns(Reference$steep128))
  expectReproduces(siteColumns(Spectra[[3L]]), siteColumns(Reference$steep256))
  # Thresholds of test-fitSaF.R of newmark 1.7.1.
  expect_lt(.relativeChange(Spectra[[1L]]$SaF, Spectra[[2L]]$SaF), 2e-3)
  expect_lt(.relativeChange(Spectra[[2L]]$SaF, Spectra[[3L]]$SaF), 5e-4)
})

test_that("refineSiteSpectrum accepts the first order within tolerance with its change", {
  Reference <- siteReference()
  Refined <- refineSiteSpectrum(Reference$uhs, VS30 = 270, Vref = 760, model = "ST20",
                                orders = c(32, 64), tolerance = 2e-3)
  expect_identical(Refined$order, 64L)
  expectReproduces(siteColumns(Refined$spectrum), siteColumns(Reference$quadrature64))
  expectReproducesChange(Refined$change, .relativeChange(Reference$quadrature32$SaF, Reference$quadrature64$SaF))
  expect_error(refineSiteSpectrum(Reference$uhs, VS30 = 270, Vref = 760, model = "ST20",
                                  orders = c(32, 64), tolerance = 1e-12), "did not converge")
})

test_that("the reference condition returns the rock spectrum", {
  Reference <- siteReference()
  OUT <- integrateSiteSpectrum(Reference$uhs, VS30 = 760, Vref = 760, model = "ST20", order = 32)
  Rock <- data.table::copy(Reference$uhs)
  data.table::setorder(Rock, Tn, p)
  expect_identical(OUT$Sa, Rock$Sa)
  expect_identical(OUT$SaF, Rock$Sa)
})

test_that("the quadrature quantile function keeps the knots and repairs folded tails", {
  # Quantiles and assertions of test-fitSaF.R of newmark 1.7.1.
  Quantiles <- data.table::data.table(
    p = c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95),
    Sa = c(0.002172658, 0.002368442, 0.002523095, 0.003518529, 0.005471477, 0.006568370, 0.007305735)
  )
  Raw <- buildQuantileSpline(Quantiles)
  Repaired <- .buildTailedQuantile(Quantiles)
  Grid <- seq(0, 1, length.out = 10001L)
  Interior <- seq(min(Quantiles$p), max(Quantiles$p), length.out = 1001L)
  expect_true(any(diff(Raw(Grid)) < 0))
  expect_gte(min(diff(Repaired(Grid))), -1e-14)
  expect_equal(Repaired(Quantiles$p), log(Quantiles$Sa), tolerance = 1e-14)
  expect_equal(Repaired(Interior), Raw(Interior), tolerance = 1e-14)
  Tied <- data.table::copy(Quantiles)
  Tied[3L, Sa := Tied$Sa[[2L]]]
  expect_equal(.buildTailedQuantile(Tied)(Tied$p), log(Tied$Sa), tolerance = 1e-14)
  expect_gte(min(diff(.buildTailedQuantile(Tied)(Grid))), -1e-14)
  Descending <- data.table::copy(Quantiles)
  Descending[3L, Sa := 0.9 * Descending$Sa[[2L]]]
  expect_error(.buildTailedQuantile(Descending), "non-decreasing", fixed = TRUE)
})

test_that("a mean-only spectrum is the conditional lognormal of evaluateSiteSa", {
  Reference <- siteReference()
  Means <- Reference$uhs[p == "mean"]
  expect_error(integrateSiteSpectrum(Means, VS30 = 270, Vref = 760, model = "ST20", order = 32), "evaluateSiteSa")
  expect_error(simulateSiteSpectrum(Means, VS30 = 270, Vref = 760, model = "ST20", NS = 10), "evaluateSiteSa")
  # The closed form that test-fitSaF.R of newmark 1.7.1 asserts for mean-only
  # input, with its tolerance.
  P <- c(0.05, 0.50, 0.95)
  PGA <- Means[Tn == 0, Sa]
  for (t in Means$Tn) {
    Sa <- Means[Tn == t, Sa]
    Site <- evaluateSiteSa(Sa = Sa, PGA = PGA, Tn = t, VS30 = 270, Vref = 760, model = "ST20")
    Factor <- evaluateSiteFactor(PGA = PGA, Tn = t, VS30 = 270, Vref = 760, model = "ST20")
    expect_equal(exp(Site$muLnSaF + Site$sigmaLnSaF * stats::qnorm(P)),
                 Sa * exp(Factor$muLnF + Factor$sigmaLnF * stats::qnorm(P)), tolerance = 1e-14)
    expect_equal(exp(Site$muLnSaF + Site$sigmaLnSaF^2 / 2), Sa * exp(Factor$muLnF + Factor$sigmaLnF^2 / 2),
                 tolerance = 1e-14)
  }
  expect_error(integrateSiteSpectrum(Reference$uhs, VS30 = 270, Vref = 760, model = "ST20", order = 3),
               "integer >= 4", fixed = TRUE)
})

# The general simulator: since the correction of finding F1 of
# validation-site-response its rock marginals are the monotone tailed
# quantile of the quadrature route, not the cubic tails of newmark 1.7.1, so
# these blocks test its mathematical contract; the measured compatibility with
# 1.7.1 that still holds is the last block. 1e-12 is the tolerance
# test-fitSaF.R of newmark 1.7.1 gives the re-centring of the quadrature route.

# The copula pair the simulator draws first at period t, as documented:
# bivariate standard normal with correlation rho_BJ(t, 0) by mvtnorm::rmvnorm().
copulaSample <- function(t, ns) {
  Rho <- correlateResiduals(Tn = t, T0 = 0)
  mvtnorm::rmvnorm(ns, sigma = matrix(c(1, Rho, Rho, 1), nrow = 2L))
}

test_that("simulateSiteSpectrum maps every copula uniform through the monotone tailed quantile", {
  # The tail fixture of the quadrature test above; its cubic extrapolation
  # (buildQuantileSpline) decreases beyond Fold.
  Quantiles <- data.table::data.table(
    p = c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95),
    Sa = c(0.002172658, 0.002368442, 0.002523095, 0.003518529, 0.005471477, 0.006568370, 0.007305735)
  )
  Grid <- seq(0, 1, length.out = 10001L)
  Tailed <- .buildTailedQuantile(Quantiles)
  expect_gt(min(diff(Tailed(Grid))), 0)
  Fold <- min(Grid[-1L][diff(buildQuantileSpline(Quantiles)(Grid)) < 0])
  # A PGA-only spectrum whose mean row is the mean of that quantile function.
  Mean <- stats::integrate(function(u) exp(Tailed(u)), lower = 0, upper = 1)$value
  Spectrum <- data.table::data.table(Tn = 0, p = c(as.character(Quantiles$p), "mean"), Sa = c(Quantiles$Sa, Mean))
  set.seed(919)
  Draws <- simulateSiteSpectrum(Spectrum, VS30 = 270, Vref = 760, model = "ST20", NS = 300)
  set.seed(919)
  U <- stats::pnorm(copulaSample(0, 300)[, 1L])
  expect_true(any(U > Fold))
  # Every draw is Q(U) plus one re-centring constant, so the copula ranks hold.
  expect_equal(log(Draws$PGA) - Tailed(U), rep(log(Mean / mean(exp(Tailed(U)))), 300L), tolerance = 1e-12)
  expect_false(is.unsorted(Draws$PGA[order(U)]))
})

test_that("simulateSiteSpectrum samples the closed-form quantiles of a single lognormal within Monte Carlo error", {
  # The single-lognormal fixture of test-siteResponse-mixture.R (ln PGA ~
  # N(ln 0.4, 0.5), an atom at 0.1 g at Tn = 1) as a rock quantile spectrum
  # at the probability levels of the 1.7.1 product, with its closed-form mean.
  P <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
  Labels <- c("0.05", "0.10", "0.16", "0.50", "0.84", "0.90", "0.95")
  Mu <- log(0.4)
  Sigma <- 0.5
  Mean <- exp(Mu + Sigma^2 / 2)
  Spectrum <- data.table::data.table(Tn = rep(c(0, 1), each = 8L), p = rep(c(Labels, "mean"), 2L),
                                     Sa = c(stats::qlnorm(P, meanlog = Mu, sdlog = Sigma), Mean, rep(0.1, 8L)))
  NS <- 300L
  set.seed(919)
  Draws <- simulateSiteSpectrum(Spectrum, VS30 = 270, Vref = 760, model = "ST20", NS = NS)
  Sampled <- Draws[Tn == 0, summarizeProfiles(.SD, levels = Labels, properties = "Sa")]$Sa
  # The sample is re-centred from the mean of the tailed quantile function to
  # the lognormal mean, so its p-quantile tends to qlnorm(p) * exp(Shift).
  Tailed <- .buildTailedQuantile(Spectrum[Tn == 0 & p != "mean"])
  Shift <- log(Mean / stats::integrate(function(u) exp(Tailed(u)), lower = 0, upper = 1)$value)
  # Standard error of the log of the re-centred sample p-quantile of NS draws
  # (Bahadur representation of the sample quantile and delta method for the
  # division by the sample mean, under the lognormal law):
  # NS Var = p (1 - p) / g^2 + (exp(Sigma^2) - 1) + 2 (Phi(z - Sigma) - p) / g,
  # with z = qnorm(p) and g = dnorm(z) / Sigma.
  z <- stats::qnorm(P)
  g <- stats::dnorm(z) / Sigma
  SE <- sqrt((P * (1 - P) / g^2 + expm1(Sigma^2) + 2 * (stats::pnorm(z - Sigma) - P) / g) / NS)
  # Two-sided normal bound at ALPHA over the seven levels (Bonferroni). ALPHA
  # is a nominal setting under the lognormal reference law behind SE: the
  # sampled law is the bounded tailed spline re-centred, so ALPHA is not its
  # exact family-wise false-alarm probability. ALPHA = 1e-3 is the owner's
  # value (2026-10-13).
  ALPHA <- 1e-3
  Deviation <- log(Sampled) - (stats::qnorm(P, mean = Mu, sd = Sigma) + Shift)
  expect_lte(max(abs(Deviation) / SE), stats::qnorm(1 - ALPHA / (2 * length(P))))
  expect_identical(Draws[Tn == 0, PGA], Draws[Tn == 0, Sa])
  # The degenerate period is its atom in every realization.
  expect_equal(Draws[Tn == 1, Sa], rep(0.1, NS), tolerance = 1e-12)
})

test_that("simulateSiteSpectrum keeps the pairing, the random stream and the supplied means", {
  Reference <- siteReference()
  set.seed(919)
  Draws <- simulateSiteSpectrum(Reference$uhs, VS30 = 270, Vref = 760, model = "ST20", NS = 300)
  Seed <- .Random.seed
  expect_identical(nrow(Draws), 900L)
  expect_identical(names(Draws), c("Tn", "realization", "PGA", "Sa", "SaF"))
  expect_identical(Draws$realization, rep(1:300, 3L))
  # Per period in increasing order: the copula pair of every realization,
  # then its factor residual. Row j carries the j-th pair: PGA follows the
  # ranks of its first coordinate and Sa those of the second.
  set.seed(919)
  for (t in sort(unique(Reference$uhs$Tn))) {
    Z <- copulaSample(t, 300)
    stats::rnorm(300)
    expect_identical(order(Draws[Tn == t, PGA]), order(Z[, 1L]))
    expect_identical(order(Draws[Tn == t, Sa]), order(Z[, 2L]))
  }
  expect_identical(.Random.seed, Seed)
  expect_identical(Draws[Tn == 0, Sa], Draws[Tn == 0, PGA])
  # The samples are re-centred on the tabulated means of PGA and of the ordinate.
  Means <- Reference$uhs[p == "mean", .(Tn, Sa)]
  expect_equal(Draws[, .(Sa = mean(Sa)), by = Tn]$Sa, Means$Sa, tolerance = 1e-12)
  expect_equal(Draws[, .(PGA = mean(PGA)), by = Tn]$PGA, rep(Means[Tn == 0, Sa], 3L), tolerance = 1e-12)
})

test_that("simulateSiteSpectrum returns the rock realizations at the reference condition for any Vref", {
  Reference <- siteReference()
  set.seed(919)
  Unit <- simulateSiteSpectrum(Reference$uhs, VS30 = 270, Vref = 270, model = "ST20", NS = 300)
  Seed <- .Random.seed
  expect_identical(Unit$SaF, Unit$Sa)
  # Only the copula pairs are drawn: the unit factor has no residual.
  set.seed(919)
  for (t in sort(unique(Reference$uhs$Tn))) copulaSample(t, 300)
  expect_identical(.Random.seed, Seed)
  set.seed(919)
  expect_identical(simulateSiteSpectrum(Reference$uhs, VS30 = 760, Vref = 760, model = "ST20", NS = 300), Unit)
})

test_that("compatibility with the fitSaF Monte Carlo of newmark 1.7.1: the random stream", {
  # Measured against fixtures/siteResponse-reference.rds (1.7.1, same seed,
  # spectrum and arguments): the draws consume the stream exactly as 1.7.1.
  # Its summaries are no longer reproduced, its cubic tails being replaced.
  Reference <- siteReference()
  set.seed(919)
  simulateSiteSpectrum(Reference$uhs, VS30 = 270, Vref = 760, model = "ST20", NS = 300)
  # R-devel adds 100000 times its binomial kind to the header; these draws use
  # no binomial, so the other header kinds and the generator state stay exact.
  expect_identical(c(.Random.seed[1L] %% 100000L, .Random.seed[-1L]),
                   c(Reference$seedAfter[1L] %% 100000L, Reference$seedAfter[-1L]))
})
