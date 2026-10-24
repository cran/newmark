# Conformity with test-fitSaFMixture.R of newmark 1.7.1: same fixtures and
# argument values; the 1.7.1 outputs come from fixtures/siteResponse-reference.rds,
# written by the installed newmark 1.7.1 (build d57b64f) in its own R process,
# never together with newmark 5.
siteReference <- function() readRDS(test_path("fixtures", "siteResponse-reference.rds"))

mixtureMean <- function(x) x[, .(Sa = sum(w * exp(muLn + sigmaLn^2 / 2))), by = Tn]

test_that("integrateSiteMixture keeps the reference-condition identity and its input", {
  GMM <- siteReference()$mixture
  Before <- data.table::copy(GMM)
  OUT <- integrateSiteMixture(mixture = GMM, VS30 = 760, Vref = 760, model = "ST20", order = 32L)
  expect_equal(OUT, GMM, tolerance = 0)
  expect_identical(GMM, Before)
})

test_that("integrateSiteMixture reproduces fitSaFMixture of newmark 1.7.1 exactly", {
  Reference <- siteReference()
  OUT <- integrateSiteMixture(Reference$mixture, VS30 = 360, Vref = 760, model = "ST20", order = 16)
  expectReproduces(as.list(OUT), as.list(Reference$mixture16))
  expect_identical(names(OUT), c("Tn", "k", "muLn", "sigmaLn", "w"))
  expect_equal(OUT[, .N, by = Tn]$N, rep(32L, 3L))
  expect_equal(OUT[, sum(w), by = Tn]$V1, rep(1, 3L), tolerance = 1e-12)
  expect_true(OUT[, .(sigmaRms = sqrt(sum(w * sigmaLn^2))), by = Tn][, data.table::uniqueN(round(sigmaRms, 12L))] > 1L)
  expect_true(mixtureMean(OUT)[, all(is.finite(Sa) & Sa > 0)])
  High <- integrateSiteMixture(Reference$mixture, VS30 = 360, Vref = 760, model = "ST20", order = 64)
  expectReproduces(as.list(High), as.list(Reference$mixture64))
  expect_true(High[, all(is.finite(w) & w > 0)])
})

test_that("period-varying participation and degenerate rock spectra match newmark 1.7.1", {
  Reference <- siteReference()
  Varying <- data.table::copy(Reference$mixture)[!(Tn == 1 & k == "b")]
  Varying[Tn == 1, w := 1]
  OUT <- integrateSiteMixture(Varying, VS30 = 360, Vref = 760, model = "ST20", order = 8)
  expectReproduces(as.list(OUT), as.list(Reference$varying8))
  expect_equal(OUT[, .N, by = Tn]$N, c(16L, 16L, 8L))
  Degenerate <- data.table::copy(Reference$mixture)[, sigmaLn := 0]
  expectReproduces(as.list(integrateSiteMixture(Degenerate, VS30 = 360, Vref = 760, model = "ST20", order = 8)),
                   as.list(Reference$degenerate8))
})

test_that("inconsistent rock components are rejected", {
  GMM <- siteReference()$mixture
  Invalid <- data.table::copy(GMM)[Tn == 1 & k == "a", w := 0.3]
  expect_error(integrateSiteMixture(mixture = Invalid, VS30 = 360, Vref = 760, model = "ST20", order = 32L), "weights must sum to one", fixed = TRUE)
  Invalid <- data.table::copy(GMM)[Tn == 1 & k == "a", sigmaLn := -0.1]
  expect_error(integrateSiteMixture(mixture = Invalid, VS30 = 360, Vref = 760, model = "ST20", order = 32L), "invalid parameters", fixed = TRUE)
  expect_error(summarizeMixture(Invalid, levels = c(0.16, 0.5, 0.84)), "mixture contains invalid parameters", fixed = TRUE)
  Invalid <- data.table::rbindlist(list(GMM, data.table::data.table(Tn = 1, k = "c", muLn = -1, sigmaLn = 0.3, w = 0.2)))
  Invalid[Tn == 1, w := w / sum(w)]
  expect_error(integrateSiteMixture(mixture = Invalid, VS30 = 360, Vref = 760, model = "ST20", order = 32L), "no matching Tn == 0 PGA component", fixed = TRUE)
})

test_that("summarizeMixture is the closed form of a single lognormal", {
  Mixture <- data.table::data.table(Tn = c(0, 1), k = "a", muLn = c(log(0.4), log(0.1)), sigmaLn = c(0.5, 0), w = 1)
  Levels <- c(0.16, 0.5, 0.84)
  OUT <- summarizeMixture(Mixture, levels = Levels)
  expect_identical(OUT$level, rep(c("0.16", "0.5", "0.84", "mean"), 2L))
  # The solver runs uniroot with tol = 1e-12 on ln Sa, which stops when the
  # function is zero or its step is below tol plus an allowance for the
  # representation error of ln Sa; 4 eps |ln Sa| is that allowance in Brent's
  # zeroin stopping rule.
  LnClosed <- stats::qnorm(Levels, mean = log(0.4), sd = 0.5)
  expect_lte(max(abs(log(OUT[Tn == 0 & level != "mean", Sa]) - LnClosed)),
             1e-12 + 4 * .Machine$double.eps * max(abs(LnClosed)))
  expect_identical(OUT[Tn == 0 & level == "mean", Sa], exp(log(0.4) + 0.5^2 / 2))
  # The zero-dispersion component is an atom: every level inside its jump
  # returns it (generalized inverse).
  expect_identical(OUT[Tn == 1, Sa], rep(exp(log(0.1)), 4L))
  expect_false(is.unsorted(OUT[Tn == 0 & level != "mean", Sa], strictly = TRUE))
})

test_that("the median of twelve equal atoms is the seventh, as added in double", {
  # Six weights of 1/12 add to 0.49999999999999994 in double, below p = 0.5,
  # so the bisection stops at the seventh atom; added in long double they
  # reach 0.5 and would select the sixth.
  Mixture <- data.table::data.table(Tn = 0, k = sprintf("k%02d", 1:12), muLn = (1:12) / 10, sigmaLn = 0, w = 1 / 12)
  expect_identical(summarizeMixture(Mixture, levels = 0.5)[level == "0.5", Sa], exp(7 / 10))
})

test_that("refineSiteMixture returns the rock mixture at the reference condition", {
  GMM <- siteReference()$mixture
  OUT <- refineSiteMixture(GMM, VS30 = 760, Vref = 760, model = "ST20", levels = c(0.16, 0.84),
                           orders = c(8, 16, 32), tolerance = 2e-3)
  expect_equal(OUT$mixture, GMM, tolerance = 0)
  expect_identical(OUT$order, NA_integer_)
  expect_error(refineSiteMixture(GMM, VS30 = 360, Vref = 760, model = "ST20", levels = 0.5,
                                 orders = c(8, 16), tolerance = 2e-3), "at least three")
})
