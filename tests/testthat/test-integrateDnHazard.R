# The PGA of the power-law hazard AEP = 0.001 PGA^-2 of these tests at the
# annual rate `rate`.
.invertHazard <- function(rate) sqrt(.001 / rate)

# Curves on unequal supports: the power-law PGA curve over the rates
# [1e-5, 1] (80 points) and twice its ordinate at Tn = 1 over [1e-4, 0.5]
# (20 points of another grid), each grid with exact end points.
UnequalCurves <- rbind(transform(data.frame(Tn = 0, AEP = 10^seq(-5, 0, length.out = 80L)), Sa = .invertHazard(AEP)),
                       transform(data.frame(Tn = 1, AEP = 10^seq(-4, log10(.5), length.out = 20L)), Sa = 2 * .invertHazard(AEP)))

# The log spectral accelerations u where the BT07 or BM17 median equals
# log(Da). Each median is the parabola m - a (u - v)^2 in u, with a and the
# vertex v of evaluateFlexibleDisplacement.R:131-132 (BT07: a = 0.244,
# v = (3.04 + 0.566 log ky) / 0.488) and :138-139 (BM17: a = 0.225,
# v = (3.060 + 0.538 log ky) / 0.450) and m its value at v, so the roots are
# v -+ sqrt((m - log Da) / a), one at a tangency and none below it.
.solveFlexibleMedian <- function(model, args) {
  a <- c(BT07 = .244, BM17 = .225)[[model]]
  Vertex <- c(BT07 = (3.04 + .566 * log(args$ky)) / .488, BM17 = (3.060 + .538 * log(args$ky)) / .450)[[model]]
  Height <- evaluateFlexibleDisplacement(model = model, ky = args$ky, Ts = args$Ts, Sa = exp(Vertex), Mw = args$Mw,
                                         PGV = 1, NFC = "D100", pulse = FALSE)$muLnD - log(args$Da)
  if (Height < 0) return(numeric())
  unique(Vertex + c(-1, 1) * sqrt(Height / a))
}

# One model's probability that D > Da at each rate of `rates` by a midpoint
# rule independent of integrateDnHazard(): `n` points per piece of the
# site-factor residual e in [-8, 8] weighted by dnorm. A rigid model splits
# [-8, 8] at the AM88 cap ky / PGA = 0.9999 and at the one intensity where
# its median displacement, monotone in log PGA, equals Da, bracketed by a
# scan of 2001 values over the intensities mu +- 8 sigma of these rates and
# refined by uniroot; a flexible model splits at the analytic roots of
# .solveFlexibleMedian(). `Sa` and `PGA` give the rock ordinates at the model
# period as functions of the rate; `args` holds the arguments of
# integrateDnHazard().
.integrateDenseProbability <- function(model, args, rates, Sa, PGA, n) {
  IsRigid <- model %in% c("AM88", "JB07", "SR08")
  .evaluateModel <- function(IM) {
    if (IsRigid) return(evaluateRigidDisplacement(model = model, PGA = IM, ky = args$ky, AI = estimateAI(IM)))
    evaluateFlexibleDisplacement(model = model, ky = args$ky, Ts = args$Ts, Sa = IM, Mw = args$Mw,
                                 PGV = rep(1, length(IM)), NFC = "D100", pulse = FALSE)
  }
  Site <- evaluateSiteSa(Sa = Sa(rates), PGA = PGA(rates), Tn = if (IsRigid) 0 else 1.5 * args$Ts,
                         VS30 = args$VS30, Vref = args$Vref, model = args$siteModel)
  if (IsRigid) {
    .measureMedian <- function(u) .evaluateModel(exp(u))$muLnD - log(args$Da)
    u <- seq(min(Site$muLnSaF - 8 * Site$sigmaLnSaF), max(Site$muLnSaF + 8 * Site$sigmaLnSaF), length.out = 2001L)
    AUX <- sign(.measureMedian(u))
    Splits <- c(log(args$ky / .9999), vapply(which(AUX[-1L] != AUX[-length(AUX)]),
                                             function(i) uniroot(.measureMedian, u[c(i, i + 1L)])$root, numeric(1L)))
  } else {
    Splits <- .solveFlexibleMedian(model, args)
  }
  Splits <- sort(Splits)
  mapply(function(mu, sigma) {
    Bounds <- c(-8, if (sigma > 0) pmin(pmax((Splits - mu) / sigma, -8), 8), 8)
    h <- diff(Bounds) / n
    e <- Bounds[-length(Bounds)] + outer(h, seq_len(n) - .5)
    Model <- .evaluateModel(exp(mu + sigma * as.vector(e)))
    sum(pnorm((log(args$Da) - Model$muLnD) / Model$sigmaLnD, lower.tail = FALSE) * dnorm(e) * h)
  }, Site$muLnSaF, Site$sigmaLnSaF)
}

# One model's annual rate of D > Da over the rates `limits`: the trapezoid
# rule on 4001 values of z = log(rate) of the probabilities of
# .integrateDenseProbability() with 1601 points per piece. Its resolution is
# established per dimension: the rate grid doubled alone and the residual
# grid doubled alone each agree with it within 1e-5 relative.
.integrateDenseAEP <- function(model, args, limits, Sa, PGA) {
  .integrateRate <- function(nz, ne) {
    z <- seq(log(limits[1L]), log(limits[2L]), length.out = nz)
    y <- .integrateDenseProbability(model, args, rates = exp(z), Sa = Sa, PGA = PGA, n = ne) * exp(z)
    sum(diff(z) * (y[-1L] + y[-length(y)]) / 2)
  }
  OUT <- .integrateRate(nz = 4001L, ne = 1601L)
  expect_lt(abs(.integrateRate(nz = 8001L, ne = 1601L) - OUT) / OUT, 1e-5)
  expect_lt(abs(.integrateRate(nz = 4001L, ne = 3202L) - OUT) / OUT, 1e-5)
  OUT
}

test_that("power-law hazard with JB07 agrees with its closed forms, with and without site dispersion", {
  skip_on_cran()
  PGA <- exp(seq(log(1e-5), log(1e4), length.out = 80L))
  Curves <- data.frame(Tn = 0, Sa = PGA, AEP = .001 * PGA^-2)
  expect_no_message(OUT <- integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 760, Vref = 760,
    Mw = 8, weights = c(JB07 = 1), Da = 100, IT = 5,
    TR = c(200, 100000), siteModel = "ST20"), class = "moduleNote")
  a <- .561 * 1.9228 + 3.833
  b <- .561 * 2.6109 - 3.833 * log(.2) - 1.474 * log(10)
  s <- .616 * log(10)
  Rate <- .001 * exp(-2 * (log(100) - b) / a + .5 * (2 * s / a)^2)
  expect_lt(abs(OUT$summary$AEP - Rate) / Rate, 1e-6)
  POE <- -expm1(-5 * OUT$summary$AEP)
  expect_lt(abs(OUT$summary$POE - POE) / POE, sqrt(.Machine$double.eps))
  expect_lte(OUT$summary$bandAEP, OUT$summary$AEP)
  expect_lt(abs(OUT$summary$rareTailAEPBound - min(Curves$AEP)) / min(Curves$AEP), sqrt(.Machine$double.eps))
  # With u = log(PGA), P(D > Da | u) = pnorm(Alpha u + Beta) and the rate
  # over [u1, u2] is antiderivative(u2) - antiderivative(u1).
  Alpha <- a / s
  Beta <- (b - log(100)) / s
  .evaluateAntiderivative <- function(u) -.001 * exp(-2 * u) * pnorm(Alpha * u + Beta) +
    .001 * exp(2 * Beta / Alpha + 2 / Alpha^2) * pnorm(Alpha * u + Beta + 2 / Alpha)
  .evaluateLogPGA <- function(rate) -log(rate / .001) / 2
  Support <- range(Curves$AEP)
  Band <- c(max(Support[1L], 1 / 100000), min(Support[2L], 1 / 200))
  ClosedForm <- .evaluateAntiderivative(.evaluateLogPGA(Support[1L])) - .evaluateAntiderivative(.evaluateLogPGA(Support[2L]))
  ClosedBand <- .evaluateAntiderivative(.evaluateLogPGA(Band[1L])) - .evaluateAntiderivative(.evaluateLogPGA(Band[2L]))
  expect_lt(abs(OUT$models$bandAEP - ClosedBand) / ClosedBand, 1e-5)
  expect_lt(abs(OUT$models$AEP - ClosedForm) / ClosedForm, 1e-5)
  expect_lt(max(abs(OUT$conditional$probability - pnorm(Alpha * .evaluateLogPGA(OUT$conditional$AEP) + Beta))), 1e-10)
  # The dense quadrature of the SR08, BM17 and ST17 test reproduces this integral.
  Dense <- .integrateDenseAEP("JB07", list(ky = .2, VS30 = 760, Vref = 760, Da = 100, siteModel = "ST20"),
                              limits = Support, Sa = .invertHazard, PGA = .invertHazard)
  expect_lt(abs(Dense - ClosedForm) / ClosedForm, 1e-4)
  # With site dispersion, ln D = a ln PGA + b of the amplified PGA and its
  # Arias intensity of estimateAI(), so the exceedance averaged over the
  # site residual is a normal probability.
  OUT <- integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 217, Vref = 760, Mw = 8, weights = c(JB07 = 1),
                           Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  Site <- evaluateSiteSa(Sa = .invertHazard(OUT$conditional$AEP), PGA = .invertHazard(OUT$conditional$AEP),
                         Tn = 0, VS30 = 217, Vref = 760, model = "ST20")
  expect_lt(max(abs(OUT$conditional$probability -
                      pnorm((a * Site$muLnSaF + b - log(100)) / sqrt(s^2 + a^2 * Site$sigmaLnSaF^2)))), 1e-8)
  # The band between the two adjacent knots where that probability crosses
  # 0.5 agrees with the dense quadrature.
  Site <- evaluateSiteSa(Sa = Curves$Sa, PGA = Curves$Sa, Tn = 0, VS30 = 217, Vref = 760, model = "ST20")
  Above <- pnorm((a * Site$muLnSaF + b - log(100)) / sqrt(s^2 + a^2 * Site$sigmaLnSaF^2)) >= .5
  TR <- 1 / c(min(Curves$AEP[!Above]), max(Curves$AEP[Above]))
  OUT <- integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 217, Vref = 760, Mw = 8, weights = c(JB07 = 1),
                           Da = 100, IT = 5, TR = TR, siteModel = "ST20")
  Dense <- .integrateDenseAEP("JB07", list(ky = .2, VS30 = 217, Vref = 760, Da = 100, siteModel = "ST20"),
                              limits = 1 / rev(TR), Sa = .invertHazard, PGA = .invertHazard)
  expect_lt(abs(OUT$models$bandAEP - Dense) / Dense, 1e-4)
})

test_that("model weighting acts on annual rates and leaves input unchanged", {
  PGA <- exp(seq(log(.001), log(3), length.out = 30L))
  Curves <- rbind(data.frame(Tn = 0, Sa = PGA, AEP = .001 * PGA^-2),
                  data.frame(Tn = 1, Sa = 2 * PGA, AEP = .001 * PGA^-2),
                  data.frame(Tn = 2, Sa = PGA, AEP = .001 * PGA^-2))
  Before <- Curves
  expect_no_message(OUT <- integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 217, Vref = 760,
    Mw = 8, weights = c(AM88 = 1, BT07 = 3), Da = 100, IT = 5,
    TR = c(200, 100000), siteModel = "ST20"), class = "moduleNote")
  Rate <- sum(OUT$models$AEP * c(.25, .75))
  expect_lt(abs(OUT$summary$AEP - Rate) / Rate, sqrt(.Machine$double.eps))
  expect_identical(Curves, Before)
  expect_true(all(OUT$conditional$probability >= 0 & OUT$conditional$probability <= 1))
  AUX <- merge(OUT$conditional, OUT$models[, .(IDn, TR.min, TR.max)], by = "IDn")
  expect_true(all(AUX$TR >= AUX$TR.min & AUX$TR <= AUX$TR.max))
  expect_error(integrateDnHazard(curves = Curves, ky = .2, Ts = 2, VS30 = 217, Vref = 760, Mw = 8,
                                 weights = c(BT07 = 1), Da = 100, IT = 5, TR = c(200, 100000),
                                 siteModel = "ST20"), "does not cover")
  expect_identical(
    integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 217, Vref = 760, Mw = 8,
      weights = c(AM88 = 1e308, BT07 = 1e308), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")$summary,
    integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 217, Vref = 760, Mw = 8,
      weights = c(AM88 = 1, BT07 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")$summary)
})

test_that("quadrature resolves the knots of a tabulated hazard curve", {
  Curves <- data.frame(Tn = 0,
    Sa = c(1e-05, 1.44e-05, 2.07e-05, 2.98e-05, 4.28e-05, 6.16e-05,
      8.87e-05, .0001276, .0001836, .0002641, .00038, .0005467,
      .0007866, .0011317, .0016282, .0023425, .0033702, .0048488,
      .0069761, .0100367, .0144401, .0207754, .0298902, .0430038,
      .0618708, .0890152, .1280686, .1842559, .265094, .381398,
      .5487279, .7894699, 1.1358321, 1.634153, 2.351101, 3.3825938),
    AEP = c(2.6294432929264, 2.59581452362498, 2.55167870340321,
      2.49378070803038, 2.41724762748248, 2.31661809662547,
      2.18813204002533, 2.03245719749849, 1.85537937242431,
      1.66572259044533, 1.47199866542859, 1.28004891712667,
      1.093085191251, .91363575943465, .744957871932129,
      .591258685395498, .456475705616138, .342762317768807,
      .25019105949596, .177224289959675, .121742978390443,
      .081216031582229, .052816443263762, .0335888522953669,
      .0208672341890493, .0125715518532524, .00723451410423333,
      .00388939489980349, .00189692502518957, .000805580893173562,
      .000279060133657404, 7.01310891273437e-05, 1.27480912565701e-05,
      1.87756676262737e-06, 2.95959743796081e-07, 1.22298400747845e-08))
  expect_no_message(OUT <- integrateDnHazard(Curves, ky = .15178479478321799, Ts = .95,
    VS30 = 217, Vref = 760, Mw = 8, weights = c(AM88 = 1),
    Da = 170, IT = 5, TR = c(200, 100000), siteModel = "ST20"), class = "moduleNote")
  expect_true(is.finite(OUT$summary$AEP))
  expect_gt(OUT$summary$AEP, 0)
  expect_lte(OUT$models$integrationError, 1e-6 * OUT$summary$AEP)
})

test_that("a repeated AEP node keeps the larger Sa and is noted, and a non-numeric curve is refused", {
  PGA <- exp(seq(log(.001), log(3), length.out = 30L))
  Curves <- data.frame(Tn = 0, Sa = PGA, AEP = .001 * PGA^-2)
  Curves$AEP[2L] <- Curves$AEP[1L]
  expect_message(OUT <- integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 217, Vref = 760,
    Mw = 8, weights = c(AM88 = 1), Da = 100, IT = 5,
    TR = c(200, 100000), siteModel = "ST20"), class = "moduleNote")
  expect_lt(abs(OUT$conditional[AEP == Curves$AEP[1L], Sa] - PGA[2L]) / PGA[2L], 1e-12)
  Curves$Sa <- as.character(Curves$Sa)
  expect_error(integrateDnHazard(Curves, ky = .2, Ts = 1, VS30 = 217, Vref = 760,
    Mw = 8, weights = c(AM88 = 1), Da = 100, IT = 5,
    TR = c(200, 100000), siteModel = "ST20"), "Tn, Sa and AEP must be numeric")
})

test_that("an intermediate period interpolates log Sa between its neighbours at equal rate", {
  PGA <- exp(seq(log(.001), log(3), length.out = 30L))
  Curves <- rbind(data.frame(Tn = 0, Sa = PGA, AEP = .001 * PGA^-2),
                  data.frame(Tn = .5, Sa = 2 * PGA, AEP = .001 * PGA^-2),
                  data.frame(Tn = 2, Sa = PGA, AEP = .001 * PGA^-2))
  OUT <- integrateDnHazard(Curves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(BT07 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  f <- (log(1 + 1e-8) - log(.5 + 1e-8)) / (log(2 + 1e-8) - log(.5 + 1e-8))
  PGA <- .invertHazard(OUT$conditional$AEP)
  Sa <- exp((1 - f) * log(2 * PGA) + f * log(PGA))
  expect_lt(max(abs(OUT$conditional$Sa - Sa) / Sa), 1e-12)
})

test_that("an intermediate period between three unequal grids interpolates at equal rate on the union of knots", {
  PGA <- exp(seq(log(.001), log(3), length.out = 30L))
  Curves <- rbind(UnequalCurves[UnequalCurves$Tn == 0, ], transform(UnequalCurves[UnequalCurves$Tn == 1, ], Tn = .5),
                  data.frame(Tn = 2, Sa = PGA, AEP = .001 * PGA^-2))
  OUT <- integrateDnHazard(Curves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(BT07 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  Support <- c(max(tapply(Curves$AEP, Curves$Tn, min)), min(tapply(Curves$AEP, Curves$Tn, max)))
  AUX <- sort(unique(c(Curves$AEP, 1 / c(200, 100000))))
  expect_identical(OUT$conditional$AEP, AUX[AUX >= Support[1L] & AUX <= Support[2L]])
  f <- (log(1 + 1e-8) - log(.5 + 1e-8)) / (log(2 + 1e-8) - log(.5 + 1e-8))
  PGA <- .invertHazard(OUT$conditional$AEP)
  Sa <- exp((1 - f) * log(2 * PGA) + f * log(PGA))
  expect_lt(max(abs(OUT$conditional$Sa - Sa) / Sa), 1e-12)
})

test_that("each model integrates over the common support of its curves and their knots", {
  OUT <- integrateDnHazard(UnequalCurves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(AM88 = 1, BT07 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  Flexible <- OUT$models[IDn == "BT07"]
  expect_identical(c(Flexible$TR.min, Flexible$TR.max, Flexible$rareTailBound), c(2, 1e4, 1e-4))
  Rigid <- OUT$models[IDn == "AM88"]
  expect_identical(c(Rigid$TR.min, Rigid$rareTailBound), c(1, 1e-5))
  # 1 / 1e-5 is not 1e5 in double precision.
  expect_identical(Rigid$TR.max, 1 / 1e-5)
  AUX <- sort(unique(c(UnequalCurves$AEP, 1e-4, .5, 1 / c(200, 100000))))
  expect_identical(OUT$conditional[IDn == "BT07", AEP], AUX[AUX >= 1e-4 & AUX <= .5])
})

test_that("the calibration range of ky comes from buildKyLimits at the largest PGA of the support", {
  Args <- list(curves = UnequalCurves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
               weights = c(AM88 = 1, BT07 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  OUT <- do.call(integrateDnHazard, Args)
  expect_identical(unlist(OUT$models[IDn == "BT07", .(kyMin, kyMax, extrapolated)], use.names = FALSE), c(.02, .40, 0))
  # The PGA support of AM88 ends at the rate 1e-5, where the PGA is 10 g.
  expect_lt(abs(OUT$models[IDn == "AM88", kyMax] - 10) / 10, sqrt(.Machine$double.eps))
  OUT <- do.call(integrateDnHazard, modifyList(Args, list(ky = .5)))
  expect_identical(OUT$models[IDn == "BT07", extrapolated], 1)
})

test_that("the band reaches the support end points and vanishes outside the support", {
  Full <- integrateDnHazard(UnequalCurves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(AM88 = 1, BT07 = 1), Da = 100, IT = 5, TR = c(2, 1e4), siteModel = "ST20")
  expect_lt(abs(Full$models[IDn == "BT07", bandAEP] - Full$models[IDn == "BT07", AEP]) / Full$models[IDn == "BT07", AEP], 1e-9)
  expect_true(all(c(.5, 1e-4) %in% Full$conditional[IDn == "AM88", AEP]))
  expect_true(all(c(.5, 1e-4) %in% Full$conditional[IDn == "BT07", AEP]))
  Empty <- integrateDnHazard(UnequalCurves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(AM88 = 1, BT07 = 1), Da = 100, IT = 5, TR = c(1e6, 1e7), siteModel = "ST20")
  expect_identical(c(Empty$summary$bandAEP, Empty$summary$bandPOE, Empty$models$bandAEP), c(0, 0, 0, 0))
  AUX <- c(Full$models$bandIntegrationError, Empty$models$bandIntegrationError)
  expect_length(AUX, 4L)
  expect_true(all(is.finite(AUX) & AUX >= 0))
})

test_that("SR08, BM17 and BT07 with ST17 agree with the dense quadrature", {
  skip_on_cran()
  Args <- list(curves = UnequalCurves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
               weights = c(SR08 = 1, BM17 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  OUT <- do.call(integrateDnHazard, Args)
  Dense <- .integrateDenseAEP("SR08", Args, limits = c(1e-5, 1), Sa = .invertHazard, PGA = .invertHazard)
  expect_lt(abs(OUT$models[IDn == "SR08", AEP] - Dense) / Dense, 1e-4)
  Dense <- .integrateDenseAEP("BM17", Args, limits = c(1e-4, .5), Sa = function(rate) 2 * .invertHazard(rate),
                              PGA = .invertHazard)
  expect_lt(abs(OUT$models[IDn == "BM17", AEP] - Dense) / Dense, 1e-4)
  Args <- modifyList(Args, list(VS30 = 300, weights = c(BT07 = 1), siteModel = "ST17"))
  OUT <- do.call(integrateDnHazard, Args)
  Dense <- .integrateDenseAEP("BT07", Args, limits = c(1e-4, .5), Sa = function(rate) 2 * .invertHazard(rate),
                              PGA = .invertHazard)
  expect_lt(abs(OUT$models$AEP - Dense) / Dense, 1e-4)
})

test_that("a keyed data.table and the random stream are left unchanged and invalid curves are refused", {
  Curves <- data.table::as.data.table(UnequalCurves)
  data.table::setkey(Curves, Tn, Sa)
  Before <- data.table::copy(Curves)
  set.seed(1)
  Seed <- .Random.seed
  integrateDnHazard(Curves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(AM88 = 1, BT07 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  expect_identical(Curves, Before)
  expect_identical(.Random.seed, Seed)
  Curves <- transform(UnequalCurves, Sa = TRUE)
  expect_error(integrateDnHazard(Curves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(AM88 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20"), "Tn, Sa and AEP must be numeric")
  Curves <- UnequalCurves
  Curves$AEP[3L] <- NA
  expect_error(integrateDnHazard(Curves, ky = .2, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
    weights = c(AM88 = 1), Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20"), "invalid rock hazard curves")
})

test_that("the former refusal case converges and agrees with the dense quadrature", {
  skip_on_cran()
  # A low threshold and a weak slope make the exceedance a steep transition in
  # the site-factor residual next to the AM88 cap; the pieces split at the cap
  # and at the median displacement resolve it.
  Args <- list(curves = UnequalCurves, ky = .05, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8,
               weights = c(AM88 = 1), Da = .1, IT = 5, TR = c(200, 100000), siteModel = "ST20")
  OUT <- do.call(integrateDnHazard, Args)
  Dense <- .integrateDenseAEP("AM88", Args, limits = c(1e-5, 1), Sa = .invertHazard, PGA = .invertHazard)
  expect_lt(abs(OUT$models$AEP - Dense) / Dense, 1e-4)
  # The acceptance criteria of the accepted order, checked as such; they are
  # not independent comparisons.
  expect_lte(OUT$models$quadratureError, 1e-10 + 1e-4 * OUT$models$AEP)
  expect_lte(OUT$models$bandQuadratureError, 1e-10 + 1e-4 * OUT$models$bandAEP)
  expect_lte(OUT$models$probabilityError, 1e-4)
})

test_that("the median crossings of BT07 and BM17 are found on both sides of the vertex as a close pair, a near-tangent pair and an end point, none when the maximum lies just below Da, one root at a maximum exactly zero, and the dense quadrature agrees across them", {
  skip_on_cran()
  PGA <- exp(seq(log(.001), log(3), length.out = 30L))
  Curves <- rbind(data.frame(Tn = 0, Sa = PGA, AEP = .001 * PGA^-2),
                  data.frame(Tn = .5, Sa = 2 * PGA, AEP = .001 * PGA^-2),
                  data.frame(Tn = 2, Sa = PGA, AEP = .001 * PGA^-2))
  f <- (log(1 + 1e-8) - log(.5 + 1e-8)) / (log(2 + 1e-8) - log(.5 + 1e-8))
  .interpolateSa <- function(rate) exp((1 - f) * log(2 * .invertHazard(rate)) + f * log(.invertHazard(rate)))
  # The BT07 median at Ts = 2/3 and Mw = 8 is concave in u = log Sa with its
  # vertex at u* = (3.04 + 0.566 log ky) / 0.488; with Da the median
  # displacement at the vertex times exp(-0.244) it equals log(Da) at
  # u* -+ sqrt((mu* - log Da) / 0.244) = u* -+ 1.
  .evaluateMedian <- function(model, ky, u) {
    evaluateFlexibleDisplacement(model = model, ky = ky, Ts = 2 / 3, Sa = exp(u), Mw = 8,
                                 PGV = rep(1, length(u)), NFC = "D100", pulse = FALSE)$muLnD
  }
  Vertex <- (3.04 + .566 * log(.02)) / .488
  Args <- list(curves = Curves, ky = .02, Ts = 2 / 3, VS30 = 217, Vref = 760, Mw = 8, weights = c(BT07 = 1),
               Da = exp(.evaluateMedian("BT07", .02, Vertex) - .244), IT = 5, TR = c(200, 100000), siteModel = "ST20")
  OUT <- do.call(integrateDnHazard, Args)
  # The enclosure of integrateDnHazard() for these nodes: both of its ends lie
  # below log(Da), and both crossings inside it.
  LnSa <- log(.interpolateSa(OUT$conditional$AEP))
  Site <- evaluateSiteSa(Sa = exp(LnSa), PGA = .invertHazard(OUT$conditional$AEP), Tn = 1, VS30 = 217, Vref = 760,
                         model = "ST20")
  MuF <- Site$muLnSaF - LnSa
  Enclosure <- c(min(LnSa) + min(MuF), max(LnSa) + max(MuF)) + c(-8, 8) * max(Site$sigmaLnSaF)
  expect_true(all(.evaluateMedian("BT07", .02, Enclosure) < log(Args$Da)))
  Roots <- .findMedianCrossings(function(u) .evaluateMedian("BT07", .02, u) - log(Args$Da), Enclosure)
  Closed <- Vertex + c(-1, 1) * sqrt((.evaluateMedian("BT07", .02, Vertex) - log(Args$Da)) / .244)
  expect_length(Roots, 2L)
  expect_lt(max(abs(Roots - Closed)), 1e-8)
  # With m the median at the vertex v and a the curvature of the parabola
  # (evaluateFlexibleDisplacement.R:131-132 and 138-139), Da = exp(m - a d^2 / 4)
  # puts the roots at v -+ d / 2: d = 0.001 keeps both inside one cell of a
  # 2001-value grid over this enclosure, so the median less log(Da) is negative
  # at every value of that grid. Da = exp(m - 5e-13) leaves the maximum 5e-13
  # above log(Da), with the roots v -+ sqrt(5e-13 / a) about 3e-6 apart, and
  # Da = exp(m + 5e-13) leaves it 5e-13 below, with no root. Whether the
  # regression's own tangent Da = exp(m) measures a maximum of exactly zero
  # depends on where stats::optimize() stops, so the zero maximum is a measure
  # exactly zero on [v - 1e-6, v + 1e-6] and negative outside: one root. An
  # enclosure that starts exactly at a root returns it once.
  for (model in c("BT07", "BM17")) {
    a <- c(BT07 = .244, BM17 = .225)[[model]]
    v <- c(BT07 = (3.04 + .566 * log(.02)) / .488, BM17 = (3.060 + .538 * log(.02)) / .450)[[model]]
    m <- .evaluateMedian(model, .02, v)
    Da <- exp(m - a * .001^2 / 4)
    Roots <- .findMedianCrossings(function(u) .evaluateMedian(model, .02, u) - log(Da), Enclosure)
    expect_length(Roots, 2L)
    expect_lt(max(abs(Roots - (v + c(-1, 1) * .0005))), 1e-8)
    Da <- exp(m - 5e-13)
    Roots <- .findMedianCrossings(function(u) .evaluateMedian(model, .02, u) - log(Da), Enclosure)
    expect_length(Roots, 2L)
    expect_lt(max(abs(Roots - .solveFlexibleMedian(model, modifyList(Args, list(Da = Da))))), 1e-8)
    Roots <- .findMedianCrossings(function(u) -pmax(abs(u - v) - 1e-6, 0)^2, Enclosure)
    expect_length(Roots, 1L)
    expect_lte(abs(Roots - v), 1e-6)
    Da <- exp(m + 5e-13)
    expect_length(.findMedianCrossings(function(u) .evaluateMedian(model, .02, u) - log(Da), Enclosure), 0L)
    Roots <- .findMedianCrossings(function(u) .evaluateMedian(model, .02, u) - .evaluateMedian(model, .02, v - 1),
                                  c(v - 1, Enclosure[2L]))
    expect_length(Roots, 2L)
    expect_identical(Roots[1L], v - 1)
    expect_lt(abs(Roots[2L] - (v + 1)), 1e-8)
  }
  Dense <- .integrateDenseAEP("BT07", Args, limits = range(Curves$AEP), Sa = .interpolateSa, PGA = .invertHazard)
  expect_lt(abs(OUT$models$AEP - Dense) / Dense, 1e-4)
  # The node probabilities, with the dense residual grid resolved by its
  # doubling alone (the rate grid does not enter them).
  Dense <- .integrateDenseProbability("BT07", Args, rates = OUT$conditional$AEP, Sa = .interpolateSa,
                                      PGA = .invertHazard, n = 1601L)
  expect_lt(max(abs(.integrateDenseProbability("BT07", Args, rates = OUT$conditional$AEP, Sa = .interpolateSa,
                                               PGA = .invertHazard, n = 3202L) - Dense)), 1e-6)
  expect_lt(max(abs(OUT$conditional$probability - Dense)), 1e-5)
  # At ky = 0.01 the vertex is lower and the site median of a node lies
  # beyond the second crossing, on the descending branch.
  Vertex <- (3.04 + .566 * log(.01)) / .488
  Args <- modifyList(Args, list(ky = .01, Da = exp(.evaluateMedian("BT07", .01, Vertex) - .244)))
  expect_gt(max(Site$muLnSaF), Vertex + 1)
  OUT <- do.call(integrateDnHazard, Args)
  Dense <- .integrateDenseAEP("BT07", Args, limits = range(Curves$AEP), Sa = .interpolateSa, PGA = .invertHazard)
  expect_lt(abs(OUT$models$AEP - Dense) / Dense, 1e-4)
})
