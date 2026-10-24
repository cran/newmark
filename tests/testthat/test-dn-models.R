# Oracle: fixtures/dn-reference.rds, written by newmark 1.7.1 in its own
# process; the model inputs are the roxygen examples of its Dn_models.R
# (SHA-256 b9533213ed5c3594...) and getKyLimits.R. The boundary checks use
# the thresholds of the published equations.
reference <- function() readRDS(test_path("fixtures", "dn-reference.rds"))

parameters <- function(x) c(muLnD = x$muLnD, sdLnD = x$sigmaLnD)

test_that("the six displacement models reproduce the 1.7.1 examples", {
  Rigid <- function(model) parameters(evaluateRigidDisplacement(model, PGA = 0.5, ky = 0.1, AI = estimateAI(0.5)))
  Flexible <- function(model, Sa, Mw, pulse, NFC = "D100") {
    parameters(evaluateFlexibleDisplacement(model, ky = 0.1, Ts = 0.3, Sa = Sa, Mw = Mw, PGV = 80, NFC = NFC, pulse = pulse))
  }
  expect_identical(list(Rigid("AM88"), Rigid("JB07"), Rigid("SR08"), Flexible("BT07", 0.8, 7, FALSE),
                        Flexible("BM17", 0.8, 8, FALSE), Flexible("BM19", 0.8, 7, FALSE),
                        Flexible("BM19", 0.8, 7, TRUE), Flexible("BM19", 0.8, 7, TRUE, NFC = "D50")),
                   reference()$models)
  expect_identical(evaluateRigidDisplacement("JB07", PGA = 0.5, ky = c(0.1, 0.2), AI = estimateAI(0.5))$IDn,
                   c("JB07", "JB07"))
})

test_that("the models switch their coefficients at the published thresholds", {
  # AM88 caps ky / PGA at 0.9999, so a block yielding at or above the peak
  # keeps the displacement of that ratio.
  expect_identical(evaluateRigidDisplacement("AM88", PGA = 0.5, ky = c(0.5, 0.7), AI = c(1, 1))$muLnD,
                   rep(evaluateRigidDisplacement("AM88", PGA = 0.5, ky = 0.49995, AI = 1)$muLnD, 2L))
  # BT07: intercept -0.22 below Ts = 0.05 s and -1.10 from there on.
  BT07 <- function(Ts) evaluateFlexibleDisplacement("BT07", ky = 0.1, Ts = Ts, Sa = 0.8, Mw = 7, PGV = 80, NFC = "D100",
                                                    pulse = FALSE)$muLnD
  expect_equal(BT07(0.05) - BT07(0.049), -1.10 + 0.22 + 1.50 * 0.001, tolerance = 1e-12)
  # BM19 pulse: c1 and c4 change above PGV = 150 cm/s.
  BM19 <- function(PGV) evaluateFlexibleDisplacement("BM19", ky = 0.1, Ts = 0.3, Sa = 0.8, Mw = 7, PGV = PGV,
                                                     NFC = "D100", pulse = TRUE)$muLnD
  expect_equal(BM19(150 + 1e-9) - BM19(150), (1.764 - -6.951) + (-0.097 - 1.547) * log(150), tolerance = 1e-6)
  expect_identical(evaluateFlexibleDisplacement("BM19", ky = 0.1, Ts = 0.3, Sa = c(0.8, 0.8), Mw = 7, PGV = c(80, 400),
                                                NFC = "D50", pulse = FALSE)$sigmaLnD, c(0.72, 0.72))
  expect_error(evaluateFlexibleDisplacement("BM19", ky = 0.1, Ts = 0.3, Sa = 0.8, Mw = 7, PGV = 80, NFC = "D84",
                                            pulse = FALSE), "NFC")
  expect_error(evaluateFlexibleDisplacement("BT07", ky = 0.1, Ts = c(0.3, 0.4), Sa = 0.8, Mw = 7, PGV = 80, NFC = "D100",
                                            pulse = FALSE), "one finite")
  expect_error(evaluateRigidDisplacement("BT07", PGA = 0.5, ky = 0.1, AI = 1), "model must be")
})

test_that("the models expand vectors of length one and reject other unequal lengths", {
  PGA <- rep(0.5, 4L)
  expect_identical(nrow(evaluateRigidDisplacement("SR08", PGA = 0.5, ky = c(0.1, 0.2), AI = estimateAI(0.5))), 2L)
  expect_identical(nrow(evaluateFlexibleDisplacement("BT07", ky = 0.1, Ts = 0.3, Sa = rep(0.8, 4L), Mw = 7,
                                                     PGV = rep(80, 4L), NFC = "D100", pulse = FALSE)), 4L)
  expect_error(evaluateRigidDisplacement("JB07", PGA = PGA, ky = c(0.1, 0.2), AI = estimateAI(PGA)),
               "PGA = 4, ky = 2, AI = 4", fixed = TRUE)
  expect_error(evaluateRigidDisplacement("AM88", PGA = PGA, ky = c(0.1, 0.2), AI = estimateAI(PGA)),
               "PGA = 4, ky = 2", fixed = TRUE)
  expect_error(evaluateFlexibleDisplacement("BT07", ky = c(0.1, 0.2), Ts = 0.3, Sa = rep(0.8, 4L), Mw = 7,
                                            PGV = rep(80, 4L), NFC = "D100", pulse = FALSE),
               "ky = 2, Sa = 4", fixed = TRUE)
})

test_that("estimatePGV and estimateAI evaluate the in-house regressions", {
  expect_equal(estimatePGV(1), 100 * exp(0.1241), tolerance = 1e-15)
  expect_equal(estimateAI(1), exp(2.6109), tolerance = 1e-15)
  expect_error(estimateAI(-0.1), "positive")
})

test_that("selectDisplacementModels keeps the flexible model of the regime", {
  Switches <- c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM17 = 1, BM19 = 1)
  expect_message(Crustal <- selectDisplacementModels(Switches, subduction = FALSE), "BM17 off, BM19 on")
  expect_identical(Crustal, c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM19 = 1))
  expect_identical(suppressMessages(selectDisplacementModels(Switches, subduction = TRUE)),
                   c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM17 = 1))
  expect_identical(selectDisplacementModels(c(AM88 = 1, BM17 = 0, BM19 = 0), subduction = TRUE), c(AM88 = 1))
  expect_error(selectDisplacementModels(c(AM88 = 0.5), subduction = TRUE), "0 or 1")
  expect_error(selectDisplacementModels(c(AM88 = 1, BM17 = 1), subduction = TRUE), "without BM19")
  expect_error(selectDisplacementModels(c(AM88 = 1, BM17 = 1, BM19 = 0), subduction = TRUE), "alike")
  expect_error(selectDisplacementModels(c(XX88 = 1), subduction = TRUE), "unsupported")
})

test_that("buildKyLimits tabulates the calibration ranges of 1.7.1", {
  expect_identical(as.list(buildKyLimits(PGA = 0.5)), stats::setNames(reference()$limits, c("IDn", "kyMin", "kyMax")))
})
