test_that("evaluateRoot returns tabulated roots exactly and marks support", {
  Prep <- buildRootModel(mode = 1, method = "nlm")
  Hit <- Prep$table[1L]
  OUT <- evaluateRoot(Prep, mo = Hit$m, lo = Hit$l)
  expect_identical(OUT$an, Hit$an)
  expect_true(OUT$exact)
  expect_true(OUT$withinSupport)
  Off <- evaluateRoot(Prep, mo = c(0.148, 1.5), lo = c(0.351, 0.351))
  expect_false(any(Off$exact))
  expect_identical(Off$withinSupport, c(TRUE, FALSE))
  # Queries below the tabulated minima, both 0, are extrapolations as well.
  expect_equal(c(Prep$loMin, Prep$moMin), c(0, 0))
  Below <- evaluateRoot(Prep, mo = c(-0.1, 0.148), lo = c(0.351, -0.05))
  expect_false(any(Below$exact))
  expect_identical(Below$withinSupport, c(FALSE, FALSE))
  # Fitted to the solved table, this surrogate gives 3.1171 for mo 0.148 and lo
  # 0.351, where the frozen ShearTable of newmark 1.7.1 records 3.1158.
  expect_equal(round(Off$an[1], 4), 3.1171)
  # One query mixing a grid node and an off-grid point keeps each row's origin.
  Mixed <- evaluateRoot(Prep, mo = c(0.148, Hit$m), lo = c(0.351, Hit$l))
  expect_identical(Mixed$exact, c(FALSE, TRUE))
  expect_identical(Mixed$an[2L], Hit$an)
  expect_equal(Mixed$an[1L], Off$an[1L])
  expect_error(buildRootModel(mode = 9, method = "nlm"), "tabulated mode")
  expect_error(buildRootModel(mode = 1, method = "rf"), "lm")
})

test_that("the exact model solves the characteristic relation of Dakoulas and Gazetas", {
  # Dakoulas and Gazetas (1985), Table 2 (m = 0) and Table 3 (m = 1/2) on the
  # printed page 170 (original PDF read on 2026-10-12; the hydrated copy is
  # faithful): modes 1 to 4 at lambda 0, 0.10, 0.25 and 0.50, printed to three
  # decimals, so a faithful root lies within the half unit 0.0005. Two cells of
  # Table 2 do not: lambda 0.50 mode 2 prints 9.605 and lambda 0.25 mode 4
  # prints 14.804, while the roots of the relation are 9.604122 and 14.803308
  # (relation residual below 1e-12 against 4e-5 at the printed values;
  # Ambraseys 1960, Table 1, prints 9.60 and 14.80). They are treated as
  # rounding errors of the source and asserted separately.
  Table2 <- rbind(c(2.405, 5.520, 8.654, 11.792), c(2.448, 5.726, 9.096, 12.510),
                  c(2.668, 6.580, 10.666, 14.804), c(3.588, 9.605, 15.818, 22.070))
  Table3 <- rbind(c(2.903, 6.033, 9.171, 12.310), c(2.974, 6.415, 10.035, 13.740),
                  c(3.301, 7.718, 12.422, 17.208), c(4.574, 11.903, 19.544, 27.245))
  Lambda <- c(0, 0.10, 0.25, 0.50)
  Solved2 <- sapply(1:4, function(Mode) evaluateRoot(buildRootModel(Mode, "exact"), mo = rep(0, 4), lo = Lambda)$an)
  Solved3 <- sapply(1:4, function(Mode) evaluateRoot(buildRootModel(Mode, "exact"), mo = rep(0.5, 4), lo = Lambda)$an)
  # Every solved root satisfies the relation to the solver's tolerance.
  Residual <- c(mapply(.characteristicRelation, Solved2, 0, rep(Lambda, 4)),
                mapply(.characteristicRelation, Solved3, 0.5, rep(Lambda, 4)))
  expect_lt(max(abs(Residual)), 1e-9)
  # Thirty cells agree within the printed half unit (absolute); the two
  # disputed cells differ by less than one unit and are reproduced to 1e-6.
  Disputed <- cbind(c(4L, 3L), c(2L, 4L))
  expect_lt(max(abs(Solved2 - Table2)[-c(8L, 15L)]), 5e-4)
  expect_lt(max(abs(Solved3 - Table3)), 5e-4)
  expect_lt(max(abs(Solved2[Disputed] - c(9.604122, 14.803308))), 1e-6)
  expect_true(all(abs(Solved2[Disputed] - Table2[Disputed]) > 5e-4))
  expect_true(all(abs(Solved2[Disputed] - Table2[Disputed]) < 1e-3))
  # Without truncation and inhomogeneity the first root is the first zero of J0.
  Exact <- evaluateRoot(buildRootModel(1, "exact"), mo = 0, lo = 0)
  expect_equal(Exact$an, 2.404825557695773, tolerance = 1e-9)
  expect_true(Exact$exact && Exact$withinSupport)
  expect_null(buildRootModel(1, "exact")$fit)
  expect_error(evaluateRoot(buildRootModel(1, "exact"), mo = 2, lo = 0.3), "mo in")
  # The table holds solved roots: at m = 0 and lambda = 0 its first root is the
  # first zero of J0, and its boundary rows m = 0 and lambda = 0 of every mode,
  # with every 500th row, satisfy the relation.
  Hit <- buildRootModel(1, "nlm")$table[l == 0 & m == 0]
  expect_lt(abs(Hit$an - 2.404825557695773), 1e-9)
  Sample <- .dataset("CylinderRoots")[m == 0 | l == 0 | seq_len(.N) %% 500L == 0L]
  expect_lt(max(abs(Sample[, mapply(.characteristicRelation, a = an, mo = m, lo = l)])), 1e-8)
})

test_that("the exact search finds the first root when the truncation ratio is near 1", {
  # The search limit is 100 / (1 - lo): with lo 0.99 the first root of mode 1
  # lies above 100 and below that limit.
  Root <- evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = 0.427, lo = 0.99)$an
  expect_true(is.finite(Root))
  expect_gt(Root, 100)
  expect_lt(abs(.characteristicRelation(Root, 0.427, 0.99)), 1e-8)
})

test_that("the general route reproduces the shear-beam periods of Table 8 of Dakoulas and Gazetas", {
  # Dakoulas and Gazetas (1985), hydrated copy, Table 8 (printed page 175):
  # dam height (crest to base, m), average S-wave velocity (m/s), truncation
  # ratio, inhomogeneity parameter and shear-beam T1 (s) of their five dams.
  # The period takes the apex-to-base height Hs / (1 - lo) (owner's ruling of
  # 2026-10-12) and must match within the 0.1 percent set for this check.
  Dams <- data.frame(Hs = c(40, 80, 120, 120, 120), Cbar = c(200, 245, 280, 280, 280),
                     lo = c(0.043, 0.027, 0.027, 0.027, 0.027), mo = c(0.57, 0.57, 0.57, 0.63, 0.63),
                     T1 = c(0.534, 0.860, 1.129, 1.129, 1.129))
  # Basal velocity from the average one through their expression for the
  # average S-wave velocity of the dam (eq-natural-average-velocity).
  VSo <- with(Dams, Cbar * (4 + mo) * (1 - lo^2) / (4 * (1 - lo^(2 + mo / 2))))
  Roots <- evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = Dams$mo, lo = Dams$lo)
  Ts <- evaluateBeamPeriod(Hs = Dams$Hs, lo = Dams$lo, mo = Dams$mo, VSo = VSo, an = Roots$an)
  expect_lt(max(abs(Ts / Dams$T1 - 1)), 0.001)
  expect_error(evaluateBeamPeriod(Hs = 10, lo = 0.3, mo = 2, VSo = 200, an = 3), "below 2")
  expect_error(evaluateBeamPeriod(Hs = 10, lo = 1, mo = 0.5, VSo = 200, an = 3), "lo in")
})
