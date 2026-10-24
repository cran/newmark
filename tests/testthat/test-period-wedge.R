# Numerical slope period of the truncated wedge (solveWedgePeriod). The
# expected values are closed forms and tables of Dakoulas and Gazetas (1985,
# hydrated copy) evaluated here, independently of the solver; "within the
# discretisation error" means that the distance to the closed form does not
# exceed the change between the profile sampled on N and on 2N layers.

# Layers of the power-law wedge G = Gb (z/H)^m of the paper, z from the apex,
# sampled at the midpoints of N equal layers below the crest (Hs = H - h).
powerLawLayers <- function(Hs, lo, m, Cb, N) {
  H <- Hs / (1 - lo)
  hs <- rep(Hs / N, N)
  zm <- cumsum(hs) - hs / 2
  data.frame(hs = hs, zm = zm, gsat = 19, VSm = Cb * ((lo * H + zm) / H)^(m / 2))
}

# Exact first period of the layered wedge: in layer k the general solution of
# eq-free-motion with constant modulus is U = A J0(w z / V_k) + B Y0(w z / V_k)
# (eq-free-solution with m = 0); U and G U' are continuous at the interfaces,
# G z U' = 0 at the crest (the regular branch, B = 0, at the apex) and
# U(H) = 0. `apex` is the crest's distance to the apex; `upper` bounds the
# search in rad/s.
wedgeExactPeriod <- function(apex, hs, VSm, gsat, upper) {
  n <- length(hs)
  z <- apex + c(0, cumsum(hs))
  G <- gsat * VSm^2
  baseDisplacement <- function(w) {
    k <- w / VSm
    AB <- if (apex == 0) c(1, 0) else c(besselY(k[1] * apex, 1), -besselJ(k[1] * apex, 1))
    for (i in seq_len(n - 1L)) {
      x <- k[i] * z[i + 1L]
      y <- k[i + 1L] * z[i + 1L]
      Continuity <- c(AB[1] * besselJ(x, 0) + AB[2] * besselY(x, 0),
                      -G[i] * k[i] * (AB[1] * besselJ(x, 1) + AB[2] * besselY(x, 1)))
      AB <- solve(matrix(c(besselJ(y, 0), -G[i + 1L] * k[i + 1L] * besselJ(y, 1),
                           besselY(y, 0), -G[i + 1L] * k[i + 1L] * besselY(y, 1)), 2L), Continuity)
    }
    x <- k[n] * z[n + 1L]
    AB[1] * besselJ(x, 0) + AB[2] * besselY(x, 0)
  }
  Grid <- seq(upper / 600, upper, length.out = 600)
  Value <- vapply(Grid, baseDisplacement, numeric(1L))
  k <- which(Value[-1L] * Value[-length(Value)] < 0)[1L]
  2 * pi / uniroot(baseDisplacement, Grid[c(k, k + 1L)], tol = upper * 1e-14)$root
}

test_that("layered power-law wedges reduce to the closed form of Dakoulas and Gazetas", {
  # eq-natural-frequency-period, T = 4 pi H / (a1 (2 - m) Cb), H = Hs / (1 - lo),
  # with a1 the first root of eq-free-characteristic; m and lambda as set by the
  # owner's approval of 2026-10-12. With mo = m and VSo = Cb the analytic route
  # of the solver is that closed form.
  Prep <- buildRootModel(mode = 1, method = "exact")
  for (m in c(0, 0.5, 0.95, 1.5)) {
    for (lo in c(0.075, 0.2, 0.351, 0.5)) {
      Closed <- 4 * pi * (10 / (1 - lo)) / (evaluateRoot(Prep, mo = m, lo = lo)$an * (2 - m) * 200)
      OUT <- lapply(c(50, 100), function(N) {
        solveWedgePeriod(powerLawLayers(10, lo, m, 200, N), Hs = 10, lo = lo, mo = m, VSo = 200,
                         method = "exact", elements = 2)
      })
      expect_equal(OUT[[2]]$TsAnalytic, Closed)
      expect_lte(abs(OUT[[2]]$Ts / Closed - 1), abs(OUT[[2]]$Ts / OUT[[1]]$Ts - 1))
    }
  }
})

test_that("the untruncated wedge reproduces Table 7 of Dakoulas and Gazetas", {
  # Table 7 (lambda = 0, printed page 172): T1 Cbar / H = 16 pi /
  # ((4 + m)(2 - m) a1), with Cbar = 4 Cb / (4 + m) (eq-natural-average-velocity
  # at lambda = 0), printed to three decimals; the printed 2.565 of m = 1/2 is
  # the exact 2.56555 truncated, so the printed values hold within one unit
  # of their last decimal.
  Table7 <- data.frame(m = c(0, 1 / 2, 2 / 3, 3 / 4), printed = c(2.613, 2.565, 2.571, 2.579))
  Prep <- buildRootModel(mode = 1, method = "exact")
  for (i in seq_len(nrow(Table7))) {
    m <- Table7$m[i]
    Closed <- 16 * pi / ((4 + m) * (2 - m) * evaluateRoot(Prep, mo = m, lo = 0)$an)
    Period <- vapply(c(50, 100), function(N) {
      solveWedgePeriod(powerLawLayers(10, 0, m, 200, N), Hs = 10, lo = 0, mo = m, VSo = 200,
                       method = "exact", elements = 2)$Ts * (4 * 200 / (4 + m)) / 10
    }, numeric(1L))
    expect_lte(abs(Period[2] / Closed - 1), abs(Period[2] / Period[1] - 1))
    expect_lte(abs(Closed - Table7$printed[i]), 1e-3)
    expect_lte(abs(Period[2] - Table7$printed[i]), 1e-3)
  }
})

test_that("the five dams of Table 8 of Dakoulas and Gazetas are reproduced within 0.2 percent", {
  # Table 8 (printed page 175): dam height (crest to base), average S-wave
  # velocity, truncation ratio, exponent and shear-beam T1. Cb from the
  # average velocity with eq-natural-average-velocity; the 0.2 percent is the
  # owner's approval of 2026-10-12.
  Dams <- data.frame(Hs = c(40, 80, 120, 120, 120), Cbar = c(200, 245, 280, 280, 280),
                     lo = c(0.043, 0.027, 0.027, 0.027, 0.027), m = c(0.57, 0.57, 0.57, 0.63, 0.63),
                     T1 = c(0.534, 0.860, 1.129, 1.129, 1.129))
  Dams$Cb <- with(Dams, Cbar * (4 + m) * (1 - lo^2) / (4 * (1 - lo^(2 + m / 2))))
  for (i in seq_len(nrow(Dams))) {
    Dam <- Dams[i, ]
    OUT <- solveWedgePeriod(powerLawLayers(Dam$Hs, Dam$lo, Dam$m, Dam$Cb, 200), Hs = Dam$Hs, lo = Dam$lo,
                            mo = Dam$m, VSo = Dam$Cb, method = "exact", elements = 1)
    expect_lte(abs(OUT$Ts / Dam$T1 - 1), 0.002)
  }
})

test_that("random layered wedges reduce to the exact solution of the layered equation of motion", {
  # Seed of the chain's validated contract (owner, 2026-10-12: semilla 20260915).
  set.seed(20260915L)
  for (Case in 1:40) {
    N <- sample(1:12, 1)
    hs <- exp(runif(N, log(0.5), log(15)))
    Layers <- data.frame(hs = hs, zm = cumsum(hs) - hs / 2, gsat = runif(N, 12, 24),
                         VSm = exp(runif(N, log(60), log(1500))))
    lo <- if (Case <= 5) 0 else runif(1, 0, 0.95)
    Hs <- sum(hs)
    OUT <- solveWedgePeriod(Layers, Hs = Hs, lo = lo, mo = 0.5, VSo = 300, method = "exact", elements = 4)
    Exact <- wedgeExactPeriod(lo * Hs / (1 - lo), hs, Layers$VSm, Layers$gsat, upper = 1.5 * 2 * pi / OUT$Ts)
    expect_lte(abs(OUT$Ts / Exact - 1), OUT$TsMeshChange)
  }
})

test_that("the analytic route carries its root diagnostics and the geometry is checked", {
  Layers <- powerLawLayers(10, 0.351, 0.5, 200, 20)
  Exact <- solveWedgePeriod(Layers, Hs = 10, lo = 0.351, mo = 0.5, VSo = 200, method = "exact", elements = 2)
  Surrogate <- solveWedgePeriod(Layers, Hs = 10, lo = 0.351, mo = 0.5, VSo = 200, method = "nlm", elements = 2)
  # The numerical period does not depend on the root model of the analytic route.
  expect_identical(Surrogate$Ts, Exact$Ts)
  expect_true(Exact$anExact && Exact$anWithinSupport)
  Root <- evaluateRoot(buildRootModel(mode = 1, method = "nlm"), mo = 0.5, lo = 0.351)
  expect_identical(Surrogate[, .(an, anExact, anWithinSupport)], Root[, .(an, anExact = exact, anWithinSupport = withinSupport)])
  expect_equal(Surrogate$TsAnalytic, evaluateBeamPeriod(Hs = 10, lo = 0.351, mo = 0.5, VSo = 200, an = Root$an))
  expect_equal(Surrogate$TsAnalyticDifference, Surrogate$TsAnalytic / Surrogate$Ts - 1)
  expect_error(solveWedgePeriod(Layers, Hs = 12, lo = 0.351, mo = 0.5, VSo = 200, method = "exact", elements = 2),
               "add up to 10 m, not to Hs = 12 m")
  for (lo in list(-0.1, 1, NA_real_, c(0.2, 0.3))) {
    expect_error(solveWedgePeriod(Layers, Hs = 10, lo = lo, mo = 0.5, VSo = 200, method = "exact", elements = 2),
                 "lo must be one truncation ratio")
  }
  expect_error(solveWedgePeriod(Layers, Hs = -10, lo = 0.3, mo = 0.5, VSo = 200, method = "exact", elements = 2),
               "Hs must be one finite positive height")
  expect_error(solveWedgePeriod(Layers, Hs = 10, lo = 0.351, mo = 2, VSo = 200, method = "exact", elements = 2),
               "mo in")
  expect_error(solveWedgePeriod(Layers, Hs = 10, lo = 0.351, mo = 0.5, VSo = 0, method = "nlm", elements = 2),
               "VSo")
  expect_error(solveWedgePeriod(Layers, Hs = 10, lo = 0.351, mo = c(0.5, 0.6), VSo = 200, method = "exact",
                                elements = 2), "one number each")
  expect_error(solveWedgePeriod(Layers, Hs = 10, lo = 0.351, mo = 0.5, VSo = 200, method = "rf", elements = 2),
               "lm")
})
