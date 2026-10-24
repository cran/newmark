# Numerical site period of the layered column (solveColumnPeriod). The
# expected values are closed forms of the hydrated sources evaluated here,
# independently of the solver; "within the discretisation error" means that
# the distance to the closed form does not exceed the solver's own mesh
# change (for a uniform profile) or the change between a profile sampled on
# N and on 2N layers (for a continuous profile represented by layers).

layerTable <- function(hs, VSm, gsat) data.frame(hs = hs, zm = cumsum(hs) - hs / 2, gsat = gsat, VSm = VSm)

# First period of the n-layer frequency equation (E) of Charara (2019),
# eq-frequency-equation-condensed: from x = 0, x <- r_k tan(t_k w + atan(x))
# over the upper layers, equal to cot(t_n w), with t_k = H_k / V_k and
# r_k = rho_k V_k / (rho_(k+1) V_(k+1)). Between poles the difference
# increases with w, so the fundamental is its first crossing from negative
# to positive; `upper` bounds the search in rad/s.
chararaPeriod <- function(hs, VSm, gsat, upper) {
  n <- length(hs)
  t <- hs / VSm
  r <- (gsat * VSm)[-n] / (gsat * VSm)[-1L]
  equationE <- function(w) {
    x <- 0
    for (k in seq_len(n - 1L)) x <- r[k] * tan(t[k] * w + atan(x))
    x - 1 / tan(t[n] * w)
  }
  Grid <- seq(upper / 20000, upper, length.out = 20000)
  Value <- vapply(Grid, equationE, numeric(1L))
  k <- which(Value[-length(Value)] < 0 & Value[-1L] > 0)[1L]
  2 * pi / uniroot(equationE, Grid[c(k, k + 1L)], tol = upper * 1e-14)$root
}

# The same exact solution without the poles of equation (E), whose root can
# fall within a grid cell of a pole when a thin stiff layer lies under a soft
# one: in each layer U = A cos(w z / V_k + phi_k) (Charara, eq-layer-solution),
# so displacement and shear G U' are carried from the free surface (1 and 0)
# across each layer and stay continuous at the interfaces
# (eq-interface-displacement, eq-interface-stress); the fundamental is the
# first zero of the displacement at the base.
columnExactPeriod <- function(hs, VSm, gsat, upper) {
  G <- gsat * VSm^2
  baseDisplacement <- function(w) {
    U <- 1
    Shear <- 0
    for (k in seq_along(hs)) {
      Phase <- w * hs[k] / VSm[k]
      Impedance <- G[k] * w / VSm[k]
      AUX <- U * cos(Phase) + Shear / Impedance * sin(Phase)
      Shear <- Shear * cos(Phase) - U * Impedance * sin(Phase)
      U <- AUX
    }
    U
  }
  Grid <- seq(upper / 600, upper, length.out = 600)
  Value <- vapply(Grid, baseDisplacement, numeric(1L))
  k <- which(Value[-1L] * Value[-length(Value)] < 0)[1L]
  2 * pi / uniroot(baseDisplacement, Grid[c(k, k + 1L)], tol = upper * 1e-14)$root
}

test_that("a uniform column converges to 4H/V from below within its mesh change", {
  # Dobry, Oweis and Urzua (1976), eq-uniform-period: T = 4H/V.
  for (Layers in list(layerTable(20, 200, 18), layerTable(rep(4, 5), rep(200, 5), rep(18, 5)))) {
    for (Elements in c(1, 4, 16)) {
      OUT <- solveColumnPeriod(Layers, elements = Elements)
      expect_gte(4 * 20 / 200 / OUT$Tr - 1, 0)
      expect_lte(4 * 20 / 200 / OUT$Tr - 1, OUT$TrMeshChange)
    }
  }
})

test_that("the power-law stratum of Dobry, Oweis and Urzua reduces to its closed form", {
  # V = VSo (z/Hs)^(p/2) from the surface, constant density: eq-power-period
  # with the height exponent (2 - p)/2, T = 4 pi Hs / ((2 - p) VSo q1), q1 the
  # first zero of J_nu, nu = (p - 1)/(2 - p); the printed exponent (2 - p/2)
  # would give 4 H^2 / V at p = 0 instead of eq-uniform-period's 4 H / V. The
  # layers sample V at their midpoints.
  Hs <- 30
  VSo <- 300
  for (p in c(0, 0.5, 1, 1.5)) {
    Nu <- (p - 1) / (2 - p)
    Grid <- seq(0.05, 10, by = 0.01)
    Value <- besselJ(Grid, Nu)
    k <- which(Value[-1L] * Value[-length(Value)] < 0)[1L]
    q1 <- uniroot(function(x) besselJ(x, Nu), Grid[c(k, k + 1L)], tol = 1e-14)$root
    Closed <- 4 * pi * Hs / ((2 - p) * VSo * q1)
    Period <- vapply(c(100, 200), function(N) {
      hs <- rep(Hs / N, N)
      solveColumnPeriod(layerTable(hs, VSo * ((cumsum(hs) - hs / 2) / Hs)^(p / 2), rep(19, N)), elements = 1)$Tr
    }, numeric(1L))
    expect_lte(abs(Period[2] / Closed - 1), abs(Period[2] / Period[1] - 1))
  }
})

test_that("the profiles of Charara reproduce his exact n-layer periods and Dobry's two-layer equation", {
  # Charara (2019), section 6.1 and Table 1 (printed page 1703): G = 0.5 Gmax with Gmax 40 (Clay 1),
  # 80 (Clay 2) and 180 MPa (sands); densities 1.8 and 2.0 t/m3; exact periods
  # printed to two decimals: 0.98, 0.79 and 0.82 s.
  Profiles <- list(list(hs = c(24, 8), G = c(20, 90), rho = c(1.8, 2.0), printed = 0.98),
                   list(hs = c(12, 12, 8), G = c(20, 40, 90), rho = c(1.8, 1.8, 2.0), printed = 0.79),
                   list(hs = c(12, 4, 12, 4), G = c(20, 90, 40, 90), rho = c(1.8, 2.0, 1.8, 2.0), printed = 0.82))
  for (Profile in Profiles) {
    VSm <- sqrt(Profile$G * 1000 / Profile$rho)
    OUT <- solveColumnPeriod(layerTable(Profile$hs, VSm, 9.81 * Profile$rho), elements = 50)
    Exact <- chararaPeriod(Profile$hs, VSm, Profile$rho, upper = 1.5 * 2 * pi / OUT$Tr)
    expect_equal(columnExactPeriod(Profile$hs, VSm, Profile$rho, upper = 1.5 * 2 * pi / OUT$Tr), Exact,
                 tolerance = 1e-10)
    expect_lte(abs(OUT$Tr / Exact - 1), OUT$TrMeshChange)
    expect_equal(round(Exact, 2), Profile$printed)
    expect_equal(round(OUT$Tr, 2), Profile$printed)
    expect_lte(OUT$TrRayleighWeightedDifference, OUT$TrMeshChange)
  }
  # Two layers: the period solves eq-twolayer-period-equation of Dobry, Oweis
  # and Urzua, tan(pi TA / 2T) tan(pi TB / 2T) = (rhoB HB / rhoA HA)(TA / TB).
  VSm <- sqrt(c(20, 90) * 1000 / c(1.8, 2.0))
  Exact <- chararaPeriod(c(24, 8), VSm, c(1.8, 2.0), upper = 10)
  TA <- 4 * 24 / VSm[1]
  TB <- 4 * 8 / VSm[2]
  expect_equal(tan(pi * TA / (2 * Exact)) * tan(pi * TB / (2 * Exact)), (2.0 * 8) / (1.8 * 24) * TA / TB,
               tolerance = 1e-9)
})

test_that("the period of Dobry's Table 4 profile is bounded below by its printed Rayleigh period", {
  # Dobry, Oweis and Urzua (1976), Fig. 11 profile and Table 4 (printed page 1310; 12 layers from
  # the bottom up, ft and ft/s, constant density): Method 7 gives 0.593 s, a
  # lower bound of the exact period (section Rayleigh bound of the manual).
  H <- c(25, 21, 7, 9, 16, 24, 65, 5, 8, 7, 5, 8)
  V <- c(2777, 1223, 1633, 1071, 1205, 1065, 1270, 951, 825, 1039, 726, 828)
  Depth <- c(187.5, 164.5, 150.5, 142.5, 130, 110, 65.5, 30.5, 24, 16.5, 10.5, 4)
  OUT <- solveColumnPeriod(data.frame(hs = H, zm = Depth, gsat = 120, VSm = V), elements = 16)
  expect_equal(round(OUT$TrRayleighWeighted, 3), 0.593)
  expect_lte(OUT$TrRayleighWeightedDifference, OUT$TrMeshChange)
  # Charara's equation (E) on the same layers, from the surface down.
  Exact <- chararaPeriod(rev(H), rev(V), rep(120, 12), upper = 1.5 * 2 * pi / OUT$Tr)
  expect_lte(abs(OUT$Tr / Exact - 1), OUT$TrMeshChange)
  expect_gt(Exact, 0.593)
})

test_that("random layered columns reduce to the exact layered solution", {
  # Seed of the chain's validated contract (owner, 2026-10-12: semilla 20260915).
  set.seed(20260915L)
  for (Case in 1:40) {
    N <- sample(1:12, 1)
    Layers <- layerTable(exp(runif(N, log(0.2), log(20))), exp(runif(N, log(50), log(2000))), runif(N, 12, 24))
    OUT <- solveColumnPeriod(Layers, elements = 4)
    Exact <- columnExactPeriod(Layers$hs, Layers$VSm, Layers$gsat, upper = 1.5 * 2 * pi / OUT$Tr)
    expect_lte(abs(OUT$Tr / Exact - 1), OUT$TrMeshChange)
  }
})

test_that("a table with Gm and no VSm uses the velocities of evaluateSoilLayers", {
  Layers <- layerTable(c(5, 10, 15), c(180, 260, 340), c(17, 19, 21))
  Moduli <- data.frame(hs = Layers$hs, zm = Layers$zm, gsat = Layers$gsat, Gm = Layers$VSm^2 * Layers$gsat / 9810)
  expect_equal(solveColumnPeriod(Moduli, elements = 4), solveColumnPeriod(Layers, elements = 4))
  # With both columns the velocities are used and the moduli ignored.
  Layers$Gm <- -1
  expect_identical(solveColumnPeriod(Layers, elements = 4), solveColumnPeriod(Layers[, -5], elements = 4))
})

test_that("degenerate layers and meshes are rejected with a diagnostic", {
  Layers <- layerTable(c(5, 10), c(180, 260), c(17, 19))
  expect_error(solveColumnPeriod(Layers[0, ], elements = 4), "columns hs, zm, gsat")
  expect_error(solveColumnPeriod(Layers[, c("hs", "zm", "gsat")], elements = 4), "VSm or Gm")
  expect_error(solveColumnPeriod(list(hs = 1, zm = 0.5, gsat = 18, VSm = 200), elements = 4), "columns")
  for (Column in c("hs", "gsat", "VSm", "zm")) {
    for (Value in list(0, -1, NA_real_, Inf)) {
      AUX <- Layers
      AUX[[Column]][2] <- Value
      expect_error(solveColumnPeriod(AUX, elements = 4), paste0("layers\\$", Column, " must be finite and positive"))
    }
  }
  Moduli <- data.frame(hs = c(5, 10), zm = c(2.5, 10), gsat = c(17, 19), Gm = c(50, 0))
  expect_error(solveColumnPeriod(Moduli, elements = 4), "layers\\$Gm must be finite and positive")
  AUX <- Layers
  AUX$zm[2] <- AUX$zm[1]
  expect_error(solveColumnPeriod(AUX, elements = 4), "distinct")
  for (Elements in list(0, 1.5, c(2, 4), NA_real_, "4")) {
    expect_error(solveColumnPeriod(Layers, elements = Elements), "positive whole number")
  }
})
