# Supplied-layer calculations: the requested difference of this family. The
# installed newmark 1.7.1 has no entry for measured layers (its physics is
# reachable only through synthetic generation), so these checks fail on the
# oracle by construction. Expected values are closed-form evaluations of the
# documented equations, not outputs of either implementation.

test_that("estimateVS30 is the thickness-weighted harmonic average, completed below 30 m", {
  expect_equal(estimateVS30(hs = c(10, 20), VSm = c(200, 400), Vref = 760), 300)
  # 30 / (10/200 + 20/760): the completion term uses the reference velocity.
  expect_equal(estimateVS30(hs = 10, VSm = 200, Vref = 760), 30 / (10 / 200 + 20 / 760))
  # Layers below 30 m do not contribute.
  expect_equal(estimateVS30(hs = c(30, 50), VSm = c(300, 900), Vref = 760), 300)
  expect_error(estimateVS30(hs = c(10, -1), VSm = c(200, 300), Vref = 760), "finite positive")
})

test_that("estimateRayleighPeriod gives 3.974 H/Vs for a uniform column and ignores layer order", {
  # Documented closed form of the static-deflection Rayleigh quotient.
  Ts <- estimateRayleighPeriod(VSm = rep(200, 20), hs = rep(1, 20), zm = seq(0.5, 19.5, by = 1))
  expect_equal(Ts, 3.974 * 20 / 200, tolerance = 1e-3)
  IDX <- c(3, 1, 2)
  expect_identical(estimateRayleighPeriod(VSm = c(220, 280, 340)[IDX], hs = c(10, 10, 10)[IDX], zm = c(5, 15, 25)[IDX]),
                   estimateRayleighPeriod(VSm = c(220, 280, 340), hs = c(10, 10, 10), zm = c(5, 15, 25)))
})

test_that("the Rayleigh periods reproduce Dobry's Table 4 profile and reduce to each other", {
  # Dobry, Oweis and Urzua (1976), Fig. 11 profile (H = 200 ft, 12 layers from the
  # bottom up) and Table 4: Method 7 gives T7 = 0.593 s (hydrated copy, Table 4 and
  # the Fig. 11 text; the figure's exact period is also printed as 0.593 s).
  H <- c(25, 21, 7, 9, 16, 24, 65, 5, 8, 7, 5, 8)
  V <- c(2777, 1223, 1633, 1071, 1205, 1065, 1270, 951, 825, 1039, 726, 828)
  Depth <- c(187.5, 164.5, 150.5, 142.5, 130, 110, 65.5, 30.5, 24, 16.5, 10.5, 4)
  Ts <- estimateRayleighPeriod(VSm = V, hs = H, zm = Depth)
  expect_equal(round(Ts, 3), 0.593)
  # Method 7 as printed: X from the layer increments, omega7^2 from the two sums.
  X <- c(0, cumsum(Depth / V^2 * H))
  Omega2 <- 4 * sum(Depth^2 / V^2 * H) / sum((X[-1L] + X[-13L])^2 * H)
  expect_equal(Ts, 2 * pi / sqrt(Omega2))
  # A uniform unit weight reduces the weighted quotient to the constant-density one.
  expect_equal(estimateWeightedRayleighPeriod(VSm = V, hs = H, zm = Depth, gsat = rep(120, 12)), Ts)
  # Fig. 11 unit weights, 105 to 140 pcf from the surface down (figure text); no
  # published value exists for the weighted estimate, so only its invariances are checked.
  Gamma <- rev(c(105, 133, 120, 120, 137, 125, 127, 119, 130, 135, 138, 140))
  Tw <- estimateWeightedRayleighPeriod(VSm = V, hs = H, zm = Depth, gsat = Gamma)
  expect_true(is.finite(Tw) && Tw > 0)
  IDX <- c(12:7, 1:6)
  expect_identical(estimateWeightedRayleighPeriod(V[IDX], H[IDX], Depth[IDX], Gamma[IDX]), Tw)
  expect_identical(estimateWeightedRayleighPeriod(V, H, Depth, Gamma * 2), Tw)
  expect_error(estimateWeightedRayleighPeriod(V, H, Depth, Gamma[-1L]), "equal length")
})

test_that("classifySite follows ASCE 7-22 Table 20.2-1 with closed lower bounds", {
  expect_identical(classifySite(c(120, 150, 209.9, 210, 300, 440, 640, 910, 1500, 3000)),
                   c("E", "DE", "DE", "D", "CD", "C", "BC", "B", "A", "A"))
  expect_error(classifySite(-1), "non-negative")
})

test_that("evaluateSoilLayers evaluates the documented stress and stiffness relations on supplied layers", {
  Layers <- data.frame(USCS = c("SM", "ML"), hs = c(5, 5), zm = c(2.5, 7.5),
                       eo = c(0.7, 0.9), gsat = c(20, 18), IP = c(NA, 10))
  OUT <- evaluateSoilLayers(Layers, Hs = 10, water = 0, POP = 100)
  expect_identical(OUT$group, c("Sands", "Fines"))
  # Dry column: effective weight is the saturated weight; Ko = 0.5 gives
  # dp = (1 + 2 * 0.5) / 3 * gamma * h = gamma * h / 1.5 ... evaluated exactly:
  Dp <- 1 / 3 * (1 + 2 * 0.5) * c(20, 18) * c(5, 5)
  expect_equal(OUT$pm, c(Dp[1] / 2, Dp[1] + Dp[2] / 2))
  expect_equal(OUT$OCR, pmax((OUT$pm + 100) / OUT$pm, 1))
  # Sands: arithmetic mean over the group's correlations of A F(e) 100^n (pm/100)^n.
  Models <- newmark::ShearModelParameters[GroupID == "Sands"]
  Fe <- ifelse(Models$ModelID == "SS75", Models$Ce - 0.7 / 1.7, (Models$Ce - 0.7)^2 / 1.7)
  expect_equal(OUT$Gm[1], mean(Models$A * Fe * 100^Models$n * (OUT$pm[1] / 100)^Models$n))
  expect_equal(OUT$VSm, sqrt(9.81 * OUT$Gm * 1000 / c(20, 18)))
  # Below the water table the effective weight drops by 10 kN/m3.
  Wet <- evaluateSoilLayers(Layers, Hs = 10, water = 1, POP = 100)
  expect_equal(Wet$gs, c(10, 8))
  expect_error(evaluateSoilLayers(Layers[, -1], Hs = 10, water = 0, POP = 100), "columns")
})

test_that("m1 follows Table 1 of Hardin and Drnevich by default and the 1.7.1 nodes on request", {
  # Hardin and Drnevich (1972), Table 1 (hydrated copy, printed page 672): K is
  # 0.48 at PI 80 and 0.50 at PI >= 100, so m1 is 0.49 at IP 90 and 0.50 at
  # IP 120. The 1.7.1 nodes end at (100, 0.48): 0.48 at IP 90 and NA at 120.
  Layers <- data.frame(USCS = c("SM", "ML"), hs = c(5, 5), zm = c(2.5, 7.5),
                       eo = c(0.7, 0.9), gsat = c(20, 18), IP = c(NA, 90))
  Legacy <- data.frame(IP = c(0, 20, 40, 60, 80, 100), m1 = c(0, 0.18, 0.30, 0.41, 0.48, 0.48))
  # Closed-form modulus of the fine layer for a given m1; no Fines correlation
  # uses the linear form of SS75.
  Models <- newmark::ShearModelParameters[GroupID == "Fines"]
  Dp <- 1 / 3 * (1 + 2 * 0.5) * c(20, 18) * c(5, 5)
  pm <- Dp[1] + Dp[2] / 2
  OCR <- (pm + 100) / pm
  modulus <- function(m1) mean(Models$A * (Models$Ce - 0.9)^2 / 1.9 * OCR^m1 * 100^Models$n * (pm / 100)^Models$n)
  expect_equal(evaluateSoilLayers(Layers, Hs = 10, water = 0, POP = 100)$Gm[2], modulus(0.49))
  expect_equal(evaluateSoilLayers(Layers, Hs = 10, water = 0, POP = 100, m1 = Legacy)$Gm[2], modulus(0.48))
  Layers$IP[2] <- 120
  expect_equal(evaluateSoilLayers(Layers, Hs = 10, water = 0, POP = 100)$Gm[2], modulus(0.50))
  expect_identical(evaluateSoilLayers(Layers, Hs = 10, water = 0, POP = 100, m1 = Legacy)$Gm[2], NA_real_)
  expect_error(evaluateSoilLayers(Layers, Hs = 10, water = 0, POP = 100, m1 = Legacy[, "IP", drop = FALSE]),
               "numeric IP and m1 nodes")
})

test_that("characterizeProfile recovers an exact power-law stiffness profile", {
  # Gm = Go (zm/Hs)^mo exactly, so the log-log fit returns Go and mo up to rounding.
  zm <- seq(0.5, 19.5, by = 1)
  Gm <- 150 * (zm / 20)^0.4
  Layers <- data.frame(USCS = "SP", hs = 1, zm = zm, gsat = 20, Gm = Gm, VSm = sqrt(9.81 * Gm * 1000 / 20))
  OUT <- characterizeProfile(Layers, Hs = 20, Vref = 760)
  expect_equal(OUT$Go, 150)
  expect_equal(OUT$mo, 0.4)
  # The general route rounds nothing and carries no threshold depth ratios.
  expect_equal(OUT$VSo, sqrt(9.81 * 150 * 1000 / 20))
  expect_false(any(c("Z500", "Z1000") %in% names(OUT)))
  expect_lt(OUT$fitSigma, 1e-12)
  expect_identical(OUT$fitLayers, 20L)
  expect_identical(OUT$SID, classifySite(OUT$VS30))
  expect_equal(OUT$Sands, 100)
  expect_identical(OUT$USCS, "SP")
})
