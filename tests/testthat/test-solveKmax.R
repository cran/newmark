# Cases of the kmax tests of newmark 1.7.1 (lib/tests/testthat at commit
# f7d11bb: test-kmax-boundary.R, test-kmax-continuation.R and
# test-kmax-mixture.R) rewritten for solveKmax(). The fixture
# kmax-closed-form.Rds holds the curves of the continuation and constant-AM88
# cases evaluated by the displacement equations Dn_* of the installed 1.7.1,
# with the exact roots its tests solve by uniroot over ky.
closedForm <- function() readRDS(test_path("fixtures", "kmax-closed-form.Rds"))

solveCase <- function(case, PGA = case$PGA) {
  solveKmax(ky = case$ky, Dn = matrix(case$Dn, nrow = 1L), PGA = PGA, model = case$model, Da = case$Da)
}

test_that("targets above the ceiling of a realization saturate at kmax zero", {
  ky <- c(0.05, 0.10, 0.20, 0.40)
  Roots <- solveKmax(ky = ky, Dn = matrix(0.02 / ky, nrow = 4L, ncol = 4L, byrow = TRUE), PGA = NULL,
                     model = "BT07", Da = c(0.1, 1))
  expect_equal(Roots[Da == 0.1, kmax], rep(0.2, 4L))
  expect_true(all(Roots[Da == 1, kmax] == 0))
  expect_identical(Roots$realization, rep(1:4, each = 2L))
  expect_false(any(Roots$beyond))
})

test_that("subnormal tails are floored and crossings resolve at the grid", {
  ky <- c(0.05, 0.10, 0.20, 0.40)
  Roots <- solveKmax(ky = ky, Dn = matrix(c(0.4, 0.02, 1e-300, 9.9e-324), nrow = 1L), PGA = NULL,
                     model = "BM17", Da = c(0.1, 0.001, 1))
  Interior <- exp(log(0.1) + (log(0.1) - log(0.02)) / (log(0.4) - log(0.02)) * (log(0.05) - log(0.1)))
  expect_equal(Roots$kmax, c(Interior, 0.2, 0))
})

test_that("roots below the drawn curve follow the closed form of each model", {
  Cases <- closedForm()
  for (model in c("BT07", "BM17", "BM19")) {
    Roots <- solveCase(Cases[[model]])
    expect_equal(Roots$kmax, Cases[[model]]$root, tolerance = 1e-6, info = model)
    expect_true(Roots$beyond, info = model)
    expect_gt(Roots$kmax, max(Cases[[model]]$ky))
  }
  for (model in c("AM88", "SR08")) {
    Roots <- solveCase(Cases[[model]])
    expect_equal(Roots$kmax, Cases[[model]]$root, tolerance = 1e-6, info = model)
    expect_lt(Roots$kmax, Cases[[model]]$PGA)
    expect_true(Roots$beyond, info = model)
  }
})

test_that("a large sampled PGA cannot produce an astronomical SR08 root", {
  Roots <- solveCase(closedForm()$SR08plateau)
  expect_lte(Roots$kmax, 36)
  expect_true(Roots$beyond)
})

test_that("roots within the drawn curve keep the grid interpolation", {
  Case <- closedForm()$BT07inside
  Roots <- solveCase(Case)
  expect_false(Roots$beyond)
  expect_true(Roots$kmax %between% range(Case$ky))
})

test_that("a constant AM88 curve is inverted exactly with its sampled PGA", {
  Case <- closedForm()$AM88flat
  Roots <- solveKmax(ky = Case$ky, Dn = matrix(Case$Dn, nrow = 1L, ncol = 3L), PGA = Case$PGA, model = "AM88",
                     Da = Case$Da)
  expect_equal(Roots$kmax, Case$root, tolerance = 1e-10)
  expect_true(all(Roots$kmax > 0 & Roots$kmax < 0.9999 * Case$PGA & Roots$beyond))
  expect_error(solveKmax(ky = Case$ky, Dn = matrix(Case$Dn, nrow = 1L, ncol = 3L), PGA = NULL, model = "AM88",
                         Da = Case$Da), "non-invertible")
  expect_error(solveKmax(ky = Case$ky, Dn = matrix(Case$Dn, nrow = 1L, ncol = 3L), PGA = 0.5, model = "AM88",
                         Da = Case$Da), "non-invertible")
  # Under another model the sub-micron constant curve is floored to zero and
  # every target saturates at null demand.
  Roots <- solveKmax(ky = Case$ky, Dn = matrix(Case$Dn, nrow = 1L, ncol = 3L), PGA = Case$PGA, model = "BT07",
                     Da = Case$Da)
  expect_true(all(Roots$kmax == 0))
})

test_that("a curve flattened by its envelope saturates above its ceiling and fails below", {
  Dn <- matrix(c(1, 2), nrow = 1L)
  expect_true(all(solveKmax(ky = c(0.01, 0.02), Dn = Dn, PGA = 0.005, model = "AM88", Da = 3)$kmax == 0))
  expect_error(solveKmax(ky = c(0.01, 0.02), Dn = Dn, PGA = 0.005, model = "AM88", Da = 1.5), "non-invertible")
  expect_error(solveKmax(ky = c(0.1, 0.2), Dn = matrix(1, nrow = 3L, ncol = 2L), PGA = NULL, model = "AM88",
                         Da = 0.5), "non-invertible")
})

test_that("invalid realizations fail explicitly", {
  ky <- c(0.05, 0.10)
  expect_error(solveKmax(ky = rev(ky), Dn = matrix(1, 1L, 2L), PGA = NULL, model = "BT07", Da = 1), "increasing")
  expect_error(solveKmax(ky = ky, Dn = matrix(-1, 1L, 2L), PGA = NULL, model = "BT07", Da = 1), "non-negative")
  expect_error(solveKmax(ky = ky, Dn = matrix(1, 1L, 2L), PGA = c(1, 1), model = "BT07", Da = 1), "per realization")
  expect_error(solveKmax(ky = ky, Dn = matrix(1, 1L, 2L), PGA = NULL, model = "XX01", Da = 1), "model must be")
  expect_error(solveKmax(ky = ky, Dn = matrix(1, 1L, 2L), PGA = NULL, model = "BT07", Da = c(1, 1)), "unique")
})
