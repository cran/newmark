# Oracle: fixtures/siteResponse-reference.rds, written by the installed
# newmark 1.7.1 (build d57b64f) in its own R process, never together with
# newmark 5: the outputs of F_ST20 and SaF_ST17 on the examples of their
# manual pages.
siteReference <- function() readRDS(test_path("fixtures", "siteResponse-reference.rds"))

test_that("evaluateSiteFactor reproduces the ST20 factor of newmark 1.7.1 exactly", {
  Reference <- siteReference()
  Factor <- evaluateSiteFactor(PGA = c(0.05, 0.2, 0.5), Tn = 0.2, VS30 = 300, Vref = 760, model = "ST20")
  expect_identical(as.list(Factor), stats::setNames(as.list(Reference$factorST20), c("muLnF", "sigmaLnF", "ID")))
  Factor <- evaluateSiteFactor(PGA = 0.1, Tn = 1.0, VS30 = 760, Vref = 3000, model = "ST20")
  expect_identical(as.list(Factor), stats::setNames(as.list(Reference$factorST20at3000), c("muLnF", "sigmaLnF", "ID")))
})

test_that("evaluateSiteSa reproduces SaF_ST17 of newmark 1.7.1 exactly", {
  Site <- evaluateSiteSa(Sa = c(0.1, 0.2, 0.3, 0.4), PGA = rep(0.1, 4), Tn = 0.5, VS30 = 300, Vref = 760,
                         model = "ST17")
  expect_identical(as.list(Site), stats::setNames(as.list(siteReference()$siteSaST17), c("muLnSaF", "sigmaLnSaF", "ID")))
})

test_that("the site ordinate is the rock ordinate shifted by the factor", {
  Factor <- evaluateSiteFactor(PGA = c(0.05, 0.3), Tn = 0.3, VS30 = 250, Vref = 760, model = "ST20")
  Site <- evaluateSiteSa(Sa = c(0.2, 0.6), PGA = c(0.05, 0.3), Tn = 0.3, VS30 = 250, Vref = 760, model = "ST20")
  expect_identical(Site$muLnSaF, log(c(0.2, 0.6)) + Factor$muLnF)
  expect_identical(Site$sigmaLnSaF, Factor$sigmaLnF)
})

test_that("the reference condition is the unit factor, one row per PGA", {
  for (Model in c("ST20", "ST17")) {
    Factor <- evaluateSiteFactor(PGA = c(0.1, 0.2, 0.4), Tn = 0.5, VS30 = 760, Vref = 760, model = Model)
    expect_identical(Factor$muLnF, c(0, 0, 0))
    expect_identical(Factor$sigmaLnF, c(0, 0, 0))
  }
})

test_that("site factor inputs outside the models fail explicitly", {
  expect_error(evaluateSiteFactor(0.1, 0.2, 300, 760, "ST99"), "model must be")
  expect_error(evaluateSiteFactor(0.1, 0.2, 300, 500, "ST20"), "760 or 3000")
  expect_error(evaluateSiteFactor(-0.1, 0.2, 300, 760, "ST20"), "non-negative")
  expect_error(evaluateSiteFactor(0.1, c(0.2, 0.3), 300, 760, "ST20"), "one finite")
  expect_warning(evaluateSiteFactor(0.1, 12, 300, 760, "ST20"), "clamped to 10 s")
  expect_error(evaluateSiteSa(c(0.1, 0.2), 0.1, 0.2, 300, 760, "ST20"), "equal length")
})
