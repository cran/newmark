# Power-law rock curves AEP = 0.001 Sa^-2 at Tn = 0, 1 and 2 s, the PGA
# curve doubled at 1 s, so that 1.5 Ts lies within them for Ts up to 4/3 s.
.powerCurves <- function() {
  PGA <- exp(seq(log(.001), log(3), length.out = 30L))
  rbind(data.frame(Tn = 0, Sa = PGA, AEP = .001 * PGA^-2), data.frame(Tn = 1, Sa = 2 * PGA, AEP = .001 * PGA^-2),
        data.frame(Tn = 2, Sa = PGA, AEP = .001 * PGA^-2))
}

# The arguments shared by integrateDnHazard() and integrateDnHazardProfiles().
.mixtureArgs <- function(weights = c(AM88 = 1, BT07 = 3)) {
  list(curves = .powerCurves(), Mw = 8, weights = weights, Da = 100, IT = 5, TR = c(200, 100000), siteModel = "ST20",
       Vref = 760)
}

.integrateOne <- function(args, ky, Ts, VS30) {
  integrateDnHazard(args$curves, ky = ky, Ts = Ts, VS30 = VS30, Vref = args$Vref, Mw = args$Mw, weights = args$weights,
                    Da = args$Da, IT = args$IT, TR = args$TR, siteModel = args$siteModel)
}

.mix <- function(args, profiles, ky) do.call(integrateDnHazardProfiles, c(args, list(profiles = profiles, ky = ky)))

SUMMARYCOLUMNS <- c("Da", "IT", "AEP", "POE", "bandAEP", "bandPOE", "rareTailAEPBound", "rareTailPOEBound")

test_that("one profile and one ky of weight one reproduce integrateDnHazard", {
  Args <- .mixtureArgs()
  Direct <- .integrateOne(Args, ky = .2, Ts = 1, VS30 = 217)
  OUT <- .mix(Args, profiles = data.frame(Ts = 1, VS30 = 217), ky = .2)
  expect_identical(names(OUT), c("summary", "components"))
  expect_identical(names(OUT$summary), c(SUMMARYCOLUMNS, "nProfiles", "nKy", "extrapolated", "mixtureError",
                                         "bandMixtureError"))
  expectReproduces(as.list(OUT$summary[, ..SUMMARYCOLUMNS]), as.list(Direct$summary[, ..SUMMARYCOLUMNS]))
  expect_equal(as.list(OUT$summary[, .(nProfiles, nKy, extrapolated, mixtureError, bandMixtureError)]),
                   list(nProfiles = 1L, nKy = 1L, extrapolated = 0, mixtureError = NA_real_, bandMixtureError = NA_real_))
  expect_identical(names(OUT$components), c("Ts", "VS30", "ky", "w", SUMMARYCOLUMNS, "extrapolated"))
  # Input copies: the component keys are the profile and ky inputs unchanged.
  expect_identical(as.list(OUT$components[, .(Ts, VS30, ky)]), list(Ts = 1, VS30 = 217, ky = .2))
  expect_equal(OUT$components$w, 1)
  expectReproduces(as.list(OUT$components[, ..SUMMARYCOLUMNS]), as.list(Direct$summary[, ..SUMMARYCOLUMNS]))
  # ky as a table of one row with any positive weight is the same mixture.
  expectReproduces(.mix(Args, profiles = data.frame(Ts = 1, VS30 = 217, w = 7), ky = data.frame(ky = .2, w = 3)), OUT)
})

test_that("equal profiles reproduce the rate of one, with no Monte Carlo spread", {
  Args <- .mixtureArgs()
  One <- .mix(Args, profiles = data.frame(Ts = 1, VS30 = 217), ky = .2)
  Two <- .mix(Args, profiles = data.frame(Ts = c(1, 1), VS30 = c(217, 217)), ky = .2)
  expectReproduces(as.list(Two$summary[, ..SUMMARYCOLUMNS]), as.list(One$summary[, ..SUMMARYCOLUMNS]))
  expect_identical(Two$summary$nProfiles, 2L)
  expect_equal(c(Two$summary$mixtureError, Two$summary$bandMixtureError), c(0, 0))
  expect_equal(Two$components$w, c(.5, .5))
})

test_that("the rigid models do not depend on Ts, so profiles that differ only in Ts give the rate of one", {
  skip_on_cran()
  Args <- .mixtureArgs(weights = c(AM88 = 1, JB07 = 1))
  One <- .mix(Args, profiles = data.frame(Ts = .3, VS30 = 217), ky = .2)
  Spread <- .mix(Args, profiles = data.frame(Ts = c(.3, 1.2), VS30 = 217), ky = .2)
  expectReproduces(Spread$components$AEP, rep(One$summary$AEP, 2L))
  expectReproduces(as.list(Spread$summary[, ..SUMMARYCOLUMNS]), as.list(One$summary[, ..SUMMARYCOLUMNS]))
})

test_that("the mixture is the weighted sum of the conditional rates, and Poisson acts after it", {
  skip_on_cran()
  Args <- .mixtureArgs()
  Profiles <- data.frame(Ts = c(.6, 1.2), VS30 = c(300, 217), w = c(1, 3))
  Ky <- data.frame(ky = c(.15, .25), w = c(2, 1))
  OUT <- .mix(Args, profiles = Profiles, ky = Ky)
  Grid <- expand.grid(k = 1:2, i = 1:2)
  Direct <- lapply(seq_len(nrow(Grid)), function(j) {
    .integrateOne(Args, ky = Ky$ky[Grid$k[j]], Ts = Profiles$Ts[Grid$i[j]], VS30 = Profiles$VS30[Grid$i[j]])$summary
  })
  w <- (Profiles$w / sum(Profiles$w))[Grid$i] * (Ky$w / sum(Ky$w))[Grid$k]
  expect_equal(OUT$components$w, w, tolerance = 1e-15)
  # Input copies: Ts and ky are the inputs of the grid, not results.
  expect_identical(OUT$components$Ts, Profiles$Ts[Grid$i])
  expect_identical(OUT$components$ky, Ky$ky[Grid$k])
  expectReproduces(OUT$components$AEP, vapply(Direct, `[[`, 0, "AEP"))
  Rate <- sum(w * vapply(Direct, `[[`, 0, "AEP"))
  Band <- sum(w * vapply(Direct, `[[`, 0, "bandAEP"))
  Tail <- sum(w * vapply(Direct, `[[`, 0, "rareTailAEPBound"))
  expect_equal(OUT$summary$AEP, Rate, tolerance = 1e-14)
  expect_equal(OUT$summary$bandAEP, Band, tolerance = 1e-14)
  expect_equal(OUT$summary$rareTailAEPBound, Tail, tolerance = 1e-14)
  expectReproduces(OUT$summary$POE, -expm1(-5 * OUT$summary$AEP))
  expectReproduces(OUT$summary$bandPOE, -expm1(-5 * OUT$summary$bandAEP))
  expectReproduces(OUT$summary$rareTailPOEBound, -expm1(-5 * OUT$summary$rareTailAEPBound))
  # The mixed probability is not the mixture of the probabilities.
  expect_gt(abs(OUT$summary$POE - sum(w * vapply(Direct, `[[`, 0, "POE"))), 0)
  # The Monte Carlo error of the profiles from their ky-mixed rates, with
  # the weights 1/4 and 3/4: effective count 1.6, divisor 1 - 10/16.
  r <- tapply(w * vapply(Direct, `[[`, 0, "AEP"), Grid$i, sum) / c(.25, .75)
  m <- sum(c(.25, .75) * r)
  expect_equal(OUT$summary$mixtureError, sqrt(sum(c(.25, .75) * (r - m)^2) / (1 - 10 / 16) * 10 / 16), tolerance = 1e-12)
})

test_that("with equal weights the mixture error is the standard error of the mean rate of the profiles", {
  skip_on_cran()
  Args <- .mixtureArgs()
  Profiles <- data.frame(Ts = c(.5, .8, 1.1), VS30 = c(320, 260, 217))
  OUT <- .mix(Args, profiles = Profiles, ky = .2)
  expect_equal(OUT$summary$mixtureError, stats::sd(OUT$components$AEP) / sqrt(3), tolerance = 1e-12)
  expect_equal(OUT$summary$bandMixtureError, stats::sd(OUT$components$bandAEP) / sqrt(3), tolerance = 1e-12)
  expect_equal(OUT$summary$AEP, mean(OUT$components$AEP), tolerance = 1e-14)
})

test_that("the weights of each set are normalized, and profiles without w weigh equally", {
  skip_on_cran()
  Args <- .mixtureArgs()
  Base <- .mix(Args, profiles = data.frame(Ts = c(.6, 1.2), VS30 = c(300, 217), w = c(.25, .75)),
               ky = data.frame(ky = c(.15, .25), w = c(2 / 3, 1 / 3)))
  Scaled <- .mix(Args, profiles = data.frame(Ts = c(.6, 1.2), VS30 = c(300, 217), w = c(2, 6)),
                 ky = data.frame(ky = c(.15, .25), w = c(10, 5)))
  expectReproduces(Scaled, Base)
  Equal <- .mix(Args, profiles = data.frame(Ts = c(.6, 1.2), VS30 = c(300, 217)), ky = .2)
  Ones <- .mix(Args, profiles = data.frame(Ts = c(.6, 1.2), VS30 = c(300, 217), w = c(5, 5)), ky = .2)
  expectReproduces(Equal, Ones)
})

test_that("extrapolated is the weighted fraction of the components where a model is outside its ky range", {
  skip_on_cran()
  Args <- .mixtureArgs(weights = c(BT07 = 1))
  OUT <- .mix(Args, profiles = data.frame(Ts = 1, VS30 = 217), ky = data.frame(ky = c(.2, .5), w = c(3, 1)))
  Flags <- vapply(c(.2, .5), function(ky) {
    as.numeric(any(.integrateOne(Args, ky = ky, Ts = 1, VS30 = 217)$models$extrapolated == 1))
  }, 0)
  expect_equal(Flags, c(0, 1))
  # The product assigns the discrete flags 0 and 1, not numerical estimates.
  expect_identical(OUT$components$extrapolated, Flags)
  expect_equal(OUT$summary$extrapolated, .25)
})

test_that("malformed profiles, ky and weights are refused before any integration", {
  Args <- .mixtureArgs()
  expect_error(.mix(Args, profiles = data.frame(Ts = 1), ky = .2), "profiles must contain rows with Ts and VS30", fixed = TRUE)
  expect_error(.mix(Args, profiles = data.frame(Ts = numeric(), VS30 = numeric()), ky = .2),
               "profiles must contain rows with Ts and VS30", fixed = TRUE)
  expect_error(.mix(Args, profiles = data.frame(Ts = "1", VS30 = 217), ky = .2), "Ts and VS30 of profiles must be numeric",
               fixed = TRUE)
  for (w in list(c(1, 0), c(1, NA), c(1, Inf), c(1, -1))) {
    expect_error(.mix(Args, profiles = data.frame(Ts = c(1, 1), VS30 = 217, w = w), ky = .2),
                 "the weights w of profiles must be positive and finite", fixed = TRUE)
  }
  for (ky in list(c(.1, .2), data.frame(ky = .2), data.frame(ky = "0.2", w = 1), "0.2", data.frame(ky = numeric(), w = numeric()))) {
    expect_error(.mix(Args, profiles = data.frame(Ts = 1, VS30 = 217), ky = ky),
                 "ky must be one number or a table of numeric ky with weights w", fixed = TRUE)
  }
  expect_error(.mix(Args, profiles = data.frame(Ts = 1, VS30 = 217), ky = data.frame(ky = .2, w = 0)),
               "the weights w of ky must be positive and finite", fixed = TRUE)
  # A matrix column of a data.frame is numeric but holds several values per
  # row, which would multiply the components and their weights.
  Profiles <- data.frame(Ts = c(.6, 1.2), VS30 = c(300, 217), w = c(1, 3))
  Ky <- data.frame(ky = c(.15, .25), w = c(2, 1))
  Messages <- c(Ts = "Ts and VS30 of profiles must be numeric, one value per row",
                VS30 = "Ts and VS30 of profiles must be numeric, one value per row",
                w = "the weights w of profiles must be positive and finite, one per row")
  for (Column in names(Messages)) {
    Matrix <- Profiles
    Matrix[[Column]] <- cbind(Profiles[[Column]], Profiles[[Column]])
    expect_error(.mix(Args, profiles = Matrix, ky = Ky), Messages[[Column]], fixed = TRUE)
  }
  Messages <- c(ky = "ky must be one number or a table of numeric ky with weights w, one value per row",
                w = "the weights w of ky must be positive and finite, one per row")
  for (Column in names(Messages)) {
    Matrix <- Ky
    Matrix[[Column]] <- cbind(Ky[[Column]], Ky[[Column]])
    expect_error(.mix(Args, profiles = Profiles, ky = Matrix), Messages[[Column]], fixed = TRUE)
  }
})

test_that("a weight near one keeps the mixture error of the reliability weights", {
  # With the weights 1 and 1e-20, 1 - sum(w^2) rounds to 0 while the divisor
  # of reliability weights is 2e-20. For two profiles the error is
  # |r1 - r2| sqrt(sum(w^2) / 2), here |r1 - r2| / sqrt(2).
  expect_equal(.measureMixtureError(c(1, 2), c(1, 1e-20)), sqrt(.5), tolerance = 1e-12)
  OUT <- .mix(.mixtureArgs(), profiles = data.frame(Ts = c(.6, 1.2), VS30 = c(300, 217), w = c(1, 1e-20)), ky = .2)
  expect_equal(OUT$summary$mixtureError, abs(diff(OUT$components$AEP)) / sqrt(2), tolerance = 1e-12)
  expect_equal(OUT$summary$bandMixtureError, abs(diff(OUT$components$bandAEP)) / sqrt(2), tolerance = 1e-12)
  expect_gt(OUT$summary$mixtureError, 0)
})
