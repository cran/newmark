test_that("column spectra share 1 m/s lookups while mean rows keep their velocity", {
  X <- readRDS(test_path("fixtures", "dn-reference.rds"))$pairs
  # S10/U1 mean row and columns 1, 4, 13 of the owner's seed-23 pilot.
  Shear <- data.table(IDg = "S10", IDm = "U1", Ts = 0.179, VS30 = 326.8, Vref = 760)
  Columns <- data.table(IDg = "S10", IDm = "U1", realization = 1:3,
                       Ts = c(0.170, 0.172, 0.186), VS30 = c(320.8, 323.9, 324.4), Vref = 760)
  Preserved <- copy(Columns)
  Mixture <- refineSiteMixture
  Spectrum <- refineSiteSpectrum
  Calls <- list()
  local_mocked_bindings(
    refineSiteMixture = function(...) {
      OUT <- Mixture(...)
      Calls[[length(Calls) + 1L]] <<- list(VS30 = list(...)$VS30, spectrum = OUT$mixture)
      OUT
    },
    refineSiteSpectrum = function(...) {
      OUT <- Spectrum(...)
      Calls[[length(Calls) + 1L]] <<- list(VS30 = list(...)$VS30,
        spectrum = OUT$spectrum[, .(Tn = as.numeric(Tn), p = as.character(p), Sa = as.numeric(SaF))])
      OUT
    })
  for (Period in c(100, NA_real_)) {
    Calls <- list()
    Rock <- if (is.na(Period)) as.data.table(X$inputs$MCE)[, .(Tn, p, Sa)] else
      as.data.table(X$inputs$rock)[TR == Period, .(Tn, k, muLn, sigmaLn, w)]
    Records <- .amplifyDemand(Rock, shear = Shear, draws = Columns, Vref = 760,
                               zTrunc = if (is.na(Period)) NA_real_ else 3,
                               demandID = if (is.na(Period)) "MCE" else "TEST", TR = Period)
    # Input copies: the velocities are the stored and the rounded inputs, and the
    # columns and the mean-row metadata pass through unchanged.
    expect_identical(vapply(Calls, `[[`, numeric(1L), "VS30"), c(326.8, 321, 324))
    expect_identical(Columns, Preserved)
    expect_identical(Records$profiles[[1L]][, !"spectrum"], Preserved)
    expect_identical(Records$Vs30, Shear$VS30)
    expect_identical(Records$Ts, Shear$Ts)
    expect_identical(Records$profiles[[1L]]$spectrum[[1L]], Calls[[2L]]$spectrum)
    expect_identical(Records$profiles[[1L]]$spectrum[[2L]], Records$profiles[[1L]]$spectrum[[3L]])
    Reference <- Calls[[1L]]$spectrum
    PGA <- if (is.na(Period)) Reference[Tn == 0 & p == "mean", Sa] else
      Reference[Tn == 0, Reduce(`+`, w * exp(muLn + sigmaLn^2 / 2))]
    expectReproduces(Records$ky[[1L]], buildKyGrid(PGA = PGA, n = 30L))
  }
})
