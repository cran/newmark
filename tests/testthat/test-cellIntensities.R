# Owner's 2026-11-08 regression: the synthetic dn-reference spectra and
# seeds, sampled with NS=200 by the frozen pre-correction public samplers.
# The fixture contains the five paired intensities and the following eps.
test_that("100 identical columns reproduce all 200 intensities and eps within 1e-12", {
  X <- readRDS(test_path("fixtures", "dn-reference.rds"))$synthetic
  Expected <- readRDS(test_path("fixtures", "intensities-200.rds"))
  Evaluate <- evaluateDisplacements
  Observed <- list()
  local_mocked_bindings(evaluateDisplacements = function(intensities, eps, ...) {
    Observed[[length(Observed) + 1L]] <<- copy(intensities)[, eps := eps]
    Evaluate(intensities, eps = eps, ...)
  })
  for (i in seq_along(X$records)) {
    Records <- copy(X$records[[i]]$records)
    Columns <- data.table(realization = seq_len(100L), Ts = Records$Ts[[1L]],
                          spectrum = rep(Records$spectrum, 100L))
    Records[, profiles := list(Columns)]
    Records[, spectrum := NULL]
    Observed <- list()
    .simulateDnDemand(Records, demandID = if (i == 1L) X$plan$model else "MCE",
                       TR = if (i == 1L) X$plan$family else NA_real_, siteID = "SUK",
                       Mw = X$plan$mw, NS = 200L, weights = X$plan$weights,
                       seed = X$seeds[[i]], workers = 1L)
    DT <- rbindlist(Observed)[order(realization)]
    # Tolerance 1e-12: across the four failing R-CMD-check jobs of run 42531685559
    # (ubuntu release/devel/oldrel-1, windows release) the fixture made on macOS
    # differs from the sampled values by a median of 3e-16 relative, a 99th
    # percentile of 2.2e-15 and a maximum of 3.0e-14 (windows, PGV), from the
    # platform's libm and BLAS; the rejected per-column sampler deviates by 0.47
    # to 0.81 in the mean-relative metric of expect_equal. 1e-12 sits 33 times
    # above the worst observed platform discrepancy and 5e11 below the defect.
    expect_equal(as.list(DT[, !"eps"]), as.list(Expected[[i]][, !"eps"]),
                 tolerance = 1e-12, info = paste("intensities", i))
    expect_equal(DT$eps, Expected[[i]]$eps, tolerance = 1e-12, info = paste("eps", i))
  }
})
