# Oracle: fixtures/dn-reference.rds, written by newmark 1.7.1 in its own
# process: the scenario of test-sampleSaCorr.R (seed 48), the dn products of
# the synthetic project of helper-syntheticProject.R (SHA-256
# 725ea79f4fb53681...) run by runProcess(steps = "dn"), with the amplified
# records of each demand, and fitDnProducts() and getDnTables() on the
# inputs of test-product-pairs.R (SHA-256 f1012729ba0324ec...) under the
# product seeds of that site.
reference <- function() readRDS(test_path("fixtures", "dn-reference.rds"))

# The object as it is once saved and read back, as the products are compared.
roundTrip <- function(x) unserialize(serialize(x, NULL))

csvLines <- function(x) {
  FILE <- tempfile(fileext = ".csv")
  on.exit(unlink(FILE), add = TRUE)
  data.table::fwrite(x, FILE, na = "NA")
  readLines(FILE)
}

# The dn products of one site; the PGATable is rebuilt by the site-response
# family's buildPGATable() from the same inputs.
expectProducts <- function(site, global, pga, expected) {
  expectReproduces(roundTrip(site$DnDraws), expected$DnDraws)
  expectReproducesLines(csvLines(pga), expected$PGATable)
  expectReproduces(roundTrip(site$DnTable), expected$DnTable)
  expectReproduces(roundTrip(global$DnTable), expected$GlobalDnTable)
  expectReproduces(roundTrip(global$DnPlotTable), expected$DnPlotTable)
}

test_that("the scenario realizations and summaries reproduce the 1.7.1 scenario", {
  R <- reference()$scenario48
  Mixture <- data.table::data.table(Tn = c(0, 0.4, 1), k = "oq", sigmaLn = c(0.20, 0.45, 0.75), w = 1)
  Mixture[, muLn := log(c(0.80, 0.55, 0.42)) - sigmaLn^2 / 2]
  set.seed(48)
  Intensities <- sampleIntensities(Mixture, Ts = 0.4, NS = 7L, zTrunc = NA)
  Eps <- stats::rnorm(7L)
  Draws <- evaluateDisplacements(Intensities, eps = Eps, ky = c(0.01, 0.10), Ts = 0.4, Mw = 7,
                                 models = "AM88", NFC = "D100", pulse = FALSE)
  expectReproduces(Draws$Dn, R$draws$Dn)
  # The caller's realization identifiers are kept, paired with their rows.
  Relabeled <- data.table::copy(Intensities)[, realization := rev(realization) + 100L]
  Kept <- evaluateDisplacements(Relabeled, eps = Eps, ky = c(0.01, 0.10), Ts = 0.4, Mw = 7,
                                models = "AM88", NFC = "D100", pulse = FALSE)
  expect_identical(Kept$realization, rep(Relabeled$realization, 2L))
  expect_identical(Kept$Dn, Draws$Dn)
  for (DT in list(data.table::copy(Intensities)[, realization := realization[[1L]]], data.table::copy(Intensities)[2L, realization := NA_integer_])) {
    expect_error(evaluateDisplacements(DT, eps = Eps, ky = 0.1, Ts = 0.4, Mw = 7, models = "AM88",
                                       NFC = "D100", pulse = FALSE), "non-missing unique realization")
  }
  expectReproduces(Intensities$PGA[Draws$realization], R$draws$PGA.s)
  Levels <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95, "mean")
  Curve <- summarizeDisplacements(Draws, levels = Levels)
  expect_identical(Curve$level, R$curve$p[R$curve$IDn == "AM88"])
  expectReproduces(Curve$Dn, R$curve$Dn[R$curve$IDn == "AM88"])
  # The ensemble of 1.7.1 is the pool of the compatibility route; the
  # general summary coincides with it here because no pool holds a tie.
  Ensemble <- R$curve$Dn[R$curve$IDn == "ensemble"]
  expect_false(any(Draws[, anyDuplicated(Dn) > 0L, by = ky]$V1))
  expectReproduces(summarizeEnsemble(Draws, weights = c(AM88 = 1), levels = Levels)$Dn, Ensemble)
  expect_error(evaluateDisplacements(Intensities, eps = stats::rnorm(7L), ky = 0.1, Ts = 0.4, Mw = 7, models = "XX99",
                                     NFC = "D100", pulse = FALSE), "models must name")
})

test_that("tied pools keep full precision in the general summary and text precision in the compatibility route", {
  # Owner's decision of 2026-10-13: summarizeEnsemble merges ties exactly;
  # the pool of buildDnTables reads every value back from its 15-digit text
  # once the pool holds a tie, as newmark 1.7.1 did; the means stay exact.
  Draws <- data.table::data.table(IDn = "AM88", ky = 0.1, Dn = c(0, 0, 0.1 + 0.2))
  General <- summarizeEnsemble(Draws, weights = c(AM88 = 1), levels = c(0.95, "mean"))$Dn
  Legacy <- .poolDisplacements(Draws, weights = c(AM88 = 1), levels = c(0.95, "mean"),
                               estimator = .textWeightedQuantile)$Dn
  expect_false(identical(General[[1L]], Legacy[[1L]]))
  expect_equal(General[[1L]], Legacy[[1L]])
  expect_identical(General[[1L]], .weightedQuantile(Draws$Dn, rep(1, 3L), 0.95))
  expect_identical(Legacy[[1L]], .weightedQuantile(c(0, 0, 0.3), rep(1, 3L), 0.95))
  expect_identical(General[[2L]], Legacy[[2L]])
})

test_that("buildDnTables and assembleDnTables reproduce the synthetic project products", {
  X <- reference()$synthetic
  Site <- buildDnTables(rock = X$inputs$rock, MCE = X$inputs$MCE, shear = X$inputs$shear, modelID = X$plan$model,
                        siteID = "SUK", TR = X$plan$family, Vref = 760, Mw = X$plan$mw, NS = X$plan$ns,
                        weights = X$plan$weights, zTrunc = X$plan$truncation, workers = 1L)
  Global <- assembleDnTables(list(SUK = Site$DnTable), shear = X$inputs$shear, weights = X$plan$weights,
                             modelID = X$plan$model, TR = X$plan$family, Vref = 760)
  PGA <- buildPGATable(X$inputs$rock, X$inputs$MCE, X$inputs$shear, modelID = X$plan$model, siteID = "SUK",
                       TR = X$plan$family)
  expectProducts(Site, Global, PGA, X$products)
})

test_that("buildDnTables reproduces the product-pairs case and leaves its inputs unchanged", {
  X <- reference()$pairs
  Inputs <- data.table::copy(X$inputs)
  Weights <- c(AM88 = 1, BM19 = 1)
  Site <- buildDnTables(rock = X$inputs$rock, MCE = X$inputs$MCE, shear = X$inputs$shear, modelID = "TEST", siteID = "TST",
                        TR = 100, Vref = 760, Mw = 7, NS = 2L, weights = Weights, zTrunc = 3, workers = 1L)
  expect_identical(X$inputs, Inputs)
  Global <- assembleDnTables(list(TST = Site$DnTable), shear = X$inputs$shear, weights = Weights, modelID = "TEST",
                             TR = 100, Vref = 760)
  PGA <- buildPGATable(X$inputs$rock, X$inputs$MCE, X$inputs$shear, modelID = "TEST", siteID = "TST", TR = 100)
  expectProducts(Site, Global, PGA, X$products)
  expect_true(all(is.na(Site$DnTable[ID == "MCE", TR])))
  expect_identical(sort(unique(Site$DnTable$IDn)), c("AM88", "BM19", "ensemble", "flexible", "rigid"))
})

test_that("the compatibility route rejects inputs outside the product contract", {
  X <- reference()$pairs
  expect_error(buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, "TEST", "TST", 100, Vref = 700, Mw = 7, NS = 2L,
                             weights = c(AM88 = 1), zTrunc = 3, workers = 1L), "760")
  expect_error(buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, "TEST", "TST", 100, Vref = 760, Mw = 7, NS = 2L,
                             weights = c(AM88 = 1, XX99 = 1), zTrunc = 3, workers = 1L), "active models")
  expect_error(buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, "TEST", "TST", 200, Vref = 760, Mw = 7, NS = 2L,
                             weights = c(AM88 = 1), zTrunc = 3, workers = 1L), "return periods")
})
