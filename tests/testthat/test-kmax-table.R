# Oracle: the stored products of project AR-M2V4D (data/newmark/M2V4D/CFD of
# 2026-09-23, kmaxTable SHA-256 521790a8e5def10a3e65c31337d66f1da0508e5c2fb0f5bece79ee81ac405ff6),
# written by newmark 1.7.1 fitKmaxProducts(). The fixture kmax-M2V4D-S27-U4.Rds
# holds its DnDraws, its PGATable read back from the CSV with TR numeric and
# its ShearTable for geometry S27 and material U4 under TR 100 and the MCE,
# with the stored kmaxTable rows of those cells at the levels of 27 m. The
# arguments are those of the project's newmark.json: ID.gmdp M2V4D, site CFD,
# NS 300, DaH.gmdp, and Dn_weights with subduction false (BM17 off).
sliceInputs <- function() readRDS(test_path("fixtures", "kmax-M2V4D-S27-U4.Rds"))

sliceTable <- function(slice) {
  buildKmaxTable(draws = slice$drawTable, PGATable = slice$PGATable, shear = slice$shearTable,
                 modelID = "M2V4D", siteID = "CFD", TR = 100, NS = 300,
                 weights = c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM19 = 1),
                 DaH = c(0.05, 0.1, 0.25, 0.4, 0.5, 1, 2.5))
}

test_that("buildKmaxTable reproduces the stored newmark 1.7.1 kmaxTable rows exactly", {
  Slice <- sliceInputs()
  expect_warning(Actual <- sliceTable(Slice), "kmax rows use roots beyond sampled support")
  Reference <- Slice$kmaxTable
  expectReproducesTable(Reference, Actual)
  expect_identical(names(Reference), names(Actual))
  expect_identical(data.table::key(Reference), data.table::key(Actual))
  expect_identical(lapply(Reference, typeof), lapply(Actual, typeof))
  expect_equal(Actual$Kh, Actual$kmax / Actual$PGA, tolerance = 0)
})

test_that("buildKmaxTable leaves the caller's tables unchanged", {
  Slice <- sliceInputs()
  Before <- lapply(Slice[c("drawTable", "PGATable", "shearTable")], data.table::copy)
  suppressWarnings(sliceTable(Slice))
  expect_identical(lapply(Slice[c("drawTable", "PGATable", "shearTable")], as.data.frame), lapply(Before, as.data.frame))
  expect_null(attr(Slice$shearTable, "index"))
})

test_that("inconsistent product inputs fail explicitly", {
  Slice <- sliceInputs()
  Args <- list(draws = Slice$drawTable, PGATable = Slice$PGATable, shear = Slice$shearTable,
               modelID = "M2V4D", siteID = "CFD", TR = 100, NS = 300,
               weights = c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM19 = 1), DaH = c(0.05, 0.1))
  expect_error(do.call(buildKmaxTable, modifyList(Args, list(NS = 200))), "NS")
  expect_error(do.call(buildKmaxTable, modifyList(Args, list(TR = 200))), "PGATable")
  expect_error(do.call(buildKmaxTable, modifyList(Args, list(weights = c(AM88 = 1, BT07 = 1)))), "IDn vocabulary")
  expect_error(do.call(buildKmaxTable, modifyList(Args, list(DaH = c(0.1, 0.1)))), "DaH")
})

# Scenario cases of test-kmax-mixture.R and test-kmax-quantiles.R of newmark
# 1.7.1 (lib/tests/testthat at commit f7d11bb): every realization of a model
# has Dn = root / ky, so its kmax at Da = 1 is the root.
scenarioDraws <- function(roots, samples) {
  rbindlist(lapply(names(roots), function(model) {
    OUT <- CJ(s = seq_len(samples), ky = c(0.05, 0.10, 0.20, 0.30, 0.40))
    OUT[, `:=`(Dn = roots[[model]] / ky, IDn = model)]
  }))
}

test_that("a scenario summarizes every model and every family with an active member", {
  OUT <- .summarizeKmaxScenario(scenarioDraws(c(AM88 = 0.10, BT07 = 0.30), 20L), weights = c(AM88 = 0.5, BT07 = 0.5),
                                Da = 1)
  expect_equal(OUT[IDn == "AM88" & level == "mean", kmax], 0.10)
  expect_equal(OUT[IDn == "BT07" & level == "mean", kmax], 0.30)
  expect_equal(OUT[IDn == "ensemble" & level == "mean", kmax], 0.20)
  expect_equal(OUT[IDn == "ensemble" & level == "0.95", kmax], 0.30)
  expect_setequal(unique(OUT$IDn), c("AM88", "BT07", "ensemble", "rigid", "flexible"))
  expect_true(OUT[, setequal(level, c("mean", "0.05", "0.10", "0.16", "0.50", "0.84", "0.90", "0.95")) && .N == 8L,
                  by = .(IDn, Da)][, all(V1)])
})

test_that("JB07 joins the ensemble and rigid families when it is active", {
  Draws <- scenarioDraws(c(AM88 = 0.10, JB07 = 100, BT07 = 0.30, SR08 = 0.20, BM17 = 0.40), 20L)
  OUT <- .summarizeKmaxScenario(Draws, weights = c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM17 = 1), Da = 1)
  expect_equal(OUT[IDn == "JB07" & level == "mean", kmax], 100)
  expect_equal(OUT[IDn == "ensemble" & level == "mean", kmax], 101 / 5)
  expect_equal(OUT[IDn == "rigid" & level == "mean", kmax], 100.3 / 3)
  expect_equal(OUT[IDn == "flexible" & level == "mean", kmax], 0.35)
})
