# Oracle: the seeded ShearTable of newmark 1.7.1 with its slope-period height
# correction, fixtures/ts-reference-ShearTable-height.Rds (SHA-256
# be042be5adac0da19839bad01c0e4aa9306b37086d017e1b58593edda57783c1), produced
# by loading the patched working tree of the newmark repository (lib/R/getTs.R
# SHA-256 cc1abac174ed1a52555ca887b3ebdc289c71243811206e99e287427ad987e3b9,
# every other lib/R and lib/data file at commit 29ec237) with pkgload in its own process
# and calling getGeometry() and getTs() as checkTs.R of the frozen
# oqt-to-newmark succession did: set.seed(20260915L), heights 10 and 120 m,
# truncation 0.35 and 0.075, slopes 2.5 and 3.5, U1 = ML MH, U2 = SP SM,
# Vref 760, nine levels. Compared with all.equal at tolerance 0 and identical
# names, key and types, except the root an and the period Ts: the surrogate is
# fitted to CylinderRoots solved at every node instead of the Mathematica
# tables of newmark 1.7.1, which moves 34 of the 36 roots by at most 0.0013 and
# 4 of the periods by at most 0.001 s; the bounds add half a unit of their
# rounding to those maxima.
referenceShearTable <- function() readRDS(test_path("fixtures", "ts-reference-ShearTable-height.Rds"))

expectReferenceRows <- function(Reference, Actual) {
  Columns <- setdiff(names(Reference), c("an", "Ts"))
  expect_true(isTRUE(all.equal(Reference[, ..Columns], Actual[, ..Columns], tolerance = 0, check.attributes = TRUE)))
  expect_lt(max(abs(Actual$an - Reference$an)), 0.00135)
  expect_lt(max(abs(Actual$Ts - Reference$Ts)), 0.0015)
}

referenceInputs <- function() {
  list(geometry = buildGeometry(Hs = c(10, 120), lo = c(0.35, 0.075), s = c(2.5, 3.5)),
       USCS = list(U1 = c("ML", "MH"), U2 = c("SP", "SM")), IDm = c("U1", "U2"), Vref = 760,
       levels = c(0.02, 0.05, 0.10, 0.16, "mean", 0.84, 0.90, 0.95, 0.98))
}

test_that("buildShearTable reproduces the seeded ShearTable of the corrected newmark 1.7.1", {
  Reference <- referenceShearTable()
  set.seed(20260915L)
  Actual <- do.call(buildShearTable, referenceInputs())
  expectReferenceRows(Reference, Actual)
  expect_identical(names(Reference), names(Actual))
  expect_identical(data.table::key(Reference), data.table::key(Actual))
  expect_identical(lapply(Reference, typeof), lapply(Actual, typeof))
})

test_that("a row outside the support of the root surrogate takes the exact root", {
  # lo is 8 / (2 * 10 * 2.5 + 8) = 0.138 for In and 30 / (2 * 10 * 1 + 30) = 0.6,
  # above the surrogate maximum 0.5, for Out.
  set.seed(20260915L)
  Shear <- buildShearTable(data.frame(IDg = c("In", "Out"), Hs = 10, b = c(8, 30), s = c(2.5, 1)),
                           USCS = list(U1 = c("ML", "MH")), IDm = "U1", Vref = 760, levels = "mean")
  Surrogate <- evaluateRoot(buildRootModel(mode = 1, method = "nlm"), mo = Shear$mo, lo = Shear$lo)
  Exact <- evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = Shear$mo, lo = Shear$lo)
  expect_identical(Surrogate$withinSupport, c(TRUE, FALSE))
  expect_identical(Shear$an, round(c(Surrogate$an[1L], Exact$an[2L]), 4))
})

test_that("inconsistent geometry or material inputs fail explicitly", {
  Geometry <- buildGeometry(Hs = c(10, 120), lo = c(0.35, 0.075), s = c(2.5, 3.5))
  expect_error(buildGeometry(Hs = c(10, 20), lo = 0.3, s = 2.5), "equally sized")
  expect_error(buildShearTable(Geometry, list(U1 = "ML"), "missing", 760, "mean"), "IDm must select")
  expect_error(buildShearTable(Geometry, list(U1 = "ML"), "U1", 760, c(0.5, NA)), "levels must contain")
})
