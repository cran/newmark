# Oracle: the products of project AR-M2V4D (data/hazard/uhs/M2V4D/CFD/UHSRock.Rds,
# data/hazard/tables/MCETable.Rds, data/newmark/ShearTable.Rds and
# data/newmark/M2V4D/CFD/PGATable.csv, SHA-256 a8b697c529af8fdc...),
# restricted to the shear pairs S135/U1 and
# S135/U2, and the newmark 1.7.1 site response of the first pair
# (.fitSiteProducts at TR 100 and .fitSaFConverged for the MCE), stored in
# fixtures/siteResponse-reference.rds by the installed newmark 1.7.1 (build
# d57b64f) in its own R process, never together with newmark 5.
projectReference <- function() readRDS(test_path("fixtures", "siteResponse-reference.rds"))$project

test_that("buildPGATable rewrites the stored PGATable.csv rows byte for byte", {
  Project <- projectReference()
  Table <- buildPGATable(rock = Project$rock, MCE = Project$mce, shear = Project$shear,
                         modelID = "M2V4D", siteID = "CFD", TR = Project$family)
  FILE <- tempfile(fileext = ".csv")
  on.exit(unlink(FILE), add = TRUE)
  data.table::fwrite(Table, FILE, na = "NA")
  expectReproducesLines(readLines(FILE), Project$csv)
  expect_identical(data.table::key(Table), c("siteID", "ID", "TR", "IDg", "IDm"))
  expect_identical(names(Table), c("IDg", "IDm", "ID", "TR", "PGA", "Vs30", "siteID", "units.PGA"))
})

test_that("the refined site demands of one pair are those of newmark 1.7.1", {
  Project <- projectReference()
  VS30 <- Project$shear[IDm == "U1", VS30]
  Mixture <- refineSiteMixture(Project$rock[TR == 100, .(Tn, k, muLn, sigmaLn, w)], VS30 = VS30, Vref = 760,
                               model = "ST20", levels = c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95),
                               orders = c(16, 32, 64, 128, 256), tolerance = 2e-3)
  expectReproduces(as.list(Mixture$mixture), as.list(Project$mixtureTR100))
  expect_identical(Mixture$order, 64L)
  Rock <- Project$mce[, .(Tn, p, Sa)]
  Spectrum <- refineSiteSpectrum(Rock, VS30 = VS30, Vref = 760, model = "ST20", orders = c(32, 64, 128, 256),
                                 tolerance = 2e-3)
  expectReproduces(as.list(Spectrum$spectrum), as.list(Project$spectrumMCE[, .(Tn, p, Sa, SaF)]))
  expect_identical(Spectrum$order, attr(Project$spectrumMCE, "quadratureOrder"))
  expectReproducesChange(Spectrum$change, attr(Project$spectrumMCE, "convergence"))
})

test_that("inconsistent product inputs fail explicitly", {
  Project <- projectReference()
  expect_error(buildPGATable(Project$rock, Project$mce, Project$shear, "M2V4D", "CFD", c(100, 200)),
               "return periods of TR")
  expect_error(buildPGATable(Project$rock, Project$mce, Project$shear, "M2V4D", "XYZ", Project$family),
               "return periods of TR")
  expect_error(buildPGATable(Project$rock, Project$mce[, TR := 1], Project$shear, "M2V4D", "CFD", Project$family),
               "must not contain TR")
})
