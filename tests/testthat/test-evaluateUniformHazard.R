# The uniform-hazard evaluation on the mean curves that newmark 1.7.1
# computed from the OpenQuake exports of project AR-SABP0 (model SABP0, site
# caldera); their provenance is in fixtures/hazard-SABP0/FIXTURES.md.
FIXTURES <- test_path("fixtures", "hazard-SABP0")

test_that("evaluateUniformHazard gives spectra that grow with the return period", {
  Curves <- readRDS(file.path(FIXTURES, "AEPTable-1.7.1.Rds"))
  TR <- c(100, 200, 500, 1000, 2000, 2500, 5000, 10000)
  UHS <- evaluateUniformHazard(Curves[p == "mean"], TR = TR, by = c("Tn", "p"))
  expect_named(UHS, c("Tn", "p", "TR", "Sa"))
  expect_identical(nrow(UHS), length(unique(Curves[p == "mean"]$Tn)) * length(TR))
  expect_true(all(UHS[, !is.unsorted(Sa), by = Tn]$V1))
  expect_error(evaluateUniformHazard(Curves[p == "mean"], TR = c(100, 1e9), by = c("Tn", "p")),
               "exceeds imported curve range")
})
