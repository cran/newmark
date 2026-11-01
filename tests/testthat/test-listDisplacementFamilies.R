# The families that the displacement and kmax tables pool, against the models
# that selectDisplacementModels() switches on in either tectonic regime.
test_that("the families are ensemble, rigid and flexible", {
  expect_identical(names(listDisplacementFamilies()), c("ensemble", "rigid", "flexible"))
})

test_that("every member of a family is a model that selectDisplacementModels() switches on", {
  Models <- unique(unlist(listDisplacementFamilies(), use.names = FALSE))
  Switches <- stats::setNames(rep(1, length(Models)), Models)
  Active <- lapply(c(TRUE, FALSE), function(x) {
    names(suppressMessages(selectDisplacementModels(switches = Switches, subduction = x)))
  })
  expect_setequal(unlist(Active), Models)
})
