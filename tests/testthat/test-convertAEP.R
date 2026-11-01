# Direct checks of the hazard-side calculations: definitions that hold
# exactly.
test_that("convertAEP applies the Poisson relations", {
  AUX <- convertAEP(AEP = c(1 / 475, 1 / 2475), IT = c(50, 1))
  expect_identical(AUX$TR, 1 / c(1 / 475, 1 / 2475))
  expect_identical(AUX$POE, -expm1(-c(50, 1) * c(1 / 475, 1 / 2475)))
  expect_error(convertAEP(AEP = 0.01, IT = 0), "finite and > 0")
})
