# The classed conditions of the AOM family that the library signals; the
# renderer of the command draws them and an R session shows their text.

test_that(".inform signals a family message whose text is one line and whose fields are named", {
  Note <- tryCatch(.inform("moduleNote", "Dry run of ", "a.json", operation = "process", phase = "plan", path = "a.json"),
                   message = identity)
  expect_s3_class(Note, c("moduleNote", "message", "condition"), exact = TRUE)
  expect_identical(conditionMessage(Note), "Dry run of a.json\n")
  expect_identical(Note[c("operation", "phase", "path")], list(operation = "process", phase = "plan", path = "a.json"))
  expect_null(Note$unit)
  Progress <- tryCatch(.inform("moduleProgress", "step ", 1L, operation = "process", count = 1L, total = 2L),
                       message = identity)
  expect_s3_class(Progress, c("moduleProgress", "message", "condition"), exact = TRUE)
  expect_identical(Progress[c("count", "total")], list(count = 1L, total = 2L))
  expect_s3_class(tryCatch(.inform("moduleDetail", "a detail", operation = "process"), message = identity),
                  c("moduleDetail", "message", "condition"), exact = TRUE)
  # Outside the command, R's default handler shows the text alone.
  expect_identical(capture_messages(.inform("moduleNote", "a note", operation = "process")), "a note\n")
})

test_that(".warn signals a family warning that names no call", {
  Warning <- tryCatch(.warn("Tn = ", 12, " s exceeds", operation = "evaluateSiteFactor", unit = "S1"), warning = identity)
  expect_s3_class(Warning, c("moduleWarning", "warning", "condition"), exact = TRUE)
  expect_identical(conditionMessage(Warning), "Tn = 12 s exceeds")
  expect_null(conditionCall(Warning))
  expect_identical(Warning[c("operation", "unit")], list(operation = "evaluateSiteFactor", unit = "S1"))
})

test_that("a failure after a publication is a partial publication of the family and a module failure", {
  Error <- tryCatch(.publishRun(function(ledger) {
    ledger$effects <- data.table::data.table(step = "ShearTable", siteID = NA_character_, product = "ShearTable",
                                             path = "ShearTable.Rds", action = "writes", rows = 1L)
    stop("a later failure")
  }), error = identity)
  expect_identical(class(Error), c("newmarkPartialPublication", "modulePartialPublication", "moduleFailure", "error",
                                   "condition"))
  expect_identical(conditionMessage(Error$cause), "a later failure")
  expect_identical(Error$effects$path, "ShearTable.Rds")
})
