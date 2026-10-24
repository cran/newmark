# newmark 1.7.1 recorded these fixtures on aarch64-apple-darwin23. Keep the
# original exact comparisons there; other platforms' libm arithmetic can
# differ by a few ulps, so only double values receive tolerance 1e-12.
isReferencePlatform <- function() identical(R.version$platform, "aarch64-apple-darwin23")

# Keep types, attributes, missing values and all non-double content exact.
# A data.table's self-reference is a process-local pointer, not stored data.
.fixtureStructure <- function(x) {
  Attributes <- attributes(x)
  if (inherits(x, "data.table")) Attributes[[".internal.selfref"]] <- NULL
  Attributes <- Attributes[sort(names(Attributes))]
  Values <- x
  if (is.double(x)) Values[is.finite(x)] <- 0
  if (is.list(x)) Values <- lapply(x, .fixtureStructure)
  list(type = typeof(x), attributes = Attributes, values = Values)
}

expectReproduces <- function(actual, expected) {
  Label <- deparse1(substitute(actual))
  Expected <- deparse1(substitute(expected))
  if (isReferencePlatform()) {
    return(expect_identical(actual, expected, label = Label, expected.label = Expected))
  }
  Structure <- .fixtureStructure(actual)
  Reference <- .fixtureStructure(expected)
  if (!identical(Structure, Reference)) {
    return(expect_identical(Structure, Reference, label = Label, expected.label = Expected))
  }
  expect_equal(actual, expected, tolerance = 1e-12, label = Label, expected.label = Expected)
}

# A convergence change is a relative difference of two close spectra, so off
# the reference platform its rounding is bounded in absolute terms.
expectReproducesChange <- function(actual, expected) {
  Label <- deparse1(substitute(actual))
  Expected <- deparse1(substitute(expected))
  if (isReferencePlatform() || !identical(.fixtureStructure(actual), .fixtureStructure(expected))) {
    return(expect_identical(actual, expected, label = Label, expected.label = Expected))
  }
  expect_lte(abs(actual - expected), 1e-12, label = paste0("|", Label, " - ", Expected, "|"))
}

# These three table oracles originally used all.equal with zero tolerance.
expectReproducesTable <- function(actual, expected) {
  if (isReferencePlatform()) {
    return(expect_true(isTRUE(all.equal(actual, expected, tolerance = 0, check.attributes = TRUE))))
  }
  expectReproduces(actual, expected)
}

expectReproducesLines <- function(actual, expected) {
  if (isReferencePlatform()) return(expect_identical(actual, expected))
  # Only complete, unquoted numeric CSV fields are numbers; quoted text,
  # identifiers containing digits, separators and line order remain exact.
  Pattern <- '"(?:[^"]|"")*"(*SKIP)(*F)|(?:^|,)\\K[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?(?=,|$)'
  Structure <- gsub(Pattern, "<number>", actual, perl = TRUE)
  Reference <- gsub(Pattern, "<number>", expected, perl = TRUE)
  if (!identical(Structure, Reference)) return(expect_identical(actual, expected))
  Values <- lapply(regmatches(actual, gregexpr(Pattern, actual, perl = TRUE)), as.numeric)
  Reference <- lapply(regmatches(expected, gregexpr(Pattern, expected, perl = TRUE)), as.numeric)
  expectReproduces(Values, Reference)
}
