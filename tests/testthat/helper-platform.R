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

# Calculated doubles use testthat edition 3's default tolerance on every host.
expectReproduces <- function(actual, expected) {
  Label <- deparse1(substitute(actual))
  Expected <- deparse1(substitute(expected))
  Structure <- .fixtureStructure(actual)
  Reference <- .fixtureStructure(expected)
  if (!identical(Structure, Reference)) {
    return(expect_identical(Structure, Reference, label = Label, expected.label = Expected))
  }
  expect_equal(actual, expected, label = Label, expected.label = Expected)
}

expectReproducesLines <- function(actual, expected) {
  # PGATable's fifth field is its calculated PGA. Mask only that unquoted
  # number: headers, quoted text, separators, row order, TR and Vs30 stay exact.
  Pattern <- '^(?:(?:"(?:[^"]|"")*"|[^",]*),){4}\\K[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?(?=,|$)'
  Structure <- sub(Pattern, "<PGA>", actual, perl = TRUE)
  Reference <- sub(Pattern, "<PGA>", expected, perl = TRUE)
  if (!identical(Structure, Reference)) return(expect_identical(actual, expected))
  Values <- data.table::fread(text = paste(actual, collapse = "\n"), colClasses = c(PGA = "numeric"))$PGA
  Reference <- data.table::fread(text = paste(expected, collapse = "\n"), colClasses = c(PGA = "numeric"))$PGA
  expectReproduces(Values, Reference)
}
