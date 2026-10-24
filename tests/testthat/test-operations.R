# The registry of operations that the command-line interface reads from the
# installed package, against the functions it names: the interface decodes
# INPUT.json with these argument kinds and refuses a registry whose schema
# it does not read.
REGISTRY <- jsonlite::read_json(system.file("api", "operations.json", package = "newmark"), simplifyVector = FALSE)

test_that("the registry has schema 1, the schema the command-line interface reads", {
  expect_identical(REGISTRY$schema, 1L)
})

test_that("every entry names an exported function and exactly its formal arguments", {
  Exports <- getNamespaceExports("newmark")
  for (Entry in REGISTRY$operations) {
    expect_true(Entry[["function"]] %in% Exports, info = Entry$name)
    Formals <- formals(getExportedValue("newmark", Entry[["function"]]))
    expect_identical(vapply(Entry$arguments, `[[`, "", "name"), as.character(names(Formals)), info = Entry$name)
    # An argument is required exactly when its formal has no default.
    expect_identical(vapply(Entry$arguments, function(x) isTRUE(x$required), logical(1L)),
                     unname(vapply(Formals, function(x) identical(x, quote(expr = )), logical(1L))), info = Entry$name)
  }
})

test_that("every entry uses the vocabulary that the command-line interface decodes", {
  Names <- vapply(REGISTRY$operations, `[[`, "", "name")
  expect_false(anyDuplicated(Names) > 0L)
  Kinds <- c("number", "integer", "text", "logical", "table", "tables", "named", "namedText", "matrix", "numberOrNull",
             "rds", "path")
  checkColumns <- function(columns, label) {
    for (Column in columns) {
      expect_true(Column$type %in% c("double", "integer", "text", "logical", "tables"), info = label)
      expect_true(is.logical(Column$null) && length(Column$null) == 1L, info = label)
      if (Column$type == "tables") checkColumns(Column$columns, paste(label, Column$name))
    }
  }
  for (Entry in REGISTRY$operations) {
    expect_true(Entry$route %in% c("general", "compatibility", "project", "meta"), info = Entry$name)
    expect_true(Entry$random %in% c("none", "stream", "products"), info = Entry$name)
    expect_true(is.logical(Entry$json) && length(Entry$json) == 1L, info = Entry$name)
    expect_true(is.character(Entry$family) && is.character(Entry$result), info = Entry$name)
    expect_identical(unlist(Entry$interfaces), c("R", "CLI"), info = Entry$name)
    for (Argument in Entry$arguments) {
      expect_true(Argument$kind %in% Kinds, info = paste(Entry$name, Argument$name))
      if (Argument$kind %in% c("table", "tables")) checkColumns(Argument$columns, paste(Entry$name, Argument$name))
    }
  }
  # Function operations and the project verb listProjectTables are named as
  # their functions; doctor, process and query are the only other names.
  Functions <- vapply(REGISTRY$operations, `[[`, "", "function")
  Other <- Names != Functions
  expect_identical(stats::setNames(Functions[Other], Names[Other]),
                   c(doctor = "describeModule", process = "processProject", query = "queryProjectTable"))
})
