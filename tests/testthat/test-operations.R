# The registry of operations that the command-line interface reads from the
# installed package, against the functions it names: the interface decodes
# INPUT.json with these argument kinds and refuses a registry whose schema
# it does not read.
REGISTRY <- jsonlite::read_json(system.file("api", "operations.json", package = "newmark"), simplifyVector = FALSE)

test_that("the registry has schema 1, the schema the command-line interface reads", {
  expect_identical(REGISTRY$schema, 1L)
})

test_that("every entry names an exported function and exactly its formal arguments, doctor's resources aside", {
  Exports <- getNamespaceExports("newmark")
  for (Entry in REGISTRY$operations) {
    expect_true(Entry[["function"]] %in% Exports, info = Entry$name)
    Formals <- formals(getExportedValue("newmark", Entry[["function"]]))
    # The command passes doctor its own resource directory; INPUT.json has no such argument.
    if (Entry$route == "meta") Formals <- Formals[names(Formals) != "resources"]
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
             "numberOrNamed", "rds", "path")
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

# The report of doctor: outside an installed product every required check
# fails with the same detail; in a resource directory laid out as the
# installer lays it out, the checks pass until a file of the receipt changes.
test_that("describeModule() reports schema, package, version, summary and the five required checks", {
  Version <- as.character(utils::packageVersion("newmark"))
  Names <- c("library", "exports", "packages", "tools", "receipt")
  Report <- describeModule()
  expect_identical(Report[c("schema", "package", "version", "summary")],
                   list(schema = 1L, package = "newmark", version = Version, summary = list(status = "fail")))
  expect_identical(vapply(Report$checks, `[[`, "", "name"), Names)
  for (Check in Report$checks) {
    expect_identical(Check[c("required", "status", "detail")],
                     list(required = TRUE, status = "fail", detail = "no installed resources: run the installed command"),
                     info = Check$name)
  }
  Prefix <- tempfile("newmarkPrefix")
  on.exit(unlink(Prefix, recursive = TRUE), add = TRUE)
  Resources <- file.path(Prefix, "libexec", "newmark")
  dir.create(file.path(Resources, "install"), recursive = TRUE)
  dir.create(file.path(Resources, "lib"))
  writeLines(c("Package: newmark", paste("Version:", Version)), file.path(Resources, "lib", "DESCRIPTION"))
  writeLines('list(exports = "describeModule", packages = "jsonlite", tools = character())',
             file.path(Resources, "install", "requirements.R"))
  Files <- c("libexec/newmark/lib/DESCRIPTION", "libexec/newmark/install/requirements.R")
  jsonlite::write_json(list(package = list(name = "newmark"), file = Files, md5 = unname(tools::md5sum(file.path(Prefix, Files)))),
                       file.path(Resources, "install.json"), auto_unbox = TRUE)
  Report <- describeModule(Resources)
  expect_identical(Report$summary, list(status = "ok"))
  expect_identical(vapply(Report$checks, `[[`, "", "name"), Names)
  expect_identical(vapply(Report$checks, `[[`, "", "status"), rep("ok", 5L))
  expect_identical(Report$checks[[5L]]$detail, "2 files keep their recorded MD5")
  writeLines(c("Package: newmark", "Version: 0.0.0"), file.path(Resources, "lib", "DESCRIPTION"))
  Report <- describeModule(Resources)
  expect_identical(Report$summary, list(status = "fail"))
  expect_identical(vapply(Report$checks, `[[`, "", "status"), c("fail", "ok", "ok", "ok", "fail"))
  expect_identical(Report$checks[[1L]]$detail, paste("newmark", Version, "loads; lib/DESCRIPTION declares 0.0.0"))
  expect_identical(Report$checks[[5L]]$detail, "changed or missing: libexec/newmark/lib/DESCRIPTION")
})
