# A contract of the ShearTable step alone, in the range form, with the seed
# of the comparisons of this chain (20260915).
SHEARONLY <- list(path = list(tables = "products"), Vref.gmdp = 760, Hs = c(10, 12), lo.min = 0.35,
                  lo.max = 0.075, s.min = 2.5, s.max = 3.5, uscs = list(U1 = c("ML", "MH")), IDm = "U1",
                  seed = 20260915)

# A product as its step wrote it, read directly.
readStored <- function(path) {
  if (grepl("[.]csv$", path)) data.table::fread(path, colClasses = c(TR = "numeric")) else readRDS(path)
}

test_that("the product files under path.tables are listed in publication order, read as their steps wrote them", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Listing <- listProjectTables("newmark.json", root = Root)
  Paths <- c("products/ShearTable.Rds", "products/M2V4D/CFD/DnDraws.Rds", "products/M2V4D/CFD/PGATable.csv",
             "products/M2V4D/CFD/kmaxTable.Rds", "products/kmaxTable.Rds")
  FILES <- file.path(Root, Paths)
  Tables <- lapply(FILES, readStored)
  # data.table() takes `key` as its own argument, hence the rename.
  Expected <- data.table::data.table(
    table = c("ShearTable", "DnDraws", "PGATable", "kmaxTable", "kmaxTable"),
    step = c("ShearTable", "DnTable", "DnTable", "kmaxTable", "kmaxTable"),
    siteID = c(NA, "CFD", "CFD", "CFD", NA), path = Paths, bytes = file.size(FILES),
    rows = vapply(Tables, nrow, integer(1L)), columns = lapply(Tables, names),
    KEY = lapply(Tables, function(x) as.character(data.table::key(x))), date = file.mtime(FILES))
  data.table::setnames(Expected, "KEY", "key")
  expect_identical(Listing, Expected)
})

test_that("only files at the destinations of processProject are listed", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Before <- listProjectTables("newmark.json", root = Root)
  # None of these is read: each would stop the listing as a file that holds
  # no table.
  writeLines("not a product", file.path(Root, "products", "notes.txt"))
  saveRDS(1, file.path(Root, "products", ".kmaxTable.Rds.a1b2"))
  for (Directory in c("M2V4D/CFD.old", "M2V4D/CFD/archive", "OTHER/CFD")) {
    dir.create(file.path(Root, "products", Directory), recursive = TRUE)
    saveRDS(1, file.path(Root, "products", Directory, "kmaxTable.Rds"))
  }
  dir.create(file.path(Root, "products", "M2V4D", "EMPTY"))
  dir.create(file.path(Root, "products", "DnTable.Rds"))
  expect_identical(listProjectTables("newmark.json", root = Root), Before)
})

test_that("without products the listing has its columns and no rows", {
  Root <- newProductRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  writeProductContract(Root, PRODUCTCONTRACT)
  Empty <- listProjectTables("newmark.json", root = Root)
  expect_identical(nrow(Empty), 0L)
  expect_identical(vapply(Empty, function(x) class(x)[[1L]], ""),
                   c(table = "character", step = "character", siteID = "character", path = "character",
                     bytes = "numeric", rows = "integer", columns = "list", key = "list", date = "POSIXct"))
  dir.create(file.path(Root, "products", "M2V4D", "CFD"), recursive = TRUE)
  expect_identical(listProjectTables("newmark.json", root = Root), Empty)
})

test_that("a contract without ID.gmdp has global products only", {
  Root <- productProject(contract = PRODUCTCONTRACT[c("path", "Vref.gmdp")])
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  expect_identical(listProjectTables("newmark.json", root = Root)$path, c("products/ShearTable.Rds", "products/kmaxTable.Rds"))
})

test_that("the listing names and counts a product as the effects record of processProject does", {
  Root <- newProductRoot()
  Products <- tempfile("newmarkProducts")
  on.exit(unlink(c(Root, Products), recursive = TRUE), add = TRUE)
  writeProductContract(Root, SHEARONLY)
  Effects <- suppressMessages(processProject(root = Root, steps = "ShearTable"))
  Listing <- listProjectTables("newmark.json", root = Root)
  expect_identical(Listing$table, Effects$product)
  for (Column in c("step", "siteID", "path", "rows")) expect_identical(Listing[[Column]], Effects[[Column]])
  # An absolute path.tables is spelled absolute by both.
  Contract <- SHEARONLY
  Contract$path$tables <- Products
  writeProductContract(Root, Contract)
  Effects <- suppressMessages(processProject(root = Root, steps = "ShearTable"))
  expect_identical(Effects$path, file.path(Products, "ShearTable.Rds"))
  expect_identical(listProjectTables("newmark.json", root = Root)$path, Effects$path)
})

test_that("the contract is read and checked as processProject reads it", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  # The message of the listing is the one the plan of processProject gives
  # for the same contract.
  planMessage <- function(contract) {
    writeProductContract(Root, contract)
    Listing <- tryCatch(listProjectTables("newmark.json", root = Root), error = conditionMessage)
    Plan <- tryCatch(suppressMessages(processProject(root = Root, steps = "kmaxTable", dryRun = TRUE)),
                     error = conditionMessage)
    expect_identical(Listing, Plan)
    Listing
  }
  PATH.contract <- file.path(Root, "newmark.json")
  expect_identical(planMessage(c(PRODUCTCONTRACT, list(sed = 1))),
                   paste0(PATH.contract, " declares keys that newmark does not read: sed"))
  expect_identical(planMessage(modifyList(PRODUCTCONTRACT, list(ID.gmdp = "M2/V4D"))),
                   paste0(PATH.contract, " must define ID.gmdp as an identifier: a letter followed by letters, digits, _ or -"))
  expect_identical(planMessage(PRODUCTCONTRACT[c("Vref.gmdp", "ID.gmdp")]), paste0(PATH.contract, " must define path.tables"))
  writeProductContract(Root, PRODUCTCONTRACT)
  expect_error(listProjectTables("newmark.json", root = file.path(Root, "absent")), "project root does not exist", fixed = TRUE)
  expect_error(listProjectTables("other.json", root = Root), "project contract not found", fixed = TRUE)
  expect_error(listProjectTables(c("newmark.json", "other.json"), root = Root), "file must be one non-empty path", fixed = TRUE)
})

test_that("a product file that cannot be read or holds no table stops the listing naming it", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  FILE <- file.path(Root, "products", "M2V4D", "CFD", "kmaxTable.Rds")
  saveRDS(list(rows = 1), FILE)
  expect_error(listProjectTables("newmark.json", root = Root), paste0(FILE, " does not hold a table"), fixed = TRUE)
  writeLines("not RDS", FILE)
  expect_error(listProjectTables("newmark.json", root = Root), paste0("could not read ", FILE, ": "), fixed = TRUE)
})

test_that("the listing writes nothing", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Tree <- productTree(Root)
  Temporary <- list.files(tempdir(), recursive = TRUE, all.files = TRUE, no.. = TRUE)
  listProjectTables("newmark.json", root = Root)
  expect_identical(productTree(Root), Tree)
  expect_identical(list.files(tempdir(), recursive = TRUE, all.files = TRUE, no.. = TRUE), Temporary)
})
