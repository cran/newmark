# Oracle: the kmax fixture of test-kmax-table.R, the stored rows of project
# AR-M2V4D for geometry S27 and material U4 under TR 100 and the MCE, with
# its PGA table read back from the CSV and its shear table; buildKmaxTable()
# reproduces that site table exactly. Another site is the same table under
# another siteID, which the product checks treat as any site.
sliceInputs <- function() readRDS(test_path("fixtures", "kmax-M2V4D-S27-U4.Rds"))

# The arguments of the project's newmark.json: Dn_weights with subduction
# false (BM17 off) and DaH.gmdp.
WEIGHTS <- c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM19 = 1)
DAH <- c(0.05, 0.1, 0.25, 0.4, 0.5, 1, 2.5)
KEY <- c("siteID", "ID", "TR", "IDg", "IDm", "IDn", "p", "Da")

siteTable <- function(slice) {
  suppressWarnings(buildKmaxTable(slice$drawTable, PGATable = slice$PGATable, shear = slice$shearTable,
                                  modelID = "M2V4D", siteID = "CFD", TR = 100, NS = 300, weights = WEIGHTS, DaH = DAH))
}

relabel <- function(x, site) data.table::copy(x)[, siteID := site][]

# assembleKmaxTables() with the arguments of the fixture, any of them
# replaced through `...`; called by name, so a condition names the function.
assemble <- function(tables, slice, ...) {
  Args <- list(shear = slice$shearTable, modelID = "M2V4D", TR = 100, weights = WEIGHTS, DaH = DAH)
  Args[names(list(...))] <- list(...)
  do.call("assembleKmaxTables", c(list(tables = tables), Args))
}

test_that("one site gives the table of that site, byte for byte once saved", {
  Slice <- sliceInputs()
  Site <- siteTable(Slice)
  expect_warning(Global <- assemble(list(CFD = Site), Slice),
                 "global: 432/896 kmax rows use roots beyond sampled support", fixed = TRUE)
  expect_identical(serialize(Global, NULL), serialize(Site, NULL))
  expect_identical(data.table::key(Global), KEY)
})

test_that("the tables of several sites are stacked, keyed and reported under the label global", {
  Slice <- sliceInputs()
  Site <- siteTable(Slice)
  Tables <- list(ZZZ = relabel(Site, "ZZZ"), CFD = Site)
  expect_warning(Global <- assemble(Tables, Slice), "global: 864/1792 kmax rows use roots beyond sampled support",
                 fixed = TRUE)
  expect_identical(data.table::key(Global), KEY)
  # Keyed by siteID first, so CFD precedes ZZZ whatever the order of tables.
  Reversed <- suppressWarnings(assemble(rev(Tables), Slice))
  expect_identical(serialize(Reversed, NULL), serialize(Global, NULL))
  expect_identical(nrow(Global), 2L * nrow(Site))
  expect_identical(unique(Global$siteID), c("CFD", "ZZZ"))
  expect_identical(as.data.frame(Global[siteID == "CFD"]), as.data.frame(Site))
  expect_identical(as.data.frame(Global[siteID == "ZZZ"]), as.data.frame(Tables$ZZZ))
})

test_that("the caller's tables keep their rows, attributes and key after the assembly", {
  Slice <- sliceInputs()
  # Unkeyed and in another row order, so that keying the caller's own table
  # would show; the shear table must not gain the index of a subset.
  Tables <- list(CFD = siteTable(Slice)[order(-kmax)], ZZZ = as.data.frame(relabel(siteTable(Slice), "ZZZ")))
  Shear <- data.table::copy(Slice$shearTable)
  Before <- lapply(list(Tables, Shear), serialize, connection = NULL)
  Global <- suppressWarnings(assemble(Tables, Slice, shear = Shear))
  expect_identical(lapply(list(Tables, Shear), serialize, connection = NULL), Before)
  expect_null(data.table::key(Tables$CFD))
  expect_null(attr(Shear, "index"))
  expect_identical(data.table::key(Global), KEY)
})

test_that("the global table takes neither the number of realizations nor the PGA tables", {
  expect_identical(names(formals(assembleKmaxTables)),
                   c("tables", "shear", "modelID", "TR", "weights", "DaH"))
})

test_that("the levels are every percentage of DaH at every height, to six significant digits", {
  expect_identical(.deriveKmaxLevels(c(0.5, 1), data.frame(Hs = c(10, 20, 20))), c(5, 10, 20))
  # 0.07 * 100 is 7.0000000000000009 in floating point.
  expect_identical(.deriveKmaxLevels(0.07, data.frame(Hs = 100)), 7)
})

test_that("buildKmaxTable and assembleKmaxTables stop with the same rules of DaH and of the heights", {
  Slice <- sliceInputs()
  Cases <- list(list(DaH = c(0.1, 0.1), shear = Slice$shearTable,
                     message = "DaH must contain distinct positive percentages of Hs"),
                list(DaH = DAH, shear = Slice$shearTable[, !"Hs"],
                     message = "shear must carry finite positive Hs heights: the levels of DaH are derived from them"))
  for (Case in Cases) {
    Built <- tryCatch(buildKmaxTable(Slice$drawTable, PGATable = Slice$PGATable, shear = Case$shear, modelID = "M2V4D",
                                     siteID = "CFD", TR = 100, NS = 300, weights = WEIGHTS, DaH = Case$DaH),
                      error = identity)
    Assembled <- tryCatch(assemble(list(CFD = Slice$kmaxTable), Slice, shear = Case$shear, DaH = Case$DaH),
                          error = identity)
    expect_identical(conditionMessage(Built), Case$message)
    expect_identical(conditionMessage(Assembled), Case$message)
    expect_identical(conditionCall(Built)[[1L]], as.name("buildKmaxTable"))
    expect_identical(conditionCall(Assembled)[[1L]], as.name("assembleKmaxTables"))
  }
})

test_that("tables outside the product contract fail naming what they break", {
  Slice <- sliceInputs()
  Site <- siteTable(Slice)
  for (Tables in list(list(Site), list(CFD = Site, CFD = Site), list(`1CFD` = Site), list())) {
    expect_error(assemble(Tables, Slice), "tables must be a named list with unique site identifiers", fixed = TRUE)
  }
  # A table built at other levels, as a site left behind by a change of
  # DaH.gmdp.
  expect_error(assemble(list(CFD = Site), Slice, DaH = c(0.05, 0.1)),
               "kmaxTable Da of siteID=CFD does not reproduce Da", fixed = TRUE)
  expect_error(assemble(list(CFD = Site), Slice, TR = 200),
               "kmaxTable has invalid siteID, ID, TR, p, or IDn vocabulary", fixed = TRUE)
  expect_error(assemble(list(CFD = Site), Slice, weights = c(AM88 = 1, BT07 = 1)),
               "kmaxTable has invalid siteID, ID, TR, p, or IDn vocabulary", fixed = TRUE)
  expect_error(assemble(list(ZZZ = Site), Slice),
               "kmaxTable has invalid siteID, ID, TR, p, or IDn vocabulary", fixed = TRUE)
  expect_error(assemble(list(CFD = data.table::data.table(x = 1)), Slice),
               "Invalid kmaxTable schema for siteID=CFD", fixed = TRUE)
})
