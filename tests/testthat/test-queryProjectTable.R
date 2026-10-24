# The rows that a direct filter of a stored table selects, as a data.frame
# with the stored types, for comparison with a query.
directRows <- function(table, rows, columns = names(table)) {
  OUT <- as.data.frame(table)[rows, columns, drop = FALSE]
  row.names(OUT) <- NULL
  OUT
}

# A query result after the JSON of the command line, numbers with 17
# significant digits and missing values as null, and back, each column
# returned to the type it left with.
roundTrip <- function(x) {
  Back <- jsonlite::fromJSON(jsonlite::toJSON(x, digits = I(17), na = "null"))
  as.data.frame(Map(function(back, stored) as.vector(back, typeof(stored)), Back, x), optional = TRUE)
}

query <- function(root, table, siteID = NULL, where = list(), columns = NULL) {
  queryProjectTable("newmark.json", root = root, table = table, siteID = siteID, where = where, columns = columns)
}

stored <- function(root, ...) {
  FILE <- file.path(root, "products", ...)
  if (grepl("[.]csv$", FILE)) data.table::fread(FILE, colClasses = c(TR = "numeric")) else readRDS(FILE)
}

# The extrapolated column of the rows of a DnTable or kmaxTable at `value`
# (their ky or kmax, g), computed from buildKyLimits() at the PGA of each
# cell in the PGA table `PGA`: a model is outside its own range, and a tree
# outside when any of its models stored for the same cell is.
expectedMarks <- function(table, value, PGA) {
  Models <- c("AM88", "JB07", "SR08", "BT07", "BM17", "BM19")
  Trees <- list(ensemble = Models, rigid = c("AM88", "JB07", "SR08"), flexible = c("BT07", "BM17", "BM19"))
  Cell <- paste(table$siteID, table$ID, table$TR, table$IDg, table$IDm)
  CellPGA <- PGA$PGA[match(Cell, paste(PGA$siteID, PGA$ID, PGA$TR, PGA$IDg, PGA$IDm))]
  Branch <- paste(Cell, table$IDn)
  OUT <- rep(NA_real_, nrow(table))
  for (Group in unique(Branch)) {
    ROWS <- which(Branch == Group)
    IDn <- table$IDn[[ROWS[[1L]]]]
    Members <- if (IDn %in% Models) IDn else intersect(Trees[[IDn]], table$IDn[Cell == Cell[[ROWS[[1L]]]]])
    Limits <- as.data.frame(buildKyLimits(PGA = CellPGA[[ROWS[[1L]]]]))
    Limits <- Limits[Limits$IDn %in% Members, ]
    OUT[ROWS] <- vapply(value[ROWS], function(v) as.numeric(any(v < Limits$kyMin | v > Limits$kyMax)), numeric(1L))
  }
  OUT
}

test_that("a query returns the rows whose where columns take any of the given values, as a direct filter", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Stored <- stored(Root, "kmaxTable.Rds")
  Query <- query(Root, "kmaxTable", where = list(IDg = "S27", IDm = "U4", ID = "MCE", IDn = "ensemble", p = "mean"),
                 columns = c("Da", "kmax", "Kh", "beyond"))
  ROWS <- which(Stored$IDg == "S27" & Stored$IDm == "U4" & Stored$ID == "MCE" & Stored$IDn == "ensemble" &
                  Stored$p == "mean")
  expect_length(ROWS, 7L)
  expect_s3_class(Query, "data.table")
  expect_identical(as.data.frame(Query), directRows(Stored, ROWS, c("Da", "kmax", "Kh", "beyond")))
  # The values of one column are alternatives; the columns combine.
  Query <- query(Root, "kmaxTable", where = list(IDn = c("AM88", "BT07"), p = c("0.05", "mean")), columns = names(Stored))
  ROWS <- which((Stored$IDn == "AM88" | Stored$IDn == "BT07") & (Stored$p == "0.05" | Stored$p == "mean"))
  expect_identical(as.data.frame(Query), directRows(Stored, ROWS))
})

test_that("where = list() and columns = NULL return every product as it is stored, kmaxTable with extrapolated last", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Cases <- list(list("ShearTable", NULL, "ShearTable.Rds"), list("kmaxTable", NULL, "kmaxTable.Rds"),
                list("DnDraws", "CFD", c("M2V4D", "CFD", "DnDraws.Rds")),
                list("PGATable", "CFD", c("M2V4D", "CFD", "PGATable.csv")),
                list("kmaxTable", "CFD", c("M2V4D", "CFD", "kmaxTable.Rds")))
  for (Case in Cases) {
    Stored <- do.call(stored, c(list(Root), as.list(Case[[3L]])))
    Query <- query(Root, Case[[1L]], siteID = Case[[2L]])
    expect_identical(names(Query), c(names(Stored), if (Case[[1L]] == "kmaxTable") "extrapolated"))
    expect_identical(as.data.frame(Query)[names(Stored)], directRows(Stored, seq_len(nrow(Stored))))
  }
})

test_that("text values are converted to the type of their column, NA to a missing number", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Stored <- stored(Root, "kmaxTable.Rds")
  # The MCE rows carry TR = NA.
  Query <- query(Root, "kmaxTable", where = list(TR = "NA"), columns = c("ID", "TR"))
  expect_identical(as.data.frame(Query), directRows(Stored, which(is.na(Stored$TR)), c("ID", "TR")))
  expect_identical(unique(Query$ID), "MCE")
  # A short decimal and its 17 significant digits name the same stored double.
  ROWS <- which(Stored$Da == 1.35)
  expect_length(ROWS, 128L)
  for (Text in c("1.35", "1.3500000000000001")) {
    expect_identical(as.data.frame(query(Root, "kmaxTable", where = list(TR = c("100", "NA"), Da = Text),
                                         columns = names(Stored))),
                     directRows(Stored, ROWS))
  }
  # The realization s of DnDraws is an integer column: "4.0" is 4.
  Draws <- stored(Root, "M2V4D", "CFD", "DnDraws.Rds")
  Query <- query(Root, "DnDraws", siteID = "CFD", where = list(s = c("3", "4.0"), IDn = "AM88", TR = "NA"),
                 columns = c("s", "ky", "Dn"))
  ROWS <- which((Draws$s == 3L | Draws$s == 4L) & Draws$IDn == "AM88" & is.na(Draws$TR))
  expect_length(ROWS, 60L)
  expect_identical(as.data.frame(Query), directRows(Draws, ROWS, c("s", "ky", "Dn")))
  expect_type(Query$s, "integer")
  # The PGA table is read back from its CSV with TR numeric.
  PGA <- stored(Root, "M2V4D", "CFD", "PGATable.csv")
  Query <- query(Root, "PGATable", siteID = "CFD", where = list(TR = "NA"))
  expect_identical(as.data.frame(Query), directRows(PGA, which(is.na(PGA$TR))))
  expect_identical(Query$ID, "MCE")
  # A text column compares the text itself: the stored levels.
  Shear <- stored(Root, "ShearTable.Rds")
  Query <- query(Root, "ShearTable", where = list(level = c("0.1", "mean")), columns = c("level", "Ts"))
  expect_identical(as.data.frame(Query), directRows(Shear, which(Shear$level == "0.1" | Shear$level == "mean"),
                                                    c("level", "Ts")))
  Query <- query(Root, "kmaxTable", where = list(p = "0.50", ID = "M2V4D"), columns = names(Stored))
  expect_identical(as.data.frame(Query), directRows(Stored, which(Stored$p == "0.50" & Stored$ID == "M2V4D")))
})

test_that("each stored number written with the 17 significant digits of the JSON selects its rows", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  # The 30 yield accelerations of the slice, written as the JSON of a query
  # writes them; as.numeric() reads some of these texts one unit in the last
  # place away where R's long double adds no precision, as on macOS for
  # aarch64, and the query would then report them as values never taken.
  Draws <- stored(Root, "M2V4D", "CFD", "DnDraws.Rds")
  Grid <- unique(Draws$ky)
  expect_length(Grid, 30L)
  Text <- strsplit(gsub("^\\[|\\]$", "", jsonlite::toJSON(Grid, digits = I(17))), ",", fixed = TRUE)[[1L]]
  expect_identical(jsonlite::fromJSON(paste0("[", paste(Text, collapse = ","), "]")), Grid)
  expect_identical(as.data.frame(query(Root, "DnDraws", siteID = "CFD", where = list(ky = Text))),
                   directRows(Draws, seq_len(nrow(Draws))))
  ROWS <- which(Draws$ky == Grid[[17L]])
  expect_identical(as.data.frame(query(Root, "DnDraws", siteID = "CFD", where = list(ky = Text[[17L]]))),
                   directRows(Draws, ROWS))
})

test_that("the values of a query come back from the 17-digit JSON of the command line unchanged", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Queries <- list(query(Root, "kmaxTable", where = list(TR = "NA", IDn = c("ensemble", "BT07"))),
                  query(Root, "DnDraws", siteID = "CFD", where = list(s = c("1", "300"), TR = c("NA", "100"))),
                  query(Root, "PGATable", siteID = "CFD"),
                  query(Root, "ShearTable"))
  for (Query in Queries) {
    expect_gt(nrow(Query), 0L)
    expect_identical(roundTrip(Query), as.data.frame(Query))
  }
})

test_that("a DnTable query marks the rows outside the ky range of their model or of a model of their tree", {
  Root <- dnProductProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Stored <- stored(Root, "TEST", "TST", "DnTable.Rds")
  PGA <- stored(Root, "TEST", "TST", "PGATable.csv")
  Expected <- expectedMarks(Stored, Stored$ky, PGA)
  # AM88 leaves its range above the PGA of the cell, and so do the trees
  # that contain it; BM19 and its tree flexible stay within theirs.
  Marked <- vapply(split(Expected, Stored$IDn), sum, numeric(1L))
  expect_gt(Marked[["AM88"]], 0)
  expect_identical(unname(Marked[c("ensemble", "rigid", "BM19", "flexible")]), c(Marked[["AM88"]], Marked[["AM88"]], 0, 0))
  Query <- query(Root, "DnTable", siteID = "TST")
  expect_identical(as.data.frame(Query), cbind(directRows(Stored, seq_len(nrow(Stored))), extrapolated = Expected))
  # The global DnTable takes the PGA of each cell from the PGATable of its site.
  Global <- stored(Root, "DnTable.Rds")
  expect_identical(as.data.frame(query(Root, "DnTable")),
                   cbind(directRows(Global, seq_len(nrow(Global))), extrapolated = expectedMarks(Global, Global$ky, PGA)))
  # Selected columns without extrapolated return the stored ones as before.
  expect_identical(as.data.frame(query(Root, "DnTable", siteID = "TST", columns = c("IDn", "ky", "Dn"))),
                   directRows(Stored, seq_len(nrow(Stored)), c("IDn", "ky", "Dn")))
})

test_that("a kmaxTable query marks the rows whose kmax other than 0 is outside the ky range of their model or tree", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Stored <- stored(Root, "kmaxTable.Rds")
  PGA <- stored(Root, "M2V4D", "CFD", "PGATable.csv")
  Expected <- expectedMarks(Stored, Stored$kmax, PGA) * (Stored$kmax > 0)
  # Both bounds are crossed, JB07 below 0.05 g and AM88 above the PGA of the
  # cell; a kmax of 0, the null demand, is never marked although it lies
  # below every range.
  expect_true(any(Expected[Stored$IDn == "JB07" & Stored$kmax > 0 & Stored$kmax < 0.05] == 1))
  expect_true(any(Expected[Stored$IDn == "AM88" & Stored$kmax > Stored$PGA] == 1))
  expect_true(any(Stored$kmax == 0))
  for (Site in list(NULL, "CFD")) {
    expect_identical(as.data.frame(query(Root, "kmaxTable", siteID = Site)),
                     cbind(directRows(Stored, seq_len(nrow(Stored))), extrapolated = Expected))
  }
  # where filters extrapolated as a numeric column, before the projection.
  ROWS <- which(Expected == 1 & Stored$IDn == "rigid")
  expect_gt(length(ROWS), 0L)
  expect_identical(as.data.frame(query(Root, "kmaxTable", where = list(extrapolated = "1", IDn = "rigid"),
                                       columns = c("p", "Da", "kmax"))),
                   directRows(Stored, ROWS, c("p", "Da", "kmax")))
  Query <- query(Root, "kmaxTable", where = list(extrapolated = c("0", "1")), columns = c("extrapolated", "kmax"))
  expect_identical(as.data.frame(Query), data.frame(extrapolated = Expected, kmax = Stored$kmax))
})

test_that("the AM88 upper bound admits the last digits of the grid top of a mean PGA above 0.80 g", {
  Root <- newProductRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Site <- file.path(Root, "products", "TEST", "TST")
  dir.create(Site, recursive = TRUE)
  # The grid of this PGA ends one unit in the last place above it; one grid
  # step further lies beyond the bound.
  PGA <- 2.721
  Ky <- buildKyGrid(PGA = PGA, n = 30L)
  expect_gt(Ky[[30L]], PGA)
  Ky <- c(Ky, exp(2 * log(Ky[[30L]]) - log(Ky[[29L]])))
  saveRDS(data.table::data.table(IDg = "G1", IDm = "U1", p = "mean", ID = "TEST", TR = 100, ky = Ky, Dn = 1, Ts = 0.3,
                                 Vs30 = 360, IDn = "AM88", siteID = "TST", units.Dn = "cm"),
          file.path(Site, "DnTable.Rds"))
  data.table::fwrite(data.table::data.table(IDg = "G1", IDm = "U1", ID = "TEST", TR = 100, PGA = PGA, Vs30 = 360,
                                            siteID = "TST", units.PGA = "g"),
                     file.path(Site, "PGATable.csv"), na = "NA")
  writeProductContract(Root, list(path = list(tables = "products"), Vref.gmdp = 760, ID.gmdp = "TEST"))
  expect_identical(query(Root, "DnTable", siteID = "TST")$extrapolated, c(rep(0, 30L), 1))
})

test_that("a DnTable query without the PGA of a cell stops naming the PGATable or the cells", {
  Root <- dnProductProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  FILE <- file.path(Root, "products", "TEST", "TST", "PGATable.csv")
  PGA <- stored(Root, "TEST", "TST", "PGATable.csv")
  data.table::fwrite(PGA[ID == "MCE"], FILE, na = "NA")
  expect_error(query(Root, "DnTable", siteID = "TST"),
               paste("products/TEST/TST/DnTable.Rds has cells without a PGA in the PGATable of their site:",
                     "siteID=TST, ID=TEST, TR=100, IDg=G1, IDm=U1; siteID=TST, ID=TEST, TR=100, IDg=G2, IDm=U2"),
               fixed = TRUE)
  unlink(FILE)
  for (Site in list("TST", NULL)) {
    expect_error(query(Root, "DnTable", siteID = Site), paste0("missing input ", FILE, " (written by the DnTable step)"),
                 fixed = TRUE)
  }
  writeProductContract(Root, list(path = list(tables = "products"), Vref.gmdp = 760))
  expect_error(query(Root, "DnTable"),
               paste0(file.path(Root, "newmark.json"),
                      " must define ID.gmdp, which locates the PGATable of the sites of products/DnTable.Rds"),
               fixed = TRUE)
})

test_that("a product saved as a data.frame is queried as a data.table", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Shear <- as.data.frame(stored(Root, "ShearTable.Rds"))
  saveRDS(Shear, file.path(Root, "products", "ShearTable.Rds"))
  Query <- query(Root, "ShearTable", where = list(level = "mean"))
  expect_s3_class(Query, "data.table")
  expect_identical(as.data.frame(Query), directRows(Shear, which(Shear$level == "mean")))
})

test_that("an unknown product, site, column or value stops the query naming it", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  expect_error(query(Root, "kmaxTables"),
               "unknown table kmaxTables: the tables are ShearTable, DnPlotTable, DnTable, kmaxTable, DnDraws, PGATable",
               fixed = TRUE)
  expect_error(query(Root, "DnDraws"), "DnDraws is a table of a site: siteID must name the site", fixed = TRUE)
  expect_error(query(Root, "ShearTable", siteID = "CFD"), "ShearTable is a global table: siteID must be NULL", fixed = TRUE)
  expect_error(query(Root, "kmaxTable", siteID = "XYZ"),
               paste0("unknown site XYZ: the sites under ", file.path(Root, "products", "M2V4D"), " are CFD"), fixed = TRUE)
  expect_error(query(Root, "kmaxTable", siteID = "CFD.old"), "siteID must be NULL for a global table or one site identifier",
               fixed = TRUE)
  # A product that its step has not published yet names its file and step.
  expect_error(query(Root, "DnTable"),
               paste0("missing input ", file.path(Root, "products", "DnTable.Rds"), " (written by the DnTable step)"),
               fixed = TRUE)
  expect_error(query(Root, "DnTable", siteID = "CFD"),
               paste0("missing input ", file.path(Root, "products", "M2V4D", "CFD", "DnTable.Rds"),
                      " (written by the DnTable step)"), fixed = TRUE)
  expect_error(query(Root, "kmaxTable", where = list(IDx = "S27")), "products/kmaxTable.Rds has no column IDx", fixed = TRUE)
  expect_error(query(Root, "kmaxTable", columns = c("Da", "kmaxx")), "products/kmaxTable.Rds has no column kmaxx",
               fixed = TRUE)
  expect_error(query(Root, "kmaxTable", where = list(IDg = c("S27", "S99"))),
               "where.IDg gives values that column IDg of products/kmaxTable.Rds never takes: S99", fixed = TRUE)
  expect_error(query(Root, "kmaxTable", where = list(TR = c("150", "NA"))),
               "where.TR gives values that column TR of products/kmaxTable.Rds never takes: 150", fixed = TRUE)
  expect_error(query(Root, "kmaxTable", where = list(p = "0.5")),
               "where.p gives values that column p of products/kmaxTable.Rds never takes: 0.5", fixed = TRUE)
  expect_error(query(Root, "ShearTable", where = list(level = "0.10")),
               "where.level gives values that column level of products/ShearTable.Rds never takes: 0.10", fixed = TRUE)
  expect_error(query(Root, "kmaxTable", where = list(TR = c("100", "abc", ""))),
               "where.TR gives text that is not a number for column TR of products/kmaxTable.Rds: abc, ", fixed = TRUE)
  # Numbers are spelled as in JSON.
  expect_error(query(Root, "kmaxTable", where = list(TR = c("0x64", " 100", "1.", "Inf", "NaN"))),
               "where.TR gives text that is not a number for column TR of products/kmaxTable.Rds: 0x64,  100, 1., Inf, NaN",
               fixed = TRUE)
  expect_error(query(Root, "DnDraws", siteID = "CFD", where = list(s = c("3.5", "1e10", "NaN"))),
               "where.s gives text that is not a whole number for column s of products/M2V4D/CFD/DnDraws.Rds: 3.5, 1e10, NaN",
               fixed = TRUE)
  Shear <- stored(Root, "ShearTable.Rds")
  Shear$flag <- TRUE
  saveRDS(Shear, file.path(Root, "products", "ShearTable.Rds"))
  expect_error(query(Root, "ShearTable", where = list(flag = "TRUE")),
               "column flag of products/ShearTable.Rds holds logical values, which where does not filter", fixed = TRUE)
})

test_that("a site product needs ID.gmdp in the contract, which a global one does not", {
  Root <- productProject(contract = PRODUCTCONTRACT[c("path", "Vref.gmdp")])
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  expect_error(query(Root, "kmaxTable", siteID = "CFD"),
               paste0(file.path(Root, "newmark.json"), " must define ID.gmdp, which locates the products of a site"),
               fixed = TRUE)
  expect_identical(nrow(query(Root, "kmaxTable")), 896L)
})

test_that("malformed arguments stop the query before the contract is read", {
  Root <- newProductRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  # The root holds no contract, so each message comes from the argument.
  for (Product in list(c("kmaxTable", "DnTable"), NA_character_, "", 1)) {
    expect_error(query(Root, Product), "table must be one table name", fixed = TRUE)
  }
  for (Where in list(NULL, list("S27"), list(IDg = "S27", IDg = "S10"), list(TR = 100), list(IDg = character()),
                     list(IDg = NA_character_), data.frame(IDg = "S27"))) {
    expect_error(query(Root, "kmaxTable", where = Where),
                 "where must be a list of non-empty character vectors named by distinct columns", fixed = TRUE)
  }
  for (Columns in list(character(), c("Da", "Da"), NA_character_, 1)) {
    expect_error(query(Root, "kmaxTable", columns = Columns), "columns must be NULL or distinct column names", fixed = TRUE)
  }
  expect_error(query(Root, "kmaxTable"), "project contract not found", fixed = TRUE)
})

test_that("a query writes nothing and leaves the product as it was stored", {
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Tree <- productTree(Root)
  Temporary <- list.files(tempdir(), recursive = TRUE, all.files = TRUE, no.. = TRUE)
  Before <- stored(Root, "kmaxTable.Rds")
  Query <- query(Root, "kmaxTable", where = list(TR = "NA"))
  query(Root, "DnDraws", siteID = "CFD", where = list(s = "1"))
  query(Root, "PGATable", siteID = "CFD", where = list(TR = "NA"))
  expect_identical(productTree(Root), Tree)
  expect_identical(list.files(tempdir(), recursive = TRUE, all.files = TRUE, no.. = TRUE), Temporary)
  expect_identical(stored(Root, "kmaxTable.Rds"), Before)
})
