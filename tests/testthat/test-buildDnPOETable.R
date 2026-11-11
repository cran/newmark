# AR-S2J2P oracle d5c8ccc99dda9ad828763324b93a6d2aa75deb4a5ff062be9e4b9ada6e33fcef:
# Doornpoort/MCE/S60/U3, 30000 stored draws, NS 200 and the contract's seven
# DaH.gmdp thresholds. The fixture includes the 210 independent oracle rows.
dnPOEInputs <- function() {
  Fixture <- readRDS(test_path("fixtures", "dnpoe-S2J2P-S60-U3.Rds"))
  setnames(Fixture$draws, old = "s", new = "realization")
  Fixture
}

test_that("pointwise exceedances reproduce the real oracle", {
  Fixture <- dnPOEInputs()
  Actual <- do.call(buildDnPOETable, Fixture[names(Fixture) != "reference"])
  expectReproduces(as.data.frame(Actual), as.data.frame(Fixture$reference))
  expect_equal(Actual[DaH == 1, POE][1:4], c(0.764, 0.740, 0.706, 0.673))
  expect_identical(unique(Actual$n), 1000L)
  expect_identical(data.table::key(Actual), DNPOEKEY)
  expect_true(all(is.na(Actual$TR)))
})

test_that("weights count realizations with strict exceedance and preserve zeros", {
  Args <- dnPOEInputs()
  Args$reference <- NULL
  Args$draws <- Args$draws[ky == min(ky) & realization <= 2 & IDn %in% c("AM88", "JB07")]
  Args$draws[, Dn := c(0, 60, 61, 120)]
  Args$NS <- 2
  Args$weights <- c(AM88 = 1, JB07 = 3)
  Args$DaH <- c(0.5, 1, 2)
  Actual <- do.call(buildDnPOETable, Args)
  expect_equal(Actual$POE, c(7 / 8, 6 / 8, 0))
  expect_identical(Actual$n, rep(4L, 3L))
  Args$draws[, Dn := 121]
  expect_equal(do.call(buildDnPOETable, Args)$POE, rep(1, 3L))
})

test_that("DnPOE uses only each geometry's exact mean height and preserves inputs", {
  Args <- dnPOEInputs()
  Args$reference <- NULL
  Args$shear[level == "mean", Hs := 60.1234567]
  Before <- lapply(Args[c("draws", "shear")], data.table::copy)
  Seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv)
  Actual <- do.call(buildDnPOETable, Args)
  expect_equal(Actual$Da, Actual$DaH * 60.1234567)
  for (Name in names(Before)) {
    expect_identical(as.list(Args[[Name]]), as.list(Before[[Name]]))
    expect_identical(data.table::key(Args[[Name]]), data.table::key(Before[[Name]]))
    expect_identical(data.table::indices(Args[[Name]]), data.table::indices(Before[[Name]]))
  }
  expect_identical(if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv), Seed)
})

test_that("incomplete, duplicated and invalid DnPOE inputs fail explicitly", {
  Args <- dnPOEInputs()
  Args$reference <- NULL
  Draws <- Args$draws
  Args$draws <- Draws[-1]
  expect_error(do.call(buildDnPOETable, Args), "NS realizations")
  Args$draws <- rbind(Draws, Draws[1])
  expect_error(do.call(buildDnPOETable, Args), "duplicate")
  Args$draws <- data.table::copy(Draws)
  Args$draws[1, Dn := NA_real_]
  expect_error(do.call(buildDnPOETable, Args), "invalid models")
  Args <- dnPOEInputs()
  Args$reference <- NULL
  Args$draws[, units.Dn := "m"]
  expect_error(do.call(buildDnPOETable, Args), "units")
  Args <- dnPOEInputs()
  Args$reference <- NULL
  Args$shear <- rbind(Args$shear, Args$shear[level == "mean"][, `:=`(IDm = "U2", Hs = 61)])
  expect_error(do.call(buildDnPOETable, Args), "one positive finite Hs")
  Args <- dnPOEInputs()
  Args$reference <- NULL
  setnames(Args$draws, old = "realization", new = "s")
  expect_error(do.call(buildDnPOETable, Args), "DnDraws is missing columns: realization", fixed = TRUE)
})

test_that("assembly validates stored sites and does not mutate caller tables", {
  Fixture <- dnPOEInputs()
  Site <- do.call(buildDnPOETable, Fixture[names(Fixture) != "reference"])
  Tables <- list(Doornpoort = Site, Other = data.table::copy(Site)[, siteID := "Other"])
  Before <- lapply(Tables, data.table::copy)
  Global <- assembleDnPOETables(Tables, shear = Fixture$shear, modelID = Fixture$modelID,
                                TR = Fixture$TR, weights = Fixture$weights, DaH = Fixture$DaH)
  expect_identical(nrow(Global), 2L * nrow(Site))
  expect_identical(unique(Global$siteID), names(Tables))
  expect_identical(lapply(Tables, as.list), lapply(Before, as.list))
  Tables$Other[, Da := Da + 1]
  expect_error(assembleDnPOETables(Tables, shear = Fixture$shear, modelID = Fixture$modelID,
                                  TR = Fixture$TR, weights = Fixture$weights, DaH = Fixture$DaH), "thresholds")
})

test_that("DnPOETable process, listing and query use the existing publication path without PGA", {
  Slice <- productSlice()
  Root <- productProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Contract <- c(PRODUCTCONTRACT, list(TR.gmdp = 100, NS = 300, subduction = FALSE,
                                     Dn_weights = list(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM17 = 1, BM19 = 1),
                                     DaH.gmdp = c(0.05, 0.1, 0.25, 0.4, 0.5, 1, 2.5)))
  Contract$path$uhs <- "uhs"
  Site <- file.path(Root, "uhs", "M2V4D", "CFD")
  dir.create(Site, recursive = TRUE)
  # Site discovery only checks the path, as in stubSiteProject().
  saveRDS(1, file.path(Site, "UHSRock.Rds"))
  writeProductContract(Root, Contract)
  unlink(file.path(Root, "products", "M2V4D", "CFD", "PGATable.csv"))
  Before <- productTree(Root)
  Plan <- suppressMessages(processProject(root = Root, steps = "DnPOETable", dryRun = TRUE))
  expect_identical(productTree(Root), Before)
  expect_identical(Plan$product, rep("DnPOETable", 2L))
  Effects <- suppressMessages(processProject(root = Root, steps = "DnPOETable"))
  expect_identical(Effects$path, c("products/M2V4D/CFD/DnPOETable.Rds", "products/DnPOETable.Rds"))
  Expected <- buildDnPOETable(Slice$drawTable, shear = Slice$shearTable, modelID = "M2V4D", siteID = "CFD",
                              TR = 100, NS = 300, weights = c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM19 = 1),
                              DaH = Contract$DaH.gmdp)
  Actual <- queryProjectTable("newmark.json", root = Root, table = "DnPOETable", siteID = NULL,
                              where = list(), columns = NULL)
  expectReproduces(as.data.frame(Actual), as.data.frame(Expected))
  expect_false("extrapolated" %in% names(Actual))
  Listing <- listProjectTables("newmark.json", root = Root)
  expect_identical(Listing[table == "DnPOETable", path], Effects$path)
  Query <- queryProjectTable("newmark.json", root = Root, table = "DnPOETable", siteID = "CFD",
                             where = list(TR = "NA"), columns = c("ID", "TR", "POE"))
  # Query copies the values already published, without another calculation.
  expect_identical(as.data.frame(Query), as.data.frame(Actual[ID == "MCE", .(ID, TR, POE)]))
  expect_equal(readRDS(file.path(Root, "products", "kmaxTable.Rds")), Slice$kmaxTable)
})
