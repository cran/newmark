writeContract <- function(root, contract) {
  writeLines(jsonlite::toJSON(contract, digits = I(17), auto_unbox = TRUE), file.path(root, "newmark.json"))
}

# Normalized, because processProject resolves its root in that spelling.
newProjectRoot <- function() {
  Root <- tempfile("newmarkProject")
  dir.create(Root, recursive = TRUE)
  normalizePath(Root, winslash = "/")
}

# The nine levels of the ShearTable step.
LEVELS <- c(.02, .05, .10, .16, "mean", .84, .90, .95, .98)

# A ShearTable contract in the range form, with the seed of the comparisons
# of this chain (20260915); a case that needs a key absent or malformed
# edits a copy.
SHEAR <- list(path = list(tables = "products"), Vref.gmdp = 760, Hs = c(10, 12),
              lo.min = 0.35, lo.max = 0.075, s.min = 2.5, s.max = 3.5,
              uscs = list(U1 = c("ML", "MH")), IDm = "U1", seed = 20260915)

# A synthetic shear table to the schema the validators of buildDnTables()
# check: one geometry and `materials` materials, each pair with a 0.84 row
# that must not reach the products.
shearTable <- function(materials = 1L) {
  data.table::data.table(IDg = "S10", IDm = rep(paste0("U", seq_len(materials)), each = 2L), Hs = 10,
                         VS30 = c(400, 330), Ts = c(0.3, 0.45), Vref = 760, level = c("mean", "0.84"))
}

# The periods of the synthetic spectra; the last covers 1.5 times the
# period of the mean rows of shearTable().
PERIODS <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1)

# The MCE quantiles and mean of one lognormal spectrum per period, for each
# of `sites`.
mceTable <- function(sites) {
  Levels <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
  OUT <- data.table::CJ(siteID = sites, Tn = PERIODS)[, list(p = c(sprintf("%.2f", Levels), "mean"),
                                                          Sa = 0.6 * exp(-Tn + 0.55 * c(qnorm(Levels), 0.55 / 2))),
                                                      by = c("siteID", "Tn")]
  OUT[, c("ID", "Vs30", "Vref") := list("MCE", 760, 760)][]
}

# The rock spectrum of one site at the return period `TR`: one lognormal
# per period whose mean is the rock Sa, scaled by `scale`.
rockTable <- function(siteID, TR, scale) {
  Sa <- 0.45 * scale * exp(-PERIODS)
  data.table::data.table(ID = "M1", siteID = siteID, TR = TR, Tn = PERIODS, Sa = Sa,
                         Vs30 = 760, Vref = 760, SaF = Sa, AF = NA_real_, k = "1", muLn = log(Sa) - 0.55^2 / 2,
                         sigmaLn = 0.55, w = 1, fitErr = 0)
}

# The site steps of a project whose shear and MCE tables are synthetic and
# whose other inputs are stubs: before the site steps calculate, only those
# two tables are read.
SITES <- list(path = list(tables = "products", data = "data", uhs = "uhs", calc = "calc.json"),
              Vref.gmdp = 760, ID.gmdp = "M1", TR.gmdp = c(500, 100, 500), NS = 30, Mw.gmdp = 7,
              subduction = FALSE, DaH.gmdp = c(0.5, 1),
              Dn_weights = list(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM17 = 1, BM19 = 1))

stubSiteProject <- function(contract = SITES) {
  Root <- newProjectRoot()
  for (Site in c("S1", "S2")) {
    dir.create(file.path(Root, "uhs", "M1", Site), recursive = TRUE)
    saveRDS(1, file.path(Root, "uhs", "M1", Site, "UHSRock.Rds"))
  }
  dir.create(file.path(Root, "uhs", "M1", "S1.site"))
  dir.create(file.path(Root, "data"))
  saveRDS(mceTable(c("S1", "S2")), file.path(Root, "data", "MCETable.Rds"))
  writeLines('{"job":{"calculation":{"truncation_level":3}}}', file.path(Root, "calc.json"))
  for (Site in c("S1", "S2")) {
    dir.create(file.path(Root, "products", "M1", Site), recursive = TRUE)
    saveRDS(1, file.path(Root, "products", "M1", Site, "DnDraws.Rds"))
    writeLines("TR", file.path(Root, "products", "M1", Site, "PGATable.csv"))
  }
  saveRDS(shearTable(), file.path(Root, "products", "ShearTable.Rds"))
  writeContract(Root, contract)
  Root
}

# A project that the DnTable step can run: synthetic rock spectra of
# `sites`, each site with its own scale, the MCE of `mceSites`, the shear
# table and a contract of one return period, five realizations and one
# model, in one process.
DN <- list(path = list(tables = "products", data = "data", uhs = "uhs", calc = "calc.json"),
           Vref.gmdp = 760, ID.gmdp = "M1", TR.gmdp = 100, NS = 5, Mw.gmdp = 7, subduction = FALSE,
           DaH.gmdp = c(0.5, 1.5), Dn_weights = list(AM88 = 1), DN_PARALLEL = FALSE)

dnProject <- function(sites = c("S1", "S2"), mceSites = sites, contract = DN) {
  Root <- newProjectRoot()
  for (Site in sites) {
    dir.create(file.path(Root, "uhs", "M1", Site), recursive = TRUE)
    saveRDS(rockTable(Site, TR = 100, scale = c(S1 = 1, S2 = 1.2, S3 = 1.4)[[Site]]),
            file.path(Root, "uhs", "M1", Site, "UHSRock.Rds"))
  }
  dir.create(file.path(Root, "data"))
  saveRDS(mceTable(mceSites), file.path(Root, "data", "MCETable.Rds"))
  writeLines('{"job":{"calculation":{"truncation_level":3}}}', file.path(Root, "calc.json"))
  dir.create(file.path(Root, "products"))
  saveRDS(shearTable(), file.path(Root, "products", "ShearTable.Rds"))
  writeContract(Root, contract)
  Root
}

planQuietly <- function(...) suppressMessages(processProject(..., dryRun = TRUE))

test_that("the ShearTable step equals buildShearTable called directly and creates path.tables", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  writeContract(Root, SHEAR)
  expect_false(dir.exists(file.path(Root, "products")))
  expect_message(Effects <- processProject(root = Root, steps = "ShearTable"),
                 "ShearTable: wrote ShearTable (18 rows, seed 20260915)", fixed = TRUE)

  set.seed(20260915L)
  Expected <- buildShearTable(buildGeometry(Hs = c(10, 12), lo = seq(0.35, 0.075, length.out = 2L),
                                            s = seq(2.5, 3.5, length.out = 2L)),
                              USCS = list(U1 = c("ML", "MH")), IDm = "U1", Vref = 760, levels = LEVELS)
  FILE <- file.path(Root, "products", "ShearTable.Rds")
  expect_equal(readRDS(FILE), Expected, tolerance = 0)
  expect_identical(Effects, data.table::data.table(step = "ShearTable", siteID = NA_character_,
                                                   product = "ShearTable", path = "products/ShearTable.Rds",
                                                   action = "writes", rows = 18L))

  # Reproducible from the contract: a second run replaces the file with the
  # same bytes.
  Bytes <- readBin(FILE, "raw", file.size(FILE))
  Effects <- suppressMessages(processProject(root = Root, steps = "ShearTable"))
  expect_identical(Effects$action, "replaces")
  expect_identical(readBin(FILE, "raw", file.size(FILE)), Bytes)
  expect_identical(list.files(file.path(Root, "products"), all.files = TRUE, no.. = TRUE), "ShearTable.Rds")
})

test_that("a geometryTable and pairs reach buildShearTable as the contract declares them", {
  Geometry <- '[{"IDg":"S10","Hs":10<D>,"b":8<D>,"s":2<D>},{"IDg":"S20","Hs":20<D>,"b":16<D>,"s":3<D>}]'
  shearFrom <- function(decimals) {
    Root <- newProjectRoot()
    on.exit(unlink(Root, recursive = TRUE), add = TRUE)
    writeLines(paste0('{"path":{"tables":"products"},"Vref.gmdp":760,"seed":20260915,',
                      '"uscs":{"U1":["ML","MH"],"U2":["SP","SM"]},"IDm":["U1","U2"],',
                      '"geometryTable":', gsub("<D>", decimals, Geometry, fixed = TRUE), ',',
                      '"pairs":[{"IDg":"S10","IDm":"U1"},{"IDg":"S20","IDm":"U2"}]}'),
               file.path(Root, "newmark.json"))
    suppressMessages(processProject(root = Root, steps = "ShearTable"))
    readRDS(file.path(Root, "products", "ShearTable.Rds"))
  }
  Integers <- shearFrom("")
  Doubles <- shearFrom(".0")
  expect_equal(Integers, Doubles, tolerance = 0)
  expect_identical(vapply(Integers[, c("Hs", "b", "s")], typeof, ""), c(Hs = "double", b = "double", s = "double"))

  set.seed(20260915L)
  Expected <- buildShearTable(data.frame(IDg = c("S10", "S20"), Hs = c(10, 20), b = c(8, 16), s = c(2, 3)),
                              USCS = list(U1 = c("ML", "MH"), U2 = c("SP", "SM")), IDm = c("U1", "U2"),
                              Vref = 760, levels = LEVELS, pairs = data.frame(IDg = c("S10", "S20"), IDm = c("U1", "U2")))
  expect_equal(Integers, Expected, tolerance = 0)
})

test_that("file and the contract paths resolve against root, absolute paths are kept", {
  Root <- newProjectRoot()
  Products <- tempfile("newmarkProducts")
  on.exit(unlink(c(Root, Products), recursive = TRUE), add = TRUE)
  Contract <- SHEAR
  Contract$path$tables <- Products
  dir.create(file.path(Root, "spec"))
  writeLines(jsonlite::toJSON(Contract, digits = I(17), auto_unbox = TRUE), file.path(Root, "spec", "other.json"))
  Effects <- planQuietly(file = file.path("spec", "other.json"), root = Root, steps = "ShearTable")
  expect_identical(Effects$path, file.path(Products, "ShearTable.Rds"))
  expect_false(dir.exists(Products))
})

test_that("a call without root stops before the contract is read", {
  expect_error(processProject(steps = "ShearTable", dryRun = TRUE), 'argument "root" is missing, with no default',
               fixed = TRUE)
})

test_that("an absent, invalid or non-object contract fails naming the file", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  PATH.contract <- file.path(Root, "newmark.json")
  expect_error(processProject(root = Root, steps = "ShearTable"), "project contract not found", fixed = TRUE)
  writeLines('{"path":', PATH.contract)
  expect_error(processProject(root = Root, steps = "ShearTable"), "is not valid JSON", fixed = TRUE)
  writeLines("[1,2,3]", PATH.contract)
  expect_error(processProject(root = Root, steps = "ShearTable"), "must be one JSON object", fixed = TRUE)
  writeLines('{"Vref.gmdp":760,"Vref.gmdp":761}', PATH.contract)
  expect_error(processProject(root = Root, steps = "ShearTable"), "uniquely named keys", fixed = TRUE)
})

test_that("keys the contract schema does not define are refused by name, nested ones too", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  for (Key in c("Hs.min", "Hs.max", "NTS", "KMAX_PARALLEL", "KMAX_WORKERS", "sed")) {
    writeContract(Root, c(SHEAR, setNames(list(1), Key)))
    expect_error(planQuietly(root = Root, steps = "ShearTable"),
                 paste0("declares keys that newmark does not read: ", Key), fixed = TRUE)
  }
  Contract <- SHEAR
  Contract$path$uhs.old <- "uhs"
  writeContract(Root, Contract)
  expect_error(planQuietly(root = Root, steps = "ShearTable"), "does not read: path.uhs.old", fixed = TRUE)
})

test_that("without seed the ShearTable step runs on the caller's stream and reports no seed", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  writeContract(Root, SHEAR[names(SHEAR) != "seed"])
  set.seed(4L)
  Stream <- get(".Random.seed", envir = globalenv())
  expect_identical(planQuietly(root = Root, steps = "ShearTable")$product, "ShearTable")
  expect_identical(get(".Random.seed", envir = globalenv()), Stream)
  expect_message(Effects <- processProject(root = Root, steps = "ShearTable"),
                 "ShearTable: wrote ShearTable (18 rows, seed none)", fixed = TRUE)
  expect_identical(Effects$rows, 18L)
  # Not reseeded: the table is the one buildShearTable() draws from the same stream.
  set.seed(4L)
  Expected <- buildShearTable(buildGeometry(Hs = c(10, 12), lo = seq(0.35, 0.075, length.out = 2L),
                                            s = seq(2.5, 3.5, length.out = 2L)),
                              USCS = list(U1 = c("ML", "MH")), IDm = "U1", Vref = 760, levels = LEVELS)
  expect_equal(readRDS(file.path(Root, "products", "ShearTable.Rds")), Expected, tolerance = 0)
})

test_that("the ShearTable step checks a seed that the contract defines, and the heights", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  for (Seed in list(-1, 1.5, "20260915", TRUE, 2147483648, c(1, 2))) {
    writeContract(Root, modifyList(SHEAR, list(seed = Seed)))
    expect_error(planQuietly(root = Root, steps = "ShearTable"),
                 "must define seed as a whole number between 0 and 2147483647", fixed = TRUE)
  }
  writeContract(Root, modifyList(SHEAR, list(seed = 0)))
  expect_identical(nrow(planQuietly(root = Root, steps = "ShearTable")), 1L)
  # No default heights: the range form needs Hs.
  writeContract(Root, SHEAR[names(SHEAR) != "Hs"])
  expect_error(planQuietly(root = Root, steps = "ShearTable"), "must define Hs", fixed = TRUE)
  writeContract(Root, modifyList(SHEAR, list(IDm = c("U1", "U9"))))
  expect_error(planQuietly(root = Root, steps = "ShearTable"),
               "names materials that uscs does not define: U9", fixed = TRUE)
  expect_identical(list.files(Root), "newmark.json")
})

test_that("a dry run rejects the ShearTable inputs that a run rejects, before seeding or writing", {
  # The geometryTable and pairs of the test above; the range form with a
  # first truncation ratio of 0.005, whose berm width rounds to 0 m.
  Shear <- list(path = list(tables = "products"), Vref.gmdp = 760, seed = 20260915,
                uscs = list(U1 = c("ML", "MH"), U2 = c("SP", "SM")), IDm = c("U1", "U2"),
                geometryTable = data.frame(IDg = c("S10", "S20"), Hs = c(10, 20), b = c(8, 16), s = c(2, 3)),
                pairs = data.frame(IDg = c("S10", "S20"), IDm = c("U1", "U2")))
  Cases <- list(
    list(replace(Shear, "geometryTable", list(transform(Shear$geometryTable, Hs = c(-1, 20)))),
         "geometry requires positive finite Hs"),
    list(replace(Shear, "geometryTable", list(transform(Shear$geometryTable, IDg = "S10"))),
         "geometry requires unique nonempty IDg values"),
    list(replace(Shear, "pairs", list(transform(Shear$pairs, IDm = c("U1", "U9")))),
         "pairs must contain unique known geometry/material identifiers"),
    list(replace(Shear, "pairs", list(data.frame(IDg = c("S10", "S10"), IDm = c("U1", "U1")))),
         "pairs must contain unique known geometry/material identifiers"),
    list(replace(Shear, "Vref.gmdp", list(-1)), "Vref must be one positive finite velocity"),
    list(modifyList(SHEAR, list(lo.min = 0.005)), "geometry requires positive finite b")
  )
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  for (Case in Cases) {
    writeContract(Root, Case[[1L]])
    set.seed(4L)
    Stream <- get(".Random.seed", envir = globalenv())
    Dry <- tryCatch(planQuietly(root = Root, steps = "ShearTable"), error = identity)
    expect_identical(get(".Random.seed", envir = globalenv()), Stream)
    expect_identical(list.files(Root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE), "newmark.json")
    Run <- tryCatch(suppressMessages(processProject(root = Root, steps = "ShearTable")), error = identity)
    expect_s3_class(Dry, "error")
    expect_identical(conditionMessage(Dry),
                     paste0(file.path(Root, "newmark.json"), " defines ShearTable inputs that buildShearTable() rejects: ",
                            Case[[2L]]))
    expect_identical(conditionMessage(Run), conditionMessage(Dry))
    expect_identical(get(".Random.seed", envir = globalenv()), Stream)
    expect_identical(list.files(Root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE), "newmark.json")
  }
})

test_that("an unknown USCS code of a selected material stops a dry run as a run, before seeding or writing", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  # The rule is sampleSoilLayers()'s, whose own condition stays as it was.
  Sampler <- tryCatch(sampleSoilLayers(Hs = 10, USCS = c("ML", "UNKNOWN"), h = 1), error = identity)
  expect_identical(conditionMessage(Sampler),
                   "USCS must contain USCS codes among SP, SM, SC, SW, GW, GP, GM, GC, MH, ML, CH, CL, OH, OL, PT")
  expect_identical(conditionCall(Sampler), quote(sampleSoilLayers(Hs = 10, USCS = c("ML", "UNKNOWN"), h = 1)))
  # U2 is not selected, so the run never samples its code.
  writeContract(Root, modifyList(SHEAR, list(uscs = list(U1 = c("ML", "UNKNOWN"), U2 = "XX"))))
  set.seed(4L)
  Stream <- get(".Random.seed", envir = globalenv())
  Dry <- tryCatch(planQuietly(root = Root, steps = "ShearTable"), error = identity)
  Run <- tryCatch(suppressMessages(processProject(root = Root, steps = "ShearTable")), error = identity)
  expect_identical(conditionMessage(Dry),
                   paste0(file.path(Root, "newmark.json"), " defines uscs.U1 with codes that sampleSoilLayers() rejects: ",
                          conditionMessage(Sampler)))
  expect_identical(conditionMessage(Run), conditionMessage(Dry))
  expect_identical(get(".Random.seed", envir = globalenv()), Stream)
  expect_identical(list.files(Root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE), "newmark.json")
  writeContract(Root, modifyList(SHEAR, list(uscs = list(U1 = c("ML", "MH"), U2 = "XX"))))
  expect_identical(nrow(planQuietly(root = Root, steps = "ShearTable")), 1L)
})

test_that("a contract without seed plans the site steps, which never read it", {
  Root <- stubSiteProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Effects <- planQuietly(root = Root, steps = c("kmaxTable", "DnTable"))
  expect_identical(Effects$step, c(rep("DnTable", 8L), rep("kmaxTable", 3L)))
  expect_identical(Effects$siteID, c(rep(c("S1", "S2"), each = 3L), NA, NA, "S1", "S2", NA))
  expect_identical(Effects$product, c(rep(c("DnDraws", "PGATable", "DnTable"), 2L), "DnPlotTable", "DnTable",
                                      "kmaxTable", "kmaxTable", "kmaxTable"))
  expect_identical(Effects$path[c(1L, 2L, 7L, 9L, 11L)],
                   c("products/M1/S1/DnDraws.Rds", "products/M1/S1/PGATable.csv", "products/DnPlotTable.Rds",
                     "products/M1/S1/kmaxTable.Rds", "products/kmaxTable.Rds"))
  expect_identical(Effects$action, c("replaces", "replaces", "writes", "replaces", "replaces", rep("writes", 6L)))
  expect_identical(Effects$rows, rep(NA_integer_, 11L))
})

test_that("Vref.gmdp must be 760 for the site steps and path.data is required by DnTable only", {
  Root <- stubSiteProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  writeContract(Root, modifyList(SITES, list(Vref.gmdp = 800)))
  for (Step in c("DnTable", "kmaxTable")) {
    expect_error(planQuietly(root = Root, steps = Step), "must define Vref.gmdp as 760 m/s", fixed = TRUE)
  }
  Contract <- SITES
  Contract$path$data <- NULL
  writeContract(Root, Contract)
  expect_identical(nrow(planQuietly(root = Root, steps = "kmaxTable")), 3L)
  expect_error(planQuietly(root = Root, steps = "DnTable"), "must define path.data", fixed = TRUE)
})

test_that("DN_WORKERS is a positive whole number and DN_PARALLEL a logical, or the plan stops", {
  Root <- stubSiteProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  for (Workers in list(2.7, "3", TRUE, 0)) {
    writeContract(Root, c(SITES, list(DN_WORKERS = Workers)))
    expect_error(planQuietly(root = Root, steps = "DnTable"),
                 "must define DN_WORKERS as a whole number between 1 and", fixed = TRUE)
  }
  for (Spelling in c("3", "3.0")) {
    writeLines(sub("}$", paste0(',"DN_WORKERS":', Spelling, "}"),
                   jsonlite::toJSON(SITES, digits = I(17), auto_unbox = TRUE)), file.path(Root, "newmark.json"))
    expect_identical(nrow(planQuietly(root = Root, steps = "DnTable")), 8L)
  }
  writeContract(Root, c(SITES, list(DN_PARALLEL = 1)))
  expect_error(planQuietly(root = Root, steps = "DnTable"), "must define DN_PARALLEL as true or false", fixed = TRUE)
})

test_that("a site key absent or malformed fails naming the contract and the key", {
  Root <- stubSiteProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Cases <- list(list(TR.gmdp = "100", "must define TR.gmdp as numbers"),
                list(NS = 2.5, "must define NS as a whole number"),
                list(ID.gmdp = "1M", "must define ID.gmdp as an identifier"),
                list(DaH.gmdp = c(1, 1), "must define DaH.gmdp as distinct positive percentages of Hs"),
                list(subduction = 0, "must define subduction as true or false"),
                list(Dn_weights = list(AM88 = 1, BM17 = 1), "Dn_weights and subduction in"))
  for (Case in Cases) {
    Contract <- SITES
    Contract[names(Case)[[1L]]] <- Case[1L]
    writeContract(Root, Contract)
    expect_error(planQuietly(root = Root, steps = c("DnTable", "kmaxTable")), Case[[2L]], fixed = TRUE)
  }
  writeContract(Root, SITES[names(SITES) != "Mw.gmdp"])
  expect_identical(nrow(planQuietly(root = Root, steps = "kmaxTable")), 3L)
  expect_error(planQuietly(root = Root, steps = "DnTable"), "must define Mw.gmdp", fixed = TRUE)
})

test_that("an input that no selected step produces names the file and its producer", {
  Root <- stubSiteProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  unlink(file.path(Root, "products", "M1", "S1", "PGATable.csv"))
  expect_error(planQuietly(root = Root, steps = "kmaxTable"), "PGATable.csv (written by the DnTable step)", fixed = TRUE)
  expect_identical(nrow(planQuietly(root = Root, steps = "kmaxTable", sites = "S2")), 2L)
  unlink(file.path(Root, "products", "M1", "S1", "DnDraws.Rds"))
  expect_error(planQuietly(root = Root, steps = "kmaxTable", sites = "S1"),
               "DnDraws.Rds (written by the DnTable step)", fixed = TRUE)
  expect_identical(nrow(planQuietly(root = Root, steps = c("DnTable", "kmaxTable"), sites = "S1")), 7L)
  unlink(file.path(Root, "products", "ShearTable.Rds"))
  expect_error(planQuietly(root = Root, steps = "DnTable"), "ShearTable.Rds (written by the ShearTable step)",
               fixed = TRUE)
  unlink(file.path(Root, "data", "MCETable.Rds"))
  writeContract(Root, c(SITES, list(seed = 1, uscs = list(U1 = "ML"), IDm = "U1", Hs = 10, lo.min = 0.3,
                                    lo.max = 0.3, s.min = 2, s.max = 2)))
  expect_error(planQuietly(root = Root, steps = c("ShearTable", "DnTable")),
               "MCETable.Rds (a hazard product of the project)", fixed = TRUE)
  unlink(file.path(Root, "calc.json"))
  expect_error(planQuietly(root = Root, steps = "DnTable"), "calc.json (the classical calculation of the project)",
               fixed = TRUE)
  unlink(file.path(Root, "uhs", "M1", c("S1", "S2"), "UHSRock.Rds"))
  expect_error(planQuietly(root = Root, steps = c("ShearTable", "kmaxTable")),
               "missing input UHSRock.Rds (a hazard product of the project)", fixed = TRUE)
})

test_that("a dry run reports every destination in publication order and changes nothing", {
  Root <- stubSiteProject(c(SITES, list(seed = 20260915, uscs = list(U1 = "ML"), IDm = "U1", Hs = 10,
                                        lo.min = 0.3, lo.max = 0.3, s.min = 2, s.max = 2)))
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Before <- list.files(Root, recursive = TRUE, all.files = TRUE)
  Directory <- getwd()
  set.seed(4L)
  Stream <- get(".Random.seed", envir = globalenv())
  Messages <- capture_messages(Effects <- processProject(root = Root, sites = "S2", dryRun = TRUE))

  expect_identical(Effects$path, c("products/ShearTable.Rds", "products/M1/S2/DnDraws.Rds",
                                   "products/M1/S2/PGATable.csv", "products/M1/S2/DnTable.Rds",
                                   "products/DnPlotTable.Rds", "products/DnTable.Rds",
                                   "products/M1/S2/kmaxTable.Rds", "products/kmaxTable.Rds"))
  expect_identical(Effects$action, c(rep("replaces", 3L), rep("writes", 5L)))
  expect_true(all(c(paste0("dry run: ", file.path(Root, "newmark.json"), "\n"), "steps: ShearTable, DnTable, kmaxTable\n",
                    "sites: S2\n", "ShearTable: replaces products/ShearTable.Rds\n",
                    "kmaxTable: writes products/kmaxTable.Rds\n") %in% Messages))
  expect_identical(getwd(), Directory)
  expect_identical(get(".Random.seed", envir = globalenv()), Stream)
  expect_identical(list.files(Root, recursive = TRUE, all.files = TRUE), Before)
  expect_error(planQuietly(root = Root, steps = "DnTable", sites = c("S2", "S9")),
               "sites names sites without UHSRock.Rds under", fixed = TRUE)
})

test_that("the step names of newmark 1.7.1 are refused before the contract is read", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  for (Steps in list("ts", c("DnTable", "kmax"))) {
    expect_error(processProject(root = Root, steps = Steps), "steps must select ShearTable, DnTable or kmaxTable",
                 fixed = TRUE)
  }
})

# A tree of three published files and the group that would replace them,
# for the failures of the publication itself.
publicationTree <- function() {
  Root <- newProjectRoot()
  Paths <- c("A.Rds", "B.Rds", "C.Rds")
  for (Path in Paths) saveRDS(paste("old", Path), file.path(Root, Path))
  list(root = Root,
       group = data.table::data.table(step = "DnTable", siteID = "S1", product = c("A", "B", "C"), path = Paths,
                                      file = file.path(Root, Paths)),
       objects = lapply(Paths, function(x) data.frame(new = rep(x, 2L))))
}

contents <- function(tree) vapply(tree$group$file, function(x) {
  if (dir.exists(x)) "directory" else as.character(unlist(readRDS(x))[[1L]])
}, "", USE.NAMES = FALSE)

# The group of one file published before the failing one.
firstGroup <- function(tree) {
  data.table::data.table(step = "ShearTable", siteID = NA_character_, product = "ShearTable",
                         path = "ShearTable.Rds", file = file.path(tree$root, "ShearTable.Rds"))
}

# A destination under `root` that cannot be replaced by a rename: a
# directory in its place.
blockProduct <- function(root, path) {
  unlink(file.path(root, path))
  dir.create(file.path(root, path), recursive = TRUE)
  writeLines("x", file.path(root, path, "keep"))
}

test_that("a failure at the first rename of a run reaches the caller unchanged and replaces nothing", {
  Tree <- publicationTree()
  on.exit(unlink(Tree$root, recursive = TRUE), add = TRUE)
  blockProduct(Tree$root, Tree$group$path[[1L]])
  Error <- tryCatch(.publishRun(function(ledger) .publishGroup(ledger, Tree$group, Tree$objects)),
                    error = identity)
  expect_false(inherits(Error, "newmarkPartialPublication"))
  expect_match(conditionMessage(Error), "could not replace A.Rds", fixed = TRUE)
  expect_identical(contents(Tree), c("directory", "old B.Rds", "old C.Rds"))
  expect_identical(sort(list.files(Tree$root, all.files = TRUE, no.. = TRUE)), c("A.Rds", "B.Rds", "C.Rds"))
})

test_that("a failure at an intermediate rename leaves the files before it replaced and says so", {
  Tree <- publicationTree()
  on.exit(unlink(Tree$root, recursive = TRUE), add = TRUE)
  blockProduct(Tree$root, Tree$group$path[[2L]])
  Error <- tryCatch(.publishRun(function(ledger) {
    ledger$step <- "DnTable"
    ledger$siteID <- "S1"
    .publishGroup(ledger, Tree$group, Tree$objects)
  }), error = identity)
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_identical(Error$effects, data.table::data.table(step = "DnTable", siteID = "S1", product = "A",
                                                         path = "A.Rds", action = "replaces", rows = 2L))
  expect_identical(Error$path, "B.Rds")
  expect_match(conditionMessage(Error$cause), "could not replace B.Rds", fixed = TRUE)
  expect_match(conditionMessage(Error),
               "DnTable, site S1 failed at B.Rds after 1 published file(s): could not replace B.Rds", fixed = TRUE)
  expect_match(conditionMessage(Error), "\npublished: A.Rds\nkeep their previous content or remain absent: B.Rds, C.Rds",
               fixed = TRUE)
  expect_identical(contents(Tree), c("A.Rds", "directory", "old C.Rds"))
  expect_identical(sort(list.files(Tree$root, all.files = TRUE, no.. = TRUE)), c("A.Rds", "B.Rds", "C.Rds"))
})

test_that("a failure writing a temporary file changes no file of its group and leaves no temporary", {
  Tree <- publicationTree()
  on.exit(unlink(Tree$root, recursive = TRUE), add = TRUE)
  Group <- Tree$group
  Group$file[[3L]] <- file.path(Tree$root, "absent", "C.Rds")
  Group$path[[3L]] <- "absent/C.Rds"
  # The connection warning of saveRDS accompanies the induced failure.
  Error <- suppressWarnings(tryCatch(.publishRun(function(ledger) {
    ledger$step <- "ShearTable"
    .publishGroup(ledger, firstGroup(Tree), list(data.frame(x = 1:3)))
    ledger$step <- "DnTable"
    ledger$siteID <- "S1"
    .publishGroup(ledger, Group, Tree$objects)
  }), error = identity))
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_identical(Error$effects$path, "ShearTable.Rds")
  expect_identical(Error$path, "absent/C.Rds")
  # The cause is the error of saveRDS itself, as it signals it for any file
  # of that directory.
  Expected <- tryCatch(suppressWarnings(saveRDS(1, file.path(Tree$root, "absent", ".probe"))), error = identity)
  expect_identical(Error$cause, Expected)
  expect_match(conditionMessage(Error),
               paste0("DnTable, site S1 failed at absent/C.Rds after 1 published file(s): ", conditionMessage(Expected)),
               fixed = TRUE)
  expect_identical(contents(Tree), c("old A.Rds", "old B.Rds", "old C.Rds"))
  expect_identical(sort(list.files(Tree$root, all.files = TRUE, no.. = TRUE)),
                   c("A.Rds", "B.Rds", "C.Rds", "ShearTable.Rds"))
})

# .publishGroup with a writer that signals `condition`, resolved in an
# environment of the test; the package's own closure stays as it is.
failingPublisher <- function(condition) {
  Publisher <- .publishGroup
  environment(Publisher) <- list2env(list(saveRDS = function(object, file) stop(condition)),
                                     parent = environment(.publishGroup))
  Publisher
}

test_that("an error of the writer reaches the caller as the same condition, or as cause after a publication", {
  Tree <- publicationTree()
  on.exit(unlink(Tree$root, recursive = TRUE), add = TRUE)
  Cause <- errorCondition("induced write failure", class = "inducedWriteError", marker = "kept")
  publishFailing <- failingPublisher(Cause)
  Error <- tryCatch(.publishRun(function(ledger) publishFailing(ledger, Tree$group, Tree$objects)),
                    error = identity)
  expect_identical(Error, Cause)
  expect_s3_class(Error, "inducedWriteError")
  expect_identical(Error$marker, "kept")

  Error <- tryCatch(.publishRun(function(ledger) {
    ledger$step <- "ShearTable"
    .publishGroup(ledger, firstGroup(Tree), list(data.frame(x = 1:3)))
    ledger$step <- "DnTable"
    ledger$siteID <- "S1"
    publishFailing(ledger, Tree$group, Tree$objects)
  }), error = identity)
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_identical(Error$cause, Cause)
  expect_s3_class(Error$cause, "inducedWriteError")
  expect_identical(Error$cause$marker, "kept")
  expect_identical(Error$path, "A.Rds")
  expect_identical(Error$effects$path, "ShearTable.Rds")
  expect_identical(conditionMessage(Error),
                   paste0("DnTable, site S1 failed at A.Rds after 1 published file(s): induced write failure",
                          "\npublished: ShearTable.Rds\nkeep their previous content or remain absent: A.Rds, B.Rds, C.Rds"))
  expect_identical(contents(Tree), c("old A.Rds", "old B.Rds", "old C.Rds"))
  expect_identical(sort(list.files(Tree$root, all.files = TRUE, no.. = TRUE)),
                   c("A.Rds", "B.Rds", "C.Rds", "ShearTable.Rds"))
})

test_that("the only file of the ShearTable group is never half replaced", {
  Root <- newProjectRoot()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  writeContract(Root, SHEAR)
  dir.create(file.path(Root, "products", "ShearTable.Rds"), recursive = TRUE)
  writeLines("x", file.path(Root, "products", "ShearTable.Rds", "keep"))
  Error <- tryCatch(suppressMessages(processProject(root = Root, steps = "ShearTable")), error = identity)
  expect_false(inherits(Error, "newmarkPartialPublication"))
  expect_match(conditionMessage(Error), "could not replace products/ShearTable.Rds", fixed = TRUE)
  expect_identical(list.files(file.path(Root, "products"), all.files = TRUE, no.. = TRUE), "ShearTable.Rds")

  unlink(file.path(Root, "products"), recursive = TRUE)
  dir.create(file.path(Root, "products"))
  Sys.chmod(file.path(Root, "products"), "0555")
  on.exit(Sys.chmod(file.path(Root, "products"), "0755"), add = TRUE, after = FALSE)
  skip_if(file.access(file.path(Root, "products"), 2L) == 0L, "directory permissions are not enforced")
  Error <- tryCatch(suppressWarnings(suppressMessages(processProject(root = Root, steps = "ShearTable"))),
                    error = identity)
  # Nothing was published, so the error of saveRDS reaches the caller as
  # saveRDS signals it for any file of that directory.
  Expected <- tryCatch(suppressWarnings(saveRDS(1, file.path(Root, "products", ".probe"))), error = identity)
  expect_identical(Error, Expected)
  expect_identical(list.files(file.path(Root, "products"), all.files = TRUE, no.. = TRUE), character())
})

# The DnTable step on projects of synthetic inputs (dnProject()).

bytes <- function(path) readBin(path, "raw", file.size(path))

# The bytes that the step writes for an object: saveRDS(), or fwrite() with
# na = "NA" for the PGA table.
savedBytes <- function(x) {
  FILE <- tempfile()
  on.exit(unlink(FILE), add = TRUE)
  saveRDS(x, FILE)
  bytes(FILE)
}

csvBytes <- function(x) {
  FILE <- tempfile(fileext = ".csv")
  on.exit(unlink(FILE), add = TRUE)
  data.table::fwrite(x, FILE, na = "NA")
  bytes(FILE)
}

# The bytes of every DnTable product under a project, named by its path
# relative to the root.
dnProducts <- function(root) {
  FILES <- list.files(root, recursive = TRUE, pattern = "^(DnTable|DnDraws|DnPlotTable)[.]Rds$|^PGATable[.]csv$")
  stats::setNames(lapply(file.path(root, FILES), bytes), FILES)
}

# The files of the group of one site, in publication order.
siteGroup <- function(siteID) file.path("products", "M1", siteID, c("DnDraws.Rds", "PGATable.csv", "DnTable.Rds"))

# The temporary files that a publication leaves behind.
temporaries <- function(root) list.files(root, recursive = TRUE, all.files = TRUE, pattern = "^[.]")

test_that("the DnTable step publishes what buildDnTables, buildPGATable and assembleDnTables return", {
  Root <- dnProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Messages <- capture_messages(Effects <- processProject(root = Root, steps = "DnTable"))
  expect_true(all(c("DnTable: wrote site S1\n", "DnTable: wrote site S2\n") %in% Messages))

  Shear <- readRDS(file.path(Root, "products", "ShearTable.Rds"))
  MCE <- readRDS(file.path(Root, "data", "MCETable.Rds"))
  for (Site in c("S1", "S2")) {
    Rock <- readRDS(file.path(Root, "uhs", "M1", Site, "UHSRock.Rds"))
    Direct <- buildDnTables(Rock, MCE, Shear, modelID = "M1", siteID = Site, TR = 100, Vref = 760, Mw = 7, NS = 5,
                            weights = c(AM88 = 1), zTrunc = 3, workers = 1)
    PGA <- buildPGATable(Rock, MCE = MCE, shear = Shear, modelID = "M1", siteID = Site, TR = 100)
    FILES <- file.path(Root, siteGroup(Site))
    expect_identical(bytes(FILES[[1L]]), savedBytes(Direct$DnDraws))
    expect_identical(bytes(FILES[[2L]]), csvBytes(PGA))
    expect_identical(bytes(FILES[[3L]]), savedBytes(Direct$DnTable))
  }
  Tables <- lapply(c(S1 = "S1", S2 = "S2"), function(x) readRDS(file.path(Root, siteGroup(x)[[3L]])))
  Global <- assembleDnTables(Tables, shear = Shear, weights = c(AM88 = 1), modelID = "M1", TR = 100, Vref = 760)
  expect_identical(bytes(file.path(Root, "products", "DnPlotTable.Rds")), savedBytes(Global$DnPlotTable))
  expect_identical(bytes(file.path(Root, "products", "DnTable.Rds")), savedBytes(Global$DnTable))

  expect_identical(Effects, data.table::data.table(
    step = "DnTable", siteID = c(rep(c("S1", "S2"), each = 3L), NA, NA),
    product = c(rep(c("DnDraws", "PGATable", "DnTable"), 2L), "DnPlotTable", "DnTable"),
    path = c(siteGroup("S1"), siteGroup("S2"), "products/DnPlotTable.Rds", "products/DnTable.Rds"),
    action = "writes", rows = c(rep(c(300L, 2L, 1440L), 2L), nrow(Global$DnPlotTable), 2880L)))

  # Running the step again replaces every file with the same bytes.
  Products <- dnProducts(Root)
  expect_length(Products, 8L)
  Effects <- suppressMessages(processProject(root = Root, steps = "DnTable"))
  expect_identical(Effects$action, rep("replaces", 8L))
  expect_identical(dnProducts(Root), Products)
  expect_identical(temporaries(Root), character())
})

test_that("dn products repeat exactly and do not depend on the sites run alongside", {
  Alone <- dnProject(sites = "S1")
  Both <- dnProject()
  on.exit(unlink(c(Alone, Both), recursive = TRUE), add = TRUE)
  suppressMessages(processProject(root = Alone, steps = "DnTable"))
  Products <- dnProducts(Alone)
  expect_length(Products, 5L)

  suppressMessages(processProject(root = Both, steps = "DnTable"))
  expect_identical(dnProducts(Both)[siteGroup("S1")], Products[siteGroup("S1")])
})

test_that("the global tables hold every site found with a DnTable.Rds, selected or not", {
  Root <- dnProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Effects <- suppressMessages(processProject(root = Root, steps = "DnTable", sites = "S1"))
  expect_identical(Effects$path, c(siteGroup("S1"), "products/DnPlotTable.Rds", "products/DnTable.Rds"))
  # S2 is found but has no product yet, so the global tables hold S1 alone.
  expect_false(dir.exists(file.path(Root, "products", "M1", "S2")))
  for (FILE in c("DnPlotTable.Rds", "DnTable.Rds")) {
    expect_identical(unique(readRDS(file.path(Root, "products", FILE))$siteID), "S1")
  }
  Effects <- suppressMessages(processProject(root = Root, steps = "DnTable", sites = "S2"))
  expect_identical(Effects$action, c(rep("writes", 3L), rep("replaces", 2L)))
})

test_that("an invalid DnTable.Rds of another site stops the global rebuild after the selected site is published", {
  Root <- dnProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  suppressMessages(processProject(root = Root, steps = "DnTable"))
  saveRDS(data.table::data.table(x = 1), file.path(Root, siteGroup("S2")[[3L]]))
  Globals <- dnProducts(Root)[c("products/DnPlotTable.Rds", "products/DnTable.Rds")]
  Error <- tryCatch(suppressMessages(processProject(root = Root, steps = "DnTable", sites = "S1")),
                    error = identity)
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_identical(Error$effects$path, siteGroup("S1"))
  expect_identical(Error$effects$action, rep("replaces", 3L))
  expect_null(Error$path)
  expect_match(conditionMessage(Error), "DnTable failed after 3 published file(s): Invalid DnTable schema for siteID=S2",
               fixed = TRUE)
  expect_identical(dnProducts(Root)[names(Globals)], Globals)
})

test_that("a failure at a rename of the DnTable step leaves published only the files before it", {
  Cases <- list(
    list(block = siteGroup("S1")[[1L]], published = character(), kept = character()),
    list(block = siteGroup("S1")[[2L]], published = siteGroup("S1")[[1L]], kept = siteGroup("S1")[2:3]),
    list(block = siteGroup("S2")[[1L]], published = siteGroup("S1"), kept = siteGroup("S2")),
    list(block = "products/DnPlotTable.Rds", published = c(siteGroup("S1"), siteGroup("S2")),
         kept = c("products/DnPlotTable.Rds", "products/DnTable.Rds"))
  )
  Roots <- character()
  on.exit(unlink(Roots, recursive = TRUE), add = TRUE)
  for (Case in Cases) {
    Root <- dnProject()
    Roots <- c(Roots, Root)
    blockProduct(Root, Case$block)
    Error <- tryCatch(suppressMessages(processProject(root = Root, steps = "DnTable")), error = identity)
    expect_match(conditionMessage(Error), paste0("could not replace ", Case$block, ": "), fixed = TRUE)
    # Nothing published: the error of the rename reaches the caller as it is.
    expect_identical(inherits(Error, "newmarkPartialPublication"), length(Case$published) > 0L)
    expect_identical(as.character(Error$effects$path), Case$published)
    expect_identical(Error$path, if (length(Case$published)) Case$block)
    expect_identical(grepl(paste0("keep their previous content or remain absent: ", paste(Case$kept, collapse = ", ")),
                           conditionMessage(Error), fixed = TRUE), length(Case$published) > 0L)
    expect_identical(sort(names(dnProducts(Root))), sort(Case$published))
    expect_identical(temporaries(Root), character())
  }
})

test_that("a failure in the calculation of a later site leaves the sites before it published", {
  Root <- dnProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  saveRDS(data.table::data.table(x = 1), file.path(Root, "uhs", "M1", "S2", "UHSRock.Rds"))
  Error <- tryCatch(suppressMessages(processProject(root = Root, steps = "DnTable")), error = identity)
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_null(Error$path)
  expect_identical(conditionMessage(Error$cause), "UHSRock has an invalid schema")
  expect_identical(conditionMessage(Error),
                   paste0("DnTable, site S2 failed after 3 published file(s): UHSRock has an invalid schema",
                          "\npublished: ", paste(siteGroup("S1"), collapse = ", ")))
  expect_identical(sort(names(dnProducts(Root))), sort(siteGroup("S1")))
})

test_that("without ShearTable, the MCE demand and the heights are checked before the site steps calculate", {
  Root <- dnProject(mceSites = "S1")
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Before <- list.files(Root, recursive = TRUE, all.files = TRUE)
  for (DryRun in c(TRUE, FALSE)) {
    expect_error(suppressMessages(processProject(root = Root, steps = "DnTable", dryRun = DryRun)),
                 "MCETable has no MCE rock rows for siteID=S2, Vs30=Vref=760", fixed = TRUE)
  }
  expect_identical(list.files(Root, recursive = TRUE, all.files = TRUE), Before)
  expect_identical(nrow(planQuietly(root = Root, steps = "DnTable", sites = "S1")), 5L)
  # The MCE grid must reach 1.5 times the largest period of the mean rows.
  Shear <- shearTable()
  Shear$Ts[[1L]] <- 0.8
  saveRDS(Shear, file.path(Root, "products", "ShearTable.Rds"))
  expect_error(planQuietly(root = Root, steps = "DnTable", sites = "S1"),
               "MCE rock period grid does not cover required Dn period 1.5*max(Ts) = 1.2", fixed = TRUE)
  # kmaxTable derives its displacements from the heights, which DnTable does
  # not read.
  Shear <- shearTable()
  Shear$Hs <- NULL
  saveRDS(Shear, file.path(Root, "products", "ShearTable.Rds"))
  expect_identical(nrow(planQuietly(root = Root, steps = "DnTable", sites = "S1")), 5L)
  expect_error(planQuietly(root = Root, steps = c("DnTable", "kmaxTable"), sites = "S1"),
               paste0(file.path(Root, "products", "ShearTable.Rds"), " holds a shear table that buildKmaxTable() rejects: ",
                      "shear must carry finite positive Hs heights: the levels of DaH are derived from them"),
               fixed = TRUE)
  # With ShearTable selected, the shear table is the one the run would write.
  writeContract(Root, c(DN, SHEAR[c("Hs", "lo.min", "lo.max", "s.min", "s.max", "uscs", "IDm", "seed")]))
  expect_identical(nrow(planQuietly(root = Root)), 12L)
})

test_that("with ShearTable selected, the DnTable checks follow the publication of the shear table", {
  Root <- dnProject(mceSites = "S1", contract = c(DN, SHEAR[c("Hs", "lo.min", "lo.max", "s.min", "s.max", "uscs",
                                                               "IDm", "seed")]))
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  expect_identical(nrow(planQuietly(root = Root, steps = c("ShearTable", "DnTable"))), 9L)
  Error <- tryCatch(suppressMessages(processProject(root = Root, steps = c("ShearTable", "DnTable"))),
                    error = identity)
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_identical(Error$effects$path, "products/ShearTable.Rds")
  expect_identical(Error$effects$action, "replaces")
  expect_null(Error$path)
  expect_identical(conditionMessage(Error),
                   paste0("DnTable failed after 1 published file(s): MCETable has no MCE rock rows for siteID=S2, ",
                          "Vs30=Vref=760\npublished: products/ShearTable.Rds"))
  expect_identical(names(dnProducts(Root)), character())
})

test_that("DnTable runs on the cores that R detects minus one unless the contract decides, never on more than its pairs", {
  Root <- dnProject(sites = "S1")
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  saveRDS(shearTable(materials = 3L), file.path(Root, "products", "ShearTable.Rds"))
  Workers <- integer()
  local_mocked_bindings(buildDnTables = function(..., workers) {
    Workers <<- c(Workers, workers)
    stop("calculation stopped by the test")
  })
  Detected <- max(1L, parallel::detectCores(logical = TRUE) - 1L, na.rm = TRUE)
  Cases <- list(list(contract = DN[names(DN) != "DN_PARALLEL"], workers = min(Detected, 3L)),
                list(contract = modifyList(DN, list(DN_PARALLEL = TRUE, DN_WORKERS = 2)), workers = 2L),
                list(contract = modifyList(DN, list(DN_PARALLEL = TRUE, DN_WORKERS = 50)), workers = 3L),
                list(contract = modifyList(DN, list(DN_WORKERS = 2)), workers = 1L))
  for (Case in Cases) {
    writeContract(Root, Case$contract)
    expect_error(processProject(root = Root, steps = "DnTable"), "calculation stopped by the test", fixed = TRUE)
  }
  expect_identical(Workers, vapply(Cases, `[[`, integer(1L), "workers"))
})

# The kmaxTable step on projects of synthetic inputs (dnProject()), whose
# DnTable products a first run publishes.

# The bytes of every kmaxTable product under a project, named by its path
# relative to the root.
kmaxProducts <- function(root) {
  FILES <- list.files(root, recursive = TRUE, pattern = "^kmaxTable[.]Rds$")
  stats::setNames(lapply(file.path(root, FILES), bytes), FILES)
}

siteKmax <- function(siteID) file.path("products", "M1", siteID, "kmaxTable.Rds")

readPGA <- function(root, siteID) {
  data.table::fread(file.path(root, "products", "M1", siteID, "PGATable.csv"), colClasses = c(TR = "numeric"))
}

dnProjectRun <- function(...) {
  Root <- dnProject(...)
  suppressMessages(processProject(root = Root, steps = "DnTable"))
  Root
}

test_that("the kmaxTable step publishes what buildKmaxTable and assembleKmaxTables return", {
  Root <- dnProjectRun()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Dn <- dnProducts(Root)
  Messages <- capture_messages(Effects <- processProject(root = Root, steps = "kmaxTable"))
  expect_true(all(c("kmaxTable: wrote site S1\n", "kmaxTable: wrote site S2\n") %in% Messages))
  # The step leaves the DnTable products as they were.
  expect_identical(dnProducts(Root), Dn)

  Shear <- readRDS(file.path(Root, "products", "ShearTable.Rds"))
  Sites <- c(S1 = "S1", S2 = "S2")
  PGA <- lapply(Sites, readPGA, root = Root)
  Tables <- lapply(Sites, function(x) {
    buildKmaxTable(readRDS(file.path(Root, "products", "M1", x, "DnDraws.Rds")), PGATable = PGA[[x]], shear = Shear,
                   modelID = "M1", siteID = x, TR = 100, NS = 5, weights = c(AM88 = 1), DaH = c(0.5, 1.5))
  })
  for (Site in Sites) expect_identical(bytes(file.path(Root, siteKmax(Site))), savedBytes(Tables[[Site]]))
  Global <- assembleKmaxTables(lapply(Sites, function(x) readRDS(file.path(Root, siteKmax(x)))), shear = Shear,
                               modelID = "M1", TR = 100, weights = c(AM88 = 1), DaH = c(0.5, 1.5))
  expect_identical(bytes(file.path(Root, "products", "kmaxTable.Rds")), savedBytes(Global))
  # 0.5 and 1.5 percent of the 10 m geometry, in cm, and the MCE at TR = NA.
  expect_identical(sort(unique(Global$Da)), c(5, 15))
  expect_true(all(is.na(Global$TR[Global$ID == "MCE"])))
  expect_identical(Effects, data.table::data.table(
    step = "kmaxTable", siteID = c("S1", "S2", NA), product = "kmaxTable",
    path = c(siteKmax("S1"), siteKmax("S2"), "products/kmaxTable.Rds"), action = "writes", rows = c(96L, 96L, 192L)))

  # Running the step again replaces every file with the same bytes.
  Products <- kmaxProducts(Root)
  Effects <- suppressMessages(processProject(root = Root, steps = "kmaxTable"))
  expect_identical(Effects$action, rep("replaces", 3L))
  expect_identical(kmaxProducts(Root), Products)
  expect_identical(temporaries(Root), character())
})

test_that("the three steps run in order, and kmaxTable reads the shear table the run has just published", {
  Root <- dnProject(contract = c(DN, SHEAR[c("Hs", "lo.min", "lo.max", "s.min", "s.max", "uscs", "IDm", "seed")]))
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Effects <- suppressMessages(processProject(root = Root))
  expect_identical(Effects$step, c("ShearTable", rep("DnTable", 8L), rep("kmaxTable", 3L)))
  expect_identical(Effects$path[10:12], c(siteKmax("S1"), siteKmax("S2"), "products/kmaxTable.Rds"))
  # The heights of the new shear table, 10 and 12 m, not the 10 m of the
  # table it replaced.
  expect_identical(sort(unique(readRDS(file.path(Root, "products", "kmaxTable.Rds"))$Da)), c(5, 6, 15, 18))
  expect_identical(temporaries(Root), character())
})

test_that("the global kmaxTable holds every site found with a kmaxTable.Rds, selected or not", {
  Root <- dnProjectRun()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Effects <- suppressMessages(processProject(root = Root, steps = "kmaxTable", sites = "S1"))
  expect_identical(Effects$path, c(siteKmax("S1"), "products/kmaxTable.Rds"))
  # S2 is found but has no kmaxTable.Rds yet, so the global table holds S1 alone.
  expect_false(file.exists(file.path(Root, siteKmax("S2"))))
  expect_identical(unique(readRDS(file.path(Root, "products", "kmaxTable.Rds"))$siteID), "S1")
  Effects <- suppressMessages(processProject(root = Root, steps = "kmaxTable", sites = "S2"))
  expect_identical(Effects$action, c("writes", "replaces"))
})

test_that("a run over one site rebuilds the global tables from every site present", {
  Root <- dnProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  suppressMessages(processProject(root = Root, steps = c("DnTable", "kmaxTable")))
  Products <- c(dnProducts(Root), kmaxProducts(Root))
  Globals <- file.path(Root, "products", c("DnTable.Rds", "DnPlotTable.Rds", "kmaxTable.Rds"))
  unlink(Globals)
  suppressMessages(processProject(root = Root, steps = c("DnTable", "kmaxTable"), sites = "S1"))
  expect_identical(c(dnProducts(Root), kmaxProducts(Root))[names(Products)], Products)
  expect_identical(vapply(Globals, function(x) toString(sort(unique(readRDS(x)$siteID))), "", USE.NAMES = FALSE),
                   rep("S1, S2", 3L))
})

test_that("the global kmaxTable reads no PGATable.csv, which a selected site still needs", {
  Root <- dnProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  suppressMessages(processProject(root = Root, steps = c("DnTable", "kmaxTable")))
  Missing <- file.path(Root, "products", "M1", "S2", "PGATable.csv")
  unlink(Missing)
  Global <- bytes(file.path(Root, "products", "kmaxTable.Rds"))
  # S2 is not selected: its kmaxTable.Rds enters the global table without its PGATable.csv.
  Effects <- suppressMessages(processProject(root = Root, steps = "kmaxTable", sites = "S1"))
  expect_identical(Effects$path, c(siteKmax("S1"), "products/kmaxTable.Rds"))
  expect_identical(bytes(file.path(Root, "products", "kmaxTable.Rds")), Global)
  # A selected site without it stops the plan, before anything is written.
  Products <- kmaxProducts(Root)
  expect_error(suppressMessages(processProject(root = Root, steps = "kmaxTable")),
               paste0("missing input ", Missing, " (written by the DnTable step)"), fixed = TRUE)
  expect_identical(kmaxProducts(Root), Products)
})

test_that("a kmaxTable.Rds of another site that the global table rejects stops it after the selected site is published", {
  Root <- dnProject()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  suppressMessages(processProject(root = Root, steps = c("DnTable", "kmaxTable")))
  Global <- bytes(file.path(Root, "products", "kmaxTable.Rds"))
  # S2 left at the levels of another DaH.gmdp, then S2 not a kmaxTable at all.
  writeContract(Root, modifyList(DN, list(DaH.gmdp = c(0.5, 2))))
  Causes <- c("kmaxTable Da of siteID=S2 does not reproduce Da", "Invalid kmaxTable schema for siteID=S2")
  for (Cause in Causes) {
    if (Cause == Causes[[2L]]) saveRDS(data.table::data.table(x = 1), file.path(Root, siteKmax("S2")))
    Error <- tryCatch(suppressMessages(processProject(root = Root, steps = "kmaxTable", sites = "S1")), error = identity)
    expect_s3_class(Error, "newmarkPartialPublication")
    expect_identical(Error$effects$path, siteKmax("S1"))
    expect_identical(Error$effects$action, "replaces")
    expect_null(Error$path)
    expect_match(conditionMessage(Error), paste0("kmaxTable failed after 1 published file(s): ", Cause), fixed = TRUE)
    expect_identical(bytes(file.path(Root, "products", "kmaxTable.Rds")), Global)
  }
})

test_that("a failure at a rename of the kmaxTable step leaves published only the files before it", {
  Cases <- list(
    list(block = siteKmax("S1"), published = character(), site = "S1"),
    list(block = siteKmax("S2"), published = siteKmax("S1"), site = "S2"),
    list(block = "products/kmaxTable.Rds", published = c(siteKmax("S1"), siteKmax("S2")), site = NA)
  )
  Roots <- character()
  on.exit(unlink(Roots, recursive = TRUE), add = TRUE)
  for (Case in Cases) {
    Root <- dnProjectRun()
    Roots <- c(Roots, Root)
    blockProduct(Root, Case$block)
    Error <- tryCatch(suppressMessages(processProject(root = Root, steps = "kmaxTable")), error = identity)
    expect_match(conditionMessage(Error), paste0("could not replace ", Case$block, ": "), fixed = TRUE)
    # Nothing published: the error of the rename reaches the caller as it is.
    expect_identical(inherits(Error, "newmarkPartialPublication"), length(Case$published) > 0L)
    expect_identical(as.character(Error$effects$path), Case$published)
    if (length(Case$published)) {
      expect_identical(Error$path, Case$block)
      expect_match(conditionMessage(Error),
                   paste0("kmaxTable", if (!is.na(Case$site)) paste0(", site ", Case$site), " failed at ", Case$block,
                          " after ", length(Case$published), " published file(s)"), fixed = TRUE)
      expect_match(conditionMessage(Error), paste0("keep their previous content or remain absent: ", Case$block),
                   fixed = TRUE)
    }
    expect_identical(sort(names(kmaxProducts(Root))), sort(Case$published))
    expect_identical(temporaries(Root), character())
  }
})

test_that("a failure in the calculation of a later site of kmaxTable leaves the sites before it published", {
  Root <- dnProjectRun()
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  saveRDS(data.table::data.table(x = 1), file.path(Root, "products", "M1", "S2", "DnDraws.Rds"))
  Error <- tryCatch(suppressMessages(processProject(root = Root, steps = "kmaxTable")), error = identity)
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_null(Error$path)
  expect_match(conditionMessage(Error$cause), "Invalid DnDraws schema for siteID=S2", fixed = TRUE)
  expect_identical(conditionMessage(Error),
                   paste0("kmaxTable, site S2 failed after 1 published file(s): ", conditionMessage(Error$cause),
                          "\npublished: ", siteKmax("S1")))
  expect_identical(names(kmaxProducts(Root)), siteKmax("S1"))
})

# The message of the height rule for the shear table of a project.
heightsMessage <- function(root) {
  paste0(file.path(root, "products", "ShearTable.Rds"), " holds a shear table that buildKmaxTable() rejects: ",
         "shear must carry finite positive Hs heights: the levels of DaH are derived from them")
}

MISSINGMCE <- "MCETable has no MCE rock rows for siteID=S2, Vs30=Vref=760"

test_that("with both site steps failing, a dry run reports the MCE demand first and a run the heights, as 1.7.1", {
  Root <- dnProject(mceSites = "S1")
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  Shear <- shearTable()
  Shear$Hs <- NULL
  saveRDS(Shear, file.path(Root, "products", "ShearTable.Rds"))
  Before <- list.files(Root, recursive = TRUE, all.files = TRUE)
  Dry <- tryCatch(suppressMessages(processProject(root = Root, steps = c("DnTable", "kmaxTable"), dryRun = TRUE)),
                  error = identity)
  Run <- tryCatch(suppressMessages(processProject(root = Root, steps = c("DnTable", "kmaxTable"))), error = identity)
  expect_identical(conditionMessage(Dry), MISSINGMCE)
  expect_identical(conditionMessage(Run), heightsMessage(Root))
  # Each mode reports the MCE demand when only DnTable is selected.
  for (DryRun in c(TRUE, FALSE)) {
    expect_error(suppressMessages(processProject(root = Root, steps = "DnTable", dryRun = DryRun)), MISSINGMCE,
                 fixed = TRUE)
  }
  expect_identical(list.files(Root, recursive = TRUE, all.files = TRUE), Before)
})

test_that("a run checks the heights right after it publishes the shear table, before anything of DnTable", {
  Root <- dnProject(mceSites = "S1", contract = c(DN, SHEAR[c("Hs", "lo.min", "lo.max", "s.min", "s.max", "uscs", "IDm",
                                                               "seed")]))
  on.exit(unlink(Root, recursive = TRUE), add = TRUE)
  # A shear table without heights can only come from another builder.
  local_mocked_bindings(buildShearTable = function(...) {
    Shear <- shearTable()
    Shear$Hs <- NULL
    Shear
  })
  Error <- tryCatch(suppressMessages(processProject(root = Root)), error = identity)
  expect_s3_class(Error, "newmarkPartialPublication")
  expect_identical(Error$effects$path, "products/ShearTable.Rds")
  expect_null(Error$path)
  expect_identical(conditionMessage(Error), paste0("kmaxTable failed after 1 published file(s): ", heightsMessage(Root),
                                                   "\npublished: products/ShearTable.Rds"))
  expect_identical(names(dnProducts(Root)), character())
})
