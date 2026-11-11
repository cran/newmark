# Oracle: fixtures/dn-reference.rds, written by newmark 1.7.1 in its own
# process: the scenario of test-sampleSaCorr.R (seed 48), the dn products of
# the synthetic project of helper-syntheticProject.R (SHA-256
# 725ea79f4fb53681...) run by runProcess(steps = "dn"), with the amplified
# records of each demand, and fitDnProducts() and getDnTables() on the
# inputs of test-product-pairs.R (SHA-256 f1012729ba0324ec...) under the
# product seeds of that site.
reference <- function() readRDS(test_path("fixtures", "dn-reference.rds"))

# The object as it is once saved and read back, as the products are compared.
roundTrip <- function(x) unserialize(serialize(x, NULL))

csvLines <- function(x) {
  FILE <- tempfile(fileext = ".csv")
  on.exit(unlink(FILE), add = TRUE)
  data.table::fwrite(x, FILE, na = "NA")
  readLines(FILE)
}

# The dn products of one site; the PGATable is rebuilt by the site-response
# family's buildPGATable() from the same inputs.
expectProducts <- function(site, global, pga, expected) {
  Draws <- data.table::copy(expected$DnDraws)
  setnames(Draws, old = "s", new = "realization")
  expect_false("s" %in% names(site$DnDraws))
  expect_type(site$DnDraws$realization, "integer")
  expectReproduces(roundTrip(site$DnDraws), Draws)
  expectReproducesLines(csvLines(pga), expected$PGATable)
  expectReproduces(roundTrip(site$DnTable), expected$DnTable)
  expectReproduces(roundTrip(global$DnTable), expected$GlobalDnTable)
  expectReproduces(roundTrip(global$DnPlotTable), expected$DnPlotTable)
}

test_that("the scenario realizations and summaries reproduce the 1.7.1 scenario", {
  R <- reference()$scenario48
  Mixture <- data.table::data.table(Tn = c(0, 0.4, 1), k = "oq", sigmaLn = c(0.20, 0.45, 0.75), w = 1)
  Mixture[, muLn := log(c(0.80, 0.55, 0.42)) - sigmaLn^2 / 2]
  set.seed(48)
  Intensities <- sampleIntensities(Mixture, Ts = 0.4, NS = 7L, zTrunc = NA)
  Eps <- stats::rnorm(7L)
  Draws <- evaluateDisplacements(Intensities, eps = Eps, ky = c(0.01, 0.10), Ts = 0.4, Mw = 7,
                                 models = "AM88", NFC = "D100", pulse = FALSE)
  expectReproduces(Draws$Dn, R$draws$Dn)
  # The caller's realization identifiers are kept, paired with their rows.
  Relabeled <- data.table::copy(Intensities)[, realization := rev(realization) + 100L]
  Kept <- evaluateDisplacements(Relabeled, eps = Eps, ky = c(0.01, 0.10), Ts = 0.4, Mw = 7,
                                models = "AM88", NFC = "D100", pulse = FALSE)
  expect_identical(Kept$realization, rep(Relabeled$realization, 2L))
  expectReproduces(Kept$Dn, Draws$Dn)
  for (DT in list(data.table::copy(Intensities)[, realization := realization[[1L]]], data.table::copy(Intensities)[2L, realization := NA_integer_])) {
    expect_error(evaluateDisplacements(DT, eps = Eps, ky = 0.1, Ts = 0.4, Mw = 7, models = "AM88",
                                       NFC = "D100", pulse = FALSE), "non-missing unique realization")
  }
  expectReproduces(Intensities$PGA[Draws$realization], R$draws$PGA.s)
  Levels <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95, "mean")
  Curve <- summarizeDisplacements(Draws, levels = Levels)
  expect_identical(Curve$level, R$curve$p[R$curve$IDn == "AM88"])
  expectReproduces(Curve$Dn, R$curve$Dn[R$curve$IDn == "AM88"])
  # The ensemble of 1.7.1 is the pool of the compatibility route; the
  # general summary coincides with it here because no pool holds a tie.
  Ensemble <- R$curve$Dn[R$curve$IDn == "ensemble"]
  expect_false(any(Draws[, anyDuplicated(Dn) > 0L, by = ky]$V1))
  expectReproduces(summarizeEnsemble(Draws, weights = c(AM88 = 1), levels = Levels)$Dn, Ensemble)
  expect_error(evaluateDisplacements(Intensities, eps = stats::rnorm(7L), ky = 0.1, Ts = 0.4, Mw = 7, models = "XX99",
                                     NFC = "D100", pulse = FALSE), "models must name")
})

test_that("tied pools keep full precision in the general summary and text precision in the compatibility route", {
  # Owner's decision of 2026-10-13: summarizeEnsemble merges ties exactly;
  # the pool of buildDnTables reads every value back from its 15-digit text
  # once the pool holds a tie, as newmark 1.7.1 did; means use the unrounded pool.
  Draws <- data.table::data.table(IDn = "AM88", ky = 0.1, Dn = c(0, 0, 0.1 + 0.2))
  General <- summarizeEnsemble(Draws, weights = c(AM88 = 1), levels = c(0.95, "mean"))$Dn
  Legacy <- .poolDisplacements(Draws, weights = c(AM88 = 1), levels = c(0.95, "mean"),
                               estimator = .textWeightedQuantile)$Dn
  # Exact distinction required by the contractual 15-digit text rounding.
  expect_false(identical(General[[1L]], Legacy[[1L]]))
  expect_equal(General[[1L]], Legacy[[1L]])
  # Keep the rounding on the compatibility route: tolerance would also accept swapped routes.
  expect_identical(General[[1L]], .weightedQuantile(Draws$Dn, rep(1, 3L), 0.95))
  expect_identical(Legacy[[1L]], .weightedQuantile(c(0, 0, 0.3), rep(1, 3L), 0.95))
  expectReproduces(General[[2L]], Legacy[[2L]])
})

test_that("buildDnTables and assembleDnTables reproduce the synthetic project products", {
  X <- reference()$synthetic
  Site <- buildDnTables(rock = X$inputs$rock, MCE = X$inputs$MCE, shear = X$inputs$shear, modelID = X$plan$model,
                        siteID = "SUK", TR = X$plan$family, Vref = 760, Mw = X$plan$mw, NS = X$plan$ns,
                        weights = X$plan$weights, zTrunc = X$plan$truncation, workers = 1L, draws = fixedShearDraws(X$inputs$shear))
  Global <- assembleDnTables(list(SUK = Site$DnTable), shear = X$inputs$shear, weights = X$plan$weights,
                             modelID = X$plan$model, TR = X$plan$family, Vref = 760)
  PGA <- buildPGATable(X$inputs$rock, X$inputs$MCE, X$inputs$shear, modelID = X$plan$model, siteID = "SUK",
                       TR = X$plan$family)
  expectProducts(Site, Global, PGA, X$products)
})

test_that("buildDnTables reproduces the product-pairs case and leaves its inputs unchanged", {
  X <- reference()$pairs
  Inputs <- data.table::copy(X$inputs)
  Weights <- c(AM88 = 1, BM19 = 1)
  Site <- buildDnTables(rock = X$inputs$rock, MCE = X$inputs$MCE, shear = X$inputs$shear, modelID = "TEST", siteID = "TST",
                        TR = 100, Vref = 760, Mw = 7, NS = 2L, weights = Weights, zTrunc = 3, workers = 1L,
                        draws = fixedShearDraws(X$inputs$shear))
  expect_identical(X$inputs, Inputs)
  Global <- assembleDnTables(list(TST = Site$DnTable), shear = X$inputs$shear, weights = Weights, modelID = "TEST",
                             TR = 100, Vref = 760)
  PGA <- buildPGATable(X$inputs$rock, X$inputs$MCE, X$inputs$shear, modelID = "TEST", siteID = "TST", TR = 100)
  expectProducts(Site, Global, PGA, X$products)
  expect_true(all(is.na(Site$DnTable[ID == "MCE", TR])))
  expect_identical(sort(unique(Site$DnTable$IDn)), c("AM88", "BM19", "ensemble", "flexible", "rigid"))
})

test_that("the compatibility route rejects inputs outside the product contract", {
  X <- reference()$pairs
  expect_error(buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, "TEST", "TST", 100, Vref = 700, Mw = 7, NS = 2L,
                             weights = c(AM88 = 1), zTrunc = 3, workers = 1L, draws = fixedShearDraws(X$inputs$shear)), "760")
  expect_error(buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, "TEST", "TST", 100, Vref = 760, Mw = 7, NS = 2L,
                             weights = c(AM88 = 1, XX99 = 1), zTrunc = 3, workers = 1L, draws = fixedShearDraws(X$inputs$shear)), "active models")
  expect_error(buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, "TEST", "TST", 200, Vref = 760, Mw = 7, NS = 2L,
                             weights = c(AM88 = 1), zTrunc = 3, workers = 1L, draws = fixedShearDraws(X$inputs$shear)), "return periods")
})

test_that("a magnitude per demand gives each demand the rows of its magnitude given for every demand", {
  X <- reference()$pairs
  build <- function(Mw) {
    buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, modelID = "TEST", siteID = "TST", TR = 100, Vref = 760,
                  Mw = Mw, NS = 2L, weights = c(AM88 = 1, BM19 = 1), zTrunc = 3, workers = 1L,
                  draws = fixedShearDraws(X$inputs$shear))
  }
  Low <- build(6)
  High <- build(7)
  Both <- build(c(MCE = 7, `100` = 6))
  for (Table in c("DnDraws", "DnTable")) {
    # BM19 depends on the magnitude, so the two magnitudes give other rows.
    expect_false(isTRUE(all.equal(Low[[Table]][ID == "TEST"], High[[Table]][ID == "TEST"])), info = Table)
    expectReproduces(Both[[Table]][ID == "TEST"], Low[[Table]][ID == "TEST"])
    expectReproduces(Both[[Table]][ID == "MCE"], High[[Table]][ID == "MCE"])
  }
})

test_that("a magnitude per demand names each return period of TR, as it prints, and the MCE once", {
  X <- reference()$pairs
  # The magnitude is checked before the rock table, so TR need not match it.
  Cases <- list(c(`100` = 6, MCE = 7), c(`100` = 6, `500` = 7), c(`100` = 6, `500` = 7, `1000` = 7, MCE = 7),
                c(`100` = 6, `100` = 6, `500` = 7, MCE = 7), c(`1e+02` = 6, `500` = 7, MCE = 7), c(6, 7, 7),
                c(`100` = 6, `500` = 0, MCE = 7), c(`100` = 6, `500` = NA, MCE = 7), "7")
  for (Mw in Cases) {
    expect_error(buildDnTables(X$inputs$rock, X$inputs$MCE, X$inputs$shear, "TEST", "TST", TR = c(500, 100), Vref = 760,
                               Mw = Mw, NS = 2L, weights = c(AM88 = 1), zTrunc = 3, workers = 1L,
                               draws = fixedShearDraws(X$inputs$shear)),
                 "Mw must be one finite positive number or one per demand (each return period of TR and MCE)", fixed = TRUE)
  }
})

test_that("movements cycle through paired column properties without changing the shared ky grid", {
  X <- reference()$pairs
  Columns <- data.table::as.data.table(X$inputs$shear)[level == "mean",
    .(realization = 1:2, Ts = c(0.2, 0.3), VS30 = c(360, 560), Vref = 760), by = .(IDg, IDm)]
  Preserved <- data.table::copy(Columns)
  Sample <- .sampleIntensities
  Evaluate <- evaluateDisplacements
  Sampling <- list()
  Evaluation <- list()
  local_mocked_bindings(
    .sampleIntensities = function(spectrum, Ts, NS, zTrunc, generator) {
      OUT <- Sample(spectrum, Ts = Ts, NS = NS, zTrunc = zTrunc, generator = generator)
      Sampling[[length(Sampling) + 1L]] <<- list(Ts = Ts, mean = OUT$mean, spectrum = spectrum)
      OUT
    },
    evaluateDisplacements = function(intensities, eps, ky, Ts, ...) {
      Evaluation[[length(Evaluation) + 1L]] <<- list(Ts = Ts, movements = intensities$realization,
                                                    intensities = copy(intensities), ky = ky)
      Evaluate(intensities, eps = eps, ky = ky, Ts = Ts, ...)
    })
  Args <- list(rock = X$inputs$rock, MCE = X$inputs$MCE, shear = X$inputs$shear,
               modelID = "TEST", siteID = "TST", TR = 100, Vref = 760, Mw = 7, NS = 5L,
               weights = c(AM88 = 1, BM19 = 1), zTrunc = 3, workers = 1L, draws = Columns)
  Site <- do.call(buildDnTables, Args)
  expect_identical(Columns, Preserved)
  # Input copies: the periods reach the sampler as the columns carry them.
  expect_identical(vapply(Sampling, `[[`, numeric(1L), "Ts"), rep(c(0.2, 0.3), 4L))
  Rock <- data.table::as.data.table(X$inputs$rock)[TR == 100, .(Tn, k, muLn, sigmaLn, w)]
  for (j in 1:2) {
    Expected <- refineSiteMixture(Rock, VS30 = c(360, 560)[[j]], Vref = 760, model = "ST20", levels = DNLEVELS,
                                  orders = c(16L, 32L, 64L, 128L, 256L), tolerance = 2e-3)$mixture
    expectReproduces(Sampling[[j]]$spectrum, Expected)
  }
  expect_identical(vapply(Evaluation, `[[`, numeric(1L), "Ts"), rep(c(0.2, 0.3), 4L))
  expect_identical(lapply(Evaluation, `[[`, "movements"), rep(list(c(1L, 3L, 5L), c(2L, 4L)), 4L))
  for (i in seq.int(1L, length(Evaluation), by = 2L)) {
    expect_identical(Evaluation[[i]]$ky, Evaluation[[i + 1L]]$ky)
    DT <- rbindlist(lapply(Evaluation[c(i, i + 1L)], `[[`, "intensities"))
    expect_equal(vapply(DT[, .(PGA, Sa13, Sa15)], mean, numeric(1L)),
                  setNames((3 * Sampling[[i]]$mean + 2 * Sampling[[i + 1L]]$mean) / 5,
                           c("PGA", "Sa13", "Sa15")), tolerance = DNTOL)
  }
  expect_true(all(Site$DnDraws[, .(N = uniqueN(realization)), by = .(ID, TR, IDg, IDm, IDn, ky)]$N == 5L))
  Fixed <- do.call(buildDnTables, c(Args[names(Args) != "draws"], list(draws = fixedShearDraws(X$inputs$shear))))
  expect_identical(Site$DnTable[, !"Dn"], Fixed$DnTable[, !"Dn"])
  expect_false(isTRUE(all.equal(Site$DnTable$Dn, Fixed$DnTable$Dn)))
  if (.Platform$OS.type != "windows") {
    Parallel <- do.call(buildDnTables, modifyList(Args, list(workers = 2L)))
    expectReproduces(Parallel, Site)
  }
  Columns[, Ts := max(X$inputs$rock$Tn)]
  expect_error(do.call(buildDnTables, Args), "period grid does not cover")
})
