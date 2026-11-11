# The realizations of the shear table: one geometry or two, the pool ML MH
# and the seed of the comparisons of this chain (20260915), with the nine
# levels of the ShearTable step.
LEVELS <- c(0.02, 0.05, 0.10, 0.16, "mean", 0.84, 0.90, 0.95, 0.98)

drawInputs <- function(Hs = 10, lo = 0.35, s = 2.5) {
  list(geometry = buildGeometry(Hs = Hs, lo = lo, s = s), USCS = list(U1 = c("ML", "MH")), IDm = "U1", Vref = 760,
       levels = LEVELS)
}

drawSeeded <- function(inputs, ...) {
  set.seed(20260915L)
  do.call(buildShearDraws, c(inputs, list(...)))
}

test_that("the shear table of buildShearDraws is the table of buildShearTable on the same stream", {
  skip_on_cran()
  Inputs <- drawInputs(Hs = c(10, 120), lo = c(0.35, 0.075), s = c(2.5, 3.5))
  Shear <- drawSeeded(Inputs)
  StreamDraws <- .Random.seed
  set.seed(20260915L)
  Table <- do.call(buildShearTable, Inputs)
  expectReproduces(Shear$ShearTable, Table)
  # The same random numbers, and no other: the stream ends where it ends.
  expect_identical(.Random.seed, StreamDraws)
  expect_identical(names(Shear), c("ShearDraws", "ShearTable"))
  expect_identical(names(Shear$ShearDraws),
                   c("Hs", "b", "s", "Vref", "IDm", "IDg", "POP", "Hw", "USCS", "Go", "mo", "Ts", "VSo", "VS30",
                     "Z500", "Z1000", "realization", "lo", "an"))
  expect_identical(typeof(Shear$ShearDraws$realization), "integer")
  expect_identical(Shear$ShearDraws[, .N, by = .(IDg, IDm)]$N, c(100L, 100L))
  expect_identical(Shear$ShearDraws$realization, rep(1:100, 2L))
  expect_identical(unique(Shear$ShearDraws[, .(IDg, IDm)]), unique(Table[, .(IDg, IDm)]))
})

test_that("the realizations summarized as buildShearTable summarizes them reproduce the rows of the shear table", {
  skip_on_cran()
  Shear <- drawSeeded(drawInputs())
  Draws <- Shear$ShearDraws
  # Newest first, as the summary of newmark 1.7.1 accumulated them.
  OUT <- summarizeProfiles(Draws[rev(seq_len(nrow(Draws)))], levels = LEVELS[c(1:4, 6:9, 5)],
                           properties = c("Go", "mo", "VSo", "VS30", "Z500", "Z1000"))
  OUT[, `:=`(Go = round(Go, 1), mo = round(mo, 3), VSo = round(VSo, 1), VS30 = round(VS30, 1), Z500 = round(Z500, 1))]
  OUT[level != "mean", Z1000 := round(Z1000, 1)]
  Columns <- c("level", "Go", "mo", "VSo", "VS30", "Z500", "Z1000")
  expectReproduces(OUT[, ..Columns], Shear$ShearTable[, ..Columns])
  # lo is rounded to three decimals by the product contract.
  expect_identical(unique(Draws$lo), unique(Shear$ShearTable$lo))
  # The period of a summary row is the beam on that row, not a statistic
  # of the periods of the realizations.
  expect_equal(Shear$ShearTable$Ts, round(evaluateBeamPeriod(Hs = 10, lo = Shear$ShearTable$lo, mo = Shear$ShearTable$mo,
                                                                 VSo = Shear$ShearTable$VSo, an = Shear$ShearTable$an), 3))
})

test_that("the period of a realization is the beam evaluated on its own row, with the exact root outside the surrogate", {
  skip_on_cran()
  # lo is 8 / (2 * 10 * 2.5 + 8) = 0.138 for In and 30 / (2 * 10 * 1 + 30) =
  # 0.6, above the surrogate maximum 0.5, for Out, whose realizations all
  # take the exact root.
  Inputs <- drawInputs()
  Inputs$geometry <- data.frame(IDg = c("In", "Out"), Hs = 10, b = c(8, 30), s = c(2.5, 1))
  Draws <- drawSeeded(Inputs)$ShearDraws
  expect_equal(Draws$lo, rep(c(0.138, 0.6), each = 100L))
  Roots <- evaluateRoot(buildRootModel(mode = 1, method = "nlm"), mo = Draws$mo, lo = Draws$lo)
  Exact <- evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = Draws$mo, lo = Draws$lo)
  expect_identical(Roots$withinSupport, rep(c(TRUE, FALSE), each = 100L))
  # Outside its support the surrogate gives other roots, so the exact one is
  # observed rather than assumed.
  # Exact by contract: the stored roots are rounded to four decimals, and the
  # surrogate and exact roots differ at that rounding outside the support.
  expect_true(all(round(Roots$an[101:200], 4) != round(Exact$an[101:200], 4)))
  expect_equal(Draws$an, round(ifelse(Roots$withinSupport, Roots$an, Exact$an), 4))
  expect_equal(Draws$Ts, round(evaluateBeamPeriod(Hs = Draws$Hs, lo = Draws$lo, mo = Draws$mo, VSo = Draws$VSo,
                                                      an = Draws$an), 3))
})

test_that("realizations beyond 100 are drawn after every pair has its first 100, which they leave unchanged", {
  skip_on_cran()
  Inputs <- drawInputs(Hs = c(10, 12), lo = c(0.35, 0.075), s = c(2.5, 3.5))
  Base <- drawSeeded(Inputs)
  Extended <- drawSeeded(Inputs, NS = 120)
  expectReproduces(Extended$ShearTable, Base$ShearTable)
  # Every numeric property of ShearDraws is assigned/copied or contract-rounded.
  expect_identical(Extended$ShearDraws[realization <= 100L], Base$ShearDraws)
  expect_identical(Extended$ShearDraws$realization, rep(1:120, 2L))
  expect_identical(Extended$ShearDraws[, .N, by = .(IDg, IDm)]$N, c(120L, 120L))
  # The second pass continues the stream where the base run of every pair
  # ends: 20 columns per pair, drawn pair by pair as buildShearTable draws.
  set.seed(20260915L)
  invisible(do.call(buildShearDraws, Inputs))
  M1 <- data.frame(IP = c(0, 20, 40, 60, 80, 100), m1 = c(0, 0.18, 0.30, 0.41, 0.48, 0.48))
  for (Pair in unique(Base$ShearTable$IDg)) {
    Hs <- Base$ShearTable[IDg == Pair, Hs][[1L]]
    Profiles <- .legacyProfile(simulateSoilProfiles(Hs = Hs, USCS = c("ML", "MH"), NS = 20, h = 1, water = 0, POP = 100,
                                                    Vref = 760, m1 = M1))
    Columns <- c("Go", "mo", "VSo", "VS30", "Z500", "Z1000")
    # These six profile properties carry the rounding of .legacyProfile.
    expect_identical(as.list(Extended$ShearDraws[IDg == Pair & realization > 100L, ..Columns]),
                     as.list(Profiles[, ..Columns]))
  }
})

test_that("NS and the domain of the root are checked", {
  Inputs <- drawInputs()
  for (NS in list(99, 100.5, c(100, 120), NA_real_, "100")) {
    expect_error(do.call(buildShearDraws, c(Inputs, list(NS = NS))), "NS must be one whole number of at least 100",
                 fixed = TRUE)
  }
  expect_error(buildShearDraws(Inputs$geometry, list(U1 = "ML"), "missing", 760, "mean"), "IDm must select")
})

# A clay column of 10 m: four of its realizations under this seed fit a
# stiffness that decreases with depth, while its mean row has a root.
clayInputs <- function(levels = "mean") {
  list(geometry = buildGeometry(Hs = 10, lo = 0.35, s = 2.5), USCS = list(U1 = "CH"), IDm = "U1", Vref = 760,
       levels = levels)
}

test_that("a realization whose fitted mo is negative is returned with mo 0, the root and period at 0, and counted", {
  skip_on_cran()
  # The realizations of the clay column as the 1.7.1 generator stores them,
  # drawn alone from the same stream.
  set.seed(20260915L)
  Fitted <- .legacyProfile(simulateSoilProfiles(Hs = 10, USCS = "CH", NS = 100, h = 1, water = 0, POP = 100, Vref = 760,
                                                m1 = data.frame(IP = c(0, 20, 40, 60, 80, 100),
                                                                m1 = c(0, 0.18, 0.30, 0.41, 0.48, 0.48))))
  Negative <- Fitted$mo < 0
  expect_identical(Fitted$realization[Negative], c(16L, 26L, 32L, 62L))
  Notes <- list()
  set.seed(20260915L)
  Shear <- withCallingHandlers(do.call(buildShearDraws, clayInputs()), moduleNote = function(m) {
    Notes[[length(Notes) + 1L]] <<- m
    invokeRestart("muffleMessage")
  })
  expect_length(Notes, 1L)
  expect_identical(Notes[[1L]]$operation, "buildShearDraws")
  expect_identical(Notes[[1L]]$count, sum(Negative))
  expect_identical(conditionMessage(Notes[[1L]]),
                   paste0("Negative mo forced to 0 in 4 realization(s): IDg=S10 IDm=U1 realization=16 mo=-0.051; ",
                          "IDg=S10 IDm=U1 realization=26 mo=-0.007; IDg=S10 IDm=U1 realization=32 mo=-0.013; ",
                          "IDg=S10 IDm=U1 realization=62 mo=-0.023\n"))
  Draws <- Shear$ShearDraws
  # mo is contract-rounded to three decimals, with 0 assigned to negative fits.
  expect_identical(Draws$mo, ifelse(Negative, 0, Fitted$mo))
  # The forced rows take the root and the period of the homogeneous beam at
  # the truncation ratio of their geometry.
  Root <- evaluateRoot(buildRootModel(mode = 1, method = "nlm"), mo = 0, lo = Draws$lo[[1L]])
  expect_true(Root$withinSupport)
  expect_equal(Draws$an[Negative], rep(round(Root$an, 4), 4L))
  expect_equal(Draws$Ts[Negative], round(evaluateBeamPeriod(Hs = 10, lo = Draws$lo[Negative], mo = 0,
                                                                VSo = Fitted$VSo[Negative], an = round(Root$an, 4)), 3))
  # The summary keeps the fitted mo: its row is the one of buildShearTable().
  set.seed(20260915L)
  expectReproduces(Shear$ShearTable, do.call(buildShearTable, clayInputs()))
  # With the nine levels of the ShearTable step, quantile rows of the summary
  # have a negative mo, which has no root, and both functions stop alike.
  set.seed(20260915L)
  Table <- tryCatch(do.call(buildShearTable, clayInputs(LEVELS)), error = conditionMessage)
  expect_identical(Table, "the exact method needs lo in [0, 1) and mo in [0, 2)")
  set.seed(20260915L)
  expect_error(do.call(buildShearDraws, clayInputs(LEVELS)), Table, fixed = TRUE)
})

test_that("a realization with mo of 2 or more stops buildShearDraws naming it, before the summary rows are evaluated", {
  skip_on_cran()
  # Realization 7 of the clay column set beyond the relation, under the nine
  # levels whose quantile rows have no root either.
  Legacy <- .legacyProfile
  local_mocked_bindings(.legacyProfile = function(profiles) {
    OUT <- Legacy(profiles)
    OUT[realization == 7L, mo := 2.1]
    OUT
  })
  set.seed(20260915L)
  expect_error(do.call(buildShearDraws, clayInputs(LEVELS)),
               "1 realization(s) have mo of 2 or more, where the first-mode root is not defined: IDg=S10 IDm=U1 realization=7 mo=2.1",
               fixed = TRUE)
})
