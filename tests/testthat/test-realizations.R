test_that("simulateSoilProfiles retains paired realizations and their layers", {
  set.seed(20260915L)
  Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 3, h = 1,
                                   water = 0, POP = 100, Vref = 760)
  expect_identical(Profiles$realization, 1:3)
  expect_true(all(c("Go", "mo", "VSo", "VS30", "TrRayleigh", "layers") %in% names(Profiles)))
  expect_identical(vapply(Profiles$layers, nrow, 1L), rep(10L, 3))
  # The row properties come from the stored layers of that same row.
  Again <- characterizeProfile(Profiles$layers[[2L]], Hs = 10, Vref = 760)
  expect_identical(Profiles[2L, .(Go, mo, VSo, VS30, TrRayleigh)], Again[, .(Go, mo, VSo, VS30, TrRayleigh)])
})

test_that("summarizeProfiles returns the requested marginal statistics in order", {
  set.seed(20260915L)
  Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 5, h = 1,
                                   water = 0, POP = 100, Vref = 760)
  OUT <- summarizeProfiles(Profiles, levels = c(0.16, "mean", 0.84), properties = c("Go", "mo", "VSo"))
  expect_identical(OUT$level, c("0.16", "mean", "0.84"))
  expect_identical(OUT$n, rep(5L, 3))
  expect_identical(names(OUT), c("level", "n", "Go", "mo", "VSo"))
  expect_equal(OUT[level == "mean", Go], mean(Profiles$Go))
  expect_equal(OUT[level == "0.84", mo], quantile(Profiles$mo, 0.84, names = FALSE))
  expect_error(summarizeProfiles(Profiles, levels = c("median"), properties = "Go"), "probabilities")
  expect_error(summarizeProfiles(Profiles, levels = "mean", properties = "Z500"), "columns")
  expect_error(summarizeProfiles(Profiles, levels = "mean", properties = "SID"), "numeric")
})

test_that("propagateBeamPeriod solves the slope period of each realization from its layers", {
  set.seed(20260915L)
  Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 4, h = 1,
                                   water = 0, POP = 100, Vref = 760)
  OUT <- propagateBeamPeriod(Profiles, lo = 0.351, method = "nlm")
  expect_identical(nrow(OUT), 4L)
  # Each row is solveWedgePeriod on the layers of that realization with two
  # elements per layer (owner's setting of 2026-10-12).
  for (i in 1:4) {
    Single <- solveWedgePeriod(Profiles$layers[[i]], Hs = 10, lo = 0.351, mo = Profiles$mo[i],
                               VSo = Profiles$VSo[i], method = "nlm", elements = 2)
    expect_identical(OUT[i, names(Single), with = FALSE], Single)
  }
  # The analytic route takes the apex-to-base height Hs / (1 - lo), owner's
  # ruling of 2026-10-12.
  expect_equal(OUT$TsAnalytic, 4 * pi * (10 / (1 - 0.351)) / (OUT$an * (2 - OUT$mo) * OUT$VSo))
  # Divergence flag at 10 percent of relative difference (owner's setting), flag only.
  expect_identical(OUT$TsDivergent, abs(OUT$TsAnalytic / OUT$Ts - 1) > 0.10)
  expect_identical(OUT$lo, rep(0.351, 4))
  expect_true(all(OUT$anWithinSupport))
  # The caller's table keeps its columns and gains none.
  expect_identical(OUT[, names(Profiles), with = FALSE], Profiles)
  expect_false(any(c("Ts", "lo") %in% names(Profiles)))
  # The slope-period population is summarized by naming it.
  Periods <- summarizeProfiles(OUT, levels = c(0.16, "mean"), properties = c("Ts", "TsAnalytic"))
  expect_equal(Periods[level == "mean", Ts], mean(OUT$Ts))
  expect_equal(Periods[level == "0.16", TsAnalytic], quantile(OUT$TsAnalytic, 0.16, names = FALSE))
  # The population mean of the analytic period is not the period of the mean profile.
  Summary <- summarizeProfiles(Profiles, levels = "mean", properties = c("mo", "VSo"))
  Mean <- evaluateBeamPeriod(Hs = 10, lo = 0.351, mo = Summary$mo, VSo = Summary$VSo,
                             an = evaluateRoot(buildRootModel(1, "nlm"), Summary$mo, 0.351)$an)
  expect_false(isTRUE(all.equal(mean(OUT$TsAnalytic), Mean)))
  expect_error(propagateBeamPeriod(Profiles[, !"layers"], lo = 0.351, method = "nlm"), "layers")
  expect_error(propagateBeamPeriod(Profiles, lo = 0, method = "nlm"), "strictly between 0 and 1")
  Taller <- copy(Profiles)[2L, Hs := 12]
  expect_error(propagateBeamPeriod(Taller, lo = 0.351, method = "nlm"), "realization row 2: .*not to Hs = 12 m")
})

test_that("propagateBeamPeriod trims the deepest sampled layer to Hs before the numerical period", {
  # Owner's decision of 2026-10-13: Hs = 10.5 m is sampled as eleven layers of
  # 1 m, and the slope period is solved on the column cut at Hs.
  set.seed(20260915L)
  Profiles <- simulateSoilProfiles(Hs = 10.5, USCS = c("ML", "MH"), NS = 3, h = 1,
                                   water = 0, POP = 100, Vref = 760)
  Sampled <- copy(Profiles$layers)
  OUT <- propagateBeamPeriod(Profiles, lo = 0.351, method = "nlm")
  expect_true(all(is.finite(OUT$Ts) & OUT$Ts > 0))
  for (i in 1:3) {
    # Only the thickness of the deepest layer changes.
    Trimmed <- copy(Profiles$layers[[i]])[layer == 11L, hs := 0.5]
    Single <- solveWedgePeriod(Trimmed, Hs = 10.5, lo = 0.351, mo = Profiles$mo[i], VSo = Profiles$VSo[i],
                               method = "nlm", elements = 2)
    expect_identical(OUT[i, names(Single), with = FALSE], Single)
  }
  # The analytic route keeps the Hs, mo and VSo of the row.
  expect_equal(OUT$TsAnalytic, 4 * pi * (10.5 / (1 - 0.351)) / (OUT$an * (2 - OUT$mo) * OUT$VSo))
  # The sampled layers stay untrimmed, in the caller's table and in the
  # result, and the site period is that of the sampled column.
  expect_identical(Profiles$layers, Sampled)
  expect_identical(OUT$layers, Sampled)
  expect_identical(vapply(OUT$layers, function(L) sum(L$hs), numeric(1L)), rep(11, 3))
  expect_identical(propagateSitePeriod(OUT)$Tr, propagateSitePeriod(Profiles)$Tr)
  # An excess not smaller than the deepest layer stops the call with its row.
  Shorter <- copy(Profiles)[2L, Hs := 10]
  expect_error(propagateBeamPeriod(Shorter, lo = 0.351, method = "nlm"),
               "realization row 2: .*exceeding Hs = 10 m by at least the deepest layer of 1 m")
})

test_that("propagateSitePeriod solves the site period of each realization from its layers", {
  set.seed(20260915L)
  Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 3, h = 1,
                                   water = 0, POP = 100, Vref = 760)
  OUT <- propagateSitePeriod(Profiles)
  # Each row is solveColumnPeriod on the layers of that realization with two
  # elements per layer (owner's setting of 2026-10-12).
  for (i in 1:3) {
    expect_identical(OUT[i, .(Tr, TrMeshChange, TrRayleighWeighted, TrRayleighWeightedDifference)],
                     solveColumnPeriod(Profiles$layers[[i]], elements = 2))
  }
  # The weighted Rayleigh estimate is bounded by the period within its mesh change.
  expect_true(all(OUT$TrRayleighWeightedDifference <= OUT$TrMeshChange))
  # The constant-density estimate of characterizeProfile stays beside them.
  expect_identical(OUT[, names(Profiles), with = FALSE], Profiles)
  expect_false("Tr" %in% names(Profiles))
  # The site and slope periods do not depend on the order of the two steps.
  Both <- propagateBeamPeriod(OUT, lo = 0.351, method = "nlm")
  expect_identical(Both[, .(Tr, Ts, TsAnalytic)],
                   propagateSitePeriod(propagateBeamPeriod(Profiles, lo = 0.351, method = "nlm"))[, .(Tr, Ts, TsAnalytic)])
  expect_error(propagateSitePeriod(Profiles[, !"layers"]), "list column layers")
  Broken <- copy(Profiles)
  set(Broken, i = 3L, j = "layers", value = list(list(copy(Profiles$layers[[3]])[2, gsat := -1])))
  expect_error(propagateSitePeriod(Broken), "realization row 3: layers\\$gsat")
})
