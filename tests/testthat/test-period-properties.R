# Properties of the numerical periods on random layered profiles (owner's
# approval of 2026-10-12: any positive modulus, density, thickness and
# truncation ratio) and the vectorised wrapper solveProfilePeriods.

# Mesh change declared for this test with 4 elements per layer: a level the
# solver must stay below on every random profile, not a threshold of the
# product (the owner sets those).
MESHLEVEL <- 0.01

randomLayers <- function(N) {
  hs <- exp(runif(N, log(0.2), log(20)))
  data.frame(hs = hs, zm = cumsum(hs) - hs / 2, gsat = runif(N, 12, 24), VSm = exp(runif(N, log(50), log(2000))))
}

test_that("random columns give finite positive periods within the declared mesh level and the Rayleigh bound", {
  # Seed of the chain's validated contract (owner, 2026-10-12: semilla 20260915).
  set.seed(20260915L)
  for (Case in 1:100) {
    Layers <- randomLayers(sample(1:25, 1))
    OUT <- solveColumnPeriod(Layers, elements = 4)
    expect_true(is.finite(OUT$Tr) && OUT$Tr > 0)
    expect_gte(OUT$TrMeshChange, 0)
    expect_lt(OUT$TrMeshChange, MESHLEVEL)
    # Rayleigh bound: TrRayleighWeighted <= T exact <= Tr (1 + TrMeshChange).
    expect_lte(OUT$TrRayleighWeightedDifference, OUT$TrMeshChange)
    # Invariance to the row order given the depths, and to the unit of gsat.
    Shuffled <- Layers[sample(nrow(Layers)), ]
    expect_identical(solveColumnPeriod(Shuffled, elements = 4), OUT)
    Heavier <- Layers
    Heavier$gsat <- 2 * Heavier$gsat
    expect_identical(solveColumnPeriod(Heavier, elements = 4), OUT)
  }
})

test_that("random wedges give finite positive periods within the declared mesh level", {
  set.seed(20260915L)
  for (Case in 1:100) {
    Layers <- randomLayers(sample(1:25, 1))
    lo <- if (Case <= 10) 0 else runif(1, 0, 0.95)
    OUT <- solveWedgePeriod(Layers, Hs = sum(Layers$hs), lo = lo, mo = 0.5, VSo = 300, method = "exact",
                            elements = 4)
    expect_true(is.finite(OUT$Ts) && OUT$Ts > 0)
    expect_gte(OUT$TsMeshChange, 0)
    expect_lt(OUT$TsMeshChange, MESHLEVEL)
    Shuffled <- Layers[sample(nrow(Layers)), ]
    expect_identical(solveWedgePeriod(Shuffled, Hs = sum(Layers$hs), lo = lo, mo = 0.5, VSo = 300,
                                      method = "exact", elements = 4), OUT)
  }
})

test_that("solveProfilePeriods pairs every realization with its own inputs", {
  set.seed(20260915L)
  Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 3, h = 1,
                                   water = 0, POP = 100, Vref = 760)
  OUT <- solveProfilePeriods(Profiles$layers, Hs = 10, lo = 0.351, mo = Profiles$mo, VSo = Profiles$VSo,
                             method = "nlm", elements = 2)
  expect_identical(nrow(OUT), 3L)
  for (i in 1:3) {
    Single <- cbind(solveWedgePeriod(Profiles$layers[[i]], Hs = 10, lo = 0.351, mo = Profiles$mo[i],
                                     VSo = Profiles$VSo[i], method = "nlm", elements = 2),
                    solveColumnPeriod(Profiles$layers[[i]], elements = 2))
    expect_identical(OUT[i], Single)
  }
  expect_true(all(OUT$TrRayleighWeightedDifference <= OUT$TrMeshChange))
  # The slope periods are those of the realization route on the same rows.
  Beam <- propagateBeamPeriod(Profiles, lo = 0.351, method = "nlm")
  expect_identical(OUT[, .(Ts, TsMeshChange, TsAnalytic, TsAnalyticDifference)],
                   Beam[, .(Ts, TsMeshChange, TsAnalytic, TsAnalyticDifference)])
  # Heights and truncation ratios may differ by realization.
  Mixed <- solveProfilePeriods(Profiles$layers, Hs = c(10, 10, 10), lo = c(0, 0.351, 0.5), mo = Profiles$mo,
                               VSo = Profiles$VSo, method = "exact", elements = 2)
  expect_identical(Mixed[3], cbind(solveWedgePeriod(Profiles$layers[[3]], Hs = 10, lo = 0.5, mo = Profiles$mo[3],
                                                    VSo = Profiles$VSo[3], method = "exact", elements = 2),
                                   solveColumnPeriod(Profiles$layers[[3]], elements = 2)))
  expect_error(solveProfilePeriods(Profiles$layers[[1]], Hs = 10, lo = 0.351, mo = 0.3, VSo = 200,
                                   method = "exact", elements = 2), "non-empty list")
  expect_error(solveProfilePeriods(Profiles$layers, Hs = c(10, 10), lo = 0.351, mo = Profiles$mo,
                                   VSo = Profiles$VSo, method = "exact", elements = 2), "one value per layer table")
  expect_error(solveProfilePeriods(Profiles$layers, Hs = 10, lo = 0.351, mo = Profiles$mo[1:2],
                                   VSo = Profiles$VSo, method = "exact", elements = 2), "mo and VSo")
  Broken <- Profiles$layers
  Broken[[2]] <- copy(Broken[[2]])[3, gsat := 0]
  expect_error(solveProfilePeriods(Broken, Hs = 10, lo = 0.351, mo = Profiles$mo, VSo = Profiles$VSo,
                                   method = "exact", elements = 2), "layer table 2: layers\\$gsat")
})
