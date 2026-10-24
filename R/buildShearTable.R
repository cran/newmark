#' Shear table of release 1.7.1 from geometries and material pools
#'
#' Compatibility route that reproduces, value for value, the `ShearTable`
#' product of the `getTs()` function of `newmark` 1.7.1 with its slope-period
#' height correction, so that the corrected library can serve as comparison
#' oracle while the calculations migrate. It
#' crosses every geometry with every selected material (or only the supplied
#' `pairs`), simulates 100 synthetic columns of 1 m layers per pair with a
#' dry column, 100 kPa of preconsolidation and a 760 m/s completion velocity,
#' summarizes their properties marginally and then evaluates the first-mode
#' truncated-beam period on each summary row with the apex-to-base height
#' \eqn{H_s/(1 - \lambda_o)} of [evaluateBeamPeriod()]. Shear tables built
#' before that correction took \eqn{H_s} itself, so their unrounded periods
#' are \eqn{1 - \lambda_o} times these. The first-mode root is the `"nlm"`
#' surrogate of [buildRootModel()] where `mo` and `lo` lie within its
#' tabulated support and the exact root of the characteristic relation
#' outside it, where `newmark` 1.7.1 extrapolated the surrogate.
#'
#' @section Semantics reproduced from `newmark` 1.7.1:
#' The period of a summary row is the beam expression evaluated on the
#' summarized `mo` and `VSo` of that row (`summarize-then-evaluate`), not a
#' quantile of the period population; [propagateBeamPeriod()] gives the
#' latter. The completion velocity of the upper 30 m is 760 m/s whatever
#' `Vref` is, because `getTs()` did not forward it; `Vref` is only recorded.
#' The layers take the plasticity exponent \eqn{m_1} of 1.7.1, interpolated
#' through `IP` 0, 20, 40, 60, 80 and 100 with 0, 0.18, 0.30, 0.41, 0.48 and
#' 0.48 and `NA` above 100, instead of the published table that
#' [evaluateSoilLayers()] uses by default. Before the summary, every realization takes the observables of the 1.7.1
#' generator: `Go` rounded to 2 decimals, `mo` and the Rayleigh period
#' (`TrRayleigh` of [characterizeProfile()], which 1.7.1 stored as `Ts`) to 3,
#' `VS30` to 1, `VSo` derived from the rounded `Go` and rounded to 1, and the
#' threshold depth ratios `Z500` and `Z1000`, the fitted normalized depth
#' \eqn{z/H_s} at which the inverse regression
#' \eqn{\log(z_m/H_s) = a + c \log V_S} reaches 500 and 1000 m/s, rounded
#' to 1 (a dimensionless ratio that 1.7.1 labelled as metres). The summaries
#' then round `Go` to 1 decimal, `mo` to 3, `Ts` to 2, `VSo`, `VS30` and
#' `Z500` to 1, and the quantile rows of `Z1000` to 1 while its mean row is
#' unrounded; `lo` is rounded to 3 decimals before the root and the period
#' use it, the root to 4 decimals and the beam period to 3.
#' The realizations are summarized newest first, as the original
#' accumulation did, so the means reproduce exactly. The general route
#' ([characterizeProfile()], [summarizeProfiles()]) rounds nothing and does
#' not compute the threshold ratios.
#'
#' @param geometry A data.frame or data.table with `IDg`, `Hs`, `b` and
#'   `s`, such as [buildGeometry()] returns.
#' @param USCS Named list of USCS code vectors, indexed by material `IDm`
#'   (the `uscs` argument of `getTs()`).
#' @param IDm Character vector of selected material identifiers.
#' @param Vref One positive reference velocity in m/s, recorded in the table.
#' @param levels Quantile levels and/or the literal `"mean"`.
#' @param pairs Optional table of unique `IDg`/`IDm` pairs to evaluate, in the
#'   supplied order; `NULL` evaluates the full geometry by material product.
#' @return A data.table with `Hs`, `b`, `s`, `Vref`, `IDm`, `IDg`, `POP`,
#'   `Hw`, `USCS`, `Go`, `mo`, `Ts` (fundamental period of the shear beam, s), `VSo`, `VS30`, `Z500`,
#'   `Z1000`, `level`, `lo` and `an`, one row per pair and level.
#' @export
#' @examples
#' \donttest{
#' Geometry <- buildGeometry(Hs = 10, lo = 0.35, s = 2.5)
#' set.seed(20260915L)
#' buildShearTable(Geometry, USCS = list(U1 = c("ML", "MH")), IDm = "U1",
#'                 Vref = 760, levels = c(0.16, "mean", 0.84))
#' }
buildShearTable <- function(geometry, USCS, IDm, Vref, levels, pairs = NULL) {
  Problem <- .diagnoseShearArguments(geometry, USCS = USCS, IDm = IDm, Vref = Vref, levels = levels, pairs = pairs)
  if (!is.null(Problem)) stop(Problem)
  # Resolved in the function scope, so a USCS column of the caller's table
  # cannot shadow the material pools inside the grouped evaluation.
  materialPool <- function(id) USCS[[id]]
  Pairs <- CJ(IDm, IDg = geometry$IDg)
  if (!is.null(pairs)) Pairs <- as.data.table(pairs)[, .(IDm, IDg)]
  # Quantile rows precede the mean row, as the original summary did.
  Levels <- c(levels[levels != "mean"], if (any(levels == "mean")) "mean")
  MaterialTable <- as.data.table(geometry)[Pairs, on = .(IDg), allow.cartesian = TRUE]
  set(MaterialTable, j = "Vref", value = Vref)
  ShearTable <- MaterialTable[, {
    Profiles <- .legacyProfile(simulateSoilProfiles(Hs = Hs, USCS = materialPool(IDm), NS = 100, h = 1,
                                                    water = 0, POP = 100, Vref = 760,
                                                    m1 = data.frame(IP = c(0, 20, 40, 60, 80, 100),
                                                                    m1 = c(0, 0.18, 0.30, 0.41, 0.48, 0.48))))
    OUT <- summarizeProfiles(Profiles[rev(seq_len(nrow(Profiles)))], levels = Levels,
                             properties = c("Go", "mo", "Ts", "VSo", "VS30", "Z500", "Z1000"))
    OUT[, `:=`(Go = round(Go, 1), mo = round(mo, 3), Ts = round(Ts, 2), VSo = round(VSo, 1),
               VS30 = round(VS30, 1), Z500 = round(Z500, 1))]
    OUT[level != "mean", Z1000 := round(Z1000, 1)]
    data.table(POP = 100, Hw = 0, USCS = paste0(materialPool(IDm), collapse = " "),
               OUT[, .(Go, mo, Ts, VSo, VS30, Z500, Z1000, level)])
  }, by = .(Hs, b, s, Vref, IDm, IDg)]
  ShearTable[, lo := round(b / (2 * Hs * s + b), 3)]
  Roots <- evaluateRoot(buildRootModel(mode = 1, method = "nlm"), mo = ShearTable$mo, lo = ShearTable$lo)
  Roots[!Roots$withinSupport, an := evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = mo, lo = lo)$an]
  ShearTable[, an := round(Roots$an, 4)]
  ShearTable[, Ts := round(evaluateBeamPeriod(Hs = Hs, lo = lo, mo = mo, VSo = VSo, an = an), 3)]
  ShearTable[]
}

# The first rule of buildShearTable() that its arguments break, as the
# message it stops with, or NULL when they keep every rule. The ShearTable
# step of processProject() applies the same rules to its contract before it
# seeds anything.
.diagnoseShearArguments <- function(geometry, USCS, IDm, Vref, levels, pairs) {
  Required <- c("IDg", "Hs", "b", "s")
  if (!is.data.frame(geometry) || !nrow(geometry) || !all(Required %in% names(geometry))) {
    return("geometry must contain rows with IDg, Hs, b and s")
  }
  if (!is.character(geometry$IDg) || anyNA(geometry$IDg) ||
      any(!nzchar(geometry$IDg)) || anyDuplicated(geometry$IDg)) {
    return("geometry requires unique nonempty IDg values")
  }
  for (Column in c("Hs", "b", "s")) {
    if (!is.numeric(geometry[[Column]]) || any(!is.finite(geometry[[Column]])) ||
        any(geometry[[Column]] <= 0)) return(paste0("geometry requires positive finite ", Column))
  }
  if (!is.character(IDm) || !length(IDm) || anyNA(IDm) || anyDuplicated(IDm) ||
      !is.list(USCS) || anyDuplicated(names(USCS)) || !all(IDm %in% names(USCS)) ||
      any(!vapply(USCS[IDm], function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x)), logical(1L)))) {
    return("IDm must select distinct named USCS material vectors")
  }
  if (!is.numeric(Vref) || length(Vref) != 1L || !is.finite(Vref) || Vref <= 0) {
    return("Vref must be one positive finite velocity")
  }
  if (!length(levels) || anyNA(levels)) return("levels must contain quantile levels and/or mean")
  if (is.null(pairs)) return(NULL)
  if (!is.data.frame(pairs) || !nrow(pairs) || !all(c("IDg", "IDm") %in% names(pairs))) {
    return("pairs must be a nonempty table with IDg and IDm")
  }
  Pairs <- as.data.table(pairs)[, .(IDm, IDg)]
  if (anyNA(Pairs) || anyDuplicated(Pairs) || !all(Pairs$IDg %in% geometry$IDg) || !all(Pairs$IDm %in% IDm)) {
    return("pairs must contain unique known geometry/material identifiers")
  }
  NULL
}

# Observables of one realization as the newmark 1.7.1 generator stored them:
# its rounding, the basal velocity derived from the rounded modulus and the
# threshold depth ratios of the inverse velocity regression.
.legacyProfile <- function(profiles) {
  OUT <- copy(as.data.table(profiles))
  Gsat <- vapply(OUT$layers, function(L) L$gsat[nrow(L)], numeric(1L))
  Ratios <- Map(.thresholdRatios, OUT$layers, OUT$Hs)
  OUT[, `:=`(Go = round(Go, 2), mo = round(mo, 3), VS30 = round(VS30, 1), Ts = round(TrRayleigh, 3))]
  OUT[, VSo := round(sqrt(9.81 * Go * 1000 / Gsat), 1)]
  OUT[, `:=`(Z500 = round(vapply(Ratios, `[`, numeric(1L), 1L), 1),
             Z1000 = round(vapply(Ratios, `[`, numeric(1L), 2L), 1))]
  OUT
}

.thresholdRatios <- function(layers, Hs) {
  AUX <- data.table(zm = layers$zm, VSm = layers$VSm)[VSm > 0][, .(X = log(zm / Hs), Z = log(VSm))]
  Inverse <- lm(X ~ Z, data = AUX)
  c(exp(unname(predict(Inverse, newdata = data.table(Z = log(500))))),
    exp(unname(predict(Inverse, newdata = data.table(Z = log(1000))))))
}
