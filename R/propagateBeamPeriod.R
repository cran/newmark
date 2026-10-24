#' Fundamental period of the truncated shear beam for every profile realization
#'
#' Solves, for each realization of [simulateSoilProfiles()], the fundamental period
#' of the truncated shear beam directly from the realization's own layers
#' with [solveWedgePeriod()], and keeps beside it the analytic period of the
#' truncated shear beam evaluated on the realization's fitted `mo` and `VSo`,
#' with the controls of both. The resulting population of periods stays
#' paired with the other properties of each realization: it can be
#' summarized with [summarizeProfiles()], naming `"Ts"` or `"TsAnalytic"`
#' among its `properties`, or propagated further. This is the
#' calculation-before-summary route; the compatibility route of
#' [buildShearTable()] evaluates the analytic period on summary rows instead.
#'
#' @section Settings:
#' Every layer of a realization is divided into two linear elements, and the
#' mesh-doubling control of [solveWedgePeriod()] repeats the solve with four,
#' so `Ts` is the period on four elements per layer. A realization is flagged
#' divergent when the analytic period differs from the numerical one by more
#' than 10 percent, \eqn{|T_{s,analytic}/T_s - 1| > 0.10}; the flag rejects
#' nothing. Both values are design settings.
#'
#' @section Height:
#' The \eqn{\lceil H_s / h \rceil} layers of [sampleSoilLayers()] reach below
#' `Hs` when `Hs` is not a multiple of `h`. Before the numerical period, the
#' deepest layer of such a realization (largest `zm`) is trimmed so that the
#' thicknesses add up to `Hs`: the trimmed geometry is the one `Ts` uses.
#' Only that thickness changes; the properties of every layer stay those
#' sampled, and `TsAnalytic` keeps the `Hs`, `mo` and `VSo` of the row. The
#' `layers` column is returned untrimmed, so the site period of
#' [propagateSitePeriod()] is that of the sampled column. A sum that already
#' meets the equality of [solveWedgePeriod()] is solved as it is; an excess
#' not smaller than the deepest layer stops the call with the row. This is a
#' design decision.
#'
#' @section Analytic route:
#' `TsAnalytic` takes the apex-to-base height \eqn{H_s/(1 - \lambda_o)} of
#' [evaluateBeamPeriod()] and the first characteristic root at `mo` and `lo`
#' from [evaluateRoot()] with `method`. Each `mo` is the exponent of the
#' profile that [characterizeProfile()] fits from the crest, used with the
#' roots of the apex-origin profile of Dakoulas and Gazetas (1985): an
#' approximation that the theory base measures on a continuous profile at
#' 1.4 to 6.0 percent for \eqn{m = 0.5} and \eqn{\lambda} from 0.075 to
#' 0.35, and whose size for each realization is `TsAnalyticDifference`,
#' flagged by `TsDivergent`. `Ts` makes no such approximation: it is the
#' period of the realization's layered shear beam.
#'
#' @param profiles The result of [simulateSoilProfiles()] or
#'   [propagateSitePeriod()]: one row per realization with columns `Hs`, `mo`,
#'   `VSo` and `layers`, the list column of layer tables whose thicknesses add
#'   up to `Hs` or exceed it by less than the deepest layer (section Height;
#'   see [solveWedgePeriod()]).
#' @param lo One truncation ratio of the geometry strictly between 0 and 1,
#'   see [buildGeometry()].
#' @param method The solution method of the root model of the analytic route:
#'   `"exact"`, or the surrogate `"lm"` or `"nlm"`; see [buildRootModel()].
#' @return A copy of `profiles` with these columns appended, one value per
#'   realization paired with the other properties of its row: `lo`; `Ts`
#'   (s), the numerical fundamental period of the realization's layers trimmed to
#'   `Hs`;
#'   `TsMeshChange`, its relative change under mesh doubling; `TsAnalytic`
#'   (s), the period of the analytic route; `TsAnalyticDifference`,
#'   `TsAnalytic / Ts - 1`; `TsDivergent`, `TRUE` when
#'   `abs(TsAnalyticDifference)` exceeds 0.10; and `an`, `anExact` and
#'   `anWithinSupport`, the root of the analytic route and its diagnostics.
#'   The population is that of `profiles`; `TsMeshChange`,
#'   `TsAnalyticDifference` and `TsDivergent` are numerical diagnostics, not
#'   physical dispersion. An invalid realization stops the call with its row.
#' @export
#' @examples
#' set.seed(20260915)
#' Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 3, h = 1,
#'                                  water = 0, POP = 100, Vref = 760)
#' Periods <- propagateBeamPeriod(Profiles, lo = 0.351, method = "nlm")
#' Periods[, .(mo, VSo, Ts, TsAnalytic, TsDivergent)]
propagateBeamPeriod <- function(profiles, lo, method) {
  Required <- c("Hs", "mo", "VSo", "layers")
  if (!is.data.frame(profiles) || !nrow(profiles) || !all(Required %in% names(profiles))) {
    stop("profiles must contain rows with columns ", paste(Required, collapse = ", "))
  }
  if (!is.numeric(lo) || length(lo) != 1L || !is.finite(lo) || lo <= 0 || lo >= 1) {
    stop("lo must be one truncation ratio strictly between 0 and 1")
  }
  OUT <- copy(as.data.table(profiles))
  Periods <- rbindlist(lapply(seq_len(nrow(OUT)), function(i) {
    tryCatch(
      solveWedgePeriod(.trimLayers(OUT$layers[[i]], Hs = OUT$Hs[[i]]), Hs = OUT$Hs[[i]], lo = lo,
                       mo = OUT$mo[[i]], VSo = OUT$VSo[[i]], method = method, elements = PERIODELEMENTS),
      error = function(e) stop("realization row ", i, ": ", conditionMessage(e), call. = FALSE)
    )
  }))
  set(Periods, j = "TsDivergent", value = abs(Periods$TsAnalyticDifference) > TSDIVERGENCE)
  setcolorder(Periods, c("Ts", "TsMeshChange", "TsAnalytic", "TsAnalyticDifference", "TsDivergent"))
  set(OUT, j = "lo", value = lo)
  set(OUT, j = names(Periods), value = Periods)
  OUT[]
}

# Design settings for the periods of the realization route: linear
# elements per layer of the declared mesh of solveWedgePeriod() and
# solveColumnPeriod(), and the relative difference between the analytic
# and the numerical slope period beyond which a realization is flagged
# (flag only, no rejection).
PERIODELEMENTS <- 2
TSDIVERGENCE <- 0.10

# The layers of one realization cut at its height (design decision): when
# the sampled thicknesses exceed Hs beyond the equality
# solveWedgePeriod() requires, the deepest layer loses the excess, its other
# properties unchanged. Any other table is returned as it is, for
# solveWedgePeriod() to solve or reject; .orderLayers() validates it first,
# so an invalid table fails as there.
.trimLayers <- function(layers, Hs) {
  Layers <- .orderLayers(layers)
  Total <- sum(Layers$hs)
  if (!isTRUE(Total > Hs) || isTRUE(all.equal(Total, Hs))) return(layers)
  n <- length(Layers$hs)
  Thickness <- Hs - sum(Layers$hs[-n])
  if (Thickness <= 0) {
    stop("the layer thicknesses add up to ", Total, " m, exceeding Hs = ", Hs,
         " m by at least the deepest layer of ", Layers$hs[[n]], " m")
  }
  Layers$hs[[n]] <- Thickness
  as.data.table(Layers)
}
