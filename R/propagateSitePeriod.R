#' Fundamental period of the soil profile for every realization
#'
#' Solves, for each realization of [simulateSoilProfiles()], the fundamental period
#' \eqn{T_r} of its layered soil column on rigid base with [solveColumnPeriod()],
#' and keeps beside it the density-weighted Rayleigh estimate and its
#' relative difference. The resulting population of site periods stays
#' paired with the other properties of each realization: it can be
#' summarized with [summarizeProfiles()], naming `"Tr"` among its
#' `properties`, or propagated further.
#'
#' @section Settings:
#' Every layer of a realization is divided into two linear elements, and the
#' mesh-doubling control of [solveColumnPeriod()] repeats the solve with
#' four, so `Tr` is the period on four elements per layer: a design
#' setting shared with [propagateBeamPeriod()].
#'
#' @section Rayleigh estimates:
#' The column is the realization's layers stacked from the surface, so its
#' height is the sum of their thicknesses, which exceeds `Hs` when `Hs` is
#' not a multiple of the layer thickness of [sampleSoilLayers()].
#' `TrRayleighWeighted` never exceeds the exact period of that column (see
#' [solveColumnPeriod()]). The `TrRayleigh` of [characterizeProfile()],
#' already in the realization table, is the constant-density estimate of
#' `newmark` 1.7.1: it ignores the variation of the unit weight between
#' layers, so it is not a bound of `Tr` when the unit weight varies.
#'
#' @param profiles The result of [simulateSoilProfiles()] or
#'   [propagateBeamPeriod()]: one row per realization with the list column
#'   `layers` of layer tables (see [solveColumnPeriod()]).
#' @return A copy of `profiles` with these columns appended, one value per
#'   realization paired with the other properties of its row: `Tr` (s), the
#'   numerical site period of the realization's column; `TrMeshChange`, its
#'   relative change under mesh doubling; `TrRayleighWeighted` (s), the
#'   density-weighted Rayleigh period of [estimateWeightedRayleighPeriod()];
#'   and `TrRayleighWeightedDifference`, `TrRayleighWeighted / Tr - 1`. The
#'   population is that of `profiles`; `TrMeshChange` and
#'   `TrRayleighWeightedDifference` are numerical diagnostics, not physical
#'   dispersion. An invalid realization stops the call with its row.
#' @export
#' @examples
#' set.seed(20260915)
#' Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 3, h = 1,
#'                                  water = 0, POP = 100, Vref = 760)
#' propagateSitePeriod(Profiles)[, .(TrRayleigh, Tr, TrRayleighWeighted, TrMeshChange)]
propagateSitePeriod <- function(profiles) {
  if (!is.data.frame(profiles) || !nrow(profiles) || !"layers" %in% names(profiles)) {
    stop("profiles must contain rows with the list column layers")
  }
  OUT <- copy(as.data.table(profiles))
  Periods <- rbindlist(lapply(seq_len(nrow(OUT)), function(i) {
    tryCatch(solveColumnPeriod(OUT$layers[[i]], elements = PERIODELEMENTS),
             error = function(e) stop("realization row ", i, ": ", conditionMessage(e), call. = FALSE))
  }))
  set(OUT, j = names(Periods), value = Periods)
  OUT[]
}
