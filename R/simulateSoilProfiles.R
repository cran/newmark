#' Simulate paired soil-profile realizations for a material pool
#'
#' Draws `NS` synthetic columns with [sampleSoilLayers()], evaluates each with
#' [evaluateSoilLayers()] and characterizes it with [characterizeProfile()].
#' Every realization is kept with its layers, so the profile quantities stay
#' paired: the `Go`, `mo`, `VSo`, `VS30` and `TrRayleigh` of one row belong
#' to the same column. [summarizeProfiles()] reduces this population to
#' marginal statistics; [propagateSitePeriod()] and [propagateBeamPeriod()]
#' solve the fundamental periods of the soil profile and truncated shear beam on each row before any reduction.
#'
#' @param Hs One positive height in metres.
#' @param USCS Character vector of USCS codes forming the material pool.
#' @param NS Number of realizations, a positive integer.
#' @param h One positive layer thickness in metres.
#' @param water Water-table height as a fraction of `Hs`.
#' @param POP Preconsolidation pressure in kPa.
#' @param Vref One positive reference velocity in m/s, see [estimateVS30()].
#' @param m1 Nodes of the plasticity exponent, see [evaluateSoilLayers()],
#'   whose published table is also the default here.
#' @return A data.table with one row per realization: `realization`, `Hs`,
#'   `Hw` (water-table height, m), `POP`, `Vref`, the columns of
#'   [characterizeProfile()] and `layers`, a list column holding the
#'   evaluated layer table of the row. The rows are the population generated
#'   by the pool under the caller's random stream.
#' @export
#' @examples
#' set.seed(20260915)
#' simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 3, h = 1,
#'                      water = 0, POP = 100, Vref = 760)
simulateSoilProfiles <- function(Hs, USCS, NS, h, water, POP, Vref, m1 = HardinExponent) {
  if (!is.numeric(NS) || length(NS) != 1L || !is.finite(NS) || NS < 1 || NS != as.integer(NS)) {
    stop("NS must be one positive integer")
  }
  LIST <- vector("list", NS)
  for (i in seq_len(NS)) {
    Layers <- evaluateSoilLayers(sampleSoilLayers(Hs = Hs, USCS = USCS, h = h),
                                 Hs = Hs, water = water, POP = POP, m1 = m1)
    Profile <- characterizeProfile(Layers, Hs = Hs, Vref = Vref)
    Profile[, `:=`(realization = i, Hs = Hs, Hw = Hs * water, POP = POP, Vref = Vref,
                   layers = list(Layers))]
    LIST[[i]] <- Profile
  }
  OUT <- rbindlist(LIST)
  setcolorder(OUT, c("realization", "Hs", "Hw", "POP", "Vref"))
  OUT[]
}
