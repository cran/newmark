#' Fundamental periods of the soil profile and truncated shear beam for every realization
#'
#' Applies [solveWedgePeriod()] and [solveColumnPeriod()] to each layer table
#' of a list, with the paired height, truncation ratio, exponent and basal
#' velocity of the same realization, so that the numerical periods and their
#' controls stay paired with the realization before any summary.
#'
#' @param layers A non-empty list of layer tables, such as the `layers`
#'   column of [simulateSoilProfiles()]; see [solveColumnPeriod()].
#' @param Hs Slope heights in metres: one for every table, or one per table.
#' @param lo Truncation ratios in \eqn{[0, 1)}: one for every table, or one
#'   per table.
#' @param mo Inhomogeneity exponents of the analytic route, one per table.
#' @param VSo Basal shear-wave velocities of the analytic route in m/s, one
#'   per table.
#' @param method The solution method of the root model: `"exact"`, or the surrogate `"lm"` or `"nlm"`; see [buildRootModel()].
#' @param elements One positive whole number of linear elements per layer of
#'   the declared mesh; the control doubles it.
#' @return A data.table with one row per layer table, in the order of the
#'   list, and the columns of [solveWedgePeriod()] followed by those of
#'   [solveColumnPeriod()]. An invalid table stops the call with its position
#'   in the list.
#' @export
#' @examples
#' set.seed(20260915)
#' Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 2, h = 1,
#'                                  water = 0, POP = 100, Vref = 760)
#' solveProfilePeriods(Profiles$layers, Hs = 10, lo = 0.351, mo = Profiles$mo,
#'                     VSo = Profiles$VSo, method = "exact", elements = 2)
solveProfilePeriods <- function(layers, Hs, lo, mo, VSo, method, elements) {
  if (!is.list(layers) || is.data.frame(layers) || !length(layers)) {
    stop("layers must be a non-empty list of layer tables")
  }
  n <- length(layers)
  if (!length(Hs) %in% c(1L, n) || !length(lo) %in% c(1L, n)) {
    stop("Hs and lo must hold one value or one value per layer table")
  }
  if (length(mo) != n || length(VSo) != n) stop("mo and VSo must hold one value per layer table")
  Hs <- rep_len(Hs, n)
  lo <- rep_len(lo, n)
  rbindlist(lapply(seq_len(n), function(i) {
    tryCatch(
      cbind(solveWedgePeriod(layers[[i]], Hs = Hs[i], lo = lo[i], mo = mo[i], VSo = VSo[i],
                             method = method, elements = elements),
            solveColumnPeriod(layers[[i]], elements = elements)),
      error = function(e) stop("layer table ", i, ": ", conditionMessage(e), call. = FALSE)
    )
  }))
}
