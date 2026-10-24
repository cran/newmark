#' Time-averaged shear-wave velocity of the upper 30 m
#'
#' @section Mathematics:
#' With layer thicknesses \eqn{h_k} from the surface and velocities
#' \eqn{V_k}, the travel time through the part of each layer inside the upper
#' 30 m is accumulated and
#' \deqn{V_{S30} = \frac{30}{\sum_k h_k^{(30)} / V_k + \max(30 - \sum_k h_k, 0)/V_{ref}},}
#' where the second term completes a column shallower than 30 m with the
#' reference velocity `Vref` and vanishes otherwise. The result is the
#' harmonic average of the velocities weighted by thickness (ASCE/SEI 7-22
#' Eq. 20.4-1), evaluated at the supplied profile; it carries no dispersion.
#'
#' @param hs Positive layer thicknesses in metres, from the surface down.
#' @param VSm Positive layer shear-wave velocities in m/s, aligned with `hs`.
#' @param Vref One positive velocity in m/s of the material below the column,
#'   used only when the column is shallower than 30 m.
#' @return One velocity in m/s, unrounded.
#' @export
#' @examples
#' estimateVS30(hs = c(10, 20), VSm = c(200, 400), Vref = 760)   # 300
#' estimateVS30(hs = 10, VSm = 200, Vref = 760)                  # padded
estimateVS30 <- function(hs, VSm, Vref) {
  if (!is.numeric(hs) || !is.numeric(VSm) || !length(hs) || length(hs) != length(VSm) ||
      any(!is.finite(c(hs, VSm))) || any(c(hs, VSm) <= 0)) {
    stop("hs and VSm must be finite positive vectors of equal length")
  }
  if (!is.numeric(Vref) || length(Vref) != 1L || !is.finite(Vref) || Vref <= 0) {
    stop("Vref must be one finite positive velocity")
  }
  ZTop <- c(0, head(cumsum(hs), -1))
  ZBot <- ZTop + hs
  t30 <- pmax(pmin(ZBot, 30) - ZTop, 0) / VSm
  if (sum(hs) < 30) t30 <- c(t30, (30 - sum(hs)) / Vref)
  30 / sum(t30)
}
