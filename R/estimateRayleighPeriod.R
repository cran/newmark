#' Fundamental period of a layered soil profile by the Rayleigh quotient
#'
#' @section Mathematics:
#' For a fixed-base, free-surface column of uniform density under a uniform
#' horizontal body force, the shear strain at depth \eqn{z} is proportional
#' to \eqn{z / V_S(z)^2}. Integrating it from the base gives the static
#' deflection, used as trial mode shape \eqn{X_k} at each layer midpoint
#' (Dobry, Oweis and Urzua 1976), and the Rayleigh quotient yields
#' \deqn{T_r = 2\pi \sqrt{\frac{\sum_k h_k X_k^2}{\sum_k h_k X_k}}.}
#' Layers are ordered from the base internally, so the input order does not
#' change the result. The midpoint approximation depends on the subdivision
#' of the profile: for a uniform column of height \eqn{H} the estimate tends
#' to \eqn{3.974 H/V_S} as the layers become thin, against the exact
#' \eqn{4 H/V_S}, while one single layer gives \eqn{\pi H/V_S}. It is a
#' conditional evaluation of the
#' supplied profile and carries no dispersion.
#'
#' @param VSm Positive layer shear-wave velocities in m/s.
#' @param hs Positive layer thicknesses in metres, aligned with `VSm`.
#' @param zm Layer midpoint depths from the surface in metres, aligned.
#' @return One estimate of the fundamental period in seconds, unrounded.
#' @export
#' @examples
#' estimateRayleighPeriod(VSm = c(200, 200), hs = c(10, 10), zm = c(5, 15))
estimateRayleighPeriod <- function(VSm, hs, zm) {
  if (!is.numeric(VSm) || !is.numeric(hs) || !is.numeric(zm) || !length(hs) ||
      length(hs) != length(VSm) || length(hs) != length(zm) ||
      any(!is.finite(c(VSm, hs, zm))) || any(c(VSm, hs, zm) <= 0)) {
    stop("VSm, hs and zm must be finite positive vectors of equal length")
  }
  IDX <- order(zm, decreasing = TRUE)
  h <- hs[IDX]
  g <- zm[IDX] / VSm[IDX]^2
  X <- cumsum(g * h) - g * h / 2
  2 * pi * sqrt(sum(h * X^2) / sum(h * X))
}
