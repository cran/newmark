#' Site's Fundamental Period
#'
#' Fundamental period of a layered soil column (fixed base, free
#' surface) by the classical Rayleigh quotient with the static-
#' deflection trial shape (e.g., Dobry, Oweis & Urzua 1976): under a
#' uniform horizontal body force the shear strain at depth z is
#' z / Vs(z)^2 (uniform density), the mode shape is its integral from
#' the base, and \deqn{T_s = 2\pi\sqrt{\sum h_i X_i^2 / \sum h_i X_i}}
#' with X_i the shape at each layer midpoint. For a uniform profile
#' this yields 3.974 H/Vs (exact: 4 H/Vs, -0.7 percent); across layered
#' test profiles the error against the exact transfer-matrix solution
#' is a few percent on the low side. Layer order (top-down or
#' bottom-up) is irrelevant: layers are sorted internally by depth.
#'
#' @param VSm Vector. Shear Wave Velocity profile in m/s (vector)
#' @param hs Vector. Layer Thickness profile in m (vector)
#' @param zm Vector. Layer mid-depth coordinates from the surface in m
#'   (vector)
#'
#' @return Double. Site's fundamental period Ts in seconds
#' @export fitModel.Ts
#'
#' @import data.table
#'
fitModel.Ts <- function(VSm, hs, zm) { # nolint: object_name_linter.
  IDX <- order(zm, decreasing = TRUE) # base first
  h <- hs[IDX]
  g <- zm[IDX] / VSm[IDX]^2 # static shear strain at layer mid
  X <- cumsum(g * h) - g * h / 2 # mode shape at layer mid, from base
  Ts <- 2 * pi * sqrt(sum(h * X^2) / sum(h * X))
  round(Ts, 3)
}
