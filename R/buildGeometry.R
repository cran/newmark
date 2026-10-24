#' Paired embankment geometries from height, truncation ratio and slope
#'
#' @section Mathematics:
#' The truncation ratio of the shear beam is \eqn{\lambda = b / (2 H_s s + b)}
#' for a berm width \eqn{b}, height \eqn{H_s} and horizontal-to-vertical
#' slope \eqn{s}; inverting it gives \eqn{b = 2 H_s s / (1/\lambda - 1)}. The
#' base width is \eqn{b + 2 H_s s} and the slope angle
#' \eqn{\beta = \arctan(1/s)} in degrees. Widths and heights are rounded to
#' metres and slopes to two decimals, as in `newmark` 1.7.1; the identifier
#' `IDg` is `S` followed by the rounded height.
#'
#' @param Hs Positive heights in metres.
#' @param lo Truncation ratios strictly between 0 and 1, aligned with `Hs`.
#' @param s Positive horizontal-to-vertical slopes, aligned with `Hs`.
#' @return A data.table with `IDg`, `b`, `bmax`, `s`, `Hs` and `beta`.
#' @export
#' @examples
#' buildGeometry(Hs = c(10, 120), lo = c(0.35, 0.075), s = c(2.5, 3.5))
buildGeometry <- function(Hs, lo, s) {
  if (!is.numeric(Hs) || !is.numeric(lo) || !is.numeric(s) || !length(Hs) ||
      length(Hs) != length(lo) || length(Hs) != length(s) ||
      any(!is.finite(c(Hs, lo, s))) || any(Hs <= 0 | s <= 0 | lo <= 0 | lo >= 1)) {
    stop("Hs, lo and s must be finite, equally sized vectors; Hs/s positive and 0 < lo < 1")
  }
  Width <- 2 * Hs * s / (1 / lo - 1)
  OUT <- data.table(
    IDg = paste0("S", round(Hs, 0)),
    b = round(Width, 0), bmax = round(Width + 2 * Hs * s, 0),
    s = round(s, 2), Hs = round(Hs, 0), beta = atan(1 / s) * 180 / pi
  )
  if (anyDuplicated(OUT$IDg) || any(OUT$Hs <= 0 | OUT$s <= 0)) {
    stop("rounded geometries must have distinct IDg and positive Hs/s")
  }
  OUT
}
