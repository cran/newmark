#' Site's Fundamental Period
#'
#' @param VSm Vector. Shear Wave Velocity profile in m/s (vector)
#' @param hs Vector. Layer Thickness profile in m (vector)
#' @param zm Vector. Layer Coordinates profile in m (vector)
#'
#' @return Double. Site's fundamental period Ts in seconds
#' @export fitModel.Ts
#'
#' @import data.table
#'
fitModel.Ts <- function(VSm, hs, zm) { # nolint: object_name_linter.
  NL <- length(VSm)
  A <- vector(mode = "double", length = NL - 1)
  B <- vector(mode = "double", length = NL - 1)
  f <- vector(mode = "double", length = NL)
  Hs <- sum(hs)
  for (j in seq(1, NL)) {
    f[j + 1] <- f[j] + hs[j] * (Hs - zm[j]) / VSm[j]
    A[j] <- (VSm[j] * (f[j + 1] - f[j]))^2 / hs[j]
    B[j] <- hs[j] * (f[j + 1] + f[j])^2
  }
  ws2 <- 4 * sum(A) / sum(B)
  Ts <- (2 * pi / sqrt(ws2)) |> round(3)
  Ts
}
