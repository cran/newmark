#' Peak ground velocity from peak ground acceleration
#'
#' Evaluates the in-house rock regression of `newmark` 1.7.1 that gives each
#' sampled peak ground acceleration its peak ground velocity, the velocity
#' the pulse equations of `BM19` in [evaluateFlexibleDisplacement()] read.
#'
#' @section Mathematics:
#' \deqn{PGV = 100\, PGA^{1.0529}\, e^{0.1241}}
#' with PGA in g and PGV in cm/s, a median relation without dispersion,
#' fitted on the rock records of Verri Kozlowski (2023), as recorded in
#' `newmark` 1.7.1. It is evaluated realization by realization, so a sampled
#' PGV stays paired with its PGA.
#'
#' @param PGA Peak ground accelerations in g, finite and positive.
#' @return Peak ground velocities in cm/s, aligned with `PGA`.
#' @export
#' @examples
#' estimatePGV(PGA = c(0.1, 0.5))
estimatePGV <- function(PGA) {
  if (!is.numeric(PGA) || !length(PGA) || any(!is.finite(PGA)) || any(PGA <= 0)) {
    stop("PGA must be finite positive accelerations")
  }
  (PGA^1.0529) * exp(0.1241) * 100
}

#' @name estimatePGV
#' @references Verri Kozlowski, A. (2023). Modelos Parametricos para la
#'   Seleccion de Sismos de Diseno en Infraestructura Minera. Seminario de
#'   Vinculacion y Transferencia (SeVyT) 2023.
NULL
