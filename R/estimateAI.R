#' Arias intensity from peak ground acceleration
#'
#' Evaluates the in-house rock regression of `newmark` 1.7.1 that gives each
#' sampled peak ground acceleration its Arias intensity, the intensity the
#' rigid models `JB07` and `SR08` of [evaluateRigidDisplacement()] read.
#'
#' @section Mathematics:
#' \deqn{I_a = PGA^{1.9228}\, e^{2.6109}}
#' with PGA in g and \eqn{I_a} in m/s, a median relation without
#' dispersion, fitted on rock records (SeVyT 2023 dataset, as recorded in
#' `newmark` 1.7.1). It is evaluated realization by realization, so a sampled
#' intensity stays paired with its PGA.
#'
#' @param PGA Peak ground accelerations in g, finite and positive.
#' @return Arias intensities in m/s, aligned with `PGA`.
#' @export
#' @examples
#' estimateAI(PGA = c(0.1, 0.5))
estimateAI <- function(PGA) {
  if (!is.numeric(PGA) || !length(PGA) || any(!is.finite(PGA)) || any(PGA <= 0)) {
    stop("PGA must be finite positive accelerations")
  }
  (PGA^1.9228) * exp(2.6109)
}
