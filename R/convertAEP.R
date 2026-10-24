#' Probability of exceedance and return period of an annual rate
#'
#' @section Mathematics:
#' Under Poisson occurrence an annual rate \eqn{AEP} gives, over an
#' investigation time \eqn{IT}, the probability of at least one exceedance
#' and the return period
#' \deqn{POE = 1 - e^{-IT \cdot AEP}, \qquad TR = 1/AEP .}
#' The conversion is exact; it carries no dispersion.
#'
#' @param AEP Annual rates of exceedance in 1/year.
#' @param IT Investigation time in years, finite and positive: one value or
#'   one per `AEP`.
#' @return A data.table with `POE` (probability within `IT`) and `TR`
#'   (years), one row per `AEP`.
#' @export
#' @examples
#' convertAEP(AEP = c(1 / 475, 1 / 2475), IT = 50)
convertAEP <- function(AEP, IT) {
  AEP <- as.numeric(AEP)
  data.table(POE = -expm1(-.recycleIT(IT, length(AEP)) * AEP), TR = 1 / AEP)
}

# One finite positive investigation time per value.
.recycleIT <- function(IT, n) {
  IT <- as.numeric(IT)
  if (!(length(IT) %in% c(1L, n))) {
    stop("Investigation time must be scalar or match the data length.")
  }
  if (any(!is.finite(IT) | IT <= 0)) {
    stop("Investigation time must be finite and > 0.")
  }
  rep_len(IT, n)
}
