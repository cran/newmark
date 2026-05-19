###############################################################################
# getKyLimits.R — Calibration ky range per Newmark displacement model
###############################################################################

#' Get calibrated ky range for each displacement model
#'
#' Returns the ky bounds within which each model's regression is calibrated.
#' For rigid models (AM88, SR08) the bounds are relative to PGA.
#' For flexible models the bounds are absolute.
#'
#' @param PGA numeric scalar, peak ground acceleration (g). Required for
#'   rigid model limits.
#'
#' @return data.table with columns IDn, ky.lo, ky.hi
#' @examples
#' getKyLimits(PGA = 0.5)
#' @export
getKyLimits <- function(PGA) {
  data.table::data.table(
    IDn   = c("AM88",  "JB07", "SR08", "BT07", "BM17", "BM19"),
    ky.lo = c(0.001,   0.05,   0.05,   0.02,   0.01,   0.01),
    ky.hi = c(PGA,     0.40,   0.30,   0.40,   0.80,   0.80)
  )
}
