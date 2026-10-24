#' Calibration range of yield acceleration of each displacement model
#'
#' Tabulates the yield accelerations within which the regression of each
#' displacement model was calibrated. The rigid model `AM88` is bounded by
#' the peak ground acceleration of the scenario; the others have absolute
#' bounds. The displacement products evaluate every model on the common
#' grid of [buildKyGrid()] and do not restrict it to these ranges.
#'
#' @param PGA One finite positive peak ground acceleration in g, the upper
#'   bound of `AM88`.
#' @return A data.table with `IDn` and the bounds `kyMin` and `kyMax` in g.
#' @export
#' @examples
#' buildKyLimits(PGA = 0.5)
buildKyLimits <- function(PGA) {
  if (!is.numeric(PGA) || length(PGA) != 1L || !is.finite(PGA) || PGA <= 0) {
    stop("PGA must be one finite positive acceleration")
  }
  data.table(IDn = c("AM88", "JB07", "SR08", "BT07", "BM17", "BM19"),
             kyMin = c(0.001, 0.05, 0.05, 0.02, 0.01, 0.01),
             kyMax = c(PGA, 0.40, 0.30, 0.40, 0.80, 0.80))
}
