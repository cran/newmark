#' Interpolate Sa(p) at arbitrary period using log-log interpolation
#'
#' @param uhs data.table with columns Tn, p, Sa.
#' @param Tn numeric scalar - target period for interpolation.
#'
#' @return data.table with interpolated Tn, p, Sa
#' @examples
#' uhs <- data.table::data.table(
#'   Tn = c(0.2, 0.5, 1.0, 0.2, 0.5, 1.0),
#'   p  = c("0.16", "0.16", "0.16", "0.84", "0.84", "0.84"),
#'   Sa = c(0.8, 0.6, 0.3, 1.2, 0.9, 0.5)
#' )
#' interpolateSaTable(uhs, Tn = 0.3)
#' @export
interpolateSaTable <- function(uhs, Tn) {
  stopifnot(all(c("Tn", "p", "Sa") %in% names(uhs)))
  stopifnot(is.numeric(Tn), length(Tn) == 1L)

  if (anyDuplicated(uhs[, .(Tn, p)])) {
    warning("interpolateSaTable(): input has duplicate (Tn, p) rows; ",
            "values will be collapsed via ties = 'mean'. Disaggregate ",
            "by (siteID, ID, Vs30, TR) before calling to silence this.",
            call. = FALSE)
  }

  p <- sort(unique(uhs[p != "mean"]$p))
  if (length(p) == 0L) {
    stop("interpolateSaTable(): missing quantiles (p != 'mean')")
  }

  SaTable <- data.table::rbindlist(lapply(p, function(pj) {
    rows <- uhs[p == pj]
    if (nrow(rows) < 2L) {
      stop("interpolateSaTable(): need >=2 Tn values for p = ", pj)
    }

    y <- suppressWarnings(stats::approx(
      x = log(rows$Tn + 1e-8),
      y = log(rows$Sa),
      xout = log(Tn),
      rule = 2
    ))$y

    data.table::data.table(Tn = Tn, p = pj, Sa = exp(y))
  }))

  MeanRow <- uhs[p == "mean"]
  if (nrow(MeanRow) >= 2L) {
    YMean <- suppressWarnings(stats::approx(
      x = log(MeanRow$Tn + 1e-8),
      y = log(MeanRow$Sa),
      xout = log(Tn),
      rule = 2
    ))$y

    SaTable <- data.table::rbindlist(list(
      SaTable,
      data.table::data.table(Tn = Tn, p = "mean", Sa = exp(YMean))
    ))
  } else if (nrow(MeanRow) == 1L) {
    SaTable <- data.table::rbindlist(list(
      SaTable,
      data.table::data.table(Tn = Tn, p = "mean", Sa = MeanRow$Sa)
    ))
  }

  data.table::setorder(SaTable, Tn, p)
  SaTable[]
}
