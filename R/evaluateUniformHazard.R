#' Uniform-hazard ordinates of hazard curves at given return periods
#'
#' Evaluates every hazard curve of a table at the requested return periods,
#' which gives the ordinates of the uniform-hazard spectra of those curves.
#'
#' @section Mathematics:
#' For one curve with points \eqn{(TR_i, S_{a,i})}, only points with finite
#' positive \eqn{TR} and \eqn{S_a} are used, and points sharing a \eqn{TR}
#' keep the largest \eqn{S_a} (a message reports it when their ordinates
#' differ beyond 12 significant digits). The ordinate at each requested
#' \eqn{TR_o} is the log-log linear interpolation
#' \deqn{S_a(TR_o) = \exp L\left(\ln TR_o;\ \ln TR_i, \ln S_{a,i}\right),}
#' with \eqn{L} piecewise linear. A requested period may lie outside the
#' curve's range by at most \eqn{\max(0.5, 10^{-4}\max(TR_i, TR_o))} years,
#' where the end ordinate is taken; beyond that margin the call stops rather
#' than extrapolate. The result evaluates the given curve (for OpenQuake, a
#' statistic curve of the logic tree) at the return period; the ordinate of
#' a mean curve is not, in general, the mean of the ordinates.
#'
#' @param curves A table with columns `TR` (years), `Sa` (g) and the columns
#'   named in `by`.
#' @param TR Finite positive return periods in years.
#' @param by Names of the columns of `curves` that identify one curve.
#' @return A data.table with the `by` columns, `TR` and `Sa` (g), one row
#'   per curve and requested period, curves in order of first appearance
#'   and periods in the order given. A curve with fewer than two usable
#'   points stops the call.
#' @export
#' @examples
#' Curves <- data.frame(ID = "A", TR = c(50, 500, 5000), Sa = c(0.05, 0.2, 0.6))
#' evaluateUniformHazard(Curves, TR = c(100, 2475), by = "ID")
evaluateUniformHazard <- function(curves, TR, by) {
  if (!is.character(by) || anyDuplicated(by) || any(by %in% c("TR", "Sa"))) {
    stop("by must name distinct columns other than TR and Sa")
  }
  if (!is.data.frame(curves) || !all(c("TR", "Sa", by) %in% names(curves))) {
    stop("curves must be a table with columns TR, Sa and ", paste(by, collapse = ", "))
  }
  if (!is.numeric(TR) || !length(TR) || any(!is.finite(TR) | TR <= 0)) {
    stop("TR must contain finite positive return periods")
  }
  Target <- TR
  as.data.table(curves)[, .evaluateCurve(TR, Sa, Target, .BY), by = by]
}

.evaluateCurve <- function(TR, Sa, Target, curve) {
  Label <- paste(vapply(names(curve), function(x) paste0(x, "=", as.character(curve[[x]])), character(1)),
                 collapse = ", ")
  DT <- data.table(TR = as.numeric(TR), Sa = as.numeric(Sa))[is.finite(TR) & TR > 0 & is.finite(Sa) & Sa > 0]
  if (nrow(DT) < 2L) {
    stop("Need at least two positive (TR, Sa) points for ", Label, ".")
  }
  setorder(DT, TR)
  if (anyDuplicated(DT$TR)) {
    DT <- DT[, .(Sa = max(Sa), nSa = uniqueN(signif(Sa, 12))), by = TR]
    if (any(DT$nSa > 1L)) {
      message(">> Collapsing duplicate TR nodes by keeping max Sa in ", Label, ".")
    }
    DT <- DT[, .(TR, Sa)]
    setorder(DT, TR)
  }
  if (nrow(DT) < 2L) {
    stop("Need at least two unique TR nodes for ", Label, ".")
  }
  Lower <- min(DT$TR)
  Upper <- max(DT$TR)
  Margin <- max(0.5, 1e-4 * max(c(abs(DT$TR), abs(Target))))
  if (min(Target) < Lower - Margin || max(Target) > Upper + Margin) {
    stop("Requested TR grid [", signif(min(Target), 8), ", ", signif(max(Target), 8),
         "] exceeds imported curve range [", signif(Lower, 8), ", ", signif(Upper, 8),
         "] for ", Label, ". Expand the imported POE grid or narrow the requested TR instead of extrapolating.")
  }
  LogSa <- stats::approx(x = log(DT$TR), y = log(DT$Sa), xout = log(pmin(pmax(Target, Lower), Upper)),
                         method = "linear", ties = "ordered")$y
  if (any(!is.finite(LogSa))) {
    stop("Interpolation failed for ", Label, ".")
  }
  data.table(TR = Target, Sa = exp(LogSa))
}
