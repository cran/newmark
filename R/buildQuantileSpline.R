#' Quantile function of log spectral acceleration, monotone between its knots
#'
#' Builds, from the tabulated quantiles of one spectral ordinate, the
#' function \eqn{Q(u) = \ln S_a(u)} that [sampleSpectra()] uses to map the
#' uniforms of its Gaussian copula to accelerations when a spectrum is
#' given by quantiles.
#'
#' @section Mathematics:
#' The quantile rows, sorted by probability, give the nodes
#' \eqn{(p_i, \ln S_{a,i})}, \eqn{0 < p_i < 1}, with finite logarithms
#' because every \eqn{S_{a,i}} is finite and positive; a row outside these
#' conditions stops the function. The ordinates are replaced by their running
#' maximum, so that no node descends, and when two consecutive ordinates
#' are then equal every ordinate is lifted by \eqn{i \times 10^{-8}}. The
#' function is the cubic spline through these nodes with the monotonicity
#' filter of Hyman (1983) (`stats::splinefun(method = "hyman")`). Its
#' argument is clamped to \eqn{[0, 1]}; between 0 and the first tabulated
#' probability, and between the last one and 1, the end pieces of the
#' spline extrapolate. The monotonicity filter constrains the spline between
#' the tabulated knots only: an end piece can descend (observed on the tail
#' fixture of the site-response tests above \eqn{u = 0.9916}, on no slice
#' the products sample). A change of the end pieces is a separate
#' scientific candidate.
#'
#' @param quantiles A table with the probability levels `p` (numbers
#'   strictly between 0 and 1, or text that converts to them) and the finite
#'   positive spectral accelerations `Sa` of one ordinate. Rows with `p`
#'   equal to `"mean"` are ignored; a missing `p` is rejected; at
#'   least three quantile rows are required. Any other column must hold a
#'   single value, so that the rows belong to one spectrum slice.
#'   Duplicate levels signal a warning and are collapsed by the spline.
#' @return A function of probabilities \eqn{u} returning \eqn{\ln S_a}
#'   (natural logarithm of the acceleration, in the units of `Sa`).
#' @export
#' @examples
#' Q <- buildQuantileSpline(data.frame(p = c(0.16, 0.50, 0.84), Sa = c(0.3, 0.5, 0.8)))
#' exp(Q(c(0.16, 0.5, 0.84)))
buildQuantileSpline <- function(quantiles) {
  if (!is.data.frame(quantiles) || !all(c("p", "Sa") %in% names(quantiles))) {
    stop("quantiles must be a table with columns p and Sa")
  }
  DT <- as.data.table(quantiles)[p != "mean"]
  # Labels that do not convert are rejected just below; the coercion warning
  # would only repeat that diagnostic. The row filter drops a missing p, so
  # it is detected on the input.
  Level <- suppressWarnings(as.numeric(DT$p))
  if (anyNA(quantiles$p) || any(!is.finite(Level) | Level <= 0 | Level >= 1)) {
    stop("quantiles must have probabilities p strictly between 0 and 1")
  }
  if (!is.numeric(DT$Sa) || any(!is.finite(DT$Sa) | DT$Sa <= 0)) {
    stop("quantiles must have finite positive accelerations Sa")
  }
  if (nrow(DT) < 3L) {
    stop("buildQuantileSpline needs 3 or more quantile rows")
  }
  if (anyDuplicated(DT$p)) {
    warning("quantiles has duplicate p rows; the spline collapses them. ",
            "Select one spectrum slice before building the quantile function.", call. = FALSE)
  }
  COLS <- setdiff(names(DT), c("p", "Sa"))
  if (length(COLS)) {
    COLS <- COLS[vapply(DT[, COLS, with = FALSE], function(x) uniqueN(x) > 1L, logical(1L))]
    if (length(COLS)) stop("quantiles has columns with several values: ", paste(COLS, collapse = ", "))
  }
  Probability <- sort(Level)
  LnSa <- cummax(log(DT$Sa[order(Level)]))
  if (any(diff(LnSa) <= 0)) {
    LnSa <- LnSa + seq_along(LnSa) * 1e-8
  }
  # The only warning splinefun raises here is the collapse of duplicate
  # levels, already signalled above with its remedy.
  Spline <- suppressWarnings(stats::splinefun(x = Probability, y = LnSa, method = "hyman"))
  function(u) Spline(pmin(pmax(u, 0), 1))
}
