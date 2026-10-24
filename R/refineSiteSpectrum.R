#' Site spectrum at the first quadrature order that agrees with the previous
#'
#' Evaluates [integrateSiteSpectrum()] at the increasing Gauss-Legendre
#' `orders` and returns the site spectrum of the first order whose ordinates
#' agree with those of the previous order, together with the order and the
#' agreement reached, so that the numerical error of the result is visible.
#'
#' @section Mathematics:
#' With \eqn{S_i} the site ordinates (`SaF`, every period and level) at
#' order \eqn{n_i} and the relative change
#' \eqn{\delta_i = \max |S_{i-1} - S_i| / \max(|S_i|, \epsilon)}
#' (\eqn{\epsilon} the smallest positive double), order \eqn{n_i},
#' \eqn{i \ge 2}, is accepted when \eqn{\delta_i} is finite and at most
#' `tolerance`. No order accepted is an error.
#'
#' @param spectrum The rock quantile spectrum, see [integrateSiteSpectrum()].
#' @param VS30 One target velocity in m/s.
#' @param Vref One reference velocity in m/s.
#' @param model Site factor model of [evaluateSiteFactor()], `"ST20"` or
#'   `"ST17"`.
#' @param orders At least two increasing integers, each at least 4.
#' @param tolerance One positive relative change.
#' @return A list with `spectrum` (the site spectrum of the accepted order,
#'   as [integrateSiteSpectrum()] returns it), `order` (the accepted order)
#'   and `change` (\eqn{\delta} of that order).
#' @export
#' @examples
#' P <- c("0.05", "0.10", "0.16", "0.50", "0.84", "0.90", "0.95", "mean")
#' Factor <- c(0.45, 0.60, 0.72, 0.95, 1.30, 1.48, 1.78, 1.00)
#' Spectrum <- data.table::data.table(Tn = rep(c(0, 0.2, 1), each = 8), p = rep(P, 3),
#'                                    Sa = rep(c(0.30, 0.45, 0.24), each = 8) * rep(Factor, 3))
#' refineSiteSpectrum(Spectrum, VS30 = 270, Vref = 760, model = "ST20",
#'                    orders = c(32, 64, 128), tolerance = 2e-3)[c("order", "change")]
refineSiteSpectrum <- function(spectrum, VS30, Vref, model, orders, tolerance) {
  if (!is.numeric(orders) || length(orders) < 2L || any(!is.finite(orders)) ||
      any(orders != as.integer(orders)) || any(orders < 4) || any(diff(orders) <= 0)) {
    stop("orders must be at least two increasing integers >= 4")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be one positive number")
  }
  Orders <- as.integer(orders)
  Previous <- integrateSiteSpectrum(spectrum = spectrum, VS30 = VS30, Vref = Vref, model = model,
                                    order = Orders[[1L]])
  for (n in Orders[-1L]) {
    Current <- integrateSiteSpectrum(spectrum = spectrum, VS30 = VS30, Vref = Vref, model = model, order = n)
    Change <- .relativeChange(Previous$SaF, Current$SaF)
    if (is.finite(Change) && Change <= tolerance) return(list(spectrum = Current, order = n, change = Change))
    Previous <- Current
  }
  stop("the site-spectrum quadrature did not converge by order ", Orders[[length(Orders)]],
       " (relative change ", format(Change, digits = 8L), ")")
}
