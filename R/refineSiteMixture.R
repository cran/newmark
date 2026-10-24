#' Site mixture at the first Gauss-Hermite order that agrees with the previous
#'
#' Evaluates [integrateSiteMixture()] at the increasing quadrature `orders`
#' and returns the site mixture of the first order whose quantiles and mean
#' agree with those of the previous order, together with the order and the
#' agreement reached, so that the numerical error of the result is visible.
#'
#' @section Mathematics:
#' At order \eqn{n_i} the report \eqn{R_i} is [summarizeMixture()] of the
#' site mixture at `levels` and the mean, at every period. With the relative
#' change \eqn{\delta(a, b) = \max |a - b| / \max(|b|, \epsilon)} over
#' corresponding entries (\eqn{\epsilon} the smallest positive double), the
#' change between \eqn{R_{i-1}} and \eqn{R_i} is the pair
#' \eqn{\Delta_i = (\delta \text{ over all entries}, \delta \text{ over the
#' means})}. The first change \eqn{\Delta_2} is the reference; order
#' \eqn{n_i}, \eqn{i \ge 3}, is accepted when
#' \eqn{\max \Delta_i \le} `tolerance` and each component of
#' \eqn{\Delta_i} does not exceed that of \eqn{\Delta_2} by more than
#' \eqn{10^{-12}}. No order accepted is an error. When `VS30` equals `Vref`
#' the rock mixture is returned without quadrature.
#'
#' @param mixture The rock mixture, see [integrateSiteMixture()].
#' @param VS30 One target velocity in m/s.
#' @param Vref One reference velocity in m/s.
#' @param model Site factor model of [evaluateSiteFactor()], `"ST20"` or
#'   `"ST17"`.
#' @param levels Probabilities in \eqn{(0, 1)} whose quantiles enter the
#'   comparison with the mean.
#' @param orders At least three increasing integers, each at least 2.
#' @param tolerance One positive relative change.
#' @return A list with `mixture` (the site mixture of the accepted order, as
#'   [integrateSiteMixture()] returns it), `order` (the accepted order) and
#'   `change` (the pair \eqn{\Delta} of that order, named `all` and `mean`);
#'   `order` and `change` are `NA` when `VS30` equals `Vref`.
#' @export
#' @examples
#' Mixture <- data.table::CJ(Tn = c(0, 0.2, 1), k = c("a", "b"))
#' Mixture[, `:=`(muLn = log(0.4) + 0.1 * (k == "b") - 0.5 * Tn, sigmaLn = 0.5, w = 0.5)]
#' refineSiteMixture(Mixture, VS30 = 360, Vref = 760, model = "ST20",
#'                   levels = c(0.16, 0.50, 0.84), orders = c(8, 16, 32, 64),
#'                   tolerance = 2e-3)[c("order", "change")]
refineSiteMixture <- function(mixture, VS30, Vref, model, levels, orders, tolerance) {
  .checkSiteCondition(VS30 = VS30, Vref = Vref, model = model)
  if (!is.numeric(orders) || length(orders) < 3L || any(!is.finite(orders)) ||
      any(orders != as.integer(orders)) || any(orders < 2) || any(diff(orders) <= 0)) {
    stop("orders must be at least three increasing integers >= 2")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("tolerance must be one positive number")
  }
  Rock <- .checkMixture(mixture)
  if (VS30 == Vref) return(list(mixture = Rock, order = NA_integer_, change = c(all = NA_real_, mean = NA_real_)))
  Previous <- NULL
  First <- NULL
  for (n in as.integer(orders)) {
    Mixture <- integrateSiteMixture(mixture = Rock, VS30 = VS30, Vref = Vref, model = model, order = n)
    Report <- summarizeMixture(Mixture, levels = levels)
    setorder(Report, Tn, level)
    if (!is.null(Previous)) {
      AUX <- merge(Previous, Report, by = c("Tn", "level"), suffixes = c(".x", ".y"))
      Change <- c(all = .relativeChange(AUX$Sa.x, AUX$Sa.y),
                  mean = .relativeChange(AUX[level == "mean", Sa.x], AUX[level == "mean", Sa.y]))
      if (!is.null(First) && max(Change) <= tolerance && all(Change <= First + 1e-12)) {
        return(list(mixture = Mixture, order = n, change = Change))
      }
      if (is.null(First)) First <- Change
    }
    Previous <- Report
  }
  stop("the site-mixture quadrature did not converge for VS30=", VS30,
       "; maximum relative change=", signif(max(Change), 4))
}

# Largest relative change of `current` against `previous`, entry by entry.
.relativeChange <- function(previous, current) {
  max(abs(previous - current) / pmax(abs(current), .Machine$double.xmin))
}
