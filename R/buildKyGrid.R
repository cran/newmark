#' Yield-acceleration grid of a displacement scenario
#'
#' Builds the log-spaced yield accelerations at which the displacement
#' curves of one scenario are evaluated, from the calibration ranges of the
#' displacement models and the mean peak ground acceleration of the
#' scenario.
#'
#' @section Mathematics:
#' For \eqn{n > 1},
#' \deqn{k_{y,i} = \exp\left(\ln 0.01 + \frac{i - 1}{n - 1} \left[\ln \max(\overline{PGA}, 0.80) - \ln 0.01\right]\right), \quad i = 1, \dots, n}
#' in g. For \eqn{n = 1} the grid is the lower bound alone,
#' \eqn{k_{y,1} = \exp(\ln 0.01)}, whatever the PGA. The lower bound
#' 0.01 g is the lowest calibration point of `BM17`;
#' the upper bound is the mean peak ground acceleration of the scenario,
#' beyond which a rigid block does not slide, extended to the 0.80 g
#' calibration ceiling of `BM17` for low-intensity scenarios. The
#' displacement products take \eqn{\overline{PGA}} as the mixture mean
#' \eqn{\sum_k w_k e^{\mu_k + \sigma_k^2/2}} at zero period of an
#' analytic spectrum, or the `"mean"` row at the smallest period of a
#' quantile spectrum.
#'
#' @param PGA One finite positive mean peak ground acceleration of the
#'   scenario in g.
#' @param n Number of grid points, one positive integer (the displacement
#'   products use 30).
#' @return A numeric vector of `n` increasing yield accelerations in g.
#' @export
#' @examples
#' buildKyGrid(PGA = 0.5, n = 30)
buildKyGrid <- function(PGA, n) {
  if (!is.numeric(PGA) || length(PGA) != 1L || !is.finite(PGA) || PGA <= 0) {
    stop("PGA must be one finite positive mean acceleration")
  }
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("n must be one positive integer")
  }
  exp(seq(from = log(0.01), to = log(max(PGA, 0.80)), length.out = n))
}
