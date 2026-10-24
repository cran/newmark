#' Quantiles and mean of a spectral mixture
#'
#' Reduces a spectrum represented as a finite lognormal mixture at every
#' oscillator period, such as the rock spectra of the hazard products or the result of
#' [integrateSiteMixture()], to the spectral accelerations at the requested
#' probabilities and to the mean, period by period.
#'
#' @section Mathematics:
#' At each period the weights are normalized, \eqn{\tilde w_k = w_k / \sum
#' w}, and the mixture distribution function is
#' \eqn{F(x) = \sum_k \tilde w_k \Phi((x - \mu_k)/\sigma_k)} for
#' \eqn{x = \ln S_a}, a component with \eqn{\sigma_k = 0} contributing the
#' step \eqn{1\{x \ge \mu_k\}}. The quantile at probability \eqn{p} is
#' \eqn{\exp(x_p)} with \eqn{x_p = \inf\{x : F(x) \ge p\}}, the smallest
#' value whose distribution function reaches \eqn{p} (generalized inverse).
#' When every dispersion is positive, \eqn{F} is continuous and strictly
#' increasing and \eqn{x_p} is the root of \eqn{F(x) = p}, found by
#' [stats::uniroot()] to \eqn{10^{-12}} on
#' \eqn{[\min(\mu - 12 s), \max(\mu + 12 s)]}, \eqn{s = \max(\sigma, 1)};
#' otherwise \eqn{F} jumps at the atoms, \eqn{F(x) = p} may have no root,
#' and \eqn{x_p} is the upper end of 100 bisections of that interval, each
#' keeping \eqn{F \ge p} at the upper end, so that a probability inside a
#' jump returns the atom. The mean is
#' \eqn{\sum_k \tilde w_k e^{\mu_k + \sigma_k^2/2}}. These are marginal
#' statistics of the ordinate at each period; they do not describe one
#' spectrum across periods.
#'
#' @param mixture A mixture as [integrateSiteMixture()] takes and returns
#'   it, validated by the same rules: one row per `(Tn, k)` with columns
#'   `Tn` (s), `k` (component label), `muLn`, `sigmaLn` (finite, `sigmaLn`
#'   non-negative) and `w` (positive weights summing to one at every period),
#'   including the PGA rows (`Tn == 0`).
#' @param levels Probabilities in \eqn{(0, 1)}.
#' @return A data.table ordered by `Tn` with `Tn`, `level` (the probability
#'   as text, then `"mean"`) and `Sa` (g), one row per period and level,
#'   unrounded.
#' @export
#' @examples
#' Mixture <- data.table::CJ(Tn = c(0, 0.2, 1), k = c("a", "b"))
#' Mixture[, `:=`(muLn = log(0.4) + 0.1 * (k == "b") - 0.5 * Tn, sigmaLn = 0.5, w = 0.5)]
#' summarizeMixture(Mixture, levels = c(0.16, 0.50, 0.84))
summarizeMixture <- function(mixture, levels) {
  Mixture <- .checkMixture(mixture)
  if (!is.numeric(levels) || !length(levels) || any(!is.finite(levels)) || any(levels <= 0 | levels >= 1)) {
    stop("levels must be probabilities in (0, 1)")
  }
  Labels <- c(as.character(levels), "mean")
  # The checked copy holds only the mixture columns, so a column of the
  # caller's table cannot shadow the levels in the grouped evaluation.
  Mixture[, {
    Weights <- w / Reduce(`+`, w)
    .(level = Labels,
      Sa = c(vapply(levels, function(x) .solveMixtureQuantile(muLn, sigmaLn, Weights, x), numeric(1L)),
             Reduce(`+`, Weights * exp(muLn + sigmaLn^2 / 2))))
  }, by = Tn]
}

# Acceleration quantile of a lognormal mixture whose components may be
# degenerate: root search on the mixture CDF when every dispersion is
# positive, bisection otherwise; its sums avoid long double.
.solveMixtureQuantile <- function(muLn, sigmaLn, w, p) {
  Weights <- w / Reduce(`+`, w)
  CDF <- function(x) {
    Terms <- ifelse(sigmaLn > 0, stats::pnorm((x - muLn) / sigmaLn), as.numeric(x >= muLn))
    Reduce(`+`, Weights * Terms)
  }
  Scale <- pmax(sigmaLn, 1)
  Lower <- min(muLn - 12 * Scale)
  Upper <- max(muLn + 12 * Scale)
  if (all(sigmaLn > 0)) {
    return(exp(stats::uniroot(function(x) CDF(x) - p, lower = Lower, upper = Upper, tol = 1e-12)$root))
  }
  for (i in seq_len(100L)) {
    Middle <- (Lower + Upper) / 2
    if (CDF(Middle) >= p) Upper <- Middle else Lower <- Middle
  }
  exp(Upper)
}
