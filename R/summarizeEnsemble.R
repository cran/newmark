#' Weighted statistics of a displacement model ensemble
#'
#' Pools the realizations of the models named in `weights` at each yield
#' acceleration and reduces them to weighted quantiles and a weighted mean,
#' the ensemble displacement of those models. Every realization of a model
#' carries the normalized weight of its model. The statistics are marginal
#' over yield acceleration, as in [summarizeDisplacements()].
#'
#' The population of each row is the `n` pooled realizations at that yield
#' acceleration. The pool inherits the conditioning of its members: the
#' realizations of the flexible models of [evaluateDisplacements()] are
#' conditional on a nonzero displacement ([evaluateFlexibleDisplacement()]),
#' so the probability of zero or negligible displacement of those published
#' models is not represented in the ensemble quantiles and mean.
#'
#' @section Mathematics:
#' With \eqn{\omega_m = w_m / \sum_j w_j} over the named models, the pooled
#' sample at one \eqn{k_y} holds \eqn{D_{n,s}^{(m)}} with weight
#' \eqn{\omega_m}. The weighted mean is
#' \eqn{\sum \omega_m D_{n,s}^{(m)} / \sum \omega_m} and the weighted
#' quantile of probability \eqn{q} is the seventh sample quantile definition
#' of Hyndman and Fan (1996), the type 7 of `quantile()`, extended to
#' cumulative weights: with the distinct sorted values \eqn{x_i}, their summed
#' weights and the total weight \eqn{W}, the order statistics of rank
#' \eqn{\lfloor h \rfloor} and \eqn{\lfloor h \rfloor + 1},
#' \eqn{h = 1 + (W - 1) q}, are the smallest \eqn{x_i} whose cumulative
#' weight reaches each rank, interpolated linearly by the fractional part of
#' \eqn{h}. With equal model weights and \eqn{N} realizations per model,
#' \eqn{W = N}. Both statistics are computed here without Hmisc and at full
#' precision: equal values are merged exactly, with the estimator of
#' [summarizeKmax()]. When the pooled sample holds equal values, Hmisc
#' merges values by their 15-significant-digit text, so every pooled value is
#' then read back from that text; the compatibility route of
#' [buildDnTables()] keeps that text merge, as `newmark` 1.7.1 did, and this
#' general summary does not (design decision).
#'
#' @param draws A table with `IDn`, `ky` and finite non-negative `Dn` (cm),
#'   as [evaluateDisplacements()] returns; rows of models not named in
#'   `weights` are ignored.
#' @param weights Named positive finite weights of the ensemble models;
#'   only their ratios matter.
#' @param levels Probabilities in (0, 1) and/or the literal `"mean"`.
#' @return A data.table with `ky` (g), `level`, `n` (pooled realizations)
#'   and `Dn` (cm), one row per yield acceleration and level, the yield
#'   accelerations in their order of appearance.
#' @references Hyndman, R. J. and Fan, Y. (1996). Sample quantiles in
#'   statistical packages. The American Statistician, 50(4), 361-365.
#'   \doi{10.1080/00031305.1996.10473566}
#' @export
#' @examples
#' Draws <- data.frame(IDn = rep(c("AM88", "BT07"), each = 4), ky = 0.1,
#'                     Dn = c(0.5, 1, 2, 4, 1, 2, 3, 5))
#' summarizeEnsemble(Draws, weights = c(AM88 = 1, BT07 = 1), levels = c(0.5, "mean"))
summarizeEnsemble <- function(draws, weights, levels) {
  .poolDisplacements(draws, weights = weights, levels = levels, estimator = .weightedQuantile)
}

# The ensemble summary with its quantile estimator as an argument: the
# exported summary passes .weightedQuantile, the compatibility route of
# buildDnTables() the text-merging variant of newmark 1.7.1.
.poolDisplacements <- function(draws, weights, levels, estimator) {
  if (!is.data.frame(draws) || !all(c("IDn", "ky", "Dn") %in% names(draws)) ||
      !is.numeric(draws$Dn) || any(!is.finite(draws$Dn)) || any(draws$Dn < 0)) {
    stop("draws must contain IDn, ky and finite non-negative Dn")
  }
  if (!is.numeric(weights) || !length(weights) || is.null(names(weights)) || anyNA(names(weights)) ||
      anyDuplicated(names(weights)) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must be named positive finite model weights")
  }
  Level <- .displacementLevels(levels)
  Weights <- weights / Reduce(`+`, weights)
  DT <- as.data.table(draws)[IDn %in% names(Weights)]
  if (!nrow(DT)) stop("draws has no realization of the weighted models")
  DT[, {
    W <- Weights[IDn]
    .(level = Level$label, n = .N,
      Dn = .levelValues(Level,
                        quantile = function(q) estimator(Dn, W, q),
                        mean = function() Reduce(`+`, W * Dn) / Reduce(`+`, W)))
  }, by = ky]
}
