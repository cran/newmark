#' Weighted mean and quantiles of pooled yield coefficients
#'
#' Pools the yield-coefficient realizations of one or more displacement
#' models with their epistemic weights and reduces the pool, at every
#' admissible displacement, to its weighted mean, its weighted quantiles and
#' the weighted share of roots that lie beyond the yield-coefficient grid.
#'
#' @section Mathematics:
#' The positive weights \eqn{W_m} of the pooled models are normalized,
#' \eqn{w_m = W_m / \sum_l W_l}. Every model contributes the same
#' realizations \eqn{s = 1, \dots, N_s}, so at one \eqn{D_a} the pool holds
#' the \eqn{N_m N_s} roots \eqn{k_{m,s}}, each weighted by its model's
#' \eqn{w_m}. The mean is
#' \deqn{\bar k_{max} = \frac{\sum_{m,s} w_m k_{m,s}}{\sum_{m,s} w_m}.}
#' For the quantile of level \eqn{p}, the pooled values are sorted and equal
#' values merged with their summed weights into \eqn{x_1 < \dots < x_K} with
#' cumulative weights \eqn{C_1 < \dots < C_K = W}; with
#' \eqn{h = 1 + (W - 1) p}, \eqn{a = \max(\lfloor h \rfloor, 1)} and
#' \eqn{b = \min(a + 1, W)}, the quantile is
#' \eqn{(1 - f)\, x(a) + f\, x(b)}, \eqn{f = h - \lfloor h \rfloor}, where
#' \eqn{x(c)} is the smallest \eqn{x_k} with \eqn{C_k \ge c} (or \eqn{x_K}
#' when none). This is the seventh sample quantile definition of Hyndman and
#' Fan (1996) extended to cumulative weights, the estimator `newmark` 1.7.1
#' obtained from `Hmisc::wtd.quantile()`; for one model of unit weight it is
#' the type-7 sample quantile of its \eqn{N_s} roots. The share beyond the
#' grid is \eqn{\sum_m w_m \bar b_m}, with \eqn{\bar b_m} the fraction of the
#' realizations of model \eqn{m} whose root lies beyond the grid.
#'
#' A root of 0 (a target above the ceiling of its realization, see
#' [solveKmax()]) stays in the pool as a realization of null demand.
#'
#' @param roots The roots of [solveKmax()] for one or more models, bound by
#'   rows: `IDn`, `realization`, `Da`, `kmax` and `beyond`, one row per
#'   model, realization and `Da`; every pooled model must hold the same
#'   realizations at every `Da`.
#' @param weights Named non-negative weights by model; the models with
#'   positive weight form the pool and must all be present in `roots` with
#'   the same realizations.
#' @param levels Unique probabilities in (0, 1), distinct at two decimals.
#' @return A data.table with one row per `Da` and level: `Da` (cm), `level`
#'   (`"mean"`, then each probability formatted `"%.2f"`), `kmax` (g), the
#'   weighted mean or quantile of the pooled realizations at that `Da`,
#'   `beyond`, the weighted share of the pool beyond the grid, and `n`, the
#'   number of pooled roots (models times realizations), repeated on every
#'   level of one `Da`.
#' @references Hyndman, R. J. and Fan, Y. (1996). Sample quantiles in
#'   statistical packages. The American Statistician, 50(4), 361-365.
#'   \doi{10.1080/00031305.1996.10473566}
#' @export
#' @examples
#' ky <- c(0.05, 0.10, 0.20, 0.30, 0.40)
#' Roots <- rbind(
#'   solveKmax(ky, matrix(0.10 / ky, 20, 5, byrow = TRUE), PGA = NULL, model = "AM88", Da = 1),
#'   solveKmax(ky, matrix(0.30 / ky, 20, 5, byrow = TRUE), PGA = NULL, model = "BT07", Da = 1))
#' summarizeKmax(Roots, weights = c(AM88 = 0.5, BT07 = 0.5), levels = c(0.05, 0.50, 0.95))
summarizeKmax <- function(roots, weights, levels) {
  .poolKmax(roots, weights = weights, levels = levels, estimator = .weightedQuantile)
}

# The pooled summary with its quantile estimator as an argument: the
# exported summary passes .weightedQuantile, the compatibility route of
# buildKmaxTable() the text-merging variant of newmark 1.7.1.
.poolKmax <- function(roots, weights, levels, estimator) {
  if (!is.numeric(levels) || !length(levels) || any(!is.finite(levels)) || any(levels <= 0 | levels >= 1) || anyDuplicated(levels)) {
    stop("levels must contain unique finite probabilities in (0, 1)")
  }
  Labels <- sprintf("%.2f", levels)
  if (anyDuplicated(Labels)) stop("levels values collide after two-decimal formatting")
  if (!is.numeric(weights) || !length(weights) || is.null(names(weights)) || anyNA(names(weights)) ||
      any(!nzchar(names(weights))) || anyDuplicated(names(weights)) || any(!is.finite(weights)) ||
      any(weights < 0) || !any(weights > 0)) {
    stop("weights must be finite, non-negative, uniquely named, and include one positive value")
  }
  Required <- c("IDn", "realization", "Da", "kmax", "beyond")
  if (!is.data.frame(roots) || !all(Required %in% names(roots))) {
    stop("roots must contain columns ", paste(Required, collapse = ", "))
  }
  Weights <- weights[weights > 0]
  Weights <- Weights / Reduce(`+`, Weights)
  Missing <- setdiff(names(Weights), unique(roots$IDn))
  if (length(Missing)) stop("kmax roots are missing positive-weight models: ", paste(Missing, collapse = ", "))
  DT <- as.data.table(roots)[IDn %in% names(Weights)]
  Models <- sort(names(Weights))
  Weights <- Weights[Models]
  Samples <- lapply(Models, function(m) sort(unique(DT$realization[DT$IDn == m])))
  if (!all(vapply(Samples, identical, logical(1L), y = Samples[[1L]]))) {
    stop("kmax roots require complete shared realization coverage")
  }
  if (anyDuplicated(DT, by = c("IDn", "realization", "Da"))) {
    stop("kmax roots must hold one row per model, realization and Da")
  }
  # The shared coverage must hold within every Da, not only over the whole table.
  Cover <- DT[, .(complete = length(unique(.SD$IDn)) == length(Models) &&
                    all(vapply(split(.SD$realization, .SD$IDn),
                               function(s) identical(sort(unique(s)), Samples[[1L]]), logical(1L)))),
              by = "Da", .SDcols = c("IDn", "realization")]
  if (!all(Cover$complete)) stop("kmax roots require complete shared realization coverage within every Da")
  # Bound outside the grouped call, where a column named levels would shadow it.
  Probs <- levels
  DT[, {
    LIST <- lapply(Models, function(m) kmax[IDn == m])
    x <- unlist(LIST)
    w <- rep(as.numeric(Weights), each = length(LIST[[1L]]))
    data.table(level = c("mean", Labels), kmax = c(Reduce(`+`, w * x) / Reduce(`+`, w), estimator(x, w, Probs)),
               beyond = Reduce(`+`, vapply(seq_along(Models), function(i) Weights[[i]] * mean(beyond[IDn == Models[[i]]]),
                                   numeric(1L))),
               n = length(x))
  }, by = "Da"]
}
