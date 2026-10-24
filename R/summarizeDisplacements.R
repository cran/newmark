#' Marginal statistics of displacement realizations by model
#'
#' Reduces the realizations of [evaluateDisplacements()] to the requested
#' quantiles and mean of the displacement of each model at each yield
#' acceleration. The statistics are marginal: a quantile row of one yield
#' acceleration and the row of the next do not belong to the same
#' realization. Keep the realizations when a later calculation needs paired
#' values.
#'
#' The population of each row is the `n` realizations given for that model
#' and yield acceleration. Realizations of the flexible models of
#' [evaluateDisplacements()] are conditional on a nonzero displacement
#' ([evaluateFlexibleDisplacement()]), so their quantiles and mean describe
#' the nonzero displacement and do not include the probability of zero or
#' negligible displacement of the published models.
#'
#' @section Mathematics:
#' For the \eqn{N} realizations \eqn{D_{n,s}} of one model at one
#' \eqn{k_y}, the quantile of probability \eqn{q} is type 7 of
#' [stats::quantile()] (linear interpolation between the order statistics
#' \eqn{\lfloor h \rfloor} and \eqn{\lceil h \rceil}, \eqn{h = 1 + (N - 1) q})
#' and the mean is \eqn{\frac{1}{N}\sum_s D_{n,s}}.
#'
#' @param draws A table with `IDn`, `ky` and finite non-negative `Dn` (cm),
#'   one row per realization, as [evaluateDisplacements()] returns.
#' @param levels Probabilities in (0, 1) and/or the literal `"mean"`, in the
#'   order the rows are wanted.
#' @return A data.table with `IDn`, `ky` (g), `level` (the probability as
#'   text, or `"mean"`), `n` (realizations summarized) and `Dn` (cm), one
#'   row per model, yield acceleration and level, the models and yield
#'   accelerations in their order of appearance in `draws`.
#' @export
#' @examples
#' Draws <- data.frame(IDn = "AM88", ky = 0.1, Dn = c(0.5, 1, 2, 4))
#' summarizeDisplacements(Draws, levels = c(0.16, 0.5, 0.84, "mean"))
summarizeDisplacements <- function(draws, levels) {
  if (!is.data.frame(draws) || !nrow(draws) || !all(c("IDn", "ky", "Dn") %in% names(draws)) ||
      !is.numeric(draws$Dn) || any(!is.finite(draws$Dn)) || any(draws$Dn < 0)) {
    stop("draws must contain rows with IDn, ky and finite non-negative Dn")
  }
  Level <- .displacementLevels(levels)
  as.data.table(draws)[, .(level = Level$label, n = .N,
                           Dn = .levelValues(Level, quantile = function(q) stats::quantile(Dn, probs = q, names = FALSE),
                                             mean = function() { s <- Reduce(`+`, Dn) / length(Dn); s + Reduce(`+`, Dn - s) / length(Dn) })),
                       by = .(IDn, ky)]
}

# Levels of a displacement summary: the labels, which of them is the mean
# and the probabilities of the others.
.displacementLevels <- function(levels) {
  if (!length(levels) || anyNA(levels)) stop("levels must contain probabilities and/or \"mean\"")
  Label <- as.character(levels)
  IsMean <- Label == "mean"
  # Tokens other than "mean" that do not convert are rejected just below; the
  # coercion warning would only repeat that diagnostic.
  Probability <- suppressWarnings(as.numeric(Label[!IsMean]))
  if (any(!is.finite(Probability) | Probability <= 0 | Probability >= 1)) {
    stop("levels must be probabilities in (0, 1) or \"mean\"")
  }
  Label[!IsMean] <- as.character(Probability)
  list(label = Label, isMean = IsMean, probability = Probability)
}

# The statistic of every level, in the order of the levels: all quantiles
# from one call, the mean where it was requested.
.levelValues <- function(levels, quantile, mean) {
  OUT <- numeric(length(levels$label))
  if (length(levels$probability)) OUT[!levels$isMean] <- quantile(levels$probability)
  if (any(levels$isMean)) OUT[levels$isMean] <- mean()
  OUT
}
