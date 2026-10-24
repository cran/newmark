#' Correlated sample of spectral accelerations
#'
#' Draws `NS` joint realizations of the spectral accelerations at the
#' oscillator periods `Tn` through a Gaussian copula: the dependence comes from the
#' correlation matrix of the log residuals and each marginal from the
#' spectrum, given either as tabulated quantiles or as an analytic
#' lognormal mixture. The rows that share `realization` form one paired realization
#' of every ordinate.
#'
#' @section Mathematics:
#' Standard normal vectors \eqn{Z_s \sim N(0, C)}, \eqn{s = 1, \dots, N_S},
#' are drawn with `mvtnorm::rmvnorm` (eigen factorization), after adding
#' \eqn{|\lambda_{min}| + 10^{-8}} to the diagonal of \eqn{C} when its
#' smallest eigenvalue is not positive, and mapped to uniforms
#' \eqn{U = \Phi(Z)}. With a truncation \eqn{z_t}, the uniforms are
#' remapped to \eqn{\Phi(-z_t) + U\,(\Phi(z_t) - \Phi(-z_t))}, which keeps
#' their ranks and restricts them to \eqn{[\Phi(-z_t), \Phi(z_t)]}. Each
#' column is then inverted through the quantile function \eqn{Q} of
#' \eqn{\ln S_a} of its marginal and rescaled to its mean,
#' \deqn{S_{a,s} = e^{Q(U_s)}\, \frac{\bar S_a}{\frac{1}{N_S} \sum_{r=1}^{N_S} e^{Q(U_r)}};}
#' this finite-sample recentering multiplies every draw of one period by one
#' positive constant: it keeps their ranks, but no fixed bound on \eqn{S_a}.
#'
#' For a quantile spectrum, \eqn{Q} is [buildQuantileSpline()] of the rows
#' at that exact period and \eqn{\bar S_a} is its `"mean"` row; a period
#' without a mean row takes the sample mean, which leaves the draws at
#' \eqn{e^{Q(U_s)}}. A quantile spectrum is sampled untruncated, so \eqn{U}
#' spans (0, 1) and reaches the end pieces that the spline extrapolates
#' beyond its first and last tabulated probabilities. For a mixture
#' \eqn{F(x) = \sum_k w_k \Phi((x - \mu_k)/\sigma_k)} of \eqn{\ln S_a},
#' \eqn{Q} is the bisection inverse of the untruncated \eqn{F}, so with a
#' truncation the draws of a tabulated period lie, before the recentering,
#' between the mixture quantiles of probability \eqn{\Phi(-z_t)} and
#' \eqn{\Phi(z_t)}. The mean is \eqn{\sum_k w_k e^{\mu_k + \sigma_k^2/2}},
#' or with truncation that of the components truncated individually at
#' \eqn{\pm z_t},
#' \deqn{\sum_k w_k e^{\mu_k + \sigma_k^2/2} \frac{\Phi(z_t - \sigma_k) - \Phi(-z_t - \sigma_k)}{\Phi(z_t) - \Phi(-z_t)},}
#' while the draws invert the untruncated mixture over the restricted
#' probabilities. A mixture period that is not tabulated interpolates
#' \eqn{Q} linearly and the mean geometrically in the coordinate
#' \eqn{\ln(T + 10^{-8})} of the tabulated periods, the requested period
#' entering as \eqn{\ln T} (\eqn{T > 0}) or \eqn{-\infty} (\eqn{T = 0}),
#' between the two neighbouring tabulated periods, and takes the end period
#' outside them.
#'
#' Validation, including the quantile and mean rows of every requested
#' period, happens before any draw, so an invalid request leaves the random
#' stream untouched; a valid one consumes `NS` multivariate normal draws of
#' the caller's stream.
#'
#' @param spectrum A quantile spectrum (`Tn`, `p`, `Sa` in g, with at least
#'   three quantile rows at every requested period, each with a probability
#'   strictly between 0 and 1 and a finite positive `Sa`, and optionally a
#'   `"mean"` row with a finite positive `Sa`), or an analytic mixture of
#'   \eqn{\ln S_a} (`Tn`, `k`, `muLn`, `sigmaLn`, `w`, with rows at
#'   `Tn == 0` and weights summing to one at every period).
#' @param Tn Requested oscillator periods in seconds, finite and non-negative, in any
#'   order; the result lists them in this order.
#' @param correlation The symmetric correlation matrix of the log residuals
#'   at `Tn`, with unit diagonal, for example from [correlateResiduals()].
#'   Correlations given only between a reference period `Tn[1]` and the
#'   others, a vector `rho`, form the matrix
#'   `C <- diag(length(Tn)); C[1, -1] <- C[-1, 1] <- rho`.
#' @param NS Number of realizations, one positive integer.
#' @param zTrunc Symmetric truncation of the copula marginals in standard
#'   deviations, one finite positive number for an analytic mixture, or
#'   `NA` for untruncated sampling.
#' @return A data.table of paired realizations with `realization` (realization,
#'   1 to `NS`), `Tn` (s) and `Sa` (g), one row per realization and
#'   requested period, ordered by the position of the period in `Tn` and
#'   then by `realization`. The rows that share `realization` form one joint draw.
#' @export
#' @examples
#' Spectrum <- data.frame(Tn = rep(c(0, 1), each = 4),
#'                        p = rep(c("0.16", "0.50", "0.84", "mean"), 2),
#'                        Sa = c(0.30, 0.48, 0.72, 0.50, 0.10, 0.19, 0.31, 0.20))
#' C <- diag(2)
#' C[1, 2] <- C[2, 1] <- correlateResiduals(Tn = 1, T0 = 0.01)
#' set.seed(1)
#' Sample <- sampleSpectra(Spectrum, Tn = c(0, 1), correlation = C, NS = 20, zTrunc = NA)
#' Sample[, .(mean = mean(Sa)), by = Tn]
sampleSpectra <- function(spectrum, Tn, correlation, NS, zTrunc) {
  OUT <- .drawSpectra(spectrum, Tn = Tn, correlation = correlation, NS = NS, zTrunc = zTrunc)
  data.table(realization = rep(seq_len(nrow(OUT)), times = length(Tn)), Tn = rep(Tn, each = nrow(OUT)),
             Sa = as.vector(OUT))
}

# The draws of sampleSpectra() as a matrix: one row per realization, one
# column per requested period.
.drawSpectra <- function(spectrum, Tn, correlation, NS, zTrunc) {
  if (!is.data.frame(spectrum)) stop("spectrum must be a table")
  if (!is.numeric(Tn) || !length(Tn) || any(!is.finite(Tn) | Tn < 0) ||
      !is.numeric(NS) || length(NS) != 1L || !is.finite(NS) || NS != as.integer(NS) || NS < 1L) {
    stop("invalid correlated spectral-sample request: Tn must be finite non-negative periods and NS one positive integer")
  }
  if (!is.matrix(correlation) || !is.numeric(correlation) || any(dim(correlation) != length(Tn)) ||
      any(!is.finite(correlation)) || !isSymmetric(unname(correlation)) || any(diag(correlation) != 1)) {
    stop("correlation must be a finite symmetric matrix with unit diagonal, one row and column per period")
  }
  if (length(zTrunc) != 1L || !(is.na(zTrunc) || (is.numeric(zTrunc) && is.finite(zTrunc) && zTrunc > 0))) {
    stop("zTrunc must be one finite positive value, or NA for untruncated sampling")
  }
  NS <- as.integer(NS)
  DT <- as.data.table(spectrum)
  IsAnalytic <- .classifySpectrum(DT) == "analytic"
  if (!is.na(zTrunc) && !IsAnalytic) {
    stop("truncated sampling is defined for the analytic mixture input")
  }
  if (IsAnalytic) {
    Mixture <- .checkMixture(DT)
  } else {
    Marginals <- lapply(Tn, function(period) {
      Rows <- DT[Tn == period]
      # A row with a missing p is kept, so that the spline rejects it.
      Quantile <- buildQuantileSpline(Rows[is.na(p) | p != "mean", .(p, Sa)])
      Mean <- Rows[p == "mean", Sa]
      if (any(!is.finite(Mean) | Mean <= 0)) {
        stop("spectrum has a mean row without a finite positive Sa at Tn = ", period)
      }
      list(quantile = Quantile, mean = Mean[1L])
    })
  }
  C <- correlation
  Eigenvalues <- eigen(x = C, symmetric = TRUE, only.values = TRUE)$values
  if (min(Eigenvalues) <= 0) {
    C <- C + diag(abs(min(Eigenvalues)) + 1e-8, nrow(C))
  }
  U <- stats::pnorm(mvtnorm::rmvnorm(n = NS, sigma = C))
  if (!is.na(zTrunc)) {
    U <- stats::pnorm(-zTrunc) + U * (stats::pnorm(zTrunc) - stats::pnorm(-zTrunc))
  }
  OUT <- matrix(NA_real_, NS, length(Tn))
  if (!IsAnalytic) {
    for (j in seq_along(Tn)) {
      Quantiles <- exp(Marginals[[j]]$quantile(U[, j]))
      Mean <- Marginals[[j]]$mean
      if (!is.finite(Mean)) Mean <- .portableMean(Quantiles)
      OUT[, j] <- Quantiles * (Mean / .portableMean(Quantiles))
    }
    return(OUT)
  }
  Periods <- sort(unique(Mixture$Tn))
  evaluate <- function(index, u) {
    Components <- Mixture[Tn == Periods[[index]]]
    Mean <- if (is.na(zTrunc)) {
      Components[, Reduce(`+`, w * exp(muLn + sigmaLn^2 / 2))]
    } else {
      Components[, Reduce(`+`, w * exp(muLn + sigmaLn^2 / 2) *
                         (stats::pnorm(zTrunc - sigmaLn) - stats::pnorm(-zTrunc - sigmaLn)) /
                         (stats::pnorm(zTrunc) - stats::pnorm(-zTrunc)))]
    }
    list(quantile = .invertMixture(p = u, muLn = Components$muLn, sigmaLn = Components$sigmaLn,
                                   w = Components$w),
         mean = Mean)
  }
  for (j in seq_along(Tn)) {
    Bracket <- .bracketPeriod(Tn[[j]], Periods)
    Lower <- evaluate(Bracket$index[[1L]], U[, j])
    QuantileLn <- Lower$quantile
    Mean <- Lower$mean
    if (Bracket$index[[1L]] != Bracket$index[[2L]]) {
      Upper <- evaluate(Bracket$index[[2L]], U[, j])
      QuantileLn <- Lower$quantile + Bracket$fraction * (Upper$quantile - Lower$quantile)
      Mean <- exp(log(Lower$mean) + Bracket$fraction * (log(Upper$mean) - log(Lower$mean)))
    }
    Quantiles <- exp(QuantileLn)
    OUT[, j] <- Quantiles * (Mean / .portableMean(Quantiles))
  }
  OUT
}

# The two passes of mean() accumulated in double, as mean() does only where the
# platform has no long double; mean() takes over when the result is not
# finite, as when the first sum overflows.
.portableMean <- function(x) {
  s <- Reduce(`+`, x) / length(x)
  s <- s + Reduce(`+`, x - s) / length(x)
  if (!is.finite(s)) return(mean(x))
  s
}

# Which of the two spectrum schemas a table follows: the analytic mixture of
# ln Sa (Tn, k, muLn, sigmaLn, w) or the quantile table (Tn, p, Sa).
.classifySpectrum <- function(spectrum) {
  IsAnalytic <- all(c("Tn", "k", "muLn", "sigmaLn", "w") %in% names(spectrum))
  IsQuantile <- all(c("Tn", "p", "Sa") %in% names(spectrum))
  if (IsAnalytic && IsQuantile) stop("spectrum input is ambiguous between analytic and quantile schemas")
  if (IsAnalytic) return("analytic")
  if (IsQuantile) return("quantile")
  if (any(c("k", "muLn", "sigmaLn", "w") %in% names(spectrum))) stop("analytic spectrum has an incomplete schema")
  stop("spectrum input has neither the analytic nor quantile schema")
}

# The tabulated periods that bracket a requested one in ln(T + 1e-8): the
# exact period, the end period outside the table, or the two neighbours and
# the fraction of the way between them.
.bracketPeriod <- function(period, periods) {
  Exact <- match(period, periods)
  if (!is.na(Exact)) return(list(index = c(Exact, Exact), fraction = 0))
  X <- log(periods + 1e-8)
  Target <- if (period == 0) -Inf else log(period)
  if (Target <= X[[1L]]) return(list(index = c(1L, 1L), fraction = 0))
  if (Target >= X[[length(X)]]) return(list(index = rep(length(X), 2L), fraction = 0))
  Lower <- findInterval(Target, X)
  list(index = c(Lower, Lower + 1L), fraction = (Target - X[[Lower]]) / (X[[Lower + 1L]] - X[[Lower]]))
}

# Quantiles of ln Sa under a lognormal mixture: the smallest x with
# sum_k w_k Phi((x - muLn_k) / sigmaLn_k) >= p, by 100 bisection steps from a
# bracket of 12 scale units around the components, widened while it does not
# contain every requested probability. A component with sigmaLn = 0 is an
# atom at muLn. The weight sums avoid long double.
.invertMixture <- function(p, muLn, sigmaLn, w) {
  if (!is.numeric(p) || !length(p) || any(!is.finite(p) | p <= 0 | p >= 1) ||
      !length(muLn) || length(muLn) != length(sigmaLn) || length(muLn) != length(w) ||
      any(!is.finite(muLn)) || any(!is.finite(sigmaLn) | sigmaLn < 0) ||
      any(!is.finite(w) | w < 0) || Reduce(`+`, w) <= 0) {
    stop("invalid analytic-mixture quantile request")
  }
  Weight <- w / Reduce(`+`, w)
  distribution <- function(x) {
    Terms <- vapply(seq_along(muLn), function(j) {
      if (sigmaLn[[j]] > 0) return(stats::pnorm(q = (x - muLn[[j]]) / sigmaLn[[j]]))
      as.numeric(x >= muLn[[j]])
    }, FUN.VALUE = numeric(length(x)))
    if (is.null(dim(Terms))) Terms <- matrix(Terms, nrow = length(x), ncol = length(muLn))
    as.numeric(Terms %*% Weight)
  }
  Scale <- pmax(sigmaLn, 1)
  Lower <- min(muLn - 12 * Scale)
  Upper <- max(muLn + 12 * Scale)
  Span <- max(Upper - Lower, 1)
  for (i in seq_len(32L)) {
    LowerProbability <- distribution(Lower)
    UpperProbability <- distribution(Upper)
    if (LowerProbability < min(p) && UpperProbability >= max(p)) break
    if (LowerProbability >= min(p)) Lower <- Lower - Span
    if (UpperProbability < max(p)) Upper <- Upper + Span
    Span <- 2 * Span
  }
  if (distribution(Lower) >= min(p) || distribution(Upper) < max(p)) {
    stop("analytic-mixture quantile could not be bracketed")
  }
  Lower <- rep(Lower, length(p))
  Upper <- rep(Upper, length(p))
  for (i in seq_len(100L)) {
    Middle <- (Lower + Upper) / 2
    MASK <- distribution(Middle) >= p
    Upper[MASK] <- Middle[MASK]
    Lower[!MASK] <- Middle[!MASK]
  }
  Upper
}

# Names resolved by data.table inside the displacement family.
utils::globalVariables(c(
  "AF", "AI", "Dn", "ID", "IDg", "IDm", "IDn", "PGA", "PGA.s", "PGV", "Sa", "Sa13", "Sa15", "SaF",
  "TR", "Tn", "Ts", "V1", "VS30", "Vref", "Vs30", "bad", "components", "fitErr", "hasPGA", "k",
  "ky", "level", "maxTn", "muLn", "nPGA", "p", "realization", "s", "sigmaLn", "siteID", "w", "weightSum",
  "weightsValid"
))
