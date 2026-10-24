#' Site spectrum of a rock quantile spectrum by quadrature
#'
#' Quantiles and mean of the spectral acceleration at a site of velocity
#' `VS30`, from a rock spectrum given at each oscillator period by its mean and its
#' quantiles, integrating the joint rock demand of PGA and the ordinate and
#' the site factor of [evaluateSiteFactor()] without random numbers. The
#' site quantiles are reported at the probability levels of the rock
#' spectrum. [refineSiteSpectrum()] chooses the order;
#' [simulateSiteSpectrum()] samples the same model; a rock spectrum given
#' only by its means is the conditional case of [evaluateSiteSa()].
#'
#' @section Mathematics:
#' At PGA and at each period the rock quantile function
#' \eqn{Q(u) = \ln S_a} is the Hyman monotone spline through the tabulated
#' knots \eqn{(p_k, \ln S_{a,k})}, continued below the first and above the
#' last knot by the adjacent secant slopes (floored at zero) and clamped to
#' \eqn{u \in [0, 1]}. The Gauss-Legendre rule of `order` nodes on
#' \eqn{(0, 1)} gives \eqn{u_i} and \eqn{w_i} (Golub-Welsch); with
#' \eqn{z_i = \Phi^{-1}(u_i)} and \eqn{\rho = \rho_{BJ}(T_n, 0)} of
#' [correlateResiduals()], the joint node \eqn{(i, j)} carries
#' \deqn{\ln PGA_i = Q_0(u_i) + c_0, \quad
#'   \ln S_{a,ij} = Q_n(\Phi(\rho z_i + \sqrt{1 - \rho^2} z_j)) + c_n,
#'   \quad W_{ij} = w_i w_j,}
#' a Gaussian copula between PGA and the ordinate, where \eqn{c_0} and
#' \eqn{c_n} re-centre each rock marginal on its tabulated mean
#' (\eqn{\sum_i w_i PGA_i = \overline{PGA}} and
#' \eqn{\sum_{ij} W_{ij} S_{a,ij} = \bar S_a}). Each joint node is one
#' lognormal component of the site ordinate,
#' \eqn{\mu_{ij} = \ln S_{a,ij} + \mu_F(PGA_i)} and
#' \eqn{\sigma_{ij} = \max(\sigma_F(PGA_i), 10^{-9})}, the site residual
#' being integrated analytically and taken independent of the rock residual
#' given PGA. The site quantile at a tabulated probability \eqn{p} solves
#' \eqn{\sum W_{ij} \Phi((x - \mu_{ij})/\sigma_{ij}) = p} (root to
#' \eqn{10^{-10}} on \eqn{[\min(\mu - 9\sigma), \max(\mu + 9\sigma)]}) and
#' the site mean is \eqn{\sum W_{ij} e^{\mu_{ij} + \sigma_{ij}^2/2}}. When
#' `VS30` equals `Vref` the site spectrum is the rock spectrum. Each period
#' is integrated with its own copula with PGA; no dependence between
#' periods is implied.
#'
#' @param spectrum A table with columns `Tn` (s), `p` (`"mean"` or a
#'   probability label) and `Sa` (g): exactly one row at `Tn == 0` and
#'   `p == "mean"`, and at PGA and every period at least three unique,
#'   positive, non-decreasing quantiles. Periods without a single mean row
#'   are skipped.
#' @param VS30 One target velocity in m/s.
#' @param Vref One reference velocity in m/s; 760 or 3000 unless it equals
#'   `VS30`.
#' @param model Site factor model of [evaluateSiteFactor()], `"ST20"` or
#'   `"ST17"`.
#' @param order One integer, at least 4: Gauss-Legendre nodes per latent
#'   dimension, \eqn{n^2} joint nodes per period.
#' @return A data.table ordered by `Tn` and `p` with `Tn`, `p` (the labels
#'   of the rock spectrum), `Sa` (the rock ordinate at that level, g) and
#'   `SaF` (the site ordinate at the same level, g), unrounded. Quantile and
#'   mean rows are marginal statistics at each period over the population
#'   of the rock spectrum and the site-factor residual; the integration
#'   error of `order` is not reported (see [refineSiteSpectrum()]). The
#'   random number stream is not used.
#' @export
#' @examples
#' P <- c("0.05", "0.10", "0.16", "0.50", "0.84", "0.90", "0.95", "mean")
#' Factor <- c(0.45, 0.60, 0.72, 0.95, 1.30, 1.48, 1.78, 1.00)
#' Spectrum <- data.table::data.table(Tn = rep(c(0, 0.2, 1), each = 8), p = rep(P, 3),
#'                                    Sa = rep(c(0.30, 0.45, 0.24), each = 8) * rep(Factor, 3))
#' integrateSiteSpectrum(Spectrum, VS30 = 270, Vref = 760, model = "ST20", order = 32)
integrateSiteSpectrum <- function(spectrum, VS30, Vref, model, order) {
  .checkSiteCondition(VS30 = VS30, Vref = Vref, model = model)
  Rock <- .checkQuantileSpectrum(spectrum)
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
      order != as.integer(order) || order < 4L) {
    stop("order must be one integer >= 4")
  }
  if (VS30 == Vref) {
    OUT <- Rock[, .(Tn, p, Sa, SaF = Sa)]
    setorder(OUT, Tn, p)
    return(OUT[])
  }
  nq <- as.integer(order)
  PGA <- Rock[Tn == 0 & p == "mean", Sa]
  Rule <- .gaussLegendreRule(nq)
  Z <- stats::qnorm(Rule$u)
  DP0 <- Rock[Tn == 0]
  P0 <- suppressWarnings(as.numeric(DP0$p))
  I0 <- which(!is.na(P0) & P0 > 0 & P0 < 1)
  LnQ0 <- .buildTailedQuantile(DP0[I0, .(p = as.numeric(p), Sa)])
  LnPGA <- LnQ0(Rule$u)
  LnPGA <- LnPGA + log(PGA / sum(Rule$w * exp(LnPGA)))
  W <- as.vector(outer(Rule$w, Rule$w))
  LIST <- lapply(sort(unique(Rock$Tn)), function(t) {
    DT <- Rock[Tn == t]
    Sa <- DT[p == "mean", Sa]
    if (length(Sa) != 1L) return(NULL)
    PNum <- suppressWarnings(as.numeric(DT$p))
    Idx <- which(!is.na(PNum) & PNum > 0 & PNum < 1)
    LnQn <- .buildTailedQuantile(DT[Idx, .(p = as.numeric(p), Sa)])
    Rho <- correlateResiduals(Tn = t, T0 = 0)
    Z1 <- outer(Rho * Z, rep(1, nq)) + outer(rep(1, nq), sqrt(pmax(1 - Rho^2, 0)) * Z)
    LnSa <- as.vector(LnQn(stats::pnorm(Z1)))
    LnSa <- LnSa + log(Sa / sum(W * exp(LnSa)))
    Factor <- evaluateSiteFactor(PGA = exp(rep(LnPGA, nq)), Tn = t, VS30 = VS30, Vref = Vref, model = model)
    Mu <- LnSa + Factor$muLnF
    Sigma <- pmax(Factor$sigmaLnF, 1e-9)
    data.table(
      Tn = t,
      p = c(DT$p[Idx], "mean"),
      Sa = c(DT$Sa[Idx], Sa),
      SaF = c(exp(.solveLogQuantile(p = PNum[Idx], muLn = Mu, sigmaLn = Sigma, w = W)),
              sum(W * exp(Mu + Sigma^2 / 2)))
    )
  })
  OUT <- rbindlist(LIST, use.names = TRUE)
  setorder(OUT, Tn, p)
  OUT[]
}

# The rock quantile spectrum shared by integrateSiteSpectrum() and
# simulateSiteSpectrum(): the three columns, one PGA mean row and some
# probability level; the "mean" label coerces to NA by design.
.checkQuantileSpectrum <- function(spectrum) {
  if (!is.data.frame(spectrum) || !all(c("Tn", "p", "Sa") %in% names(spectrum))) {
    stop("spectrum must be a table with columns Tn, p and Sa")
  }
  DT <- as.data.table(spectrum)[, .(Tn, p, Sa)]
  if (length(DT[Tn == 0 & p == "mean", Sa]) != 1L) {
    stop("spectrum needs exactly one Tn == 0, p == \"mean\" row for PGA")
  }
  PNum <- suppressWarnings(as.numeric(DT$p))
  if (!any(!is.na(PNum) & PNum > 0 & PNum < 1)) {
    stop("spectrum holds only means; its site ordinates are the conditional lognormal of evaluateSiteSa()")
  }
  DT
}

# Gauss-Legendre nodes and weights on (0, 1). For the n-point rule the
# Golub-Welsch Jacobi matrix has off-diagonal i / sqrt(4 i^2 - 1).
# Transforming its [-1, 1] weights to [0, 1] leaves the squared first
# eigenvector components normalized to one.
.gaussLegendreRule <- function(n) {
  Index <- seq_len(n - 1L)
  Jacobi <- matrix(0, n, n)
  OffDiagonal <- Index / sqrt(4 * Index^2 - 1)
  Jacobi[cbind(Index, Index + 1L)] <- OffDiagonal
  Jacobi[cbind(Index + 1L, Index)] <- OffDiagonal
  Eigen <- eigen(Jacobi, symmetric = TRUE)
  Order <- order(Eigen$values)
  list(u = (Eigen$values[Order] + 1) / 2, w = Eigen$vectors[1L, Order]^2)
}

# Monotone log-quantile function of a rock quantile spectrum, shared by
# integrateSiteSpectrum() and simulateSiteSpectrum(). Hyman interpolation
# preserves increasing and flat segments between knots; its cubic
# extrapolation can fold outside them, so the unobserved tails continue the
# adjacent secant slopes in probability space.
.buildTailedQuantile <- function(quantiles) {
  DT <- copy(as.data.table(quantiles))
  DT[, p := suppressWarnings(as.numeric(p))]
  DT <- DT[is.finite(p) & p > 0 & p < 1]
  setorder(DT, p)
  if (nrow(DT) < 3L || anyDuplicated(DT$p) || DT[, any(!is.finite(Sa) | Sa <= 0)] || any(diff(DT$Sa) < 0)) {
    stop("the rock quantile spectrum needs at least three unique, positive, non-decreasing quantiles at PGA and every period")
  }
  P <- DT$p
  Y <- log(DT$Sa)
  # splinefun warns when the Hyman filter adjusts slopes; the filtered spline
  # is the intended interpolant.
  Q <- suppressWarnings(stats::splinefun(x = P, y = Y, method = "hyman"))
  LeftSlope <- max((Y[[2L]] - Y[[1L]]) / (P[[2L]] - P[[1L]]), 0)
  N <- length(P)
  RightSlope <- max((Y[[N]] - Y[[N - 1L]]) / (P[[N]] - P[[N - 1L]]), 0)
  function(u) {
    u <- pmin(pmax(u, 0), 1)
    Value <- Q(pmin(pmax(u, P[[1L]]), P[[N]]))
    Left <- u < P[[1L]]
    Right <- u > P[[N]]
    Value[Left] <- Y[[1L]] + LeftSlope * (u[Left] - P[[1L]])
    Value[Right] <- Y[[N]] + RightSlope * (u[Right] - P[[N]])
    Value
  }
}

# Log-quantiles of a lognormal mixture with positive dispersions: the root of
# sum(w * pnorm((x - muLn) / sigmaLn)) = p for each probability.
.solveLogQuantile <- function(p, muLn, sigmaLn, w) {
  Lower <- min(muLn - 9 * sigmaLn)
  Upper <- max(muLn + 9 * sigmaLn)
  vapply(p, function(probability) {
    stats::uniroot(function(x) sum(w * stats::pnorm(q = (x - muLn) / sigmaLn)) - probability,
                   lower = Lower, upper = Upper, tol = 1e-10)$root
  }, numeric(1L))
}
