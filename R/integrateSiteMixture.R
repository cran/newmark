#' Site mixture of an analytic rock mixture by Gauss-Hermite quadrature
#'
#' Transforms a rock spectrum represented at each oscillator period as a finite
#' lognormal mixture into the corresponding site spectrum, integrating the
#' site factor of [evaluateSiteFactor()] over the rock PGA of each component
#' with Gauss-Hermite quadrature. Every quadrature node becomes one component
#' of the returned mixture, so the result is again an analytic mixture:
#' [summarizeMixture()] reduces it to quantiles and mean, and
#' [refineSiteMixture()] chooses the order.
#'
#' @section Mathematics:
#' At period \eqn{T_n} the rock component \eqn{k} gives
#' \eqn{\ln S_a \sim N(\mu_k, \sigma_k)} with weight \eqn{w_k}, and the
#' component of the same label at \eqn{T_n = 0} gives
#' \eqn{\ln PGA \sim N(\mu_{0k}, \sigma_{0k})}. With the Gauss-Hermite rule
#' of `order` nodes \eqn{x_m} and weights \eqn{\omega_m} normalized to one
#' (Golub-Welsch, weight \eqn{e^{-x^2}}), the PGA nodes are
#' \eqn{\ln PGA_{km} = \mu_{0k} + \sqrt{2} \sigma_{0k} x_m}. With
#' \eqn{\rho = \rho_{BJ}(\max(T_n, 0.01), 0.01)} of [correlateResiduals()],
#' \eqn{\rho_k = \rho} and \eqn{r_k = \sigma_k / \sigma_{0k}} when
#' \eqn{\sigma_{0k} > 0} and both zero otherwise, the ordinate conditional on
#' the PGA node is
#' \deqn{\ln S_a \mid PGA_{km} \sim N(\mu_k + \rho_k r_k (\ln PGA_{km} -
#'   \mu_{0k}), \sqrt{1 - \rho_k^2}\, \sigma_k),}
#' and with \eqn{(\mu_F, \sigma_F)} of the factor at \eqn{PGA_{km}} the site
#' component \eqn{(k, m)}, labelled `"k:m"`, is
#' \deqn{\mu_{km} = \mu_k + \rho_k r_k (\ln PGA_{km} - \mu_{0k}) + \mu_F,
#'   \quad \sigma_{km} = \sqrt{(1 - \rho_k^2) \sigma_k^2 + \sigma_F^2},
#'   \quad w_{km} \propto w_k \omega_m,}
#' the weights renormalized to one at every period. The site residual is
#' taken independent of the rock residual given PGA. When `VS30` equals
#' `Vref` the site mixture is the rock mixture. No random numbers are drawn;
#' the error of the order is not measured here (see [refineSiteMixture()]).
#'
#' @param mixture A table with one row per `(Tn, k)` and columns `Tn` (oscillator period, s),
#'   `k` (component label), `muLn`, `sigmaLn` (natural-log location and
#'   scale, `sigmaLn` non-negative) and `w` (positive weights summing to one
#'   at every period). It must hold the PGA rows (`Tn == 0`), and every
#'   component used at another period must have a PGA row of the same label.
#' @param VS30 One target velocity in m/s.
#' @param Vref One reference velocity in m/s; 760 or 3000 unless it equals
#'   `VS30`.
#' @param model Site factor model of [evaluateSiteFactor()], `"ST20"` or
#'   `"ST17"`.
#' @param order One integer, at least 2: Gauss-Hermite nodes per component.
#' @return A data.table with `Tn` (s), `k`, `muLn`, `sigmaLn` and `w`,
#'   ordered by `Tn` and `k`: the site spectrum as a lognormal mixture at
#'   every period, conditional on the target velocity and the model, its
#'   population being that of the rock mixture. `mixture` is not modified.
#' @export
#' @examples
#' Mixture <- data.table::CJ(Tn = c(0, 0.2, 1), k = c("a", "b"))
#' Mixture[, `:=`(muLn = log(0.4) + 0.1 * (k == "b") - 0.5 * Tn, sigmaLn = 0.5, w = 0.5)]
#' integrateSiteMixture(Mixture, VS30 = 360, Vref = 760, model = "ST20", order = 16)
integrateSiteMixture <- function(mixture, VS30, Vref, model, order) {
  .checkSiteCondition(VS30 = VS30, Vref = Vref, model = model)
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
      order != as.integer(order) || order < 2L) {
    stop("order must be one integer >= 2")
  }
  nq <- as.integer(order)
  Input <- .checkMixture(mixture)
  if (VS30 == Vref) return(Input)
  PGA <- Input[Tn == 0, .(k, muLn.PGA = muLn, sigmaLn.PGA = sigmaLn)]
  Rule <- .gaussHermiteRule(nq)
  LIST <- lapply(sort(unique(Input$Tn)), function(period) {
    Components <- copy(Input[Tn == period])
    Components[PGA, on = "k", `:=`(muLn.PGA = i.muLn.PGA, sigmaLn.PGA = i.sigmaLn.PGA)]
    if (anyNA(Components$muLn.PGA) || anyNA(Components$sigmaLn.PGA)) {
      stop("component at Tn=", period, " has no matching Tn == 0 PGA component")
    }
    LnPgaNodes <- outer(X = sqrt(2) * Components$sigmaLn.PGA, Y = Rule$x) + Components$muLn.PGA
    Site <- evaluateSiteFactor(PGA = exp(as.vector(LnPgaNodes)), Tn = period, VS30 = VS30,
                               Vref = Vref, model = model)
    Rho <- correlateResiduals(Tn = max(period, 0.01), T0 = 0.01)
    RhoComponent <- ifelse(Components$sigmaLn.PGA > 0, Rho, 0)
    Ratio <- ifelse(Components$sigmaLn.PGA > 0, Components$sigmaLn / Components$sigmaLn.PGA, 0)
    MuConditional <- as.vector(Components$muLn + RhoComponent * Ratio * (LnPgaNodes - Components$muLn.PGA))
    SigmaConditional <- as.vector(outer(X = sqrt((1 - RhoComponent^2) * Components$sigmaLn^2),
                                        Y = rep(x = 1, times = nq)))
    data.table(
      Tn = period,
      k = as.vector(outer(X = Components$k, Y = seq_len(nq),
                          FUN = function(component, node) paste0(component, ":", node))),
      muLn = MuConditional + Site$muLnF,
      sigmaLn = sqrt(SigmaConditional^2 + Site$sigmaLnF^2),
      w = as.vector(outer(X = Components$w, Y = Rule$w))
    )
  })
  OUT <- rbindlist(l = LIST, use.names = TRUE)
  OUT[, w := w / sum(w), by = Tn]
  .checkMixture(OUT)
}

# Gauss-Hermite rule of n nodes for the weight exp(-x^2): nodes are the
# eigenvalues of the Golub-Welsch Jacobi matrix (off-diagonal sqrt(i / 2));
# weights are 1 / sum_j h_j(x)^2 over the orthonormal Hermite polynomials,
# normalized to sum one, so sum(w * f(mu + sqrt(2) * sigma * x)) estimates
# E[f(Y)] for Y ~ N(mu, sigma).
.gaussHermiteRule <- function(n) {
  i <- seq_len(n - 1L)
  Jacobi <- matrix(data = 0, nrow = n, ncol = n)
  Jacobi[cbind(i, i + 1L)] <- sqrt(i / 2)
  Jacobi[cbind(i + 1L, i)] <- Jacobi[cbind(i, i + 1L)]
  Eigen <- eigen(x = Jacobi, symmetric = TRUE)
  Order <- order(Eigen$values)
  Nodes <- Eigen$values[Order]
  Polynomial.previous <- rep(x = 1, times = n)
  SumSquares <- Polynomial.previous^2
  Polynomial.current <- sqrt(2) * Nodes
  SumSquares <- SumSquares + Polynomial.current^2
  for (j in seq_len(max(n - 2L, 0L))) {
    Polynomial.next <- (Nodes * Polynomial.current - sqrt(j / 2) * Polynomial.previous) / sqrt((j + 1) / 2)
    SumSquares <- SumSquares + Polynomial.next^2
    Polynomial.previous <- Polynomial.current
    Polynomial.current <- Polynomial.next
  }
  Weights <- 1 / SumSquares
  if (any(!is.finite(Weights) | Weights <= 0) || sum(Weights) <= 0) {
    stop("Gauss-Hermite weights of order ", n, " are not representable")
  }
  list(x = Nodes, w = Weights / sum(Weights))
}

# The analytic mixture contract: the five columns, one row per (Tn, k),
# valid parameters, PGA rows and unit weight at every period; returns a copy
# with character labels ordered by Tn and k.
.checkMixture <- function(mixture) {
  Required <- c("Tn", "k", "muLn", "sigmaLn", "w")
  if (!is.data.frame(mixture) || !all(Required %in% names(mixture))) {
    stop("mixture is missing columns: ", paste(setdiff(Required, names(mixture)), collapse = ", "))
  }
  OUT <- copy(as.data.table(mixture)[, Required, with = FALSE])
  OUT[, k := as.character(k)]
  if (OUT[, any(!is.finite(Tn) | Tn < 0 | is.na(k) | !nzchar(k) | !is.finite(muLn) |
                  !is.finite(sigmaLn) | sigmaLn < 0 | !is.finite(w) | w <= 0)]) {
    stop("mixture contains invalid parameters")
  }
  if (!any(OUT$Tn == 0)) stop("mixture must include PGA rows with Tn == 0")
  if (anyDuplicated(OUT[, .(Tn, k)])) stop("mixture contains duplicate (Tn, k) rows")
  if (OUT[, .(weightSum = sum(w)), by = Tn][, any(abs(weightSum - 1) > 1e-12)]) {
    stop("mixture weights must sum to one at every Tn")
  }
  setorder(OUT, Tn, k)
  OUT[]
}
