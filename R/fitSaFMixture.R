#' Site-amplified spectra from an analytic ground-motion mixture
#'
#' Exact (sampling-free) composition of a scenario rock-hazard mixture
#' with the lognormal site-response dispersion. The rock demand per
#' period is a finite mixture of lognormal components indexed by
#' rupture-variant x GMPE-branch (the DSHA scenario engine output);
#' conditioning on the component and on its same-realization rock PGA,
#' \eqn{\ln SaF} is exactly normal:
#' \deqn{\mu_k(q) = m_k(T_n) + \rho\,\frac{s_k(T_n)}{s_k(0)}(q - m_k(0))
#'   + \mu_{\ln F}(T_n, V_{s30}, e^q)}
#' \deqn{\sigma_k(q)^2 = (1-\rho^2)\,s_k(T_n)^2
#'   + \sigma_{\ln F}(T_n, V_{s30}, e^q)^2}
#' so the pooled CDF is a mixture of one-dimensional Gaussian integrals
#' over the PGA channel, evaluated here by Gauss--Hermite quadrature
#' (\code{nq} nodes per component) and inverted monotonically. The same
#' inverter produces the rock quantiles \code{Sa}, so \code{Sa} and
#' \code{SaF} are quantiles of the same class of random object and
#' \code{AF = SaF/Sa} is the amplification of the p-quantile.
#' The within-component correlation \eqn{\rho(PGA, Sa(T_n))} is the
#' Baker--Jayaram (2008) within-event model (\code{\link{rhoBJ}}).
#' No random numbers are drawn: outputs are byte-stable.
#'
#' If \code{vs30 == vref} the site factor is deterministically 1 and
#' \code{SaF = Sa} exactly.
#'
#' @param gmm data.table with one row per (period, mixture component):
#'   columns \code{Tn} (s; must include \code{Tn == 0} rows for the PGA
#'   channel), \code{k} (component key, constant set across periods),
#'   \code{mu} and \code{sigma} (natural-log mean and sd of the
#'   component lognormal, sigma > 0), \code{w} (component weight,
#'   normalized internally per period).
#' @param vs30 numeric scalar - target Vs30 (m/s).
#' @param vref numeric scalar - reference rock Vs30 (m/s), default 760.
#' @param models character. Site-response model id: \code{"ST20"}
#'   (default) or \code{"ST17"} - see \code{\link{F_ST20}},
#'   \code{\link{F_ST17}}.
#' @param p_TARGET numeric probabilities reported alongside "mean".
#' @param nq integer, Gauss--Hermite nodes per component (default 32).
#'
#' @return data.table with columns \code{Tn, p, Sa, SaF, AF}; one row
#'   per (Tn, p) with p = \code{p_TARGET} plus \code{"mean"}.
#' @examples
#' # 6-component scenario mixture (2 rupture variants x 3 GMPE branches)
#' gmm <- data.table::CJ(Tn = c(0, 0.2, 1), k = 1:6)
#' gmm[, `:=`(mu = log(0.2) + 0.1 * k - 0.6 * Tn, sigma = 0.6, w = 1)]
#' fitSaFMixture(gmm, vs30 = 360, nq = 16)
#' @export
fitSaFMixture <- function(gmm, vs30, vref = 760, models = "ST20",
                          p_TARGET = c(0.05, 0.10, 0.16, 0.50, 0.84,
                                       0.90, 0.95),
                          nq = 32) {
  stopifnot(data.table::is.data.table(gmm))
  if (!all(c("Tn", "k", "mu", "sigma", "w") %in% names(gmm))) {
    stop("fitSaFMixture: 'gmm' must have columns Tn, k, mu, sigma, w")
  }
  if (gmm[, any(sigma <= 0)]) stop("fitSaFMixture: sigma must be > 0")
  FSite <- switch(models, ST20 = F_ST20, ST17 = F_ST17,
                  stop("fitSaFMixture: unknown site-response model: ",
                       models))

  P0 <- gmm[Tn == 0, .(k, mu0 = mu, s0 = sigma)]
  if (!nrow(P0)) stop("fitSaFMixture: 'gmm' must include Tn == 0 rows")
  GH <- .gaussHermite(nq)

  LIST <- lapply(sort(unique(gmm$Tn)), function(t) {
    DT <- gmm[Tn == t][P0, on = "k", `:=`(mu0 = i.mu0, s0 = i.s0)]
    if (anyNA(DT$mu0) || nrow(DT) != nrow(P0)) {
      stop("fitSaFMixture: component set at Tn = ", t,
           " does not match the Tn == 0 set")
    }
    W <- DT$w / sum(DT$w)

    ## rock mixture: quantiles and closed-form mean
    Sa <- c(exp(.mixQ(p_TARGET, DT$mu, DT$sigma, W)),
            sum(W * exp(DT$mu + DT$sigma^2 / 2)))

    if (vs30 == vref) {
      SaF <- Sa                       # deterministic unit factor (V1)
    } else {
      rho <- rhoBJ(Tn = max(t, 0.01), T0 = 0.01)
      Q   <- outer(sqrt(2) * DT$s0, GH$x) + DT$mu0     # K x nq, ln PGA
      AUX <- FSite(PGA = exp(as.vector(Q)), Tn = t, vs30 = vs30,
                   vref = vref)
      MU  <- as.vector(DT$mu + rho * (DT$sigma / DT$s0) * (Q - DT$mu0)) +
        AUX$muLnF
      SIG <- pmax(sqrt((1 - rho^2) * DT$sigma^2 + AUX$sdLnF^2), 1e-9)
      WQ  <- as.vector(outer(W, GH$w))
      SaF <- c(exp(.mixQ(p_TARGET, MU, SIG, WQ)),
               sum(WQ * exp(MU + SIG^2 / 2)))
    }

    data.table::data.table(
      Tn = t, p = c(sprintf("%.2f", p_TARGET), "mean"), Sa = Sa, SaF = SaF
    )
  })

  OUT <- data.table::rbindlist(LIST, use.names = TRUE)
  OUT[, AF := SaF / Sa]
  data.table::setorder(OUT, Tn, p)
  OUT[]
}


# Gauss-Hermite nodes/weights (physicists'), Golub-Welsch on the Jacobi
# matrix; weights returned already normalized to sum 1 (probabilist use:
# E[g(Z)] = sum(w * g(sqrt(2) * x)) for Z ~ N(0, 1)).
.gaussHermite <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- sqrt(i / 2)
  J[cbind(i + 1L, i)] <- sqrt(i / 2)
  E <- eigen(J, symmetric = TRUE)
  list(x = E$values, w = E$vectors[1L, ]^2)
}


# Quantiles of a finite normal mixture in ln-space: monotone CDF
# inversion by uniroot. mu/sig/w are flattened component vectors.
.mixQ <- function(p, mu, sig, w) {
  LO <- min(mu - 9 * sig)
  HI <- max(mu + 9 * sig)
  vapply(p, function(pi) {
    stats::uniroot(function(x) sum(w * stats::pnorm((x - mu) / sig)) - pi,
                   lower = LO, upper = HI, tol = 1e-10)$root
  }, numeric(1))
}
