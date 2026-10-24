#' Paired realizations of a site spectrum from a rock quantile spectrum
#'
#' Draws, at every oscillator period, `NS` paired realizations of the rock PGA, the
#' rock ordinate and the site ordinate from a rock spectrum given by its
#' mean and quantiles and from the site factor of [evaluateSiteFactor()].
#' The realizations are returned before any summary; reduce them per period
#' with [summarizeProfiles()] (see the example). [integrateSiteSpectrum()]
#' evaluates the same model without random numbers.
#'
#' @section Mathematics:
#' At PGA and at each period the rock quantile function
#' \eqn{Q(u) = \ln S_a} is that of [integrateSiteSpectrum()]: the Hyman
#' monotone spline through the tabulated knots \eqn{(p_k, \ln S_{a,k})},
#' continued below the first and above the last knot by the adjacent secant
#' slopes (floored at zero) and clamped to \eqn{u \in [0, 1]}. It is
#' non-decreasing, so every copula uniform maps to the rock quantile of its
#' rank. At period \eqn{T_n}, \eqn{(Z_{0j}, Z_{nj})}, \eqn{j = 1,
#' \ldots, n_s}, are drawn from the bivariate normal with unit variances and
#' correlation \eqn{\rho_{BJ}(T_n, 0)} of [correlateResiduals()]
#' (`mvtnorm::rmvnorm()`), a Gaussian copula between PGA and the ordinate,
#' and each marginal sample is re-centred on its tabulated mean:
#' \deqn{PGA_j = e^{Q_0(\Phi(Z_{0j}))} \frac{\overline{PGA}}{\langle
#'   e^{Q_0(\Phi(Z_0))} \rangle}, \quad S_{a,j} = e^{Q_n(\Phi(Z_{nj}))}
#'   \frac{\bar S_a}{\langle e^{Q_n(\Phi(Z_n))} \rangle},}
#' with \eqn{\langle \cdot \rangle} the sample mean. Then
#' \eqn{\ln F_j \sim N(\mu_F(PGA_j), \sigma_F(PGA_j))} independently and
#' \eqn{S_{aF,j} = S_{a,j} F_j}. When `VS30` equals `Vref`, for any
#' reference velocity, \eqn{F} is the unit factor (\eqn{\mu_F = 0},
#' \eqn{\sigma_F = 0}) and \eqn{S_{aF,j} = S_{a,j}}.
#'
#' The triple \eqn{(PGA_j, S_{a,j}, S_{aF,j})} is paired within one period
#' only: every period draws its own copula sample, so rows with the same
#' `realization` at different periods do not form one spectrum. Periods are drawn
#' in increasing order, each consuming from the caller's random stream the
#' `NS` copula draws and then, where the factor dispersion is positive, the
#' `NS` factor residuals (`stats::rnorm` draws nothing where the dispersion
#' is zero, as ST20 gives from 3000 m/s against a 760 m/s reference); when
#' `VS30` equals `Vref` no factor residual is drawn. Sample statistics carry Monte Carlo
#' error of order \eqn{n_s^{-1/2}}.
#'
#' @param spectrum The rock quantile spectrum, see [integrateSiteSpectrum()].
#' @param VS30 One target velocity in m/s.
#' @param Vref One reference velocity in m/s; 760 or 3000 unless it equals
#'   `VS30`.
#' @param model Site factor model of [evaluateSiteFactor()], `"ST20"` or
#'   `"ST17"`.
#' @param NS One positive integer: realizations per period.
#' @return A data.table with one row per period and realization: `Tn` (s),
#'   `realization` (1 to `NS`), `PGA`, `Sa` and `SaF` (g; rock PGA, rock ordinate
#'   and site ordinate of the realization), in the order drawn. Advances the
#'   random number stream.
#' @export
#' @examples
#' P <- c("0.05", "0.10", "0.16", "0.50", "0.84", "0.90", "0.95", "mean")
#' Factor <- c(0.45, 0.60, 0.72, 0.95, 1.30, 1.48, 1.78, 1.00)
#' Spectrum <- data.table::data.table(Tn = rep(c(0, 0.2, 1), each = 8), p = rep(P, 3),
#'                                    Sa = rep(c(0.30, 0.45, 0.24), each = 8) * rep(Factor, 3))
#' set.seed(919)
#' Draws <- simulateSiteSpectrum(Spectrum, VS30 = 270, Vref = 760, model = "ST20", NS = 300)
#' Draws[, summarizeProfiles(.SD, levels = c(0.16, 0.5, 0.84, "mean"), properties = "SaF"), by = Tn]
simulateSiteSpectrum <- function(spectrum, VS30, Vref, model, NS) {
  .checkSiteCondition(VS30 = VS30, Vref = Vref, model = model)
  Rock <- .checkQuantileSpectrum(spectrum)
  if (!is.numeric(NS) || length(NS) != 1L || !is.finite(NS) || NS < 1 || NS != as.integer(NS)) {
    stop("NS must be one positive integer")
  }
  PGA <- Rock[Tn == 0 & p == "mean", Sa]
  DP0 <- Rock[Tn == 0]
  P0 <- suppressWarnings(as.numeric(DP0$p))
  I0 <- which(!is.na(P0) & P0 > 0 & P0 < 1)
  LnQ0 <- .buildTailedQuantile(DP0[I0, .(p = as.numeric(p), Sa)])
  LIST <- lapply(sort(unique(Rock$Tn)), function(t) {
    DT <- Rock[Tn == t]
    Sa <- DT[p == "mean", Sa]
    if (length(Sa) != 1L) return(NULL)
    PNum <- suppressWarnings(as.numeric(DT$p))
    Idx <- which(!is.na(PNum) & PNum > 0 & PNum < 1)
    LnQn <- .buildTailedQuantile(DT[Idx, .(p = as.numeric(p), Sa)])
    Rho <- correlateResiduals(Tn = t, T0 = 0)
    Z2 <- mvtnorm::rmvnorm(NS, sigma = matrix(c(1, Rho, Rho, 1), nrow = 2L))
    LnPGA <- LnQ0(stats::pnorm(Z2[, 1L]))
    LnPGA <- log(exp(LnPGA) * (PGA / mean(exp(LnPGA))))
    LnSa <- LnQn(stats::pnorm(Z2[, 2L]))
    LnSa <- log(exp(LnSa) * (Sa / mean(exp(LnSa))))
    LnSaF <- LnSa
    # The reference condition is the unit factor for any Vref, whereas
    # evaluateSiteFactor() takes only the model references 760 and 3000.
    if (VS30 != Vref) {
      Factor <- evaluateSiteFactor(PGA = exp(LnPGA), Tn = t, VS30 = VS30, Vref = Vref, model = model)
      LnSaF <- LnSa + stats::rnorm(NS, mean = Factor$muLnF, sd = Factor$sigmaLnF)
    }
    data.table(Tn = t, realization = seq_len(NS), PGA = exp(LnPGA), Sa = exp(LnSa), SaF = exp(LnSaF))
  })
  rbindlist(LIST)
}
