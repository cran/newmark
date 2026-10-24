#' Correlation of spectral residuals between two oscillator periods
#'
#' Evaluates the Baker and Jayaram (2008) correlation between the logarithmic
#' residuals of spectral acceleration at two oscillator periods, the
#' coefficient that couples the ordinates of a correlated spectral sample
#' ([sampleSpectra()]). It depends on the pair of periods only, not on
#' magnitude or distance.
#'
#' @section Mathematics:
#' With \eqn{T_{min} = \min(T_n, T_0)} and \eqn{T_{max} = \max(T_n, T_0)},
#' equations 5 and 6 of Baker and Jayaram (2008) are
#' \deqn{C_1 = 1 - \cos\left(\frac{\pi}{2} - 0.366 \ln \frac{T_{max}}{\max(T_{min}, 0.109)}\right)}
#' \deqn{C_2 = 1 - 0.105 \left(1 - \frac{1}{1 + e^{100 T_{max} - 5}}\right) \frac{T_{max} - T_{min}}{T_{max} - 0.0099}}
#' when \eqn{T_{max} < 0.2} and \eqn{C_2 = 0} otherwise,
#' \eqn{C_3 = C_2} when \eqn{T_{max} < 0.109} and \eqn{C_3 = C_1} otherwise,
#' \deqn{C_4 = C_1 + \frac{1}{2} \left(\sqrt{C_3} - C_3\right) \left(1 + \cos \frac{\pi T_{min}}{0.109}\right)}
#' and the correlation is \eqn{C_2} when \eqn{T_{max} < 0.109}, \eqn{C_1}
#' when \eqn{T_{min} > 0.109}, \eqn{\min(C_2, C_4)} when
#' \eqn{T_{max} < 0.2} and \eqn{C_4} otherwise, tried in that order. The
#' coefficient is symmetric in the two periods. The fit is calibrated for
#' periods from 0.01 to 10 s; a pair outside that range is still evaluated
#' and a warning is signalled. As in `newmark` 1.7.1, a period of zero or less
#' stands for peak ground acceleration and is evaluated at 0.01 s, a
#' convention of the package and not of the paper.
#'
#' @param Tn,T0 Oscillator periods in seconds: finite numeric vectors of
#'   equal length, or one of them of length one.
#' @return Dimensionless correlation coefficients, one per period pair.
#' @references Baker, J. W. and Jayaram, N. (2008). Correlation of spectral
#'   acceleration values from NGA ground motion models. *Earthquake
#'   Spectra*, 24(1), 299-317. \doi{10.1193/1.2857544}
#' @export
#' @examples
#' correlateResiduals(Tn = c(0, 0.3, 1.5), T0 = 0.01)
correlateResiduals <- function(Tn, T0) {
  if (!is.numeric(Tn) || !is.numeric(T0) || !length(Tn) || !length(T0) ||
      any(!is.finite(Tn)) || any(!is.finite(T0)) ||
      (length(Tn) != length(T0) && min(length(Tn), length(T0)) != 1L)) {
    stop("Tn and T0 must be finite periods of equal length, or one of them a single period")
  }
  Tn[Tn <= 0] <- 0.01
  T0[T0 <= 0] <- 0.01
  if (any(c(Tn, T0) < 0.01 | c(Tn, T0) > 10)) {
    warning("Baker and Jayaram (2008) is calibrated for periods from 0.01 to 10 s; ",
            "the correlation outside that range is an extrapolation")
  }
  Tmin <- pmin(T0, Tn)
  Tmax <- pmax(T0, Tn)
  C1 <- 1 - cos(pi / 2 - 0.366 * log(Tmax / pmax(Tmin, 0.109)))
  C2 <- ifelse(Tmax < 0.2, 1 - 0.105 * (1 - 1 / (1 + exp(100 * Tmax - 5))) * (Tmax - Tmin) / (Tmax - 0.0099), 0)
  C3 <- ifelse(Tmax < 0.109, C2, C1)
  C4 <- C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi * Tmin / 0.109))
  fcase(Tmax < 0.109, C2,
        Tmin > 0.109, C1,
        Tmax < 0.2, pmin(C2, C4),
        Tmax >= 0.2, C4)
}
