#' Baker and Jayaram (2008) inter-period correlation model
#'
#' Computes the Pearson correlation coefficient between the residuals of
#' ln(Sa) at two oscillator periods, using the canonical predictive equation
#' of Baker and Jayaram (2008), Earthquake Spectra 24(1), 299-317
#' (DOI 10.1193/1.2857544), Eqs. 5 and 6 verbatim.
#'
#' The model has no dependence on rupture distance or magnitude.
#' \code{Rrup} is retained for signature compatibility and is silently ignored.
#' Periods of zero are treated as 0.01 s (standard PGA proxy in B&J 2008).
#'
#' Form (Eqs. 5-6 of B&J 2008):
#'   C1 = 1 - cos(pi/2 - 0.366 * ln(Tmax / max(Tmin, 0.109)))
#'   C2 = (1 - 0.105 * (1 - 1/(1 + exp(100*Tmax - 5))) * (Tmax-Tmin)/(Tmax-0.0099))
#'         if Tmax < 0.2; else 0
#'   C3 = C2 if Tmax < 0.109; else C1
#'   C4 = C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi * Tmin / 0.109))
#'
#'   rho = C2          if Tmax < 0.109
#'       = C1          else if Tmin > 0.109
#'       = min(C2, C4) else if Tmax < 0.2
#'       = C4          otherwise
#'
#' Validity: 0.01 <= T1, T2 <= 10 s.
#'
#' @param Tn   Numeric scalar. Target spectral period in seconds.
#' @param T0   Numeric scalar. Reference period in seconds. Default 0.01
#'   (used as a numerical proxy for PGA). Values <= 0 are silently
#'   converted to 0.01.
#' @param Rrup Numeric scalar. Not used; retained for backward compatibility.
#'
#' @return Numeric scalar correlation coefficient rho(T0, Tn).
#'
#' @examples
#' rhoBJ(Tn = 0.3, T0 = 0.01)   # rho(PGA, 0.3 s)
#' rhoBJ(Tn = 1.0, T0 = 0.5)
#'
#' @references
#' Baker, J. W. & Jayaram, N. (2008). Correlation of Spectral Acceleration
#' Values from NGA Ground Motion Models. Earthquake Spectra, 24(1), 299-317.
#' DOI: 10.1193/1.2857544.
#'
#' @export
rhoBJ <- function(Tn, T0 = 0.01, Rrup = 100) {
  if (T0 <= 0) T0 <- 0.01
  if (Tn <= 0) Tn <- 0.01
  stopifnot(length(Tn) == 1L, length(T0) == 1L)
  stopifnot(is.finite(Tn), is.finite(T0))

  if (any(c(T0, Tn) < 0.01) || any(c(T0, Tn) > 10)) {
    warning("rhoBJ: B&J 2008 is calibrated for T in [0.01, 10] s. ",
            "Extrapolation outside this range is not supported by the source paper.")
  }

  Tmin <- min(T0, Tn)
  Tmax <- max(T0, Tn)

  C1 <- 1 - cos(pi / 2 - 0.366 * log(Tmax / max(Tmin, 0.109)))

  if (Tmax < 0.2) {
    C2 <- 1 - 0.105 * (1 - 1 / (1 + exp(100 * Tmax - 5))) *
      (Tmax - Tmin) / (Tmax - 0.0099)
  } else {
    C2 <- 0
  }

  C3 <- if (Tmax < 0.109) C2 else C1
  C4 <- C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi * Tmin / 0.109))

  if (Tmax < 0.109) {
    rho <- C2
  } else if (Tmin > 0.109) {
    rho <- C1
  } else if (Tmax < 0.2) {
    rho <- min(C2, C4)
  } else {
    rho <- C4
  }

  rho
}
