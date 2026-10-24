#' Deformable sliding-mass displacement
#'
#' Evaluates the lognormal distribution of the nonzero seismic
#' displacement of a deformable sliding mass with yield coefficient
#' \eqn{k_y} and initial fundamental period \eqn{T_s} under an earthquake of moment
#' magnitude \eqn{M_w}, with one of the three flexible models of the
#' displacement ensemble, named by its author-year identifier `model`. The
#' spectral acceleration is read at the degraded period of the model:
#' \eqn{1.5 T_s} for `BT07` and `BM17`, \eqn{1.3 T_s} for `BM19`.
#'
#' Each published model treats the displacement as a mixed random variable:
#' a probability of negligible ("zero") displacement, with a threshold of
#' 1 cm in `BT07` and 0.5 cm in `BM17` and `BM19`, and a lognormal
#' distribution of the nonzero displacement. This function evaluates only
#' the nonzero branch, equation 14.5 of Bray (2007, p. 344), equation 4 of
#' Bray, Macedo and Travasarou (2018, p. 5) and equations 6, 8 and 10 of
#' Bray and Macedo (2023, pp. 6-7), and omits the probability of zero
#' displacement of equation 14.4 (p. 343), equations 2 and 3 (p. 5) and
#' equations 5, 7 and 9 (pp. 6-7) of the same sources. Its result is the
#' distribution of \eqn{D_n} conditional on the inputs and on a nonzero
#' displacement.
#'
#' `ky` and `Sa` are evaluated element by element: one of length one is
#' repeated for every element of the other, otherwise both must have the
#' same length, and `PGV` has the length of `Sa`; the function stops and
#' names the lengths otherwise.
#'
#' @section Mathematics:
#' With \eqn{D_n} in cm, \eqn{S_a} in g, \eqn{T_s} in s and
#' \eqn{x = \ln k_y}, \eqn{y = \ln S_a}:
#'
#' `BT07`, Bray and Travasarou (2007), equation 5, with
#' \eqn{a_0 = -1.10}, replaced by \eqn{-0.22} for nearly rigid masses
#' (\eqn{T_s < 0.05} s):
#' \deqn{\ln D_n = a_0 - 2.83 x - 0.333 x^2 + 0.566 x y + 3.04 y - 0.244 y^2 + 1.50 T_s + 0.278 (M_w - 7), \qquad \sigma = 0.66}
#'
#' `BM17`, Bray, Macedo and Travasarou (2018), equation 4, subduction
#' interface earthquakes, with
#' \eqn{(c_1, c_2, c_3) = (-5.864, -9.421, 0)} for \eqn{T_s < 0.1} s and
#' \eqn{(-6.896, 3.081, -0.803)} otherwise:
#' \deqn{\ln D_n = c_1 + 0.550 M_w + c_2 T_s + c_3 T_s^2 - 3.353 x - 0.390 x^2 + (3.060 + 0.538 x) y - 0.225 y^2, \qquad \sigma = 0.73}
#'
#' `BM19`, Bray and Macedo (2019) as corrected in 2023, shallow crustal
#' earthquakes. Ordinary motions (`pulse = FALSE`, equation 10):
#' \deqn{\ln D_n = a_0 + (2.649 + 0.344 x) y - 0.090 y^2, \qquad \sigma = 0.72}
#' with \eqn{a_0 = -4.684 + 0.603 M_w - 9.471 T_s - 2.482 x - 0.244 x^2}
#' for \eqn{T_s < 0.1} s and
#' \eqn{a_0 = -5.981 + 0.603 M_w + 3.223 T_s - 0.945 T_s^2 - 2.482 x - 0.244 x^2}
#' otherwise. Near-fault forward-directivity pulse motions
#' (`pulse = TRUE`), a property of the scenario declared by the analyst:
#' \deqn{\ln D_n = c_1 - b_1 x - b_2 x^2 + b_3 x y + b_4 y - b_5 y^2 + c_2 T_s + c_3 T_s^2 + b_m M_w + c_4 \ln PGV}
#' with PGV in cm/s and the coefficients of equation 6 for the maximum
#' component (`NFC = "D100"`: \eqn{b = (2.632, 0.278, 0.527, 1.978, 0.233)},
#' \eqn{b_m = 0.01}, \eqn{\sigma = 0.56}) or equation 8 for the median
#' component (`"D50"`: \eqn{b = (2.931, 0.319, 0.584, 2.261, 0.241)},
#' \eqn{b_m = 0.05}, \eqn{\sigma = 0.54}). \eqn{c_1} and \eqn{c_4} change at
#' PGV = 150 cm/s, \eqn{c_1, c_2, c_3} at \eqn{T_s = 0.1} s: D100
#' \eqn{c_1 = -6.724 / 1.991} and \eqn{(c_2, c_3) = (-2.744, 0)} below
#' 0.1 s, \eqn{c_1 = -6.951 / 1.764} and \eqn{(1.069, -0.498)} above,
#' \eqn{c_4 = 1.547 / -0.097}; D50 \eqn{c_1 = -7.497 / 2.480} and
#' \eqn{(-2.731, 0)}, \eqn{c_1 = -7.718 / -0.369} and \eqn{(1.031, -0.480)},
#' \eqn{c_4 = 1.458 / 0.025} (PGV up to / above 150 cm/s). For D50 above
#' 150 cm/s and \eqn{T_s < 0.1} s the paper's \eqn{c_1 = 2.480} differs from
#' the authors' spreadsheet (-0.148); the paper is followed.
#'
#' @param model One flexible model identifier: `"BT07"`, `"BM17"` or
#'   `"BM19"`.
#' @param ky Yield coefficients in g, finite and positive.
#' @param Ts The initial fundamental period of the sliding mass in seconds, one
#'   finite non-negative value.
#' @param Sa Spectral accelerations in g at the degraded period of the
#'   model, finite and positive.
#' @param Mw One finite positive moment magnitude.
#' @param PGV Peak ground velocities in cm/s, finite and positive, aligned
#'   with `Sa`, for example [estimatePGV()] of each PGA; only the pulse
#'   equations of `BM19` read them.
#' @param NFC Near-fault component of the `BM19` pulse equations, `"D100"`
#'   or `"D50"`.
#' @param pulse One logical: `TRUE` for a near-fault pulse motion of `BM19`,
#'   `FALSE` for ordinary motions.
#' @return A data.table with `IDn` (the model), `muLnD` and `sigmaLnD`, the
#'   mean and standard deviation of \eqn{\ln D_n} (\eqn{D_n} in cm) given the
#'   inputs and a nonzero displacement: \eqn{e^{\mu_{\ln D}}} is the median
#'   nonzero displacement and \eqn{\sigma_{\ln D}} its dispersion in natural
#'   logarithms. One row per element of the longer of `ky` and `Sa`.
#' @references Bray, J. D. and Travasarou, T. (2007). Simplified procedure
#'   for estimating earthquake-induced deviatoric slope displacements.
#'   *Journal of Geotechnical and Geoenvironmental Engineering*, 133(4),
#'   381-392. \doi{10.1061/(ASCE)1090-0241(2007)133:4(381)}
#'
#'   Bray, J. D. (2007). Simplified seismic slope displacement procedures. In
#'   K. D. Pitilakis (ed.), *Earthquake Geotechnical Engineering*, 327-353.
#'   Springer.
#'
#'   Bray, J. D., Macedo, J. and Travasarou, T. (2018). Simplified procedure
#'   for estimating seismic slope displacements for subduction zone
#'   earthquakes. *Journal of Geotechnical and Geoenvironmental
#'   Engineering*, 144(3), 04017124.
#'   \doi{10.1061/(ASCE)GT.1943-5606.0001833}
#'
#'   Bray, J. D. and Macedo, J. (2019). Procedure for estimating
#'   shear-induced seismic slope displacement for shallow crustal
#'   earthquakes. *Journal of Geotechnical and Geoenvironmental
#'   Engineering*, 145(12), 04019106.
#'
#'   Bray, J. D. and Macedo, J. (2023). Performance-based seismic assessment
#'   of slope systems. *Soil Dynamics and Earthquake Engineering*, 168,
#'   107835. \doi{10.1016/j.soildyn.2023.107835} The `BM19` coefficients
#'   implemented here are those printed in its equations 6, 8 and 10.
#' @export
#' @examples
#' evaluateFlexibleDisplacement(model = "BM19", ky = 0.1, Ts = 0.3, Sa = 0.8, Mw = 7,
#'                              PGV = estimatePGV(0.5), NFC = "D100", pulse = FALSE)
evaluateFlexibleDisplacement <- function(model, ky, Ts, Sa, Mw, PGV, NFC, pulse) {
  if (!is.character(model) || length(model) != 1L || !model %in% c("BT07", "BM17", "BM19")) {
    stop("model must be one of BT07, BM17 and BM19")
  }
  if (!is.numeric(ky) || !is.numeric(Sa) || !is.numeric(Ts) || !is.numeric(Mw) || !is.numeric(PGV) ||
      !length(ky) || !length(Sa) || length(Ts) != 1L || length(Mw) != 1L ||
      any(!is.finite(c(ky, Sa, Ts, Mw, PGV))) || any(c(ky, Sa, Mw, PGV) <= 0) || Ts < 0) {
    stop("ky, Sa, Mw and PGV must be finite and positive, Ts one finite non-negative period and Mw one magnitude")
  }
  if (length(PGV) != length(Sa)) stop("Sa and PGV must have the same length")
  .checkModelLengths(ky = ky, Sa = Sa)
  if (!is.character(NFC) || length(NFC) != 1L || !NFC %in% c("D100", "D50")) stop("NFC must be 'D100' or 'D50'")
  if (!is.logical(pulse) || length(pulse) != 1L || is.na(pulse)) stop("pulse must be one TRUE or FALSE")
  y <- log(Sa)
  x <- log(ky)
  if (model == "BT07") {
    a0 <- if (Ts < 0.05) -0.22 else -1.10
    return(data.table(IDn = model, muLnD = a0 - 2.83 * x - 0.333 * x^2 + 0.566 * x * y + 3.04 * y -
                        0.244 * y^2 + 1.50 * Ts + 0.278 * (Mw - 7), sigmaLnD = 0.66))
  }
  if (model == "BM17") {
    c1 <- if (Ts < 0.1) -5.864 else -6.896
    c2 <- if (Ts < 0.1) -9.421 else 3.081
    c3 <- if (Ts < 0.1) 0.0 else -0.803
    a0 <- c1 + 0.550 * Mw + c2 * Ts + c3 * Ts^2 - 3.353 * x - 0.390 * x^2
    return(data.table(IDn = model, muLnD = a0 + (3.060 + 0.538 * x) * y - 0.225 * y^2, sigmaLnD = 0.73))
  }
  if (!pulse) {
    a0 <- if (Ts < 0.1) {
      -4.684 + 0.603 * Mw - 9.471 * Ts - 2.482 * x - 0.244 * x^2
    } else {
      -5.981 + 0.603 * Mw + 3.223 * Ts - 0.945 * Ts^2 - 2.482 * x - 0.244 * x^2
    }
    return(data.table(IDn = model, muLnD = a0 + (2.649 + 0.344 * x) * y - 0.090 * y^2, sigmaLnD = 0.72))
  }
  # Pulse coefficients by component; c1 and c4 hold the values for PGV up to
  # 150 cm/s and above it.
  K <- list(D100 = list(c1 = if (Ts < 0.1) c(-6.724, 1.991) else c(-6.951, 1.764),
                        c2 = if (Ts < 0.1) -2.744 else 1.069, c3 = if (Ts < 0.1) 0 else -0.498,
                        c4 = c(1.547, -0.097), b = c(2.632, 0.278, 0.527, 1.978, 0.233), bm = 0.01, sd = 0.56),
            D50 = list(c1 = if (Ts < 0.1) c(-7.497, 2.480) else c(-7.718, -0.369),
                       c2 = if (Ts < 0.1) -2.731 else 1.031, c3 = if (Ts < 0.1) 0 else -0.480,
                       c4 = c(1.458, 0.025), b = c(2.931, 0.319, 0.584, 2.261, 0.241), bm = 0.05, sd = 0.54))[[NFC]]
  Regime <- 1L + (PGV > 150)
  data.table(IDn = model,
             muLnD = K$c1[Regime] - K$b[[1L]] * x - K$b[[2L]] * x^2 + K$b[[3L]] * x * y + K$b[[4L]] * y -
               K$b[[5L]] * y^2 + K$c2 * Ts + K$c3 * Ts^2 + K$bm * Mw + K$c4[Regime] * log(PGV),
             sigmaLnD = K$sd)
}
