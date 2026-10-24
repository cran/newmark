#' Conditional site amplification factor of the NGA-East models
#'
#' Natural-log mean and standard deviation of the site amplification factor
#' \eqn{F} that converts a spectral ordinate at the reference condition
#' `Vref` into the ordinate at a site of time-averaged velocity `VS30`, given
#' the peak ground acceleration at the reference condition and the oscillator
#' period. Model `"ST20"` is the final NGA-East ergodic model (linear part of
#' Stewart et al. 2020, nonlinear part of Hashash et al. 2020, coefficients of
#' the USGS nshmp-lib table carried in `SiteAmp2020`); model `"ST17"` is the
#' 2017 generation (PEER reports 2017/04 and 2017/05). The pair defines the
#' lognormal \eqn{\ln F \sim N(\mu, \sigma)} conditional on PGA, period,
#' `VS30` and `Vref`; it is the factor, not an amplified spectrum (see
#' [evaluateSiteSa()], [integrateSiteMixture()], [integrateSiteSpectrum()]
#' and [simulateSiteSpectrum()], which integrate it over the rock demand).
#'
#' @section Mathematics:
#' *Model `"ST20"`.* The coefficients \eqn{c, V_1, V_2, V_f, V_c, f_3, f_4,
#' f_5, \sigma_{vc}, \sigma_l, \sigma_u, \sigma_c, f_{760i}, f_{760g},
#' f_{760is}, f_{760gs}} of the period are read from `SiteAmp2020`: the PGA
#' row for \eqn{T_n \le 0.01} s, the tabulated row at a tabulated period, and
#' otherwise the linear interpolation in \eqn{\ln T_n} between the bracketing
#' rows; a period above 10 s is clamped to 10 s with a warning. With
#' \eqn{w(v) = \min(\max(0.667 \ln(v/400)/\ln 1.5 + 0.1, 0.1), 0.767)} and
#' \eqn{F_{760}(v) = w f_{760i} + (1 - w) f_{760g}}, the linear term is
#' \eqn{f_{lin}(v) = f_v(v) + F_{760}(v)}, where \eqn{f_v(v)} is
#' \eqn{c \ln(V_1/760)} for \eqn{v \le V_1}, \eqn{c \ln(v/760)} up to
#' \eqn{V_2}, \eqn{c \ln(V_2/760)} up to 2000 m/s,
#' \eqn{c \ln(V_2/760) - (c \ln(V_2/760) + F_{760}(v)) \ln(v/2000)/\ln 1.5}
#' below 3000 m/s and \eqn{-F_{760}(v)} from 3000 m/s (Stewart et al. 2020,
#' Eqs. 3 and 5-7). The nonlinear slope is
#' \eqn{f_2(v) = f_4 [e^{f_5(\min(v, V_{NL}) - 360)} - e^{f_5(V_{NL} - 360)}]}
#' for \eqn{v < V_c} and 0 otherwise, with \eqn{V_{NL} = 760} m/s for
#' \eqn{T_n < 0.4} s and 3000 m/s above (Hashash et al. 2020, Eq. 3), and it
#' multiplies \eqn{F_{Rk} = \ln((PGA_{3000} + f_3)/f_3)}, where
#' \eqn{PGA_{3000} = PGA/2.275} when `Vref` is 760 m/s (PEER 2017/05,
#' Eq. 2.4) and \eqn{PGA} when it is 3000 m/s. The linear dispersion is
#' \eqn{\sigma_l - 2 D x + D x^2} with \eqn{D = \sigma_l - \sigma_{vc}} and
#' \eqn{x = (v - 200)/(V_f - 200)} below \eqn{V_f}, \eqn{\sigma_{vc}} up to
#' \eqn{V_2}, \eqn{\sigma_{vc} + (\sigma_u - \sigma_{vc})((v - V_2)/(2000 -
#' V_2))^2} below 2000 m/s, \eqn{\sigma_u (1 - \ln(v/2000)/\ln 1.5)} up to
#' 3000 m/s and 0 above (Stewart et al. 2020, Eq. 4); the nonlinear
#' dispersion is \eqn{\sigma_{f2} F_{Rk}} with \eqn{\sigma_{f2} = \sigma_c}
#' below 300 m/s, \eqn{\sigma_c (1 - \ln(v/300)/\ln(1000/300))} below
#' 1000 m/s and 0 above. For `Vref` = 760 m/s
#' \deqn{\mu = f_{lin}(V_{S30}) - f_{lin}(760) + (f_2(V_{S30}) - f_2(760))
#'   F_{Rk}, \quad \sigma = \sqrt{\sigma_{lin}^2(V_{S30}) +
#'   (\sigma_{f2}(V_{S30}) F_{Rk})^2},}
#' the reference side adding no dispersion; for `Vref` = 3000 m/s
#' \eqn{\mu = f_{lin}(V_{S30}) + f_2(V_{S30}) F_{Rk}} and the weighted
#' \eqn{F_{760}} dispersion \eqn{w f_{760is} + (1 - w) f_{760gs}} enters the
#' root sum of squares.
#'
#' *Model `"ST17"`.* The coefficients are interpolated linearly in
#' \eqn{T_n} from the tables of PEER 2017/04 (linear part, \eqn{F_{760}}) and
#' 2017/05 (nonlinear part, Table 2.1). The linear mean is
#' \eqn{c \ln(V_1/760)}, \eqn{c \ln(V_{S30}/760)} or
#' \eqn{c \ln(V_2/760) + (c/2) \ln(V_{S30}/V_2)} below \eqn{V_1}, up to
#' \eqn{V_2} and above; the nonlinear mean is
#' \eqn{f_2 \ln(PGA/b + 1)} for \eqn{V_{S30} < V_c} with
#' \eqn{f_2 = f_4 [e^{f_5(\min(V_{S30}, 3000) - 360)} - e^{f_5(3000 - 360)}]}
#' and \eqn{b = 2.275 f_3} (`Vref` 760) or \eqn{f_3} (`Vref` 3000), and
#' \eqn{F_{760}(T_n)} is added for `Vref` = 3000. The dispersion combines the
#' linear term (quadratic below \eqn{V_f} between 200 m/s and \eqn{V_f},
#' \eqn{\sigma_{vc}} up to \eqn{V_2}, quadratic towards 2000 m/s above), the
#' tabulated \eqn{F_{760}} dispersion for `Vref` = 3000 and the nonlinear
#' term \eqn{\sigma_{c} \ln(PGA/b + 1)} tapered from 300 to 1000 m/s.
#'
#' When `VS30` equals `Vref` the factor is the deterministic unit factor,
#' \eqn{\mu = 0} and \eqn{\sigma = 0}. Every result is a conditional evaluation at the supplied PGA; the
#' dispersion combines the epistemic uncertainties of the published median
#' model, not the uncertainty of the rock demand.
#'
#' @param PGA Peak ground accelerations in g at the reference condition,
#'   finite and non-negative.
#' @param Tn One oscillator period in seconds, finite and non-negative; 0 is
#'   PGA.
#' @param VS30 One target time-averaged shear-wave velocity in m/s.
#' @param Vref One reference velocity in m/s, 760 or 3000.
#' @param model `"ST20"` or `"ST17"`.
#' @return A data.table with `muLnF` and `sigmaLnF` (natural-log mean and
#'   standard deviation of \eqn{F}, one row per `PGA`) and `ID` (the model),
#'   unrounded.
#' @references Stewart, J.P. et al. (2020). Ergodic site amplification model
#'   for central and eastern North America. Earthquake Spectra 36(1).
#'   Hashash, Y.M.A. et al. (2020). Nonlinear site amplification model for
#'   ergodic seismic hazard analysis in central and eastern North America.
#'   Earthquake Spectra 36(1). Stewart, J.P. et al. (2017). PEER Report
#'   2017/04. Hashash, Y.M.A. et al. (2017). PEER Report 2017/05.
#' @export
#' @examples
#' evaluateSiteFactor(PGA = c(0.05, 0.2, 0.5), Tn = 0.2, VS30 = 300, Vref = 760,
#'                    model = "ST20")
#' evaluateSiteFactor(PGA = 0.1, Tn = 1.0, VS30 = 760, Vref = 3000, model = "ST20")
evaluateSiteFactor <- function(PGA, Tn, VS30, Vref, model) {
  .checkSiteCondition(VS30 = VS30, Vref = Vref, model = model)
  if (!(Vref %in% c(760, 3000))) stop("Vref must be 760 or 3000 m/s")
  if (!is.numeric(PGA) || !length(PGA) || any(!is.finite(PGA)) || any(PGA < 0)) {
    stop("PGA must be finite and non-negative")
  }
  if (!is.numeric(Tn) || length(Tn) != 1L || !is.finite(Tn) || Tn < 0) {
    stop("Tn must be one finite non-negative period")
  }
  if (VS30 == Vref) return(data.table(muLnF = rep(0, length(PGA)), sigmaLnF = 0, ID = model))
  if (model == "ST20") return(.siteFactorST20(PGA, Tn, VS30, Vref))
  if (model == "ST17") return(.siteFactorST17(PGA, Tn, VS30, Vref))
  stop("model must be \"ST20\" or \"ST17\"")
}

# The model and the two velocities shared by every site-response
# calculation; the reference velocities a model supports are checked by
# evaluateSiteFactor(), so that an identity (VS30 == Vref) holds for any reference.
.checkSiteCondition <- function(VS30, Vref, model) {
  if (!is.character(model) || length(model) != 1L || !(model %in% c("ST20", "ST17"))) {
    stop("model must be \"ST20\" or \"ST17\"")
  }
  if (!is.numeric(VS30) || length(VS30) != 1L || !is.finite(VS30) || VS30 <= 0) {
    stop("VS30 must be one finite positive velocity")
  }
  if (!is.numeric(Vref) || length(Vref) != 1L || !is.finite(Vref) || Vref <= 0) {
    stop("Vref must be one finite positive velocity")
  }
  invisible(TRUE)
}

.siteFactorST20 <- function(PGA, Tn, VS30, Vref) {
  C <- .coefficientsST20(Tn)
  VrefNL <- if (Tn < 0.4) 760 else 3000
  PGA3000 <- if (Vref == 760) PGA / 2.275 else PGA
  FRk <- log((PGA3000 + C$f3) / C$f3)
  weight760 <- function(v) min(max(
    (0.767 - 0.1) * log(v / 400) / log(600 / 400) + 0.1, 0.1), 0.767)
  f760 <- function(v) weight760(v) * C$f760i + (1 - weight760(v)) * C$f760g
  fv <- function(v) {
    if (v >= 3000) return(-f760(v))
    if (v <= C$V1) return(C$c * log(C$V1 / 760))
    if (v <= C$V2) return(C$c * log(v / 760))
    if (v <= 2000) return(C$c * log(C$V2 / 760))
    C$c * log(C$V2 / 760) - (C$c * log(C$V2 / 760) + f760(v)) * log(v / 2000) / log(3000 / 2000)
  }
  fLin <- function(v) fv(v) + f760(v)
  f2 <- function(v) {
    if (v >= C$Vc) return(0)
    C$f4 * (exp(C$f5 * (min(v, VrefNL) - 360)) - exp(C$f5 * (VrefNL - 360)))
  }
  sigmaLin <- function(v) {
    if (v < C$Vf) {
      D <- C$sig_l - C$sig_vc
      x <- (v - 200) / (C$Vf - 200)
      return(C$sig_l - 2 * D * x + D * x^2)
    }
    if (v <= C$V2) return(C$sig_vc)
    if (v < 2000) return(C$sig_vc + (C$sig_u - C$sig_vc) * ((v - C$V2) / (2000 - C$V2))^2)
    if (v <= 3000) return(C$sig_u * (1 - log(v / 2000) / log(3000 / 2000)))
    0
  }
  sigmaF2 <- function(v) {
    if (v < 300) return(C$sig_c)
    if (v < 1000) return(C$sig_c - C$sig_c * log(v / 300) / log(1000 / 300))
    0
  }
  if (Vref == 760) {
    muLnF <- (fLin(VS30) - fLin(760)) + (f2(VS30) - f2(760)) * FRk
    sigmaLnF <- sqrt(sigmaLin(VS30)^2 + (sigmaF2(VS30) * FRk)^2)
    return(data.table(muLnF = muLnF, sigmaLnF = sigmaLnF, ID = "ST20"))
  }
  Sigma760 <- weight760(VS30) * C$f760is + (1 - weight760(VS30)) * C$f760gs
  muLnF <- fLin(VS30) + f2(VS30) * FRk
  sigmaLnF <- sqrt(Sigma760^2 + sigmaLin(VS30)^2 + (sigmaF2(VS30) * FRk)^2)
  data.table(muLnF = muLnF, sigmaLnF = sigmaLnF, ID = "ST20")
}

# SiteAmp2020 row at period Tn: linear in ln(Tn) between tabulated periods
# (OpenQuake CoeffsTable rule); Tn <= 0.01 s is the PGA row; Tn > 10 s
# clamps to the last row.
.coefficientsST20 <- function(Tn) {
  DT <- SiteAmp2020
  if (Tn > 10) {
    warning("evaluateSiteFactor: Tn = ", Tn, " s exceeds the published 10 s table of ST20; clamped to 10 s.",
            call. = FALSE)
    Tn <- 10
  }
  if (Tn <= 0.01) return(as.list(DT[1L]))
  i <- findInterval(Tn, DT$Tn)
  if (DT$Tn[i] == Tn) return(as.list(DT[i]))
  r <- (log(Tn) - log(DT$Tn[i])) / (log(DT$Tn[i + 1L]) - log(DT$Tn[i]))
  as.list(DT[i] + r * (DT[i + 1L] - DT[i]))
}

.siteFactorST17 <- function(PGA, Tn, VS30, Vref) {
  # Coefficient tables of PEER 2017/04 (linear part, F760) and 2017/05
  # (nonlinear part, Table 2.1), interpolated linearly in Tn.
  at <- function(x, y) stats::approxfun(x = x, y = y, ties = mean)(Tn)
  TnLinear <- c(0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80, 1.00, 2.00, 3.00, 4.00, 5.00, 20.00)
  TnF760 <- c(0.00, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017,
              0.019, 0.02, 0.022, 0.023, 0.025, 0.027, 0.028, 0.031,
              0.033, 0.035, 0.038, 0.04, 0.043, 0.046, 0.05, 0.053,
              0.057, 0.061, 0.066, 0.071, 0.076, 0.081, 0.087, 0.093,
              0.1, 0.107, 0.115, 0.123, 0.132, 0.142, 0.152, 0.163,
              0.175, 0.187, 0.201, 0.215, 0.231, 0.248, 0.266, 0.285,
              0.305, 0.327, 0.351, 0.376, 0.404, 0.433, 0.464, 0.498,
              0.534, 0.572, 0.614, 0.658, 0.705, 0.756, 0.811, 0.87,
              0.933, 1, 1.072, 1.15, 1.233, 1.322, 1.417, 1.52, 1.63,
              1.748, 1.874, 2.009, 2.154, 2.31, 2.477, 2.656, 2.848,
              3.054, 3.275, 3.511, 3.765, 4.037, 4.329, 4.642, 4.977,
              5.337, 5.722, 6.136, 6.579, 7.055, 7.565, 8.111, 8.697,
              9.326, 10, 20)
  TnNonlinear <- c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.)
  C <- list(
    c = at(TnLinear, c(-0.359, -0.359, -0.344, -0.408, -0.416, -0.418, -0.446, -0.527,
                       -0.554, -0.580, -0.593, -0.597, -0.582, -0.582)),
    V1 = at(TnLinear, c(333, 333, 319, 314, 200, 200, 200, 230, 278, 286, 313, 322, 325, 325)),
    V2 = at(TnLinear, c(760, 760, 760, 760, 760, 760, 764, 942, 1103, 1201, 876, 881, 855, 855)),
    Vf = at(TnLinear, c(333, 333, 319, 314, 250, 250, 280, 280, 300, 300, 313, 322, 325, 325)),
    sigmaVc = at(TnLinear, c(0.268, 0.268, 0.270, 0.251, 0.225, 0.225, 0.225, 0.225,
                             0.225, 0.259, 0.306, 0.340, 0.340, 0.340)),
    sigmaL = at(TnLinear, c(0.281, 0.281, 0.263, 0.306, 0.276, 0.275, 0.311, 0.334,
                            0.377, 0.548, 0.538, 0.435, 0.400, 0.400)),
    sigmaU = at(TnLinear, c(0.472, 0.472, 0.470, 0.334, 0.381, 0.381, 0.323, 0.308,
                            0.361, 0.388, 0.551, 0.585, 0.587, 0.587)),
    F760 = at(TnF760, c(0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185,
                        0.185, 0.185, 0.185, 0.185, 0.189, 0.195, 0.203, 0.212,
                        0.224, 0.238, 0.252, 0.267, 0.283, 0.3, 0.318, 0.337,
                        0.356, 0.377, 0.4, 0.425, 0.454, 0.475, 0.512, 0.558,
                        0.613, 0.674, 0.73, 0.76, 0.759, 0.714, 0.647, 0.586,
                        0.534, 0.488, 0.449, 0.419, 0.39, 0.362, 0.332, 0.301,
                        0.278, 0.27, 0.262, 0.242, 0.224, 0.209, 0.197, 0.186,
                        0.175, 0.166, 0.157, 0.15, 0.142, 0.135, 0.127, 0.12,
                        0.111, 0.103, 0.095, 0.088, 0.083, 0.08, 0.078, 0.078,
                        0.077, 0.077, 0.078, 0.079, 0.079, 0.078, 0.076, 0.075,
                        0.074, 0.073, 0.073, 0.072, 0.07, 0.068, 0.066, 0.065,
                        0.065, 0.064, 0.063, 0.061, 0.06, 0.058, 0.057, 0.056,
                        0.056, 0.055, 0.055, 0.053, 0.053)),
    sigma760 = at(TnF760, c(0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434,
                            0.434, 0.434, 0.434, 0.432, 0.429, 0.422, 0.414, 0.404, 0.393,
                            0.387, 0.387, 0.39, 0.39, 0.381, 0.363, 0.34, 0.32, 0.308, 0.306,
                            0.312, 0.322, 0.335, 0.346, 0.357, 0.366, 0.372, 0.37, 0.351,
                            0.323, 0.284, 0.253, 0.234, 0.222, 0.214, 0.214, 0.207, 0.195,
                            0.177, 0.156, 0.141, 0.131, 0.124, 0.117, 0.115, 0.112, 0.113,
                            0.111, 0.105, 0.104, 0.113, 0.123, 0.132, 0.139, 0.138, 0.133,
                            0.13, 0.128, 0.124, 0.118, 0.112, 0.108, 0.11, 0.114, 0.12,
                            0.122, 0.124, 0.124, 0.118, 0.113, 0.112, 0.111, 0.109, 0.108,
                            0.111, 0.116, 0.12, 0.122, 0.12, 0.116, 0.112, 0.108, 0.104,
                            0.1, 0.096, 0.091, 0.087, 0.082, 0.078, 0.074, 0.07, 0.069, 0.069)),
    f3 = at(TnNonlinear, c(0.16249, 0.16249, 0.15083, 0.12815, 0.1307, 0.09414, 0.09888, 0.07357,
                           0.04367, 0.00164, 0.00746, 0.00269, 0.00242, 0.05329, 0.05329)),
    f4 = at(TnNonlinear, c(-0.50667, -0.50667, -0.44661, -0.30481, -0.22825, -0.11591, -0.07793,
                           -0.01592, -0.00478, -0.00236, -0.00626, -0.00331, -0.00256, -0.00631,
                           -0.00631)),
    f5 = at(TnNonlinear, c(-0.00273, -0.00273, -0.00335, -0.00488, -0.00655, -0.00872, -0.01028,
                           -0.01515, -0.01823, -0.01296, -0.01043, -0.01215, -0.01325,
                           -0.01403, -0.01403)),
    Vc = at(TnNonlinear, c(2990, 2990, 2990, 1533, 1152, 1018, 938, 832, 951, 879, 894, 875,
                           856, 837, 837)),
    sigmaC = at(TnNonlinear, c(0.12, 0.12, 0.12, 0.12, 0.15, 0.15, 0.15, 0.1, 0.06, 0.04, 0.04, 0.03,
                               0.02, 0.02, 0.02))
  )
  # Bounds of the quadratic interpolation of the linear dispersion.
  Vl <- 200
  Vu <- 2000
  muLnPGA <- log(PGA)
  # PGA at 760 m/s over PGA at 3000 m/s (PEER 2017/05, Eq. 2.4).
  RatioPGA <- if (Vref == 760) 2.275 else 1.0
  muLinear <- function(v) {
    if (v <= C$V1) return(C$c * log(C$V1 / 760))
    if (v <= C$V2) return(C$c * log(v / 760))
    C$c * log(C$V2 / 760) + C$c / 2 * log(v / C$V2)
  }
  a <- C$f4 * (exp(C$f5 * (min(VS30, 3000) - 360)) - exp(C$f5 * (3000 - 360)))
  b <- RatioPGA * C$f3
  muNL <- if (VS30 < C$Vc) a * log(exp(muLnPGA) / b + 1) else 0
  muI <- if (Vref == 3000) C$F760 else 0
  muLnF <- muLinear(VS30) + muI + muNL
  sigmaLinear <- function(v) {
    if (v <= C$Vf) {
      return(C$sigmaL - 2 * (C$sigmaL - C$sigmaVc) * (v - Vl) / (C$Vf - Vl) +
               (C$sigmaL - C$sigmaVc) * (v - Vl)^2 / (C$Vf - Vl)^2)
    }
    if (v <= C$V2) return(C$sigmaVc)
    C$sigmaVc + (C$sigmaU - C$sigmaVc) * (v - C$V2)^2 / (Vu - C$V2)^2
  }
  sigmaNL <- function(v) {
    if (v < 300) return(C$sigmaC)
    if (v < 1000) return(-C$sigmaC * log(v / 300) / log(1000 / 300) + C$sigmaC)
    0
  }
  # The nonlinear dispersion takes the PGA-dependent multiplier of the mean
  # (PEER 2017/05, Eq. 2.5).
  FRk <- log(exp(muLnPGA) / b + 1)
  sigmaI <- if (Vref == 3000) C$sigma760 else 0
  sigmaLnF <- sqrt(sigmaLinear(VS30)^2 + sigmaI^2 + (sigmaNL(VS30) * FRk)^2)
  data.table(muLnF = muLnF, sigmaLnF = sigmaLnF, ID = "ST17")
}

# Names resolved by data.table inside the site-response family.
utils::globalVariables(c(
  "Tn", "p", "Sa", "SaF", "k", "muLn", "sigmaLn", "w", "level",
  "muLn.PGA", "sigmaLn.PGA", "i.muLn.PGA", "i.sigmaLn.PGA", "weightSum", "Sa.x", "Sa.y",
  "ID", "siteID", "TR", "Vs30", "Vref", "IDg", "IDm", "VS30", "PGA"
))
