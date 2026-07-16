# nolint start


# ============================================================
# Site-factor model SaF_ST17  (canonical)
# ============================================================
#' Site-amplified spectral acceleration (2017 NGA-East generation)
#'
#' Computes site-amplified spectral acceleration Sa_F using the
#' 2017-generation NGA-East non-linear site factor model (PEER Report
#' 2017/05 lineage; see \code{\link{F_ST17}}). This wrapper calls the
#' underlying F_ST17() function.
#'
#' @param Sa numeric vector - rock spectral accelerations at period Tn (g)
#' @param pga numeric vector - peak ground acceleration on rock (g), same length as Sa
#' @param Tn numeric scalar - oscillator period (s)
#' @param vs30 numeric scalar - time-averaged shear-wave velocity to 30 m (m/s)
#' @param vref numeric scalar - reference velocity (m/s), default 760
#'
#' @return A \link[data.table]{data.table} with columns:
#'   \describe{
#'     \item{\code{muLnSaF}}{natural-log mean of amplified Sa}
#'     \item{\code{sdLnSaF}}{natural-log standard deviation of amplified Sa}
#'     \item{\code{ID}}{character string "ST17"}}
#'
#' @section References:
#' Hashash, Y.M.A. et al. (2017) "Recommendation for Ergodic Nonlinear
#' Site Amplification in Central and Eastern North America." PEER
#' Report 2017/05, Pacific Earthquake Engineering Research Center.
#'
#' @seealso \code{\link{F_ST17}} for the underlying site factor model
#' @seealso \code{\link{fitSaF}} for Monte-Carlo site amplification analysis
#'
#' @export
#'
#' @examples
#' library(data.table)
#'
#' # Site amplification for soft soil site
#' Sa_rock <- c(0.1, 0.2, 0.3, 0.4)
#' PGA_rock <- c(0.1, 0.1, 0.1, 0.1)
#'
#' result <- SaF_ST17(
#'   Sa = Sa_rock,
#'   pga = PGA_rock,
#'   Tn = 0.5,      # 0.5s period
#'   vs30 = 300,    # Soft soil
#'   vref = 760     # Rock reference
#' )
#'
#' print(result)
SaF_ST17 <- function(Sa, pga, Tn, vs30, vref = 760) {
    stopifnot(length(Sa) == length(pga))

    ## deterministic (rock) case --------------------------------------------
    if (vs30 == vref) {
        return(data.table(
            muLnSaF = log(Sa),
            sdLnSaF = 0,
            ID      = "ST17"
        ))
    }

    ## call the F-model - note the argument name **PGA**
    Fmodel <- F_ST17(PGA = pga, Tn = Tn, vs30 = vs30, vref = vref)

    data.table(
        muLnSaF = log(Sa) + Fmodel$muLnF,
        sdLnSaF = Fmodel$sdLnF,
        ID      = "ST17"
    )
}



# ============================================================
#  NGA-East 2017-generation site-factor model (PEER 2017/05)  -  F_ST17()
# ============================================================
# nolint start
#' Non-linear site amplification factor **F<sub>ST17</sub>**
#'
#' Implements the period-dependent, Vs30-dependent NGA-East site
#' amplification recommendations of the 2017 PEER-report generation
#' (nonlinear block transcribed from PEER Report 2017/05, Table 2.1;
#' the historical "ST17" id is kept for backward compatibility) for the
#' natural-log mean and standard deviation of the amplification factor
#' \eqn{F(T_n,PGA,V_{s30})}. The nonlinear sigma applies the
#' PGA-dependent multiplier of the source model,
#' \eqn{\sigma_{F_{nl}} = \sigma_{f2}\,\ln((PGA_{3000}+f_3)/f_3)}
#' (corrected 2026-07; earlier versions used bare \eqn{\sigma_{f2}},
#' making \code{sdLnF} amplitude-independent). For the final 2020
#' model generation see \code{\link{F_ST20}}.
#'
#' *If \code{vs30 == vref}, the site factor equals 1 deterministically, so
#'   the function returns \code{muLnF = 0} and \code{sdLnF = 0} (zero
#'   dispersion; \code{rnorm(n, 0, 0)} returns zeros and \code{ln F = 0}).*
#'
#' @param PGA   numeric **vector** - peak ground acceleration (in g).
#' @param Tn    numeric **scalar** - oscillator period (s).
#' @param vs30  numeric **scalar** - time-averaged shear-wave velocity to 30 m (m/s).
#' @param vref  numeric **scalar** - reference velocity (760 m/s or 3000 m/s), default 760.
#' @param Vl    numeric scalar - lower Vs30 bound for NL sigma interpolation (default 200).
#' @param Vu    numeric scalar - upper Vs30 bound for NL sigma interpolation (default 2000).
#'
#' @return A \link[data.table]{data.table} with columns:
#'   \describe{
#'     \item{\code{muLnF}}{natural-log mean of \eqn{F}}
#'     \item{\code{sdLnF}}{natural-log standard deviation of \eqn{F}
#'                         (set to \code{0} when \code{vs30 == vref})}
#'     \item{\code{ID}}{character string \code{"ST17"}}}
#'
#' @section References:
#' Stewart, J.P., Parker, G.A., Harmon, J.A., Atkinson, G.M., Boore,
#' D.M., Darragh, R.B., Silva, W.J. & Hashash, Y.M.A. (2017) "Expert
#' Panel Recommendations for Ergodic Site Amplification in Central and
#' Eastern North America." PEER Report 2017/04 (linear model and F760
#' tables, transcribed from its electronic supplement).
#'
#' Hashash, Y.M.A., Harmon, J., Ilhan, O., Stewart, J.P., Rathje, E.M.,
#' Campbell, K.W., Silva, W.J. & Xu, B. (2017) "Recommendation for
#' Ergodic Nonlinear Site Amplification in Central and Eastern North
#' America." PEER Report 2017/05 (nonlinear model, Table 2.1).
#'
#' @export

F_ST17 <- function(PGA, Tn, vs30, vref = 760, Vl = 200, Vu = 2000) {
    # nolint end
    ## -- basic input checks -----------------------------------------------
    OK <- length(vref) == 1L && length(vs30) == 1L && vref %in% c(760, 3000)
    stopifnot(OK)

    ## -- deterministic case (vs30 == vref) --------------------------------
    if (vs30 == vref) {
        return(data.table(
            muLnF = 0,
            sdLnF = 0,
            ID    = "ST17"
        ))
    }

    ## ---------- interpolated coefficient tables (unchanged) ---------------
    cI <- stats::approxfun(
        x = c(
            0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
            1.00, 2.00, 3.00, 4.00, 5.00, 20.00
        ),
        y = c(
            -0.359, -0.359, -0.344, -0.408, -0.416, -0.418, -0.446, -0.527,
            -0.554, -0.580, -0.593, -0.597, -0.582, -0.582
        ),
        ties = mean
    )
    V1I <- stats::approxfun(
        x = c(
            0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
            1.00, 2.00, 3.00, 4.00, 5.00, 20.00
        ),
        y = c(
            333, 333, 319, 314, 200, 200, 200, 230,
            278, 286, 313, 322, 325, 325
        ),
        ties = mean
    )
    V2I <- stats::approxfun(
        x = c(
            0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
            1.00, 2.00, 3.00, 4.00, 5.00, 20.00
        ),
        y = c(
            760, 760, 760, 760, 760, 760, 764, 942,
            1103, 1201, 876, 881, 855, 855
        ),
        ties = mean
    )
    VfI <- stats::approxfun(
        x = c(
            0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
            1.00, 2.00, 3.00, 4.00, 5.00, 20.00
        ),
        y = c(
            333, 333, 319, 314, 250, 250, 280, 280,
            300, 300, 313, 322, 325, 325
        ),
        ties = mean
    )
    sdVcI <- stats::approxfun(
        x = c(
            0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
            1.00, 2.00, 3.00, 4.00, 5.00, 20.00
        ),
        y = c(
            0.268, 0.268, 0.270, 0.251, 0.225, 0.225, 0.225, 0.225,
            0.225, 0.259, 0.306, 0.340, 0.340, 0.340
        ),
        ties = mean
    )
    sdLI <- stats::approxfun(
        x = c(
            0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
            1.00, 2.00, 3.00, 4.00, 5.00, 20.00
        ),
        y = c(
            0.281, 0.281, 0.263, 0.306, 0.276, 0.275, 0.311, 0.334,
            0.377, 0.548, 0.538, 0.435, 0.400, 0.400
        ),
        ties = mean
    )
    sdUI <- stats::approxfun(
        x = c(
            0.00, 0.08, 0.10, 0.20, 0.30, 0.40, 0.50, 0.80,
            1.00, 2.00, 3.00, 4.00, 5.00, 20.00
        ),
        y = c(
            0.472, 0.472, 0.470, 0.334, 0.381, 0.381, 0.323, 0.308,
            0.361, 0.388, 0.551, 0.585, 0.587, 0.587
        ),
        ties = mean
    )
    F760I <- stats::approxfun(
        x = c(
            0.00, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017,
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
            9.326, 10, 20
        ),
        y = c(
            0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185, 0.185,
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
            0.056, 0.055, 0.055, 0.053, 0.053
        ),
        ties = mean
    )
    sd760I <- stats::approxfun(
        x = c(
            0.00, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017,
            0.019, 0.02, 0.022, 0.023, 0.025, 0.027, 0.028, 0.031, 0.033,
            0.035, 0.038, 0.04, 0.043, 0.046, 0.05, 0.053, 0.057, 0.061,
            0.066, 0.071, 0.076, 0.081, 0.087, 0.093, 0.1, 0.107, 0.115,
            0.123, 0.132, 0.142, 0.152, 0.163, 0.175, 0.187, 0.201, 0.215,
            0.231, 0.248, 0.266, 0.285, 0.305, 0.327, 0.351, 0.376, 0.404,
            0.433, 0.464, 0.498, 0.534, 0.572, 0.614, 0.658, 0.705, 0.756,
            0.811, 0.87, 0.933, 1, 1.072, 1.15, 1.233, 1.322, 1.417, 1.52,
            1.63, 1.748, 1.874, 2.009, 2.154, 2.31, 2.477, 2.656, 2.848,
            3.054, 3.275, 3.511, 3.765, 4.037, 4.329, 4.642, 4.977, 5.337,
            5.722, 6.136, 6.579, 7.055, 7.565, 8.111, 8.697, 9.326, 10, 20
        ),
        y = c(
            0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434, 0.434,
            0.434, 0.434, 0.434, 0.432, 0.429, 0.422, 0.414, 0.404, 0.393,
            0.387, 0.387, 0.39, 0.39, 0.381, 0.363, 0.34, 0.32, 0.308, 0.306,
            0.312, 0.322, 0.335, 0.346, 0.357, 0.366, 0.372, 0.37, 0.351,
            0.323, 0.284, 0.253, 0.234, 0.222, 0.214, 0.214, 0.207, 0.195,
            0.177, 0.156, 0.141, 0.131, 0.124, 0.117, 0.115, 0.112, 0.113,
            0.111, 0.105, 0.104, 0.113, 0.123, 0.132, 0.139, 0.138, 0.133,
            0.13, 0.128, 0.124, 0.118, 0.112, 0.108, 0.11, 0.114, 0.12,
            0.122, 0.124, 0.124, 0.118, 0.113, 0.112, 0.111, 0.109, 0.108,
            0.111, 0.116, 0.12, 0.122, 0.12, 0.116, 0.112, 0.108, 0.104,
            0.1, 0.096, 0.091, 0.087, 0.082, 0.078, 0.074, 0.07, 0.069, 0.069
        ),
        ties = mean
    )
    f3I <- stats::approxfun(
        x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
        y = c(
            0.16249, 0.16249, 0.15083, 0.12815, 0.1307, 0.09414, 0.09888, 0.07357,
            0.04367, 0.00164, 0.00746, 0.00269, 0.00242, 0.05329, 0.05329
        ),
        ties = mean
    )
    f4I <- stats::approxfun(
        x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
        y = c(
            -0.50667, -0.50667, -0.44661, -0.30481, -0.22825, -0.11591, -0.07793,
            -0.01592, -0.00478, -0.00236, -0.00626, -0.00331, -0.00256, -0.00631,
            -0.00631
        ),
        ties = mean
    )
    f5I <- stats::approxfun(
        x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
        y = c(
            -0.00273, -0.00273, -0.00335, -0.00488, -0.00655, -0.00872, -0.01028,
            -0.01515, -0.01823, -0.01296, -0.01043, -0.01215, -0.01325,
            -0.01403, -0.01403
        ),
        ties = mean
    )
    VcI <- stats::approxfun(
        x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
        y = c(
            2990, 2990, 2990, 1533, 1152, 1018, 938, 832, 951, 879, 894, 875,
            856, 837, 837
        ),
        ties = mean
    )
    sdcI <- stats::approxfun(
        x = c(0.0, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1., 2., 3., 4., 5., 10., 20.),
        y = c(
            0.12, 0.12, 0.12, 0.12, 0.15, 0.15, 0.15, 0.1, 0.06, 0.04, 0.04, 0.03,
            0.02, 0.02, 0.02
        ),
        ties = mean
    )

    ## ---------- mean term --------------------------------------------------
    muLnPGA <- log(PGA)
    C7603000 <- if (vref == 760) 2.275 else 1.0

    if (vs30 <= V1I(Tn)) {
        muL <- cI(Tn) * log(V1I(Tn) / 760)
    } else if (vs30 <= V2I(Tn)) {
        muL <- cI(Tn) * log(vs30 / 760)
    } else {
        muL <- cI(Tn) * log(V2I(Tn) / 760) + cI(Tn) / 2 * log(vs30 / V2I(Tn))
    }

    a <- f4I(Tn) * (exp(f5I(Tn) * (min(vs30, 3000) - 360)) -
        exp(f5I(Tn) * (3000 - 360)))
    b <- C7603000 * f3I(Tn)

    muNL <- if (vs30 < VcI(Tn)) a * log(exp(muLnPGA) / b + 1) else 0
    muI <- if (vref == 3000) F760I(Tn) else 0

    muLnF <- muL + muI + muNL

    ## ---------- sigma term -------------------------------------------------
    if (vs30 <= VfI(Tn)) {
        sdL <- sdLI(Tn) - 2 * (sdLI(Tn) - sdVcI(Tn)) *
            (vs30 - Vl) / (VfI(Tn) - Vl) +
            (sdLI(Tn) - sdVcI(Tn)) * (vs30 - Vl)^2 / (VfI(Tn) - Vl)^2
    } else if (vs30 <= V2I(Tn)) {
        sdL <- sdVcI(Tn)
    } else {
        sdL <- sdVcI(Tn) + (sdUI(Tn) - sdVcI(Tn)) *
            (vs30 - V2I(Tn))^2 / (Vu - V2I(Tn))^2
    }

    if (vs30 < 300) {
        sdNL <- sdcI(Tn)
    } else if (vs30 < 1000) {
        sdNL <- -sdcI(Tn) * log(vs30 / 300) / log(1000 / 300) + sdcI(Tn)
    } else {
        sdNL <- 0
    }

    ## sigma_Fnl = sigma_f2 * ln((PGA_3000 + f3)/f3): PEER 2017/05 Eq. 2.5
    ## gives sigma_f2 (slope sigma); the Fnl sigma requires the same
    ## PGA-dependent multiplier as muNL (cf. USGS nshmp-lib, OQ nga_east.py).
    FRk <- log(exp(muLnPGA) / b + 1)
    sdI <- if (vref == 3000) sd760I(Tn) else 0
    sdLnF <- sqrt(sdL^2 + sdI^2 + (sdNL * FRk)^2)

    ## ---------- output -----------------------------------------------------
    data.table(
        muLnF = muLnF,
        sdLnF = sdLnF,
        ID    = "ST17"
    )
}

# nolint end


# ============================================================
#  NGA-East ergodic site amplification (2020) - F_ST20()
# ============================================================
#' Non-linear site amplification factor **F<sub>ST20</sub>**
#'
#' Final-generation NGA-East ergodic site-amplification model: linear
#' component and F760 impedance/gradient weighting per Stewart et al.
#' (2020), nonlinear component per Hashash et al. (2020), with the
#' period-dependent nonlinear reference velocity (760 m/s for
#' Tn < 0.4 s, 3000 m/s above) and the PGA-dependent nonlinear sigma
#' \eqn{\sigma_{F_{nl}} = \sigma_{f2}\,\ln((PGA_{3000}+f_3)/f_3)}.
#' Coefficients are the USGS nshmp-lib table (see
#' \code{data-raw/SiteAmp2020.R} for provenance); reference
#' implementation cross-checked function-by-function against OpenQuake
#' \code{gsim/nga_east.py}.
#'
#' Conventions:
#' * The published model is referenced to 3000 m/s. For
#'   \code{vref = 760} the mean is the exact difference
#'   \eqn{F(V_{s30}) - F(760)} of the composite model (including the
#'   F760 weighting difference below 600 m/s and the small
#'   \eqn{F_{nl}(760)} term at Tn >= 0.4 s), and the input \code{PGA}
#'   (at 760 m/s) is converted to the 3000 m/s anchor via
#'   \eqn{PGA_{3000} = PGA_{760}/2.275} (PEER 2017/05 Eq. 2.4).
#' * \code{sdLnF} is evaluated at the target \code{vs30} only; for
#'   \code{vref = 760} no reference-side sigma is added (the rock
#'   hazard carries its own dispersion). For \code{vref = 3000} the
#'   weighted F760 sigma enters the SRSS.
#' * \code{vs30 == vref} returns \code{muLnF = 0, sdLnF = 0}
#'   (deterministic unit factor, same contract as \code{F_ST17}).
#' * Coefficients interpolate linearly in \eqn{\ln T_n} (OpenQuake
#'   CoeffsTable rule); \code{Tn = 0} uses the PGA row; \code{Tn} above
#'   the published 10 s limit is clamped to 10 s with a warning.
#'
#' @param PGA   numeric vector - peak ground acceleration (g) at the
#'   \code{vref} condition.
#' @param Tn    numeric scalar - oscillator period (s).
#' @param vs30  numeric scalar - target Vs30 (m/s).
#' @param vref  numeric scalar - reference velocity, 760 or 3000 m/s.
#'
#' @return A \link[data.table]{data.table} with columns \code{muLnF},
#'   \code{sdLnF}, \code{ID = "ST20"}.
#'
#' @section References:
#' Stewart, J.P. et al. (2020) "Ergodic site amplification model for
#' central and eastern North America." *Earthquake Spectra* **36**(1).
#' Hashash, Y.M.A. et al. (2020) "Nonlinear site amplification model
#' for ergodic seismic hazard analysis in central and eastern North
#' America." *Earthquake Spectra* **36**(1).
#'
#' @examples
#' # soft-soil amplification of the 0.2 s ordinate, hazard at 760 m/s
#' F_ST20(PGA = c(0.05, 0.2, 0.5), Tn = 0.2, vs30 = 300)
#'
#' # very-hard-rock reference (3000 m/s)
#' F_ST20(PGA = 0.1, Tn = 1.0, vs30 = 760, vref = 3000)
#' @export
F_ST20 <- function(PGA, Tn, vs30, vref = 760) {
    OK <- length(Tn) == 1L && length(vs30) == 1L && length(vref) == 1L &&
        vref %in% c(760, 3000)
    stopifnot(OK)
    if (vs30 == vref) {
        return(data.table(muLnF = 0, sdLnF = 0, ID = "ST20"))
    }

    C <- .st20Coef(Tn)
    VREFNL <- if (Tn < 0.4) 760 else 3000
    PGA3 <- if (vref == 760) PGA / 2.275 else PGA
    FRk <- log((PGA3 + C$f3) / C$f3)

    ## F760 impedance/gradient weight, Stewart et al. (2020) eqs 5-7
    .w <- function(v) min(max(
        (0.767 - 0.1) * log(v / 400) / log(600 / 400) + 0.1, 0.1), 0.767)
    .f760 <- function(v) .w(v) * C$f760i + (1 - .w(v)) * C$f760g

    ## linear Vs30 scaling, Stewart et al. (2020) eq 3 (3000-referenced
    ## composite flin = fv + f760; fv branches as in OQ _get_fv)
    .flin <- function(v) {
        X <- if (v <= C$V1) C$c * log(C$V1 / 760)
        else if (v <= C$V2) C$c * log(v / 760)
        else if (v <= 2000) C$c * log(C$V2 / 760)
        else C$c * log(C$V2 / 760) -
            (C$c * log(C$V2 / 760) + .f760(v)) * log(v / 2000) / log(3000 / 2000)
        if (v >= 3000) X <- -.f760(v)
        X + .f760(v)
    }

    ## nonlinear slope f2, Hashash et al. (2020) eq 3
    .f2 <- function(v) {
        if (v >= C$Vc) return(0)
        C$f4 * (exp(C$f5 * (min(v, VREFNL) - 360)) -
                exp(C$f5 * (VREFNL - 360)))
    }

    ## linear sigma, Stewart et al. (2020) eq 4
    .sdLin <- function(v) {
        if (v < C$Vf) {
            D <- C$sig_l - C$sig_vc
            x <- (v - 200) / (C$Vf - 200)
            C$sig_l - 2 * D * x + D * x^2
        } else if (v <= C$V2) C$sig_vc
        else if (v < 2000) C$sig_vc +
            (C$sig_u - C$sig_vc) * ((v - C$V2) / (2000 - C$V2))^2
        else if (v <= 3000) C$sig_u * (1 - log(v / 2000) / log(3000 / 2000))
        else 0
    }

    ## nonlinear slope sigma, Hashash et al. (2020) eq 2.5 taper
    .sdF2 <- function(v) {
        if (v < 300) C$sig_c
        else if (v < 1000) C$sig_c - C$sig_c * log(v / 300) / log(1000 / 300)
        else 0
    }

    if (vref == 760) {
        muLnF <- (.flin(vs30) - .flin(760)) + (.f2(vs30) - .f2(760)) * FRk
        sdLnF <- sqrt(.sdLin(vs30)^2 + (.sdF2(vs30) * FRk)^2)
    } else {
        SD760 <- .w(vs30) * C$f760is + (1 - .w(vs30)) * C$f760gs
        muLnF <- .flin(vs30) + .f2(vs30) * FRk
        sdLnF <- sqrt(SD760^2 + .sdLin(vs30)^2 + (.sdF2(vs30) * FRk)^2)
    }

    data.table(muLnF = muLnF, sdLnF = sdLnF, ID = "ST20")
}


# Interpolated SiteAmp2020 coefficient row at period Tn: linear in
# ln(Tn) between tabulated periods (OpenQuake CoeffsTable logratio
# rule); Tn = 0 is the PGA row; Tn > 10 s clamps to the last row.
.st20Coef <- function(Tn) {
    DT <- SiteAmp2020
    if (Tn > 10) {
        warning("F_ST20: Tn = ", Tn,
                " s exceeds the published 10 s table; clamped to 10 s.",
                call. = FALSE)
        Tn <- 10
    }
    if (Tn <= 0.01) return(as.list(DT[1L]))
    i <- findInterval(Tn, DT$Tn)
    if (DT$Tn[i] == Tn) return(as.list(DT[i]))
    r <- (log(Tn) - log(DT$Tn[i])) / (log(DT$Tn[i + 1L]) - log(DT$Tn[i]))
    as.list(DT[i] + r * (DT[i + 1L] - DT[i]))
}
