# nolint start


# ============================================================
# Site–factor model SaF_ST17  (canonical)
# ============================================================
#' Site-amplified spectral acceleration using Seyhan & Stewart (2014)
#'
#' Computes site-amplified spectral acceleration Sa_F using the Seyhan & Stewart (2014)
#' non-linear site factor model. This is the canonical implementation that calls
#' the underlying F_ST17() function.
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
#' Seyhan, E. & Stewart, J.P. (2014) "Semi-empirical nonlinear site amplification
#' from NGA-West2 data and simulations." *Earthquake Spectra* **30**(3):1241--1256.
#' \doi{10.1193/063013EQS181M}
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

    ## call the F-model — note the argument name **PGA**
    Fmodel <- F_ST17(PGA = pga, Tn = Tn, vs30 = vs30, vref = vref)

    data.table(
        muLnSaF = log(Sa) + Fmodel$muLnF,
        sdLnSaF = Fmodel$sdLnF,
        ID      = "ST17"
    )
}



# ============================================================
#  Seyhan & Stewart (2014) site-factor model  —  F_ST17()
# ============================================================
# nolint start
#' Non-linear site amplification factor **F<sub>ST17</sub>**
#'
#' Implements the period-dependent, Vs30-dependent model of Seyhan & Stewart (2014)
#' for the natural-log mean and standard deviation of the amplification factor
#' \eqn{F(T_n,PGA,V_{s30})}.
#'
#' *If \code{vs30 == vref}, the site factor equals 1 deterministically, so
#'   the function returns \code{muLnF = 0} and \code{sdLnF = 0} (zero
#'   dispersion; \code{rnorm(n, 0, 0)} returns zeros and \code{ln F = 0}).*
#'
#' @param PGA   numeric **vector** – peak ground acceleration (in g).
#' @param Tn    numeric **scalar** – oscillator period (s).
#' @param vs30  numeric **scalar** – time-averaged shear-wave velocity to 30 m (m/s).
#' @param vref  numeric **scalar** – reference velocity (760 m/s or 3000 m/s), default 760.
#' @param Vl    numeric scalar – lower Vs30 bound for NL sigma interpolation (default 200).
#' @param Vu    numeric scalar – upper Vs30 bound for NL sigma interpolation (default 2000).
#'
#' @return A \link[data.table]{data.table} with columns:
#'   \describe{
#'     \item{\code{muLnF}}{natural-log mean of \eqn{F}}
#'     \item{\code{sdLnF}}{natural-log standard deviation of \eqn{F}
#'                         (set to \code{0} when \code{vs30 == vref})}
#'     \item{\code{ID}}{character string \code{"ST17"}}}
#'
#' @section References:
#' Stewart, J.P. & Seyhan, E. (2017) “Semi-empirical nonlinear site
#' amplification model for global application.” *Earthquake Spectra* **33**(1).
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

    sdI <- if (vref == 3000) sd760I(Tn) else 0
    sdLnF <- sqrt(sdL^2 + sdI^2 + sdNL^2)

    ## ---------- output -----------------------------------------------------
    data.table(
        muLnF = muLnF,
        sdLnF = sdLnF,
        ID    = "ST17"
    )
}

# nolint end
