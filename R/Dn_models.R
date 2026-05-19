# nolint start
#  Displacement models (empirical equations) --------------------------------
#' Ambraseys & Menu (1988) rigid sliding-block model
#' @param PGA Peak ground acceleration (g).
#' @param ky  Yield acceleration (g).
#' @return data.table(muLnD, sdLnD, ID = "AM88")
#' @examples
#' Dn_AM88(PGA = 0.5, ky = 0.1)
#' @export
Dn_AM88 <- function(PGA, ky) {
    r <- ky / PGA
    r <- pmin(r, 0.9999) # keep (1 - r) > 0
    muLog10D <- 0.90 + log10((1 - r)^2.53 * r^-1.09)
    sdLog10D <- 0.30

    data.table(
        muLnD = muLog10D * log(10),
        sdLnD = sdLog10D * log(10),
        ID = "AM88"
    )
}


#' Jibson (2007) empirical model
#' @inheritParams Dn_AM88
#' @param AI Arias Intensity (m/s); if NULL, back-calculated from PGA.
#' @return data.table(muLnD, sdLnD, ID = "JB07")
#' @examples
#' Dn_JB07(PGA = 0.5, ky = 0.1)
#' @export
Dn_JB07 <- function(PGA, ky, AI = NULL) {
    if (is.null(AI)) AI <- (PGA^1.9228) * exp(2.6109)
    r <- ky / PGA
    muLog10D <- 0.561 * log10(AI) - 3.833 * log10(r) - 1.474
    sdLog10D <- 0.616

    data.table(
        muLnD = muLog10D * log(10),
        sdLnD = sdLog10D * log(10),
        ID = "JB07"
    )
}

#' Saygili & Rathje (2008) model
#' @inheritParams Dn_JB07
#' @return data.table(muLnD, sdLnD, ID = "SR08")
#' @examples
#' Dn_SR08(PGA = 0.5, ky = 0.1)
#' @export
Dn_SR08 <- function(PGA, ky, AI = NULL) {
    if (is.null(AI)) AI <- (PGA^1.9228) * exp(2.6109)
    r <- ky / PGA
    muLnD <- 2.39 - 5.24 * r - 18.78 * r^2 + 42.01 * r^3 - 29.15 * r^4 -
        1.56 * log(PGA) + 1.38 * log(AI)
    sdLnD <- 0.46 + 0.56 * r

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "SR08"
    )
}

#' Bray & Travasarou (2007) flexible sliding-block model
#' @param ky Yield acceleration (g).
#' @param Sa Spectral acceleration at 1.5 * Ts (g).
#' @param Ts Fundamental period (s).
#' @param Mw Moment magnitude.
#' @return data.table(muLnD, sdLnD, ID = "BT07")
#' @examples
#' Dn_BT07(ky = 0.1, Sa = 0.8, Ts = 0.3, Mw = 7.0)
#' @export
Dn_BT07 <- function(ky, Sa, Ts, Mw = 6.5) {
    lnSa <- log(Sa)
    lnky <- log(ky)

    muLnD <- -1.10 - 2.83 * lnky - 0.333 * lnky^2 +
        0.566 * lnky * lnSa + 3.04 * lnSa -
        0.244 * lnSa^2 + 0.278 * (Mw - 7)

    sdLnD <- 0.66

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "BT07"
    )
}

#' Bray & Macedo (2017) - subduction/interface events
#' @param ky Yield acceleration (g).
#' @param Sa Spectral acceleration at 1.5 * Ts (g).
#' @param Ts Fundamental period (s).
#' @param Mw Moment magnitude (default 7.5).
#' @return data.table(muLnD, sdLnD, ID = "BM17")
#' @examples
#' Dn_BM17(ky = 0.1, Sa = 0.8, Ts = 0.3, Mw = 8.0)
#' @export
Dn_BM17 <- function(ky, Sa, Ts, Mw = 7.5) {
    lnSa <- log(Sa)

    if (Ts < 0.1) {
        c1 <- -5.864
        c2 <- -9.421
        c3 <- 0.0
    } else {
        c1 <- -6.896
        c2 <- 3.081
        c3 <- -0.803
    }

    a0 <- c1 + 0.550 * Mw + c2 * Ts + c3 * Ts^2 -
        3.353 * log(ky) - 0.390 * log(ky)^2
    a1 <- 3.060 + 0.538 * log(ky)
    a2 <- -0.225

    muLnD <- a0 + a1 * lnSa + a2 * lnSa^2
    sdLnD <- 0.73

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "BM17"
    )
}

#' Bray & Macedo (2019, corrected 2023) - shallow-crustal events
#'
#' Multi-equation structure per BM23 (Soil Dynamics and Earthquake Engineering
#' 168, 107835, DOI 10.1016/j.soildyn.2023.107835):
#'   Eq (10): ordinary motions, PGV <= 115 cm/s
#'   Eq (6):  D100 max-component, PGV > 115 cm/s
#'   Eq (7):  D50  median-component, PGV > 115 cm/s
#' with a sub-regime split at PGV = 150 cm/s within Eqs (6) and (7).
#'
#' @param ky   Yield acceleration (g).
#' @param Ts   Fundamental period (s).
#' @param Sa   Spectral acceleration at 1.3 * Ts (g).
#' @param PGA  Peak ground acceleration (g).
#' @param PGV  Peak ground velocity (cm/s).
#' @param Mw   Moment magnitude (default 6.5).
#' @param NFC  Near-fault component: "D100" (default; max-component, slope
#'   within 45 deg of fault-normal) or "D50" (median-component, other
#'   orientations). Used only when PGV > 115 cm/s.
#' @return data.table(muLnD, sdLnD, ID = "BM19")
#' @examples
#' Dn_BM19(ky = 0.1, Ts = 0.3, Sa = 0.8, PGA = 0.5, PGV = 80, Mw = 7.0)
#' @export
Dn_BM19 <- function(ky, Ts, Sa, PGA, PGV, Mw = 6.5, NFC = "D100") {
    n <- length(Sa)
    if (!(length(PGA) == n && length(PGV) == n)) {
        stop("Sa, PGA, and PGV must have the same length")
    }
    if (!NFC %in% c("D100", "D50")) {
        stop("NFC must be 'D100' or 'D50'")
    }

    lnSa  <- log(Sa)
    lnky  <- log(ky)
    lnPGV <- log(PGV)

    Ordinary  <- PGV <= 115
    PulseLow  <- PGV >  115 & PGV <= 150
    PulseHigh <- PGV >  150

    # Eq (10): ordinary motions — identical to the pre-fix implementation
    A0Ord <- if (Ts < 0.1) {
        -4.684 + 0.603 * Mw - 9.471 * Ts - 2.482 * lnky - 0.244 * lnky^2
    } else {
        -5.981 + 0.603 * Mw + 3.223 * Ts - 0.945 * Ts^2 - 2.482 * lnky - 0.244 * lnky^2
    }
    MuOrd <- A0Ord + (2.649 + 0.344 * lnky) * lnSa - 0.090 * lnSa^2

    # Pulse equation coefficients (Eq 6 D100 or Eq 7 D50)
    if (NFC == "D100") {
        # Eq (6) — BM23 coefficients (verified from paper text)
        if (Ts < 0.1) {
            C1Low <- -6.724; C1High <-  1.991; c2 <- -2.744; c3 <-  0.000
        } else {
            C1Low <- -6.951; C1High <-  1.764; c2 <-  1.069; c3 <- -0.498
        }
        C4Low <- 1.547; C4High <- -0.097
        b1 <- 2.632; b2 <- 0.278; b3 <- 0.527; b4 <- 1.978; b5 <- 0.233
        Bm <- 0.01; SdPulse <- 0.56
    } else {
        # Eq (7) — BM23 coefficients (verified from paper text)
        if (Ts < 0.1) {
            C1Low <- -7.497; C1High <-  2.480; c2 <- -2.731; c3 <-  0.000
        } else {
            C1Low <- -7.718; C1High <- -0.369; c2 <-  1.031; c3 <- -0.480
        }
        C4Low <- 1.458; C4High <- 0.025
        b1 <- 2.931; b2 <- 0.319; b3 <- 0.584; b4 <- 2.261; b5 <- 0.241
        Bm <- 0.05; SdPulse <- 0.54
    }

    C1Vec <- ifelse(PulseLow, C1Low, ifelse(PulseHigh, C1High, NA_real_))
    C4Vec <- ifelse(PulseLow, C4Low, ifelse(PulseHigh, C4High, NA_real_))

    MuPulse <- C1Vec - b1 * lnky - b2 * lnky^2 + b3 * lnky * lnSa +
               b4 * lnSa - b5 * lnSa^2 + c2 * Ts + c3 * Ts^2 +
               Bm * Mw + C4Vec * lnPGV

    muLnD <- ifelse(Ordinary, MuOrd, MuPulse)
    sdLnD <- ifelse(Ordinary, 0.72, SdPulse)

    data.table(
        muLnD = muLnD,
        sdLnD = sdLnD,
        ID = "BM19"
    )
}
