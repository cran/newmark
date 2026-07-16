#' MCER / Design elastic spectrum for a Vs30/Vref slice (ASCE 7-22 shape)
#'
#' Builds the two-period design/MCER spectrum shape of ASCE 7-22
#' Sec. 11.4.5.2 (Eqs. 11.4-3/-4/-5, DESIGN = 2/3 x MCER), with the
#' spectral parameters derived from the input site-specific spectrum
#' per ASCE 7-22 Sec. 21.4: at MCER level, SMS = 0.9 x max SaM(T) over
#' T from 0.2 to 5 s, and SM1 = max(0.9 x max of T x SaM(T), SaM(1 s))
#' over T from 1 to 2 s when vs30 > 442 m/s, else 1 to 5 s (the 1.5 and 2/3
#' factors of Eqs. 21.4/11.4-1/-2 cancel at MCER level). Evaluation
#' windows are clamped to the period range of the input. The
#' Sec. 21.2.3 80-percent floor and the Sec. 21.2.2 deterministic cap
#' are NOT applied (both require mapped/deterministic inputs this
#' function does not receive). See NEWS.md for the audit summary.
#'
#' The input `UHS` slice must contain at least `Tn` and `SaF`.
#' Typical usage passes a single TR and `p == "mean"`.
#'
#' @param UHS data.table slice from UHSTable (single TR, p == "mean").
#' @param TL  numeric, long-period transition (seconds). Default 8.
#' @param vs30 numeric scalar or NULL. Site Vs30 (m/s) selecting the
#'   Sec. 21.4 SD1 window: 1 to 2 s when vs30 > 442 m/s, else 1 to 5 s.
#'   NULL (default) uses the conservative 1 to 5 s window.
#' @param spectrum character, "MCER" (default) or "DESIGN" (2/3 * MCER).
#'
#' @return data.table with columns `Tn, SaF`.
#' @examples
#' uhs <- data.table::data.table(Tn = c(0, 0.2, 0.5, 1, 2), SaF = c(0.4, 1.0, 0.9, 0.6, 0.3))
#' designUHS(uhs)
#' @export
designUHS <- function(UHS, TL = 8, spectrum = c("MCER", "DESIGN"),
                      vs30 = NULL) {
  spectrum <- match.arg(spectrum)

  # Sort and extract raw vectors
  UHS <- UHS[order(Tn)]
  Tn0 <- as.numeric(UHS$Tn)
  Sa0 <- as.numeric(UHS$SaF)

  # Drop NA/Inf and collapse duplicates in Tn using the maximum SaF
  ok <- is.finite(Tn0) & is.finite(Sa0)
  Tn0 <- Tn0[ok]
  Sa0 <- Sa0[ok]
  if (length(Tn0) < 1L) {
    stop("designUHS(): no valid (Tn, SaF) after cleaning")
  }
  if (any(duplicated(Tn0))) {
    DTu <- data.table::data.table(Tn = Tn0, Sa = Sa0)[, .(Sa = max(Sa, na.rm = TRUE)), by = Tn]
    Tn0 <- DTu$Tn
    Sa0 <- DTu$Sa
  }
  if (length(Tn0) < 2L) {
    stop("designUHS(): need >= 2 (Tn, SaF) points after cleaning")
  }

  # Safe interpolator (linear in T; ASCE uses linear T-spacing)
  fSa <- stats::approxfun(Tn0, Sa0, rule = 2, ties = "mean")

  # ASCE 7-22 Sec. 21.4 parameters at MCER level (1.5 and 2/3 cancel):
  # SMS = 0.9 max SaM over [0.2, 5] s; SM1 = max(0.9 max T*SaM, SaM(1 s))
  # over [1, 2] s (vs30 > 442 m/s) or [1, 5] s; windows clamped to the
  # input period range. fSa is piecewise linear, so maxima are searched
  # on the knots plus a dense grid (T*Sa is quadratic within segments).
  .gmax <- function(lo, hi, f) {
    hi <- min(hi, max(Tn0))
    if (hi <= lo) return(f(lo))
    G <- sort(unique(c(lo, hi, Tn0[Tn0 >= lo & Tn0 <= hi],
                       exp(seq(log(lo), log(hi), length.out = 256L)))))
    max(f(G))
  }
  SMS <- 0.9 * .gmax(0.2, 5, fSa)
  TSTAR <- if (!is.null(vs30) && is.finite(vs30) && vs30 > 442) 2 else 5
  SM1 <- max(0.9 * .gmax(1, TSTAR, function(t) t * fSa(t)), fSa(1.0))
  if (!is.finite(SMS) || !is.finite(SM1)) {
    stop("designUHS(): invalid SMS or SM1 from the Sec. 21.4 derivation")
  }

  # Corner periods and evaluation grid
  Ts <- SM1 / SMS
  if (!is.finite(Ts) || Ts <= 0) stop("designUHS(): Ts must be positive")
  if (!is.finite(TL) || TL <= 0) stop("designUHS(): TL must be positive")

  T0 <- 0.2 * Ts
  TnDecay1 <- exp(seq(log(max(Ts, 1e-6)), log(TL),          length.out = 20))
  # sample the long-period branch (Eq. 11.4-5) one octave beyond TL;
  # the previous bound max(TL, 5) collapsed to TL whenever TL >= 5,
  # leaving the T > TL branch unreachable in the output
  TnDecay2 <- exp(seq(log(TL),               log(2 * TL),    length.out = 20))
  Tn <- sort(unique(c(0, T0, Ts, TL, TnDecay1, TnDecay2)))

  # MCER or DESIGN spectrum
  if (tolower(spectrum) == "mcer") {
    SaF <- piecewiseSa(t = Tn, sds = SMS, sd1 = SM1, TL = TL)
  } else {
    SaF <- piecewiseSa(t = Tn, sds = (2 / 3) * SMS, sd1 = (2 / 3) * SM1, TL = TL)
  }

  data.table::data.table(Tn = Tn, SaF = SaF)
}

# Piecewise ASCE 7-22 spectrum generator (SDS/SD1 form)
piecewiseSa <- function(t, sds, sd1, TL = 8) {
  Ts <- sd1 / sds
  T0 <- 0.2 * Ts
  t  <- pmax(t, 1e-6)
  ifelse(
    t < T0,  sds * (0.4 + 0.6 * t / T0),
    ifelse(
      t <= Ts, sds,
      ifelse(
        t <= TL, sd1 / t,
        sd1 * TL / t^2
      )
    )
  )
}
