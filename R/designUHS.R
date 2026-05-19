#' MCER / Design elastic spectrum for a Vs30/Vref slice (ASCE 7-22)
#'
#' The input `UHS` slice must contain at least `Tn` and `SaF`.
#' Typical usage passes a single TR and `p == "mean"`.
#'
#' @param UHS data.table slice from UHSTable (single TR, p == "mean").
#' @param TL  numeric, long-period transition (seconds). Default 8.
#' @param spectrum character, "MCER" (default) or "DESIGN" (2/3 * MCER).
#'
#' @return data.table with columns `Tn, SaF`.
#' @examples
#' uhs <- data.table::data.table(Tn = c(0, 0.2, 0.5, 1, 2), SaF = c(0.4, 1.0, 0.9, 0.6, 0.3))
#' designUHS(uhs)
#' @export
designUHS <- function(UHS, TL = 8, spectrum = c("MCER", "DESIGN")) {
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

  # Rock ordinates (already site-adjusted if SaF came from sdLnSaFC1)
  Ss <- fSa(0.2)  # S_s
  S1 <- fSa(1.0)  # S_1
  if (!is.finite(Ss) || !is.finite(S1)) {
    stop("designUHS(): invalid Ss or S1 from interpolation")
  }

  # Site factors (set to 1 here; apply Fa/Fv upstream if desired)
  Fa  <- 1
  Fv  <- 1
  SMS <- Fa * Ss
  SM1 <- Fv * S1

  # Corner periods and evaluation grid
  Ts <- SM1 / SMS
  if (!is.finite(Ts) || Ts <= 0) stop("designUHS(): Ts must be positive")
  if (!is.finite(TL) || TL <= 0) stop("designUHS(): TL must be positive")

  T0 <- 0.2 * Ts
  TnDecay1 <- exp(seq(log(max(Ts, 1e-6)), log(TL),          length.out = 20))
  TnDecay2 <- exp(seq(log(TL),               log(max(TL, 5)), length.out = 20))
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
