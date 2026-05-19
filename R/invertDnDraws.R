###############################################################################
# invertDnDraws.R — Direct kmax inversion from Dn realisations
###############################################################################

#' Invert Dn draws to kmax(Da)
#'
#' For each realisation s, inverts each model's Dn(ky) curve to kmax(Da),
#' builds weighted ensemble mean, then averages over realisations.
#'
#' Log-log linear extrapolation is used when Da falls outside the support of
#' the computed Dn(ky) curve. No calibration-range clamping is applied — the
#' natural support of the curve bounds the inversion.
#'
#' All inputs must be in the same displacement units (cm).
#'
#' @param draws   data.table(ky, s, Dn, IDn) — Dn in cm
#' @param Da      numeric vector — displacement targets in cm
#' @param weights named numeric vector, ensemble weights by IDn
#' @param p       numeric quantiles in (0,1) to report alongside the mean.
#'                Default c(0.05, 0.10, 0.16, 0.84, 0.90, 0.95). The output
#'                always includes the ensemble mean as the first p row.
#'
#' @return data.table(Da, p, kmax). Column `p` is character: "mean" or the
#'         numeric quantile formatted as "0.84", "0.90", etc.
#' @examples
#' \dontrun{
#' uhs <- data.table::fread(
#'   system.file("extdata", "uhs.csv", package = "newmark"))
#' out <- fitDnCurve(uhs, ky = getDnKy(uhs, Ts = 0.3), Ts = 0.3, Mw = 7.5,
#'                   NS = 100,
#'                   weights = c(AM88=1, BT07=1, SR08=1, BM17=1))
#' invertDnDraws(out$draws, Da = c(2.5, 5, 10),
#'               weights = c(AM88=1, BT07=1, SR08=1, BM17=1))
#' }
#' @export
invertDnDraws <- function(draws, Da, weights,
                          p = c(0.05, 0.10, 0.16, 0.84, 0.90, 0.95)) {

  stopifnot(is.numeric(p), all(p > 0 & p < 1))

  w   <- as.numeric(weights)
  w   <- w / sum(w)
  ids <- names(weights)
  S   <- sort(unique(draws$s))

  OUT <- vector("list", length(S))

  for (j in seq_along(S)) {
    si <- S[j]
    km <- numeric(length(Da))

    for (k in seq_along(Da)) {
      x <- numeric(length(ids))
      for (i in seq_along(ids)) {
        DT <- draws[s == si & IDn == ids[i]]
        if (!nrow(DT)) {
          x[i] <- NA_real_
          next
        }
        data.table::setorder(DT, ky)
        # force Dn non-increasing in ky (monotone projection)
        DT[, Dn := rev(cummax(rev(Dn)))]
        # invert in log-log: Dn decreasing -> ky increasing
        # approx needs x increasing, so use Dn reversed
        lDn <- rev(log(DT$Dn + .Machine$double.eps))
        lky <- rev(log(DT$ky))
        lDa <- log(Da[k])
        n   <- length(lDn)
        if (lDa >= lDn[1L] && lDa <= lDn[n]) {
          x[i] <- exp(stats::approx(lDn, lky, xout = lDa,
                                    rule = 1L, ties = "ordered")$y)
        } else if (lDa > lDn[n]) {
          t <- (lky[n] - lky[n - 1L]) / (lDn[n] - lDn[n - 1L])
          x[i] <- exp(lky[n] + t * (lDa - lDn[n]))
        } else {
          t <- (lky[2L] - lky[1L]) / (lDn[2L] - lDn[1L])
          x[i] <- exp(lky[1L] + t * (lDa - lDn[1L]))
        }
      }
      km[k] <- weighted.mean(x, w, na.rm = TRUE)
    }

    OUT[[j]] <- data.table::data.table(Da = Da, kmax = km)
  }

  DT <- data.table::rbindlist(OUT, use.names = TRUE)
  DT[, list(
    p    = c("mean", sprintf("%.2f", p)),
    kmax = c(mean(kmax, na.rm = TRUE),
             stats::quantile(kmax, p, na.rm = TRUE, names = FALSE))
  ), by = .(Da)]
}
