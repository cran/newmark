#' Paired intensity measures of a displacement scenario
#'
#' Draws `NS` joint realizations of the intensity measures that the six
#' displacement models read for a sliding mass of initial fundamental period `Ts`:
#' peak ground acceleration, the spectral accelerations at the degraded
#' periods \eqn{1.3 T_s} and \eqn{1.5 T_s}, and the peak ground velocity
#' and Arias intensity derived from each sampled PGA. One row is one
#' realization; every column of a row belongs to the same draw.
#'
#' @section Mathematics:
#' The periods \eqn{(0,\ 1.3 T_s,\ 1.5 T_s)} are sampled jointly with
#' [sampleSpectra()] under the correlation matrix
#' \eqn{C_{jk} = \rho(T_j, T_k)} of [correlateResiduals()], with
#' \eqn{C_{jj} = 1} and the zero period evaluated at 0.01 s. A quantile
#' spectrum that lacks \eqn{1.3 T_s} or \eqn{1.5 T_s} is first interpolated
#' there with [interpolateSpectrum()]; an analytic mixture is interpolated by
#' [sampleSpectra()] itself. Then \eqn{PGV} = [estimatePGV()] and
#' \eqn{I_a} = [estimateAI()] of each sampled PGA. The draws consume `NS`
#' trivariate normal deviates of the caller's random stream.
#'
#' @param spectrum The spectrum of the scenario as [sampleSpectra()] accepts
#'   it: a quantile table (`Tn`, `p`, `Sa` in g) or an analytic mixture of
#'   \eqn{\ln S_a} (`Tn`, `k`, `muLn`, `sigmaLn`, `w`).
#' @param Ts One finite positive initial fundamental period of the sliding mass in
#'   seconds.
#' @param NS Number of realizations, one positive integer.
#' @param zTrunc Truncation of the copula marginals in standard deviations
#'   for an analytic mixture, or `NA`, see [sampleSpectra()].
#' @return A data.table with one row per realization: `realization`, `PGA` (g),
#'   `Sa13` and `Sa15` (spectral accelerations at \eqn{1.3 T_s} and
#'   \eqn{1.5 T_s}, g), `PGV` (cm/s) and `AI` (m/s).
#' @export
#' @examples
#' Sigma <- c(0.2, 0.45, 0.75)
#' Spectrum <- data.frame(Tn = c(0, 0.4, 1), k = "oq", muLn = log(c(0.8, 0.55, 0.42)) - Sigma^2 / 2,
#'                        sigmaLn = Sigma, w = 1)
#' set.seed(46)
#' sampleIntensities(Spectrum, Ts = 0.4, NS = 5, zTrunc = NA)
sampleIntensities <- function(spectrum, Ts, NS, zTrunc) {
  if (!is.data.frame(spectrum)) stop("spectrum must be a table")
  if (!is.numeric(Ts) || length(Ts) != 1L || !is.finite(Ts) || Ts <= 0) {
    stop("Ts must be one finite positive period")
  }
  Periods <- c(0, 1.3 * Ts, 1.5 * Ts)
  DT <- as.data.table(spectrum)
  if (.classifySpectrum(DT) == "quantile") {
    DT <- DT[, .(Tn, p = as.character(p), Sa)]
    DT <- rbindlist(lapply(Periods, function(period) {
      if (period %in% DT$Tn) return(DT[Tn == period])
      interpolateSpectrum(DT, Tn = period)
    }), use.names = TRUE)
  }
  C <- outer(Periods, Periods, function(a, b) correlateResiduals(Tn = a, T0 = b))
  diag(C) <- 1
  OUT <- .drawSpectra(DT, Tn = Periods, correlation = C, NS = NS, zTrunc = zTrunc)
  data.table(realization = seq_len(nrow(OUT)), PGA = OUT[, 1L], Sa13 = OUT[, 2L], Sa15 = OUT[, 3L],
             PGV = estimatePGV(OUT[, 1L]), AI = estimateAI(OUT[, 1L]))
}
