#' Site spectral acceleration given the rock ordinate and PGA
#'
#' Natural-log mean and standard deviation of the spectral acceleration at a
#' site of velocity `VS30`, conditional on the rock ordinate `Sa` and on the
#' rock PGA of the same realization at the reference condition, under a site
#' factor of [evaluateSiteFactor()]. Applied to the mean ordinates of a
#' spectrum with the mean PGA, it is the site spectrum of a deterministic
#' rock spectrum.
#'
#' @section Mathematics:
#' With \eqn{\ln F \mid PGA \sim N(\mu_F, \sigma_F)} from
#' [evaluateSiteFactor()],
#' \deqn{\ln S_{aF} \mid S_a, PGA \sim N(\ln S_a + \mu_F, \sigma_F),}
#' row by row; its quantile at probability \eqn{p} is
#' \eqn{\exp(\ln S_a + \mu_F + \sigma_F \Phi^{-1}(p))} and its mean
#' \eqn{\exp(\ln S_a + \mu_F + \sigma_F^2/2)}. The dispersion is that of the
#' site factor alone: `Sa` and `PGA` are conditioning values, not
#' distributions, so the result does not carry the uncertainty of the rock
#' demand. When `VS30` equals `Vref` the result is \eqn{\ln S_a} with zero
#' dispersion.
#'
#' @param Sa Rock spectral accelerations in g at period `Tn`, finite and
#'   positive.
#' @param PGA Rock peak ground accelerations in g, aligned with `Sa`.
#' @param Tn One oscillator period in seconds.
#' @param VS30 One target velocity in m/s.
#' @param Vref One reference velocity in m/s, 760 or 3000.
#' @param model Site factor model of [evaluateSiteFactor()], `"ST20"` or
#'   `"ST17"`.
#' @return A data.table with `muLnSaF` and `sigmaLnSaF` (natural-log mean and
#'   standard deviation of the site ordinate, one row per `Sa`) and `ID`
#'   (the model), unrounded.
#' @export
#' @examples
#' evaluateSiteSa(Sa = c(0.1, 0.2, 0.3, 0.4), PGA = rep(0.1, 4), Tn = 0.5,
#'                VS30 = 300, Vref = 760, model = "ST17")
evaluateSiteSa <- function(Sa, PGA, Tn, VS30, Vref, model) {
  if (!is.numeric(Sa) || !length(Sa) || any(!is.finite(Sa)) || any(Sa <= 0)) {
    stop("Sa must be finite and positive")
  }
  if (length(PGA) != length(Sa)) stop("Sa and PGA must have equal length")
  Factor <- evaluateSiteFactor(PGA = PGA, Tn = Tn, VS30 = VS30, Vref = Vref, model = model)
  data.table(muLnSaF = log(Sa) + Factor$muLnF, sigmaLnSaF = Factor$sigmaLnF, ID = model)
}
