###############################################################################
# getDnKy.R — Scenario-specific ky grid for Dn(ky) evaluation
###############################################################################

#' Build a log-spaced ky grid for one scenario
#'
#' Derives the effective intensities PGA, Sa(1.3Ts), Sa(1.5Ts) from the
#' scenario UHS and constructs a log-spaced ky grid covering the calibrated
#' range of the displacement models.
#'
#' Default bounds come from the calibration ranges of the models in the
#' ensemble:
#' - lower: 0.01g (BM17 lowest calibration point)
#' - upper: max(PGA, 0.80g). The physical limit is ky < PGA (no slip below
#'   the yield acceleration); the 0.80g floor extends the grid for
#'   low-PGA scenarios up to BM17's calibration ceiling.
#'
#' @param uhs  data.table with columns Tn, p, Sa (one scenario)
#' @param Ts   numeric scalar, fundamental period (s)
#' @param kyN  integer, number of grid points (default 30)
#' @param pRef character, probability level for intensity anchor (default "mean")
#'
#' @return numeric vector of ky values (log-spaced)
#' @examples
#' uhs <- data.table::data.table(Tn = 0, p = "mean", Sa = 0.5)
#' getDnKy(uhs, Ts = 0.3)
#' @export
getDnKy <- function(uhs, Ts, kyN = 30L, pRef = "mean") {

  uhs <- data.table::as.data.table(uhs)[, .(Tn, p, Sa)]
  uhs[, p := as.character(p)]

  PGA <- uhs[Tn == min(Tn) & p == pRef, Sa]

  kyMin <- 0.01
  kyMax <- max(PGA, 0.80)

  exp(seq(log(kyMin), log(kyMax), length.out = kyN))
}
