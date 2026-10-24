#' ASCE 7-22 site class of a shear-wave velocity
#'
#' Classifies \eqn{V_{S30}} by the intervals of ASCE/SEI 7-22 Table 20.2-1:
#' A at or above 1500 m/s, B from 910, BC from 640, C from 440, CD from 300,
#' D from 210, DE from 150 and E below 150 m/s. The lower bound of each
#' interval belongs to it. Class F requires site-specific criteria and is not
#' assigned from a velocity.
#'
#' @param VS30 Finite non-negative velocities in m/s.
#' @return A character vector of class identifiers aligned with `VS30`.
#' @export
#' @examples
#' classifySite(c(120, 150, 790, 910, 3000))
classifySite <- function(VS30) {
  if (!is.numeric(VS30) || any(!is.finite(VS30)) || any(VS30 < 0)) {
    stop("VS30 must be finite non-negative velocities")
  }
  c("E", "DE", "D", "CD", "C", "BC", "B", "A")[findInterval(VS30, c(0, 150, 210, 300, 440, 640, 910, 1500))]
}
