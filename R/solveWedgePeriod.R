#' Fundamental period of the truncated shear beam as a Galerkin eigenvalue problem
#'
#' Solves the fundamental period \eqn{T_s} of the truncated inhomogeneous shear beam
#' directly from the layers of one realization,
#' and returns it with its controls: the change under mesh doubling and the
#' period of the analytic route beside it, with their relative difference.
#' [solveColumnPeriod()] documents the shared numerical method, and
#' [solveProfilePeriods()] applies both to a list of realizations.
#'
#' @section Mathematics:
#' Dakoulas and Gazetas (1985) measure the depth \eqn{z} from the apex where
#' the extended slopes meet; the crest lies at \eqn{z = h}, the base at
#' \eqn{z = H} and the truncation ratio is \eqn{\lambda = h/H} (hydrated
#' copy, Fig. 1 on printed page 167 and `eq-free-truncation-ratio`). The width of
#' the dam is proportional to \eqn{z}, and the free lateral vibration obeys
#' (`eq-free-motion`)
#' \deqn{\rho \frac{\partial^2 u}{\partial t^2} = \frac{1}{z}
#'       \frac{\partial}{\partial z}\left[G(z)\, z\, \frac{\partial u}{\partial z}\right],}
#' with zero displacement at the base and zero shear at the crest
#' (`eq-free-base-condition` and `eq-free-crest-condition`):
#' \eqn{u(H) = 0} and \eqn{G z\, \partial u/\partial z = 0} at \eqn{z = h}.
#' For the slope of [buildGeometry()], of height \eqn{H_s} from crest to base
#' and truncation ratio \eqn{\lambda_o}, the extended slopes meet
#' \eqn{h = \lambda_o H_s / (1 - \lambda_o)} above the crest and
#' \eqn{H = H_s / (1 - \lambda_o)} (see [evaluateBeamPeriod()]). For harmonic
#' motion \eqn{u = U(z) e^{i \omega t}} (`eq-free-harmonic`), multiplying by
#' \eqn{z v} with \eqn{v(H) = 0} and integrating by parts, the crest term
#' vanishes by the free condition and
#' \deqn{\int_h^H G\, z\, U' v'\, dz = \omega^2 \int_h^H \rho\, z\, U v\, dz,}
#' a regular Sturm-Liouville problem with weight \eqn{z} whose smallest
#' eigenvalue gives \eqn{T_s = 2\pi/\omega_1}. The layers of the realization
#' are stacked from the crest, \eqn{z = h} plus the depth below it, in the
#' order of their midpoint depths, with \eqn{G_k = \rho_k V_k^2} and
#' \eqn{\rho_k = \gamma_k / g} constant in each; their thicknesses must add up
#' to \eqn{H_s}.
#'
#' With \eqn{\lambda_o = 0} the crest is the apex, where the width and the
#' crest term vanish with \eqn{z} itself. The problem stays well posed in the
#' energy space: the stiffness and mass of the first element remain positive,
#' and the branch \eqn{Y_q} of the general solution (`eq-free-solution`) has
#' infinite energy \eqn{\int G z U'^2 dz} at \eqn{z = 0}, so the Galerkin
#' solution converges to the regular branch, the paper's untruncated case.
#' The mesh starts at the apex, with no offset; the package tests reproduce
#' the untruncated periods of Table 7 of the paper (printed page 172).
#'
#' @section Analytic route:
#' `TsAnalytic` is the closed form of the paper (`eq-natural-frequency-period`)
#' evaluated by [evaluateBeamPeriod()] with the apex-to-base
#' height \eqn{H_s/(1 - \lambda_o)}, the basal velocity `VSo` and the first
#' root of the characteristic relation at `mo` and `lo` from
#' [evaluateRoot()] with the solution method `method`. It is exact for the
#' power law \eqn{G = G_b (z/H)^m} measured from the apex; with the exponent
#' that [characterizeProfile()] fits from the crest it is an approximation.
#' `TsAnalyticDifference` (`TsAnalytic / Ts - 1`) is returned so that a
#' consumer can flag a divergence under a threshold of its own; no
#' threshold is applied here. The package tests reduce the numerical period
#' to the closed form on layered power-law profiles with \eqn{m} = 0, 0.5,
#' 0.95 and 1.5 and \eqn{\lambda} = 0.075, 0.2, 0.351 and 0.5, and reproduce
#' the shear-beam periods of the five dams of Table 8 of the paper (printed
#' page 175) within 0.2 percent.
#'
#' @inheritSection solveColumnPeriod Numerical method
#'
#' @param layers One layer table, as for [solveColumnPeriod()], whose
#'   thicknesses add up to `Hs`.
#' @param Hs One positive slope height from crest to base, in metres.
#' @param lo One truncation ratio \eqn{\lambda_o} in \eqn{[0, 1)}; see
#'   [buildGeometry()].
#' @param mo One inhomogeneity exponent of the analytic route, the `mo` that
#'   [characterizeProfile()] fitted on the same layers.
#' @param VSo One basal shear-wave velocity of the analytic route in m/s, the
#'   `VSo` of the same fit.
#' @param method The solution method of the root model of the analytic route:
#'   `"exact"`, or the surrogate `"lm"` or `"nlm"`; see [buildRootModel()].
#' @param elements One positive whole number of linear elements per layer of
#'   the declared mesh; the control doubles it.
#' @return A one-row data.table: `Ts` (s), the numerical fundamental period on
#'   `2 * elements` elements per layer; `TsMeshChange`, the relative change
#'   \eqn{1 - T_e/T_{2e}}; `TsAnalytic` (s), the period of the analytic route;
#'   `TsAnalyticDifference`, `TsAnalytic / Ts - 1`; and `an`, `anExact` and
#'   `anWithinSupport`, the root of the analytic route and its diagnostics
#'   from [evaluateRoot()]. Each value is a conditional evaluation of the
#'   supplied layers and paired inputs; the controls are numerical
#'   diagnostics, not physical dispersion.
#' @export
#' @examples
#' Layers <- data.frame(hs = rep(1, 10), zm = seq(0.5, 9.5, by = 1), gsat = 19, VSm = 200)
#' solveWedgePeriod(Layers, Hs = 10, lo = 0.35, mo = 0, VSo = 200, method = "exact", elements = 4)
solveWedgePeriod <- function(layers, Hs, lo, mo, VSo, method, elements) {
  if (!is.numeric(Hs) || length(Hs) != 1L || !is.finite(Hs) || Hs <= 0) {
    stop("Hs must be one finite positive height in metres")
  }
  if (!is.numeric(lo) || length(lo) != 1L || !is.finite(lo) || lo < 0 || lo >= 1) {
    stop("lo must be one truncation ratio in [0, 1)")
  }
  if (!is.numeric(mo) || length(mo) != 1L || !is.numeric(VSo) || length(VSo) != 1L) {
    stop("mo and VSo must be one number each")
  }
  Layers <- .orderLayers(layers)
  if (!isTRUE(all.equal(sum(Layers$hs), Hs))) {
    stop("the layer thicknesses add up to ", sum(Layers$hs), " m, not to Hs = ", Hs, " m")
  }
  # The analytic route validates mo and VSo before the numerical solve.
  Root <- evaluateRoot(buildRootModel(mode = 1, method = method), mo = mo, lo = lo)
  TsAnalytic <- evaluateBeamPeriod(Hs = Hs, lo = lo, mo = mo, VSo = VSo, an = Root$an)
  Apex <- lo * Hs / (1 - lo)
  Beam <- .solveShearBeam(Layers, elements, width = function(depth) Apex + depth)
  data.table(Ts = Beam[["period"]], TsMeshChange = Beam[["meshChange"]], TsAnalytic = TsAnalytic,
             TsAnalyticDifference = TsAnalytic / Beam[["period"]] - 1,
             an = Root$an, anExact = Root$exact, anWithinSupport = Root$withinSupport)
}
