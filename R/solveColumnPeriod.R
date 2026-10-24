#' Fundamental period of a layered soil profile as a Galerkin eigenvalue problem
#'
#' Solves the fundamental period of the soil column on rigid base, the site
#' period \eqn{T_r}, directly from the layers of one realization, and returns
#' it with its controls: the change under mesh doubling and the
#' density-weighted Rayleigh period beside it. [solveWedgePeriod()] solves
#' the fundamental period of the truncated shear beam with the same numerical method,
#' and [solveProfilePeriods()] applies both to a list of realizations.
#'
#' @section Mathematics:
#' Dobry, Oweis and Urzua (1976) model the profile as a shear beam on rigid
#' rock (hydrated copy, `eq-beam-wave-equation`,
#' `eq-beam-rock-condition` and `eq-beam-surface-condition`):
#' \deqn{\rho \frac{\partial^2 u}{\partial t^2} = \frac{\partial}{\partial z}
#'       \left(G \frac{\partial u}{\partial z}\right),\qquad u(H, t) = 0,\qquad
#'       \left(G \frac{\partial u}{\partial z}\right)_{z=0} = 0,}
#' with \eqn{z} the depth below the surface, \eqn{H} the depth of the rigid
#' base and \eqn{\rho = \gamma / g} (section Basic Equation). For harmonic motion
#' \eqn{u = X(z) e^{i \omega t}} (`eq-beam-harmonic-motion`) it is the regular
#' Sturm-Liouville problem \eqn{-(G X')' = \omega^2 \rho X} on \eqn{[0, H]},
#' whose smallest eigenvalue gives \eqn{T_r = 2\pi/\omega_1}. Multiplying by a
#' test function \eqn{v} with \eqn{v(H) = 0} and integrating by parts, the
#' surface term vanishes by the free condition and
#' \deqn{\int_0^H G X' v'\, dz = \omega^2 \int_0^H \rho X v\, dz.}
#' On layer \eqn{k}, of thickness \eqn{h_k}, \eqn{G_k = \rho_k V_k^2} and
#' \eqn{\rho_k = \gamma_k / g} are constant; the layers are stacked from the
#' surface in the order of their midpoint depths `zm`, and \eqn{H} is the sum
#' of their thicknesses.
#'
#' @section Numerical method:
#' The weak form is discretised with continuous piecewise-linear functions on
#' a mesh that follows the layers: every layer is divided into `elements`
#' equal elements, so every interface is a node, and the node at the bottom
#' of the last layer is fixed. With a width \eqn{w(z)}, 1 for the column and
#' the distance to the apex for the wedge of [solveWedgePeriod()], the
#' integrals of an element of length \eqn{L} between an upper node \eqn{a}
#' and a lower node \eqn{b} are exact, because \eqn{G} and \eqn{\rho} are
#' constant within a layer and \eqn{w} is linear: the element stiffness is
#' \deqn{k_e = G \frac{w_a + w_b}{2L},}
#' assembled as \eqn{k_e (\delta_a - \delta_b)(\delta_a - \delta_b)^T}, and
#' the consistent mass has diagonal entries \eqn{\rho L (3 w_a + w_b)/12} and
#' \eqn{\rho L (w_a + 3 w_b)/12} and coupling \eqn{\rho L (w_a + w_b)/12}.
#' The generalized eigenvalue problem \eqn{K \phi = \omega^2 M \phi} has
#' \eqn{K} and \eqn{M} symmetric positive definite for positive \eqn{G} and
#' \eqn{\rho} with the base fixed. The stiffness factors as
#' \eqn{K = U^T U} with \eqn{U = D^{1/2} C}, where \eqn{D} is the diagonal of
#' the \eqn{k_e} and \eqn{C} the bidiagonal difference matrix (1 on the
#' diagonal, \eqn{-1} above it: element \eqn{e} joins nodes \eqn{e} and
#' \eqn{e + 1}, the base node removed); \eqn{U} is the Cholesky factor of
#' \eqn{K}, and its inverse is explicit because \eqn{C^{-1}} is the upper
#' triangle of ones. The Cholesky reduction gives the symmetric matrix
#' \deqn{B = U^{-T} M U^{-1} = D^{-1/2} S D^{-1/2},\qquad
#'       S_{ij} = c_{\min(i, j)}\ (i \neq j),\qquad S_{ii} = c_i - M_{i, i+1},}
#' with \eqn{c_i} the consistent mass from the top to node \eqn{i} (the
#' cumulative row sums of \eqn{M}). The eigenvalues of \eqn{B} are
#' \eqn{1/\omega^2}, so the fundamental mode is its largest eigenvalue
#' \eqn{\mu_1}, computed by the direct symmetric solver
#' `eigen(symmetric = TRUE)`, and \eqn{T = 2\pi\sqrt{\mu_1}}; no iterative
#' solver or convergence loop is involved. Every entry of \eqn{S} is a sum
#' of positive masses, and the largest eigenvalue of a symmetric matrix is
#' resolved to relative machine precision.
#'
#' The discretisation is a Rayleigh-Ritz method: the numerical period never
#' exceeds the exact period of the layered beam, and it does not decrease
#' when every element is halved, since the trial spaces are nested. The
#' period returned is the one with `2 * elements` elements per layer, and the
#' mesh change \eqn{1 - T_{e}/T_{2e} \ge 0} between the declared and the
#' doubled mesh is its convergence certificate. Within a layer the solution
#' is smooth and the interfaces are nodes, so the error decreases with the
#' square of the element length and is asymptotically a third of the mesh
#' change. A profile whose stiffness tends to zero at the surface converges
#' more slowly. No threshold is applied to the mesh change: it is returned
#' for the consumer.
#'
#' @section Rayleigh bound:
#' `TrRayleighWeighted` is [estimateWeightedRayleighPeriod()] on the same
#' layers. It evaluates the Rayleigh quotient of Dobry, Oweis and Urzua
#' (`eq-rayleigh-frequency`),
#' \eqn{\omega^2 = \int \rho V^2 X'^2 dz / \int \rho X^2 dz}, on the
#' first-iteration deflection of `eq-rayleigh-mode-iteration`. Let \eqn{u} be
#' the piecewise-linear shape whose values \eqn{U_k} at the interfaces are
#' that deflection accumulated from the base with the slopes
#' \eqn{s_k = \sigma_k / (\gamma_k V_k^2)}, \eqn{\sigma_k} the overburden at
#' the midpoint of layer \eqn{k}, and \eqn{X_k = (U_{k-1} + U_k)/2}. Summation
#' by parts gives exactly \eqn{\int \rho V^2 u'^2 dz = g^{-1} \sum_k \gamma_k
#' h_k X_k}, while \eqn{\int \rho u^2 dz = g^{-1} \sum_k \gamma_k h_k
#' (U_{k-1}^2 + U_{k-1} U_k + U_k^2)/3 \ge g^{-1} \sum_k \gamma_k h_k X_k^2}.
#' The estimate's \eqn{\omega^2} is therefore at least the Rayleigh quotient
#' of \eqn{u}, which is at least \eqn{\omega_1^2} by Rayleigh's principle:
#' `TrRayleighWeighted` never exceeds the exact period of the layered
#' column, and `TrRayleighWeightedDifference` cannot exceed the
#' discretisation error of `Tr`, asymptotically `TrMeshChange / 3`. With a
#' uniform unit weight it is the
#' constant-density Method 7 of [estimateRayleighPeriod()].
#'
#' The package tests reproduce the uniform layer \eqn{T = 4H/V}
#' (`eq-uniform-period`), the power-law stratum of Dobry, Oweis and Urzua
#' (`eq-power-period`, with the height exponent \eqn{(2 - p)/2}), their
#' two-layer equation (`eq-twolayer-period-equation`) and the n-layer
#' frequency equation of Charara (2019, `eq-frequency-equation-condensed`).
#'
#' @param layers One layer table, such as an element of the `layers` column
#'   of [simulateSoilProfiles()]: a data.frame or data.table with one row per
#'   layer and numeric columns `hs` (thickness, m), `zm` (midpoint depth, m;
#'   distinct, it orders the layers), `gsat` (unit weight, kN/m3; any unit,
#'   only the ratios between layers matter) and `VSm` (m/s) or, when `VSm` is
#'   absent, `Gm` (MPa), converted with
#'   \eqn{V_S = \sqrt{9.81 \cdot 1000\, G_m / \gamma_{sat}}} as in
#'   [evaluateSoilLayers()]. Every value must be finite and positive; other
#'   columns are ignored.
#' @param elements One positive whole number of linear elements per layer of
#'   the declared mesh; the control doubles it.
#' @return A one-row data.table: `Tr` (s), the numerical site period on
#'   `2 * elements` elements per layer; `TrMeshChange`, the relative change
#'   \eqn{1 - T_e/T_{2e}}; `TrRayleighWeighted` (s), the density-weighted
#'   Rayleigh period; and `TrRayleighWeightedDifference`,
#'   `TrRayleighWeighted / Tr - 1`. Each value is
#'   a conditional evaluation of the supplied layers; the two controls are
#'   numerical diagnostics, not physical dispersion.
#' @export
#' @examples
#' Layers <- data.frame(hs = c(24, 8), zm = c(12, 28), gsat = c(17.66, 19.62),
#'                      VSm = sqrt(c(20e6 / 1800, 90e6 / 2000)))
#' solveColumnPeriod(Layers, elements = 8)
solveColumnPeriod <- function(layers, elements) {
  Layers <- .orderLayers(layers)
  Beam <- .solveShearBeam(Layers, elements, width = function(depth) rep(1, length(depth)))
  TrRayleighWeighted <- estimateWeightedRayleighPeriod(VSm = Layers$VSm, hs = Layers$hs, zm = Layers$zm,
                                                       gsat = Layers$gsat)
  data.table(Tr = Beam[["period"]], TrMeshChange = Beam[["meshChange"]], TrRayleighWeighted = TrRayleighWeighted,
             TrRayleighWeightedDifference = TrRayleighWeighted / Beam[["period"]] - 1)
}
