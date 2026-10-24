#' Characteristic roots of the truncated inhomogeneous shear beam
#'
#' Roots \eqn{a_n} of the characteristic equation of Dakoulas and Gazetas
#' (1985), a shear beam whose modulus grows with depth as a power law of
#' exponent \eqn{m} and whose crest truncates the full wedge at the ratio
#' \eqn{\lambda}. Every root is solved in R at its node of the grid by the
#' exact method of [buildRootModel()] (`data-raw/CylinderRoots.R`, recorded
#' in `data-raw/README.md`). The grid, modes 1 to 8 with \eqn{m} from 0 to
#' 0.95 in steps of 0.005 and \eqn{\lambda} from 0 to 0.5 in steps of
#' 0.0025, is that of the tables of `newmark` 1.7.1, computed with
#' Mathematica. The surrogates of [buildRootModel()] are fitted to this
#' table.
#'
#' @format A data.table with 307128 rows and 4 variables, ordered by `n`,
#'   `m` and `l`:
#' \describe{
#'   \item{n}{Integer: mode number, 1 to 8.}
#'   \item{m}{Double: inhomogeneity exponent, 0 to 0.95.}
#'   \item{l}{Double: truncation ratio, 0 to 0.5.}
#'   \item{an}{Double: characteristic root of mode `n`.}
#' }
#' @source Dakoulas, P. and Gazetas, G. (1985). A class of inhomogeneous shear
#'   models for seismic response of dams and embankments. Soil Dynamics and
#'   Earthquake Engineering 4(4), 166-182.
"CylinderRoots"

#' Small-strain shear-modulus correlations by soil group
#'
#' Coefficients of the empirical correlations
#' \eqn{G_0 = A F(e) (\sigma'_0)^n} tabulated by Ishihara (1996), grouped by
#' soil family. Ishihara's Table 6.1 states \eqn{\sigma'_0} and \eqn{G_0} in
#' kPa; this table keeps the stress in kPa and expresses `A` for the
#' modulus in MPa (7.00 for the round-grained Ottawa sand of Hardin and
#' Richart (1963), printed 7000 for the modulus in kPa).
#' [evaluateSoilLayers()] applies each row at the octahedral mean effective
#' stress \eqn{p_m} (kPa) of a layer as
#' \deqn{G = A F(e)\, OCR^{m_1}\, 100^{n} (p_m/100)^{n} = A F(e)\, OCR^{m_1}\, p_m^{n}}
#' in MPa, with the overconsolidation factor \eqn{OCR^{m_1}} of Hardin and
#' Drnevich (1972) for fine soils (\eqn{m_1 = 0} otherwise). The factor
#' \eqn{100^{n} (p_m/100)^{n}} equals \eqn{p_m^{n}}; it does not normalize
#' `A`, which is not the modulus at 100 kPa. Carried
#' byte-identical from `newmark` 1.7.1 (`data-raw/MaterialProperties.xlsx`).
#'
#' @format A data.table with 20 rows and 10 variables:
#' \describe{
#'   \item{ModelID}{Character: correlation identifier; `SS75` uses the linear
#'     void-ratio function, every other row the quadratic one.}
#'   \item{GroupID}{Character: `Gravels`, `Sands` or `Fines`.}
#'   \item{NameID, AuthorID}{Character: material and authors.}
#'   \item{emin, emax}{Double: void-ratio interval of the correlation.}
#'   \item{A}{Double: coefficient of \eqn{F(e)\, p^n} for the modulus in MPa
#'     with the stress \eqn{p} in kPa.}
#'   \item{Ce}{Double: void-ratio coefficient of \eqn{F(e)}.}
#'   \item{n}{Double: stress exponent.}
#'   \item{UN}{Character: units of the modulus.}
#' }
#' @source Ishihara, K. (1996). Soil Behaviour in Earthquake Geotechnics.
#'   Oxford University Press, Table 6.1.
"ShearModelParameters"

#' ASCE 7-22 site classes
#'
#' Site classes of ASCE/SEI 7-22 Table 20.2-1 with their shear-wave velocity
#' intervals in m/s and ft/s. Carried byte-identical from `newmark` 1.7.1.
#'
#' @format A data.table with 9 rows and 4 variables: `SC`, `Description`,
#'   `Vs30 (m/s)` and `Vs30 (ft/s)`.
"SiteClass"
