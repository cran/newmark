#' newmark: uncertainty analysis in dynamic site and slope response
#'
#' The library is built around explicit scientific calculations
#' whose results carry their representative value, their applicable quantiles
#' or dispersion and the information needed to propagate uncertainty. The
#' soil and dynamics family separates synthetic layer generation
#' ([sampleSoilLayers()]) from the deterministic layer physics
#' ([evaluateSoilLayers()]), the profile characterization
#' ([characterizeProfile()], [estimateVS30()], [estimateRayleighPeriod()],
#' [classifySite()]), the characteristic roots of the truncated shear beam
#' ([buildRootModel()], [evaluateRoot()]), the analytic modal periods of the
#' truncated shear beam ([evaluateBeamPeriod()]) and the numerical fundamental
#' periods of a layered soil profile and a truncated shear beam ([solveColumnPeriod()],
#' [solveWedgePeriod()], [solveProfilePeriods()]). [simulateSoilProfiles()] retains every
#' realization with its paired properties and layers; [propagateSitePeriod()]
#' and [propagateBeamPeriod()] solve the fundamental periods `Tr` and `Ts` of each
#' realization before any summary, and [summarizeProfiles()]
#' reduces them to marginal statistics. [buildShearTable()] reproduces the
#' product of `newmark` 1.7.1 for comparison with the installed oracle.
#'
#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @importFrom stats approx coef glm lm predict quantile runif sigma uniroot
#' @importFrom utils head
NULL

# Names resolved by data.table inside package code.
utils::globalVariables(c(
  "USCS", "GroupID", "an", "l", "m", "n", "Gm", "VSm", "zm", "Go", "mo", "Ts",
  "VSo", "VS30", "Z500", "Z1000", "level", "IDg", "IDm", "Hs", "b", "s", "lo",
  "X", "Y", "Z", ".", "Dr", "TrRayleigh", "emax", "emin", "eo", "group", "gsat",
  "ModelID", "Ce", "A", "Fe", "G", "layer", "OCR", "pm", "m1", "i.query", "x.an"
))

# Lazy-loaded datasets of this package, resolved whether the package is
# attached, only loaded, or loaded from source with pkgload.
.dataset <- function(name) {
  getExportedValue("newmark", name)
}
