#' Uncertainty analysis in dynamic site and slope response
#'
#' Estimates the seismic response of sites and slopes from the results of a
#' seismic hazard analysis. Uncertainty in the soil profile and in the ground
#' motion is represented by Monte Carlo samples, and every result is reported
#' with its representative value, its applicable quantiles or dispersion, and
#' the population or conditioning to which it refers. The exported functions
#' form the families below.
#'
#' Soil profiles and their periods: [buildGeometry()], [sampleSoilLayers()],
#' [evaluateSoilLayers()], [characterizeProfile()], [estimateVS30()],
#' [classifySite()], [estimateRayleighPeriod()],
#' [estimateWeightedRayleighPeriod()], [buildRootModel()], [evaluateRoot()],
#' [evaluateBeamPeriod()], [solveColumnPeriod()], [solveWedgePeriod()],
#' [solveProfilePeriods()], [simulateSoilProfiles()], [propagateBeamPeriod()],
#' [propagateSitePeriod()], [summarizeProfiles()], [buildShearTable()].
#'
#' Hazard inputs: [readHazardWorkbook()], [convertAEP()],
#' [evaluateUniformHazard()].
#'
#' Site response: [evaluateSiteFactor()], [evaluateSiteSa()],
#' [integrateSiteMixture()], [summarizeMixture()], [refineSiteMixture()],
#' [integrateSiteSpectrum()], [refineSiteSpectrum()], [simulateSiteSpectrum()],
#' [buildPGATable()].
#'
#' Earthquake-induced displacement: [correlateResiduals()],
#' [buildQuantileSpline()], [interpolateSpectrum()], [sampleSpectra()],
#' [sampleIntensities()], [estimatePGV()], [estimateAI()], [listDisplacementFamilies()],
#' [selectDisplacementModels()], [buildKyGrid()], [buildKyLimits()],
#' [evaluateRigidDisplacement()], [evaluateFlexibleDisplacement()],
#' [evaluateDisplacements()], [summarizeDisplacements()], [summarizeEnsemble()],
#' [integrateDnHazard()], [buildDnTables()], [assembleDnTables()].
#'
#' Pseudostatic coefficient: [solveKmax()], [summarizeKmax()],
#' [buildKmaxTable()], [assembleKmaxTables()].
#'
#' Projects: [processProject()], [listProjectTables()], [queryProjectTable()],
#' [describeModule()].
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
