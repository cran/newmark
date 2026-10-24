# newmark 1.9.0

This release replaces the interface of newmark 1.3.0, the previous CRAN
release. The package was rebuilt around functions named by what they
compute, one family at a time: soil profiles, site response, hazard
inputs, displacements and the pseudostatic coefficient.

## Breaking changes

* The 34 exports of 1.3.0 are gone. The two that had survived as
  compatibility adapters, `buildGMDP()` and `fitSaFMixture()`, are
  retired in this release: `fitSaFMixture()` is replaced by
  `integrateSiteMixture()`, and the hazard inputs that `buildGMDP()`
  staged are read with `readHazardWorkbook()`, `convertAEP()` and
  `evaluateUniformHazard()`.
* Removed exports of 1.3.0 with no function of the same name:
  `approx.spline`, `buildQSpline`, `checkUHS`, `designUHS`, `Dn_AM88`,
  `Dn_BM17`, `Dn_BM19`, `Dn_BT07`, `Dn_JB07`, `Dn_SR08`, `F_ST17`,
  `F_ST20`, `fitDnCurve`, `fitDnModel`, `fitDnScenario`,
  `fitDnScenarios`, `fitKmaxScenario`, `fitKmaxScenarios`,
  `fitModel.Ts`, `fitSaF`, `geSiteTable`, `getCylinderRoots`, `getDnKy`,
  `getKyLimits`, `getSiteProperties`, `interpolateSaTable`,
  `invertDnDraws`, `rhoBJ`, `SaF_ST17`, `sampleSaCorr`, `SIDtoVs30`,
  `Vs30toSID`.
* The site amplification models ST17 and ST20 are selected through the
  `model` argument of `evaluateSiteFactor()` and the site response
  functions. The displacement models AM88, BM17, BM19, BT07, JB07 and
  SR08 are selected through `selectDisplacementModels()` and evaluated
  by `evaluateRigidDisplacement()`, `evaluateFlexibleDisplacement()`
  and `solveKmax()`.
* The dataset `SiteTable` is no longer shipped. `CylinderRoots`,
  `ShearModelParameters` and `SiteClass` remain.
* `Depends` is `R (>= 4.2)`. `digest`, `Hmisc`, `randomForest`, `rpart`,
  `stringr` and `triangle` are no longer imported; `readxl` moves to
  Suggests; `jsonlite` and `parallel` are imported.

## New interface

* Soil profiles: `simulateSoilProfiles()`, `sampleSoilLayers()`,
  `evaluateSoilLayers()`, `characterizeProfile()`, `summarizeProfiles()`,
  `estimateVS30()`, `classifySite()`, `estimateRayleighPeriod()`,
  `estimateWeightedRayleighPeriod()`, `evaluateBeamPeriod()`,
  `propagateBeamPeriod()`, `propagateSitePeriod()`, `solveColumnPeriod()`,
  `solveProfilePeriods()`, `solveWedgePeriod()`, `buildRootModel()`,
  `evaluateRoot()`, `buildGeometry()`, `buildShearTable()`.
* Site response: `evaluateSiteFactor()`, `evaluateSiteSa()`,
  `integrateSiteMixture()`, `integrateSiteSpectrum()`,
  `refineSiteMixture()`, `refineSiteSpectrum()`, `simulateSiteSpectrum()`,
  `summarizeMixture()`, `interpolateSpectrum()`, `buildQuantileSpline()`,
  `correlateResiduals()`, `sampleSpectra()`.
* Hazard inputs: `readHazardWorkbook()`, `convertAEP()`,
  `evaluateUniformHazard()`, `buildPGATable()`, `sampleIntensities()`,
  `estimateAI()`, `estimatePGV()`.
* Displacements and pseudostatic coefficient: `selectDisplacementModels()`,
  `evaluateRigidDisplacement()`, `evaluateFlexibleDisplacement()`,
  `evaluateDisplacements()`, `summarizeDisplacements()`,
  `summarizeEnsemble()`, `buildDnTables()`, `assembleDnTables()`,
  `buildKyGrid()`, `buildKyLimits()`, `solveKmax()`, `summarizeKmax()`,
  `buildKmaxTable()`, `assembleKmaxTables()`.
* Projects: `processProject()`, `listProjectTables()`,
  `queryProjectTable()`, `describeModule()`.
