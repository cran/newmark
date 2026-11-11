# newmark 1.14.0

* The `BM19` pulse-motion `D50` equation now uses the intercept -0.148
  from the official erratum when `Ts < 0.1` s and `PGV > 150` cm/s,
  replacing 2.480. Other regimes and the ordinary-motion equations used
  by the project tables and hazard integration are unchanged.
* Licence changed from MIT to AGPL (>= 3). Files in `R/` carry SPDX
  licence headers. `LICENSE.note` has been removed.
* Column site spectra are reused at VS30 rounded to the nearest 1 m/s.
  `ShearDraws` retains its stored VS30 values. Mean-row site response uses
  the row's stored VS30 for PGATable, the common ky grid and reference metadata.
* `buildDnTables()` requires `draws`, the paired `Ts` and `VS30` of
  `ShearDraws`. The `NS` movements of each pair cycle through its columns
  in realization order. Each column supplies its site response, spectral
  sampling periods and displacement period. The movement count, model
  weights, table schemas and empirical exceedance calculation are unchanged.
* PGA, Sa13 and Sa15 are centered once over the complete cell sample on
  the column means weighted by their assigned movements. Column subsets
  are not centered separately. Identical columns reproduce the intensities
  and model deviates of a single-column sample with the same seed and NS.
  Public `sampleIntensities()` and `sampleSpectra()` calls are unchanged.
* The `DnTable` project step reads `ShearDraws.Rds`. Generate it explicitly
  with the `ShearDraws` step if absent. `ShearTable` mean rows continue to
  define the common ky grid, PGATable and reference metadata; `PGA.s` holds
  each movement's sampled PGA. The hazard-integration functions are unchanged.

# newmark 1.13.0

* `integrateDnHazard()` and `integrateDnHazardProfiles()` accept `BM19`
  among the `weights`, with the spectral acceleration at its degraded
  period `1.3 Ts` and its ordinary-motion equation, as the `DnTable`
  products evaluate it. The project route keeps `BM19` for shallow
  crustal regimes (`subduction` false), so the hazard integration can now
  take the model set of such a project. Every other model integrates as
  before.

# newmark 1.12.0

* New function `buildShearDraws()`: the realizations of the soil columns
  that `buildShearTable()` summarizes, 100 per geometry and material, each
  with the period of the truncated shear beam evaluated on its own `mo` and
  `VSo` with the rounding of the summary rows, together with the shear
  table summarized from the same draws. Realizations beyond 100 (`NS`) are
  drawn in a second pass after every pair has its first 100, so the first
  100 and the shear table do not depend on `NS`. As a rule of the product,
  a realization whose fitted `mo` is negative is returned with `mo` forced
  to 0 and with `an` and `Ts` evaluated at that 0, and a message counts
  them, while the shear table summarizes the fitted `mo`; a realization
  with `mo` of 2 or more, where the first-mode root is not defined, stops
  the function naming it. `buildShearTable()` keeps its arguments, its
  random numbers and its values.
* The new `ShearDraws` step of `process` follows `ShearTable` and runs only
  when selected. It requires `seed`, draws the realizations again from the
  contract and the seed, and publishes `ShearDraws.Rds` at the root of
  `path.tables` only when their summary equals the published
  `ShearTable.Rds` value for value; it never writes `ShearTable.Rds`. Its
  message reports how many realizations it publishes with a negative `mo`
  forced to 0. The
  contract key `NS.shear` (a whole number of at least 100, 100 when absent)
  sets the realizations per pair and is read by this step alone. Listing
  and query know the table.
* New function `integrateDnHazardProfiles()`: the annual displacement
  exceedance rate of `integrateDnHazard()` mixed over site profiles, each a
  pair of `Ts` and `VS30` with an optional weight, and over yield
  coefficients with weights, under a declared independence of profile and
  yield coefficient. The Poisson probabilities, the band and the rare-tail
  bound are formed on the mixed rates. It reports each component, the
  weighted fraction of components with `ky` outside a calibration range,
  and the Monte Carlo standard error of the mixture over the profiles
  (`mixtureError`, `bandMixtureError`), which is not an error bound. It is
  an operation of the displacement family, where the command takes `ky` as
  a table of `ky` and `w`. `integrateDnHazard()` does not change.
* The exports that the installer verifies and `newmark doctor` checks now
  include `integrateDnHazard()` and `listDisplacementFamilies()`, which the
  command already called but which that list had left out.
* The products `ShearTable`, `DnTable`, `kmaxTable` and `DnPOETable` and
  their columns do not change.

# newmark 1.11.0

* `buildDnPOETable()` counts the weighted empirical fraction of stored
  pointwise displacement realizations exceeding `DaH * Hs`, by site,
  demand, geometry, material, yield acceleration and threshold.
  `assembleDnPOETables()` combines the site tables. The fraction is
  conditional on the demand (return period or MCE) and inherits the
  positive lognormal branch represented by the flexible-model draws;
  stored zeros, including SR08 underflow, count as non-exceedance.
* The new `DnPOETable` step follows `kmaxTable` in `process` and can run
  alone from the accepted ShearTable and DnDraws, without PGATable or a
  new simulation. It publishes per-site and global DnPOETable products,
  available through listing and query. Both functions are CLI operations.
  This product does not integrate displacement hazard or add an
  `extrapolated` column.

# newmark 1.10.0

* New function `integrateDnHazard()`: the annual rate at which the Newmark
  displacement of a sliding mass with a fixed yield coefficient exceeds a
  threshold, integrated over one site's rock hazard curves with the site
  factor, the displacement regressions and their weights; it reports the
  Poisson probability over the exposure, the contribution of a band of
  return periods, a bound of the rare tail beyond the curves, the
  calibration range of the yield coefficient of each model and the
  numerical diagnostics of its quadratures. The command exposes it as the
  other operations of the displacement family.
* A `tables` argument given as a JSON array whose members are all file
  references is decoded like the object form; it was refused as if it were
  neither an array nor an object of tables.
* A file that an argument of the command refers to and that does not exist
  is a usage error, exit 2 with the `newmark:` line alone, as an unreadable
  referenced `.json` file already was; it exited 1 as a failed operation.
* `sampleSoilLayers()` counts the layers of a column as the ceiling of its
  height over the layer thickness with a stated rounding tolerance. The
  quotient of two decimal lengths can land a few units of double precision
  above an integer, as with 2.1 m divided by 0.3 m, and the column got one
  layer more than its height holds. Only a quotient that exceeds an integer
  by at most four times the double-precision epsilon, relative to its value,
  changes its count.
* The tests of the numerical periods and roots that compare two equal
  evaluations accept a difference in the last digit: under the Intel MKL
  library the two evaluations differ there, as CRAN's additional check of
  1.9.0 showed.
* The tests that reproduce the products of one hundred realizations byte
  for byte, the dense quadratures of `integrateDnHazard()` and the steep
  upper tail of `integrateSiteSpectrum()` run outside CRAN, in continuous
  integration and locally; on CRAN's Windows server the whole check took
  more than ten minutes with them.

# newmark 1.9.1

* The library signals its notes, progress and warnings as classed
  conditions of the AOM family (`moduleNote`, `moduleProgress` and
  `moduleWarning`) with the fields `operation`, `phase`, `unit`, `path`,
  `count` and `total`. Their texts lose the prefixes `dry run:`, `<step>:`,
  `evaluateSiteFactor:` and `>>`, and the kmax warnings name their table in
  parentheses. A failure after a publication of `processProject()` has the
  classes `newmarkPartialPublication`, `modulePartialPublication` and
  `moduleFailure`. The command reads its arguments from the launcher, whose
  renderer writes these conditions on stderr and takes `--verbose` and
  `--quiet`.
* The private predicates `.isName()` and `.isKmaxIdentifier()` become one
  predicate, `.isIdentifier()`, and `.dnDemandsOK()` and
  `.hasKmaxDemands()` become `.hasDemands()`. They are internal; no
  exported function changes its behavior.
* The test files are renamed to the `test-<function>.R` form; for
  example, `test-kmax-table.R` is now `test-buildKmaxTable.R`.
* The package help page, `?newmark`, describes all the families of
  exported functions.
* New dataset `VancouverSpectra`: the uniform hazard spectra of the 2020
  National Building Code of Canada for Vancouver on rock, from the
  seismic hazard tool of Natural Resources Canada. The examples of
  `convertAEP()`, `evaluateUniformHazard()` and `evaluateSiteSa()` read it.
* New function `listDisplacementFamilies()`: the families of displacement
  models (`ensemble`, `rigid` and `flexible`) that the displacement and
  kmax tables summarize as trees, each with the identifiers of its models.
  The command line runs it as the operation of the same name.
* `?newmark` loses the package name from its title. The help pages of
  `estimatePGV()` and `estimateAI()` cite Verri (2023), the
  source of the rock records their relations were fitted on. On CRAN, the
  tests that repeat a step of `processProject()` skip; one run of each
  step and the dry run still execute there.
* `Mw.gmdp` takes one magnitude for every demand, as before, or an object
  that names each return period of `TR.gmdp` and `MCE` with its own
  magnitude; `buildDnTables()` takes the same as a named vector in `Mw`.
  The DnTable step applies the magnitude of each demand. A contract with
  one magnitude gives the same tables as before.
* `newmark doctor` prints one JSON report of the installation, as every
  module of the family does: `schema`, `package`, `version`, `summary` and
  the checks `library`, `exports`, `packages`, `tools` and `receipt`, each
  with `name`, `required`, `status` and `detail`; it exits 1 when a
  required check fails, and `doctor --help` lists the checks.
  `describeModule(resources)` returns that report. `newmark --version`
  prints two lines, the package with its version and the build commit.
* A `.json` file that an argument of the command refers to and that cannot
  be read as JSON is a usage error, exit 2 with the `newmark:` line alone,
  as an invalid `INPUT.json` already was; it exited 1 as a failed
  operation.
* The accumulations that feed a discrete choice, a rounding or a cut in
  the site mixture, the site spectrum, the soil profile and the PGA and
  displacement tables sum in plain double on every platform, as the kmax,
  displacement and mixture summaries have since 1.9.0; the tables are
  unchanged on macOS, where the reference products were recorded.

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
