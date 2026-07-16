# newmark 1.3.0

## Model corrections (audit-driven; change numerical results)

* `designUHS()` now derives the spectral parameters per ASCE 7-22
  Sec. 21.4 (site-specific input): at MCER level,
  `SMS = 0.9 * max SaM(T)` over `T` in [0.2, 5] s and
  `SM1 = max(0.9 * max[T * SaM(T)], SaM(1 s))` over [1, 2] s when
  `vs30 > 442` m/s, else [1, 5] s (new optional `vs30` argument;
  `NULL` uses the conservative wide window). The previous seeding
  (`Ss = Sa(0.2 s)`, `S1 = Sa(1.0 s)`, `Fa = Fv = 1`) matched no
  ASCE 7-22 pathway for site-specific spectra. ASCE-table products
  shift on the order of 10 percent. The Sec. 21.2.3 floor and the
  Sec. 21.2.2 deterministic cap remain out of scope (documented).

## Fixes

* `designUHS()`: the output grid beyond `TL` was degenerate whenever
  `TL >= 5` (including the default 8 s), so the long-period branch
  (ASCE 7-22 Eq. 11.4-5) never appeared in the returned table; the
  grid now samples one octave beyond `TL`. Documentation now states
  precisely which parts of ASCE 7-22 the function implements (the
  Sec. 21.4 derivation and Sec. 11.4.5.2 two-period shape) and which
  it does not (the Sec. 21.2 floor and deterministic cap).

## Breaking changes

* The dead `Rrup` argument has been REMOVED from `rhoBJ()`, `fitSaF()`,
  `fitDnCurve()`, and `fitDnModel()`, and from the `fitDnScenario()`
  records contract (a legacy `Rrup` column in `records` is tolerated
  and ignored). The Baker & Jayaram (2008) inter-period correlation
  model has no distance dependence, so the argument never influenced
  any result: outputs are identical draw-for-draw under the same seed.
  Callers passing `Rrup` by position or by name must drop it.

# newmark 1.2.1

## Model corrections (audit-driven; change numerical results)

* `Dn_BT07()`: restored the published `+ 1.50 * Ts` term of
  Bray & Travasarou (2007) Eq. 5, previously missing (the function
  received `Ts` but never used it), and the `-0.22` intercept for
  nearly rigid sliding masses (`Ts < 0.05` s). The implementation now
  reproduces the worked example published by Bray (2007). Median
  displacements increase by `exp(1.5 * Ts)` (about x2.1 at
  `Ts = 0.5` s).

* `geSiteTable()` / `getSiteProperties()`: the Shibata-Soelarno (1975)
  row of the Ishihara (1996) shear-modulus table now uses its
  published LINEAR void-ratio form `F(e) = Ce - e/(1 + e)`; the
  quadratic form of the other models had been applied to it,
  under-predicting that model by a factor of about 65 and dragging the
  Sands-group mean modulus down (about -20 percent in `Go` for typical
  sand profiles).

* `fitModel.Ts()`: rebuilt on the classical Rayleigh quotient with the
  static-deflection trial shape. The previous trial shape violated the
  base boundary condition under package usage and overestimated the
  fundamental period of normal (stiffness-increasing) profiles by tens
  of percent (e.g. +44 percent on a 200/400 m/s two-layer column,
  +27 percent on the vignette profile). The rebuilt estimator is
  within -0.4 to -2.8 percent of the exact layered solution across
  test profiles and is invariant to the layer ordering of the input.
  All Ts-dependent products (site properties, flexible-block
  displacements) shift on regeneration.

* `Vs30toSID()`: the class-B lower bound is now 910 m/s (3000 ft/s,
  ASCE 7-22 Table 20.2-1), matching the packaged `SiteClass` table;
  it was 900. The class-E row in `SiteClass` now prints the
  correct inequality (`<150` m/s).

## Behaviour changes

* `Dn_BM19()` (and `fitDnCurve()`, `fitDnModel()`): near-fault pulse
  behaviour is now selected with the explicit `pulse` argument
  (default `FALSE` = ordinary motion), matching the papers' definition
  of pulse as a declared property of the record and the authors'
  spreadsheet. The automatic `PGV > 115` cm/s switch is retired; PGV
  remains a required predictor inside the pulse equations, whose
  published 150 cm/s sub-regime is unchanged.

# newmark 1.2.0

## New functions

* `F_ST20()` implements the final 2020 NGA-East ergodic site-amplification
  generation (Stewart et al. 2020 linear model and F760 impedance/gradient
  weighting; Hashash et al. 2020 nonlinear model with period-dependent
  reference velocity and PGA-dependent nonlinear sigma). Coefficients are
  regenerated from the USGS nshmp-lib machine-readable table
  (`data-raw/SiteAmp2020.R`); the implementation is validated at machine
  precision against the OpenQuake `nga_east.py` reference functions.

* `fitSaFMixture()` computes site-amplified spectra from an analytic
  scenario mixture (one lognormal component per rupture-variant x GMPE
  branch, as produced by the DSHA scenario engine). The composition of the
  rock mixture with the lognormal site dispersion is exact: Gauss-Hermite
  quadrature over the PGA channel plus monotone mixture-CDF inversion. No
  random numbers are drawn - outputs are byte-stable; `Sa` and `SaF` are
  quantiles of the same conditional demand distribution and
  `AF = SaF / Sa` is the amplification of the p-quantile.

## Behaviour changes

* `fitSaF()` gains a working `models` menu: `"ST20"` (new default) and
  `"ST17"`. Pass `models = "ST17"` to reproduce the 2017-generation
  amplification.

* `F_ST17()` sigma corrected: the nonlinear component now applies the
  PGA-dependent multiplier of its source model
  (`sigma_Fnl = sigma_f2 * ln((PGA_3000 + f3)/f3)`, PEER 2017/05 Eq. 2.5
  as implemented by USGS nshmp-lib and OpenQuake). Earlier versions used
  bare `sigma_f2`, overstating dispersion at low PGA and understating it
  at high PGA. `muLnF` is unchanged; `sdLnF` now varies with PGA.

## Documentation

* Corrected provenance: the `F_ST17()` family implements the PEER Report
  2017/05 NGA-East recommendations (not Seyhan & Stewart 2014 as
  previously documented). The inter-period correlation model is
  Baker & Jayaram (2008); it has no distance dependence and the `Rrup`
  arguments are retained for backward compatibility only.

# newmark 1.1.1

## New functions

* `fitDnScenario()` fits the Newmark displacement curve and Monte Carlo draws
  for a single scenario row (one `(TR, IDg)` combination). Returns a list with
  `curve` and `draws` data.tables matching the contract of `fitDnCurve()`.

* `fitDnScenarios()` applies `fitDnScenario()` over a `records` data.table of
  scenario rows and row-binds the results. Accepts `useParallel` and `workers`
  arguments for fork-based parallel execution on multi-core systems.

* `fitKmaxScenario()` fits the seismic coefficient kmax for a single scenario
  row given pre-computed displacement draws and model weights. Returns a
  data.table with `Da`, `p`, `kmax`, `PGA`, `Kh`, `TR`, `IDg`, and `IDn`.

* `fitKmaxScenarios()` applies `fitKmaxScenario()` over a `records`
  data.table and row-binds the results. Supports the same `useParallel` /
  `workers` interface as `fitDnScenarios()`.

## Performance

* `fitDnScenarios()` now runs Dn scenario rows in parallel when
  `useParallel = TRUE` and `workers > 1L`. On systems where fork-based
  parallelism is available (Linux, macOS), wall time scales nearly linearly
  with `workers` for large scenario grids. Sequential execution (the default)
  is unchanged.

## Internal tooling

* `inst/cran-release.R` simplified to a two-phase `source()` from the
  garage canonical scripts (`cran-check.R` then `cran-submit.R`). The
  prior inline copy of those steps is removed; the garage copy is the
  single source of truth.

# newmark 1.1.0

## Bug fixes

* `inst/CITATION`: replace `utils::packageVersion("newmark")` with
  `meta$Version` so the citation file evaluates correctly when the
  package is not yet installed (CRAN pre-test failure).
* `inst/WORDLIST`: add 135 domain-specific terms (author surnames,
  technical vocabulary, package-internal identifiers) so
  `spelling::spell_check_package()` exits clean. Add `Language: en-US`
  to DESCRIPTION.



## Documentation

* Vignettes rebuilt for fidelity to the implemented code and the
  underlying paper (Verri Kozlowski 2026, working). Hallucinated
  cross-references to deleted articles, `oqt`-pipeline terminology
  leaking into `newmark` prose, and inconsistent function-name
  references were removed.

* New vignette `dynamic-site-response` documents the dsra-side
  methodology merged into `newmark` 1.1.0: Ishihara (1996)
  small-strain shear-modulus model, Gazetas & Dakoulas (1985)
  inhomogeneous truncated shear-beam fundamental period,
  synthetic-profile generation via `getSiteProperties()` /
  `geSiteTable()`, and the ASCE 7-22 site-class utilities
  `Vs30toSID()` / `SIDtoVs30()`.

* `pipeline` vignette rewritten as a function-level overview of
  the four library workflows (dynamic site response, hazard import,
  displacement curves, seismic coefficient), with a single
  top-of-file note that orchestration on a hazard grid (looping
  scenarios, persisting tables) is the consumer's responsibility.

* `ensemble-formulation` vignette now states the inverse-problem
  formulation explicitly, documents the BM19 near-fault PGV >
  115 cm s⁻¹ correction and the `NFC` D100/D50 selector, and lists
  the operational defaults (Monte Carlo size, k_y grid, reported
  quantiles) that the package and downstream consumers share.

* `newmark-quickstart` vignette: stale links to deleted articles
  removed; `Next steps` rewritten as in-package `vignette()` and
  `?fnname` references.

## Improvements (no behavioural change)

* `interpolateSaTable()` and `buildQSpline()` now emit a descriptive
  `warning()` when called on a UHS table with duplicate `(Tn, p)` rows,
  pointing the user to the standard disaggregation key
  `(siteID, ID, Vs30, TR)`. R's cryptic `regularize.values: collapsing
  to unique 'x' values` warning is suppressed in favour of the new
  message — it is replaced, not stacked on top. Numerical behaviour
  is unchanged: duplicate rows are still silently collapsed via
  `stats::approx` / `stats::splinefun` `ties = "mean"`. SoT-validated
  bit-identical on clean fixtures (5/5 direct + 36/36 pipeline) and
  byte-equal output on duplicated input.

## Internal renames (no public API change)

* Private helpers `.summarizeDn()`, `.summarizeDnDraws()`,
  `.aggregateWeighted()` parameter `p.vals` renamed to `pVals`.
  Producer `.sampleDnScenario()` returns the list element under the
  new key `pVals` (was `p.vals`). Consumers `fitDnCurve()` and
  `fitDnModel()` updated to read `.sample$pVals`. All references are
  intra-package; no exported symbol changes signature. SoT-validated
  (Stage B helpers direct identity 3/3; Stage C swap-in-place via
  `assignInNamespace` over the 6 affected functions running the full
  fitDnCurve grid 216/216, both bit-identical).

## Breaking changes

* `getDnKy()` parameters renamed for STYLE.md compliance: `ky.n` →
  `kyN`, `p.ref` → `pRef`. Output is bit-identical for any call that
  uses the defaults or passes the args by position. Callers that
  passed by name must update the keyword. No live caller in the
  ecosystem (`oqt`, research scripts, MCP server, skill templates)
  was found to pass these by name; the rename is a clean keyword
  change. SoT-validated (Stage A direct identity 27/27 across 3
  fixtures × 3 Ts × 3 kyN values; Stage C swap-in-place via
  `assignInNamespace` 216/216 across the full fitDnCurve pipeline).

* `fitDnCurve()` no longer accepts `NS.draws`. The parameter was a
  partial-retention knob (filter `draws` to the first `NS.draws`
  realisations) that no caller used; every observed call site set
  `NS.draws = NS`. The default `NS.draws = 0L` produced
  `draws = NULL`, which silently broke `invertDnDraws()` for callers
  who followed the README pattern. `draws` is now always returned
  as a populated `data.table` (or an empty `data.table` with the
  four canonical columns when no model produces output, replacing
  the prior `NULL` corner-case return).

  Migration: drop `NS.draws = ...` from any `fitDnCurve()` call.
  Output for `NS.draws = NS` callers (the universal pattern) is
  bit-identical (verified across vignette + 2 real UHSTable fixtures
  × 3 seeds × 2 NS × 3 Ts × 2 NFC × 2 weight sets).

# newmark 1.0.0

* First CRAN release.
