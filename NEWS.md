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
