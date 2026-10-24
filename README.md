# newmark

`newmark` is an R package for uncertainty analysis in dynamic site and slope
response. The fundamental period of a layered soil column, which describes its
free, undamped horizontal shear vibration on rigid rock, is approximated by the
first Rayleigh iteration in layered form and by a finite element eigenvalue
problem. Uniform hazard spectra are evaluated from seismic hazard curves, which
give the mean annual frequency with which a ground-motion intensity level is
exceeded at a site. The site-response model transforms spectral acceleration at
a reference condition into spectral acceleration at a target site condition
described by `VS30`, through a random amplification factor that is lognormal at
fixed period, velocities and reference peak ground acceleration. Permanent
seismic displacement measures the accumulated movement of a sliding mass during
an earthquake. For fixed model inputs, rigid-block and flexible-block
regressions represent the positive displacement by a lognormal distribution,
and at a fixed yield coefficient an ensemble pools the realizations of its
models with normalized weights. The yield coefficient `kmax` for an admissible
displacement `Da` is the infimum of the positive yield coefficients at which
the probability that the projected, non-increasing displacement curve exceeds
`Da` is at most a target exceedance probability. A project stores
soil-property summaries, seismic displacement realizations and their
summaries, and yield coefficients associated with admissible displacements in
the tables `ShearTable`, `DnDraws`, `DnTable` and `kmaxTable`.

| Family | Article |
| --- | --- |
| `soil` | [Fundamental period of a layered soil column: Rayleigh approximations and numerical solution](https://kanameishi.github.io/newmark/articles/soil.html) |
| `hazard` | [Uniform hazard spectra from seismic hazard curves](https://kanameishi.github.io/newmark/articles/hazard.html) |
| `siteResponse` | [Conditional spectral acceleration under site amplification](https://kanameishi.github.io/newmark/articles/site-response.html) |
| `displacement` | [Distribution of permanent seismic displacement in a model ensemble](https://kanameishi.github.io/newmark/articles/displacement.html) |
| `kmax` | [Yield coefficient for an admissible seismic displacement](https://kanameishi.github.io/newmark/articles/kmax.html) |
| `project` | [Generation and querying of project tables](https://kanameishi.github.io/newmark/articles/project.html) |
| `meta` | [Installed module identity](https://kanameishi.github.io/newmark/articles/module-identity.html) |

The reference, [Functions and data](https://kanameishi.github.io/newmark/reference/index.html), groups the functions
by the same families and lists the data sets.

The library and the `newmark` command are installed together from the
repository root with the product installer.
