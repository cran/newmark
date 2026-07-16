# newmark

<!-- badges: start -->
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.2.0-blue)](https://www.r-project.org/)
[![CRAN status](https://www.r-pkg.org/badges/version/newmark)](https://CRAN.R-project.org/package=newmark)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/newmark)](https://CRAN.R-project.org/package=newmark)
<!-- badges: end -->

## Overview

The package implements probabilistic propagation of seismic hazard
uncertainty through the Newmark sliding-block displacement framework.
Two stages of variability — the joint distribution of rock-level
spectral demand from PSHA quantiles, and the conditional dispersion
of empirical site-amplification and sliding-block models — are
sampled via Gaussian copulas with mean-preserving recentering, and
propagated to the mean, dispersion and design quantiles of the
residual displacement Dn(ky) and the performance-based seismic
coefficient kmax(d*).

Inputs are uniform-hazard spectra (UHS) at multiple return periods,
produced by any external PSHA engine.
Outputs are tabulated kmax(d*) and the normalised pseudostatic
coefficient Kh = kmax / PGA at user-specified displacement targets,
intended for performance-based seismic stability of slopes,
embankments, tailings storage facilities and waste rock dumps.

## Installation

```r
# Stable release
install.packages("newmark")

# Development version
devtools::install_github("averriK/newmark")
```

## Getting started

Four vignettes ship with the package:

```r
vignette("newmark-quickstart",     package = "newmark")    # runnable hello world
vignette("pipeline",               package = "newmark")    # function-level overview
vignette("dynamic-site-response",  package = "newmark")    # Ts, stiffness profile, site amplification
vignette("ensemble-formulation",   package = "newmark")    # mathematical derivation
```

Function-level reference: `?fitSaF`, `?getDnKy`, `?fitDnCurve`,
`?invertDnDraws`, `?getSiteProperties`, `?geSiteTable`,
`?getCylinderRoots`, `?fitModel.Ts`, `?buildGMDP`, `?Vs30toSID`,
`?SIDtoVs30`, or `help(package = "newmark")`.

Rendered documentation: <https://averriK.github.io/newmark/>.

## Getting help

Open an issue at <https://averriK.github.io/newmark/>.

## Citation

If you use this package in a publication, please cite the companion paper:

```bibtex
@unpublished{VerriKozlowski2026,
  author = {Verri Kozlowski, Alejandro},
  title  = {Probabilistic estimation of {Newmark} displacements and seismic
            coefficients under hazard uncertainty},
  year   = {2026},
  note   = {Working paper}
}
```

To cite the software package itself:

```bibtex
@Manual{Verri2026,
  title  = {{newmark}: Uncertainty Analysis in Dynamic Site and Slope Response},
  author = {Verri Kozlowski, Alejandro},
  year   = {2026},
  note   = {R package version 1.1.0},
  url    = {https://averriK.github.io/newmark/}
}
```

## License

MIT — see [LICENSE](LICENSE)

## Author

**Alejandro Verri Kozlowski** — averri@fi.uba.ar  
[ORCID 0000-0002-8535-1170](https://orcid.org/0000-0002-8535-1170)  
Universidad de Buenos Aires, Facultad de Ingeniería
