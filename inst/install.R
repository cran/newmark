## ── newmark: developer workflow ──────────────────────────────────────────────
## Run sections manually as needed.  Do NOT source() this file end-to-end.
## For CRAN release workflow see inst/cran-release.R.

## ── Install ───────────────────────────────────────────────────────────────────
## From CRAN (stable release — update only when a new version is accepted):
# install.packages("newmark")

## From GitHub (latest commit on main — picks up every push immediately):
# remotes::install_github("averriK/newmark")
# pak::pak("averriK/newmark")   # alternative if pak is available

## ── Development cycle ─────────────────────────────────────────────────────────
## 1. Regenerate documentation (roxygen2 → NAMESPACE + man/)
devtools::document()

## 2. Local CRAN check
devtools::check(document = FALSE, cran = TRUE, force_suggests = TRUE)

## 3. Local install with vignettes
devtools::install(build_vignettes = TRUE, dependencies = TRUE,
                  build = TRUE, upgrade = TRUE)

## 4. Render pkgdown site locally (optional sanity check before pushing)
# pkgdown::build_site()
