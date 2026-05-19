## ── newmark: CRAN release workflow ───────────────────────────────────────────
## Run each section manually and in order. Do NOT source() end-to-end.
## Every step must pass before advancing to the next.
##
## This script reflects the flow that produced the 1.1.0 submission on
## 2026-05-09. Update any version-specific defaults (Step 1, Step 5)
## as you go, but the structure should hold for future patch / minor
## releases.

## ── Step 0: pre-flight ────────────────────────────────────────────────────────
## Confirm you are on main and the working tree is clean.
system("git checkout main && git pull --ff-only")
system("git status")
system("git log --oneline -5")

## ── Step 1: bump version ──────────────────────────────────────────────────────
## Edit DESCRIPTION manually before running anything else.
## Versioning rules:
##   bug fix only           → patch:  1.0.x
##   new exported function  → minor:  1.x.0
##   removed/renamed param  → major:  x.0.0
## CRAN policy: no more than one submission every ~2 months.

## ── Step 2: update NEWS.md ────────────────────────────────────────────────────
## Prepend a new top-level entry above the previous version:
##   # newmark 1.x.y
##   ## Bug fixes
##   * <change 1>
##   ## Documentation
##   * <change 2>
## Keep entries terse — one line per change.

## ── Step 3: regenerate documentation ─────────────────────────────────────────
devtools::document()
## Verify: man/ files updated, NAMESPACE correct.

## ── Step 4: local CRAN check ──────────────────────────────────────────────────
devtools::check(document = FALSE, cran = TRUE, force_suggests = TRUE)
## Must exit with 0 ERRORs, 0 WARNINGs, 0 NOTEs (or justify each NOTE).
## Common NOTE to justify on a first submission: "New submission" or
## "Maintainer field".

## ── Step 5: edit cran-comments.md ─────────────────────────────────────────────
## The file already exists at the repo root and is in .Rbuildignore
## (does not ship inside the tarball). Open it and update:
##
##   ## Test environments
##   * local: <macOS / R version of the day>
##   * GitHub Actions: ubuntu-latest, macos-latest, windows-latest (R-release)
##   * R-hub v2: linux, windows, macos-arm64        # if Step 6 was run
##
##   ## R CMD check results
##   0 errors | 0 warnings | 0 notes               # confirm from Step 4
##
##   ## Reverse dependencies
##   None on CRAN.                                 # or list affected
##
##   ## Notes for reviewer
##   <ONLY on resubmissions or unusual situations — explain what
##   changed in response to prior reviewer feedback. Delete on the
##   next clean submission.>

## ── Step 6: cross-platform check (optional) ──────────────────────────────────
## (a) R-hub v2 — triggers the .github/workflows/rhub.yaml workflow.
##     Requires RHUB_TOKEN secret in the repo.
# rhub::rhub_check(platforms = c("linux", "windows", "macos-arm64"))
##
## (b) win-builder — accepts a tarball at win-builder.r-project.org.
# devtools::check_win_devel()
# devtools::check_win_release()
##
## Monitor R-hub at: https://github.com/averriK/newmark/actions

## ── Step 7: commit, tag, push ─────────────────────────────────────────────────
## Stage: DESCRIPTION + NEWS.md + cran-comments.md + any man/, NAMESPACE
## changes from Step 3. Commit as a release commit, tag, push main + tag,
## then ff-sync dev.
# system('git add DESCRIPTION NEWS.md cran-comments.md man/ NAMESPACE')
# system('git commit -m "chore: release vX.Y.Z"')
# system('git tag -a vX.Y.Z -m "newmark X.Y.Z"')
# system('git push origin main vX.Y.Z')
# system('git checkout dev && git merge --ff-only main && git push origin dev && git checkout main')

## ── Step 8: submit to CRAN ────────────────────────────────────────────────────
## MUST be run from an INTERACTIVE R session (RStudio, or `R` in a
## terminal — not `R -e` / `Rscript`). devtools::submit_cran() prompts
## once to confirm the maintainer email; in non-interactive contexts
## it aborts with "Called from non-interactive context".
##
## What it does:
## 1. Builds the tarball.
## 2. Uploads it to https://xmpalantir.wu.ac.at/cransubmit/.
## 3. CRAN-SUBMISSION at the repo root is updated with the new
##    Date / SHA stamp (commit this small change after the run).
## 4. Sends a confirmation email to the maintainer.
devtools::submit_cran()

## ── Step 9: confirm submission via email ─────────────────────────────────────
## Open the email from "CRAN <CRAN@R-project.org>" with subject
## "[CRAN-pretest] newmark X.Y.Z" or similar, and follow the link to
##   https://xmpalantir.wu.ac.at/cransubmit/conf_mail.php?code=...
##
## The page is "Step 3 of 3" — a confirmation form, NOT a re-upload.
## Tick the two checkboxes:
##   [ ] I have read the CRAN policies and would like to submit this
##       package to CRAN.
##   [ ] I have checked the submission using R CMD check --as-cran and
##       a current version of r-devel, as mandated by the CRAN
##       Repository Policy.
## Then click submit. The package is now in CRAN's review queue.

## ── Step 10: wait for CRAN response ──────────────────────────────────────────
## CRAN typically responds within 1-5 business days by email:
##   * accepted   → no further action; the package appears at
##                  https://cran.r-project.org/package=newmark
##   * pre-test   → automated checks failed; address and resubmit
##                  (incrementing patch is usually NOT required if you
##                  are still pre-acceptance).
##   * manual     → reviewer requested changes; address, edit
##                  cran-comments.md "Notes for reviewer" with what you
##                  changed, and resubmit.

## ── Step 11: post-acceptance ──────────────────────────────────────────────────
## 1. Verify badge resolves: https://www.r-pkg.org/badges/version/newmark
## 2. The pkgdown workflow auto-deploys https://averriK.github.io/newmark/
##    on every push to main; nothing to do unless docs need a rebuild.
## 3. Commit the CRAN-SUBMISSION update from Step 8.3 if not yet committed.
