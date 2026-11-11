# SPDX-FileCopyrightText: 2026 Alejandro Verri Kozlowski <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of newmark: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
#' Shear-table realizations of release 1.7.1 with their beam periods
#'
#' Draws the population of soil columns that [buildShearTable()] summarizes
#' and keeps it: for every geometry crossed with every selected material (or
#' only the supplied `pairs`), the realizations that [buildShearTable()]
#' simulates, each with the first-mode period of the truncated shear beam
#' evaluated on that realization, together with the shear table summarized
#' from the same draws.
#'
#' Realizations 1 to 100 of every pair are drawn pair by pair, in the order
#' of the table, with the calls, the rounding and the random numbers of
#' [buildShearTable()], and the shear table summarizes exactly those 100: under
#' the same random stream it is the table that [buildShearTable()] returns.
#' When `NS` exceeds 100, realizations 101 to `NS` of every pair are drawn
#' after all the pairs have their first 100, in a second pass over the pairs
#' in the same order, so the first 100 and the shear table do not depend on
#' `NS`.
#'
#' The period `Ts` of a realization is the beam expression of
#' [evaluateBeamPeriod()] evaluated on its own `mo` and `VSo` with the
#' truncation ratio `lo` of its geometry, rounded as on a summary row of the
#' shear table: `lo` to 3 decimals, the first-mode root `an` (the `"nlm"`
#' surrogate of [buildRootModel()] within its tabulated support, the exact
#' root outside it) to 4 and `Ts` to 3. It is neither the Rayleigh period of
#' the realization nor a statistic of the summary rows. The other
#' observables of a realization are the ones [buildShearTable()] summarizes,
#' with the rounding of the `newmark` 1.7.1 generator: `Go` to 2 decimals,
#' `mo` to 3, and `VSo`, `VS30`, `Z500` and `Z1000` to 1. As a rule of this
#' product, a realization whose fitted `mo` is negative, a column whose
#' stiffness fit decreases with depth, is returned with `mo` forced to 0, the
#' homogeneous beam, and with `an` and `Ts` evaluated at that 0, and a
#' message reports how many realizations were forced, while the shear table
#' summarizes the fitted `mo`. The root and the period are defined for `mo`
#' below 2: a realization at or above it stops the function naming it before
#' any period is evaluated.
#'
#' @inheritParams buildShearTable
#' @param NS Realizations per pair, a whole number of at least 100.
#' @return A list of two data.tables. `ShearDraws` has one row per pair and
#'   realization, ordered by pair as the shear table and then by
#'   realization, with the columns `Hs`, `b`, `s`, `Vref`, `IDm`, `IDg`,
#'   `POP`, `Hw`, `USCS`, `Go`, `mo`, `Ts` (fundamental period of the shear
#'   beam of the realization, s), `VSo`, `VS30`, `Z500`, `Z1000`,
#'   `realization` (an integer from 1), `lo` and `an`. `ShearTable` is the
#'   table that [buildShearTable()] describes, summarized from realizations 1
#'   to 100 of every pair.
#' @export
#' @examples
#' \donttest{
#' Geometry <- buildGeometry(Hs = 10, lo = 0.35, s = 2.5)
#' set.seed(20260915L)
#' Shear <- buildShearDraws(Geometry, USCS = list(U1 = c("ML", "MH")), IDm = "U1",
#'                          Vref = 760, levels = c(0.16, "mean", 0.84))
#' Shear$ShearDraws[, .(realization, mo, VSo, VS30, Ts)]
#' Shear$ShearTable
#' }
buildShearDraws <- function(geometry, USCS, IDm, Vref, levels, pairs = NULL, NS = 100) {
  if (!is.numeric(NS) || length(NS) != 1L || !is.finite(NS) || NS < 100 || NS != trunc(NS)) {
    stop("NS must be one whole number of at least 100")
  }
  Problem <- .diagnoseShearArguments(geometry, USCS = USCS, IDm = IDm, Vref = Vref, levels = levels, pairs = pairs)
  if (!is.null(Problem)) stop(Problem)
  Shear <- .drawShearTable(geometry, USCS = USCS, IDm = IDm, Vref = Vref, levels = levels, pairs = pairs, NS = NS)
  # The first five realizations of `rows`, as the messages name them.
  nameRealizations <- function(rows) {
    paste(head(rows[, paste0("IDg=", IDg, " IDm=", IDm, " realization=", realization, " mo=", mo)], 5L),
          collapse = "; ")
  }
  # The root of the characteristic relation and the beam period are defined
  # for mo below 2; a realization beyond it is named before any period of
  # the summary rows is evaluated.
  Beyond <- Shear$ShearDraws[mo >= 2]
  if (nrow(Beyond)) {
    stop(nrow(Beyond), " realization(s) have mo of 2 or more, where the first-mode root is not defined: ",
         nameRealizations(Beyond))
  }
  .evaluateLegacyPeriod(Shear$ShearTable)
  # A rule of this product: a negative fit, a stiffness that decreases with
  # depth, is published as the homogeneous beam, mo = 0; the summary rows
  # keep the fitted mo.
  Negative <- Shear$ShearDraws[mo < 0]
  Shear$ShearDraws[mo < 0, mo := 0]
  .evaluateLegacyPeriod(Shear$ShearDraws)
  if (nrow(Negative)) {
    .inform("moduleNote", "Negative mo forced to 0 in ", nrow(Negative), " realization(s): ", nameRealizations(Negative),
            operation = "buildShearDraws", count = nrow(Negative))
  }
  Shear
}

# The realizations of every pair and the shear table summarized from the
# first 100, as buildShearTable() and buildShearDraws() describe them, for
# arguments that their caller has checked, before any root is evaluated:
# the column Ts of both tables is still the rounded Rayleigh period, which
# the summary takes as newmark 1.7.1 did and the caller replaces with the
# beam period of .evaluateLegacyPeriod().
.drawShearTable <- function(geometry, USCS, IDm, Vref, levels, pairs, NS) {
  # Resolved in the function scope, so a USCS column of the caller's table
  # cannot shadow the material pools inside the grouped evaluation.
  materialPool <- function(id) USCS[[id]]
  drawProfiles <- function(Hs, IDm, NS) {
    Profiles <- .legacyProfile(simulateSoilProfiles(Hs = Hs, USCS = materialPool(IDm), NS = NS, h = 1,
                                                    water = 0, POP = 100, Vref = 760,
                                                    m1 = data.frame(IP = c(0, 20, 40, 60, 80, 100),
                                                                    m1 = c(0, 0.18, 0.30, 0.41, 0.48, 0.48))))
    data.table(POP = 100, Hw = 0, USCS = paste0(materialPool(IDm), collapse = " "),
               Profiles[, .(Go, mo, Ts, VSo, VS30, Z500, Z1000, realization)])
  }
  Pairs <- CJ(IDm, IDg = geometry$IDg)
  if (!is.null(pairs)) Pairs <- as.data.table(pairs)[, .(IDm, IDg)]
  # Quantile rows precede the mean row, as the original summary did.
  Levels <- c(levels[levels != "mean"], if (any(levels == "mean")) "mean")
  MaterialTable <- as.data.table(geometry)[Pairs, on = .(IDg), allow.cartesian = TRUE]
  set(MaterialTable, j = "Vref", value = Vref)
  # Each pair is summarized before the next one is drawn, as buildShearTable()
  # summarized them before the draws were kept, so that an error stops the
  # random stream at the same number; the realizations of a pair ride on its
  # summary rows as the list column draws until every pair is summarized.
  ShearTable <- MaterialTable[, {
    Profiles <- drawProfiles(Hs, IDm, NS = 100L)
    OUT <- summarizeProfiles(Profiles[rev(seq_len(nrow(Profiles)))], levels = Levels,
                             properties = c("Go", "mo", "Ts", "VSo", "VS30", "Z500", "Z1000"))
    OUT[, `:=`(Go = round(Go, 1), mo = round(mo, 3), Ts = round(Ts, 2), VSo = round(VSo, 1),
               VS30 = round(VS30, 1), Z500 = round(Z500, 1))]
    OUT[level != "mean", Z1000 := round(Z1000, 1)]
    data.table(Profiles[1L, .(POP, Hw, USCS)], OUT[, .(Go, mo, Ts, VSo, VS30, Z500, Z1000, level)],
               draws = list(Profiles))
  }, by = .(Hs, b, s, Vref, IDm, IDg)]
  ShearDraws <- ShearTable[, draws[[1L]], by = .(Hs, b, s, Vref, IDm, IDg)]
  set(ShearTable, j = "draws", value = NULL)
  if (NS > 100) {
    Extension <- MaterialTable[, drawProfiles(Hs, IDm, NS = NS - 100L), by = .(Hs, b, s, Vref, IDm, IDg)]
    set(Extension, j = "realization", value = Extension$realization + 100L)
    ShearDraws <- rbind(ShearDraws, Extension)[, .SD[order(realization)], by = .(Hs, b, s, Vref, IDm, IDg)]
  }
  list(ShearDraws = ShearDraws[], ShearTable = ShearTable[])
}

# The beam period of every row of a shear table or of its realizations, set
# by reference with the rounding of newmark 1.7.1: the truncation ratio lo
# of the geometry to 3 decimals, the first-mode root an to 4 (the nlm
# surrogate within its support, the exact root outside it, where 1.7.1
# extrapolated the surrogate) and the period Ts to 3.
.evaluateLegacyPeriod <- function(shear) {
  shear[, lo := round(b / (2 * Hs * s + b), 3)]
  Roots <- evaluateRoot(buildRootModel(mode = 1, method = "nlm"), mo = shear$mo, lo = shear$lo)
  Roots[!Roots$withinSupport, an := evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = mo, lo = lo)$an]
  shear[, an := round(Roots$an, 4)]
  shear[, Ts := round(evaluateBeamPeriod(Hs = Hs, lo = lo, mo = mo, VSo = VSo, an = an), 3)]
  invisible(shear)
}

# The list column of .drawShearTable(), resolved by data.table.
utils::globalVariables("draws")
