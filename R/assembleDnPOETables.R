# SPDX-FileCopyrightText: 2026 Alejandro Verri Kozlowski <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of newmark: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
#' Assemble the displacement exceedance tables of sites
#'
#' Validates and stacks the site products of [buildDnPOETable()]. Each site
#' must have the declared demands, geometry-material pairs and thresholds.
#' The weighted counts are already in the site tables and are not recomputed.
#'
#' @param tables Named list of site DnPOETable products, one per site identifier.
#' @param shear,modelID,TR,weights,DaH As for [buildDnPOETable()]. Weights
#'   define the active models of the contract; stored tables do not retain
#'   individual model counts or weights, so their values cannot be recovered
#'   or checked from the aggregate POE.
#' @return A data.table with the columns, key and meanings of
#'   [buildDnPOETable()], holding every supplied site. Inputs and the random
#'   stream are unchanged; no files are written.
#' @export
assembleDnPOETables <- function(tables, shear, modelID, TR, weights, DaH) {
  IDs <- names(tables)
  if (!is.list(tables) || !length(tables) || is.null(IDs) || anyDuplicated(IDs) ||
      !all(vapply(IDs, .isIdentifier, logical(1L)))) {
    stop("tables must be a named list with unique site identifiers")
  }
  Shear <- .checkDnPOEArguments(shear, modelID = modelID, TR = TR, weights = weights, DaH = DaH)
  DT <- rbindlist(lapply(IDs, function(siteID) {
    .checkDnPOETable(copy(as.data.table(tables[[siteID]])), shear = Shear, siteID = siteID,
                      modelID = modelID, TR = TR, DaH = DaH)
  }), use.names = TRUE)
  setkeyv(DT, DNPOEKEY)
  DT
}
