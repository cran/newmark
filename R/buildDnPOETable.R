# SPDX-FileCopyrightText: 2026 Alejandro Verri Kozlowski <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of newmark: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
#' Empirical displacement exceedance table of a site
#'
#' Counts the stored, pointwise `DnDraws` realizations above each admissible
#' displacement. For each demand, geometry, material and `ky`, each draw
#' receives its model's weight and `POE` is `sum(w[Dn > Da]) / sum(w)`.
#' The threshold is `DaH * Hs` in cm, with `DaH` in percent and the mean
#' shear-table height `Hs` in metres. No displacement or threshold is rounded.
#'
#' This fraction is conditional on the demand (return period or MCE) and
#' inherits the population represented by the draws: the positive lognormal
#' branch for the flexible models. Stored zeros, including SR08 underflow,
#' count as non-exceedance. It is not an integrated displacement hazard.
#'
#' @param draws Site `DnDraws`, with `siteID`, `ID`, `TR`, `IDg`, `IDm`,
#'   `IDn`, `ky`, `realization`, `Dn` and `units.Dn`. Every cell and ky must contain
#'   the `NS` realizations of every active model, with no duplicate identities.
#' @param shear ShearTable with `IDg`, `IDm`, `Hs` and `level`. Its mean
#'   rows define the geometry-material pairs and one height per geometry.
#' @param modelID Identifier of the return-period demand model.
#' @param siteID Site identifier.
#' @param TR Distinct positive return periods in years; `numeric()` selects
#'   MCE alone. MCE is always included, with `ID = "MCE"` and `TR = NA`.
#' @param NS Positive number of realizations per active model and cell.
#' @param weights Positive finite weights named by the active displacement
#'   models, as returned by [selectDisplacementModels()].
#' @param DaH Distinct positive displacement thresholds in percent of Hs.
#' @return A data.table with columns `siteID`, `ID`, `TR`, `IDg`, `IDm`,
#'   `IDn` (`"ensemble"`), `ky` (g), `DaH` (%), `Da` (cm), `POE`, `n`
#'   (integer count of all draws in the cell and ky) and `units.Da` (`"cm"`),
#'   keyed by `siteID`, `ID`, `TR`, `IDg`, `IDm`, `ky`, `DaH`.
#'   The inputs and random stream are unchanged; no files are written.
#' @export
buildDnPOETable <- function(draws, shear, modelID, siteID, TR, NS, weights, DaH) {
  Shear <- .checkDnPOEArguments(shear, modelID = modelID, TR = TR, weights = weights, DaH = DaH)
  if (!is.numeric(NS) || length(NS) != 1L || !is.finite(NS) || NS < 1 ||
      NS > .Machine$integer.max || NS != as.integer(NS)) {
    stop("NS must be a positive integer")
  }
  Draws <- as.data.table(draws)
  Missing <- setdiff(c("siteID", "ID", "TR", "IDg", "IDm", "IDn", "ky", "realization", "Dn", "units.Dn"), names(Draws))
  if (length(Missing)) stop("DnDraws is missing columns: ", toString(Missing))
  .checkDnPOECells(Draws, shear = Shear, siteID = siteID, modelID = modelID, TR = TR)
  if (!is.character(Draws$IDn) || anyNA(Draws$IDn) || !setequal(Draws$IDn, names(weights)) ||
      !is.numeric(Draws$Dn) || any(!is.finite(Draws$Dn)) || any(Draws$Dn < 0) ||
      anyNA(Draws$units.Dn) || any(Draws$units.Dn != "cm") ||
      !is.double(Draws$ky) || any(!is.finite(Draws$ky)) || any(Draws$ky <= 0) ||
      !is.numeric(Draws$realization) || any(!is.finite(Draws$realization)) ||
      any(Draws$realization < 1 | Draws$realization > NS | Draws$realization != trunc(Draws$realization))) {
    stop("DnDraws has invalid models, realizations, displacements, ky, or units")
  }
  if (anyDuplicated(Draws, by = c("siteID", "ID", "TR", "IDg", "IDm", "ky", "IDn", "realization"))) {
    stop("DnDraws has duplicate realization identities")
  }
  AUX <- Draws[, .N, by = .(ID, TR, IDg, IDm, ky, IDn)]
  if (any(AUX$N != NS) || any(AUX[, .N, by = .(ID, TR, IDg, IDm, ky)]$N != length(weights))) {
    stop("DnDraws must contain NS realizations of every active model in each cell and ky")
  }
  Thresholds <- sort(as.numeric(DaH))
  DT <- Draws[, {
    Targets <- Thresholds * Shear$Hs[match(.BY$IDg, Shear$IDg)]
    w <- weights[IDn]
    .(IDn = "ensemble", DaH = Thresholds, Da = Targets,
      POE = vapply(Targets, function(x) sum(w[Dn > x]) / sum(w), numeric(1L)), n = .N)
  }, by = .(siteID, ID, TR, IDg, IDm, ky)]
  DT[, units.Da := "cm"]
  setcolorder(DT, DNPOECOLS)
  .checkDnPOETable(DT, shear = Shear, siteID = siteID, modelID = modelID, TR = TR, DaH = DaH)
  DT
}

DNPOECOLS <- c("siteID", "ID", "TR", "IDg", "IDm", "IDn", "ky", "DaH", "Da", "POE", "n", "units.Da")
DNPOEKEY <- c("siteID", "ID", "TR", "IDg", "IDm", "ky", "DaH")

# Shared product context; unlike kmax, MCE alone is valid and each geometry
# contributes only its own unrounded heights from mean rows.
.checkDnPOEArguments <- function(shear, modelID, TR, weights, DaH) {
  if (!.isIdentifier(modelID) || modelID == "MCE") stop("modelID must be a demand identifier other than MCE")
  if (!is.numeric(TR) || any(!is.finite(TR)) || any(TR <= 0) || anyDuplicated(TR)) {
    stop("TR must contain distinct positive finite return periods, or be empty for MCE alone")
  }
  IDs <- names(weights)
  if (!is.numeric(weights) || !length(weights) || is.null(IDs) || anyNA(IDs) || anyDuplicated(IDs) ||
      any(!IDs %in% DNFAMILIES$ensemble) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must name active models with positive finite values")
  }
  if (!is.numeric(DaH) || !length(DaH) || any(!is.finite(DaH)) || any(DaH <= 0) || anyDuplicated(DaH)) {
    stop("DaH must contain distinct positive finite percentages of Hs")
  }
  if (!is.data.frame(shear) || !all(c("IDg", "IDm", "Hs", "level") %in% names(shear))) {
    stop("ShearTable must carry IDg, IDm, Hs and level")
  }
  # A positional subset avoids adding an automatic index to the caller.
  Shear <- as.data.table(shear)[which(shear$level == "mean"), .(IDg, IDm, Hs)]
  if (!nrow(Shear) || !is.character(Shear$IDg) || !is.character(Shear$IDm) || anyNA(Shear) ||
      any(!nzchar(Shear$IDg)) || any(!nzchar(Shear$IDm)) || !is.numeric(Shear$Hs) ||
      any(!is.finite(Shear$Hs)) || any(Shear$Hs <= 0) || anyDuplicated(Shear, by = c("IDg", "IDm")) ||
      anyDuplicated(unique(Shear[, .(IDg, Hs)]), by = "IDg")) {
    stop("ShearTable mean rows must have unique pairs and one positive finite Hs per geometry")
  }
  Shear
}

# Coverage shared by the draws and stored tables; TR = NA belongs to MCE.
.checkDnPOECells <- function(DT, shear, siteID, modelID, TR) {
  if (!.isIdentifier(siteID) || !identical(unique(DT$siteID), siteID) ||
      !all(vapply(DT[, .(ID, IDg, IDm)], is.character, logical(1L))) ||
      anyNA(DT[, .(ID, IDg, IDm)]) || !is.double(DT$TR)) {
    stop("DnPOE inputs have invalid site or cell identifiers/types")
  }
  Cells <- unique(DT[, .(ID, TR, IDg, IDm)])
  Demands <- data.table(ID = c(rep(modelID, length(TR)), "MCE"), TR = c(TR, NA_real_))
  if (nrow(Cells[!Demands, on = .(ID, TR)]) || nrow(Cells[!shear, on = .(IDg, IDm)]) ||
      nrow(Cells) != nrow(Demands) * nrow(shear)) {
    stop("DnPOE inputs do not cover the declared demands and geometry-material pairs")
  }
  invisible(NULL)
}

# Checks and keys an owned result. Assemblers pass a copy of caller tables.
.checkDnPOETable <- function(DT, shear, siteID, modelID, TR, DaH) {
  if (!identical(names(DT), DNPOECOLS)) stop("Invalid DnPOETable schema for siteID=", siteID)
  .checkDnPOECells(DT, shear = shear, siteID = siteID, modelID = modelID, TR = TR)
  if (!all(vapply(DT[, .(IDn, units.Da)], is.character, logical(1L))) ||
      !all(vapply(DT[, .(ky, DaH, Da, POE)], is.double, logical(1L))) ||
      anyNA(DT[, .(IDn, ky, DaH, Da, POE, n, units.Da)]) ||
      any(!is.finite(DT$ky) | DT$ky <= 0) || any(!is.finite(DT$Da) | DT$Da <= 0) ||
      any(!is.finite(DT$POE) | DT$POE < 0 | DT$POE > 1) ||
      !is.integer(DT$n) || any(DT$n <= 0) || any(DT$IDn != "ensemble") || any(DT$units.Da != "cm") ||
      !setequal(DT$DaH, DaH) || anyDuplicated(DT, by = DNPOEKEY) ||
      !identical(DT$Da, DT$DaH * shear$Hs[match(DT$IDg, shear$IDg)])) {
    stop("DnPOETable has invalid values, thresholds, units, or keys for siteID=", siteID)
  }
  AUX <- DT[, .(thresholds = .N, sizes = uniqueN(n)), by = .(ID, TR, IDg, IDm, ky)]
  if (any(AUX$thresholds != length(DaH) | AUX$sizes != 1L)) {
    stop("DnPOETable has incomplete thresholds or inconsistent n within a cell and ky")
  }
  setkeyv(DT, DNPOEKEY)
  invisible(DT)
}

utils::globalVariables(c("DaH", "POE", "units.Da"))
