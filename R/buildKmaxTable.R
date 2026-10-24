#' Yield-coefficient table of release 1.7.1 from displacement realizations
#'
#' Compatibility route that reproduces, byte for byte once saved with
#' [saveRDS()], the site `kmaxTable` product of the `fitKmaxProducts()`
#' function of `newmark` 1.7.1, so that the stored project tables can serve
#' as comparison oracle. For every demand (each return period of `TR` and
#' the MCE), geometry and material of the site, the realizations of every
#' active model are inverted with [solveKmax()] and summarized, model by
#' model and for each model family, with the pooled summary of
#' [summarizeKmax()]. With one site, as in the current projects, this table
#' is also the global `kmaxTable` of 1.7.1.
#'
#' @section Semantics reproduced from `newmark` 1.7.1:
#' Every geometry is inverted at the same levels, the union over all the
#' heights of the shear table of \eqn{D_a = D_{aH} H_s} cm rounded to six
#' significant digits, so a geometry's own levels are the rows whose `Da`
#' equals `signif(DaH * Hs, 6)` for its height. Each model is summarized
#' alone with weight 1; the families `ensemble`, `rigid` and `flexible` of
#' the displacement models pool their active members with their weights
#' normalized, and normalized once more by the summary.
#' The levels are the mean and the probabilities 0.05, 0.10, 0.16, 0.50,
#' 0.84, 0.90 and 0.95. The quantiles are those of `Hmisc::wtd.quantile()`
#' as 1.7.1 called it: once the pool of one `Da` holds a tie (a saturated
#' root of 0 or a root at the grid resolution), every pooled value is read
#' back from its 15-significant-digit text before the quantile, while the
#' mean uses the exact roots; [summarizeKmax()] has no such text step.
#' `PGA` is the amplified peak ground acceleration of the cell as read back
#' from the CSV of the PGA table, `Kh = kmax / PGA` is a dimensionless ratio
#' (the exported `fitKmaxScenario()` of 1.7.1 returned `100 kmax / PGA`,
#' which its product replaced by the ratio), and `Ts` and `Vs30` are the
#' fundamental period of the shear beam and `VS30` of the mean row of the geometry and material in
#' the shear table. The inputs are checked as 1.7.1 checked them (schema, 30
#' distinct `ky` per curve spanning 0.01 g to `max(PGA, 0.8)` g, one sampled
#' PGA per realization, `NS` realizations of every active model, coverage of
#' the PGA table); the result is checked and keyed by `siteID`, `ID`, `TR`,
#' `IDg`, `IDm`, `IDn`, `p` and `Da`, and a warning reports the rows that
#' rely on roots beyond the stored grid and the largest such `kmax`.
#'
#' @param draws The `DnDraws` product of the site: `IDg`, `IDm`, `ID`,
#'   `TR`, `ky`, `s`, `PGA.s`, `IDn`, `Dn`, `siteID` and `units.Dn`, in this
#'   order, the realizations of one cell sharing one `ky` grid.
#' @param PGATable The `PGATable` product of the site as read back from its CSV:
#'   `IDg`, `IDm`, `ID`, `TR`, `PGA`, `Vs30`, `siteID` and `units.PGA`, in this order, with `TR`
#'   numeric and `NA` on the MCE rows. It is named after the product because `PGA` names the numeric acceleration elsewhere.
#' @param shear The shear table of [buildShearTable()]; its heights
#'   `Hs` fix the levels and its `mean` rows supply `Ts` and `VS30`.
#' @param modelID Identifier of the hazard model, the `ID` of the
#'   return-period demands.
#' @param siteID Identifier of the site.
#' @param TR Positive return periods of the model demands, in years.
#' @param NS Number of realizations of every curve.
#' @param weights Positive weights named by the active displacement models.
#' @param DaH Distinct positive admissible displacements in percent of the
#'   height.
#' @return A data.table keyed by `siteID`, `ID`, `TR`, `IDg`, `IDm`, `IDn`,
#'   `p` and `Da`, with columns `IDg`, `IDm`, `p` (`"mean"` or the
#'   probability as `"%.2f"`), `ID` (`modelID` or `"MCE"`), `TR` (years,
#'   `NA` for the MCE), `Da` (cm), `kmax` (g), the mean or quantile of the
#'   yield coefficient that limits the displacement to `Da` over the pooled
#'   realizations of the model or family `IDn` under that demand, geometry
#'   and material, `PGA` (g), `Kh` (ratio), `beyond` (weighted share of the
#'   pool beyond the stored grid), `Ts` (s), `Vs30` (m/s), `IDn`, `siteID`
#'   and the unit columns `units.Da` (`"cm"`), `units.kmax` and `units.PGA`
#'   (`"g"`) and `units.Kh` (`"ratio"`). The caller's tables are not
#'   modified.
#' @export
buildKmaxTable <- function(draws, PGATable, shear, modelID, siteID, TR, NS, weights, DaH) {
  .checkKmaxArguments(modelID = modelID, siteID = siteID, TR = TR, NS = NS, weights = weights)
  Problem <- .diagnoseKmaxLevels(DaH, shear)
  if (!is.null(Problem)) stop(Problem)
  Da <- .deriveKmaxLevels(DaH, shear)
  Periods <- sort(unique(as.numeric(TR)))
  Shear <- .selectKmaxShear(shear)
  PGA <- as.data.table(PGATable)
  Draws <- as.data.table(draws)
  .checkKmaxPGA(PGA, shear = Shear, siteID = siteID, modelID = modelID, periods = Periods)
  .checkKmaxDraws(Draws, shear = Shear, PGA = PGA, siteID = siteID, modelID = modelID, periods = Periods,
                  weights = weights, NS = NS)
  Rows <- Draws[, .summarizeKmaxScenario(.SD, weights = weights, Da = Da), by = .(ID, TR, IDg, IDm),
                .SDcols = c("ky", "s", "PGA.s", "Dn", "IDn")]
  Rows[PGA, on = .(ID, TR, IDg, IDm), PGA := i.PGA]
  Rows[Shear, on = .(IDg, IDm), `:=`(Ts = i.Ts, Vs30 = i.VS30)]
  # Built by one j call, so the attributes serialize in the order of the
  # stored product: the RDS bytes, not only the values, are compared.
  DT <- Rows[, .(IDg, IDm, p = level, ID, TR, Da, kmax, PGA, Kh = kmax / PGA, beyond, Ts, Vs30, IDn, siteID = siteID,
                 units.Da = "cm", units.kmax = "g", units.PGA = "g", units.Kh = "ratio")]
  .checkKmaxTable(DT, shear = Shear, PGA = PGA, siteID = siteID, modelID = modelID, periods = Periods,
                  weights = weights, Da = Da)
  .warnKmaxBeyond(DT, label = paste0("siteID=", siteID))
  DT
}

# Columns and key of the kmaxTable product of newmark 1.7.1.
KMAXCOLS <- c("IDg", "IDm", "p", "ID", "TR", "Da", "kmax", "PGA", "Kh", "beyond", "Ts", "Vs30", "IDn",
              "siteID", "units.Da", "units.kmax", "units.PGA", "units.Kh")
KMAXKEY <- c("siteID", "ID", "TR", "IDg", "IDm", "IDn", "p", "Da")

# Probability levels of the kmax product (PTARGET of newmark 1.7.1).
KMAXLEVELS <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)

# One demand, geometry and material: the realizations of each active model
# as a matrix over the shared ky grid, their roots, and the summaries of
# every model alone and of every model family with an active member. The
# family membership is DNFAMILIES of the displacement family; the family
# weights are normalized before the pool normalizes them again, as 1.7.1
# did, which changes the bits when all six models are active.
.summarizeKmaxScenario <- function(draws, weights, Da) {
  Missing <- setdiff(names(weights), unique(as.character(draws$IDn)))
  if (length(Missing)) stop("kmax draws are missing positive-weight models: ", paste(Missing, collapse = ", "))
  Roots <- rbindlist(lapply(names(weights), function(m) {
    Curves <- draws[IDn == m]
    S <- sort(unique(Curves$s))
    Ky <- sort(unique(Curves$ky))
    Dn <- matrix(NA_real_, nrow = length(S), ncol = length(Ky))
    Dn[cbind(match(Curves$s, S), match(Curves$ky, Ky))] <- Curves$Dn
    if (anyNA(Dn)) stop("the realizations of model ", m, " in one cell do not share one ky grid")
    solveKmax(ky = Ky, Dn = Dn, PGA = Curves$PGA.s[match(S, Curves$s)], model = m, Da = Da)
  }))
  Models <- lapply(names(weights), function(m) {
    .poolKmax(Roots[IDn == m], weights = stats::setNames(1, m), levels = KMAXLEVELS,
              estimator = .textWeightedQuantile)[, IDn := m]
  })
  Pools <- lapply(names(DNFAMILIES), function(family) {
    Members <- intersect(DNFAMILIES[[family]], names(weights))
    if (!length(Members)) return(NULL)
    .poolKmax(Roots[IDn %in% Members], weights = weights[Members] / Reduce(`+`, weights[Members]), levels = KMAXLEVELS,
              estimator = .textWeightedQuantile)[, IDn := family]
  })
  rbindlist(c(Models, Pools))
}

.isKmaxIdentifier <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[A-Za-z][A-Za-z0-9_-]*$", x)
}

.checkKmaxArguments <- function(modelID, siteID, TR, NS, weights) {
  if (!.isKmaxIdentifier(modelID) || !.isKmaxIdentifier(siteID)) stop("modelID and siteID must be valid identifiers")
  if (!is.numeric(TR) || !length(TR) || any(!is.finite(TR)) || any(TR <= 0)) {
    stop("TR must contain positive finite return periods")
  }
  if (!is.numeric(NS) || length(NS) != 1L || !is.finite(NS) || NS < 1 || NS != as.integer(NS)) {
    stop("NS must be a positive integer")
  }
  IDs <- names(weights)
  if (!is.numeric(weights) || !length(weights) || is.null(IDs) || anyNA(IDs) || anyDuplicated(IDs) ||
      any(!IDs %in% DNFAMILIES$ensemble) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must name active models with positive finite values")
  }
  invisible(NULL)
}

# The mean rows of the shear table, one per geometry and material. The copy
# keeps the caller's table free of the index a subset would attach to it.
.selectKmaxShear <- function(shear) {
  DT <- copy(as.data.table(shear))
  Missing <- setdiff(c("IDg", "IDm", "VS30", "Ts", "level"), names(DT))
  if (length(Missing)) stop("ShearTable is missing columns: ", paste(Missing, collapse = ", "))
  Shear <- DT[level == "mean", .(IDg = as.character(IDg), IDm = as.character(IDm), VS30 = as.numeric(VS30),
                                 Ts = as.numeric(Ts))]
  if (!nrow(Shear)) stop("ShearTable has no rows with level == 'mean'")
  if (anyNA(Shear) || any(!nzchar(Shear$IDg)) || any(!nzchar(Shear$IDm))) {
    stop("ShearTable mean rows contain missing or empty identifiers/values")
  }
  if (any(!is.finite(Shear$VS30)) || any(Shear$VS30 <= 0) || any(!is.finite(Shear$Ts)) || any(Shear$Ts <= 0)) {
    stop("ShearTable mean rows require finite positive VS30 and Ts")
  }
  if (anyDuplicated(Shear, by = c("IDg", "IDm"))) stop("ShearTable mean rows are not unique by IDg x IDm")
  Shear
}

# The demands of a product: the model with exactly `periods` and the MCE
# with a missing return period.
.hasKmaxDemands <- function(DT, modelID, periods) {
  setequal(unique(DT$ID), c(modelID, "MCE")) && setequal(unique(DT$TR[DT$ID == modelID]), periods) &&
    all(is.na(DT$TR[DT$ID == "MCE"]))
}

# Relative agreement within 1e-10, the product tolerance of newmark 1.7.1.
.isSameKmaxNumber <- function(x, y) {
  length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
    all(abs(x - y) <= 1e-10 * pmax(1, abs(x), abs(y)))
}

# Whether the distinct values of `cols` in `product` reproduce, cell for
# cell, the columns of `reference` they name.
.reproducesKmax <- function(product, reference, by, cols) {
  Observed <- unique(product[, c(by, names(cols)), with = FALSE])
  if (nrow(Observed) != nrow(reference) || anyDuplicated(Observed, by = by)) return(FALSE)
  Expected <- reference[, c(by, unname(cols)), with = FALSE]
  setnames(Expected, unname(cols), paste0(names(cols), ".reference"))
  AUX <- merge(Observed, Expected, by = by, sort = FALSE)
  nrow(AUX) == nrow(reference) &&
    all(vapply(names(cols), function(x) .isSameKmaxNumber(AUX[[x]], AUX[[paste0(x, ".reference")]]), logical(1L)))
}

.checkKmaxPGA <- function(PGA, shear, siteID, modelID, periods) {
  Columns <- c("IDg", "IDm", "ID", "TR", "PGA", "Vs30", "siteID", "units.PGA")
  if (!identical(names(PGA), Columns)) {
    stop("Invalid PGATable schema for siteID=", siteID, ". Expected: ", paste(Columns, collapse = ", "),
         ". Found: ", paste(names(PGA), collapse = ", "))
  }
  if (!identical(unique(PGA$siteID), siteID) || !.hasKmaxDemands(PGA, modelID, periods) ||
      anyDuplicated(PGA, by = c("siteID", "ID", "TR", "IDg", "IDm")) ||
      nrow(PGA) != (length(periods) + 1L) * nrow(shear)) {
    stop("PGATable has invalid key or demand x IDg x IDm coverage")
  }
  if (anyNA(PGA[, .(IDg, IDm, ID, siteID)]) || any(!nzchar(PGA$IDg)) || any(!nzchar(PGA$IDm)) ||
      any(!is.finite(PGA$PGA)) || any(PGA$PGA <= 0) || any(!is.finite(PGA$Vs30)) || any(PGA$Vs30 <= 0) ||
      any(PGA$units.PGA != "g")) {
    stop("PGATable has invalid identifiers, numerical values, or units")
  }
  if (!.reproducesKmax(PGA, shear, by = c("IDg", "IDm"), cols = c(Vs30 = "VS30"))) {
    stop("PGATable Vs30 does not reproduce ShearTable for siteID=", siteID)
  }
  invisible(NULL)
}

# The DnDraws contract of newmark 1.7.1: 30 distinct ky per curve from
# 0.01 g to max(PGA, 0.8) g of its cell, one sampled PGA per realization and
# NS realizations of every active model in every cell.
.checkKmaxDraws <- function(DT, shear, PGA, siteID, modelID, periods, weights, NS) {
  Columns <- c("IDg", "IDm", "ID", "TR", "ky", "s", "PGA.s", "IDn", "Dn", "siteID", "units.Dn")
  if (!identical(names(DT), Columns)) {
    stop("Invalid DnDraws schema for siteID=", siteID, ". Expected: ", paste(Columns, collapse = ", "),
         ". Found: ", paste(names(DT), collapse = ", "))
  }
  Cells <- (length(periods) + 1L) * nrow(shear) * length(weights)
  if (!identical(unique(DT$siteID), siteID) || !.hasKmaxDemands(DT, modelID, periods) ||
      !setequal(unique(DT$IDn), names(weights))) {
    stop("DnDraws has invalid siteID, ID, TR, or IDn vocabulary")
  }
  if (anyNA(DT$IDg) || anyNA(DT$IDm) || anyNA(DT$ID) || anyNA(DT$IDn) || anyNA(DT$siteID) ||
      any(!nzchar(DT$IDg)) || any(!nzchar(DT$IDm)) || any(!is.finite(DT$ky)) || any(DT$ky <= 0) ||
      any(!is.finite(DT$s)) || any(DT$s != as.integer(DT$s)) || any(!is.finite(DT$PGA.s)) ||
      any(DT$PGA.s <= 0) || any(!is.finite(DT$Dn)) || any(DT$Dn < 0) || any(DT$units.Dn != "cm")) {
    stop("DnDraws has invalid identifiers, numerical values, or units")
  }
  if (!setequal(unique(DT$s), seq_len(NS))) stop("DnDraws simulation indices do not match NS")
  if (any(DT[, uniqueN(PGA.s), by = .(ID, TR, IDg, IDm, s)]$V1 != 1L)) {
    stop("DnDraws PGA.s must be constant within each sampled demand")
  }
  if (anyDuplicated(DT, by = c("siteID", "ID", "TR", "IDg", "IDm", "IDn", "s", "ky"))) {
    stop("DnDraws key is duplicated for siteID=", siteID)
  }
  if (nrow(DT) != Cells * NS * 30L) {
    stop("DnDraws row count mismatch for siteID=", siteID, ": expected ", Cells * NS * 30L, ", found ", nrow(DT))
  }
  Curves <- DT[, .(n = .N, nky = uniqueN(ky), ky.min = min(ky), ky.max = max(ky)), by = .(ID, TR, IDg, IDm, IDn, s)]
  if (any(Curves$n != 30L) || any(Curves$nky != 30L)) {
    stop("DnDraws does not have 30 unique ky values per simulation curve")
  }
  Counts <- Curves[, .(ns = uniqueN(s)), by = .(ID, TR, IDg, IDm, IDn)]
  if (any(Counts$ns != NS) || nrow(Counts) != Cells) {
    stop("DnDraws does not have NS curves for every demand x IDg x IDm x IDn")
  }
  Expected <- PGA[, .(ID, TR, IDg, IDm, ky.min = 0.01, ky.max = pmax(PGA, 0.8))]
  if (!.reproducesKmax(Curves, Expected, by = c("ID", "TR", "IDg", "IDm"), cols = c(ky.min = "ky.min", ky.max = "ky.max"))) {
    stop("DnDraws ky range does not match PGATable for siteID=", siteID)
  }
  invisible(NULL)
}

# Checks one site's kmaxTable against its shear table and, unless PGA is NULL (the global table), its PGA table,
# and keys it by reference; used on the table just built and on stored site tables.
.checkKmaxTable <- function(DT, shear, PGA, siteID, modelID, periods, weights, Da) {
  if (!identical(names(DT), KMAXCOLS)) {
    stop("Invalid kmaxTable schema for siteID=", siteID, ". Expected: ", paste(KMAXCOLS, collapse = ", "),
         ". Found: ", paste(names(DT), collapse = ", "))
  }
  Levels <- c("mean", sprintf("%.2f", KMAXLEVELS))
  Families <- names(DNFAMILIES)[vapply(DNFAMILIES, function(ids) any(names(weights) %in% ids), logical(1L))]
  Models <- c(names(weights), Families)
  if (!identical(unique(DT$siteID), siteID) || !.hasKmaxDemands(DT, modelID, periods) ||
      !setequal(unique(DT$p), Levels) || !setequal(unique(DT$IDn), Models)) {
    stop("kmaxTable has invalid siteID, ID, TR, p, or IDn vocabulary")
  }
  if (anyNA(DT[, .(IDg, IDm, ID, p, IDn, siteID)]) || any(!nzchar(DT$IDg)) || any(!nzchar(DT$IDm)) ||
      any(!is.finite(DT$Da)) || any(DT$Da <= 0) || any(!is.finite(DT$kmax)) || any(DT$kmax < 0) ||
      any(!is.finite(DT$PGA)) || any(DT$PGA <= 0) || any(!is.finite(DT$Kh)) || any(DT$Kh < 0) ||
      any(!is.finite(DT$beyond)) || any(DT$beyond < 0) || any(DT$beyond > 1) ||
      any(!is.finite(DT$Ts)) || any(DT$Ts <= 0) || any(!is.finite(DT$Vs30)) || any(DT$Vs30 <= 0) ||
      any(DT$units.Da != "cm") || any(DT$units.kmax != "g") || any(DT$units.PGA != "g") ||
      any(DT$units.Kh != "ratio")) {
    stop("kmaxTable has invalid identifiers, numerical values, or units")
  }
  if (anyDuplicated(DT, by = KMAXKEY)) stop("kmaxTable key is duplicated for siteID=", siteID)
  if (!.isSameKmaxNumber(sort(unique(DT$Da)), sort(Da))) stop("kmaxTable Da of siteID=", siteID, " does not reproduce Da")
  Rows <- (length(periods) + 1L) * nrow(shear) * length(Models) * length(Levels) * length(Da)
  if (nrow(DT) != Rows) stop("kmaxTable row count mismatch for siteID=", siteID, ": expected ", Rows, ", found ", nrow(DT))
  Groups <- DT[, .(n = .N, np = uniqueN(p)), by = .(ID, TR, IDg, IDm, IDn, Da)]
  if (any(Groups$n != length(Levels)) || any(Groups$np != length(Levels))) {
    stop("kmaxTable does not have the complete p family per demand cell")
  }
  if (!.reproducesKmax(DT, shear, by = c("IDg", "IDm"), cols = c(Ts = "Ts", Vs30 = "VS30"))) {
    stop("kmaxTable Ts/Vs30 does not reproduce ShearTable")
  }
  if (!is.null(PGA) && !.reproducesKmax(DT, PGA, by = c("ID", "TR", "IDg", "IDm"), cols = c(PGA = "PGA"))) {
    stop("kmaxTable PGA does not reproduce PGATable")
  }
  if (!.isSameKmaxNumber(DT$Kh, DT$kmax / DT$PGA)) stop("kmaxTable Kh must equal kmax/PGA as a dimensionless ratio")
  setkeyv(DT, KMAXKEY)
  invisible(DT)
}

.warnKmaxBeyond <- function(x, label) {
  Flagged <- x[beyond > 0]
  if (!nrow(Flagged)) return(invisible(NULL))
  Worst <- Flagged[order(-kmax)][1L]
  warning(label, ": ", nrow(Flagged), "/", nrow(x),
          " kmax rows use roots beyond sampled support (max beyond=", format(max(Flagged$beyond), digits = 4),
          "; max kmax=", format(Worst$kmax, digits = 4), " g at siteID=", Worst$siteID,
          " ID=", Worst$ID, " IDg=", Worst$IDg, " IDm=", Worst$IDm,
          " IDn=", Worst$IDn, " p=", Worst$p, " Da=", Worst$Da, ")", call. = FALSE)
  invisible(NULL)
}

# The rule of the heights from which buildKmaxTable() derives the levels of
# DaH, as the message it stops with, or NULL when the shear table keeps it.
# processProject() applies it to the shear table that the kmaxTable step
# reads before the site steps calculate anything.
.diagnoseKmaxHeights <- function(shear) {
  if (!is.data.frame(shear) || !is.numeric(shear$Hs) || any(!is.finite(shear$Hs)) ||
      any(shear$Hs <= 0)) {
    return("shear must carry finite positive Hs heights: the levels of DaH are derived from them")
  }
  NULL
}

# The rules of DaH and of the heights from which the kmax family derives its
# levels, as the message the caller stops with, or NULL when both hold.
# buildKmaxTable() and assembleKmaxTables() apply them before they derive
# the levels with .deriveKmaxLevels().
.diagnoseKmaxLevels <- function(DaH, shear) {
  if (!is.numeric(DaH) || !length(DaH) || any(!is.finite(DaH)) || any(DaH <= 0) || anyDuplicated(DaH)) {
    return("DaH must contain distinct positive percentages of Hs")
  }
  .diagnoseKmaxHeights(shear)
}

# The levels of the kmax family in cm: every percentage of DaH at every
# height of the shear table, rounded to six significant digits as newmark
# 1.7.1 rounded them, so a geometry's own levels are the rows whose Da is
# signif(DaH * Hs, 6) for its height.
.deriveKmaxLevels <- function(DaH, shear) {
  sort(unique(signif(as.vector(outer(DaH, unique(shear$Hs))), 6)))
}

# Names resolved by data.table inside the kmax family.
utils::globalVariables(c(
  ".", "ky", "s", "Dn", "IDn", "PGA.s", "kmax", "beyond", "p", "ID", "TR", "IDg", "IDm", "Da",
  "PGA", "Ts", "Vs30", "VS30", "level", "i.PGA", "i.Ts", "i.VS30"
))

#' @name buildKmaxTable
#' @examples
#' # The inputs of the product test of newmark 1.7.1: one rock lognormal per
#' # period, its MCE quantiles and two geometry-material pairs,
#' # here with the heights Hs from which the levels of DaH are derived.
#' Rock <- data.table::data.table(ID = "TEST", siteID = "TST", TR = 100, Tn = c(0, 0.2, 1),
#'                                Sa = c(0.3, 0.45, 0.24), Vs30 = 760, Vref = 760)
#' Rock[, `:=`(SaF = Sa, AF = NA_real_, k = "oq", muLn = log(Sa) - 0.4^2 / 2, sigmaLn = 0.4,
#'             w = 1, fitErr = 0)]
#' Levels <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
#' MCE <- Rock[, .(p = c(sprintf("%.2f", Levels), "mean"),
#'                 Sa = c(exp(muLn + 0.4 * qnorm(Levels)), Sa)), by = Tn]
#' MCE[, `:=`(ID = "MCE", siteID = "TST", Vs30 = 760, Vref = 760)]
#' Shear <- data.table::data.table(IDg = c("G1", "G2"), IDm = c("U1", "U2"), Hs = c(10, 120),
#'                                 VS30 = c(360, 560), Ts = 0.3, Vref = 760, level = "mean")
#' Weights <- c(AM88 = 1, BM19 = 1)
#' Site <- buildDnTables(Rock, MCE, Shear, modelID = "TEST", siteID = "TST", TR = 100, Vref = 760,
#'                       Mw = 7, NS = 2, weights = Weights, zTrunc = 3, workers = 1)
#' PGATable <- buildPGATable(Rock, MCE, Shear, modelID = "TEST", siteID = "TST", TR = 100)
#' Kmax <- buildKmaxTable(Site$DnDraws, PGATable = PGATable, shear = Shear, modelID = "TEST",
#'                        siteID = "TST", TR = 100, NS = 2, weights = Weights, DaH = c(0.05, 0.1))
#' Kmax[IDn == "ensemble" & p == "0.50"]
NULL
