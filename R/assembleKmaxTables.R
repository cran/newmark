#' Global yield-coefficient table of release 1.7.1
#'
#' Compatibility route that reproduces, byte for byte once saved with
#' [saveRDS()], the global `kmaxTable` product that the kmax step of `newmark`
#' 1.7.1 assembles with its `getKmaxTable()` from the site `kmaxTable`
#' products of [buildKmaxTable()], so that the stored project tables can
#' serve as comparison oracle. With one site, as in the current projects,
#' the global table equals the table of that site.
#'
#' @section Semantics reproduced from `newmark` 1.7.1:
#' Every site table is checked as [buildKmaxTable()] checks its own, except
#' against its PGA table: the global assembly does not read the per-site PGA
#' tables, which [buildKmaxTable()] validates, with each site table against
#' its own, when it builds that table. A site table is checked against the
#' mean rows of the shear table, the demands of `modelID` at `TR` and of the
#' MCE, the models of `weights` and every family with an active member, the
#' eight levels of `p` and the levels \eqn{D_a = D_{aH} H_s} cm rounded to
#' six significant digits over every height of the shear table. The tables are stacked in
#' the order of `tables` and keyed by `siteID`, `ID`, `TR`, `IDg`, `IDm`,
#' `IDn`, `p` and `Da`, and a warning labelled `global` reports the rows
#' that rely on roots beyond the stored grid and the largest such `kmax`.
#' The number of realizations is not an argument, and of `weights` only the
#' names are used: the table records neither the number of realizations nor
#' the numerical weights that generated it, so a stored table cannot be
#' checked against them.
#'
#' @param tables Named list of site `kmaxTable` products, one per site
#'   identifier.
#' @param shear,modelID,weights,DaH As for [buildKmaxTable()].
#' @param TR Distinct positive return periods of the model demands, in
#'   years. A repeated period stops the call at the row count of the site
#'   tables, whereas [buildKmaxTable()] drops the repeat.
#' @return A data.table of the rows of every site, keyed by `siteID`, `ID`,
#'   `TR`, `IDg`, `IDm`, `IDn`, `p` and `Da`, with the columns of
#'   [buildKmaxTable()]. The caller's tables are not modified.
#' @export
assembleKmaxTables <- function(tables, shear, modelID, TR, weights, DaH) {
  IDs <- names(tables)
  if (!is.list(tables) || !length(tables) || is.null(IDs) || anyDuplicated(IDs) ||
      !all(vapply(IDs, .isKmaxIdentifier, logical(1L)))) {
    stop("tables must be a named list with unique site identifiers")
  }
  Problem <- .diagnoseKmaxLevels(DaH, shear)
  if (!is.null(Problem)) stop(Problem)
  Da <- .deriveKmaxLevels(DaH, shear)
  Shear <- .selectKmaxShear(shear)
  # Each site table is checked on a copy, because the check keys it by
  # reference.
  DT <- rbindlist(lapply(IDs, function(siteID) {
    .checkKmaxTable(copy(as.data.table(tables[[siteID]])), shear = Shear, PGA = NULL, siteID = siteID,
                    modelID = modelID, periods = TR, weights = weights, Da = Da)
  }))
  setkeyv(DT, KMAXKEY)
  .warnKmaxBeyond(DT, label = "global")
  DT
}

#' @name assembleKmaxTables
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
#' Global <- assembleKmaxTables(list(TST = Kmax), shear = Shear, modelID = "TEST", TR = 100,
#'                              weights = Weights, DaH = c(0.05, 0.1))
#' # With one site, the global table equals the table of that site.
#' identical(Global, Kmax)
NULL
