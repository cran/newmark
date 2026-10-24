#' PGA table of release 1.7.1 from rock demands and a shear table
#'
#' Compatibility route that reproduces, value for value, the `PGATable`
#' product that the `dn` step of `newmark` 1.7.1 writes for one site, so that
#' the stored products serve as comparison oracle while the calculations
#' migrate. For every return period of `TR` and for the MCE, and for
#' every geometry-material pair of the shear table, it evaluates the mean
#' site PGA at the pair's VS30 with the ST20 factor and the numerical
#' controls of that product.
#'
#' @section Semantics reproduced from `newmark` 1.7.1:
#' The pairs are the `level == "mean"` rows of `shear`, whose `VS30` is the
#' site velocity; the rock condition is 760 m/s. A return period's rock
#' demand is the analytic mixture of `rock` for `modelID` and `siteID` at
#' that period, and its PGA is the mean
#' \eqn{\sum w e^{\mu + \sigma^2/2}} of the rows at \eqn{|T_n| \le
#' 10^{-10}} of the site mixture of [refineSiteMixture()] with the levels
#' 0.05, 0.10, 0.16, 0.50, 0.84, 0.90 and 0.95, the orders 16, 32, 64, 128
#' and 256 and the tolerance \eqn{2 \times 10^{-3}}. The MCE demand is the
#' quantile spectrum of `MCE` (`ID == "MCE"`) for the site, and its PGA is the
#' `"mean"` site ordinate at \eqn{|T_n| \le 10^{-10}} of
#' [refineSiteSpectrum()] with the orders 32, 64, 128 and 256 and the same
#' tolerance. Pairs with equal VS30 share one evaluation, which is
#' deterministic. The table is keyed by `siteID`, `ID`, `TR`, `IDg` and
#' `IDm`, as the product was written; `data.table::fwrite(x, na = "NA")`
#' gives its CSV file.
#'
#' @param rock The `UHSRock` table of the site: columns `ID`, `siteID`, `TR`
#'   (years), `Tn` (s), `k`, `muLn`, `sigmaLn`, `w`, `Vs30` and `Vref` (m/s).
#' @param MCE The MCE table: columns `ID`, `siteID`, `Tn` (s), `p`, `Sa`
#'   (g), `Vs30` and `Vref` (m/s), without `TR`.
#' @param shear The shear table, see [buildShearTable()]: columns `IDg`,
#'   `IDm`, `VS30` (m/s), `Vref` (760 m/s) and `level`.
#' @param modelID One identifier of the probabilistic demand model (`ID` of
#'   `rock`).
#' @param siteID One site identifier.
#' @param TR Positive return periods in years; the return periods of
#'   `rock` for the model and site must be exactly these.
#' @return A data.table with `IDg`, `IDm`, `ID` (`modelID` or `"MCE"`), `TR`
#'   (years, `NA` for the MCE), `PGA` (g, the mean site PGA, unrounded),
#'   `Vs30` (m/s), `siteID` and `units.PGA` (`"g"`), one row per demand and
#'   pair, keyed as the product.
#' @export
buildPGATable <- function(rock, MCE, shear, modelID, siteID, TR) {
  for (Table in list(list(rock, c("ID", "siteID", "TR", "Tn", "k", "muLn", "sigmaLn", "w", "Vs30", "Vref"), "rock"),
                     list(MCE, c("ID", "siteID", "Tn", "p", "Sa", "Vs30", "Vref"), "MCE"),
                     list(shear, c("IDg", "IDm", "VS30", "Vref", "level"), "shear"))) {
    if (!is.data.frame(Table[[1L]]) || !all(Table[[2L]] %in% names(Table[[1L]]))) {
      stop(Table[[3L]], " must be a table with columns ", paste(Table[[2L]], collapse = ", "))
    }
  }
  if ("TR" %in% names(MCE)) stop("MCE must not contain TR")
  for (Name in list(modelID, siteID)) {
    if (!is.character(Name) || length(Name) != 1L || is.na(Name) || !nzchar(Name)) {
      stop("modelID and siteID must be single nonempty identifiers")
    }
  }
  if (!is.numeric(TR) || !length(TR) || any(!is.finite(TR)) || any(TR <= 0)) {
    stop("TR must contain positive finite return periods")
  }
  # Numerical controls of the newmark 1.7.1 dn product (siteConvergence.R and
  # productTables.R), reproduced as its contract.
  Levels <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
  Tolerance <- 2e-3
  ZeroPeriod <- 1e-10
  Model <- modelID
  Site <- siteID
  Family <- sort(unique(as.numeric(TR)))
  Shear <- as.data.table(shear)[level == "mean", .(IDg = as.character(IDg), IDm = as.character(IDm),
                                                   VS30 = as.numeric(VS30), Vref = as.numeric(Vref))]
  if (!nrow(Shear) || anyNA(Shear) || any(!nzchar(Shear$IDg)) || any(!nzchar(Shear$IDm)) ||
      any(Shear$VS30 <= 0) || any(Shear$Vref != 760) || anyDuplicated(Shear, by = c("IDg", "IDm"))) {
    stop("shear needs unique mean rows by IDg and IDm with positive VS30 and Vref = 760")
  }
  Rock <- as.data.table(rock)[ID == Model & siteID == Site & Vs30 == 760 & Vref == 760,
                              .(TR, Tn, k, muLn, sigmaLn, w)]
  if (!nrow(Rock) || !setequal(unique(Rock$TR), Family)) {
    stop("rock must hold the return periods of TR at Vs30 = Vref = 760 for ID=", Model, ", siteID=", Site)
  }
  Mce <- as.data.table(MCE)[ID == "MCE" & siteID == Site & Vs30 == 760 & Vref == 760,
                            .(Tn = as.numeric(Tn), p = as.character(p), Sa = as.numeric(Sa))]
  if (!nrow(Mce)) stop("MCE has no MCE rows at Vs30 = Vref = 760 for siteID=", Site)
  setorder(Mce, Tn, p)
  LIST <- lapply(sort(unique(Shear$VS30)), function(v) {
    PGA.TR <- vapply(Family, function(tr) {
      Mixture <- refineSiteMixture(Rock[TR == tr, .(Tn, k, muLn, sigmaLn, w)], VS30 = v, Vref = 760,
                                   model = "ST20", levels = Levels, orders = c(16L, 32L, 64L, 128L, 256L),
                                   tolerance = Tolerance)$mixture
      Mixture[abs(Tn) <= ZeroPeriod, sum(w * exp(muLn + sigmaLn^2 / 2))]
    }, numeric(1L))
    Spectrum <- refineSiteSpectrum(Mce, VS30 = v, Vref = 760, model = "ST20",
                                   orders = c(32L, 64L, 128L, 256L), tolerance = Tolerance)$spectrum
    PGA.MCE <- Spectrum[abs(Tn) <= ZeroPeriod & p == "mean", SaF]
    data.table(VS30 = v, ID = c(rep(Model, length(Family)), "MCE"), TR = c(Family, NA_real_),
               PGA = c(PGA.TR, PGA.MCE))
  })
  OUT <- Shear[, .(IDg, IDm, VS30)][rbindlist(LIST), on = "VS30", allow.cartesian = TRUE]
  OUT <- OUT[, .(IDg, IDm, ID, TR, PGA, Vs30 = VS30, siteID = Site, units.PGA = "g")]
  setkeyv(OUT, c("siteID", "ID", "TR", "IDg", "IDm"))
  OUT[]
}

#' @name buildPGATable
#' @examples
#' # The inputs of the product test of newmark 1.7.1: one rock lognormal per
#' # period, its MCE quantiles and two geometry-material pairs.
#' Rock <- data.table::data.table(ID = "TEST", siteID = "TST", TR = 100, Tn = c(0, 0.2, 1),
#'                                Sa = c(0.3, 0.45, 0.24), Vs30 = 760, Vref = 760)
#' Rock[, `:=`(SaF = Sa, AF = NA_real_, k = "oq", muLn = log(Sa) - 0.4^2 / 2, sigmaLn = 0.4,
#'             w = 1, fitErr = 0)]
#' Levels <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
#' MCE <- Rock[, .(p = c(sprintf("%.2f", Levels), "mean"),
#'                 Sa = c(exp(muLn + 0.4 * qnorm(Levels)), Sa)), by = Tn]
#' MCE[, `:=`(ID = "MCE", siteID = "TST", Vs30 = 760, Vref = 760)]
#' Shear <- data.table::data.table(IDg = c("G1", "G2"), IDm = c("U1", "U2"), VS30 = c(360, 560),
#'                                 Ts = 0.3, Vref = 760, level = "mean")
#' buildPGATable(Rock, MCE, Shear, modelID = "TEST", siteID = "TST", TR = 100)
NULL
