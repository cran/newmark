###############################################################################
# fitDnModel.R — Single-submodel Newmark displacement + internal helpers
###############################################################################

#' Newmark displacement for a single submodel and scenario
#'
#' @param uhs   data.table with columns Tn, p, Sa (one scenario)
#' @param ky    numeric scalar, yield acceleration (g)
#' @param Ts    numeric scalar, fundamental period (s)
#' @param Mw    numeric scalar, moment magnitude
#' @param NS    integer, Monte Carlo samples
#' @param Rrup  numeric scalar, rupture distance (km)
#' @param IDn   character, submodel identifier
#' @param NFC   character, near-fault component selector for model BM19:
#'   \code{"D100"} (maximum component, default) or \code{"D50"}.
#'   Ignored by all other submodels.
#'
#' @return data.table with columns p, Dn, IDn. Dn in cm.
#' @examples
#' \dontrun{
#' uhs <- data.table::fread(
#'   system.file("extdata", "uhs.csv", package = "newmark"))
#' fitDnModel(uhs, ky = 0.1, Ts = 0.3, Mw = 7.5, NS = 100, IDn = "AM88")
#' }
#' @export
fitDnModel <- function(uhs, ky, Ts, Mw = 6.5, NS = 30L, Rrup = 100,
                       IDn, NFC = "D100") {

  .sample <- .sampleDnScenario(uhs, Ts, NS, Rrup)
  Dn      <- .drawDnModel(.sample, ky, Ts, Mw, IDn, NFC = NFC)
  if (is.null(Dn)) {
    return(data.table::data.table(p = character(), Dn = numeric(), IDn = character()))
  }
  DT <- .summarizeDn(Dn, .sample$pVals)
  DT[, IDn := IDn]
  DT
}


# --- Internal helpers (not exported, shared with fitDnCurve) ----------------

.sampleDnScenario <- function(uhs, Ts, NS, Rrup) {
  uhs <- data.table::as.data.table(uhs)[, .(Tn, p, Sa)]
  uhs[, p := as.character(p)]

  periods <- c(0, 1.3 * Ts, 1.5 * Ts)

  LIST <- lapply(periods, function(Ti) {
    if (Ti %in% uhs$Tn) uhs[Tn == Ti]
    else interpolateSaTable(uhs, Tn = Ti)
  })
  SaTable <- data.table::rbindlist(LIST, use.names = TRUE)

  J <- length(periods)
  C <- diag(J)
  for (j in seq_len(J)) for (k in seq_len(J)) {
    if (j != k) {
      Tj <- if (periods[j] == 0) 0.01 else periods[j]
      Tk <- if (periods[k] == 0) 0.01 else periods[k]
      C[j, k] <- rhoBJ(Tj, T0 = Tk)
    }
  }

  SaMatrix <- .sampleSaCorrFull(SaTable, Tn = periods, rho = C, NS = NS)

  PGA  <- SaMatrix[, 1L]
  Sa13 <- SaMatrix[, 2L]
  Sa15 <- SaMatrix[, 3L]
  PGV  <- (PGA^1.0529) * exp(0.1241) * 100
  AI   <- (PGA^1.9228) * exp(2.6109)

  x <- suppressWarnings(as.numeric(uhs[p != "mean", p]))
  pVals <- sort(unique(x[is.finite(x)]))

  list(PGA = PGA, Sa13 = Sa13, Sa15 = Sa15, PGV = PGV, AI = AI,
       pVals = pVals)
}


.drawDnModel <- function(.sample, ky, Ts, Mw, IDn, eps = NULL, NFC = "D100") {
  spec <- tryCatch({
    switch(IDn,
           AM88 = Dn_AM88(.sample$PGA, ky),
           JB07 = Dn_JB07(.sample$PGA, ky, AI = .sample$AI),
           SR08 = Dn_SR08(.sample$PGA, ky, AI = .sample$AI),
           BT07 = Dn_BT07(ky = ky, Sa = .sample$Sa15, Ts = Ts, Mw = Mw),
           BM17 = Dn_BM17(ky = ky, Sa = .sample$Sa15, Ts = Ts, Mw = Mw),
           BM19 = Dn_BM19(ky = ky, Ts = Ts, Sa = .sample$Sa13,
                          PGA = .sample$PGA, PGV = .sample$PGV, Mw = Mw,
                          NFC = NFC),
           stop("Unsupported IDn: ", IDn))
  }, error = function(e) NULL)

  if (is.null(spec)) return(NULL)
  if (is.null(eps)) eps <- stats::rnorm(length(.sample$PGA))
  exp(spec$muLnD + spec$sdLnD * eps)
}


# Draw coherent Dn(ky) curve for one model — same epsilon across all ky.
# If eps is supplied (shared across models), all models use the same
# standard-normal vector; each still scales by its own sigma_lnD.
.drawDnModelCurve <- function(.sample, ky, Mw, Ts, IDn, NS, eps = NULL,
                              NFC = "D100") {
  if (is.null(eps)) eps <- stats::rnorm(NS)
  LIST <- lapply(seq_along(ky), function(j) {
    Dn <- .drawDnModel(.sample, ky[j], Ts, Mw, IDn, eps = eps, NFC = NFC)
    if (is.null(Dn)) return(NULL)
    data.table::data.table(ky = ky[j], s = seq_len(NS), Dn = Dn, IDn = IDn)
  })
  data.table::rbindlist(Filter(Negate(is.null), LIST), use.names = TRUE)
}


# Summarize a draw table (ky, s, Dn) to quantiles + mean per ky.
.summarizeDnDraws <- function(draws, pVals) {
  draws[, .summarizeDn(Dn, pVals), by = .(ky)]
}


.summarizeDn <- function(Dn, pVals) {
  Q <- data.table::data.table()
  if (length(pVals)) {
    q <- stats::quantile(Dn, probs = pVals, na.rm = TRUE, names = FALSE)
    Q <- data.table::data.table(p = as.character(pVals), Dn = as.numeric(q))
  }
  M <- data.table::data.table(p = "mean", Dn = mean(Dn, na.rm = TRUE))
  data.table::rbindlist(list(Q, M), use.names = TRUE)
}


.aggregateWeighted <- function(LIST, w, pVals) {
  x  <- unlist(LIST)
  wt <- rep(w, each = length(LIST[[1L]]))

  Q <- data.table::data.table()
  if (length(pVals)) {
    q <- Hmisc::wtd.quantile(x, weights = wt,
                             probs = pVals, type = "quantile",
                             na.rm = TRUE)
    Q <- data.table::data.table(p = as.character(pVals), Dn = as.numeric(q))
  }
  M <- data.table::data.table(
    p  = "mean",
    Dn = Hmisc::wtd.mean(x, weights = wt, na.rm = TRUE)
  )
  data.table::rbindlist(list(Q, M), use.names = TRUE)
}
