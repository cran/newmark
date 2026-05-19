###############################################################################
# fitDnCurve.R — Full Dn(ky) curve for one scenario
###############################################################################

#' Newmark displacement curve Dn(ky) for one scenario
#'
#' Evaluates ensemble and per-submodel Dn over a vector of ky values.
#' The scenario is sampled once and the same epsilon is reused across all ky
#' so that each realisation s represents a coherent physical curve Dn_s(ky).
#'
#' @param uhs      data.table with columns Tn, p, Sa (one scenario)
#' @param ky       numeric vector, yield accelerations (g)
#' @param Ts       numeric scalar, fundamental period (s)
#' @param Mw       numeric scalar, moment magnitude
#' @param NS       integer, Monte Carlo samples for curve quantiles and
#'   per-realisation draws.
#' @param Rrup     numeric scalar, rupture distance (km)
#' @param weights  named numeric vector, ensemble weights by IDn
#' @param NFC      character, near-fault component selector for model BM19:
#'   \code{"D100"} (maximum component, default) or \code{"D50"}.
#'   Ignored by all other submodels.
#'
#' @return list(curve, draws). curve = data.table(ky, p, Dn, IDn, w).
#'   draws = data.table(ky, s, Dn, IDn) — always populated (possibly
#'   empty when no model produces output). Dn in cm.
#' @examples
#' \dontrun{
#' uhs <- data.table::fread(
#'   system.file("extdata", "uhs.csv", package = "newmark"))
#' out <- fitDnCurve(uhs, ky = getDnKy(uhs, Ts = 0.3), Ts = 0.3, Mw = 7.5,
#'                   NS = 100, weights = c(AM88=1, BT07=1, SR08=1, BM17=1))
#' }
#' @export
fitDnCurve <- function(uhs, ky, Ts, Mw = 6.5, NS = 30L,
                       Rrup = 100, weights, NFC = "D100") {

  w <- as.numeric(weights)
  w <- w / sum(w)

  # Sample scenario once — shared across all ky and all models
  .sample <- .sampleDnScenario(uhs, Ts, NS, Rrup)

  # Generate coherent draws for every model across the full ky grid.
  # Each model gets its own pre-generated epsilon so realisations are
  # physical curves, not independent per-ky draws.
  Eps <- stats::rnorm(NS)
  AllDraws <- stats::setNames(
    lapply(names(weights), function(id) {
      .drawDnModelCurve(.sample, ky, Mw, Ts, id, NS, eps = Eps, NFC = NFC)
    }),
    names(weights)
  )

  ok <- names(weights)[vapply(
    AllDraws, function(d) !is.null(d) && nrow(d) > 0L, logical(1)
  )]

  if (!length(ok)) {
    EmptyCurve <- data.table::data.table(
      ky = numeric(), p = character(), Dn = numeric(),
      IDn = character(), w = numeric()
    )
    EmptyDraws <- data.table::data.table(
      ky = numeric(), s = integer(), Dn = numeric(), IDn = character()
    )
    return(list(curve = EmptyCurve, draws = EmptyDraws))
  }

  WOk <- w[match(ok, names(weights))]

  OUT <- vector("list", length(ky))

  for (j in seq_along(ky)) {
    k <- ky[j]

    LIST <- lapply(ok, function(id) {
      d <- AllDraws[[id]][ky == k, Dn]
      if (!length(d)) NULL else d
    })
    OkMask <- !vapply(LIST, is.null, logical(1))
    if (!any(OkMask)) next

    ListJ  <- LIST[OkMask]
    WOkJ   <- WOk[OkMask]
    IdsOkJ <- ok[OkMask]

    # Per-model curve: quantiles + mean
    AUX <- data.table::rbindlist(lapply(seq_along(IdsOkJ), function(i) {
      DT <- .summarizeDn(ListJ[[i]], .sample$pVals)
      DT[, `:=`(IDn = IdsOkJ[i], w = WOkJ[i])]
      DT
    }), use.names = TRUE)

    # Ensemble curve at this ky
    L <- .aggregateWeighted(ListJ, WOkJ, .sample$pVals)
    L[, `:=`(IDn = "ensemble", w = NA_real_)]

    DT <- data.table::rbindlist(list(L, AUX), use.names = TRUE)
    DT[, ky := k]
    OUT[[j]] <- DT
  }

  Curve <- data.table::rbindlist(Filter(Negate(is.null), OUT), use.names = TRUE)
  data.table::setcolorder(Curve, c("ky", "p", "Dn", "IDn", "w"))

  Draws <- data.table::rbindlist(AllDraws[ok], use.names = TRUE)

  list(curve = Curve, draws = Draws)
}
