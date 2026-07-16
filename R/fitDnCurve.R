###############################################################################
# fitDnCurve.R - Full Dn(ky) curve for one scenario
###############################################################################

.fitDnCurvePrepared <- function(dnSample, eps, ky, Ts, Mw, weights, NFC,
                                pulse = FALSE) {
  w <- as.numeric(weights)
  w <- w / sum(w)
  NS <- length(eps)

  AllDraws <- stats::setNames(
    lapply(names(weights), function(id) {
      .drawDnModelCurve(dnSample, ky, Mw, Ts, id, NS, eps = eps, NFC = NFC,
                        pulse = pulse)
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

    AUX <- data.table::rbindlist(lapply(seq_along(IdsOkJ), function(i) {
      DT <- .summarizeDn(ListJ[[i]], dnSample$pVals)
      DT[, `:=`(IDn = IdsOkJ[i], w = WOkJ[i])]
      DT
    }), use.names = TRUE)

    L <- .aggregateWeighted(ListJ, WOkJ, dnSample$pVals)
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
#' @param weights  named numeric vector, ensemble weights by IDn
#' @param NFC      character, near-fault component selector for model BM19:
#'   \code{"D100"} (maximum component, default) or \code{"D50"}.
#'   Ignored by all other submodels.
#' @param pulse    logical, near-fault forward-directivity pulse motion
#'   declared by the analyst (BM19 only; see \code{?Dn_BM19}). Default
#'   \code{FALSE} = ordinary-motion equation.
#'
#' @return list(curve, draws). curve = data.table(ky, p, Dn, IDn, w).
#'   draws = data.table(ky, s, Dn, IDn) - always populated (possibly
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
                       weights, NFC = "D100", pulse = FALSE) {
  DnSample <- .sampleDnScenario(uhs, Ts, NS)
  Eps <- stats::rnorm(NS)
  .fitDnCurvePrepared(
    dnSample = DnSample, eps = Eps, ky = ky, Ts = Ts, Mw = Mw,
    weights = weights, NFC = NFC, pulse = pulse
  )
}
