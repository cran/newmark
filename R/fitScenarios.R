###############################################################################
# fitScenarios.R - Scenario-row helpers for Dn and kmax pipelines
###############################################################################

.scenarioFamilies <- function() {
  list(
    ensemble = c("AM88", "JB07", "SR08", "BT07", "BM17", "BM19"),
    rigid = c("AM88", "JB07", "SR08"),
    flexible = c("BT07", "BM17", "BM19")
  )
}

.scenarioWeights <- function(weights, family, ids) {
  OK <- intersect(intersect(family, names(weights)), ids)
  OK <- OK[weights[OK] > 0]
  OUT <- weights[OK]
  if (!length(OUT)) return(OUT)
  OUT / sum(OUT)
}

.aggregateWeights <- function(weights) {
  OUT <- weights
  OUT["JB07"] <- 0
  OUT
}

.aggregateDnDraws <- function(draws, weights, pVals, idn) {
  draws <- data.table::as.data.table(draws)[IDn %in% names(weights)]
  if (!nrow(draws)) return(data.table::data.table())
  OUT <- draws[
    ,
    {
      w <- weights[IDn]
      AUX <- data.table::data.table()
      if (length(pVals)) {
        q <- Hmisc::wtd.quantile(Dn, weights = w, probs = pVals,
                                 type = "quantile", na.rm = TRUE)
        AUX <- data.table::data.table(p = as.character(pVals),
                                      Dn = as.numeric(q))
      }
      DT <- data.table::data.table(
        p = "mean",
        Dn = Hmisc::wtd.mean(Dn, weights = w, na.rm = TRUE)
      )
      data.table::rbindlist(list(AUX, DT), use.names = TRUE)
    },
    by = ky
  ]
  OUT[, IDn := idn]
  OUT
}

.buildDnScenarioTable <- function(uhs, ky, Ts, Mw, NS, weights,
                                  dnSample = NULL, eps = NULL) {
  Families <- .scenarioFamilies()
  PVals <- sort(unique(as.numeric(uhs[p != "mean", p])))
  PVals <- PVals[is.finite(PVals)]
  AUX <- if (is.null(dnSample) && is.null(eps)) {
    fitDnCurve(
      uhs = uhs, ky = ky, Ts = Ts, Mw = Mw,
      NS = NS, weights = weights
    )
  } else {
    .fitDnCurvePrepared(
      dnSample = dnSample, eps = eps, ky = ky, Ts = Ts, Mw = Mw,
      weights = weights, NFC = "D100", pulse = FALSE
    )
  }
  if (is.null(AUX$curve) || !nrow(AUX$curve)) {
    return(list(curve = data.table::data.table(),
                draws = data.table::data.table()))
  }
  Draws <- AUX$draws
  AggregateWeights <- .aggregateWeights(weights)
  LIST <- lapply(names(Families), function(s) {
    AUX <- .scenarioWeights(weights = AggregateWeights,
                            family = Families[[s]],
                            ids = unique(Draws$IDn))
    if (!length(AUX)) return(NULL)
    .aggregateDnDraws(Draws, AUX, PVals, s)
  })
  names(LIST) <- names(Families)
  LIST <- Filter(Negate(is.null), LIST)
  OUT <- data.table::rbindlist(
    c(list(AUX$curve[IDn %in% names(weights), .(ky, p, Dn, IDn)]), LIST),
    use.names = TRUE
  )
  list(curve = OUT, draws = Draws)
}

.prepareDnScenario <- function(row) {
  row <- data.table::as.data.table(row)
  if (nrow(row) != 1L) stop("Dn scenario row must have one row")
  list(
    row = row,
    dnSample = .sampleDnScenario(
      row$UHS[[1L]], row$Ts[[1L]], row$NS[[1L]]
    ),
    eps = stats::rnorm(row$NS[[1L]])
  )
}

.fitDnScenarioPrepared <- function(x, weights) {
  row <- x$row
  AUX <- .buildDnScenarioTable(
    uhs = row$UHS[[1L]],
    ky = row$ky[[1L]],
    Ts = row$Ts[[1L]],
    Mw = row$Mw[[1L]],
    NS = row$NS[[1L]],
    weights = weights,
    dnSample = x$dnSample,
    eps = x$eps
  )
  DT <- AUX$curve
  if (nrow(DT)) {
    DT[, `:=`(TR = row$TR[[1L]], Ts = row$Ts[[1L]],
              Vs30 = row$Vs30[[1L]], IDg = row$IDg[[1L]],
              IDm = row$IDm[[1L]], ID = row$ID[[1L]])]
  }
  OUT <- AUX$draws
  if (nrow(OUT)) {
    OUT[, `:=`(TR = row$TR[[1L]], IDg = row$IDg[[1L]],
               IDm = row$IDm[[1L]])]
  }
  list(curve = DT, draws = OUT)
}

#' Fit Newmark displacement for one scenario row
#'
#' @param row one-row `data.table` containing `UHS`, `ky`, `Ts`, `Mw`, `NS`,
#'   `TR`, `IDg`, `IDm`, `Vs30`, and `ID`. A legacy `Rrup` column is
#'   tolerated and ignored (the correlation model has no distance
#'   dependence).
#' @param weights named numeric vector of active Newmark model weights by `IDn`.
#'
#' @return list with `curve` and `draws` data.tables. `curve` includes
#'   individual model rows and aggregate `IDn` rows: `ensemble`, `rigid`, and
#'   `flexible`. `draws` contains raw per-realisation model draws.
#' @export
fitDnScenario <- function(row, weights) {
  .fitDnScenarioPrepared(.prepareDnScenario(row), weights = weights)
}

#' Fit Newmark displacement for multiple scenario rows
#'
#' @param records `data.table` of scenario rows. Each row has the contract
#'   required by `fitDnScenario()`.
#' @param weights named numeric vector of active Newmark model weights by `IDn`.
#' @param useParallel logical. Use fork parallelism when available.
#' @param workers integer worker count when `useParallel = TRUE`.
#'
#' @return list with row-bound `curve` and `draws` data.tables.
#' @export
fitDnScenarios <- function(records, weights, useParallel = FALSE,
                           workers = 1L) {
  records <- data.table::as.data.table(records)
  n <- nrow(records)
  if (!n) return(list(curve = data.table::data.table(),
                      draws = data.table::data.table()))
  LIST <- lapply(seq_len(n), function(i) .prepareDnScenario(records[i]))
  processOne <- function(i) .fitDnScenarioPrepared(LIST[[i]],
                                                   weights = weights)
  if (isTRUE(useParallel) && workers > 1L && n > 1L && .forkOK()) {
    Threads <- data.table::getDTthreads()
    data.table::setDTthreads(1L)
    on.exit(data.table::setDTthreads(Threads), add = TRUE)
    OUT <- parallel::mclapply(seq_len(n), processOne, mc.cores = workers,
                              mc.preschedule = FALSE)
  } else {
    OUT <- lapply(seq_len(n), processOne)
  }
  list(
    curve = data.table::rbindlist(lapply(OUT, `[[`, "curve"),
                                  use.names = TRUE),
    draws = data.table::rbindlist(lapply(OUT, `[[`, "draws"),
                                  use.names = TRUE)
  )
}

.invertCurve <- function(ky, dn, Da) {
  if (length(ky) < 2L) stop("Need at least two ky points for kmax inversion")
  IDX <- order(ky)
  x <- ky[IDX]
  y <- rev(cummax(rev(dn[IDX])))
  DnLog <- rev(log(y + .Machine$double.eps))
  KyLog <- rev(log(x))
  DaLog <- log(Da)
  n <- length(DnLog)
  OUT <- numeric(length(Da))
  OK <- DaLog >= DnLog[1L] & DaLog <= DnLog[n]
  if (any(OK)) {
    OUT[OK] <- exp(stats::approx(DnLog, KyLog, xout = DaLog[OK],
                                 rule = 1L, ties = "ordered")$y)
  }
  OK <- DaLog > DnLog[n]
  if (any(OK)) {
    z <- (KyLog[n] - KyLog[n - 1L]) / (DnLog[n] - DnLog[n - 1L])
    OUT[OK] <- exp(KyLog[n] + z * (DaLog[OK] - DnLog[n]))
  }
  OK <- DaLog < DnLog[1L]
  if (any(OK)) {
    z <- (KyLog[2L] - KyLog[1L]) / (DnLog[2L] - DnLog[1L])
    OUT[OK] <- exp(KyLog[1L] + z * (DaLog[OK] - DnLog[1L]))
  }
  OUT
}

.invertDnDrawsFast <- function(draws, Da, weights,
                               p = c(0.05, 0.10, 0.16, 0.84, 0.90, 0.95)) {
  stopifnot(is.numeric(p), all(p > 0 & p < 1))
  Models <- names(weights)
  draws <- data.table::as.data.table(draws)
  weights <- as.numeric(weights)
  weights <- weights / sum(weights)
  LIST <- lapply(Models, function(model) {
    DT <- draws[
      IDn == model,
      .(K = list(.invertCurve(ky, dn = Dn, Da = Da))),
      by = s
    ]
    if (!nrow(DT)) return(NULL)
    DT[, IDn := model]
    DT
  })
  DT <- data.table::rbindlist(Filter(Negate(is.null), LIST), use.names = TRUE)
  OUT <- DT[
    ,
    {
      IDX <- match(Models, IDn)
      OK <- !is.na(IDX)
      AUX <- do.call(rbind, K[IDX[OK]])
      data.table::data.table(
        Da = Da,
        kmax = vapply(seq_along(Da), function(i) {
          stats::weighted.mean(AUX[, i], weights[OK], na.rm = TRUE)
        }, numeric(1))
      )
    },
    by = s
  ]
  OUT[
    ,
    .(p = c("mean", sprintf("%.2f", p)),
      kmax = c(mean(kmax, na.rm = TRUE),
               stats::quantile(kmax, p, na.rm = TRUE, names = FALSE))),
    by = Da
  ]
}

.buildKmaxScenarioTable <- function(draws, weights, Da) {
  draws <- data.table::as.data.table(draws)
  Families <- .scenarioFamilies()
  AggregateWeights <- .aggregateWeights(weights)
  LIST <- lapply(names(weights), function(model) {
    Weight <- stats::setNames(1, model)
    DT <- draws[IDn == model]
    if (!nrow(DT)) return(NULL)
    AUX <- .invertDnDrawsFast(
      draws = DT[, .(ky, s, Dn, IDn)],
      Da = Da,
      weights = Weight
    )
    AUX[, IDn := model]
    AUX
  })
  AUX <- lapply(names(Families), function(scenario) {
    Weight <- .scenarioWeights(weights = AggregateWeights,
                               family = Families[[scenario]],
                               ids = unique(draws$IDn))
    if (!length(Weight)) return(NULL)
    DT <- draws[IDn %in% names(Weight)]
    if (!nrow(DT)) return(NULL)
    AUX <- .invertDnDrawsFast(
      draws = DT[, .(ky, s, Dn, IDn)],
      Da = Da,
      weights = Weight
    )
    AUX[, IDn := scenario]
    AUX
  })
  data.table::rbindlist(
    c(Filter(Negate(is.null), LIST), Filter(Negate(is.null), AUX)),
    use.names = TRUE
  )
}

#' Fit kmax for one scenario row
#'
#' @param row one-row `data.table` containing `TR`, `IDg`, and `PGA_rock`.
#' @param draws `data.table` with `TR`, `IDg`, `ky`, `s`, `Dn`, and `IDn`.
#' @param weights named numeric vector of active Newmark model weights by `IDn`.
#' @param Da numeric displacement targets in cm.
#'
#' @return `data.table` with `Da`, `p`, `kmax`, `PGA`, `Kh`, `TR`, `IDg`,
#'   and `IDn`.
#' @export
fitKmaxScenario <- function(row, draws, weights, Da) {
  row <- data.table::as.data.table(row)
  if (nrow(row) != 1L) stop("kmax scenario row must have one row")
  DT <- draws[TR == row$TR[[1L]] & IDg == row$IDg[[1L]],
              .(ky, s, Dn, IDn)]
  if (!nrow(DT)) stop("Missing DnDraws for TR=", row$TR[[1L]],
                      " IDg=", row$IDg[[1L]])
  OUT <- .buildKmaxScenarioTable(draws = DT, weights = weights, Da = Da)
  OUT[, `:=`(PGA = row$PGA_rock[[1L]], Kh = 100 * kmax / row$PGA_rock[[1L]],
             TR = row$TR[[1L]], IDg = row$IDg[[1L]])]
  OUT
}

.forkOK <- function() {
  !identical(.Platform$OS.type, "windows")
}

#' Fit kmax for multiple scenario rows
#'
#' @param records `data.table` of kmax scenario rows. Each row has the contract
#'   required by `fitKmaxScenario()`.
#' @param draws `data.table` with `TR`, `IDg`, `ky`, `s`, `Dn`, and `IDn`.
#' @param weights named numeric vector of active Newmark model weights by `IDn`.
#' @param Da numeric displacement targets in cm.
#' @param useParallel logical. Use fork parallelism when available.
#' @param workers integer worker count when `useParallel = TRUE`.
#'
#' @return row-bound `data.table` of kmax results.
#' @export
fitKmaxScenarios <- function(records, draws, weights, Da,
                             useParallel = FALSE, workers = 1L) {
  records <- data.table::as.data.table(records)
  n <- nrow(records)
  if (!n) return(data.table::data.table())
  processOne <- function(i) {
    fitKmaxScenario(row = records[i], draws = draws, weights = weights,
                    Da = Da)
  }
  if (isTRUE(useParallel) && workers > 1L && n > 1L && .forkOK()) {
    Threads <- data.table::getDTthreads()
    data.table::setDTthreads(1L)
    on.exit(data.table::setDTthreads(Threads), add = TRUE)
    LIST <- parallel::mclapply(seq_len(n), processOne, mc.cores = workers,
                               mc.preschedule = FALSE)
  } else {
    LIST <- lapply(seq_len(n), processOne)
  }
  data.table::rbindlist(LIST, use.names = TRUE)
}
