#' Displacement exceedance integrated over rock hazard curves
#'
#' Integrates the existing displacement regressions and ST20/ST17 site
#' uncertainty over annual rock exceedance rates, at a fixed yield coefficient.
#' Rigid models use the PGA curve. Flexible models use Sa at 1.5 Ts and,
#' for the site factor only, the PGA at the same return period. This last
#' pairing is an equal-hazard approximation, not an event joint distribution.
#'
#' For BT07 and BM17 ([evaluateFlexibleDisplacement()]), the integrand uses
#' the lognormal distribution conditional on nonzero displacement; it is not
#' multiplied by the probability of nonzero displacement, so these rates do
#' not represent the unconditional exceedance hazard of the published mixed
#' models. For JB07 and SR08 ([evaluateRigidDisplacement()]), the Arias
#' intensity is [estimateAI()] of each site-amplified PGA, without additional
#' dispersion of the Arias intensity.
#'
#' Repeated AEP nodes retain the largest Sa. Curves are interpolated linearly
#' in log rate and log acceleration; an intermediate oscillator period uses
#' linear interpolation of log acceleration in `log(Tn + 1e-8)` between its two
#' neighbours at equal annual rate. Integration is restricted to their common
#' observed support and the PGA support. No tail is extrapolated.
#' The integral is split at the native hazard knots before adaptive quadrature.
#' The rare-tail bound assumes exceedance probability one beyond that support.
#' Model rates are combined with normalized weights before the Poisson transform.
#' The calibration range of `ky` of each model comes from [buildKyLimits()] with
#' the largest rock PGA of the model's common support; `extrapolated` is 1 when
#' `ky` lies outside it, and the rate is integrated all the same.
#'
#' The standard-normal residual of the site factor is integrated between -8
#' and 8 with its density, in pieces split at each rate where the
#' site-amplified PGA reaches the AM88 cap `ky / PGA = 0.9999` and at the
#' crossings of the median displacement with `Da` on each monotone side of the
#' regression's maximum, by Gauss-Legendre rules of 32, 64 and 128 points per
#' piece. An order is accepted when, against the
#' previous one, the rate and the band rate each agree within 1e-4 of the
#' value plus 1e-10 and the probabilities at the nodes within 1e-4; the
#' function refuses the site with `site quadrature did not converge` when the
#' last order fails. Each interval between hazard knots is integrated
#' adaptively with a relative tolerance of 1e-6, an absolute tolerance of
#' 1e-12 shared by the intervals and at most 1000 subdivisions.
#' `quadratureError`, `bandQuadratureError` and `probabilityError` report the
#' disagreements of the accepted order, and `integrationError` and
#' `bandIntegrationError` the adaptive quadrature's own estimates; none of
#' them is an error bound.
#'
#' @param curves One site's mean rock hazard: a data.frame or data.table with
#'   finite numeric columns `Tn` (s, non-negative, with a PGA curve at
#'   `Tn = 0`), `Sa` (g, positive) and `AEP` (1/year, positive), each curve
#'   decreasing in `Sa` with at least two distinct rates.
#' @param ky Fixed yield coefficient in g, one finite positive value.
#' @param Ts Initial sliding-mass period in seconds, one finite positive value.
#' @param VS30 One target time-averaged shear-wave velocity in m/s.
#' @param Vref One reference velocity in m/s, 760 or 3000.
#' @param Mw Fixed moment magnitude, one finite positive value.
#' @param weights Positive named dimensionless weights of AM88, JB07, SR08,
#'   BT07 or BM17.
#' @param Da Displacement threshold in cm, one finite positive value.
#' @param IT Exposure in years, one finite positive value.
#' @param TR Two increasing finite positive return periods in years; the band
#'   contribution covers their intersection with each model's observed support
#'   and is zero for an empty intersection.
#' @param siteModel Site-factor model, ST20 or ST17.
#' @return A list of three data.tables. `summary` has one row: `ky` (g), `Da`
#'   (cm), `IT` (years), the weighted annual exceedance rate `AEP` (1/year) and
#'   its Poisson probability of exceedance in `IT` years `POE`, the band
#'   contribution `bandAEP` (1/year) and `bandPOE`, the Poisson probability of
#'   exceedance in `IT` years of the band rate, not a fraction of `POE`, and
#'   the rare-tail upper bounds `rareTailAEPBound` (1/year) and
#'   `rareTailPOEBound`, which bound only the omitted events rarer than the
#'   common support; the omitted frequent side is not bounded. `models` has
#'   one row per displacement model, `IDn`: the normalized `weight`, the model
#'   period `Tn` (s), `ky`, `Da`, the bounds `kyMin` and `kyMax` (g) of the
#'   calibration range of `ky` of the model, `extrapolated` (1 when `ky` lies
#'   outside that range, else 0; the diagnostic at the largest rock PGA of the
#'   common support, not a certification of calibration along the integral),
#'   `AEP`, `POE`, `bandAEP`, `rareTailBound` (1/year), the common support
#'   `TR.min` and `TR.max` (years), `quadratureError` (1/year),
#'   `bandQuadratureError` (1/year) and `probabilityError`, the disagreement
#'   between the last two Gauss-Legendre orders of the rate, of the band rate
#'   and of the node probabilities (the largest over the nodes), each not an
#'   error bound, `integrationError` (1/year), the sum of the adaptive
#'   quadrature's own estimates over the intervals, and `bandIntegrationError`
#'   (1/year), the same sum over the intervals of the band. `conditional` has
#'   one row per model and rate node, `IDn` and `AEP` (1/year): `TR` (years),
#'   the rock `Sa` at the model period and the rock `PGA` (g), and the
#'   `probability` that the displacement exceeds `Da` at that node.
#'   Probabilities are dimensionless.
#' @export
#' @examples
#' Curves <- data.table::as.data.table(VancouverSpectra)[, .(Tn, Sa, AEP = 1 / TR)]
#' OUT <- integrateDnHazard(Curves, ky = 0.1, Ts = 0.3, VS30 = 300, Vref = 760, Mw = 7,
#'                          weights = c(BT07 = 1, BM17 = 1), Da = 15, IT = 50,
#'                          TR = c(475, 2475), siteModel = "ST20")
#' OUT$summary
#' OUT$models[, .(IDn, Tn, AEP, POE, extrapolated)]
integrateDnHazard <- function(curves, ky, Ts, VS30, Vref, Mw, weights, Da, IT, TR, siteModel) {
  .checkSiteCondition(VS30 = VS30, Vref = Vref, model = siteModel)
  for (x in list(ky, Ts, Mw, Da, IT)) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop("ky, Ts, Mw, Da and IT must be positive finite scalars")
    }
  }
  Allowed <- c("AM88", "JB07", "SR08", "BT07", "BM17")
  if (!is.numeric(weights) || !length(weights) || is.null(names(weights)) ||
      anyDuplicated(names(weights)) || any(!names(weights) %in% Allowed) ||
      any(!is.finite(weights) | weights <= 0)) stop("weights must name supported models with positive values")
  if (!is.numeric(TR) || length(TR) != 2L || any(!is.finite(TR) | TR <= 0) || TR[1L] >= TR[2L]) {
    stop("TR must contain two increasing positive return periods")
  }
  if (!is.data.frame(curves) || !all(c("Tn", "Sa", "AEP") %in% names(curves))) {
    stop("curves must contain Tn, Sa and AEP")
  }
  DT <- as.data.table(curves)[, .(Tn, Sa, AEP)]
  if (!all(vapply(DT, is.numeric, logical(1L)))) stop("Tn, Sa and AEP must be numeric")
  if (any(!is.finite(as.matrix(DT))) || any(DT$Tn < 0 | DT$Sa <= 0 | DT$AEP <= 0) ||
      !any(DT$Tn == 0) || anyDuplicated(DT, by = c("Tn", "Sa"))) stop("invalid rock hazard curves")
  Periods <- sort(unique(DT$Tn))
  Curves <- lapply(Periods, function(period) {
    AUX <- DT[Tn == period][order(Sa)]
    if (nrow(AUX) < 2L || any(diff(AUX$AEP) > 0)) stop("each hazard curve must decrease with acceleration")
    AUX <- AUX[, .(Sa = max(Sa), nSa = uniqueN(signif(Sa, 12))), by = AEP][order(AEP)]
    if (any(AUX$nSa > 1L)) {
      .inform("moduleNote", "Collapsing duplicate AEP nodes by keeping max Sa at Tn = ", period, ".", operation = "integrateDnHazard")
    }
    if (nrow(AUX) < 2L) stop("each curve needs distinct rates")
    list(rate = AUX$AEP, at = stats::approxfun(x = log(AUX$AEP), y = log(AUX$Sa), ties = "ordered"))
  })
  names(Curves) <- as.character(Periods)
  Weights <- weights / max(weights)
  Weights <- Weights / sum(Weights)
  Rules <- lapply(c(32L, 64L, 128L), .gaussLegendreRule)
  OUT <- lapply(names(Weights), function(model) {
    IsRigid <- model %in% c("AM88", "JB07", "SR08")
    Period <- if (IsRigid) 0 else 1.5 * Ts
    if (Period > max(Periods)) stop("rock hazard does not cover the model period")
    Bracket <- Period
    if (!Period %in% Periods) Bracket <- c(max(Periods[Periods < Period]), min(Periods[Periods > Period]))
    Used <- Curves[as.character(unique(c(0, Bracket)))]
    Limits <- c(max(vapply(Used, function(x) min(x$rate), numeric(1L))),
                min(vapply(Used, function(x) max(x$rate), numeric(1L))))
    if (Limits[1L] >= Limits[2L]) stop("hazard curves have no common rate support")
    if (length(Bracket) == 2L) Fraction <- (log(Period + 1e-8) - log(Bracket[1L] + 1e-8)) / diff(log(Bracket + 1e-8))
    .interpolateSa <- function(z) {
      if (length(Bracket) == 1L) return(exp(Curves[[as.character(Bracket)]]$at(z)))
      exp((1 - Fraction) * Curves[[as.character(Bracket[1L])]]$at(z) +
            Fraction * Curves[[as.character(Bracket[2L])]]$at(z))
    }
    .evaluateSite <- function(z) {
      evaluateSiteSa(Sa = .interpolateSa(z), PGA = exp(Curves[["0"]]$at(z)), Tn = Period,
                     VS30 = VS30, Vref = Vref, model = siteModel)
    }
    .evaluateModel <- function(IM) {
      if (IsRigid) return(evaluateRigidDisplacement(model = model, PGA = IM, ky = ky, AI = estimateAI(IM)))
      evaluateFlexibleDisplacement(model = model, ky = ky, Ts = Ts, Sa = IM, Mw = Mw,
                                   PGV = rep(1, length(IM)), NFC = "D100", pulse = FALSE)
    }
    Nodes <- sort(unique(c(unlist(lapply(Used, `[[`, "rate")), Limits, 1 / TR)))
    Nodes <- Nodes[Nodes >= Limits[1L] & Nodes <= Limits[2L]]
    # The log site intensities where the integrand over the site residual
    # bends: the AM88 cap of ky / PGA at 0.9999 of evaluateRigidDisplacement(),
    # and the crossings of the median displacement with Da on each monotone
    # side of the regression's maximum, within the enclosure of the intensities
    # the rule evaluates. Over the support, log Sa and log PGA are monotone in
    # the rate, and the site-factor mean `muLnF` and dispersion of ST20 and
    # ST17 are monotone in PGA, so the extremes of log Sa, of the site-factor
    # mean `muLnF` and of the dispersion sit at the end nodes; the residual
    # pieces live in [-8, 8]. A site model added later must also be monotone
    # in PGA.
    .measureMedian <- function(u) .evaluateModel(exp(u))$muLnD - log(Da)
    AUX <- .evaluateSite(log(Nodes))
    LnSa <- log(.interpolateSa(log(Nodes)))
    MuF <- AUX$muLnSaF - LnSa
    Enclosure <- c(min(LnSa) + min(MuF), max(LnSa) + max(MuF)) + c(-8, 8) * max(AUX$sigmaLnSaF)
    Splits <- sort(c(if (model == "AM88") log(ky / 0.9999), .findMedianCrossings(.measureMedian, Enclosure)))
    .evaluateExceedance <- function(z, rule) {
      Site <- .evaluateSite(z)
      Residual <- pmin(pmax(outer(-Site$muLnSaF, Splits, "+") / Site$sigmaLnSaF, -8), 8)
      # A site without dispersion integrates one constant piece.
      Residual[Site$sigmaLnSaF == 0, ] <- -8
      Bounds <- cbind(-8, Residual, 8)
      Lower <- as.vector(Bounds[, -ncol(Bounds)])
      Width <- as.vector(Bounds[, -1L]) - Lower
      e <- Lower + outer(Width, rule$u)
      Model <- .evaluateModel(exp(as.vector(Site$muLnSaF + Site$sigmaLnSaF * e)))
      p <- stats::pnorm((log(Da) - Model$muLnD) / Model$sigmaLnD, lower.tail = FALSE)
      rowSums(matrix(p * stats::dnorm(e) * outer(Width, rule$w), nrow = length(z)))
    }
    .integrateRate <- function(rule, limits) {
      if (limits[1L] >= limits[2L]) return(list(value = 0, abs.error = 0))
      Breaks <- sort(unique(log(c(limits, Nodes[Nodes > limits[1L] & Nodes < limits[2L]]))))
      LIST <- lapply(seq_len(length(Breaks) - 1L), function(i) {
        stats::integrate(function(z) .evaluateExceedance(z, rule) * exp(z),
                         lower = Breaks[i], upper = Breaks[i + 1L],
                         rel.tol = 1e-6, abs.tol = 1e-12 / (length(Breaks) - 1L),
                         subdivisions = 1000L)
      })
      list(value = sum(vapply(LIST, `[[`, numeric(1L), "value")),
           abs.error = sum(vapply(LIST, `[[`, numeric(1L), "abs.error")))
    }
    Band <- c(max(Limits[1L], 1 / TR[2L]), min(Limits[2L], 1 / TR[1L]))
    .integrateOrder <- function(rule) {
      list(rate = .integrateRate(rule, Limits), band = .integrateRate(rule, Band),
           probability = .evaluateExceedance(log(Nodes), rule))
    }
    Previous <- .integrateOrder(Rules[[1L]])
    Converged <- FALSE
    for (Rule in Rules[-1L]) {
      Current <- .integrateOrder(Rule)
      Error <- c(rate = abs(Current$rate$value - Previous$rate$value),
                 band = abs(Current$band$value - Previous$band$value),
                 probability = max(abs(Current$probability - Previous$probability)))
      if (Error[["rate"]] <= 1e-10 + 1e-4 * abs(Current$rate$value) &&
          Error[["band"]] <= 1e-10 + 1e-4 * abs(Current$band$value) && Error[["probability"]] <= 1e-4) {
        Converged <- TRUE
        break
      }
      Previous <- Current
    }
    if (!Converged) stop("site quadrature did not converge for ", model)
    Calibration <- buildKyLimits(PGA = exp(Curves[["0"]]$at(log(Limits[1L]))))[IDn == model]
    list(model = data.table(IDn = model, weight = unname(Weights[[model]]), Tn = Period,
                            ky = ky, Da = Da, kyMin = Calibration$kyMin, kyMax = Calibration$kyMax,
                            extrapolated = as.numeric(ky < Calibration$kyMin | ky > Calibration$kyMax),
                            AEP = Current$rate$value,
                            POE = -expm1(-IT * Current$rate$value),
                            bandAEP = Current$band$value, rareTailBound = Limits[1L],
                            TR.min = 1 / Limits[2L], TR.max = 1 / Limits[1L],
                            quadratureError = Error[["rate"]], bandQuadratureError = Error[["band"]],
                            probabilityError = Error[["probability"]],
                            integrationError = Current$rate$abs.error,
                            bandIntegrationError = Current$band$abs.error),
         conditional = data.table(IDn = model, TR = 1 / Nodes, AEP = Nodes,
                                    Sa = .interpolateSa(log(Nodes)), PGA = exp(Curves[["0"]]$at(log(Nodes))),
                                    probability = Current$probability))
  })
  Models <- rbindlist(lapply(OUT, `[[`, "model"))
  Rate <- sum(Models$weight * Models$AEP)
  Band <- sum(Models$weight * Models$bandAEP)
  Tail <- sum(Models$weight * Models$rareTailBound)
  list(summary = data.table(ky = ky, Da = Da, IT = IT,
                            AEP = Rate, POE = -expm1(-IT * Rate),
                            bandAEP = Band, bandPOE = -expm1(-IT * Band),
                            rareTailAEPBound = Tail, rareTailPOEBound = -expm1(-IT * Tail)),
       models = Models, conditional = rbindlist(lapply(OUT, `[[`, "conditional")))
}

# The intensities u of `enclosure` where `measure(u)`, the median log
# displacement less log(Da), is zero. The flexible regressions are concave in
# log intensity and the rigid ones nondecreasing, so `measure` is monotone on
# each side of its maximum, which stats::optimize() locates. The measured
# height of that maximum, not a tolerance, decides the crossings: a positive
# maximum brackets one root on each side whose end value is negative, which
# stats::uniroot() refines; a maximum exactly zero is the one root; a
# negative maximum gives none. An end value exactly zero is a root, and roots
# within 1e-9 of each other are kept once. A regression added later must be
# unimodal in log intensity.
.findMedianCrossings <- function(measure, enclosure) {
  u <- c(enclosure[1L], stats::optimize(measure, enclosure, maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum,
         enclosure[2L])
  AUX <- measure(u)
  Roots <- sort(c(u[AUX == 0], vapply(which(AUX[-3L] * AUX[-1L] < 0), function(i) {
    stats::uniroot(measure, u[c(i, i + 1L)], tol = .Machine$double.eps)$root
  }, numeric(1L))))
  Roots[diff(c(-Inf, Roots)) > 1e-9]
}
