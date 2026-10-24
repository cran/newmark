#' Displacement realizations of the models of a scenario
#'
#' Evaluates, for every requested displacement model, yield acceleration and
#' realization of a scenario, the Newmark displacement of that realization.
#' One standard normal deviate per realization is shared by every model and
#' every yield acceleration, so each realization traces one coherent
#' displacement curve \eqn{D_{n,s}(k_y)} per model and the models stay
#' paired realization by realization.
#'
#' @section Mathematics:
#' \deqn{D_{n,s}(k_y) = \exp\left(\mu_{\ln D}(k_y, IM_s) + \sigma_{\ln D}(k_y, IM_s)\, \varepsilon_s\right)}
#' in cm, where \eqn{IM_s} are the intensity measures of realization
#' \eqn{s} ([sampleIntensities()]), \eqn{\varepsilon_s} its deviate, and
#' \eqn{\mu_{\ln D}}, \eqn{\sigma_{\ln D}} the conditional lognormal
#' parameters of each model: [evaluateRigidDisplacement()] for `AM88`
#' (PGA), `JB07` and `SR08` (PGA, \eqn{I_a}), and
#' [evaluateFlexibleDisplacement()] for `BT07` and `BM17`
#' (\eqn{S_a(1.5 T_s)}) and `BM19` (\eqn{S_a(1.3 T_s)}, PGV).
#' The calculation is deterministic; the caller supplies the deviates.
#' For `BT07`, `BM17` and `BM19` these are the parameters of the nonzero
#' branch of the published models ([evaluateFlexibleDisplacement()]): their
#' probability of zero or negligible displacement is omitted, so every
#' realization of those models is a positive displacement drawn conditional
#' on a nonzero displacement.
#'
#' @param intensities Paired realizations of one scenario with columns `realization`,
#'   `PGA`, `Sa13`, `Sa15`, `PGV` and `AI`, as [sampleIntensities()]
#'   returns; `realization` identifies each realization and must be non-missing and
#'   unique.
#' @param eps Standard normal deviates of the displacement models, one per
#'   row of `intensities`, finite.
#' @param ky Yield accelerations in g, finite and positive.
#' @param Ts One initial fundamental period of the sliding mass in seconds, the one
#'   the intensities were sampled for.
#' @param Mw One moment magnitude of the scenario, read by the flexible
#'   models.
#' @param models Distinct model identifiers among `AM88`, `JB07`, `SR08`,
#'   `BT07`, `BM17` and `BM19`.
#' @param NFC,pulse Near-fault component and pulse declaration of the
#'   `BM19` model, see [evaluateFlexibleDisplacement()]; read only when
#'   `models` includes `BM19`.
#' @return A data.table with `IDn`, `ky` (g), `realization` and `Dn` (cm), one row per
#'   model, yield acceleration and realization, ordered by the order of
#'   `models`, then `ky`, then the rows of `intensities`. Each `Dn` is one
#'   realization, not a statistic: the population is the realizations of
#'   the scenario, identified by the caller's `realization` and paired across models
#'   and yield accelerations, conditional on a nonzero displacement for the
#'   flexible models; [summarizeDisplacements()] and [summarizeEnsemble()]
#'   reduce it to quantiles and means.
#' @export
#' @examples
#' Sigma <- c(0.2, 0.45, 0.75)
#' Spectrum <- data.frame(Tn = c(0, 0.4, 1), k = "oq", muLn = log(c(0.8, 0.55, 0.42)) - Sigma^2 / 2,
#'                        sigmaLn = Sigma, w = 1)
#' set.seed(46)
#' Intensities <- sampleIntensities(Spectrum, Ts = 0.4, NS = 5, zTrunc = NA)
#' evaluateDisplacements(Intensities, eps = stats::rnorm(5), ky = c(0.05, 0.1), Ts = 0.4, Mw = 7,
#'                       models = c("AM88", "BM19"), NFC = "D100", pulse = FALSE)
evaluateDisplacements <- function(intensities, eps, ky, Ts, Mw, models, NFC, pulse) {
  Required <- c("realization", "PGA", "Sa13", "Sa15", "PGV", "AI")
  if (!is.data.frame(intensities) || !nrow(intensities) || !all(Required %in% names(intensities))) {
    stop("intensities must contain rows with columns ", paste(Required, collapse = ", "))
  }
  if (anyNA(intensities$realization) || anyDuplicated(intensities$realization)) {
    stop("intensities must identify every realization by a non-missing unique realization")
  }
  if (!is.numeric(eps) || length(eps) != nrow(intensities) || any(!is.finite(eps))) {
    stop("eps must hold one finite deviate per realization")
  }
  if (!is.numeric(ky) || !length(ky) || any(!is.finite(ky)) || any(ky <= 0)) {
    stop("ky must be finite positive accelerations")
  }
  if (!is.numeric(Mw) || length(Mw) != 1L || !is.finite(Mw) || Mw <= 0) {
    stop("Mw must be one finite positive magnitude")
  }
  if (!is.character(models) || !length(models) || anyNA(models) || anyDuplicated(models) ||
      !all(models %in% c("AM88", "JB07", "SR08", "BT07", "BM17", "BM19"))) {
    stop("models must name distinct models among AM88, JB07, SR08, BT07, BM17 and BM19")
  }
  # Every realization at every yield acceleration, ky varying slowest.
  IM <- as.data.table(intensities)[rep(seq_len(nrow(intensities)), times = length(ky))]
  Ky <- rep(ky, each = nrow(intensities))
  E <- rep(eps, times = length(ky))
  rbindlist(lapply(models, function(model) {
    Model <- if (model %in% c("AM88", "JB07", "SR08")) {
      evaluateRigidDisplacement(model, PGA = IM$PGA, ky = Ky, AI = IM$AI)
    } else {
      evaluateFlexibleDisplacement(model, ky = Ky, Ts = Ts, Sa = if (model == "BM19") IM$Sa13 else IM$Sa15,
                                   Mw = Mw, PGV = IM$PGV, NFC = NFC, pulse = pulse)
    }
    data.table(IDn = model, ky = Ky, realization = IM$realization, Dn = exp(Model$muLnD + Model$sigmaLnD * E))
  }))
}
