# SPDX-FileCopyrightText: 2026 Alejandro Verri Kozlowski <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of newmark: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
#' Displacement exceedance mixed over site profiles and yield coefficients
#'
#' Mixes the annual displacement exceedance rates of [integrateDnHazard()]
#' over a discrete distribution of the site profile, given by pairs of the
#' period `Ts` and velocity `VS30` of each profile, and over yield
#' coefficients `ky` with weights, under the declared hypothesis that the
#' profile and the yield coefficient are independent.
#'
#' @section Mathematics:
#' By the total probability theorem, with the normalized weights \eqn{w_i}
#' of the profiles and \eqn{v_k} of the yield coefficients, the mixed
#' annual rate is
#' \deqn{\lambda_D(D_a) = \sum_i \sum_k w_i v_k \lambda_D(D_a \mid T_{s,i}, V_{S30,i}, k_{y,k}),}
#' with each conditional rate the `AEP` of [integrateDnHazard()]. The band
#' rate and the rare-tail bound are mixed with the same weights, and the
#' Poisson probabilities in `IT` years are formed after the mixture, on the
#' mixed rates. Each weight set is normalized to sum one; without a column
#' `w` the profiles weigh equally. For the rigid models the conditional rate
#' does not depend on `Ts`, so the mixture acts on them through `VS30` and
#' `ky` only.
#'
#' When the profiles are realizations of a Monte Carlo population, such as
#' the rows of one pair of the `ShearDraws` product of [processProject()],
#' the mixed rate is an estimate. `mixtureError` is the standard error of
#' the weighted mean over the profiles of \eqn{r_i = \sum_k v_k \lambda_D(D_a \mid T_{s,i}, V_{S30,i}, k_{y,k})},
#' the weighted standard deviation of the \eqn{r_i} (with the divisor
#' \eqn{1 - \sum_i w_i^2} of reliability weights) over the square root of the
#' effective number of profiles \eqn{1 / \sum_i w_i^2}; with equal weights
#' it is the sample standard deviation over \eqn{\sqrt{n}}.
#' `bandMixtureError` is the same quantity for the band rate. Neither is an
#' error bound, and both are `NA` for one profile.
#'
#' @inheritParams integrateDnHazard
#' @param profiles A data.frame or data.table with one row per profile and
#'   the columns `Ts` (s) and `VS30` (m/s), each value as
#'   [integrateDnHazard()] takes it, and optionally `w`, positive finite
#'   weights.
#' @param ky One yield coefficient in g, or a data.frame or data.table with
#'   the columns `ky` (g) and `w`, positive finite weights.
#' @return A list of two data.tables. `summary` has one row: `Da` (cm), `IT`
#'   (years), the mixed annual rate `AEP` (1/year) and its Poisson
#'   probability of exceedance in `IT` years `POE`, the mixed band rate
#'   `bandAEP` (1/year) and its probability `bandPOE`, the mixed rare-tail
#'   bound `rareTailAEPBound` (1/year) and its probability
#'   `rareTailPOEBound`, the counts `nProfiles` and `nKy`, `extrapolated`,
#'   the weighted fraction of the components in which some model has its
#'   `ky` outside its calibration range, and `mixtureError` and
#'   `bandMixtureError` (1/year). `components` has one row per profile and
#'   yield coefficient, profiles first and in their order: `Ts`, `VS30`,
#'   `ky`, `w` (the product of the normalized weights of the profile and of
#'   `ky`), the columns of the `summary` of [integrateDnHazard()] after `ky`
#'   (`Da`, `IT`, `AEP`, `POE`, `bandAEP`, `bandPOE`, `rareTailAEPBound`,
#'   `rareTailPOEBound`) and `extrapolated` (1 when some model of
#'   [integrateDnHazard()] reports its `ky` outside its calibration range,
#'   else 0).
#' @export
#' @examples
#' \donttest{
#' Curves <- data.table::as.data.table(VancouverSpectra)[, .(Tn, Sa, AEP = 1 / TR)]
#' Profiles <- data.frame(Ts = c(0.25, 0.3, 0.4), VS30 = c(320, 300, 270))
#' Ky <- data.frame(ky = c(0.08, 0.1), w = c(1, 3))
#' OUT <- integrateDnHazardProfiles(Curves, profiles = Profiles, ky = Ky, Mw = 7,
#'                                  weights = c(BT07 = 1, BM17 = 1), Da = 15, IT = 50,
#'                                  TR = c(475, 2475), siteModel = "ST20", Vref = 760)
#' OUT$summary
#' OUT$components[, .(Ts, VS30, ky, w, AEP, extrapolated)]
#' }
integrateDnHazardProfiles <- function(curves, profiles, ky, Mw, weights, Da, IT, TR, siteModel, Vref) {
  if (!is.data.frame(profiles) || !nrow(profiles) || !all(c("Ts", "VS30") %in% names(profiles))) {
    stop("profiles must contain rows with Ts and VS30")
  }
  if (!.isNumericColumn(profiles$Ts, profiles) || !.isNumericColumn(profiles$VS30, profiles)) {
    stop("Ts and VS30 of profiles must be numeric, one value per row")
  }
  W <- .normalizeMixtureWeights(if ("w" %in% names(profiles)) profiles$w else rep(1, nrow(profiles)), table = profiles,
                                label = "profiles")
  if (is.numeric(ky) && length(ky) == 1L) ky <- data.frame(ky = ky, w = 1)
  if (!is.data.frame(ky) || !nrow(ky) || !all(c("ky", "w") %in% names(ky)) || !.isNumericColumn(ky$ky, ky)) {
    stop("ky must be one number or a table of numeric ky with weights w, one value per row")
  }
  V <- .normalizeMixtureWeights(ky$w, table = ky, label = "ky")
  # Profiles first, each with every ky in the order given.
  Components <- data.table(Ts = rep(as.numeric(profiles$Ts), each = length(V)),
                           VS30 = rep(as.numeric(profiles$VS30), each = length(V)),
                           ky = rep(as.numeric(ky$ky), length(W)), w = rep(W, each = length(V)) * rep(V, length(W)))
  LIST <- lapply(seq_len(nrow(Components)), function(i) {
    OUT <- integrateDnHazard(curves, ky = Components$ky[[i]], Ts = Components$Ts[[i]], VS30 = Components$VS30[[i]],
                             Vref = Vref, Mw = Mw, weights = weights, Da = Da, IT = IT, TR = TR, siteModel = siteModel)
    OUT$summary[, `:=`(ky = NULL, extrapolated = as.numeric(any(OUT$models$extrapolated == 1)))]
  })
  Components <- cbind(Components, rbindlist(LIST))
  Rate <- sum(Components$w * Components$AEP)
  Band <- sum(Components$w * Components$bandAEP)
  Tail <- sum(Components$w * Components$rareTailAEPBound)
  # The ky-mixed rates of each profile, one column per profile, whose spread
  # over the profiles measures the Monte Carlo error of the mixture.
  Rates <- matrix(Components$AEP, nrow = length(V))
  BandRates <- matrix(Components$bandAEP, nrow = length(V))
  list(summary = data.table(Da = Da, IT = IT, AEP = Rate, POE = -expm1(-IT * Rate),
                            bandAEP = Band, bandPOE = -expm1(-IT * Band),
                            rareTailAEPBound = Tail, rareTailPOEBound = -expm1(-IT * Tail),
                            nProfiles = length(W), nKy = length(V),
                            extrapolated = sum(Components$w * Components$extrapolated),
                            mixtureError = .measureMixtureError(colSums(V * Rates), W),
                            bandMixtureError = .measureMixtureError(colSums(V * BandRates), W)),
       components = Components)
}

# A numeric column with one value per row of its table. is.numeric() alone
# also accepts a matrix column of a data.frame, whose values would spread
# over more components than the table has rows.
.isNumericColumn <- function(x, table) is.numeric(x) && is.null(dim(x)) && length(x) == nrow(table)

# Positive finite weights, one per row of their table, scaled to sum one;
# dividing by the largest first keeps weights near the largest double finite.
.normalizeMixtureWeights <- function(w, table, label) {
  if (!.isNumericColumn(w, table) || any(!is.finite(w) | w <= 0)) {
    stop("the weights w of ", label, " must be positive and finite, one per row")
  }
  w <- w / max(w)
  w / sum(w)
}

# The standard error of the weighted mean of `x` with normalized weights
# `w`: the weighted variance with the divisor of reliability weights,
# 1 - sum(w^2), over the effective count 1 / sum(w^2); NA for one value.
# With weights that sum one the divisor is twice the sum of the products of
# distinct weights, a sum of positive terms that keeps its precision when
# one weight is near one and 1 - sum(w^2) would cancel to zero.
.measureMixtureError <- function(x, w) {
  if (length(w) == 1L) return(NA_real_)
  Divisor <- 2 * sum(w[-1L] * cumsum(w)[-length(w)])
  Concentration <- sum(w^2)
  Mean <- sum(w * x)
  sqrt(sum(w * (x - Mean)^2) / Divisor * Concentration)
}
