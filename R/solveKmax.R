#' Yield coefficient that limits each realized displacement curve
#'
#' For every realization of one displacement model, finds the yield
#' coefficient \eqn{k_{max}} at which the realized Newmark displacement
#' equals each admissible displacement \eqn{D_a}, or the censored value
#' that the rules below assign when no such equality exists on the grid
#' (null demand, grid floor, rigid-block cap). The realizations stay
#' separate and paired with their row of `Dn`; [summarizeKmax()] reduces
#' them, alone or pooled with other models, to a mean and quantiles.
#'
#' @section Mathematics:
#' A realization gives the displacements \eqn{D_j} (cm) at the yield
#' coefficients \eqn{k_{y,1} < \dots < k_{y,J}} (g). Displacements below
#' \eqn{10^{-4}} cm (one micron) are set to zero and the curve is replaced
#' by its upper monotone envelope \eqn{\bar D_j = \max_{i \ge j} D_i},
#' non-increasing in \eqn{k_y}. The coefficient \eqn{k_{max}} solves
#' \eqn{\bar D(k_{max}) = D_a}:
#'
#' - when \eqn{D_a > \bar D_1}, the ceiling of the realization, every
#'   coefficient meets the target and \eqn{k_{max} = 0} (null demand);
#' - between two positive envelope values, \eqn{\log k_{max}} is linear in
#'   \eqn{\log \bar D};
#' - below the smallest positive envelope value of a curve that reaches the
#'   zero floor, \eqn{k_{max}} is the first \eqn{k_y} of that floor (the
#'   resolution of the grid);
#' - below a curve that stays positive on the whole grid, the root is
#'   continued beyond the grid with the closed form of the model, whose free
#'   constants are recovered from the top of the drawn curve. BT07, BM17 and
#'   BM19 are quadratic in \eqn{L = \log k_y}: the parabola
#'   \eqn{\log D = c_1 L^2 + c_2 L + c_3} through the three largest positive
#'   points is solved for its larger root, and a curvature
#'   \eqn{-10^{-8} < c_1 \le 10^{-8}}, the linear limit of the same closed
#'   form, takes the line through the top two points (a curvature above
#'   \eqn{10^{-8}} stops with an error). JB07
#'   is a power law in \eqn{k_y}, inverted along the line through the top
#'   two points. AM88 and SR08 are written in \eqn{r = k_y / PGA}, with
#'   \eqn{PGA} the sampled peak ground acceleration of the realization:
#'   \eqn{\log D = C + 2.53 \log(1 - r) - 1.09 \log r} for AM88 and
#'   \eqn{\log D = C - 5.24 r - 18.78 r^2 + 42.01 r^3 - 29.15 r^4 + 0.56 \varepsilon r}
#'   for SR08, with \eqn{\varepsilon} recovered from the top two points and
#'   \eqn{C} from the top point. The root is sought for \eqn{r \le 0.9999},
#'   the rigid-block limit; a target below the curve at that limit gives
#'   \eqn{k_{max} = 0.9999\, PGA}.
#'
#' An AM88 realization whose displacements are one positive constant lies
#' wholly at the rigid-block limit of the model; its exact inverse solves
#' \eqn{\mu(r) = \log D_a - \log D + \mu(0.9999)}, with
#' \eqn{\mu(r) = 0.90 \log 10 + 2.53 \log(1 - r) - 1.09 \log r}, and
#' \eqn{k_{max} = PGA\, r}. Roots found by continuation or by this inverse
#' lie outside the stored grid and are flagged with `beyond`; null demand,
#' interpolated roots and the grid floor are not.
#'
#' No random number is drawn and nothing is rounded.
#'
#' @section Sources of the closed forms:
#' BT07: Bray and Travasarou (2007), Eq. 5 (p. 387). BM17: Bray, Macedo
#' and Travasarou (2018), Eq. 4 (p. 04017124-5). BM19: Bray and Macedo
#' (2019), Eqs. 3a, 5 and 7, with the corrections of 2023. JB07: Jibson
#' (2007), Eq. 10 (p. 215). SR08: Saygili and Rathje (2008), Eq. 6 and
#' Tables 1-2 (pp. 797-800). AM88: Ambraseys and Menu (1988) as `newmark`
#' 1.7.1 implemented it; the original paper was not verified locally. The
#' full entries are in [evaluateRigidDisplacement()] and
#' [evaluateFlexibleDisplacement()]. The displacement floor of
#' \eqn{10^{-4}} cm, the monotone envelope, the rigid-block cap
#' \eqn{r \le 0.9999} and the first-zero-grid rule are policies of this
#' implementation, not published model terms.
#'
#' @param ky Strictly increasing positive yield coefficients in g, at least
#'   two; a continuation beyond the grid needs three positive displacements
#'   at the top of the curve.
#' @param Dn Matrix of finite non-negative displacements in cm, one row per
#'   realization and one column per `ky`.
#' @param PGA Sampled peak ground acceleration in g of each realization
#'   (one per row of `Dn`), or `NULL` when the realizations carry none; the
#'   continuation of AM88 and SR08 and the inverse of a constant AM88 curve
#'   require it.
#' @param model The displacement model of the realizations: `"AM88"`,
#'   `"JB07"`, `"SR08"`, `"BT07"`, `"BM17"` or `"BM19"`.
#' @param Da Unique positive admissible displacements in cm.
#' @return A data.table with one row per realization and `Da`: `IDn`,
#'   `realization` (the row of `Dn`), `Da` (cm), `kmax` (g), the yield
#'   coefficient at which that realization of the model reaches `Da`, and
#'   `beyond` (`TRUE` when the root comes from the model's continuation
#'   beyond the grid or from the constant-AM88 inverse; `FALSE` for null
#'   demand, interpolated roots and the grid floor).
#' @export
#' @examples
#' ky <- c(0.05, 0.10, 0.20, 0.40)
#' Dn <- matrix(0.02 / ky, nrow = 4, ncol = 4, byrow = TRUE)
#' solveKmax(ky = ky, Dn = Dn, PGA = NULL, model = "BT07", Da = c(0.1, 1))
solveKmax <- function(ky, Dn, PGA, model, Da) {
  Models <- c("AM88", "JB07", "SR08", "BT07", "BM17", "BM19")
  if (!is.character(model) || length(model) != 1L || !model %in% Models) {
    stop("model must be one of ", paste(Models, collapse = ", "))
  }
  if (!is.numeric(ky) || length(ky) < 2L || any(!is.finite(ky)) || any(ky <= 0) || is.unsorted(ky, strictly = TRUE)) {
    stop("ky must hold at least two strictly increasing positive yield coefficients")
  }
  if (!is.matrix(Dn) || !is.numeric(Dn) || !nrow(Dn) || ncol(Dn) != length(ky) || any(!is.finite(Dn)) || any(Dn < 0)) {
    stop("Dn must be a matrix of finite non-negative displacements with one column per ky")
  }
  if (!is.null(PGA) && (!is.numeric(PGA) || length(PGA) != nrow(Dn) || any(!is.finite(PGA)) || any(PGA <= 0))) {
    stop("PGA must be NULL or one finite positive acceleration per realization")
  }
  if (!is.numeric(Da) || !length(Da) || any(!is.finite(Da)) || any(Da <= 0) || anyDuplicated(Da)) {
    stop("Da must contain unique finite positive displacements")
  }
  # Displacements below one micron are numerical underflow of the models,
  # censored before any logarithm.
  FLOOR <- 1e-4
  N <- nrow(Dn)
  kmax <- matrix(0, nrow = length(Da), ncol = N)
  beyond <- matrix(FALSE, nrow = length(Da), ncol = N)
  for (i in seq_len(N)) {
    Curve <- Dn[i, ]
    # A constant positive AM88 curve is resolved on its stored value, before
    # the floor: such curves can sit orders of magnitude below one micron.
    if (model == "AM88" && Curve[[1L]] > 0 && length(unique(Curve)) == 1L) {
      if (is.null(PGA) || ky[[1L]] / PGA[[i]] < RMAX) {
        stop("kmax inversion has non-invertible Dn support for model AM88 in realization ", i)
      }
      kmax[, i] <- .invertFlatAM88(PGA = PGA[[i]], dn = Curve[[1L]], Da = Da)
      beyond[, i] <- TRUE
      next
    }
    Curve[Curve < FLOOR] <- 0
    Root <- .invertDnCurve(ky = ky, dn = Curve, Da = Da, model = model, PGA = if (is.null(PGA)) NULL else PGA[[i]])
    kmax[, i] <- Root$kmax
    beyond[, i] <- Root$beyond
  }
  if (any(!is.finite(kmax)) || any(kmax < 0)) stop("kmax inversion produced invalid roots")
  data.table(IDn = model, realization = rep(seq_len(N), each = length(Da)), Da = rep(Da, N),
             kmax = as.vector(kmax), beyond = as.vector(beyond))
}

# The rigid-block limit of r = ky / PGA: the AM88 model caps r there, and
# the AM88 and SR08 roots are sought below it.
RMAX <- 0.9999

# One floored curve on ascending ky: 0 above its ceiling, the log-log
# interpolation of its monotone envelope, the first floored ky below its
# positive part, and the model continuation below a curve with no floored
# tail.
.invertDnCurve <- function(ky, dn, Da, model, PGA) {
  Envelope <- rev(cummax(rev(dn)))
  kmax <- numeric(length(Da))
  beyond <- logical(length(Da))
  Within <- Da <= Envelope[[1L]]
  if (!any(Within)) return(list(kmax = kmax, beyond = beyond))
  if (length(unique(Envelope)) < 2L) stop("kmax inversion has non-invertible Dn support for model ", model)
  Positive <- Envelope > 0
  Tail <- which(!Positive)
  if (length(unique(Envelope[Positive])) < 2L) {
    kmax[Within] <- ky[Tail[[1L]]]
    return(list(kmax = kmax, beyond = beyond))
  }
  DnLog <- rev(log(Envelope[Positive]))
  KyLog <- rev(log(ky[Positive]))
  DaLog <- log(Da)
  Inside <- Within & DaLog >= DnLog[[1L]] & DaLog <= DnLog[[length(DnLog)]]
  if (any(Inside)) {
    kmax[Inside] <- exp(approx(DnLog, KyLog, xout = DaLog[Inside], rule = 1L, ties = "ordered")$y)
  }
  Below <- Within & DaLog < DnLog[[1L]]
  if (any(Below) && length(Tail)) kmax[Below] <- ky[Tail[[1L]]]
  if (any(Below) && !length(Tail)) {
    kmax[Below] <- .continueDnCurve(ky = ky, dn = dn, target = DaLog[Below], model = model, PGA = PGA)
    beyond[Below] <- TRUE
  }
  list(kmax = kmax, beyond = beyond)
}

# Strong-demand roots beyond the grid from the closed form of each model.
# The constants come from the top of the drawn curve, so the continuation is
# the drawn curve itself rather than a log-log line through its last two
# points (which turned near-flat SR08 plateaus into e17 roots).
.continueDnCurve <- function(ky, dn, target, model, PGA) {
  OK <- dn > 0
  if (sum(OK) < 3L) stop("model continuation needs three positive Dn points for model ", model)
  IDX <- order(ky[OK], decreasing = TRUE)[seq_len(3L)]
  Ky <- ky[OK][IDX]
  Top <- dn[OK][IDX]
  # Along decreasing ky a drawn curve rises strictly; ties or inversions at
  # the top leave the constants without a unique solution.
  if (is.unsorted(Top, strictly = TRUE)) stop("model continuation needs a strictly decreasing curve top for model ", model)
  Y <- log(Top)
  if (model %in% c("BT07", "BM17", "BM19")) {
    L <- log(Ky)
    Coef <- solve(cbind(L^2, L, 1), Y)
    if (!is.finite(Coef[[1L]]) || Coef[[1L]] > 1e-8) stop("quadratic continuation lost its curvature for model ", model)
    # A zero measured curvature is the linear limit of the same closed form;
    # the published curvatures (0.244 to 0.390) lie seven orders above.
    if (Coef[[1L]] > -1e-8) {
      Slope <- (Y[[1L]] - Y[[2L]]) / (L[[1L]] - L[[2L]])
      return(exp(L[[1L]] + (target - Y[[1L]]) / Slope))
    }
    Disc <- Coef[[2L]]^2 - 4 * Coef[[1L]] * (Coef[[3L]] - target)
    if (any(Disc < 0)) stop("quadratic continuation has no real root for model ", model)
    return(exp(pmax((-Coef[[2L]] + sqrt(Disc)) / (2 * Coef[[1L]]), (-Coef[[2L]] - sqrt(Disc)) / (2 * Coef[[1L]]))))
  }
  if (model == "JB07") {
    # The r parameterization cancels in the power law, so no PGA is needed.
    Slope <- (Y[[1L]] - Y[[2L]]) / (log(Ky[[1L]]) - log(Ky[[2L]]))
    return(exp(log(Ky[[1L]]) + (target - Y[[1L]]) / Slope))
  }
  if (!is.numeric(PGA) || length(PGA) != 1L) stop("model continuation requires the sampled PGA for model ", model)
  R <- Ky / PGA
  if (R[[1L]] >= RMAX) stop("model continuation received a saturated curve top for model ", model)
  polySR08 <- function(r) -5.24 * r - 18.78 * r^2 + 42.01 * r^3 - 29.15 * r^4
  shape <- if (model == "AM88") {
    function(r) 2.53 * log1p(-r) - 1.09 * log(r)
  } else {
    Eps <- ((Y[[1L]] - Y[[2L]]) - (polySR08(R[[1L]]) - polySR08(R[[2L]]))) / (0.56 * (R[[1L]] - R[[2L]]))
    function(r) polySR08(r) + 0.56 * Eps * r
  }
  C <- Y[[1L]] - shape(R[[1L]])
  vapply(target, function(t) {
    if (C + shape(RMAX) - t > 0) return(PGA * RMAX)
    PGA * uniroot(function(r) C + shape(r) - t, lower = R[[1L]], upper = RMAX, tol = 1e-12, maxiter = 1000L)$root
  }, numeric(1L))
}

# Exact inverse of an AM88 realization that is one positive constant on the
# whole grid (every ky / PGA at the rigid-block limit).
.invertFlatAM88 <- function(PGA, dn, Da) {
  if (any(Da <= dn)) stop("flat AM88 inversion cannot identify roots for Da <= Dn support")
  Mu <- function(r) 0.90 * log(10) + 2.53 * log1p(-r) - 1.09 * log(r)
  MuMax <- Mu(RMAX)
  vapply(Da, function(target) {
    Target <- log(target) - log(dn) + MuMax
    PGA * uniroot(function(r) Mu(r) - Target, lower = .Machine$double.xmin, upper = RMAX,
                  tol = 1e-12, maxiter = 1000L)$root
  }, numeric(1L))
}
