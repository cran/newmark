#' Build the characteristic-root model of one mode
#'
#' The roots \eqn{a_n} of the truncated inhomogeneous shear beam are tabulated
#' in [CylinderRoots] over the inhomogeneity exponent `m` and the truncation
#' ratio `l`. Method `"exact"` solves the characteristic relation of the beam
#' for every query; methods `"lm"` and `"nlm"` predict off-grid roots with a
#' polynomial surrogate fitted once to the whole table of the mode. The
#' root model is reused by [evaluateRoot()] for every query.
#'
#' @section Mathematics:
#' For the shear beam of height \eqn{H} truncated at \eqn{h = \lambda H}
#' with modulus \eqn{G \propto z^m}, Dakoulas and Gazetas (1985) reduce the
#' free vibration to a Bessel equation of order \eqn{q = m/(2-m)} and
#' obtain the characteristic relation (hydrated copy, `eq-free-characteristic`)
#' \deqn{J_{q+1}(a \lambda^{1-m/2})\, Y_q(a) - Y_{q+1}(a \lambda^{1-m/2})\, J_q(a) = 0,}
#' whose ordered positive roots are \eqn{a_n}, \eqn{n = 1, 2, \ldots}; for
#' \eqn{\lambda = 0} the relation reduces to \eqn{J_q(a) = 0}. Method
#' `"exact"` brackets the \eqn{n}-th sign change of the left-hand side on a
#' grid of step 0.02 between \eqn{a = 0.05} and \eqn{a = 100 / (1 - \lambda)}
#' and refines it with [stats::uniroot()] to a tolerance of \eqn{10^{-10}},
#' using R's fractional-order Bessel functions. The search is finite: a root
#' above the limit (for example mode 6 at \eqn{m = 1.8}, \eqn{\lambda = 0.5}, near 258.8)
#' is not found and [evaluateRoot()] stops with an error that names the limit.
#'
#' With \eqn{l_2 = l^2} and \eqn{m_2 = m^2}, the `"lm"` surrogate fits the
#' Gaussian linear model \eqn{a_n \sim l + m + l_2 + m_2} and the `"nlm"`
#' surrogate adds the interactions \eqn{l m} and \eqn{l_2 m_2}. These fits
#' are numerical interpolants of the table [CylinderRoots], whose roots the
#' `"exact"` method solves in R at every node of its grid; the fits
#' themselves do not solve the characteristic relation. The decision-tree
#' and random-forest surrogates of `newmark` 1.7.1 are not carried into this
#' version.
#'
#' @param mode One mode number, 1 to 8.
#' @param method The solution method of the root model: `"exact"`, or the surrogate `"lm"` or `"nlm"`.
#' @return The root model, a list with `mode`, `method`, `table` (the unique `(l, m, an)` rows
#'   of the mode), `fit` (the surrogate, `NULL` for `"exact"`), and the
#'   support `loMin`, `loMax`, `moMin` and `moMax`: for a surrogate the
#'   minimum and maximum `l` and `m` of the table, outside which predictions
#'   are extrapolations; for `"exact"` the domain of the relation,
#'   \eqn{0 \le \lambda < 1} and \eqn{0 \le m < 2}.
#' @export
#' @examples
#' RootModel <- buildRootModel(mode = 1, method = "nlm")
#' evaluateRoot(RootModel, mo = 0.45, lo = 0.20)
buildRootModel <- function(mode, method) {
  Roots <- .dataset("CylinderRoots")
  if (!is.numeric(mode) || length(mode) != 1L || !(mode %in% unique(Roots$n))) {
    stop("mode must be one tabulated mode number: ", paste(sort(unique(Roots$n)), collapse = ", "))
  }
  if (!is.character(method) || length(method) != 1L || !(method %in% c("exact", "lm", "nlm"))) {
    stop("method must be \"exact\", \"lm\" or \"nlm\"")
  }
  DATA <- unique(Roots[n == mode, -c("n")])
  if (method == "exact") {
    return(list(mode = mode, method = method, table = DATA, fit = NULL,
                loMin = 0, loMax = 1, moMin = 0, moMax = 2))
  }
  Fit <- if (method == "lm") {
    glm(an ~ l + m + l2 + m2, data = DATA[, .(an, l, m, l2 = l^2, m2 = m^2)])
  } else {
    glm(an ~ l + m + l2 + m2 + l * m + l2 * m2, data = DATA[, .(an, l, m, l2 = l^2, m2 = m^2)])
  }
  list(mode = mode, method = method, table = DATA, fit = Fit,
       loMin = min(DATA$l), loMax = max(DATA$l), moMin = min(DATA$m), moMax = max(DATA$m))
}

# Left-hand side of the characteristic relation of Dakoulas and Gazetas
# (1985) for inhomogeneity mo and truncation lo, evaluated at the
# dimensionless frequencies a; J_q(a) when the beam is not truncated.
.characteristicRelation <- function(a, mo, lo) {
  q <- mo / (2 - mo)
  if (lo == 0) return(besselJ(a, q))
  x <- a * lo^(1 - mo / 2)
  besselJ(x, q + 1) * besselY(a, q) - besselY(x, q + 1) * besselJ(a, q)
}

# The mode-th positive root of the characteristic relation: the mode-th sign
# change on a grid of step 0.02, refined by uniroot.
.solveCharacteristicRoot <- function(mode, mo, lo, amax = 100 / (1 - lo)) {
  Grid <- seq(0.05, amax, by = 0.02)
  Value <- .characteristicRelation(Grid, mo, lo)
  Change <- which(Value[-1L] * Value[-length(Value)] < 0)
  if (length(Change) < mode) stop("no root of mode ", mode, " below a = ", amax, " for mo=", mo, " lo=", lo)
  k <- Change[mode]
  uniroot(.characteristicRelation, lower = Grid[k], upper = Grid[k + 1L], mo = mo, lo = lo, tol = 1e-10)$root
}
