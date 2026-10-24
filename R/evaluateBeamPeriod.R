#' Modal period of the inhomogeneous truncated shear beam
#'
#' @section Mathematics:
#' Dakoulas and Gazetas (1985) measure the depth \eqn{z} from the apex where
#' the extended slopes meet: the crest lies at \eqn{z = h}, the base at
#' \eqn{z = H}, the modulus is \eqn{G = G_b (z/H)^m} and the truncation ratio
#' is \eqn{\lambda = h/H} (hydrated copy, Fig. 1 on printed page 167 and
#' equations `eq-modulus-depth-power` and `eq-free-truncation-ratio`). The
#' slope of [buildGeometry()] has height \eqn{H_s = H - h} from crest to
#' base, and its truncation ratio \eqn{\lambda_o = b/b_{max}} is that
#' \eqn{h/H}, because the extended slopes meet \eqn{h = b/(2s)} above the
#' crest; hence \eqn{H_s = H (1 - \lambda_o)}. The period of the mode whose
#' characteristic root is \eqn{a_n(m_o, \lambda_o)}, with basal velocity
#' \eqn{V_{So}}, is (`eq-natural-frequency-period`)
#' \deqn{T_n = \frac{4 \pi H}{a_n (2 - m_o) V_{So}}, \qquad H = \frac{H_s}{1 - \lambda_o}.}
#' With exact roots it reproduces the shear-beam periods of the five dams of
#' Table 8 of the paper (printed page 175) within 0.1 percent from their
#' printed height, truncation ratio, exponent and average velocity.
#'
#' The exponent \eqn{m_o} and velocity \eqn{V_{So}} that this package
#' estimates come from the fit of [characterizeProfile()], whose origin is
#' the crest (\eqn{G = G_o (z/H_s)^{m_o}}, \eqn{z} below the crest), while the
#' roots belong to the apex-origin profile of the paper. Evaluating the
#' period of a fitted profile this way is an approximation. The theory base
#' measures it on a continuous profile at 1.4, 3.7 and 6.0 percent for
#' \eqn{m = 0.5} and \eqn{\lambda} = 0.075, 0.20 and 0.35; for a layered
#' realization the measure is the difference that [propagateBeamPeriod()]
#' reports against the numerical period, with its divergence flag.
#'
#' The result is a conditional evaluation of its five inputs; evaluating it
#' on summary statistics of \eqn{m_o} and \eqn{V_{So}} does not give the
#' statistics of the period, which [propagateBeamPeriod()] obtains
#' realization by realization.
#'
#' @param Hs Slope heights from crest to base, in metres.
#' @param lo Truncation ratios \eqn{\lambda_o} in \eqn{[0, 1)}, aligned
#'   with `Hs`; see [buildGeometry()].
#' @param mo Inhomogeneity exponents, below 2.
#' @param VSo Basal shear-wave velocities in m/s.
#' @param an Characteristic roots of the requested mode at `mo` and `lo`.
#' @return Periods in seconds, aligned with the inputs, unrounded.
#' @export
#' @examples
#' # Dam A of Table 8 of Dakoulas and Gazetas (1985): 0.534 s, with the basal
#' # velocity derived from its average velocity of 200 m/s and the exact root.
#' evaluateBeamPeriod(Hs = 40, lo = 0.043, mo = 0.57, VSo = 228.25, an = 3.0123)
evaluateBeamPeriod <- function(Hs, lo, mo, VSo, an) {
  if (!is.numeric(Hs) || !is.numeric(lo) || !is.numeric(mo) || !is.numeric(VSo) || !is.numeric(an) ||
      any(!is.finite(c(Hs, lo, mo, VSo, an))) || any(Hs <= 0) || any(lo < 0 | lo >= 1) ||
      any(VSo <= 0) || any(an <= 0) || any(mo >= 2)) {
    stop("Hs, VSo and an must be finite and positive, lo in [0, 1) and mo finite and below 2")
  }
  H <- Hs / (1 - lo)
  4 * pi * H / (an * (2 - mo) * VSo)
}
