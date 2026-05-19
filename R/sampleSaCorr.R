#' Draw a correlated sample of spectral accelerations
#'
#' Draw a correlated sample of spectral accelerations
#'
#' @param uhs data.table with columns Tn, p, Sa (must include p == "mean").
#' @param Tn  numeric vector; the first element is the reference period \code{Tn[1]},
#'            e.g., 0.01 or 0.00. It does not have to be sorted.
#' @param rho numeric vector of length \code{length(Tn) - 1}; correlations between the
#'            reference period \code{Tn[1]} and each remaining period.
#' @param NS  integer, Monte Carlo sample size.
#' @return matrix of dimension NS x length(Tn) with Sa draws (g).
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats    pnorm
#' @examples
#' \dontrun{
#' uhs <- data.table::data.table(
#'   Tn = c(0, 0, 0, 1, 1, 1),
#'   p  = c("0.16", "mean", "0.84", "0.16", "mean", "0.84"),
#'   Sa = c(0.3, 0.5, 0.7, 0.1, 0.2, 0.3)
#' )
#' rho <- rhoBJ(Tn = 1.0, T0 = 0.01)
#' mat <- sampleSaCorr(uhs, Tn = c(0, 1), rho = rho, NS = 100)
#' }
#' @export

sampleSaCorr <- function(uhs, Tn, rho, NS) {
  if (length(rho) != length(Tn) - 1L) {
    stop("length(rho) must be length(Tn) - 1")
  }

  ## 1) Correlation matrix (reference period is Tn[1])
  C <- diag(length(Tn))
  C[1, -1] <- C[-1, 1] <- as.numeric(rho)

  .sampleFromCorMatrix(C, uhs, Tn, NS)
}


.sampleSaCorrFull <- function(uhs, Tn, rho, NS) {
  stopifnot(is.matrix(rho), nrow(rho) == length(Tn), ncol(rho) == length(Tn))
  .sampleFromCorMatrix(rho, uhs, Tn, NS)
}


# Shared body of sampleSaCorr() and .sampleSaCorrFull(): given a fully
# constructed correlation matrix C, draw NS correlated Sa samples per
# Tn via Gaussian copula + empirical inverse + recentering to tabulated
# mean. Both wrappers above only differ in how they build C.
.sampleFromCorMatrix <- function(C, uhs, Tn, NS) {
  ## Ensure positive-definite numerically
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    C <- C + diag(abs(min(ev)) + 1e-8, nrow(C))
  }

  ## Gaussian copula -> uniforms
  Z <- mvtnorm::rmvnorm(NS, sigma = C)
  U <- stats::pnorm(Z)

  ## Map to Sa via empirical inverse; re-center to tabulated mean at each Tn
  out <- matrix(NA_real_, NS, length(Tn))
  for (j in seq_along(Tn)) {
    Dt <- uhs[Tn == Tn[j] & p != "mean", .(p, Sa)]
    if (nrow(Dt) < 2L) stop("Not enough quantiles at Tn = ", Tn[j])
    Dt[, p := as.numeric(p)]

    Qln <- buildQSpline(Dt)
    Q   <- exp(Qln(U[, j]))

    mu <- uhs[Tn == Tn[j] & p == "mean", Sa][1]
    if (!is.finite(mu)) mu <- mean(Q)
    out[, j] <- Q * (mu / mean(Q))
  }

  out
}
