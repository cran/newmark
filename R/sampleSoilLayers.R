#' Sample one synthetic soil column from a USCS material pool
#'
#' Discretizes a column of height `Hs` into layers of thickness `h` and draws
#' the index properties of each layer from the ranges tabulated by USCS
#' class. This is the stochastic step of the synthetic route; the physics of
#' the sampled layers is evaluated afterwards by [evaluateSoilLayers()].
#'
#' @section Mathematics:
#' The column has \eqn{N_L = \lceil H_s / h \rceil} layers of thickness `h`,
#' so its depth \eqn{N_L h} can exceed `Hs` when `Hs` is not a multiple of
#' `h`; layer `k` spans \eqn{[(k-1)h, kh]} with midpoint \eqn{z_m}. The class
#' of every layer is drawn uniformly, with replacement, from `USCS`. For a
#' layer of class `c`, with the intervals of `VoidRatiosUSCS` and
#' `UnitWeightRanges`,
#' \deqn{e_{min} \sim U(e_{min}^{lo}, e_{min}^{hi}),\quad
#'       e_{max} \sim U(e_{max}^{lo}, e_{max}^{hi}),\quad
#'       e_o \sim U(e_{min}, e_{max}),\quad
#'       \gamma_{sat} \sim U(\gamma_{sat}^{lo}, \gamma_{sat}^{hi}).}
#' Fine soils also draw their Atterberg limits: the liquid limit
#' \eqn{LL \sim U(50, 100)} for CH, MH and OH and \eqn{LL \sim U(8, 50)} for
#' CL, ML and OL; peat (PT) draws no liquid limit. With the Casagrande lines
#' \eqn{U = 0.9 (LL - 8)} and \eqn{A = 0.73 (LL - 20)} for \eqn{LL \ge 20}
#' (else \eqn{A = 0}), the plasticity index is \eqn{IP \sim U(A, U)} for
#' clays and \eqn{IP \sim U(0, A)} for silts and organic soils; when the A
#' line is null (peat, or a liquid limit below 20) that draw is the
#' degenerate \eqn{U(0, 0)}, so `IP` is 0 and the random stream is not
#' consumed, as in `newmark` 1.7.1.
#'
#' @section Random stream:
#' The draws consume the caller's random stream in the order class vector,
#' then for each layer in depth order \eqn{e_{min}}, \eqn{e_{max}}, \eqn{e_o},
#' \eqn{\gamma_{sat}}, and for fine soils `LL` then `IP`. This is the order
#' of the `newmark` 1.7.1 generator, so a seeded stream reproduces its columns.
#'
#' @param Hs One positive height of the column in metres.
#' @param USCS Character vector of USCS codes forming the material pool; each
#'   layer draws its class from this vector.
#' @param h One positive layer thickness in metres.
#' @return A data.table with one row per layer, from the surface down:
#'   `layer`, `USCS`, `hs` (thickness, m), `zo`, `zi`, `zm` (top, bottom and
#'   midpoint depth, m), `emin`, `eo`, `emax` (void ratios), `gsat`
#'   (saturated unit weight, kN/m3), `LL` and `IP` (percent; `NA` for coarse
#'   classes, `LL` also `NA` for peat, `IP` 0 where the A line is null).
#'   Every value is one realization, not a statistic.
#' @export
#' @examples
#' set.seed(1)
#' sampleSoilLayers(Hs = 5, USCS = c("ML", "MH"), h = 1)
sampleSoilLayers <- function(Hs, USCS, h) {
  if (!is.numeric(Hs) || length(Hs) != 1L || !is.finite(Hs) || Hs <= 0) {
    stop("Hs must be one finite positive height in metres")
  }
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("h must be one finite positive layer thickness in metres")
  }
  Problem <- .diagnoseUSCSPool(USCS)
  if (!is.null(Problem)) stop(Problem)
  NL <- ceiling(Hs / h)
  Class <- sample(USCS, size = NL, replace = TRUE)
  zo <- zi <- emin <- emax <- eo <- gsat <- LL <- IP <- rep(NA_real_, NL)
  for (k in seq_len(NL)) {
    # Depth coordinates accumulate layer by layer, as the generator did.
    zo[k] <- if (k > 1L) zi[k - 1L] else 0
    zi[k] <- zo[k] + h
    RANGE <- VoidRatiosUSCS[USCS == Class[k]]
    emin[k] <- runif(1, min = RANGE$eminMin, max = RANGE$eminMax)
    emax[k] <- runif(1, min = RANGE$emaxMin, max = RANGE$emaxMax)
    eo[k] <- runif(1, min = emin[k], max = emax[k])
    RANGE <- UnitWeightRanges[USCS == Class[k]]
    gsat[k] <- runif(n = 1, min = RANGE$gsatMin, max = RANGE$gsatMax)
    if (Class[k] %in% c("CH", "MH", "OH")) LL[k] <- runif(n = 1, min = 50, max = 100)
    if (Class[k] %in% c("CL", "ML", "OL")) LL[k] <- runif(n = 1, min = 8, max = 50)
    if (Class[k] %in% ValidFines) {
      ALine <- if (is.na(LL[k]) || LL[k] < 20) 0 else 0.73 * (LL[k] - 20)
      if (Class[k] %in% ValidClays) IP[k] <- runif(n = 1, min = ALine, max = 0.9 * (LL[k] - 8))
      if (Class[k] %in% c(ValidSilts, ValidOrganic)) IP[k] <- runif(n = 1, min = 0, max = ALine)
    }
  }
  data.table(layer = seq_len(NL), USCS = Class, hs = rep(h, NL), zo = zo, zi = zi,
             zm = (zi + zo) / 2, emin = emin, eo = eo, emax = emax, gsat = gsat,
             LL = LL, IP = IP)
}

# The rule of a material pool that sampleSoilLayers() applies before it
# draws, as the message it stops with, or NULL when `USCS` keeps it. The
# ShearTable step of processProject() applies it to the pool of every
# selected material before it seeds anything.
.diagnoseUSCSPool <- function(USCS) {
  if (!is.character(USCS) || !length(USCS) || anyNA(USCS) || !all(USCS %in% ValidUSCS)) {
    return(paste0("USCS must contain USCS codes among ", paste(ValidUSCS, collapse = ", ")))
  }
  NULL
}
