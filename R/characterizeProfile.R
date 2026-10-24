#' Characterize one evaluated soil column
#'
#' Reduces the layers of one column to the profile quantities used by the
#' shear-beam vibration and site-response calculations: the travel-time velocity of
#' the upper 30 m, its site class, the fitted power law of the stiffness
#' profile, the basal velocity, the constant-density Rayleigh estimate of the
#' site period and the thickness composition.
#'
#' @section Mathematics:
#' The stiffness profile is fitted by least squares in logarithmic
#' coordinates over the layers with positive modulus,
#' \deqn{\log G_m = \log G_o + m_o \log (z_m / H_s),}
#' so that \eqn{G_o} is the modulus at the base \eqn{z = H_s} and \eqn{m_o}
#' the exponent of a profile whose origin is the crest, \eqn{z} being the
#' depth below it. [evaluateBeamPeriod()] uses \eqn{m_o} with the roots of
#' the truncated shear beam, whose profile has its origin at the apex of the
#' extended slopes, as an approximation. The basal
#' velocity uses the saturated unit weight of the deepest layer,
#' \eqn{V_{So} = \sqrt{9.81 \cdot 1000 G_o / \gamma_{sat,N_L}}}. The
#' composition columns are thickness percentages.
#'
#' No value is rounded: every quantity is returned at the precision of its
#' calculation. The rounding of `newmark` 1.7.1, its basal velocity derived
#' from a rounded modulus and its threshold depth ratios `Z500` and `Z1000`
#' belong to the compatibility route of [buildShearTable()]. `fitSigma` is
#' the residual standard error of the stiffness fit in log units, a numerical
#' diagnostic separate from any physical dispersion.
#'
#' @param layers The result of [evaluateSoilLayers()], or a table with
#'   columns `USCS`, `hs`, `zm`, `gsat`, `Gm` and `VSm`.
#' @param Hs One positive column height in metres.
#' @param Vref One positive reference velocity in m/s completing the upper
#'   30 m when the column is shallower; see [estimateVS30()].
#' @return A one-row data.table: `NL` (layers), `Go` (MPa), `mo`, `VSo` and
#'   `VS30` (m/s), `SID`, `TrRayleigh` (s, the Rayleigh estimate of the site
#'   period with constant density of [estimateRayleighPeriod()], as `newmark`
#'   1.7.1 computed it per realization), `Gravels`, `Sands`,
#'   `Fines`, `Clays`, `Silts`, `Organic` (percent of thickness), `USCS`
#'   (sorted classes present, joined by `.`), `fitSigma` and `fitLayers`.
#'   Every value is a conditional evaluation of the supplied layers.
#' @export
#' @examples
#' set.seed(1)
#' Layers <- evaluateSoilLayers(sampleSoilLayers(Hs = 10, USCS = c("ML", "MH"), h = 1),
#'                              Hs = 10, water = 0, POP = 100)
#' characterizeProfile(Layers, Hs = 10, Vref = 760)
characterizeProfile <- function(layers, Hs, Vref) {
  Required <- c("USCS", "hs", "zm", "gsat", "Gm", "VSm")
  if (!is.data.frame(layers) || !nrow(layers) || !all(Required %in% names(layers))) {
    stop("layers must contain rows with columns ", paste(Required, collapse = ", "))
  }
  if (!is.numeric(Hs) || length(Hs) != 1L || !is.finite(Hs) || Hs <= 0) {
    stop("Hs must be one finite positive height in metres")
  }
  DT <- as.data.table(layers)
  NL <- nrow(DT)
  VS30 <- estimateVS30(hs = DT$hs, VSm = DT$VSm, Vref = Vref)
  FIT <- data.table(zm = DT$zm, Gm = DT$Gm)[Gm > 0][, .(X = log(zm / Hs), Y = log(Gm))]
  LM <- lm(Y ~ X, data = FIT)
  Go <- exp(unname(coef(LM)[1L]))
  mo <- unname(coef(LM)[2L])
  VSo <- sqrt(9.81 * Go * 1000 / DT$gsat[NL])
  TrRayleigh <- estimateRayleighPeriod(VSm = DT$VSm, hs = DT$hs, zm = DT$zm)
  share <- function(codes) sum(DT$hs[DT$USCS %in% codes]) / sum(DT$hs) * 100
  data.table(
    NL = NL, Go = Go, mo = mo, VSo = VSo, VS30 = VS30, SID = classifySite(VS30), TrRayleigh = TrRayleigh,
    Gravels = share(ValidGravels), Sands = share(ValidSands), Fines = share(ValidFines),
    Clays = share(ValidClays), Silts = share(ValidSilts), Organic = share(ValidOrganic),
    USCS = paste0(sort(unique(DT$USCS)), collapse = "."),
    fitSigma = sigma(LM), fitLayers = nrow(FIT)
  )
}
