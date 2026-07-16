#' Site-amplified spectral acceleration (NGA-East site-response models)
#'
#' Computes site-amplified spectral acceleration using two modes:
#' (A) If the input UHS contains only p == "mean", rock Sa is treated as
#'     deterministic and the dispersion in SaF comes solely from the
#'     site-response model (sdLnF), sampled on the natural-log scale.
#' (B) If the input UHS contains numeric quantiles (e.g., "0.10", "0.50", ...),
#'     ln Sa(Tn) is sampled from those quantiles using an empirical inverse
#'     quantile function (buildQSpline). Dependence between PGA and Sa(Tn)
#'     is represented by a per-period Gaussian "star" copula with the
#'     Baker & Jayaram (2008) correlation rhoBJ(Tn, T0 = 0), which has
#'     no distance dependence. For each Monte Carlo draw, the
#'     site-response model is evaluated with the paired PGA draw, and
#'     ln SaF = ln Sa + ln F with ln F ~ Normal(muLnF, sdLnF).
#'     Independence between ln Sa and ln F is assumed. Each sampled
#'     vector is re-centered to the tabulated mean of its marginal.
#'
#' For scenario (DSHA) tables with analytic mixture components available,
#' prefer \code{\link{fitSaFMixture}} (exact, sampling-free).
#'
#' Requirements:
#' - The input table must have columns Tn, p, Sa.
#' - It must include exactly one row (Tn == 0, p == "mean") providing PGA on rock.
#'
#' Output columns: Tn, p, Sa, SaF, AF, where AF = SaF / Sa.
#' In mode (B), the same numeric p values found in the input are returned plus
#' "mean". In mode (A), the set of p values is taken from p_TARGET plus "mean".
#'
#' @param uhs       data.table with columns Tn, p, Sa.
#' @param vs30      target Vs30 (m/s).
#' @param vref      reference Vs30 (m/s), default 760.
#' @param ns        Monte Carlo size (default 1000).
#' @param models    character. Site-response model id: \code{"ST20"}
#'                  (default; Stewart/Hashash 2020 NGA-East generation, see
#'                  \code{\link{F_ST20}}) or \code{"ST17"} (PEER 2017/05
#'                  generation, see \code{\link{F_ST17}}).
#' @param p_TARGET  probabilities to report when only 'mean' is provided,
#'                  default c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95).
#'
#' @return data.table with columns Tn, p, Sa, SaF, AF.
#' @examples
#' \dontrun{
#' uhs <- data.table::fread(
#'   system.file("extdata", "uhs.csv", package = "newmark"))
#' fitSaF(uhs, vs30 = 360, vref = 760, ns = 500)
#' }
#' @export
fitSaF <- function(uhs,
                   vs30,
                   vref = 760,
                   ns = 1000,
                   models = "ST20",
                   p_TARGET = c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)) {

  FSite <- switch(models, ST20 = F_ST20, ST17 = F_ST17,
                  stop("fitSaF: unknown site-response model: ", models))
  stopifnot(data.table::is.data.table(uhs))
  if (!all(c("Tn", "p", "Sa") %in% names(uhs))) {
    stop("fitSaF: 'uhs' must have columns Tn, p, Sa")
  }

  ## PGA (roca, media unica)
  PGA <- uhs[Tn == 0 & p == "mean", Sa]
  if (length(PGA) != 1L) {
    stop("fitSaF: need exactly one (Tn == 0, p == 'mean') row for PGA")
  }

  ## Deteccion del tipo de entrada (ramificacion unica)
  PNumAll <- suppressWarnings(as.numeric(uhs$p))
  HaveQuantiles <- any(!is.na(PNumAll) & PNumAll > 0 & PNumAll < 1)

  Tn <- sort(unique(uhs$Tn))
  OUT <- data.table::data.table(Tn = numeric(0), p = character(0),
                                Sa = numeric(0), SaF = numeric(0))

  if (!HaveQuantiles) {
    # ---------------------- (A) solo p == "mean" -----------------------
    for (t in Tn) {
      Sa <- uhs[Tn == t & p == "mean", Sa]
      if (length(Sa) != 1L) {
        next  # Manejar Tns sin datos
      }

      AUX  <- FSite(PGA = PGA, Tn = t, vs30 = vs30, vref = vref)
      muLnF <- AUX$muLnF
      sdLnF <- AUX$sdLnF

      X <- stats::rnorm(ns, mean = muLnF, sd = sdLnF)  # ln F
      Z <- log(Sa) + X                                 # ln SaF

      OUT <- data.table::rbindlist(list(
        OUT,
        data.table::data.table(
          Tn  = t,
          p   = sprintf("%.2f", p_TARGET),
          Sa  = Sa,
          SaF = as.numeric(stats::quantile(exp(Z), p_TARGET, type = 7))
        ),
        data.table::data.table(
          Tn  = t,
          p   = "mean",
          Sa  = Sa,
          SaF = mean(exp(Z))
        )
      ), use.names = TRUE)
    }

  } else {
    # ---------------------- (B) hay cuantiles numericos -----------------
    ## Preparar el muestreo de PGA una unica vez
    DP0 <- uhs[Tn == 0]
    p0  <- suppressWarnings(as.numeric(DP0$p))
    id0 <- which(!is.na(p0) & p0 > 0 & p0 < 1)
    DP0q <- DP0[id0, .(p = as.numeric(p), Sa)]
    DP0q <- DP0q[order(p)]
    LnQ0 <- buildQSpline(DP0q)  # u -> ln PGA

    for (t in Tn) {
      DT <- uhs[Tn == t]
      Sa <- DT[p == "mean", Sa]
      if (length(Sa) != 1L) {
        next  # Manejar Tns sin datos
      }

      ## inverso empirico en ln-espacio para Sa(Tn)
      PNum <- suppressWarnings(as.numeric(DT$p))
      idx  <- which(!is.na(PNum) & PNum > 0 & PNum < 1)
      DTn <- DT[idx, .(p = as.numeric(p), Sa)]
      DTn <- DTn[order(p)]
      LnQn <- buildQSpline(DTn)  # u in (0,1) -> ln Sa(Tn)

      ## copula gaussiana "estrella" PGA-Sa(Tn)
      rho <- rhoBJ(Tn = t, T0 = 0)
      Sig <- matrix(c(1, rho, rho, 1), nrow = 2L)
      Z2  <- mvtnorm::rmvnorm(ns, sigma = Sig)  # (Z0, Zn)
      U0  <- stats::pnorm(Z2[, 1L])
      U1  <- stats::pnorm(Z2[, 2L])

      ## mapeo a ln-marginales y recentering
      LnPGASample <- LnQ0(U0)
      LnPGASample <- log(exp(LnPGASample) * (PGA / mean(exp(LnPGASample))))

      LnSaSample  <- LnQn(U1)
      LnSaSample  <- log(exp(LnSaSample) * (Sa / mean(exp(LnSaSample))))

      ## modelo de sitio con PGA por muestreo
      AUX <- FSite(PGA = exp(LnPGASample), Tn = t, vs30 = vs30, vref = vref)
      LnFSample <- stats::rnorm(ns, mean = AUX$muLnF, sd = AUX$sdLnF)
      lnSaF <- LnSaSample + LnFSample

      ## mismos p que trae el hazard en este Tn + media
      OUT <- data.table::rbindlist(list(
        OUT,
        data.table::data.table(
          Tn  = t,
          p   = DT$p[idx],
          Sa  = DT$Sa[idx],
          SaF = as.numeric(stats::quantile(exp(lnSaF), PNum[idx], type = 7))
        ),
        data.table::data.table(
          Tn  = t,
          p   = "mean",
          Sa  = Sa,
          SaF = mean(exp(lnSaF))
        )
      ), use.names = TRUE)
    }
  }

  OUT[, AF := SaF / Sa]
  data.table::setorder(OUT, Tn, p)
  OUT[]
}
