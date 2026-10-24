#' Rigid sliding-block displacement
#'
#' Evaluates the conditional lognormal distribution of the Newmark
#' displacement of a rigid block with yield acceleration \eqn{k_y} under a
#' ground motion of peak acceleration PGA and Arias intensity \eqn{I_a},
#' with one of the three rigid models of the displacement ensemble, named
#' by its author-year identifier `model`.
#'
#' @section Mathematics:
#' With \eqn{r = k_y / PGA}, \eqn{D_n} in cm, PGA in g and \eqn{I_a} in m/s:
#'
#' `AM88`, Ambraseys and Menu (1988), with the ratio capped at 0.9999 so
#' that \eqn{1 - r > 0} when the yield acceleration reaches the peak:
#' \deqn{\log_{10} D_n = 0.90 + \log_{10}\left[(1 - r)^{2.53} r^{-1.09}\right], \qquad \sigma_{\log_{10} D_n} = 0.30}
#'
#' `JB07`, Jibson (2007), equation 10:
#' \deqn{\log_{10} D_n = 0.561 \log_{10} I_a - 3.833 \log_{10} r - 1.474, \qquad \sigma_{\log_{10} D_n} = 0.616}
#'
#' `SR08`, Saygili and Rathje (2008), equation 6 with the (PGA,
#' \eqn{I_a}) coefficients of tables 1 and 2 (the authors recommend the
#' (PGA, PGV) vector; this published variant keeps the rigid models on Arias
#' intensity):
#' \deqn{\ln D_n = 2.39 - 5.24 r - 18.78 r^2 + 42.01 r^3 - 29.15 r^4 - 1.56 \ln PGA + 1.38 \ln I_a, \qquad \sigma_{\ln D_n} = 0.46 + 0.56 r}
#'
#' The base-10 results are returned in natural logarithms
#' (\eqn{\times \ln 10}).
#'
#' The vectors a model reads, `PGA` and `ky` and, for `JB07` and `SR08`,
#' `AI`, are evaluated element by element: a vector of length one is
#' repeated for every element of the others, and the vectors longer than one
#' must have one common length; otherwise the function stops and names the
#' lengths.
#'
#' @param model One rigid model identifier: `"AM88"`, `"JB07"` or `"SR08"`.
#' @param PGA Peak ground accelerations in g, finite and positive.
#' @param ky Yield accelerations in g, finite and positive.
#' @param AI Arias intensities in m/s, finite and positive, for example
#'   [estimateAI()] of each PGA; `AM88` does not read them.
#' @return A data.table with `IDn` (the model), `muLnD` and `sigmaLnD`, the
#'   mean and standard deviation of \eqn{\ln D_n} (\eqn{D_n} in cm), one
#'   row per element of the longest vector the model reads.
#' @references Ambraseys, N. N. and Menu, J. M. (1988). Earthquake-induced
#'   ground displacements. *Earthquake Engineering & Structural Dynamics*,
#'   16(7), 985-1006. \doi{10.1002/eqe.4290160704}
#'
#'   Jibson, R. W. (2007). Regression models for estimating coseismic
#'   landslide displacement. *Engineering Geology*, 91(2-4), 209-218.
#'   \doi{10.1016/j.enggeo.2007.01.013}
#'
#'   Saygili, G. and Rathje, E. M. (2008). Empirical predictive models for
#'   earthquake-induced sliding displacements of slopes. *Journal of
#'   Geotechnical and Geoenvironmental Engineering*, 134(6), 790-803.
#'   \doi{10.1061/(ASCE)1090-0241(2008)134:6(790)}
#' @export
#' @examples
#' evaluateRigidDisplacement(model = "SR08", PGA = 0.5, ky = c(0.05, 0.1), AI = estimateAI(0.5))
evaluateRigidDisplacement <- function(model, PGA, ky, AI) {
  if (!is.character(model) || length(model) != 1L || !model %in% c("AM88", "JB07", "SR08")) {
    stop("model must be one of AM88, JB07 and SR08")
  }
  if (!is.numeric(PGA) || !is.numeric(ky) || !is.numeric(AI) || !length(PGA) || !length(ky) ||
      !length(AI) || any(!is.finite(c(PGA, ky, AI))) || any(c(PGA, ky, AI) <= 0)) {
    stop("PGA, ky and AI must be finite and positive")
  }
  if (model == "AM88") {
    .checkModelLengths(PGA = PGA, ky = ky)
    r <- pmin(ky / PGA, 0.9999)
    return(data.table(IDn = model, muLnD = (0.90 + log10((1 - r)^2.53 * r^-1.09)) * log(10),
                      sigmaLnD = 0.30 * log(10)))
  }
  .checkModelLengths(PGA = PGA, ky = ky, AI = AI)
  if (model == "JB07") {
    return(data.table(IDn = model, muLnD = (0.561 * log10(AI) - 3.833 * log10(ky / PGA) - 1.474) * log(10),
                      sigmaLnD = 0.616 * log(10)))
  }
  r <- ky / PGA
  data.table(IDn = model,
             muLnD = 2.39 - 5.24 * r - 18.78 * r^2 + 42.01 * r^3 - 29.15 * r^4 - 1.56 * log(PGA) + 1.38 * log(AI),
             sigmaLnD = 0.46 + 0.56 * r)
}

# The vectors one displacement model reads, named: each of length one, which
# is repeated, or of the common length of the others, so that no vector is
# recycled against another.
.checkModelLengths <- function(...) {
  Lengths <- lengths(list(...))
  if (length(unique(Lengths[Lengths != 1L])) > 1L) {
    stop("the vectors of one model must have length one or one common length: ",
         paste(names(Lengths), Lengths, sep = " = ", collapse = ", "))
  }
  invisible(NULL)
}
