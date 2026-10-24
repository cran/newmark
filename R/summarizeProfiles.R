#' Marginal statistics of a profile population
#'
#' Reduces the realizations of [simulateSoilProfiles()], of
#' [propagateSitePeriod()], of [propagateBeamPeriod()] or of any table with
#' one row per realization to
#' the requested quantiles and mean of the named numeric properties. The
#' statistics are marginal: a quantile row does not describe one column, and
#' the pairing between properties is lost in this reduction. Keep the
#' realization table when a later calculation needs paired values.
#'
#' @param profiles A table with one row per realization.
#' @param levels Probabilities in (0, 1) and/or the literal `"mean"`, in the
#'   order the rows are wanted.
#' @param properties Names of the numeric columns of `profiles` to summarize,
#'   for example `c("Go", "mo", "VSo")`, `"Tr"` or `"Ts"`.
#' @return A data.table with one row per requested level: `level` (the
#'   probability as text, or `"mean"`), `n` (realizations summarized) and the
#'   requested properties. Quantiles use type 7 of [stats::quantile()] on the
#'   stored values; the mean is arithmetic.
#' @export
#' @examples
#' set.seed(20260915)
#' Profiles <- simulateSoilProfiles(Hs = 10, USCS = c("ML", "MH"), NS = 5, h = 1,
#'                                  water = 0, POP = 100, Vref = 760)
#' summarizeProfiles(Profiles, levels = c(0.16, "mean", 0.84), properties = c("Go", "mo", "VSo"))
summarizeProfiles <- function(profiles, levels, properties) {
  if (!is.character(properties) || !length(properties) || anyNA(properties) || anyDuplicated(properties)) {
    stop("properties must name distinct columns to summarize")
  }
  if (!is.data.frame(profiles) || !nrow(profiles) || !all(properties %in% names(profiles))) {
    stop("profiles must contain rows with columns ", paste(properties, collapse = ", "))
  }
  if (!all(vapply(properties, function(p) is.numeric(profiles[[p]]), logical(1L)))) {
    stop("every property to summarize must be numeric")
  }
  Properties <- properties
  if (!length(levels) || anyNA(levels)) stop("levels must contain probabilities and/or \"mean\"")
  Level <- as.character(levels)
  IsMean <- Level == "mean"
  # Non-numeric tokens other than "mean" are rejected below; the coercion
  # warning would only duplicate that diagnostic.
  Prob <- suppressWarnings(as.numeric(Level))
  if (any(!IsMean & (!is.finite(Prob) | Prob <= 0 | Prob >= 1))) {
    stop("levels must be probabilities in (0, 1) or \"mean\"")
  }
  Level[!IsMean] <- as.character(Prob[!IsMean])
  DT <- as.data.table(profiles)
  LIST <- lapply(seq_along(Level), function(i) {
    ROW <- if (IsMean[i]) {
      DT[, lapply(.SD, mean), .SDcols = Properties]
    } else {
      DT[, lapply(.SD, quantile, probs = Prob[i], names = FALSE), .SDcols = Properties]
    }
    ROW[, `:=`(level = Level[i], n = nrow(DT))]
    ROW
  })
  OUT <- rbindlist(LIST)
  setcolorder(OUT, c("level", "n", Properties))
  OUT[]
}
