#' Spectral quantiles at one oscillator period by log-log interpolation
#'
#' Interpolates every quantile curve of one spectrum slice (one site,
#' model, Vs30 and return period) at a target oscillator period, as the
#' displacement sampling does when the spectrum lacks the periods
#' \eqn{1.3 T_s} or \eqn{1.5 T_s}. Each curve is interpolated
#' independently; the slice is not grouped.
#'
#' @section Mathematics:
#' For the rows \eqn{(T_i, S_{a,i})} of one probability label,
#' \deqn{S_a(T) = \exp\, L\left(\ln T;\ \ln(T_i + 10^{-8}),\ \ln S_{a,i}\right)}
#' where \eqn{L} is linear interpolation (`stats::approx`) with constant
#' extension beyond the tabulated periods (`rule = 2`). The offset
#' \eqn{10^{-8}} admits a zero-period node; a target of zero falls below
#' every node and takes the lowest-period ordinate. Duplicate
#' \eqn{(T, p)} rows signal a warning and are averaged in
#' \eqn{\ln S_a}, a geometric mean. The `"mean"` rows follow the same rule
#' when there are two or more; a single mean row is copied to the target
#' period and a slice without one returns no mean row.
#'
#' @param spectrum A table of one spectrum slice with `Tn` (s), `p`
#'   (probability labels and optionally `"mean"`) and `Sa` (g). Every
#'   quantile label needs two or more rows; other columns are ignored.
#' @param Tn One finite non-negative target oscillator period in seconds.
#' @return A data.table with `Tn`, `p` and `Sa`, one row per quantile
#'   label plus the mean row when the slice has one, ordered by `Tn` and
#'   `p`.
#' @export
#' @examples
#' Spectrum <- data.frame(Tn = rep(c(0.2, 0.8), each = 3L),
#'                        p = rep(c("0.16", "0.84", "mean"), 2L),
#'                        Sa = c(1, 2, 1.5, 4, 8, 6))
#' interpolateSpectrum(Spectrum, Tn = 0.4)
interpolateSpectrum <- function(spectrum, Tn) {
  if (!is.data.frame(spectrum) || !all(c("Tn", "p", "Sa") %in% names(spectrum))) {
    stop("spectrum must be a table with columns Tn, p and Sa")
  }
  if (!is.numeric(Tn) || length(Tn) != 1L || !is.finite(Tn) || Tn < 0) {
    stop("Tn must be one finite non-negative period")
  }
  Target <- Tn
  DT <- as.data.table(spectrum)
  if (anyDuplicated(DT[, .(Tn, p)])) {
    warning("spectrum has duplicate (Tn, p) rows; they are averaged in log space. ",
            "Select one spectrum slice before interpolating.", call. = FALSE)
  }
  Labels <- sort(unique(DT[p != "mean"]$p))
  if (!length(Labels)) {
    stop("spectrum has no quantile rows (p other than mean)")
  }
  # approx() warns only about the tied periods already reported above.
  interpolate <- function(rows) {
    exp(suppressWarnings(stats::approx(x = log(rows$Tn + 1e-8), y = log(rows$Sa),
                                       xout = log(Target), rule = 2))$y)
  }
  OUT <- rbindlist(lapply(Labels, function(label) {
    Rows <- DT[p == label]
    if (nrow(Rows) < 2L) stop("spectrum needs two or more periods for p = ", label)
    data.table(Tn = Target, p = label, Sa = interpolate(Rows))
  }))
  Rows <- DT[p == "mean"]
  if (nrow(Rows) >= 2L) OUT <- rbindlist(list(OUT, data.table(Tn = Target, p = "mean", Sa = interpolate(Rows))))
  if (nrow(Rows) == 1L) OUT <- rbindlist(list(OUT, data.table(Tn = Target, p = "mean", Sa = Rows$Sa)))
  setorder(OUT, Tn, p)
  OUT[]
}
