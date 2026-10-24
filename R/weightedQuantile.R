# The two weighted quantile estimators of the pooled summaries, shared by the
# displacement and kmax families (design decision): the general
# summaries, summarizeEnsemble() and summarizeKmax(), merge ties at full
# precision; the compatibility routes of newmark 1.7.1, buildDnTables() and
# buildKmaxTable(), keep its 15-significant-digit text merge.

# The seventh sample quantile definition of Hyndman and Fan (1996), The American
# Statistician 50(4), 361-365 (type 7 of quantile() in R), extended to cumulative
# weights as the ensemble chapter states; sums avoid long double.
.weightedQuantile <- function(x, w, p) {
  IDX <- order(x)
  x <- x[IDX]
  Values <- unique(x)
  Cumulative <- Reduce(`+`, rowsum(w[IDX], match(x, Values))[, 1L], accumulate = TRUE)
  W <- Cumulative[length(Cumulative)]
  h <- 1 + (W - 1) * p
  a <- pmax(floor(h), 1)
  b <- a + 1
  b[b > W] <- W
  position <- function(rank) pmin(findInterval(rank, Cumulative, left.open = TRUE) + 1L, length(Values))
  f <- h - floor(h)
  (1 - f) * Values[position(a)] + f * Values[position(b)]
}

# The quantile of Hmisc::wtd.quantile() as newmark 1.7.1 called it: once the
# sorted pool holds a tie, Hmisc merges values through their text, so every
# value is read back from 15 significant digits.
.textWeightedQuantile <- function(x, w, p) {
  IDX <- order(x)
  x <- x[IDX]
  if (anyDuplicated(x)) x <- as.numeric(as.character(x))
  .weightedQuantile(x, w[IDX], p)
}
