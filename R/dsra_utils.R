# Internal utilities from dsra

#' @keywords internal
#' @noRd
.mapLayers <- function(UID, hs, USCS = ValidUSCS) {
  DT <- data.table::data.table(UID, hs)
  Hs <- sum(hs)
  DT <- DT[, list(hs = round(sum(hs) / Hs * 100)), by = UID]
  DT <- data.table::dcast(DT, formula = "1 ~ UID", value.var = "hs", fill = 0)
  DT[, (USCS) := lapply(USCS, function(u) ifelse(u %in% names(DT), DT[[u]], 0))]
  DT[, .SD, .SDcols = USCS]
}

#' @keywords internal
#' @noRd
.find <- function(V, X) {
  V <- sort(unique(V))
  X <- sort(unique(X))
  results <- numeric(length(X))
  for (i in seq_along(X)) {
    x <- X[i]
    if (x %in% V) {
      results[i] <- x
    } else {
      idxBefore <- max(which(V < x), na.rm = TRUE)
      idxAfter  <- min(which(V > x), na.rm = TRUE)
      if (length(idxBefore) == 0) idxBefore <- idxAfter
      if (length(idxAfter)  == 0) idxAfter  <- idxBefore
      if (abs(V[idxBefore] - x) < abs(V[idxAfter] - x)) {
        results[i] <- V[idxBefore]
      } else if (abs(V[idxBefore] - x) > abs(V[idxAfter] - x)) {
        results[i] <- V[idxAfter]
      } else {
        results[i] <- mean(c(V[idxBefore], V[idxAfter]))
      }
    }
  }
  results
}
