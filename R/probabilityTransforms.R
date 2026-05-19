#' Internal: Normalize investigation time vectors
#' @keywords internal
#' @noRd
.normalizeIt <- function(it, n) {
  IT <- as.numeric(it)
  if (!(length(IT) %in% c(1L, n))) {
    stop("Investigation time must be scalar or match the data length.")
  }
  if (length(IT) == 1L) {
    IT <- rep(IT, n)
  }
  if (any(!is.finite(IT) | IT <= 0, na.rm = TRUE)) {
    stop("Investigation time must be finite and > 0.")
  }
  IT
}

#' Internal: Convert POE to annual exceedance probability
#' @keywords internal
#' @noRd
poeToAep <- function(poe, it) {
  POE <- as.numeric(poe)
  IT <- .normalizeIt(it, length(POE))
  -log1p(-POE) / IT
}

#' Internal: Convert POE to return period
#' @keywords internal
#' @noRd
poeToTr <- function(poe, it) {
  POE <- as.numeric(poe)
  IT <- .normalizeIt(it, length(POE))
  IT / -log1p(-POE)
}

#' Internal: Convert annual exceedance probability to POE
#' @keywords internal
#' @noRd
aepToPoe <- function(aep, it) {
  AEP <- as.numeric(aep)
  IT <- .normalizeIt(it, length(AEP))
  -expm1(-IT * AEP)
}

#' Internal: Convert annual exceedance probability to return period
#' @keywords internal
#' @noRd
aepToTr <- function(aep) {
  1 / as.numeric(aep)
}
