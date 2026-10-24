#' Characteristic root for given inhomogeneity and truncation
#'
#' With a surrogate root model, returns the tabulated root when the query
#' falls exactly on the grid of [CylinderRoots] and the surrogate prediction
#' otherwise; with the `"exact"` method, solves the characteristic
#' relation of Dakoulas and Gazetas (1985) for every query (see
#' [buildRootModel()]). The diagnostics tell a consumer where each value
#' comes from.
#'
#' @param rootModel The root model that [buildRootModel()] returns.
#' @param mo Inhomogeneity exponents, aligned with `lo`; below 2 for `"exact"`.
#' @param lo Truncation ratios \eqn{\lambda = b / (2 H_s s + b)}, aligned;
#'   in \eqn{[0, 1)} for `"exact"`.
#' @return A data.table with `mo`, `lo`, `an` (the root), `exact` (`TRUE`
#'   when read from the table or solved from the relation) and
#'   `withinSupport` (`FALSE` when `lo` or `mo` lies outside the tabulated
#'   range of a surrogate, below its minimum or above its maximum, so the
#'   value is an extrapolation; both minima are 0 in [CylinderRoots]). A surrogate root is a
#'   numerical interpolation whose error is not characterized; an `"exact"`
#'   root is converged to \eqn{10^{-10}}.
#' @export
evaluateRoot <- function(rootModel, mo, lo) {
  if (!is.list(rootModel) || !all(c("table", "fit", "loMin", "loMax", "moMin", "moMax") %in% names(rootModel))) {
    stop("rootModel must be the result of buildRootModel()")
  }
  if (!is.numeric(mo) || !is.numeric(lo) || !length(mo) || length(mo) != length(lo) ||
      any(!is.finite(c(mo, lo)))) {
    stop("mo and lo must be finite vectors of equal length")
  }
  if (identical(rootModel$method, "exact")) {
    if (any(lo < 0 | lo >= 1 | mo < 0 | mo >= 2)) stop("the exact method needs lo in [0, 1) and mo in [0, 2)")
    an <- vapply(seq_along(mo), function(i) .solveCharacteristicRoot(rootModel$mode, mo[i], lo[i]), numeric(1L))
    return(data.table(mo = mo, lo = lo, an = an, exact = TRUE, withinSupport = TRUE))
  }
  # Exact grid hits are read by an equality join on the tabulated (l, m);
  # the remaining queries are predicted together by the surrogate.
  Query <- data.table(query = seq_along(mo), lo = lo, mo = mo)
  Hits <- rootModel$table[Query, on = .(l = lo, m = mo), .(query = i.query, an = x.an), nomatch = NA]
  if (anyDuplicated(Hits$query)) {
    Dup <- Hits$query[duplicated(Hits$query)][1L]
    stop("the root table holds several roots for lo=", lo[Dup], " mo=", mo[Dup])
  }
  an <- Hits$an
  Exact <- !is.na(an)
  if (any(!Exact)) {
    an[!Exact] <- unname(predict(rootModel$fit, newdata = data.table(l = lo[!Exact], m = mo[!Exact],
                                                                 l2 = lo[!Exact]^2, m2 = mo[!Exact]^2)))
  }
  data.table(mo = mo, lo = lo, an = an, exact = Exact,
             withinSupport = lo >= rootModel$loMin & lo <= rootModel$loMax & mo >= rootModel$moMin & mo <= rootModel$moMax)
}

#' @name evaluateRoot
#' @examples
#' RootModel <- buildRootModel(mode = 1, method = "nlm")
#' # A node of CylinderRoots, a query between its nodes and one beyond its largest m.
#' Query <- data.frame(mo = c(0, 0.148, 1.5), lo = c(0, 0.351, 0.351))
#' evaluateRoot(RootModel, mo = Query$mo, lo = Query$lo)
#' # The exact method solves the characteristic relation at the same queries.
#' evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = Query$mo, lo = Query$lo)
NULL
