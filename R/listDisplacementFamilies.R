#' Families of the displacement models
#'
#' Lists the families of displacement models that the displacement and
#' kmax tables summarize as trees: `ensemble` holds the six models, `rigid`
#' the models of [evaluateRigidDisplacement()] and `flexible` those of
#' [evaluateFlexibleDisplacement()]. [buildDnTables()] and
#' [buildKmaxTable()] pool the active members of every family that has one,
#' with their weights normalized, and write the pool with the family name
#' as `IDn`; [queryProjectTable()] marks a tree as extrapolated from the
#' calibration ranges of its members.
#'
#' @return A named list of three character vectors, `ensemble`, `rigid` and
#'   `flexible`, in this order: the identifiers of the models of each
#'   family, in the order in which [buildDnTables()] and [buildKmaxTable()]
#'   pool them.
#' @export
#' @examples
#' listDisplacementFamilies()
listDisplacementFamilies <- function() DNFAMILIES
