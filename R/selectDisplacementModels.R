#' Displacement models active in a tectonic regime
#'
#' Turns the on/off switches of the six displacement models of a project
#' into the named weights of the active models. The flexible model comes in
#' two fits for two populations: `BM17`, Bray, Macedo and Travasarou (2018),
#' for subduction interface earthquakes, and `BM19`, Bray and Macedo (2019),
#' for shallow crustal ones; both are declared with the same switch and the
#' regime keeps the one fitted for it.
#'
#' @param switches Named switches, each 0 or 1, of distinct models among
#'   `AM88`, `JB07`, `BT07`, `SR08`, `BM17` and `BM19`, with `BM17` and
#'   `BM19` declared together and alike. A non-binary value would change
#'   the ensemble mean and quantiles without any product recording it, so it
#'   is rejected.
#' @param subduction One logical: `TRUE` keeps `BM17`, `FALSE` keeps `BM19`.
#' @return The named weights (all 1) of the active models, in the order of
#'   `switches`. A message reports which half of the flexible pair the regime
#'   switched off.
#' @export
#' @examples
#' selectDisplacementModels(switches = c(AM88 = 1, JB07 = 1, BT07 = 1, SR08 = 1, BM17 = 1, BM19 = 1),
#'                          subduction = FALSE)
selectDisplacementModels <- function(switches, subduction) {
  Pair <- c("BM17", "BM19")
  if (!is.logical(subduction) || length(subduction) != 1L || is.na(subduction)) {
    stop("subduction must be one non-missing logical")
  }
  if (!is.numeric(switches) || !length(switches)) {
    stop("switches must be a non-empty named numeric vector")
  }
  IDs <- names(switches)
  if (is.null(IDs) || anyNA(IDs) || any(!nzchar(IDs)) || anyDuplicated(IDs)) {
    stop("switches must have unique non-empty names")
  }
  Unknown <- setdiff(IDs, c("AM88", "JB07", "BT07", "SR08", "BM17", "BM19"))
  if (length(Unknown)) {
    stop("switches has unsupported models: ", paste(Unknown, collapse = ", "))
  }
  if (any(!is.finite(switches)) || !all(switches %in% c(0, 1))) {
    stop("every switch must be 0 or 1. Found: ",
         paste(paste(IDs, switches, sep = "="), collapse = ", "))
  }
  Have <- intersect(Pair, IDs)
  if (length(Have) == 1L) {
    stop("switches declares ", Have, " without ", setdiff(Pair, Have),
         ". Declare both alike; subduction picks which one runs.")
  }
  if (length(Have) == 2L && switches[["BM17"]] != switches[["BM19"]]) {
    stop("switches must declare BM17 and BM19 alike: they are the two fits ",
         "of the flexible model, subduction and crustal, and `subduction` selects which runs.",
         "\n  Set both to 1 to use the flexible pair, or both to 0 to skip it.")
  }
  if (length(Have) == 2L && switches[["BM17"]] == 1) {
    Off <- if (subduction) "BM19" else "BM17"
    switches[[Off]] <- 0
    message(">> Dn models: subduction = ", subduction, " -> ", Off,
            " off, ", setdiff(Pair, Off), " on")
  }
  Weights <- switches[switches > 0]
  if (!length(Weights)) {
    stop("switches must switch on at least one model")
  }
  Weights
}
