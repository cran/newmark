#' Installed module identity
#'
#' @return A list with the package name and installed version; consumed by the
#'   `doctor` operation of the command-line interface.
#' @export
describeModule <- function() {
  list(package = "newmark", version = as.character(utils::packageVersion("newmark")))
}

#' @name describeModule
#' @examples
#' describeModule()
NULL
