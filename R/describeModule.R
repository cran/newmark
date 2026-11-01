#' Describe the installed module
#'
#' Diagnoses the installation of this module without repairs; the installed
#' command prints this report for `doctor`.
#'
#' @param resources The resource directory of the installed command, or `NULL`
#'   outside an installed product.
#' @return A list with `schema`, `package`, `version`, `summary` and `checks`.
#'   Each check has `name`, `required`, `status` (`ok` or `fail`) and `detail`.
#'   The checks `library`, `exports`, `packages`, `tools` and `receipt` are
#'   required: the library loads with the version of `lib/DESCRIPTION`, the
#'   declared exports exist, the declared packages load, the declared tools are
#'   on `PATH` and the files of the installation receipt keep their recorded
#'   MD5.
#' @export
describeModule <- function(resources = NULL) {
  Checks <- list()
  # A failing check reports its condition as the detail, so the report is always complete.
  check <- function(name, required, expr) {
    Checks[[length(Checks) + 1L]] <<- tryCatch(
      list(name = name, required = required, status = "ok", detail = force(expr)),
      error = function(e) list(name = name, required = required, status = "fail", detail = conditionMessage(e)))
  }
  # The installed command supplies its resource directory; an R call outside it has none.
  runtime <- function(path) {
    if (is.null(resources)) stop("no installed resources: run the installed command")
    file.path(resources, path)
  }
  listed <- function(x) if (length(x)) paste(x, collapse = ", ") else "(none declared)"
  declared <- function(field) as.character(source(runtime("install/requirements.R"), local = new.env())$value[[field]])
  Version <- as.character(utils::packageVersion("newmark"))
  check("library", TRUE, {
    Declared <- unname(read.dcf(runtime("lib/DESCRIPTION"), fields = "Version")[1L, 1L])
    if (!identical(Version, Declared)) stop("newmark ", Version, " loads; lib/DESCRIPTION declares ", Declared)
    paste("newmark", Version, "loads")
  })
  check("exports", TRUE, {
    Exports <- declared("exports")
    Missing <- setdiff(Exports, getNamespaceExports("newmark"))
    if (length(Missing)) stop("missing: ", listed(Missing))
    paste("present:", listed(Exports))
  })
  check("packages", TRUE, {
    Packages <- declared("packages")
    Missing <- Packages[!vapply(Packages, requireNamespace, logical(1L), quietly = TRUE)]
    if (length(Missing)) stop("do not load: ", listed(Missing))
    paste("load:", listed(Packages))
  })
  check("tools", TRUE, {
    Tools <- declared("tools")
    Missing <- Tools[!nzchar(Sys.which(Tools))]
    if (length(Missing)) stop("not on PATH: ", listed(Missing))
    paste("on PATH:", listed(Tools))
  })
  check("receipt", TRUE, {
    Record <- jsonlite::read_json(runtime("install.json"), simplifyVector = TRUE)
    if (!identical(Record$package$name, "newmark")) stop("the receipt names ", Record$package$name)
    Sums <- unname(tools::md5sum(file.path(dirname(dirname(resources)), Record$file)))
    Changed <- Record$file[is.na(Sums) | Sums != Record$md5]
    if (length(Changed)) stop("changed or missing: ", listed(Changed))
    paste(length(Record$file), "files keep their recorded MD5")
  })
  Failed <- vapply(Checks, function(x) x$required && x$status == "fail", logical(1L))
  list(schema = 1L, package = "newmark", version = Version,
       summary = list(status = if (any(Failed)) "fail" else "ok"), checks = Checks)
}

#' @name describeModule
#' @examples
#' describeModule()
NULL
