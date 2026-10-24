#' Hazard curves of a user workbook
#'
#' Reads hazard curves tabulated by hand in an Excel workbook: one sheet per
#' curve kind named `p=<kind>` (for example `p=mean` or `p=0.84`), a column
#' `Tn` with oscillator periods in seconds and one column per intensity level
#' in g, named by the level, holding annual rates of exceedance. Reading
#' needs the suggested package readxl.
#'
#' @section Mathematics:
#' Each cell is an annual rate \eqn{AEP} at the intensity level of its
#' column; [convertAEP()] gives the probability of exceedance over the
#' investigation time `IT` and the return period \eqn{1/AEP}. Cells with a
#' rate that is not positive are dropped. The curves are taken as supplied:
#' the workbook states which statistic of the hazard each sheet holds.
#'
#' @param file Path of the `.xlsx` workbook.
#' @param IT Investigation time in years over which the probabilities of
#'   exceedance are stated, finite and positive.
#' @return A data.table with one row per positive cell: `Tn` (s), `p` (the
#'   sheet's kind; a numeric kind written as a number, text otherwise), `IT`
#'   (years), `Sa` (g), `POE`, `AEP` (1/year) and `TR` (years). Sheets
#'   without intensity columns are skipped with a warning.
#' @export
#' @examples
#' # The mean hazard curves of the package tests: five periods by 39 levels, sheet p=mean.
#' Workbook <- system.file("extdata", "hazardWorkbook.xlsx", package = "newmark")
#' if (requireNamespace("readxl", quietly = TRUE)) readHazardWorkbook(Workbook, IT = 50)
readHazardWorkbook <- function(file, IT) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  if (!is.numeric(IT) || length(IT) != 1L || !is.finite(IT) || IT <= 0) {
    stop("IT must be one finite positive investigation time")
  }
  Sheets <- grep("^p=", readxl::excel_sheets(file), value = TRUE)
  if (!length(Sheets)) {
    stop("No sheets named 'p=...' found in: ", file)
  }
  rbindlist(lapply(Sheets, function(Sheet) {
    DT <- as.data.table(readxl::read_xlsx(file, sheet = Sheet))
    if (!"Tn" %in% names(DT)) {
      stop("Missing column 'Tn' in sheet '", Sheet, "'.")
    }
    COLS <- setdiff(names(DT), "Tn")
    if (!length(COLS)) {
      warning("No measure columns in sheet '", Sheet, "'. Skipping.")
      return(NULL)
    }
    DT <- melt(DT, id.vars = "Tn", measure.vars = COLS, variable.name = "Sa", value.name = "AEP",
               variable.factor = FALSE)[AEP > 0]
    Kind <- sub("p=", "", Sheet, fixed = TRUE)
    # Kinds such as "mean" are not numbers; the coercion warning carries no information.
    Number <- suppressWarnings(as.numeric(Kind))
    if (Kind != "mean" && !is.na(Number)) Kind <- as.character(Number)
    # Column names that are not numbers give NA levels, kept visible.
    DT[, `:=`(p = Kind, IT = IT, Sa = suppressWarnings(as.numeric(Sa)))]
    DT[, c("POE", "TR") := convertAEP(AEP, IT)]
    setcolorder(DT, c("Tn", "p", "IT", "Sa", "POE", "AEP", "TR"))
    DT
  }), use.names = TRUE, fill = TRUE)
}

# Names resolved by data.table inside the hazard-side functions.
utils::globalVariables("AEP")
