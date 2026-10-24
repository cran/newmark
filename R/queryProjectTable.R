#' Query a table of a project
#'
#' Reads one table that [processProject()] published under
#' `path.tables` and returns its rows whose `where` columns take any of the
#' given values, with the selected `columns`. It writes nothing.
#'
#' The table is named as in the effects record of [processProject()] and
#' in [listProjectTables()]: with `siteID = NULL`, one of the global tables
#' `ShearTable`, `DnPlotTable`, `DnTable` and `kmaxTable` at the root of
#' `path.tables`; with a site identifier, one of the tables `DnDraws`,
#' `PGATable`, `DnTable` and `kmaxTable` of that site under
#' `<ID.gmdp>/<siteID>/`, which locates them with the `ID.gmdp` of the
#' contract. The contract is read and checked, and its paths resolved
#' against `root`, as [processProject()] does, and the table is read as
#' [listProjectTables()] reads it.
#'
#' Each element of `where` names a column and gives, as text, the values it
#' may take; a row is returned when every named column takes one of its
#' values. The text is converted to the type of the column, as the values
#' arrive from the command line: a number as JSON spells it (`"100"`,
#' `"1.35"`, `"1e-3"`) for a double column and a whole number for an integer
#' column, `"NA"` being a missing value in either, such as the `TR` of the
#' MCE rows; a text column compares the text itself, so its levels are
#' written as stored, `"mean"` or `"0.50"` for `p`. Numbers match exactly: a
#' decimal names its nearest double, so the 17 significant digits of a
#' stored value, as the JSON of a query writes it, select that value. An
#' unknown table or site, a table whose file is absent, a column that
#' the table does not have, a text that does not convert and a value that
#' its column never takes stop the query naming them.
#'
#' A `DnTable` or `kmaxTable` gains the column `extrapolated`, derived as it
#' is read and never written: 1 when the `ky` of the row, or its `kmax`
#' other than 0, lies outside the calibration range that [buildKyLimits()]
#' gives its model at the mean PGA of the scenario (for `DnTable`, that of
#' the `PGATable` of its site, whose absence stops the query naming the
#' file) or, for a tree (`ensemble`, `rigid`, `flexible`), the range of any
#' of its models stored for the same site, scenario, geometry and material,
#' and 0 otherwise. The upper bound of `AM88`, the mean PGA, is compared with
#' a relative tolerance of `1e-9`: when that PGA exceeds 0.80 g the grid of
#' [buildKyGrid()] ends at it, and its last point can exceed the PGA, as
#' computed or as rounded in the CSV of `PGATable`, in the last significant
#' digits. `where` filters the column and `columns` selects it as a stored
#' column.
#'
#' @param file,root As for [listProjectTables()].
#' @param table Name of the table.
#' @param siteID `NULL` for a global table, or the identifier of the site
#'   whose table is read.
#' @param where Named list of non-empty character vectors, one per filtered
#'   column; `list()` keeps every row.
#' @param columns Names of the columns returned, in that order; `NULL`
#'   returns every column, `extrapolated` last in a `DnTable` or `kmaxTable`.
#' @return A data.table with the `columns` of the rows of the table that
#'   pass `where`, in their stored order and with their stored types.
#' @export
#' @examples
#' Root <- file.path(tempdir(), "newmarkQuery")
#' dir.create(Root)
#' writeLines('{"path": {"tables": "data/newmark"}, "Vref.gmdp": 760,
#'   "Hs": [10, 12], "lo.min": 0.35, "lo.max": 0.075, "s.min": 2.5, "s.max": 3.5,
#'   "uscs": {"U1": ["ML", "MH"]}, "IDm": ["U1"], "seed": 20260915}',
#'   file.path(Root, "newmark.json"))
#' \donttest{
#' processProject(root = Root, steps = "ShearTable")
#' queryProjectTable("newmark.json", root = Root, table = "ShearTable", siteID = NULL,
#'                   where = list(level = c("0.16", "mean", "0.84")),
#'                   columns = c("IDg", "level", "Ts", "VS30"))
#' }
#' unlink(Root, recursive = TRUE)
queryProjectTable <- function(file, root, table, siteID, where, columns) {
  if (!is.character(table) || length(table) != 1L || is.na(table) || !nzchar(table)) {
    stop("table must be one table name")
  }
  if (!is.null(siteID) && !.isName(siteID)) stop("siteID must be NULL for a global table or one site identifier")
  GlobalProducts <- unique(unlist(lapply(PRODUCTFILES, function(x) names(x$global))))
  SiteProducts <- unique(unlist(lapply(PRODUCTFILES, function(x) names(x$site))))
  if (!table %in% c(GlobalProducts, SiteProducts)) {
    stop("unknown table ", table, ": the tables are ", toString(union(GlobalProducts, SiteProducts)))
  }
  if (is.null(siteID) && !table %in% GlobalProducts) {
    stop(table, " is a table of a site: siteID must name the site")
  }
  if (!is.null(siteID) && !table %in% SiteProducts) stop(table, " is a global table: siteID must be NULL")
  if (!is.list(where) || is.data.frame(where) || (length(where) && !.isObject(where)) ||
      !all(vapply(where, function(x) is.character(x) && length(x) > 0L && !anyNA(x), logical(1L)))) {
    stop("where must be a list of non-empty character vectors named by distinct columns")
  }
  if (!is.null(columns) && (!is.character(columns) || !length(columns) || anyNA(columns) || anyDuplicated(columns))) {
    stop("columns must be NULL or distinct column names")
  }
  Plan <- .planProducts(file, root)
  if (!is.null(siteID)) {
    if (is.null(Plan$model)) stop(Plan$contract, " must define ID.gmdp, which locates the products of a site")
    Sites <- .listProductSites(Plan)
    if (!siteID %in% Sites) {
      stop("unknown site ", siteID, ": ",
           if (length(Sites)) paste0("the sites under ", Plan$modelDir, " are ", toString(Sites))
           else paste0("no site directory under ", Plan$modelDir))
    }
    Plan$selected <- siteID
  }
  DT <- .planDestinations(Plan)
  # Outside DT[...], where `siteID` would name its column.
  OK <- DT$product == table & is.na(DT$siteID) == is.null(siteID)
  DT <- DT[OK]
  .checkInput(DT$file, paste0("written by the ", DT$step, " step"))
  Table <- .readProduct(DT$file)
  if (!is.data.table(Table)) setDT(Table)
  # The cells are bound before the joins, where data.table would report an
  # error of their evaluation as a name not found. A kmax of 0 is the null
  # demand of solveKmax(), not a root at a yield acceleration.
  if (table == "kmaxTable") {
    Cells <- unique(Table[, c(CELLKEY, "PGA"), with = FALSE])
    set(Table, j = "extrapolated",
        value = as.numeric(Table$kmax > 0 & .isExtrapolated(Table, ky = Table$kmax, cells = Cells, path = DT$path)))
  }
  if (table == "DnTable") {
    Cells <- .readCellPGA(Plan, Table, path = DT$path)
    set(Table, j = "extrapolated", value = as.numeric(.isExtrapolated(Table, ky = Table$ky, cells = Cells, path = DT$path)))
  }
  Unknown <- setdiff(c(names(where), columns), names(Table))
  if (length(Unknown)) stop(DT$path, " has no column ", toString(Unknown))
  MASK <- rep_len(TRUE, nrow(Table))
  for (Name in names(where)) MASK <- MASK & .matchWhere(Table[[Name]], where[[Name]], name = Name, path = DT$path)
  COLS <- if (is.null(columns)) names(Table) else columns
  Table[MASK, COLS, with = FALSE]
}

# The scenario, geometry and material of a row of DnTable or kmaxTable: the
# cell whose mean PGA bounds AM88.
CELLKEY <- c("siteID", "ID", "TR", "IDg", "IDm")

# Whether the yield acceleration `ky` (g) of each row of a DnTable or
# kmaxTable lies outside the calibration range of buildKyLimits() at the
# mean PGA of its cell (`cells`: CELLKEY and PGA). A model has its own
# range; a tree of DNFAMILIES is outside when any of its models stored for
# the cell is, that is, outside the intersection of their ranges.
.isExtrapolated <- function(table, ky, cells, path) {
  Branches <- unique(table[, c(CELLKEY, "IDn"), with = FALSE])
  Members <- rbindlist(c(list(data.table(IDn = DNFAMILIES$ensemble, member = DNFAMILIES$ensemble)),
                         lapply(names(DNFAMILIES), function(x) data.table(IDn = x, member = DNFAMILIES[[x]]))))
  # Each branch of a cell with those of its members that the cell stores.
  Ranges <- Members[Branches, on = "IDn", nomatch = NULL, allow.cartesian = TRUE]
  Ranges <- Ranges[Branches, on = c(CELLKEY, member = "IDn"), nomatch = NULL]
  Ranges[cells, on = CELLKEY, PGA := i.PGA]
  Limits <- unique(Ranges[, .(PGA)])[, buildKyLimits(PGA), by = PGA]
  # The upper bound of AM88 is the mean PGA, at which the grid of
  # buildKyGrid() ends when it exceeds 0.80 g; that last point can exceed the
  # PGA, as computed or as rounded in PGATable.csv, in the last significant
  # digits, so the bound admits a relative excess of 1e-9. The other bounds
  # are compared exactly.
  Limits[IDn == "AM88", kyMax := kyMax * (1 + 1e-9)]
  Ranges <- Ranges[Limits, on = c("PGA", member = "IDn"), nomatch = NULL][
    , .(kyMin = max(kyMin), kyMax = min(kyMax)), by = c(CELLKEY, "IDn")]
  Bounds <- Ranges[table, on = c(CELLKEY, "IDn"), .(kyMin, kyMax)]
  if (anyNA(Bounds$kyMin)) {
    stop(path, " has IDn without a calibration range of ky: ", toString(unique(table$IDn[is.na(Bounds$kyMin)])))
  }
  ky < Bounds$kyMin | ky > Bounds$kyMax
}

# The cells of a DnTable with their mean PGA, read from the PGATable that
# the DnTable step publishes beside the DnTable of each of their sites.
.readCellPGA <- function(plan, table, path) {
  if (is.null(plan$model)) stop(plan$contract, " must define ID.gmdp, which locates the PGATable of the sites of ", path)
  plan$selected <- unique(table$siteID)
  DT <- .planDestinations(plan)
  DT <- DT[DT$product == "PGATable"]
  Cells <- rbindlist(lapply(DT$file, .readPGATable))[, c(CELLKEY, "PGA"), with = FALSE]
  Missing <- unique(table[, CELLKEY, with = FALSE])[!Cells, on = CELLKEY]
  if (nrow(Missing)) {
    stop(path, " has cells without a PGA in the PGATable of their site: ",
         paste(Missing[, paste0("siteID=", siteID, ", ID=", ID, ", TR=", TR, ", IDg=", IDg, ", IDm=", IDm)],
               collapse = "; "))
  }
  Cells
}

# The rows of a product column that take any of the values given as text,
# converted to the type of the column; a value that the column never takes
# stops the query naming it, as the text gave it.
.matchWhere <- function(column, text, name, path) {
  Values <- if (is.character(column)) text else .convertWhere(text, column, name = name, path = path)
  MASK <- column %in% Values
  Unknown <- !Values %in% column[MASK]
  if (any(Unknown)) {
    stop("where.", name, " gives values that column ", name, " of ", path, " never takes: ", toString(text[Unknown]))
  }
  MASK
}

# Text as the numbers of a numeric column: "NA" is a missing value and any
# other text must be a number as JSON spells it, a whole number within R's
# integer range for an integer column. jsonlite parses the numbers because
# it rounds a decimal to the nearest double, which as.numeric() misses by a
# unit in the last place for some 17-digit texts where R's long double adds
# no precision, as on macOS for aarch64.
.convertWhere <- function(text, column, name, path) {
  if (!is.numeric(column)) {
    stop("column ", name, " of ", path, " holds ", class(column)[[1L]], " values, which where does not filter")
  }
  Number <- text != "NA"
  Bad <- Number & !grepl(JSONNUMBER, text)
  OK <- Number & !Bad
  Values <- rep(NA_real_, length(text))
  if (any(OK)) {
    Values[OK] <- jsonlite::parse_json(paste0("[", paste(text[OK], collapse = ","), "]"), simplifyVector = TRUE)
  }
  if (is.integer(column)) {
    Bad <- Bad | OK & !(is.finite(Values) & Values == trunc(Values) & abs(Values) <= .Machine$integer.max)
  }
  if (any(Bad)) {
    stop("where.", name, " gives text that is not ", if (is.integer(column)) "a whole number" else "a number",
         " for column ", name, " of ", path, ": ", toString(text[Bad]))
  }
  if (is.integer(column)) as.integer(Values) else Values
}

# A number as JSON spells it.
JSONNUMBER <- "^-?(0|[1-9][0-9]*)([.][0-9]+)?([eE][+-]?[0-9]+)?$"

# Names resolved by data.table inside the query.
utils::globalVariables(c("kyMin", "kyMax"))
