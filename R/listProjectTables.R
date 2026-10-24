#' List the tables of a project
#'
#' Reads the project contract and lists the table files of
#' [processProject()] that exist under `path.tables`, reading each one as
#' its step wrote it to report its rows, columns and key. It writes nothing.
#'
#' The tables are the files at the destinations of [processProject()]:
#' `ShearTable.Rds`, `DnPlotTable.Rds`, `DnTable.Rds` and `kmaxTable.Rds` at
#' the root of `path.tables`, and `DnDraws.Rds`, `PGATable.csv`,
#' `DnTable.Rds` and `kmaxTable.Rds` under `<ID.gmdp>/<siteID>/` for every
#' directory there whose name is a site identifier (a letter followed by
#' letters, digits, `_` or `-`). A destination without its file is not
#' listed, and neither is any other file. A contract without `ID.gmdp`, as
#' one for the `ShearTable` step alone may be, has global tables only. The
#' contract is read and checked, and its paths resolved against `root`, as
#' [processProject()] does. The RDS files are read with [readRDS()] and
#' `PGATable.csv` with [data.table::fread()] and `TR` numeric, as the
#' `kmaxTable` step reads it back; a table file that cannot be read, or
#' that holds no table, stops the listing naming it.
#'
#' @param file Path of the JSON contract, resolved against `root` when it is
#'   relative.
#' @param root Existing project root against which the contract and its
#'   paths are resolved.
#' @return A data.table with one row per table file, in the order in which
#'   [processProject()] publishes them, and the columns `table`, `step`
#'   (the step that publishes it), `siteID` (`NA` for a global table),
#'   `path` (the file spelled from `path.tables` as the contract declares
#'   it, as in the effects record of [processProject()]), `bytes` (the size
#'   of the file), `rows` (the rows of the table it holds), `columns` and
#'   `key` (list columns: the column names of that table, and its key
#'   columns, empty when it has none) and `date` (the modification time of
#'   the file). Without tables, the table has these columns and no rows.
#' @export
#' @examples
#' Root <- file.path(tempdir(), "newmarkProducts")
#' dir.create(Root)
#' writeLines('{"path": {"tables": "data/newmark"}, "Vref.gmdp": 760,
#'   "Hs": [10, 12], "lo.min": 0.35, "lo.max": 0.075, "s.min": 2.5, "s.max": 3.5,
#'   "uscs": {"U1": ["ML", "MH"]}, "IDm": ["U1"], "seed": 20260915}',
#'   file.path(Root, "newmark.json"))
#' listProjectTables("newmark.json", root = Root)
#' \donttest{
#' processProject(root = Root, steps = "ShearTable")
#' listProjectTables("newmark.json", root = Root)
#' }
#' unlink(Root, recursive = TRUE)
listProjectTables <- function(file, root) {
  Plan <- .planProducts(file, root)
  Plan$selected <- .listProductSites(Plan)
  DT <- .planDestinations(Plan)
  OK <- file.exists(DT$file) & !dir.exists(DT$file)
  DT <- DT[OK]
  Info <- file.info(DT$file, extra_cols = FALSE)
  Tables <- lapply(DT$file, function(path) {
    Table <- .readProduct(path)
    list(rows = nrow(Table), columns = names(Table), key = as.character(key(Table)))
  })
  # key is an argument of data.table(): the column takes its name afterwards.
  OUT <- data.table(table = DT$product, step = DT$step, siteID = DT$siteID, path = DT$path, bytes = Info$size,
                    rows = vapply(Tables, `[[`, integer(1L), "rows"), columns = lapply(Tables, `[[`, "columns"),
                    KEY = lapply(Tables, `[[`, "key"), date = Info$mtime)
  setnames(OUT, "KEY", "key")
  OUT
}

# The products of a project as a plan of every step that processProject()
# lays out: the root, the path of the contract, path.tables as the
# contract spells it and, when the contract declares ID.gmdp, the model and
# the directory of the site products; `selected` names the sites whose
# products are laid out, none until the caller names them.
.planProducts <- function(file, root) {
  Root <- .resolveRoot(root)
  PATH.contract <- .resolveContract(file, Root)
  Contract <- .readContract(PATH.contract)
  Products <- .readPath(Contract, "tables", PATH.contract)
  Model <- if (.hasKey(Contract, "ID.gmdp")) .readModel(Contract, PATH.contract)
  list(steps = STEPS, root = Root, contract = PATH.contract, products = Products, model = Model,
       modelDir = if (!is.null(Model)) file.path(.resolvePath(Products, Root), Model), selected = character())
}

# The sites of the products: the directories under <path.tables>/<ID.gmdp>/
# whose names are site identifiers, in alphabetical order, as the sites of
# path.uhs are; none without a model.
.listProductSites <- function(plan) {
  if (is.null(plan$model)) return(character())
  DIRS <- list.dirs(plan$modelDir, recursive = FALSE, full.names = FALSE)
  sort(DIRS[vapply(DIRS, .isName, logical(1L))])
}

# A product as its step wrote it: an RDS file read with readRDS(), the CSV
# of the PGA table read back with TR numeric. A file that cannot be read, or
# that holds no table, stops the caller naming the file.
.readProduct <- function(path) {
  Table <- tryCatch(switch(sub("^.*[.]", "", path), Rds = readRDS(path), csv = .readPGATable(path),
                           stop("no product format for this name")),
                    error = function(e) stop("could not read ", path, ": ", conditionMessage(e)))
  if (!is.data.frame(Table)) stop(path, " does not hold a table")
  Table
}
