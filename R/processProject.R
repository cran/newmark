#' Run the steps of a project
#'
#' Reads the project contract, plans the selected steps and publishes their
#' tables under `path.tables`, file by file, returning the record of
#' what it published. The steps are named after the table each one
#' publishes and always run in the order `ShearTable`, `DnTable`,
#' `kmaxTable`.
#'
#' Everything the run decides is resolved before its first calculation, and
#' a dry run makes the same checks: the contract keys of the selected steps,
#' the paths, the sites and the destinations. For `ShearTable`, the plan
#' also applies to the geometries, the materials, the pairs and `Vref.gmdp`
#' the rules that [buildShearTable()] applies to its arguments, and to the
#' USCS codes of every selected material the rule of [sampleSoilLayers()].
#' For `DnTable` and `kmaxTable`, the inputs that no selected step produces
#' must exist, and what the shear table and the MCE table hold is checked
#' before those steps calculate anything: for `DnTable`, the mean rows of
#' the shear table and the MCE demand of every selected site over the
#' period `1.5 * max(Ts)` of those rows, as [buildDnTables()] requires them;
#' for `kmaxTable`, the heights `Hs` from which [buildKmaxTable()] derives
#' its displacements. A dry run makes the checks of `DnTable` first and a
#' run the check of the heights first, as each phase of `newmark` 1.7.1 does,
#' so that when both fail the two report different failures. With
#' `ShearTable` selected as well, the shear table is the one that step
#' publishes: a run makes these checks after publishing it, and a dry run,
#' which publishes nothing, does not make them. A contract or an input that
#' these checks reject stops the run before the step it belongs to writes
#' anything, in a dry run as in a run.
#'
#' @section Project contract:
#' The contract is one JSON object with the keys of `newmark` 1.7.1 (schema
#' 1), validated strictly. Keys are looked up exactly, and a key outside the
#' list below stops the run naming it; that includes `Hs.min`, `Hs.max`,
#' `NTS`, `KMAX_PARALLEL` and `KMAX_WORKERS`, which 1.7.1 read. Each key is
#' required only by the steps named with it. `file` and the paths of the
#' `path` block are resolved against `root` when they are relative; the
#' working directory is never changed and reading evaluates no code. Numbers
#' are read as doubles whatever their JSON spelling, so `760` and `760.0`
#' are the same value.
#'
#' \describe{
#'   \item{`path.tables`}{Directory that `newmark` writes, and the only one:
#'     the global tables at its root and the tables of each site under
#'     `<ID.gmdp>/<siteID>/`. Every step; `ShearTable` creates it when it is
#'     absent.}
#'   \item{`path.data`}{Directory of the hazard table `MCETable.Rds`, read
#'     only. `DnTable`.}
#'   \item{`path.uhs`}{Directory holding `<ID.gmdp>/<siteID>/UHSRock.Rds`,
#'     read only; the sites are its directories that hold one. `DnTable`
#'     and `kmaxTable`.}
#'   \item{`path.calc`}{JSON of the classical hazard calculation, whose
#'     finite positive `job.calculation.truncation_level` truncates the
#'     sampling of the rock demand. `DnTable`.}
#'   \item{`Vref.gmdp`}{Reference rock velocity in m/s, one number. Every
#'     step; `DnTable` and `kmaxTable` require 760.}
#'   \item{`ID.gmdp`, `TR.gmdp`, `NS`, `subduction`, `Dn_weights`}{Identifier
#'     of the hazard model (a letter followed by letters, digits, `_` or
#'     `-`), return periods in years (sorted, duplicates dropped),
#'     realizations per cell (a positive whole number), the regime (`true`
#'     for subduction) and the 0/1 switches of the displacement models, which
#'     [selectDisplacementModels()] turns into the weights of the active
#'     models. `DnTable` and `kmaxTable`.}
#'   \item{`Mw.gmdp`}{Moment magnitude of the scenarios. `DnTable`.}
#'   \item{`DaH.gmdp`}{Admissible displacements, distinct positive
#'     percentages of the height `Hs`. `kmaxTable`.}
#'   \item{`DN_PARALLEL`, `DN_WORKERS`}{Optional parallel controls of
#'     `DnTable`: `true` (the default) or `false`, and a positive whole
#'     number of workers, by default the cores that R detects minus one, at
#'     least one. A string, a logical or a fraction is rejected. The workers
#'     evaluate the geometry and material pairs of a demand, so the step
#'     never uses more of them than the mean rows of the shear table, and
#'     `false` keeps it in one process; the tables do not depend on
#'     them.}
#'   \item{`uscs`, `IDm`}{Object of USCS code arrays indexed by material, and
#'     the distinct selected materials, each one a key of `uscs`.
#'     `ShearTable`.}
#'   \item{`Hs`, `lo.min`, `lo.max`, `s.min`, `s.max`}{Range form of the
#'     geometries: the heights in metres, and the first and last truncation
#'     ratio and slope, spread evenly over the heights in their declared
#'     order (they are not a minimum and a maximum) and passed to
#'     [buildGeometry()]. `ShearTable`, unless `geometryTable` is declared;
#'     the heights have no default.}
#'   \item{`geometryTable`}{Optional array of `IDg`, `Hs`, `b` and `s` row
#'     objects, used instead of the range form.}
#'   \item{`pairs`}{Optional array of `IDg` and `IDm` row objects restricting
#'     the product of geometries by materials.}
#'   \item{`seed`}{Optional whole number between 0 and 2147483647 that seeds
#'     the `ShearTable` step. A contract without it, such as one written for
#'     `newmark` 1.7.1, runs that step on the caller's random stream, as 1.7.1
#'     does without a seed, and the contract alone does not make its table
#'     reproducible.}
#' }
#'
#' @section Steps:
#' `ShearTable` builds the geometries, calls [set.seed()] with `seed` when
#' the contract defines it, which replaces the caller's random stream, calls
#' [buildShearTable()] with the levels 0.02, 0.05, 0.10, 0.16, mean, 0.84,
#' 0.90, 0.95 and 0.98, and publishes `<path.tables>/ShearTable.Rds` with
#' [saveRDS()]. Its message reports the seed, or `seed none` without one.
#'
#' `DnTable` reads `ShearTable.Rds` under `path.tables` (the one the
#' `ShearTable` step has just published, when both run), `MCETable.Rds` and
#' the `UHSRock.Rds` of each selected site, which it processes in
#' alphabetical order. For a site it calls [buildDnTables()] with `ID.gmdp`,
#' `TR.gmdp`, `Vref.gmdp`, `Mw.gmdp`, `NS`, the weights that
#' [selectDisplacementModels()] makes of `Dn_weights` and `subduction`, the
#' truncation level of `path.calc` and the workers, and [buildPGATable()],
#' and publishes under `<path.tables>/<ID.gmdp>/<siteID>/`, in this order,
#' `DnDraws.Rds` with [saveRDS()], `PGATable.csv` with
#' [data.table::fwrite()] and `na = "NA"`, and `DnTable.Rds` with
#' [saveRDS()]. It then rebuilds the global `DnPlotTable.Rds` and
#' `DnTable.Rds` with [assembleDnTables()] from the `DnTable.Rds` of every
#' site found under `path.uhs` that has one, selected or not, and publishes
#' them, in that order, at the root of `path.tables`; a site found
#' without that file stays out of them. [buildDnTables()] seeds each site
#' and demand with its own fixed seed, so the step replaces the caller's
#' random stream and the tables of a site do not depend on the sites run
#' with it.
#'
#' `kmaxTable` reads `ShearTable.Rds` under `path.tables` and, for each
#' selected site in alphabetical order, its `DnDraws.Rds` and its
#' `PGATable.csv`, read back with [data.table::fread()] and `TR` numeric.
#' For a site it calls [buildKmaxTable()] with `ID.gmdp`, `TR.gmdp`, `NS`,
#' the weights that [selectDisplacementModels()] makes of `Dn_weights` and
#' `subduction`, and `DaH.gmdp`, and publishes
#' `<path.tables>/<ID.gmdp>/<siteID>/kmaxTable.Rds` with [saveRDS()]. It
#' then assembles the global `kmaxTable.Rds` with [assembleKmaxTables()]
#' from the `kmaxTable.Rds` of every site found under `path.uhs` that has
#' one, selected or not, without reading their `PGATable.csv`, and publishes
#' it at the root of `path.tables`; a site found without `kmaxTable.Rds`
#' stays out of it. The warnings of both functions about roots beyond the
#' stored grid reach the caller, and the step draws no random number.
#'
#' The inputs of the site steps that no selected step produces must exist
#' when they are planned: `ShearTable.Rds` under `path.tables`,
#' `MCETable.Rds` for `DnTable`, and the `DnDraws.Rds` and `PGATable.csv`
#' of each selected site for `kmaxTable` without `DnTable`.
#'
#' @section Publication:
#' The tables are published in groups: the shear table; the three files
#' of a site and the two global files of `DnTable`; the site and the global
#' `kmaxTable`. Every file of a group is first written to a temporary file
#' `.<name>.<suffix>` beside its destination, in the format its name
#' declares, and then each destination is replaced by renaming its
#' temporary file, in the order of the group. Each rename is atomic; the
#' group is not, and nothing is rolled back. A failure while writing the
#' temporary files changes no file of the group; a failure at the rename of
#' a file leaves the files before it replaced and the others with their
#' previous content, or still absent; a failure in a later calculation or
#' check leaves everything published before it: the checks of the site
#' steps that follow a published shear table, the calculation of a later
#' site, or a table of another site that the global tables cannot take (a
#' `DnTable.Rds` that [assembleDnTables()] rejects, or a `kmaxTable.Rds`
#' that [assembleKmaxTables()] rejects), which stops the step with the
#' selected sites published and the
#' global files unchanged. The temporary files are removed when the
#' function exits, by return or by error, but not when the R process ends
#' between two renames. Running a step again with the same inputs publishes
#' the same tables, except `ShearTable` without `seed`.
#'
#' An error raised before the first file is published reaches the caller
#' unchanged, as the same condition object; that includes the error of the
#' writer of a temporary file. An error raised after it, while publishing,
#' calculating or validating, is signalled as a condition of class
#' `newmarkPartialPublication`. Its message names the step and site that
#' failed, the destination when the failure was a publication, the original
#' message, the files already published and the files of the failing group
#' that keep their previous content; its field `cause` holds the original
#' condition, `effects` the effects record of the files already published,
#' and `path` the destination whose temporary file or rename failed (`NULL`
#' when the failure was not a publication). A rename that fails returns
#' `FALSE`, with at most a warning, so it is reported as an error naming
#' its destination and the reason.
#'
#' @param file Path of the JSON contract, resolved against `root` when it is
#'   relative.
#' @param root Existing project root against which the contract and its
#'   relative paths are resolved; required, with no default. The CLI passes
#'   the working directory when `--root` is absent.
#' @param steps Any of `"ShearTable"`, `"DnTable"` and `"kmaxTable"`; they
#'   run in that order.
#' @param sites Identifiers of the sites to process; empty selects every
#'   site found under `path.uhs`. Only `DnTable` and `kmaxTable` process
#'   sites, and an unknown site stops them before anything is written. The
#'   global tables are rebuilt from the tables of every site found,
#'   selected or not.
#' @param dryRun `TRUE` makes the checks that a run makes before its first
#'   calculation (see Details), reports with [message()] the contract, the
#'   steps, the sites and every destination it would write or replace, and
#'   returns without calculating, seeding, creating a directory or writing
#'   anything.
#' @return Invisibly, the effects record: a data.table with one row per
#'   published file in publication order (with `dryRun`, one per
#'   destination) and the columns `step`, `siteID` (`NA` for a global
#'   table), `product`, `path` (the destination spelled from
#'   `path.tables` as the contract declares it, so relative to `root`
#'   when that path is relative), `action` (`"writes"` when the destination
#'   did not exist, `"replaces"` when it did) and `rows` (the rows of the
#'   published table; `NA` in a dry run).
#' @export
#' @examples
#' Root <- file.path(tempdir(), "newmarkProject")
#' dir.create(Root)
#' writeLines('{"path": {"tables": "data/newmark"}, "Vref.gmdp": 760,
#'   "Hs": [10, 12], "lo.min": 0.35, "lo.max": 0.075, "s.min": 2.5, "s.max": 3.5,
#'   "uscs": {"U1": ["ML", "MH"]}, "IDm": ["U1"], "seed": 20260915}',
#'   file.path(Root, "newmark.json"))
#' print(processProject(root = Root, steps = "ShearTable", dryRun = TRUE))
#' \donttest{
#' Effects <- processProject(root = Root, steps = "ShearTable")
#' Effects
#' readRDS(file.path(Root, Effects$path))[level == "mean"]
#' }
#' unlink(Root, recursive = TRUE)
processProject <- function(file = "newmark.json", root,
                           steps = c("ShearTable", "DnTable", "kmaxTable"),
                           sites = character(), dryRun = FALSE) {
  if (!is.logical(dryRun) || length(dryRun) != 1L || is.na(dryRun)) stop("dryRun must be TRUE or FALSE")
  Plan <- .planProject(file = file, root = root, steps = steps, sites = sites)
  if (dryRun) {
    # The checks of the site steps in the order of the dry run of newmark
    # 1.7.1: DnTable first, then the heights of kmaxTable.
    if (any(c("DnTable", "kmaxTable") %in% Plan$steps) && !"ShearTable" %in% Plan$steps) {
      Shear <- readRDS(Plan$shear)
      if ("DnTable" %in% Plan$steps) .readDnInputs(Plan, Shear)
      if ("kmaxTable" %in% Plan$steps) .checkKmaxHeights(Plan, Shear)
    }
    return(invisible(.describePlan(Plan)))
  }
  invisible(.publishRun(function(ledger) {
    if ("ShearTable" %in% Plan$steps) .runShearTable(Plan, ledger)
    if (!any(c("DnTable", "kmaxTable") %in% Plan$steps)) return()
    # The shear table of the site steps, the one that ShearTable has just
    # published when it runs, checked in the order of the run of newmark
    # 1.7.1: the heights of kmaxTable before anything of DnTable.
    ledger$step <- intersect(c("kmaxTable", "DnTable"), Plan$steps)[[1L]]
    Shear <- readRDS(Plan$shear)
    if ("kmaxTable" %in% Plan$steps) .checkKmaxHeights(Plan, Shear)
    if ("DnTable" %in% Plan$steps) .runDnTable(Plan, Shear, ledger)
    if ("kmaxTable" %in% Plan$steps) .runKmaxTable(Plan, Shear, ledger)
  }))
}

# Steps in their order of execution, and the products each one publishes:
# per site under <path.tables>/<ID.gmdp>/<siteID>/, then the globals at the
# root of path.tables, each list in publication order.
STEPS <- c("ShearTable", "DnTable", "kmaxTable")
PRODUCTFILES <- list(
  ShearTable = list(site = character(), global = c(ShearTable = "ShearTable.Rds")),
  DnTable = list(site = c(DnDraws = "DnDraws.Rds", PGATable = "PGATable.csv", DnTable = "DnTable.Rds"),
                 global = c(DnPlotTable = "DnPlotTable.Rds", DnTable = "DnTable.Rds")),
  kmaxTable = list(site = c(kmaxTable = "kmaxTable.Rds"), global = c(kmaxTable = "kmaxTable.Rds"))
)

# The quantile levels of the ShearTable product of newmark 1.7.1.
SHEARLEVELS <- c(.02, .05, .10, .16, "mean", .84, .90, .95, .98)

# Schema 1 of the contract: every key a step reads. Anything else is refused.
CONTRACTKEYS <- c("path", "Vref.gmdp", "ID.gmdp", "TR.gmdp", "NS", "Mw.gmdp", "subduction", "Dn_weights",
                  "DaH.gmdp", "uscs", "IDm", "geometryTable", "Hs", "lo.min", "lo.max", "s.min", "s.max",
                  "pairs", "seed", "DN_PARALLEL", "DN_WORKERS")
PATHKEYS <- c("tables", "data", "uhs", "calc")

# Everything the run decides before it calculates or writes: the contract of
# the selected steps, the paths it resolves, the sites, the inputs no
# selected step produces and every destination. It creates nothing.
.planProject <- function(file, root, steps, sites) {
  if (!is.character(steps) || !length(steps) || anyNA(steps) || any(!steps %in% STEPS)) {
    stop("steps must select ShearTable, DnTable or kmaxTable")
  }
  if (!is.character(sites) || anyNA(sites)) stop("sites must be a character vector of site identifiers")
  Root <- .resolveRoot(root)
  PATH.contract <- .resolveContract(file, Root)
  Contract <- .readContract(PATH.contract)
  Plan <- list(steps = STEPS[STEPS %in% steps], root = Root, contract = PATH.contract,
               products = .readPath(Contract, "tables", PATH.contract),
               Vref = .readNumber(Contract, "Vref.gmdp", PATH.contract), selected = character())
  if ("ShearTable" %in% Plan$steps) Plan <- c(Plan, .planShearTable(Contract, PATH.contract, Vref = Plan$Vref))
  if (any(c("DnTable", "kmaxTable") %in% Plan$steps)) Plan <- .planSites(Plan, Contract, sites)
  Plan$destinations <- .planDestinations(Plan)
  Plan
}

.planShearTable <- function(contract, path, Vref) {
  Seed <- if (.hasKey(contract, "seed")) .readWhole(contract, "seed", path, minimum = 0)
  Geometry <- if (.hasKey(contract, "geometryTable")) .readGeometry(contract, path) else {
    Hs <- .readNumbers(contract, "Hs", path)
    buildGeometry(Hs = Hs,
                  lo = seq(.readNumber(contract, "lo.min", path), .readNumber(contract, "lo.max", path),
                           length.out = length(Hs)),
                  s = seq(.readNumber(contract, "s.min", path), .readNumber(contract, "s.max", path),
                          length.out = length(Hs)))
  }
  USCS <- .readValue(contract, "uscs", path)
  if (!.isObject(USCS) || !length(USCS) ||
      !all(vapply(USCS, function(x) is.character(x) && length(x) > 0L && !anyNA(x) && all(nzchar(x)), logical(1L)))) {
    stop(path, " must define uscs as an object of non-empty arrays of USCS codes")
  }
  IDm <- .readValue(contract, "IDm", path)
  if (!is.character(IDm) || !length(IDm) || anyNA(IDm) || anyDuplicated(IDm)) {
    stop(path, " must define IDm as an array of distinct material identifiers")
  }
  if (length(setdiff(IDm, names(USCS)))) {
    stop("IDm in ", path, " names materials that uscs does not define: ",
         paste(setdiff(IDm, names(USCS)), collapse = ", "))
  }
  Pairs <- if (.hasKey(contract, "pairs")) .readValue(contract, "pairs", path)
  if (!is.null(Pairs) && (!is.data.frame(Pairs) || !is.character(Pairs[["IDg"]]) || !is.character(Pairs[["IDm"]]))) {
    stop(path, " must define pairs as an array of IDg and IDm row objects")
  }
  Problem <- .diagnoseShearArguments(Geometry, USCS = USCS, IDm = IDm, Vref = Vref, levels = SHEARLEVELS, pairs = Pairs)
  if (!is.null(Problem)) stop(path, " defines ShearTable inputs that buildShearTable() rejects: ", Problem)
  for (Material in IDm) {
    Problem <- .diagnoseUSCSPool(USCS[[Material]])
    if (!is.null(Problem)) stop(path, " defines uscs.", Material, " with codes that sampleSoilLayers() rejects: ", Problem)
  }
  list(seed = Seed, geometry = Geometry, USCS = USCS, IDm = IDm, pairs = Pairs)
}

# The keys, paths, sites and inputs of DnTable and kmaxTable: the run reads
# the shear table, the hazard tables and the site products named here.
.planSites <- function(plan, contract, sites) {
  PATH.contract <- plan$contract
  if (plan$Vref != 760) stop(PATH.contract, " must define Vref.gmdp as 760 m/s for the DnTable and kmaxTable steps")
  plan$model <- .readModel(contract, PATH.contract)
  plan$TR <- sort(unique(.readNumbers(contract, "TR.gmdp", PATH.contract)))
  if (any(plan$TR <= 0)) stop(PATH.contract, " must define TR.gmdp as positive return periods")
  plan$NS <- .readWhole(contract, "NS", PATH.contract, minimum = 1)
  Switches <- .readValue(contract, "Dn_weights", PATH.contract)
  if (!.isObject(Switches) || !length(Switches) ||
      !all(vapply(Switches, function(x) is.numeric(x) && length(x) == 1L && !is.na(x), logical(1L)))) {
    stop(PATH.contract, " must define Dn_weights as an object of model names and 0/1 switches")
  }
  Subduction <- .readFlag(contract, "subduction", PATH.contract)
  plan$weights <- tryCatch(
    selectDisplacementModels(switches = vapply(Switches, as.numeric, numeric(1L)), subduction = Subduction),
    error = function(e) stop("Dn_weights and subduction in ", PATH.contract, ": ", conditionMessage(e))
  )
  if ("DnTable" %in% plan$steps) {
    plan$Mw <- .readNumber(contract, "Mw.gmdp", PATH.contract)
    plan$truncation <- .readTruncation(.resolvePath(.readPath(contract, "calc", PATH.contract), plan$root))
    plan$MCE <- file.path(.resolvePath(.readPath(contract, "data", PATH.contract), plan$root), "MCETable.Rds")
    plan$parallel <- if (.hasKey(contract, "DN_PARALLEL")) .readFlag(contract, "DN_PARALLEL", PATH.contract) else TRUE
    # Without DN_WORKERS, the cores that R detects minus one, at least one,
    # as newmark 1.7.1 chose them; the step bounds them by its pairs.
    plan$workers <- if (.hasKey(contract, "DN_WORKERS")) {
      .readWhole(contract, "DN_WORKERS", PATH.contract, minimum = 1)
    } else {
      max(1L, parallel::detectCores(logical = TRUE) - 1L, na.rm = TRUE)
    }
  }
  if ("kmaxTable" %in% plan$steps) {
    plan$DaH <- .readNumbers(contract, "DaH.gmdp", PATH.contract)
    if (any(plan$DaH <= 0) || anyDuplicated(plan$DaH)) {
      stop(PATH.contract, " must define DaH.gmdp as distinct positive percentages of Hs")
    }
  }
  plan$rockDir <- file.path(.resolvePath(.readPath(contract, "uhs", PATH.contract), plan$root), plan$model)
  ProductsDir <- .resolvePath(plan$products, plan$root)
  plan$shear <- file.path(ProductsDir, "ShearTable.Rds")
  plan$modelDir <- file.path(ProductsDir, plan$model)
  if (!"ShearTable" %in% plan$steps) .checkInput(plan$shear, "written by the ShearTable step")
  if ("DnTable" %in% plan$steps) .checkInput(plan$MCE, "a hazard product of the project")
  plan$sites <- .listRockSites(plan$rockDir)
  if (length(setdiff(sites, plan$sites))) {
    stop("sites names sites without UHSRock.Rds under ", plan$rockDir, ": ",
         paste(setdiff(sites, plan$sites), collapse = ", "))
  }
  plan$selected <- if (length(sites)) plan$sites[plan$sites %in% sites] else plan$sites
  if ("kmaxTable" %in% plan$steps && !"DnTable" %in% plan$steps) {
    for (Site in plan$selected) {
      .checkInput(file.path(plan$modelDir, Site, "DnDraws.Rds"), "written by the DnTable step")
      .checkInput(file.path(plan$modelDir, Site, "PGATable.csv"), "written by the DnTable step")
    }
  }
  plan
}

# What DnTable reads and checks with the validators of its family before it
# calculates anything: the mean rows of the shear table and the MCE demand
# of every selected site over the period 1.5 * max(Ts) of those rows.
.readDnInputs <- function(plan, shear) {
  Pairs <- .checkDnShear(shear, plan$Vref)
  MCE <- readRDS(plan$MCE)
  for (Site in plan$selected) .checkDnMCE(MCE, Site, plan$Vref, 1.5 * max(Pairs$Ts))
  list(MCE = MCE, pairs = nrow(Pairs))
}

# The heights from which buildKmaxTable() derives its displacements,
# checked with the rule of its family before kmaxTable calculates anything.
.checkKmaxHeights <- function(plan, shear) {
  Problem <- .diagnoseKmaxHeights(shear)
  if (!is.null(Problem)) stop(plan$shear, " holds a shear table that buildKmaxTable() rejects: ", Problem)
}

# One row per destination of the selected steps, in publication order: the
# path as the record spells it and the file it resolves to.
.planDestinations <- function(plan) {
  DT <- rbindlist(lapply(plan$steps, function(step) {
    Files <- PRODUCTFILES[[step]]
    Site <- if (length(Files$site)) lapply(plan$selected, function(siteID) {
      data.table(step = step, siteID = siteID, product = names(Files$site),
                 path = file.path(plan$products, plan$model, siteID, Files$site))
    })
    rbindlist(c(Site, list(data.table(step = step, siteID = NA_character_, product = names(Files$global),
                                      path = file.path(plan$products, Files$global)))))
  }))
  set(DT, j = "file", value = .resolvePath(DT$path, plan$root))
  DT
}

# What a dry run reports. Whether a destination already exists is the only
# observation it makes of the tree it would replace.
.describePlan <- function(plan) {
  message("dry run: ", plan$contract)
  message("steps: ", paste(plan$steps, collapse = ", "))
  if (length(plan$selected)) message("sites: ", paste(plan$selected, collapse = ", "))
  Effects <- plan$destinations[, !"file"]
  set(Effects, j = "action", value = ifelse(file.exists(plan$destinations$file), "replaces", "writes"))
  set(Effects, j = "rows", value = NA_integer_)
  for (Line in paste0(Effects$step, ": ", Effects$action, " ", Effects$path)) message(Line)
  Effects
}

.runShearTable <- function(plan, ledger) {
  ledger$step <- "ShearTable"
  if (!is.null(plan$seed)) set.seed(plan$seed)
  Shear <- buildShearTable(plan$geometry, USCS = plan$USCS, IDm = plan$IDm, Vref = plan$Vref, levels = SHEARLEVELS,
                           pairs = plan$pairs)
  Destinations <- plan$destinations[plan$destinations$step == "ShearTable"]
  .createDirectory(dirname(Destinations$file))
  .publishGroup(ledger, Destinations, list(Shear))
  message("ShearTable: wrote ShearTable (", nrow(Shear), " rows, seed ",
          if (is.null(plan$seed)) "none" else plan$seed, ")")
}

# The products of each selected site, published group by group in the order
# of discovery, and then the global ones, rebuilt from the DnTable.Rds that
# every site found holds once the selected ones are published.
.runDnTable <- function(plan, shear, ledger) {
  ledger$step <- "DnTable"
  Inputs <- .readDnInputs(plan, shear)
  Workers <- if (plan$parallel) min(plan$workers, Inputs$pairs) else 1L
  Destinations <- plan$destinations[plan$destinations$step == "DnTable"]
  for (Site in plan$selected) {
    ledger$siteID <- Site
    Rock <- readRDS(file.path(plan$rockDir, Site, "UHSRock.Rds"))
    Products <- buildDnTables(Rock, MCE = Inputs$MCE, shear = shear, modelID = plan$model, siteID = Site,
                              TR = plan$TR, Vref = plan$Vref, Mw = plan$Mw, NS = plan$NS, weights = plan$weights,
                              zTrunc = plan$truncation, workers = Workers)
    Products$PGATable <- buildPGATable(Rock, MCE = Inputs$MCE, shear = shear, modelID = plan$model,
                                       siteID = Site, TR = plan$TR)
    Group <- Destinations[Destinations$siteID %in% Site]
    .createDirectory(dirname(Group$file[[1L]]))
    .publishGroup(ledger, Group, Products[Group$product])
    message("DnTable: wrote site ", Site)
  }
  ledger$siteID <- NA_character_
  FILES <- file.path(plan$modelDir, plan$sites, "DnTable.Rds")
  OK <- file.exists(FILES)
  Global <- assembleDnTables(stats::setNames(lapply(FILES[OK], readRDS), plan$sites[OK]), shear = shear,
                             weights = plan$weights, modelID = plan$model, TR = plan$TR, Vref = plan$Vref)
  Group <- Destinations[is.na(Destinations$siteID)]
  .publishGroup(ledger, Group, Global[Group$product])
}

# The kmaxTable of each selected site, published in the order of discovery,
# and then the global one, assembled from the kmaxTable.Rds that every site
# found holds once the selected ones are published.
.runKmaxTable <- function(plan, shear, ledger) {
  ledger$step <- "kmaxTable"
  Destinations <- plan$destinations[plan$destinations$step == "kmaxTable"]
  for (Site in plan$selected) {
    ledger$siteID <- Site
    Draws <- readRDS(file.path(plan$modelDir, Site, "DnDraws.Rds"))
    PGA <- .readPGATable(file.path(plan$modelDir, Site, "PGATable.csv"))
    Kmax <- buildKmaxTable(Draws, PGATable = PGA, shear = shear, modelID = plan$model, siteID = Site,
                           TR = plan$TR, NS = plan$NS, weights = plan$weights, DaH = plan$DaH)
    .publishGroup(ledger, Destinations[Destinations$siteID %in% Site], list(Kmax))
    message("kmaxTable: wrote site ", Site)
  }
  ledger$siteID <- NA_character_
  FILES <- file.path(plan$modelDir, plan$sites, "kmaxTable.Rds")
  OK <- file.exists(FILES)
  Global <- assembleKmaxTables(stats::setNames(lapply(FILES[OK], readRDS), plan$sites[OK]),
                               shear = shear, modelID = plan$model, TR = plan$TR, weights = plan$weights,
                               DaH = plan$DaH)
  .publishGroup(ledger, Destinations[is.na(Destinations$siteID)], list(Global))
}

# The PGA table of a site as newmark 1.7.1 read it back from its CSV, with
# TR numeric so that the MCE rows keep a missing number; an absent file
# names itself and the step that writes it.
.readPGATable <- function(path) {
  fread(.checkInput(path, "written by the DnTable step"), colClasses = c(TR = "numeric"))
}

# The boundary of a run's effects. `work` receives the ledger in which every
# published file is recorded; an error before the first publication reaches
# the caller as the same condition, one after it becomes a
# newmarkPartialPublication condition whose cause is that condition, and
# the record is the result.
.publishRun <- function(work) {
  Ledger <- new.env(parent = emptyenv())
  Ledger$effects <- data.table(step = character(), siteID = character(), product = character(),
                               path = character(), action = character(), rows = integer())
  Ledger$step <- NA_character_
  Ledger$siteID <- NA_character_
  tryCatch(work(Ledger), error = function(e) {
    Effects <- Ledger$effects
    if (!nrow(Effects)) stop(e)
    Kept <- if (!is.null(Ledger$group)) setdiff(Ledger$group$path, Effects$path)
    stop(errorCondition(
      paste0(Ledger$step, if (!is.na(Ledger$siteID)) paste0(", site ", Ledger$siteID), " failed",
             if (!is.null(Ledger$path)) paste0(" at ", Ledger$path), " after ", nrow(Effects),
             " published file(s): ", conditionMessage(e),
             "\npublished: ", paste(Effects$path, collapse = ", "),
             if (length(Kept)) paste0("\nkeep their previous content or remain absent: ", paste(Kept, collapse = ", "))),
      class = "newmarkPartialPublication", cause = e, effects = Effects, path = Ledger$path
    ))
  })
  Ledger$effects
}

# One group of products: every temporary file is written beside its
# destination first, in the format of the destination's name as newmark
# 1.7.1 wrote it (RDS with saveRDS, CSV with fwrite and na = "NA"), then
# each destination is replaced by a rename in the order of the group. The
# ledger records each file as it lands, and while the group is open it
# names the group and the destination being written, which is what a
# failure leaves behind. An error of the writer passes through as it was
# signalled; a rename that fails signals no error, so its failure is raised
# here naming the destination.
.publishGroup <- function(ledger, destinations, objects) {
  ledger$group <- destinations
  Temps <- tempfile(pattern = paste0(".", basename(destinations$file), "."), tmpdir = dirname(destinations$file))
  on.exit(unlink(Temps[file.exists(Temps)]), add = TRUE)
  for (i in seq_along(Temps)) {
    ledger$path <- destinations$path[[i]]
    switch(sub("^.*[.]", "", destinations$file[[i]]),
           Rds = saveRDS(objects[[i]], Temps[[i]]),
           csv = fwrite(objects[[i]], Temps[[i]], na = "NA"),
           stop("no product format for ", destinations$path[[i]]))
  }
  for (i in seq_along(Temps)) {
    ledger$path <- destinations$path[[i]]
    Action <- if (file.exists(destinations$file[[i]])) "replaces" else "writes"
    Failure <- tryCatch(if (!file.rename(Temps[[i]], destinations$file[[i]])) "the rename failed",
                        warning = conditionMessage)
    if (!is.null(Failure)) stop("could not replace ", destinations$path[[i]], ": ", Failure)
    ledger$effects <- rbind(ledger$effects,
                            cbind(destinations[i, c("step", "siteID", "product", "path")],
                                  action = Action, rows = nrow(objects[[i]])))
  }
  ledger$group <- NULL
  ledger$path <- NULL
  invisible(NULL)
}

.readJSON <- function(path) {
  tryCatch(jsonlite::read_json(path, simplifyVector = TRUE),
           error = function(e) stop(path, " is not valid JSON: ", conditionMessage(e)))
}

.isObject <- function(x) {
  is.list(x) && !is.data.frame(x) && !is.null(names(x)) && all(nzchar(names(x))) && !anyDuplicated(names(x))
}

.readContract <- function(path) {
  if (!file.exists(path) || dir.exists(path)) stop("project contract not found: ", path)
  Contract <- .readJSON(path)
  if (!.isObject(Contract)) stop(path, " must be one JSON object with uniquely named keys")
  Unknown <- setdiff(names(Contract), CONTRACTKEYS)
  if ("path" %in% names(Contract)) {
    if (!.isObject(Contract[["path"]])) stop(path, " must define path as an object with uniquely named keys")
    Unknown <- c(Unknown, paste0("path.", setdiff(names(Contract[["path"]]), PATHKEYS), recycle0 = TRUE))
  }
  if (length(Unknown)) stop(path, " declares keys that newmark does not read: ", paste(Unknown, collapse = ", "))
  Contract
}

# `key` is one name or the chain of names of a nested block.
.hasKey <- function(contract, key) {
  for (Name in key) {
    if (!is.list(contract) || !Name %in% names(contract)) return(FALSE)
    contract <- contract[[Name]]
  }
  TRUE
}

.readValue <- function(contract, key, path) {
  if (!.hasKey(contract, key)) stop(path, " must define ", paste(key, collapse = "."))
  for (Name in key) contract <- contract[[Name]]
  contract
}

.readNumbers <- function(contract, key, path) {
  x <- .readValue(contract, key, path)
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop(path, " must define ", paste(key, collapse = "."), " as numbers")
  }
  as.numeric(x)
}

.readNumber <- function(contract, key, path) {
  x <- .readNumbers(contract, key, path)
  if (length(x) != 1L) stop(path, " must define ", paste(key, collapse = "."), " as one number")
  x
}

# A count or a seed: a whole number however JSON spells it, within R's
# integer range.
.readWhole <- function(contract, key, path, minimum) {
  x <- .readValue(contract, key, path)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum || x > .Machine$integer.max || x != trunc(x)) {
    stop(path, " must define ", key, " as a whole number between ", minimum, " and ", .Machine$integer.max)
  }
  as.integer(x)
}

.readText <- function(contract, key, path) {
  x <- .readValue(contract, key, path)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(path, " must define ", paste(key, collapse = "."), " as one non-empty text")
  }
  x
}

.readFlag <- function(contract, key, path) {
  x <- .readValue(contract, key, path)
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(path, " must define ", key, " as true or false")
  x
}

# A path of the path block as the contract spells it, a home-relative one
# expanded.
.readPath <- function(contract, name, path) {
  x <- .readText(contract, c("path", name), path)
  if (.isAbsolutePath(x)) path.expand(x) else x
}

# ID.gmdp, the identifier of the hazard model, which names the directory of
# the site products under path.tables.
.readModel <- function(contract, path) {
  x <- .readText(contract, "ID.gmdp", path)
  if (!.isName(x)) stop(path, " must define ID.gmdp as an identifier: a letter followed by letters, digits, _ or -")
  x
}

# The numeric columns as doubles, so a geometry spelled without decimal
# points reaches the calculation with the values it declares.
.readGeometry <- function(contract, path) {
  Table <- .readValue(contract, "geometryTable", path)
  if (!is.data.frame(Table) || !is.character(Table[["IDg"]]) ||
      !all(vapply(c("Hs", "b", "s"), function(x) is.numeric(Table[[x]]), logical(1L)))) {
    stop(path, " must define geometryTable as an array of IDg, Hs, b and s row objects")
  }
  for (Column in names(Table)[vapply(Table, is.numeric, logical(1L))]) Table[[Column]] <- as.numeric(Table[[Column]])
  Table
}

.readTruncation <- function(path) {
  .checkInput(path, "the classical calculation of the project")
  Spec <- .readJSON(path)
  x <- .readNumber(Spec, c("job", "calculation", "truncation_level"), path)
  if (x <= 0) stop(path, " must define a positive job.calculation.truncation_level")
  x
}

.isAbsolutePath <- function(path) {
  grepl("^(/|~)", path) |
    (.Platform$OS.type == "windows" & grepl("^([A-Za-z]:[/\\\\]|\\\\\\\\)", path))
}

# The path of the contract that `file` names, resolved against the root
# when it is relative.
.resolveContract <- function(file, root) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) stop("file must be one non-empty path")
  .resolvePath(file, root)
}

.resolvePath <- function(path, root) {
  ifelse(.isAbsolutePath(path), path.expand(path), file.path(root, path))
}

.resolveRoot <- function(root) {
  if (!is.character(root) || length(root) != 1L || is.na(root) || !nzchar(root)) stop("root must be one non-empty path")
  if (!dir.exists(root)) stop("project root does not exist: ", root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

# An absent input names the file and what produces it, instead of R's
# connection error.
.checkInput <- function(path, producer) {
  if (!file.exists(path)) stop("missing input ", path, " (", producer, ")")
  invisible(path)
}

.createDirectory <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) stop("could not create ", path)
  invisible(path)
}

.listRockSites <- function(rockDir) {
  DIRS <- list.dirs(rockDir, recursive = FALSE, full.names = TRUE)
  Sites <- sort(unique(basename(DIRS[file.exists(file.path(DIRS, "UHSRock.Rds"))])))
  if (!length(Sites)) {
    stop("missing input UHSRock.Rds (a hazard product of the project): no site directory under ", rockDir,
         " holds one")
  }
  Bad <- Sites[!vapply(Sites, .isName, logical(1L))]
  if (length(Bad)) stop("invalid site directory names under ", rockDir, ": ", paste(Bad, collapse = ", "))
  Sites
}
