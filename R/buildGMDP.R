# nolint start
#' Build Ground-Motion Design Parameters (GMDP)
#'
#' Assembles an annual-exceedance-probability (AEP) table, uniform-hazard-spectrum
#' (UHS) table, and--where available--a magnitude-distance disaggregation table,
#' using OpenQuake or user-supplied hazard files.
#'
#' @param IDo    Character. Identifier for the GMDP. Default `"gmdp"`.
#' @param path   Character. Directory with hazard data files.
#' @param engine Character. `"openquake"` (default) or `"user"`.
#' @param vref   Numeric. Reference Vs30 (m/s). Default 760.
#' @param TRo    Numeric vector. Target return periods (years).
#'
#' @import data.table

#' @return A named list with `AEPTable`, `UHSTable`, and `RMwTable`
#'   (`NULL` if disaggregation is unavailable).
#'
#' @examples
#' \dontrun{
#' # Requires an 'OpenQuake' Engine output directory containing classical
#' # PSHA + disaggregation ZIPs; not bundled with the package.
#' out <- buildGMDP(path = "path/to/openquake/output", vref = 760)
#' }
#' @export
buildGMDP <- function(path,
                      IDo = "GEM",
                      engine = "openquake",
                      vref = 760,
                      TRo = c(100, 200, 500, 1000, 2000, 2500, 5000, 10000)) {
    OUT <- importAEPTable(path = path, engine = engine, vref = vref)
    AEPTable <- OUT$AEPTable
    ITo <- unique(AEPTable$ITo)[1]
    RMwTable <- NULL

    if (engine == "openquake") {
        message("> Building Disaggregation Hazard Table...")
        RMwTable <- tryCatch(
            importModel.oqRMw(path = path, ITo = ITo, vref = vref),
            error = function(e) {
                message(">> Skipping disagg (", e$message, ")")
                NULL
            }
        )
        if (!is.null(RMwTable)) {
            RMwTable <- RMwTable[, `:=`(SID = Vs30toSID(vref), Vs30 = vref, SM = engine, IT = ITo)]
        } else {
            message("> Disaggregation data not available.")
        }
    }

    message("> Re-mesh hazard data on a uniform TR grid ...")
    COLS <- setdiff(names(AEPTable), c("Sa", "POE", "AEP", "TR"))
    # Include all output kinds (mean, std, quantiles) as separate p-levels.
    # This ensures downstream stages can propagate non-mean curves without breaking.
    UHSTable <- AEPTable[Tn != -1, remeshGroup(.SD, TRo, .BY), by = COLS]

    UHSTable[, `:=`(Vref = vref, Vs30 = vref, AF = 1, SaF = Sa)]
    UHSTable[, ID := IDo]
    AEPTable[, ID := IDo]
    if (!any(UHSTable$Tn == 0)) {
        stop("Missing PGA (Tn==0). Check openQuake job.ini file.")
    }
    list(
        AEPTable = AEPTable,
        UHSTable = UHSTable,
        RMwTable = RMwTable
    )
}

#' Import OpenQuake or user AEP Data (Internal)
#'
#' @keywords internal
#' @noRd
importAEPTable <- function(path, engine = "openquake", vref = 760) {
    message("> Build AEP Table...")
    AEPTable <- NULL

    if (engine == "openquake") {
        message("> Unzip OQ data ...")

        FILES <- list.files(path, pattern = "\\.zip$", full.names = TRUE)
        if (!length(FILES)) {
            stop("No .zip files found in ", path)
        }

        for (z in FILES) {
            message("> Import AEP data from: ", z)
            tryCatch(
                {
                    TempDir <- tempdir()
                    unlink(TempDir, recursive = TRUE)
                    dir.create(TempDir, showWarnings = FALSE)
                    utils::unzip(z, junkpaths = TRUE, exdir = TempDir)

                    AUX <- importModel.oqAEP(path = TempDir, vref = vref)
                    AEPTable <- if (is.null(AEPTable)) {
                        AUX
                    } else {
                        data.table::rbindlist(list(AEPTable, AUX), use.names = TRUE, fill = TRUE)
                    }
                },
                error = function(e) {
                    message(">> Skipping file: ", z, " (", e$message, ")")
                }
            )
        }
    } else if (engine == "user") {
        message("> Import user-supplied AEP.xlsx ...")
        AEPTable <- tryCatch(
            importModel.userAEP(path, filename = "AEP.xlsx"),
            error = function(e) stop("Error in user AEP data: ", e$message)
        )
    } else {
        stop("Unknown engine: ", engine)
    }

    if (is.null(AEPTable) || !nrow(AEPTable)) {
        stop("No valid AEP data found in path: ", path)
    }
    AEPTable[, `:=`(Vref = vref, Vs30 = vref)]
    return(list(AEPTable = AEPTable[]))
}

#' Import OpenQuake AEP File (Internal)
#'
#' @keywords internal
#' @noRd
importModel.oqAEP <- function(path, vref) {
    DT <- data.table::data.table()

    FILES <- c(
        list.files(path, pattern = "quantile_curve"),
        list.files(path, pattern = "hazard_curve"),
        list.files(path, pattern = "hcurves-csv")
    )
    if (!length(FILES)) {
        stop("No hazard_curve, quantile_curve, or hcurves-csv files found in ", path)
    }

    for (.files in FILES) {
        HEADER <- readLines(file.path(path, .files), n = 1)
        p <- .extractQuantileFromHeader(HEADER)
        ITo <- .extractInvestigationTime(HEADER)
        Tn <- .extractTnFromHeader(HEADER)

        if (is.na(p) || is.na(ITo)) {
            stop("Malformed header in file: ", file.path(path, .files))
        }

        AUX <- data.table::fread(file.path(path, .files), skip = 1, header = FALSE, blank.lines.skip = TRUE)
        if (!nrow(AUX)) {
            warning("Empty data in file: ", file.path(path, .files))
            next
        }

        COLS <- unlist(AUX[1]) |> as.vector()
        data.table::setnames(AUX, COLS)
        AUX <- AUX[-1]

        poeCols <- grep("poe-", COLS, value = TRUE)
        if (!length(poeCols)) {
            stop("No 'poe-' columns in file: ", file.path(path, .files))
        }

        DATA <- data.table::melt(
            data = AUX,
            id.vars = c("lon", "lat", "depth"),
            measure.vars = poeCols,
            variable.name = "Sa",
            value.name = "POE"
        )
        DATA[, POE := as.numeric(POE)]
        DATA[, Sa := as.numeric(sub("^poe-", "", Sa))]
        DATA[, AEP := poeToAep(POE, ITo)]
        DATA[, TR := poeToTr(POE, ITo)]
        DATA[, `:=`(Tn = Tn, p = p, ITo = ITo)]

        DT <- data.table::rbindlist(list(DT, DATA), use.names = TRUE)
    }

    # remove infinite or negative TR
    DT <- DT[is.finite(TR) & TR > 0 & !is.na(TR)]

    if (!nrow(DT)) {
        stop("No valid hazard data in: ", path)
    }
    ## ---- Check: debe existir PGA ----------------------------------------
    if (!any(DT$Tn == 0)) {
        stop("PGA row not found after import; Check input files. Check job.ini.")
    }

    return(DT[])
}

#' Import OpenQuake Mag-Dist Disaggregation (Internal)
#'
#' @keywords internal
#' @noRd
importModel.oqRMw <- function(path, ITo, vref) {
    if (!dir.exists(path)) {
        stop("Path does not exist: ", path)
    }

    DIR <- file.path(path, paste0(".temp_oqRMw_", as.integer(Sys.time())))
    if (dir.exists(DIR)) {
        unlink(DIR, recursive = TRUE, force = TRUE)
    }
    dir.create(DIR, showWarnings = FALSE)

    FILES <- list.files(path, pattern = "\\.zip$", full.names = TRUE)
    if (length(FILES)) {
        for (.files in FILES) {
            utils::unzip(.files, exdir = DIR, junkpaths = TRUE)
        }
    }
    searchDir <- if (length(FILES)) DIR else path

    FILES <- list.files(searchDir, pattern = "Mag_Dist", full.names = TRUE)
    FILES <- FILES[!grepl("TRT", FILES)]
    if (!length(FILES)) {
        unlink(DIR, recursive = TRUE, force = TRUE)
        message("> No 'Mag_Dist' files found in path. Skipping...")
        return(NULL)
    }
    meanFiles <- grep("Mag_Dist-mean", FILES, value = TRUE)
    if (length(meanFiles)) {
        FILES <- meanFiles
    }

    DHT <- data.table::data.table()
    for (.files in FILES) {
        DT <- tryCatch(
            data.table::fread(.files, skip = 1, header = TRUE, blank.lines.skip = TRUE),
            error = function(e) NULL
        )
        if (is.null(DT) || !nrow(DT)) {
            next
        }

        requiredCols <- c("mag", "dist", "poe", "imt")
        missingCols <- setdiff(requiredCols, names(DT))
        if (length(missingCols)) {
            next
        }

        data.table::setnames(DT, old = c("mag", "dist", "poe"), new = c("Mw", "R", "POE"), skip_absent = TRUE)
        if ("iml" %in% names(DT)) {
            DT[, iml := NULL]
        }

        rlzCol <- grep("rlz|mean", names(DT), value = TRUE)
        if (length(rlzCol) == 1) {
            data.table::setnames(DT, old = rlzCol, new = "p")
        } else {
            next
        }

        DT[imt == "PGA", imt := "Sa(0.0)"]
        DT[, Tn := stringr::str_extract(imt, "(?<=\\()\\d+\\.*\\d*(?=\\))")]
        DT[, Tn := as.numeric(Tn)]
        DT[is.na(Tn), Tn := 0]
        DT[, IT := ITo]
        DT[, `:=`(AEP = poeToAep(POE, IT), TR = poeToTr(POE, IT))]
        DT <- DT[is.finite(TR) & TR > 0 & is.finite(AEP) & AEP > 0]

        DHT <- data.table::rbindlist(list(DHT, DT), fill = TRUE)
    }

    unlink(DIR, recursive = TRUE, force = TRUE)
    if (!nrow(DHT)) {
        message("No valid disagg data found.")
        return(NULL)
    }
    return(unique(DHT))
}

#' Internal: Extract Quantile from OpenQuake Header
#' @keywords internal
#' @noRd
.extractQuantileFromHeader <- function(line) {
    if (grepl("kind='mean'", line)) {
        return("mean")
    }
    if (grepl("kind='std'", line)) {
        # Standard deviation curve exported by OpenQuake
        return("std")
    }
    mt <- regexpr("kind='quantile-([0-9\\.]+)'", line)
    if (mt > 0) {
        val <- regmatches(line, mt)
        qval <- sub("kind='quantile-", "", val)
        qval <- sub("'", "", qval)
        return(as.numeric(qval))
    }
    return(NA)
}



#' @keywords internal
#' @noRd
.extractInvestigationTime <- function(line) {
    mt <- regexpr("investigation_time=([0-9\\.]+)", line)
    if (mt > 0) {
        val <- regmatches(line, mt)
        num <- sub("investigation_time=", "", val)
        return(as.numeric(num))
    }
    return(NA_real_)
}

#' @keywords internal
#' @noRd
.extractTnFromHeader <- function(line) {
    l <- tolower(line) # ignora mayusculas
    # --- PGA --------------------------------------------------------------
    if (grepl("imt=[\"']?pga[\"']?", l)) {
        return(0)
    }
    # --- PGV --------------------------------------------------------------
    if (grepl("imt=[\"']?pgv[\"']?", l)) {
        return(-1)
    }
    # --- SA(T) ------------------------------------------------------------
    mt <- regexpr("imt=[\"']?sa\\(([0-9\\.]+)\\)[\"']?", l)
    if (mt > 0) {
        val <- regmatches(l, mt)
        num <- sub("imt=[\"']?sa\\(", "", val)
        num <- sub("\\)[\"']?", "", num)
        return(as.numeric(num))
    }
    return(NA_real_)
}

#' Internal: Interpolate Hazard Curve to Uniform TR Grid
#' @keywords internal
#' @noRd
.formatRemeshContext <- function(groupBy = NULL) {
    if (is.null(groupBy) || !length(groupBy)) {
        return("hazard curve")
    }
    AUX <- vapply(names(groupBy), function(x) paste0(x, "=", as.character(groupBy[[x]])), character(1))
    paste(AUX, collapse = ", ")
}

.prepareCurveForRemesh <- function(TRi, Sai, groupBy = NULL) {
    DT <- data.table::data.table(TR = as.numeric(TRi), Sa = as.numeric(Sai))
    DT <- DT[is.finite(TR) & TR > 0 & is.finite(Sa) & Sa > 0]
    if (nrow(DT) < 2L) {
        stop("Need at least two positive (TR, Sa) points for ", .formatRemeshContext(groupBy), ".")
    }

    data.table::setorder(DT, TR)
    if (anyDuplicated(DT$TR)) {
        AUX <- DT[, .(
            Sa = max(Sa),
            nSa = data.table::uniqueN(signif(Sa, 12))
        ), by = TR]
        if (any(AUX$nSa > 1L)) {
            message(">> Collapsing duplicate TR nodes by keeping max Sa in ", .formatRemeshContext(groupBy), ".")
        }
        DT <- AUX[, .(TR, Sa)]
        data.table::setorder(DT, TR)
    }
    if (nrow(DT) < 2L) {
        stop("Need at least two unique TR nodes for ", .formatRemeshContext(groupBy), ".")
    }
    DT
}

.checkTrCoverage <- function(TRi, TRo, groupBy = NULL) {
    trMin <- min(TRi)
    trMax <- max(TRi)
    targetMin <- min(TRo)
    targetMax <- max(TRo)
    tol <- max(0.5, 1e-4 * max(c(abs(TRi), abs(TRo))))
    if (targetMin < trMin - tol || targetMax > trMax + tol) {
        stop(
            "Requested TR grid [", signif(targetMin, 8), ", ", signif(targetMax, 8),
            "] exceeds imported curve range [", signif(trMin, 8), ", ", signif(trMax, 8),
            "] for ", .formatRemeshContext(groupBy),
            ". Expand the oqt POE grid or narrow TRo instead of extrapolating."
        )
    }
    tol
}

reMeshCurve <- function(TRi, Sai, TRo, groupBy = NULL) {
    if (!length(TRo) || any(!is.finite(TRo) | TRo <= 0)) {
        stop("TRo must contain finite positive return periods.")
    }

    DT <- .prepareCurveForRemesh(TRi, Sai, groupBy)
    tol <- .checkTrCoverage(DT$TR, TRo, groupBy)
    TRfit <- TRo
    TRfit[TRfit < min(DT$TR)] <- min(DT$TR)
    TRfit[TRfit > max(DT$TR)] <- max(DT$TR)
    if (any(TRo < min(DT$TR) - tol | TRo > max(DT$TR) + tol)) {
        stop("Requested TR grid exceeds the tolerated interpolation support.")
    }

    logTarget <- stats::approx(
        x = log(DT$TR),
        y = log(DT$Sa),
        xout = log(TRfit),
        method = "linear",
        ties = "ordered"
    )$y
    if (any(!is.finite(logTarget))) {
        stop("Interpolation failed for ", .formatRemeshContext(groupBy), ".")
    }
    exp(logTarget)
}


#' Internal: Remesh a Group to Uniform TR Grid
#' @keywords internal
#' @noRd
remeshGroup <- function(.SD, TRo, groupBy = NULL) {
    if (!all(c("TR", "Sa") %in% colnames(.SD))) {
        stop("'.SD' must contain columns 'TR' and 'Sa'.")
    }

    data.table::setorder(.SD, TR)
    .SD <- .SD[is.finite(TR) & TR > 0]
    if (nrow(.SD) < 2) {
        return(NULL)
    }

    TRi <- .SD$TR
    Sai <- .SD$Sa
    SaStar <- reMeshCurve(TRi, Sai, TRo, groupBy)
    data.table(TR = TRo, Sa = SaStar)
}

#' Internal: Remesh a Group to Uniform TR Grid
#' @keywords internal
#' @noRd
importModel.userAEP <- function(path = NULL, filename = "AEP.xlsx") {
    if (is.null(path) || !dir.exists(path)) {
        stop("Invalid `path`: must be a directory.")
    }

    inputFile <- file.path(path, filename)
    if (!file.exists(inputFile)) {
        stop("File not found: ", inputFile)
    }

    sheetsAll <- readxl::excel_sheets(inputFile)
    sheetsP <- grep(pattern = "^p=", sheetsAll, value = TRUE)
    if (length(sheetsP) == 0) {
        stop("No sheets named 'p=...' found in: ", inputFile)
    }

    AT <- data.table()
    for (SHEET in sheetsP) {
        DT <- data.table::as.data.table(readxl::read_xlsx(inputFile, sheet = SHEET))
        if (!("Tn" %in% names(DT))) {
            stop("Missing column 'Tn' in sheet '", SHEET, "'.")
        }
        idVar <- "Tn"
        measureVars <- setdiff(names(DT), idVar)
        if (length(measureVars) == 0) {
            warning("No measure columns in sheet '", SHEET, "'. Skipping.")
            next
        }

        aux <- melt(
            DT,
            id.vars = idVar,
            measure.vars = measureVars,
            variable.name = "Sa",
            value.name = "AEP",
            variable.factor = FALSE
        )
        aux[, Sa := as.character(Sa)]

        po <- stringr::str_remove(SHEET, "p=")
        if (po != "mean") {
            poNum <- suppressWarnings(as.numeric(po))
            if (!is.na(poNum)) {
                po <- poNum
            }
        }

        aux <- aux[AEP > 0]
        if (!is.numeric(aux$Sa)) {
            aux[, Sa := suppressWarnings(as.numeric(as.character(Sa)))]
        }

        aux[, p := po]
        aux[, IT := 50]
        aux[, POE := aepToPoe(AEP, IT)]
        aux[, TR := aepToTr(AEP)]

        AT <- rbind(AT, aux, use.names = TRUE, fill = TRUE)
    }

    colOrder <- c("Tn", "Sa", "AEP", "p", "POE", "TR", "IT")
    colOrder <- intersect(colOrder, names(AT))
    setcolorder(AT, colOrder)

    return(AT[])
}
