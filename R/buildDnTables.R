# SPDX-FileCopyrightText: 2026 Alejandro Verri Kozlowski <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of newmark: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
#' Site displacement tables with soil-column uncertainty
#'
#' Builds `DnDraws` and `DnTable` from the paired Ts and VS30 values of
#' `ShearDraws`. `PGATable` is built separately by [buildPGATable()] on the
#' mean rows of `ShearTable`. For every return period
#' of the probabilistic demand and for the MCE, and for every geometry and
#' material pair, it amplifies the rock demand to the VS30 of each column
#' rounded to the nearest 1 m/s for spectrum lookup,
#' samples its scenario intensities at that column's Ts with
#' [sampleIntensities()], evaluates the active models with
#' [evaluateDisplacements()] and summarizes them with
#' [summarizeDisplacements()] and, for the families `ensemble`, `rigid` and
#' `flexible`, the weighted summary of [summarizeEnsemble()] with the
#' quantile of `newmark` 1.7.1 (see below).
#'
#' @section Population and summaries:
#' Pairs are ordered by `IDg` and `IDm`. Within each pair, the `NS`
#' movements cycle through the columns in increasing `realization` order;
#' 200 movements and 100 columns give two movements per column. Each
#' movement uses the same column for site response, intensity sampling and
#' displacement evaluation, and keeps its identity over models and ky.
#' Each measure (PGA, Sa13 and Sa15) is centered once over the complete
#' `NS` sample on the column means weighted by their assigned movements.
#' The column subsets are not centered separately. Innovations are drawn in
#' movement order, so identical columns reproduce a single-column sample.
#' The `mean` rows of `ShearTable` supply the reference `Ts` and `Vs30`
#' metadata in `DnTable` and the PGA from which its common ky grid is built.
#' Site response for a mean row uses that row's stored VS30. Column spectrum
#' lookup uses `round(VS30)`, with a 1 m/s step; `ShearDraws` and the table
#' reference metadata retain their stored values.
#' The probabilistic demands sample the analytic mixture of the amplified
#' spectrum, truncated at `zTrunc` like the hazard calculation; the MCE
#' samples its amplified quantile spectrum untruncated. The yield grid is
#' [buildKyGrid()] with 30 points on the mean PGA of the amplified spectrum.
#' Each demand reseeds the generator with a fixed seed of the site, the
#' demand and the return period (a base-31 hash of
#' `"<siteID>/<demand>/<TR>"` modulo \eqn{2^{31} - 1}), after the
#' deterministic site response and before its draws, so a site's products do
#' not depend on the sites or demands run with it (the session's random
#' stream is replaced, as the dn step of 1.7.1 replaced it); within a demand the pairs
#' draw in order, each its `NS` intensity realizations and then its `NS`
#' model deviates. Quantile levels are 0.05, 0.10, 0.16, 0.50, 0.84, 0.90 and
#' 0.95 plus the mean, labelled `"%.2f"`; `BM19` takes ordinary motions
#' (`NFC = "D100"`, `pulse = FALSE`). The MCE rows carry `ID = "MCE"` and
#' `TR = NA`; PGA and ky are in g, Dn in cm. The family quantiles are those
#' of `Hmisc::wtd.quantile()` as 1.7.1 called it: once the pool of one yield
#' acceleration holds a tie, every pooled value is read back from its
#' 15-significant-digit text before the quantile, while the mean uses the
#' exact values; [summarizeEnsemble()], the general summary, merges ties at
#' full precision.
#'
#' @param rock The `UHSRock` table of the site: one lognormal distribution at each oscillator period
#'   per return period at Vs30 = Vref = 760 m/s (`ID`, `siteID`, `TR`, `Tn`,
#'   `Sa`, `Vs30`, `Vref`, `SaF`, `AF`, `k`, `muLn`, `sigmaLn`, `w`,
#'   `fitErr`).
#' @param MCE The `MCETable`: rock quantiles and mean of the MCE (`ID`,
#'   `siteID`, `Tn`, `p`, `Sa`, `Vs30`, `Vref`), without `TR`.
#' @param shear The `ShearTable`, whose `mean` rows define the cell
#'   reference properties and common ky grid.
#' @param modelID Identifier of the probabilistic demand model.
#' @param siteID Identifier of the site.
#' @param TR Return periods of the probabilistic demand, in years.
#' @param Vref Reference rock velocity in m/s; the product contract admits
#'   only 760.
#' @param Mw Moment magnitude of the scenarios: one number for every demand, or a numeric vector
#'   named by the return periods of `TR`, written as they print, and `MCE`.
#' @param NS Realizations per geometry, material and demand.
#' @param weights Named positive weights of the active models, for example
#'   from [selectDisplacementModels()].
#' @param zTrunc Truncation level of the probabilistic rock demand, in
#'   standard deviations.
#' @param workers Number of forked processes that evaluate the pairs of a
#'   demand after their draws; 1 evaluates them in this process. The
#'   products do not depend on it.
#' @param draws The `ShearDraws` table with `IDg`, `IDm`, `realization`,
#'   `Ts`, `VS30` and `Vref`. Each pair of `shear` must be represented,
#'   with realization identifiers from 1 without gaps or duplicates.
#'   Its columns are sampled in order and reused cyclically until all `NS`
#'   movements are assigned; `NS` remains the total per pair and demand.
#' @return A list with the keyed data.tables `DnDraws` (one row per demand,
#'   pair, model, realization and yield acceleration, with the integer draw
#'   identifier `realization` and the sampled PGA `PGA.s`) and `DnTable`
#'   (quantiles and mean of every model and family). Stored draws whose
#'   identifier is named `s` require that column to be renamed to
#'   `realization` before passing them to the draw consumers.
#'   The rows of the flexible models (`BT07`, `BM17`, `BM19`) are conditional
#'   on a nonzero displacement, as [evaluateFlexibleDisplacement()]
#'   documents; the probability of zero displacement is not represented.
#'   The inputs are not modified.
#' @export
#' @examples
#' \donttest{
#' # The inputs of the product test of newmark 1.7.1: one rock lognormal per
#' # period, its MCE quantiles and two geometry-material pairs.
#' Rock <- data.table::data.table(ID = "TEST", siteID = "TST", TR = 100, Tn = c(0, 0.2, 1),
#'                                Sa = c(0.3, 0.45, 0.24), Vs30 = 760, Vref = 760)
#' Rock[, `:=`(SaF = Sa, AF = NA_real_, k = "oq", muLn = log(Sa) - 0.4^2 / 2, sigmaLn = 0.4,
#'             w = 1, fitErr = 0)]
#' Levels <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
#' MCE <- Rock[, .(p = c(sprintf("%.2f", Levels), "mean"),
#'                 Sa = c(exp(muLn + 0.4 * qnorm(Levels)), Sa)), by = Tn]
#' MCE[, `:=`(ID = "MCE", siteID = "TST", Vs30 = 760, Vref = 760)]
#' Shear <- data.table::data.table(IDg = c("G1", "G2"), IDm = c("U1", "U2"), VS30 = c(360, 560),
#'                                 Ts = 0.3, Vref = 760, level = "mean")
#' Draws <- Shear[, .(IDg, IDm, realization = 1L, Ts, VS30, Vref)]
#' Site <- buildDnTables(Rock, MCE, Shear, modelID = "TEST", siteID = "TST", TR = 100, Vref = 760,
#'                       Mw = 7, NS = 2, weights = c(AM88 = 1, BM19 = 1), zTrunc = 3,
#'                       workers = 1, draws = Draws)
#' Site$DnTable[IDn == "ensemble" & p == "0.50" & ky < 0.02]
#' }
buildDnTables <- function(rock, MCE, shear, modelID, siteID, TR, Vref, Mw, NS, weights, zTrunc, workers, draws) {
  .checkDnArguments(modelID = modelID, siteID = siteID, TR = TR, Vref = Vref, NS = NS, weights = weights)
  if (!.isMagnitude(Mw, TR)) stop("Mw must be one finite positive number or one per demand (each return period of TR and MCE)")
  if (!is.numeric(zTrunc) || length(zTrunc) != 1L || !is.finite(zTrunc) || zTrunc <= 0) {
    stop("zTrunc must be a finite positive scalar")
  }
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) || workers < 1 ||
      workers != as.integer(workers)) {
    stop("workers must be one positive integer")
  }
  Family <- sort(unique(as.numeric(TR)))
  Shear <- .checkDnShear(shear, Vref)
  Profiles <- .checkDnProfiles(draws, shear = Shear, Vref = Vref)[realization <= NS]
  Rock <- .checkDnRock(rock, modelID, siteID, Vref, Family, 1.5 * max(Profiles$Ts))
  Mce <- .checkDnMCE(MCE, siteID, Vref, 1.5 * max(Profiles$Ts))
  LIST <- lapply(c(Family, NA_real_), function(tr) {
    Demand <- if (is.na(tr)) "MCE" else modelID
    Records <- .amplifyDemand(if (is.na(tr)) Mce else Rock[TR == tr, .(Tn, k, muLn, sigmaLn, w)],
                              shear = Shear, draws = Profiles, Vref = Vref, zTrunc = if (is.na(tr)) NA_real_ else zTrunc,
                              demandID = Demand, TR = tr)
    .simulateDnDemand(Records, demandID = Demand, TR = tr, siteID = siteID, Mw = .selectMagnitude(Mw, tr), NS = as.integer(NS),
                      weights = weights, seed = .seedDemand(siteID, Demand, tr), workers = as.integer(workers))
  })
  OUT <- list(DnDraws = rbindlist(lapply(LIST, `[[`, "DnDraws")),
              DnTable = rbindlist(lapply(LIST, `[[`, "DnTable")))
  .checkDnTable(OUT$DnTable, Shear, siteID, weights, modelID, Family)
  .checkDnDraws(OUT$DnDraws, Shear, siteID, weights, NS, modelID, Family)
  OUT
}

# The column spectra and common mean-profile ky grid of each pair. Site
# response depends on velocity only and is computed once per distinct VS30.
.amplifyDemand <- function(rockSpectrum, shear, draws, Vref, zTrunc, demandID, TR) {
  IsAnalytic <- !is.na(TR)
  Velocities <- unique(c(shear$VS30, round(draws$VS30)))
  Spectra <- lapply(Velocities, function(v) {
    if (IsAnalytic) {
      return(refineSiteMixture(rockSpectrum, VS30 = v, Vref = Vref, model = "ST20", levels = DNLEVELS,
                               orders = c(16L, 32L, 64L, 128L, 256L), tolerance = 2e-3)$mixture)
    }
    refineSiteSpectrum(rockSpectrum, VS30 = v, Vref = Vref, model = "ST20", orders = c(32L, 64L, 128L, 256L),
                       tolerance = 2e-3)$spectrum[, .(Tn = as.numeric(Tn), p = as.character(p), Sa = as.numeric(SaF))]
  })
  LIST <- lapply(seq_len(nrow(shear)), function(i) {
    Spectrum <- Spectra[[match(shear$VS30[[i]], Velocities)]]
    Cell <- paste0("ID=", demandID, ", TR=", format(TR, scientific = FALSE), ", IDg=", shear$IDg[[i]],
                   ", IDm=", shear$IDm[[i]])
    .checkAmplified(Spectrum, IsAnalytic, Cell)
    PGA <- if (IsAnalytic) .checkMixture(Spectrum)[Tn == 0, Reduce(`+`, w * exp(muLn + sigmaLn^2 / 2))] else Spectrum[Tn == min(Tn) & p == "mean", Sa]
    if (length(PGA) != 1L) stop("site response has no single mean PGA for ", Cell)
    Ky <- buildKyGrid(PGA = PGA, n = 30L)
    if (anyDuplicated(Ky)) stop("the yield grid has repeated values for ", Cell)
    Profiles <- draws[.(shear$IDg[[i]], shear$IDm[[i]]), on = .(IDg, IDm)]
    Profiles[, spectrum := Spectra[match(round(VS30), Velocities)]]
    for (Spectrum in Profiles$spectrum) .checkAmplified(Spectrum, IsAnalytic, Cell)
    data.table(IDg = shear$IDg[[i]], IDm = shear$IDm[[i]], Vs30 = shear$VS30[[i]], Ts = shear$Ts[[i]],
               zTrunc = zTrunc, ky = list(Ky), profiles = list(Profiles))
  })
  rbindlist(LIST)
}

# A site spectrum the displacement sampling can use: a valid mixture, or the
# complete quantile family of the product levels with positive ordinates.
.checkAmplified <- function(spectrum, analytic, cell) {
  if (analytic) {
    if (!setequal(names(spectrum), c("Tn", "k", "muLn", "sigmaLn", "w"))) stop("site response has an invalid analytic schema for ", cell)
    .checkMixture(spectrum)
    return(invisible(TRUE))
  }
  if (!setequal(names(spectrum), c("Tn", "p", "Sa"))) stop("site response has an invalid quantile schema for ", cell)
  if (anyDuplicated(spectrum, by = c("Tn", "p")) || !setequal(unique(spectrum$p), DNLABELS) ||
      nrow(spectrum) != uniqueN(spectrum$Tn) * length(DNLABELS)) {
    stop("site response returned an incomplete Tn x p family for ", cell)
  }
  if (any(!is.finite(spectrum$Sa)) || any(spectrum$Sa <= 0)) stop("site response returned non-finite or non-positive Sa for ", cell)
  invisible(TRUE)
}

# The fixed seed of one site and demand: the base-31 hash of
# "<siteID>/<demandID>/<TR>" (no TR for the MCE) modulo 2^31 - 1.
.seedDemand <- function(siteID, demandID, TR) {
  Label <- if (is.na(TR)) character() else format(TR, scientific = FALSE, trim = TRUE)
  Seed <- 0
  for (x in utf8ToInt(paste(c(siteID, demandID, Label), collapse = "/"))) {
    Seed <- (Seed * 31 + x) %% 2147483647
  }
  as.integer(Seed)
}

# Movements cycle through each pair's column records. All random draws
# precede evaluation, which may run in forked workers without changing them.
.simulateDnDemand <- function(records, demandID, TR, siteID, Mw, NS, weights, seed, workers) {
  # Local names, so that no column of the records can shadow them in j.
  Demand <- demandID
  Period <- TR
  Site <- siteID
  set.seed(seed)
  Draws <- lapply(seq_len(nrow(records)), function(i) {
    Intensities <- .sampleDnIntensities(records$profiles[[i]], NS = NS, zTrunc = records$zTrunc[[i]])
    list(intensities = Intensities, eps = stats::rnorm(NS))
  })
  evaluate <- function(i) {
    Draw <- Draws[[i]]
    Profiles <- records$profiles[[i]]
    Realizations <- rbindlist(lapply(seq_len(nrow(Profiles)), function(j) {
      Intensities <- Draw$intensities[seq.int(j, NS, by = nrow(Profiles))]
      DT <- evaluateDisplacements(Intensities, eps = Draw$eps[Intensities$realization], ky = records$ky[[i]],
                                   Ts = Profiles$Ts[[j]], Mw = Mw, models = names(weights),
                                   NFC = "D100", pulse = FALSE)
      DT[, PGA.s := Intensities$PGA[match(realization, Intensities$realization)]]
      DT
    }))
    Levels <- c(.dnLevels(Profiles$spectrum[[1L]]), "mean")
    Families <- lapply(names(DNFAMILIES), function(family) {
      Members <- .familyWeights(weights, family)
      if (!length(Members)) return(NULL)
      .poolDisplacements(Realizations, weights = Members, levels = Levels,
                         estimator = .textWeightedQuantile)[, IDn := family]
    })
    Curve <- rbindlist(c(list(summarizeDisplacements(Realizations, levels = Levels)), Families), use.names = TRUE)
    Curve[, `:=`(IDg = records$IDg[[i]], IDm = records$IDm[[i]], Ts = records$Ts[[i]], Vs30 = records$Vs30[[i]])]
    Realizations[, `:=`(IDg = records$IDg[[i]], IDm = records$IDm[[i]])]
    list(curve = Curve, draws = Realizations)
  }
  OUT <- if (workers > 1L && nrow(records) > 1L && !identical(.Platform$OS.type, "windows")) {
    Threads <- getDTthreads()
    setDTthreads(1L)
    on.exit(setDTthreads(Threads), add = TRUE)
    parallel::mclapply(seq_len(nrow(records)), evaluate, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(records)), evaluate)
  }
  Failed <- vapply(OUT, function(x) inherits(x, "try-error") || inherits(x, "condition"), logical(1L))
  if (any(Failed)) {
    stop("dn evaluation failed: ", paste(unique(vapply(OUT[Failed], function(x) {
      conditionMessage(if (inherits(x, "try-error")) attr(x, "condition") else x)
    }, character(1L))), collapse = "; "))
  }
  Curve <- rbindlist(lapply(OUT, `[[`, "curve"))
  Curve[, `:=`(p = .canonicalLevel(level), ID = Demand, TR = Period, siteID = Site, units.Dn = "cm")]
  Realizations <- rbindlist(lapply(OUT, `[[`, "draws"))
  Realizations[, `:=`(ID = Demand, TR = Period, siteID = Site, units.Dn = "cm")]
  list(DnTable = Curve[, DNCOLS, with = FALSE], DnDraws = Realizations[, DRAWCOLS, with = FALSE])
}

# Complete-cell sampling: the random stream is independent of column grouping,
# and the target mean weights each column by its assigned movements.
.sampleDnIntensities <- function(profiles, NS, zTrunc) {
  Innovations <- matrix(stats::rnorm(3L * NS), ncol = 3L, byrow = TRUE)
  Sample <- Means <- matrix(NA_real_, nrow = NS, ncol = 3L)
  for (j in seq_len(nrow(profiles))) {
    IDX <- seq.int(j, NS, by = nrow(profiles))
    Draw <- .sampleIntensities(profiles$spectrum[[j]], Ts = profiles$Ts[[j]], NS = length(IDX),
                               zTrunc = zTrunc, generator = function(n) as.vector(t(Innovations[IDX, , drop = FALSE])))
    Sample[IDX, ] <- Draw$sample
    Means[IDX, ] <- rep(Draw$mean, each = length(IDX))
  }
  .intensityTable(.centerSpectra(Sample, mean = apply(Means, 2L, .portableMean)))
}

# The active models of one family of DNFAMILIES with their weights, in the
# order of the family; empty when the family has no active model. The
# consumer normalizes them.
.familyWeights <- function(weights, family) {
  weights[intersect(DNFAMILIES[[family]], names(weights))]
}

# Probability levels of a scenario summary: those of the product for an
# analytic spectrum, the quantile labels of a quantile spectrum.
.dnLevels <- function(spectrum) {
  if (.classifySpectrum(spectrum) == "analytic") return(DNLEVELS)
  # Labels that do not convert, such as "mean", are dropped by the
  # finiteness filter; the coercion warning carries no information.
  x <- suppressWarnings(as.numeric(as.character(spectrum$p[spectrum$p != "mean"])))
  sort(unique(x[is.finite(x)]))
}

# Product labels of the summary levels: "mean" or the probability as "%.2f",
# one of the product levels.
.canonicalLevel <- function(x) {
  OUT <- trimws(as.character(x))
  Mean <- !is.na(OUT) & tolower(OUT) == "mean"
  # Unconvertible labels become NA and are rejected below.
  Probability <- suppressWarnings(as.numeric(OUT[!Mean]))
  if (any(!is.finite(Probability) | Probability <= 0 | Probability >= 1)) {
    stop("displacement summary contains an invalid probability label")
  }
  Index <- match(Probability, DNLEVELS)
  if (anyNA(Index)) stop("displacement summary contains an unsupported probability level")
  OUT[Mean] <- "mean"
  OUT[!Mean] <- sprintf("%.2f", DNLEVELS[Index])
  OUT
}

# Product contract of the dn tables of newmark 1.7.1: probability levels,
# their labels, model families, columns, keys and the relative tolerance of
# the reproduction checks.
DNLEVELS <- c(0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95)
DNLABELS <- c("mean", sprintf("%.2f", DNLEVELS))
DNFAMILIES <- list(ensemble = c("AM88", "JB07", "SR08", "BT07", "BM17", "BM19"),
                   rigid = c("AM88", "JB07", "SR08"),
                   flexible = c("BT07", "BM17", "BM19"))
DNCOLS <- c("IDg", "IDm", "p", "ID", "TR", "ky", "Dn", "Ts", "Vs30", "IDn", "siteID", "units.Dn")
DNKEY <- c("siteID", "ID", "TR", "IDg", "IDm", "IDn", "p", "ky")
DRAWCOLS <- c("IDg", "IDm", "ID", "TR", "ky", "realization", "PGA.s", "IDn", "Dn", "siteID", "units.Dn")
DRAWKEY <- c("siteID", "ID", "TR", "IDg", "IDm", "IDn", "realization", "ky")
DNTOL <- 1e-10

.isIdentifier <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[A-Za-z][A-Za-z0-9_-]*$", x)
}

.checkDnArguments <- function(modelID, siteID, TR, Vref, NS, weights) {
  if (!.isIdentifier(modelID) || !.isIdentifier(siteID)) stop("modelID and siteID must be valid identifiers")
  if (!is.numeric(Vref) || length(Vref) != 1L || is.na(Vref) || Vref != 760) stop("Vref must be 760 m/s")
  if (!is.numeric(TR) || !length(TR) || any(!is.finite(TR)) || any(TR <= 0)) {
    stop("TR must contain positive finite return periods")
  }
  if (!is.numeric(NS) || length(NS) != 1L || !is.finite(NS) || NS < 1 || NS != as.integer(NS)) {
    stop("NS must be a positive integer")
  }
  IDs <- names(weights)
  if (!is.numeric(weights) || !length(weights) || is.null(IDs) || anyNA(IDs) || anyDuplicated(IDs) ||
      any(!IDs %in% DNFAMILIES$ensemble) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must name active models with positive finite values")
  }
  invisible(NULL)
}

# The mean rows of the shear table, one per IDg x IDm, ordered.
.checkDnShear <- function(shear, Vref) {
  Required <- c("IDg", "IDm", "VS30", "Ts", "Vref", "level")
  if (!is.data.frame(shear) || !all(Required %in% names(shear))) {
    stop("ShearTable is missing columns: ", paste(setdiff(Required, names(shear)), collapse = ", "))
  }
  Shear <- as.data.table(shear)[level == "mean", .(IDg = as.character(IDg), IDm = as.character(IDm),
                                                   VS30 = as.numeric(VS30), Ts = as.numeric(Ts),
                                                   Vref = as.numeric(Vref))]
  if (!nrow(Shear)) stop("ShearTable has no rows with level == 'mean'")
  if (anyNA(Shear) || any(!nzchar(Shear$IDg)) || any(!nzchar(Shear$IDm))) {
    stop("ShearTable mean rows contain missing or empty identifiers/values")
  }
  if (any(!is.finite(Shear$VS30)) || any(Shear$VS30 <= 0) || any(!is.finite(Shear$Ts)) || any(Shear$Ts <= 0)) {
    stop("ShearTable mean rows require finite positive VS30 and Ts")
  }
  Reference <- Vref
  if (any(Shear$Vref != Reference)) stop("ShearTable mean rows must all have Vref = ", format(Reference))
  if (anyDuplicated(Shear, by = c("IDg", "IDm"))) stop("ShearTable mean rows are not unique by IDg x IDm")
  Shear[, Vref := NULL]
  setorder(Shear, IDg, IDm)
  Shear
}

# Paired column properties in their stored realization order. The movement
# identity cycles over this order without adding a second identity to DnDraws.
.checkDnProfiles <- function(draws, shear, Vref) {
  Required <- c("IDg", "IDm", "realization", "Ts", "VS30", "Vref")
  if (!is.data.frame(draws) || !nrow(draws) || !all(Required %in% names(draws))) {
    stop("ShearDraws must contain IDg, IDm, realization, Ts, VS30 and Vref")
  }
  DT <- as.data.table(draws)[, Required, with = FALSE]
  if (!is.character(DT$IDg) || !is.character(DT$IDm) || anyNA(DT) ||
      !all(vapply(DT[, .(realization, Ts, VS30, Vref)], is.numeric, logical(1L))) ||
      any(!is.finite(DT$Ts) | DT$Ts <= 0) || any(!is.finite(DT$VS30) | DT$VS30 <= 0) ||
      any(!is.finite(DT$realization) | DT$realization < 1 | DT$realization != trunc(DT$realization)) ||
      any(DT$Vref != Vref)) {
    stop("ShearDraws requires paired positive Ts/VS30, positive integer realizations and the reference Vref")
  }
  if (!fsetequal(unique(DT[, .(IDg, IDm)]), shear[, .(IDg, IDm)])) {
    stop("ShearDraws must cover exactly the geometry-material pairs of ShearTable")
  }
  setorder(DT, IDg, IDm, realization)
  if (any(DT[, .(complete = identical(as.numeric(realization), as.numeric(seq_len(.N)))), by = .(IDg, IDm)]$complete == FALSE)) {
    stop("ShearDraws realizations must run from 1 without duplicates or gaps within each pair")
  }
  DT
}

# The analytic rock mixtures of one site and model at Vs30 = Vref, checked
# against the UHSRock lineage, the return periods and the period coverage.
.checkDnRock <- function(rock, modelID, siteID, Vref, TR, requiredTn) {
  COLS <- c("ID", "siteID", "TR", "Tn", "Sa", "Vs30", "Vref", "SaF", "AF", "k", "muLn", "sigmaLn", "w", "fitErr")
  DT <- copy(as.data.table(rock))
  if (!identical(names(DT), COLS)) stop("UHSRock has an invalid schema")
  if (DT[, any(is.na(ID) | !nzchar(ID) | is.na(siteID) | !nzchar(siteID) | !is.finite(TR) | TR <= 0 |
                 !is.finite(Tn) | Tn < 0 | !is.finite(Sa) | Sa <= 0 | !is.finite(SaF) | SaF <= 0 |
                 abs(SaF - Sa) > 1e-12 * pmax(1, abs(SaF)) | !is.na(AF) | !is.finite(Vs30) | Vs30 != 760 |
                 !is.finite(Vref) | Vref != 760 | is.na(k) | !nzchar(as.character(k)) | !is.finite(muLn) |
                 !is.finite(sigmaLn) | sigmaLn < 0 | !is.finite(w) | w <= 0 | !is.finite(fitErr) | fitErr < 0)]) {
    stop("UHSRock has invalid analytic 760-m/s rock lineage")
  }
  if (anyDuplicated(DT, by = c("ID", "siteID", "TR", "Tn", "k"))) {
    stop("UHSRock is not unique by ID, siteID, TR, Tn, k")
  }
  if (nrow(DT[, .(components = uniqueN(k), weightsValid = all(abs(w - 1) <= 1e-12), hasPGA = any(Tn == 0)),
              by = .(ID, siteID, TR)][components != 1L | !weightsValid | !hasPGA])) {
    stop("UHSRock must contain one period-wise lognormal per spectrum")
  }
  if (DT[, any(!is.finite(abs(exp(muLn + sigmaLn^2 / 2) / Sa - 1)) | abs(exp(muLn + sigmaLn^2 / 2) / Sa - 1) > 1e-12)]) {
    stop("UHSRock analytic parameters do not preserve Sa")
  }
  Model <- modelID
  Site <- siteID
  Reference <- Vref
  Rock <- DT[ID == Model & siteID == Site & Vs30 == Reference & Vref == Reference]
  if (!nrow(Rock)) stop("UHSRock has no rock rows for ID=", modelID, ", siteID=", siteID, ", Vs30=Vref=", Vref)
  if (!setequal(unique(Rock$TR), TR)) {
    stop("UHSRock return periods differ from TR for siteID=", siteID,
         ". UHSRock: ", paste(sort(unique(Rock$TR)), collapse = ", "), ". TR: ", paste(TR, collapse = ", "))
  }
  Coverage <- Rock[, .(hasPGA = any(abs(Tn) <= DNTOL), maxTn = max(Tn)), by = TR]
  if (any(!Coverage$hasPGA)) stop("Rock spectrum has no Tn = 0 PGA row at TR=", paste(Coverage[!(hasPGA), TR], collapse = ", "))
  if (any(Coverage$maxTn + DNTOL < requiredTn)) {
    stop("Rock period grid does not cover required Dn period 1.5*max(Ts) = ", format(requiredTn),
         " at TR=", paste(Coverage[maxTn + DNTOL < requiredTn, TR], collapse = ", "))
  }
  Rock <- Rock[, .(TR, Tn, k, muLn, sigmaLn, w)]
  setorder(Rock, TR, Tn, k)
  Rock
}

# The MCE rock quantiles and mean of one site at Vs30 = Vref: complete,
# monotone and covering the required period.
.checkDnMCE <- function(MCE, siteID, Vref, requiredTn) {
  Required <- c("ID", "siteID", "Tn", "p", "Sa", "Vs30", "Vref")
  if (!is.data.frame(MCE) || !all(Required %in% names(MCE))) {
    stop("MCETable is missing columns: ", paste(setdiff(Required, names(MCE)), collapse = ", "))
  }
  if ("TR" %in% names(MCE)) stop("MCETable must not contain TR")
  Site <- siteID
  Reference <- Vref
  Rock <- as.data.table(MCE)[ID == "MCE" & siteID == Site & Vs30 == Reference & Vref == Reference,
                             .(Tn = as.numeric(Tn), p = as.character(p), Sa = as.numeric(Sa))]
  if (!nrow(Rock)) stop("MCETable has no MCE rock rows for siteID=", siteID, ", Vs30=Vref=", Vref)
  if (anyDuplicated(Rock, by = c("Tn", "p"))) stop("MCE rock rows are duplicated by Tn x p for siteID=", siteID)
  if (!setequal(unique(Rock$p), DNLABELS)) {
    stop("MCE rock p levels differ from the required family: ", paste(DNLABELS, collapse = ", "))
  }
  if (uniqueN(Rock[, .(n = uniqueN(Tn)), by = p]$n) != 1L || nrow(Rock) != uniqueN(Rock$Tn) * length(DNLABELS)) {
    stop("MCE rock p levels do not share one complete Tn grid")
  }
  if (any(!is.finite(Rock$Tn)) || any(Rock$Tn < 0) || any(!is.finite(Rock$Sa)) || any(Rock$Sa <= 0)) {
    stop("MCE rock rows require finite non-negative Tn and positive Sa")
  }
  if (!any(abs(Rock$Tn) <= DNTOL & Rock$p == "mean")) stop("MCE rock spectrum has no Tn = 0, p = mean PGA row")
  if (max(Rock$Tn) + DNTOL < requiredTn) {
    stop("MCE rock period grid does not cover required Dn period 1.5*max(Ts) = ", format(requiredTn))
  }
  Descending <- Rock[p != "mean"][order(as.numeric(p)), .(bad = any(diff(Sa) < -DNTOL)), by = Tn]
  if (any(Descending$bad)) {
    stop("MCE rock Sa quantiles are non-monotone at ", sum(Descending$bad), " Tn values for siteID=", siteID)
  }
  setorder(Rock, Tn, p)
  Rock
}

# Whether every observed (by, value) combination of a product reproduces the
# reference values within a relative 1e-10.
.reproducesDn <- function(product, reference, by, cols) {
  Observed <- unique(product[, c(by, names(cols)), with = FALSE])
  if (nrow(Observed) != nrow(reference) || anyDuplicated(Observed, by = by)) return(FALSE)
  Expected <- reference[, c(by, unname(cols)), with = FALSE]
  setnames(Expected, unname(cols), paste0(names(cols), ".reference"))
  AUX <- merge(Observed, Expected, by = by, sort = FALSE)
  nrow(AUX) == nrow(reference) && all(vapply(names(cols), function(x) {
    a <- AUX[[x]]
    b <- AUX[[paste0(x, ".reference")]]
    all(is.finite(a)) && all(is.finite(b)) && all(abs(a - b) <= DNTOL * pmax(1, abs(a), abs(b)))
  }, logical(1L)))
}

# The demands of a product: the model with exactly `TR` and the MCE
# with a missing return period.
.hasDemands <- function(DT, modelID, TR) {
  setequal(unique(DT$ID), c(modelID, "MCE")) && setequal(unique(DT$TR[DT$ID == modelID]), TR) &&
    all(is.na(DT$TR[DT$ID == "MCE"]))
}

.checkDnTable <- function(DT, shear, siteID, weights, modelID, TR) {
  if (!identical(names(DT), DNCOLS)) {
    stop("Invalid DnTable schema for siteID=", siteID, ". Expected: ", paste(DNCOLS, collapse = ", "),
         ". Found: ", paste(names(DT), collapse = ", "))
  }
  if (!identical(unique(DT$siteID), siteID) || !.hasDemands(DT, modelID, TR) || !setequal(unique(DT$p), DNLABELS)) {
    stop("DnTable has invalid siteID, ID, TR, or p vocabulary for siteID=", siteID)
  }
  Families <- names(DNFAMILIES)[vapply(names(DNFAMILIES), function(x) length(.familyWeights(weights, x)) > 0L, logical(1L))]
  ExpectedIDn <- c(names(weights), Families)
  if (!setequal(unique(DT$IDn), ExpectedIDn)) stop("DnTable has an unexpected IDn family for siteID=", siteID)
  if (anyDuplicated(DT, by = DNKEY)) stop("DnTable key is duplicated for siteID=", siteID)
  if (any(!is.finite(DT$ky)) || any(DT$ky <= 0) || any(!is.finite(DT$Dn)) || any(DT$Dn < 0) ||
      any(!is.finite(DT$Ts)) || any(DT$Ts <= 0) || any(!is.finite(DT$Vs30)) || any(DT$Vs30 <= 0) ||
      any(DT$units.Dn != "cm")) {
    stop("DnTable has invalid numerical values or units for siteID=", siteID)
  }
  Groups <- DT[, .(n = .N, nky = uniqueN(ky)), by = .(ID, TR, IDg, IDm, IDn, p)]
  if (any(Groups$n != 30L) || any(Groups$nky != 30L)) {
    stop("DnTable does not have 30 unique ky per curve for siteID=", siteID)
  }
  Demands <- length(TR) + 1L
  Expected <- Demands * nrow(shear) * length(ExpectedIDn) * length(DNLABELS) * 30L
  if (nrow(DT) != Expected) {
    stop("DnTable row count mismatch for siteID=", siteID, ": expected ", Expected, ", found ", nrow(DT))
  }
  if (nrow(unique(DT[, .(ID, TR, IDg, IDm)])) != Demands * nrow(shear)) {
    stop("DnTable does not cover every demand x IDg x IDm cell")
  }
  if (!.reproducesDn(DT, shear, by = c("IDg", "IDm"), cols = c(Ts = "Ts", Vs30 = "VS30"))) {
    stop("DnTable Ts/Vs30 does not reproduce ShearTable for siteID=", siteID)
  }
  setkeyv(DT, DNKEY)
  invisible(DT)
}

.checkDnDraws <- function(DT, shear, siteID, weights, NS, modelID, TR) {
  if (!identical(names(DT), DRAWCOLS)) stop("Invalid DnDraws schema for siteID=", siteID)
  if (!identical(unique(DT$siteID), siteID) || !.hasDemands(DT, modelID, TR) ||
      !setequal(unique(DT$IDn), names(weights))) {
    stop("DnDraws has invalid siteID, ID, TR, or IDn vocabulary")
  }
  if (anyDuplicated(DT, by = DRAWKEY)) stop("DnDraws key is duplicated for siteID=", siteID)
  if (any(!is.finite(DT$ky)) || any(DT$ky <= 0) || any(!is.finite(DT$PGA.s)) || any(DT$PGA.s <= 0) ||
      any(!is.finite(DT$Dn)) || any(DT$Dn < 0) || any(DT$units.Dn != "cm")) {
    stop("DnDraws has invalid numerical values or units")
  }
  if (any(DT[, .(nPGA = uniqueN(PGA.s)), by = .(ID, TR, IDg, IDm, realization)]$nPGA != 1L)) {
    stop("DnDraws PGA.s must be constant within each sampled demand")
  }
  Expected <- (length(TR) + 1L) * nrow(shear) * length(weights) * NS * 30L
  if (nrow(DT) != Expected) {
    stop("DnDraws row count mismatch for siteID=", siteID, ": expected ", Expected, ", found ", nrow(DT))
  }
  Groups <- DT[, .(n = .N, ns = uniqueN(realization), nky = uniqueN(ky)), by = .(ID, TR, IDg, IDm, IDn)]
  if (any(Groups$n != NS * 30L) || any(Groups$ns != NS) || any(Groups$nky != 30L)) {
    stop("DnDraws does not have NS x 30 values per simulation family")
  }
  setkeyv(DT, DRAWKEY)
  invisible(DT)
}

# A demand is named by its return period as it prints, without exponent
# ("10000"), and the MCE, whose TR is missing, by "MCE": the names of a
# magnitude per demand, in Mw and in the Mw.gmdp of a contract.
.nameDemands <- function(TR) {
  vapply(TR, function(x) if (is.na(x)) "MCE" else format(x, scientific = FALSE, trim = TRUE), character(1L))
}

# One finite positive magnitude for every demand, or one for each demand of
# TR and the MCE, named once each.
.isMagnitude <- function(Mw, TR) {
  is.numeric(Mw) && all(is.finite(Mw)) && all(Mw > 0) &&
    (length(Mw) == 1L || (!anyDuplicated(names(Mw)) && setequal(names(Mw), .nameDemands(c(TR, NA)))))
}

# The magnitude of the demand of return period TR (NA for the MCE).
.selectMagnitude <- function(Mw, TR) {
  if (length(Mw) == 1L) return(Mw)
  Mw[[.nameDemands(TR)]]
}
