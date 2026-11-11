# buildShearTable() of newmark 1.11.0 (lib/R/buildShearTable.R, lines 65-96, at
# commit 7199e40), copied without its documentation: each pair drawn and
# summarized before the next one is drawn. test-buildShearTable.R sources it
# into an environment whose parent is the package namespace, where it finds
# the helpers that it calls.
buildShearTable <- function(geometry, USCS, IDm, Vref, levels, pairs = NULL) {
  Problem <- .diagnoseShearArguments(geometry, USCS = USCS, IDm = IDm, Vref = Vref, levels = levels, pairs = pairs)
  if (!is.null(Problem)) stop(Problem)
  # Resolved in the function scope, so a USCS column of the caller's table
  # cannot shadow the material pools inside the grouped evaluation.
  materialPool <- function(id) USCS[[id]]
  Pairs <- CJ(IDm, IDg = geometry$IDg)
  if (!is.null(pairs)) Pairs <- as.data.table(pairs)[, .(IDm, IDg)]
  # Quantile rows precede the mean row, as the original summary did.
  Levels <- c(levels[levels != "mean"], if (any(levels == "mean")) "mean")
  MaterialTable <- as.data.table(geometry)[Pairs, on = .(IDg), allow.cartesian = TRUE]
  set(MaterialTable, j = "Vref", value = Vref)
  ShearTable <- MaterialTable[, {
    Profiles <- .legacyProfile(simulateSoilProfiles(Hs = Hs, USCS = materialPool(IDm), NS = 100, h = 1,
                                                    water = 0, POP = 100, Vref = 760,
                                                    m1 = data.frame(IP = c(0, 20, 40, 60, 80, 100),
                                                                    m1 = c(0, 0.18, 0.30, 0.41, 0.48, 0.48))))
    OUT <- summarizeProfiles(Profiles[rev(seq_len(nrow(Profiles)))], levels = Levels,
                             properties = c("Go", "mo", "Ts", "VSo", "VS30", "Z500", "Z1000"))
    OUT[, `:=`(Go = round(Go, 1), mo = round(mo, 3), Ts = round(Ts, 2), VSo = round(VSo, 1),
               VS30 = round(VS30, 1), Z500 = round(Z500, 1))]
    OUT[level != "mean", Z1000 := round(Z1000, 1)]
    data.table(POP = 100, Hw = 0, USCS = paste0(materialPool(IDm), collapse = " "),
               OUT[, .(Go, mo, Ts, VSo, VS30, Z500, Z1000, level)])
  }, by = .(Hs, b, s, Vref, IDm, IDg)]
  ShearTable[, lo := round(b / (2 * Hs * s + b), 3)]
  Roots <- evaluateRoot(buildRootModel(mode = 1, method = "nlm"), mo = ShearTable$mo, lo = ShearTable$lo)
  Roots[!Roots$withinSupport, an := evaluateRoot(buildRootModel(mode = 1, method = "exact"), mo = mo, lo = lo)$an]
  ShearTable[, an := round(Roots$an, 4)]
  ShearTable[, Ts := round(evaluateBeamPeriod(Hs = Hs, lo = lo, mo = mo, VSo = VSo, an = an), 3)]
  ShearTable[]
}
