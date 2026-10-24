#' Evaluate stress, stiffness and shear-wave velocity of ordered layers
#'
#' Deterministic physics of a layered column whose index properties are
#' supplied: effective unit weight under the water table, octahedral mean
#' effective stress, overconsolidation ratio, small-strain shear modulus from
#' the correlations of the soil group and the shear-wave velocity of each
#' layer. The layers may come from [sampleSoilLayers()] or from measured or
#' specified data with the same columns.
#'
#' @section Mathematics:
#' With \eqn{K_o = 0.5} and the effective unit weight \eqn{\gamma'} of each
#' layer, assigned by its midpoint depth \eqn{z_m}: \eqn{\max(\gamma_{sat} -
#' 10, 0)} kN/m3 when \eqn{z_m} lies below the water table \eqn{z_w = H_s (1 -
#' w)}, else \eqn{\gamma_{sat}}; the octahedral stress increases by
#' \eqn{\Delta p = (1 + 2K_o) \gamma' h / 3} across each layer; the layer mean
#' \eqn{p_m} is the average of its top and bottom values, and
#' \eqn{OCR = \max((p_m + POP)/p_m, 1)} with the preconsolidation `POP` in kPa.
#' For every correlation `j` of the layer's soil group in
#' [ShearModelParameters] the modulus at the mean stress is
#' \deqn{G_j = A_j F_j(e_o)\, OCR^{m_1}\, 100^{n_j} (p_m/100)^{n_j},}
#' with \eqn{F_j(e) = (C_e - e)^2/(1+e)}, except the linear
#' \eqn{F(e) = C_e - e/(1+e)} of `SS75`, and \eqn{m_1} the exponent \eqn{K}
#' of the overconsolidation ratio of Hardin and Drnevich (1972) for fine
#' soils and zero otherwise (hydrated copy, equation `eq-design-gmax`,
#' printed page 671). \eqn{m_1} is interpolated linearly in `IP` through the
#' nodes of `m1` and is `NA` outside them, so a fine layer whose `IP` lies
#' below the first node or above the last has an `NA` modulus whenever
#' \eqn{OCR > 1}. The default nodes are Table 1 of that paper (printed page
#' 672): `IP` 0, 20, 40, 60 and 80 with 0, 0.18, 0.30, 0.41 and 0.48, and
#' the row "PI \eqn{\ge} 100: 0.50", encoded as the nodes (100, 0.50) and
#' (`Inf`, 0.50) so that \eqn{m_1} is 0.50 at and above 100. The `Inf` node
#' is that encoding of the published row, not a further choice. The layer
#' modulus is the arithmetic mean \eqn{G_m = \overline{G_j}} over the group's
#' correlations, in MPa, and
#' \deqn{V_{S} = \sqrt{9.81 \cdot 1000\, G_m / \gamma_{sat}}}
#' in m/s uses the saturated unit weight. These are the relations of `newmark`
#' 1.7.1 evaluated layer by layer, except the default nodes: 1.7.1 ended
#' them at (100, 0.48), with `NA` above, and [buildShearTable()] passes
#' those nodes to reproduce it. The spread of \eqn{G_j} across correlations
#' is not retained by this version.
#'
#' @param layers A data.frame or data.table with one row per layer from the
#'   surface down and columns `USCS`, `hs` (m), `zm` (m), `eo`, `gsat`
#'   (kN/m3) and `IP` (percent, `NA` for coarse soils); `emin` and `emax`
#'   are used for the relative density when present.
#' @param Hs One positive column height in metres.
#' @param water Water-table height as a fraction of `Hs` (0 for a dry column).
#' @param POP Preconsolidation pressure in kPa, added to the mean stress in
#'   the overconsolidation ratio.
#' @param m1 Nodes of the interpolation of \eqn{m_1} in `IP`: a data.frame
#'   or data.table with numeric columns `IP` (percent) and `m1`. The default
#'   is Table 1 of Hardin and Drnevich (1972): `IP` 0, 20, 40, 60, 80, 100
#'   and `Inf` with 0, 0.18, 0.30, 0.41, 0.48, 0.50 and 0.50.
#' @return A data.table with the input columns and `group` (`Gravels`,
#'   `Sands` or `Fines`), `Dr` (relative density, when the void-ratio limits
#'   are supplied), `gs` (effective unit weight, kN/m3), `pm` (kPa), `OCR`,
#'   `Gm` (MPa) and `VSm` (m/s). Each row is a conditional evaluation of its
#'   inputs; it carries no dispersion of its own.
#' @export
#' @examples
#' Layers <- data.frame(USCS = c("SM", "ML"), hs = c(5, 5), zm = c(2.5, 7.5),
#'                      eo = c(0.7, 0.9), gsat = c(20, 18), IP = c(NA, 10))
#' evaluateSoilLayers(Layers, Hs = 10, water = 0, POP = 100)
evaluateSoilLayers <- function(layers, Hs, water, POP, m1 = HardinExponent) {
  Required <- c("USCS", "hs", "zm", "eo", "gsat", "IP")
  if (!is.data.frame(layers) || !nrow(layers) || !all(Required %in% names(layers))) {
    stop("layers must contain rows with columns ", paste(Required, collapse = ", "))
  }
  if (!is.numeric(Hs) || length(Hs) != 1L || !is.finite(Hs) || Hs <= 0) {
    stop("Hs must be one finite positive height in metres")
  }
  if (!is.numeric(water) || length(water) != 1L || !is.finite(water) || water < 0 || water > 1) {
    stop("water must be one fraction of Hs between 0 and 1")
  }
  if (!is.numeric(POP) || length(POP) != 1L || !is.finite(POP) || POP < 0) {
    stop("POP must be one finite non-negative pressure in kPa")
  }
  if (!is.data.frame(m1) || !all(c("IP", "m1") %in% names(m1)) || !is.numeric(m1$IP) ||
      !is.numeric(m1$m1) || anyNA(m1$IP) || anyNA(m1$m1)) {
    stop("m1 must be a table of numeric IP and m1 nodes without NA")
  }
  OUT <- as.data.table(layers)
  if (anyNA(OUT$USCS) || !all(OUT$USCS %in% ValidUSCS)) {
    stop("layers$USCS must contain USCS codes among ", paste(ValidUSCS, collapse = ", "))
  }
  for (Column in c("hs", "zm", "eo", "gsat")) {
    if (!is.numeric(OUT[[Column]]) || any(!is.finite(OUT[[Column]])) || any(OUT[[Column]] <= 0)) {
      stop("layers$", Column, " must be finite and positive")
    }
  }
  Group <- .soilGroup(OUT$USCS)
  if (any(Group == "Fines" & !is.finite(OUT$IP))) stop("layers$IP must be finite for fine soils")
  NL <- nrow(OUT)
  zw <- Hs - Hs * water
  gs <- ifelse(OUT$zm > zw, pmax(0, OUT$gsat - 10), OUT$gsat)
  Ko <- 0.5
  # The octahedral stress accumulates from the surface, one double-precision
  # addition per layer, as the layer loop of newmark 1.7.1 did.
  Dp <- 1 / 3 * (1 + 2 * Ko) * gs * OUT$hs
  StressBottom <- Reduce(`+`, Dp, accumulate = TRUE)
  StressTop <- c(0, StressBottom[-NL])
  pm <- (StressBottom + StressTop) / 2
  OCR <- pmax((pm + POP) / pm, 1)
  Exponent <- numeric(NL)
  Fines <- Group == "Fines"
  if (any(Fines)) {
    Exponent[Fines] <- approx(x = m1$IP, y = m1$m1, xout = OUT$IP[Fines])$y
  }
  # Every correlation of the layer's group is evaluated on one row; the layer
  # modulus is the base-R mean over its group, in the order of the parameter
  # table, so the value equals the per-layer loop of newmark 1.7.1 bit for bit.
  Layer <- data.table(layer = seq_len(NL), group = Group, eo = OUT$eo, OCR = OCR, pm = pm, m1 = Exponent)
  Pairs <- .dataset("ShearModelParameters")[Layer, on = .(GroupID = group), allow.cartesian = TRUE]
  Pairs[, Fe := ifelse(ModelID == "SS75", Ce - eo / (1 + eo), (Ce - eo)^2 / (1 + eo))]
  Pairs[, G := A * Fe * (OCR^m1) * 100^(n) * (pm / 100)^(n)]
  Gm <- vapply(split(Pairs$G, Pairs$layer), mean, numeric(1L), USE.NAMES = FALSE)
  OUT[, group := Group]
  if (all(c("emin", "emax") %in% names(OUT))) OUT[, Dr := (emax - eo) / (emax - emin)]
  OUT[, `:=`(gs = gs, pm = pm, OCR = OCR, Gm = Gm, VSm = sqrt(9.81 * Gm * 1000 / gsat))]
  OUT[]
}

# Table 1 of Hardin and Drnevich (1972), printed page 672: the default nodes
# of m1 for the general route. Its last row, "PI >= 100: 0.50", is the pair of
# nodes (100, 0.50) and (Inf, 0.50) under approx's NA-outside rule.
HardinExponent <- data.frame(IP = c(0, 20, 40, 60, 80, 100, Inf),
                             m1 = c(0, 0.18, 0.30, 0.41, 0.48, 0.50, 0.50))

# Soil group of each USCS code: the three families of the stiffness table.
.soilGroup <- function(USCS) {
  OUT <- rep(NA_character_, length(USCS))
  OUT[USCS %in% ValidGravels] <- "Gravels"
  OUT[USCS %in% ValidSands] <- "Sands"
  OUT[USCS %in% ValidFines] <- "Fines"
  OUT
}
