# Layers of one column ordered from the top down and validated for the
# numerical periods: positive thickness, unit weight and stiffness in every
# layer, and distinct positive midpoint depths, which order the layers while
# the thicknesses place their interfaces. The stiffness is the shear-wave
# velocity `VSm` when the table has it, otherwise the velocity that `Gm` gives
# through the relation of evaluateSoilLayers(), sqrt(9.81 * 1000 * Gm / gsat).
.orderLayers <- function(layers) {
  if (!is.data.frame(layers) || !nrow(layers) || !all(c("hs", "zm", "gsat") %in% names(layers)) ||
      !any(c("VSm", "Gm") %in% names(layers))) {
    stop("layers must contain rows with columns hs, zm, gsat and VSm or Gm")
  }
  Stiffness <- intersect(c("VSm", "Gm"), names(layers))[1L]
  for (Column in c("hs", "zm", "gsat", Stiffness)) {
    if (!is.numeric(layers[[Column]]) || any(!is.finite(layers[[Column]])) || any(layers[[Column]] <= 0)) {
      stop("layers$", Column, " must be finite and positive in every layer")
    }
  }
  if (anyDuplicated(layers$zm)) stop("layers$zm must be distinct, since it orders the layers from the top down")
  IDX <- order(layers$zm)
  VSm <- switch(Stiffness, VSm = layers$VSm, Gm = sqrt(9.81 * 1000 * layers$Gm / layers$gsat))
  list(hs = layers$hs[IDX], zm = layers$zm[IDX], gsat = layers$gsat[IDX], VSm = VSm[IDX])
}

# Fundamental period of the layered shear beam and its mesh-doubling control:
# the Galerkin problem with `elements` linear elements per layer and again
# with every element halved, each solved by the Cholesky reduction on the
# stiffness matrix documented in solveColumnPeriod(). `width` gives the width
# of the beam, up to a constant factor, at depths below its top; the top is
# free and the node at the bottom of the last layer is fixed.
.solveShearBeam <- function(Layers, elements, width) {
  if (!is.numeric(elements) || length(elements) != 1L || !is.finite(elements) || elements < 1 ||
      elements != round(elements)) {
    stop("elements must be one positive whole number of elements per layer")
  }
  period <- function(count) {
    Layer <- rep(seq_along(Layers$hs), each = count)
    Length <- Layers$hs[Layer] / count
    Width <- width(c(0, cumsum(Length)))
    n <- length(Length)
    Wa <- Width[-(n + 1L)]
    Wb <- Width[-1L]
    # Unit weights stand for the densities and unit weight times squared
    # velocity for the moduli: the common factor 1/g cancels in omega^2.
    Rho <- Layers$gsat[Layer]
    k <- Rho * Layers$VSm[Layer]^2 * (Wa + Wb) / (2 * Length)
    Md <- (c(Rho * Length * (3 * Wa + Wb), 0) + c(0, Rho * Length * (Wa + 3 * Wb)))[seq_len(n)] / 12
    Mo <- (Rho * Length * (Wa + Wb) / 12)[-n]
    # Mass from the top to each free node, c_i, and the matrix S of the reduction.
    Mass <- cumsum(Md + c(0, Mo) + c(Mo, 0))
    S <- outer(Mass, Mass, pmin)
    diag(S) <- Mass - c(Mo, 0)
    2 * pi * sqrt(eigen(S / sqrt(outer(k, k)), symmetric = TRUE, only.values = TRUE)$values[1L])
  }
  Coarse <- period(elements)
  Fine <- period(2 * elements)
  c(period = Fine, meshChange = 1 - Coarse / Fine)
}
