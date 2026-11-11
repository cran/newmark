# SPDX-FileCopyrightText: 2026 Alejandro Verri Kozlowski <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of newmark: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
#' Peak ground velocity from peak ground acceleration
#'
#' Evaluates the in-house rock regression of `newmark` 1.7.1 that gives each
#' sampled peak ground acceleration its peak ground velocity, the velocity
#' the pulse equations of `BM19` in [evaluateFlexibleDisplacement()] read.
#'
#' @section Mathematics:
#' \deqn{PGV = 100\, PGA^{1.0529}\, e^{0.1241}}
#' with PGA in g and PGV in cm/s, a median relation without dispersion,
#' fitted on the rock records of Verri (2023), as recorded in
#' `newmark` 1.7.1. It is evaluated realization by realization, so a sampled
#' PGV stays paired with its PGA.
#'
#' @param PGA Peak ground accelerations in g, finite and positive.
#' @return Peak ground velocities in cm/s, aligned with `PGA`.
#' @export
#' @examples
#' estimatePGV(PGA = c(0.1, 0.5))
estimatePGV <- function(PGA) {
  if (!is.numeric(PGA) || !length(PGA) || any(!is.finite(PGA)) || any(PGA <= 0)) {
    stop("PGA must be finite positive accelerations")
  }
  (PGA^1.0529) * exp(0.1241) * 100
}

#' @name estimatePGV
#' @references Verri, A. D. (2023). Modelos Parametricos para la
#'   Seleccion de Sismos de Diseno en Infraestructura Minera. Seminario de
#'   Vinculacion y Transferencia (SeVyT) 2023.
NULL
