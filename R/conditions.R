# SPDX-FileCopyrightText: 2026 Alejandro Verri Kozlowski <averri@fi.uba.ar>
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of newmark: free software under the GNU Affero General Public License, version 3 or later, see <https://www.gnu.org/licenses/>.
# Text and data mining reserved for any use not licensed under the AGPL-3: see <https://kanameishi.github.io/tdm-policy.html>.
# The library signals the classed conditions of the AOM family; the newmark
# command draws them on stderr through its renderer (install/cli/console.R),
# and an R session shows them as ordinary messages, warnings and errors.

# level is moduleNote, moduleProgress or moduleDetail; operation, phase, unit,
# path, count and total are the fields the renderer adds with --verbose. The
# final newline keeps the text a complete line under R's default message
# handler. The condition is built as messageCondition() of R 4.6 builds it,
# because the package runs on R 4.2 and later.
.inform <- function(level, ..., operation, phase = NULL, unit = NULL, path = NULL, count = NULL, total = NULL) {
  message(structure(class = c(level, "message", "condition"),
                    list(message = .makeMessage(..., appendLF = TRUE), call = NULL, operation = operation, phase = phase,
                         unit = unit, path = path, count = count, total = total)))
}

# A recoverable anomaly, with the same fields. Like warning() with
# call. = FALSE, it names no call.
.warn <- function(..., operation, phase = NULL, unit = NULL, path = NULL, count = NULL, total = NULL) {
  warning(warningCondition(.makeMessage(...), class = "moduleWarning", operation = operation, phase = phase,
                           unit = unit, path = path, count = count, total = total))
}

# The classes of a failure after files were published: the family's partial
# publication, a failure of the module, so the renderer writes its message
# alone and the command exits 1.
PARTIALPUBLICATION <- c("newmarkPartialPublication", "modulePartialPublication", "moduleFailure")
