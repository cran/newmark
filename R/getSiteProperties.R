#' Monte Carlo site properties from synthetic soil profiles
#'
#' Generates `NR` realisations of synthetic soil profiles for a given
#' total height `Hs` and USCS classification, then summarises the
#' resulting Ts, mo, Go, VSo and VS30 by user-requested quantile
#' levels (and/or the mean). Each realisation is built by
#' \code{geSiteTable()}.
#'
#' @param Hs Numeric vector. Total height(s) to hard ground (m).
#' @param USCS Character vector. USCS soil classification codes
#'   (e.g. \code{c("GC", "CL", "ML")}).
#' @param POP Numeric. Pre-consolidation pressure (kPa). Default 100.
#' @param Water Numeric. Water-table depth as fraction of \code{Hs}.
#'   Default 0 (no water).
#' @param NR Integer. Number of Monte Carlo realisations. Default 1.
#' @param h Numeric. Layer thickness (m) in the discretisation.
#'   Default 1.00.
#' @param levels Vector. Quantile levels to report; entries can be
#'   numeric in (0, 1) and/or the literal \code{"mean"}.
#'   Default \code{c(0.05, 0.5, "mean", 0.95)}.
#' @param Vref Numeric. Reference shear-wave velocity (m/s) for
#'   site-class assignment. Default 760.
#'
#' @return data.table with one row per (Hs, POP, Hw, level)
#'   combination and columns USCS, Go, mo, Ts, VSo, VS30, Z500,
#'   Z1000, level. Numeric properties are rounded for readability.
#' @export getSiteProperties
#'
#' @import data.table
#' @importFrom stats predict
#'
#' @examples
#' \donttest{
#' getSiteProperties(
#'   Hs   = 30,
#'   USCS = c("GC", "CL", "ML"),
#'   NR   = 50,
#'   levels = c(0.16, "mean", 0.84)
#' )
#' }
getSiteProperties <- function(Hs, USCS, POP = 100, Water = 0, NR = 1, h = 1.00,
                              levels = c(0.05, 0.5, "mean", 0.95),
                              Vref = 760) {
  . <- NULL

  # Build dataset
  SiteTable <- data.table()


  for (pop in POP) {
    for (hs in Hs) {
      for (w in Water) {
        for (i in 1:NR) {
          if (interactive()) cat(sprintf("\rBuilding sample %d/%d..", i, NR))
          DT <- geSiteTable(Hs = hs, h = h, Water = w, USCS = USCS, POP = pop)
          SiteTable <- rbindlist(list(DT, SiteTable), use.names = TRUE)
        }
      }
    }
  }
  if (interactive()) cat("Done!\n")


  DT <- data.table()
  if (any(levels == "mean")) {
    DT <- SiteTable[, .(
      USCS = paste0(USCS, collapse = " "),
      Go = mean(Go) |> round(1),
      mo = mean(mo) |> round(3),
      Ts = mean(Ts) |> round(2),
      VSo = mean(VSo) |> round(1),
      VS30 = mean(VS30) |> round(1),
      Z500 = mean(Z500) |> round(1),
      Z1000 = mean(Z1000),
      level = "mean"
    ), by = .(Hs, POP, Hw)] |>
      rbind(DT)
  }
  if (any(levels != "mean")) {
    probs <- levels[levels != "mean"] |> as.numeric()
    DT <- SiteTable[, .(
      USCS = paste0(USCS, collapse = " "),
      Go = quantile(Go, probs = probs) |> round(1),
      mo = quantile(mo, probs = probs) |> round(3),
      Ts = quantile(Ts, probs = probs) |> round(2),
      VSo = quantile(VSo, probs = probs) |> round(1),
      VS30 = quantile(VS30, probs = probs) |> round(1),
      Z500 = quantile(Z500, probs = probs) |> round(1),
      Z1000 = quantile(Z1000, probs = probs) |> round(1),
      level = probs
    ), by = .(Hs, POP, Hw)] |>
      rbind(DT)
  }
  DT
}
