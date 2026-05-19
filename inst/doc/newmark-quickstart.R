## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----load---------------------------------------------------------------------
library(newmark)
library(data.table)

uhs <- fread(system.file("extdata", "uhs.csv", package = "newmark"))
uhs[p %in% c("0.16", "mean", "0.84") & Tn <= 0.2, .(Tn, p, Sa)]

## ----params-------------------------------------------------------------------
# Ts: fundamental period of the sliding mass (s).
# In production, derived from getSiteProperties() + getCylinderRoots()
# using the soil USCS profile and slope geometry (Ishihara 1996,
# Gazetas & Dakoulas 1985). Here set to a representative value.
Ts <- 0.60

# Mw: moment magnitude from PSHA disaggregation.
Mw <- 6.8

# Ensemble weights (0 = model inactive).
weights <- c(AM88 = 1, JB07 = 0, BT07 = 1, SR08 = 1, BM17 = 0, BM19 = 1)

# Displacement targets (cm).
Da <- c(0.5, 2.5, 5.0, 25.0)

## ----fitDnCurve---------------------------------------------------------------
ky     <- getDnKy(uhs, Ts = Ts)
result <- fitDnCurve(
  uhs     = uhs,
  ky      = ky,
  Ts      = Ts,
  Mw      = Mw,
  NS      = 200,
  weights = weights
)

result$curve[IDn == "ensemble" & p == "mean", .(ky, Dn)]

## ----invertDnDraws------------------------------------------------------------
kmax <- invertDnDraws(result$draws, Da = Da, weights = weights)
kmax[p %in% c("0.16", "mean", "0.84")]

