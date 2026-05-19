# nolint start
#' Single-realisation synthetic soil profile and site properties
#'
#' Builds one Monte Carlo realisation of a layered soil column from a
#' total height \code{Hs} and a USCS classification, samples random
#' geotechnical properties (void ratio, unit weight, plasticity,
#' shear modulus parameters), computes the fundamental period
#' \code{Ts} via Rayleigh's method (\code{fitModel.Ts}) and the
#' Vs30-equivalent shear-wave velocity. \code{getSiteProperties}
#' calls this function repeatedly to assemble the Monte Carlo
#' ensemble; call directly when one realisation suffices.
#'
#' @param Hs Numeric scalar. Total height to hard ground (m).
#' @param Water Numeric scalar. Water-table depth as fraction of
#'   \code{Hs}: \code{zw = Hs - Water * Hs}.
#' @param USCS Character vector. Unified Soil Classification System
#'   codes (e.g. \code{c("GC", "CL", "ML")}).
#' @param Group Character vector or \code{NULL}. Soil groups,
#'   \code{c("Gravels","Fines","Sands")}; if \code{NULL}, derived from
#'   \code{USCS}.
#' @param h Numeric scalar. Layer thickness in the discretisation
#'   (m). Default 0.50.
#' @param DrID Character or \code{NULL}. Relative-density label,
#'   \code{c("Very Dense","Dense","Compact","Loose","Very Loose")};
#'   \code{NULL} samples relative density from the standard prior.
#' @param UniformDistribution Logical. If \code{TRUE} (default),
#'   sample soil properties from uniform distributions; if
#'   \code{FALSE}, use normal distributions.
#' @param POP Numeric. Pre-consolidation pressure (kPa); enters the
#'   over-consolidation ratio \code{OCR = (pm + POP) / pm}.
#'   Default 100.
#' @param IgnoreModelIntervals Logical. If \code{TRUE} (default),
#'   ignore model-specific layer-thickness intervals and use \code{h}
#'   as the discretisation step.
#' @param getSiteLayers Logical. If \code{TRUE}, return the per-layer
#'   detail in addition to the per-site summary.
#' @param Vref Numeric. Reference shear-wave velocity (m/s) for the
#'   site-class assignment. Default 760.
#'
#' @return data.table with one row per (Hs, POP, Hw, sample)
#'   combination and columns Hs, Hw, NL, Z500, Z1000, SID, Go, mo,
#'   Ts, VSo, VS30, plus USCS-fraction columns.
#' @export geSiteTable
#'
#' @import data.table
#' @import digest
#' @import stats
#' @importFrom triangle rtriangle
#' @import stringr
#' @import utils
#'
geSiteTable <- function(Hs, Water = 0, USCS, Group = NULL, h = 0.50, DrID = NULL, Vref = 760, UniformDistribution = TRUE, POP = 100, IgnoreModelIntervals = TRUE, getSiteLayers = FALSE) {
    . <- NULL

    # Ensure the lazy-loaded dataset is bound in this scope. Without
    # this call, bare-name references below resolve via the package
    # namespace under `library(newmark)` but can fail to resolve when
    # the package is loaded via `pkgload::load_all()` (the path pkgdown
    # takes when install = FALSE), producing
    #   "object 'ShearModelParameters' not found".
    utils::data("ShearModelParameters", package = "newmark",
                envir = environment())



    NL <- ceiling(Hs / h)
    # Dimensions  ----
    OV <- vector(mode = "double", length = NL)
    OS <- vector(mode = "character", length = NL)
    # Set vertical coordinates ----
    hs <- rep(h, NL)
    zo <- OV
    zi <- OV
    zm <- OV
    for (k in seq(1, NL)) {
        if (k > 1) {
            zo[k] <- zi[k - 1]
        }
        zi[k] <- zo[k] + h
        zm[k] <- (zi[k] + zo[k]) / 2
    }

    # Set group ID ----
    if (is.null(Group) & is.null(USCS)) {
        USCS <- ValidUSCS
    }

    # Case 2: Sampling soils from the same gropu(s)

    if (!is.null(Group) & is.null(USCS)) {
        OK <- tolower(Group) %in% c("gravels", "sands", "fines", "silts", "clays", "organic")
        stopifnot(OK)
        USCS <- switch(tolower(Group),
            gravels = ValidGravels,
            sands = ValidSands,
            fines = ValidFines,
            silts = ValidSilts,
            clays = ValidClays,
            silts = ValidSilts
        )
    }

    stopifnot(!is.null(USCS) && all(USCS %in% ValidUSCS))
    UID <- sample(USCS, size = NL, replace = TRUE)
    GID <- OV
    gs <- OV # Sumerged,Saturated, depending on zw
    gd <- OV # Dry
    gsat <- OV # Sat
    eo <- OV
    emax <- OV
    emin <- OV
    gsmax <- OV
    gsmin <- OV
    Dr <- OV
    LL <- OV
    IP <- OV
    po <- OV
    pi <- OV
    pm <- OV
    Gm <- OV
    VSm <- OV
    OCR <- OV
    gw <- OV
    # rho <- OV
    Hw <- Hs * Water

    zw <- Hs - Hw

    for (k in seq(1, NL)) {
        # Groupid
        if (UID[k] %in% c(ValidGravels)) {
            GID[k] <- "Gravels"
        }
        if (UID[k] %in% ValidSands) {
            GID[k] <- "Sands"
        }
        if (UID[k] %in% ValidFines) {
            GID[k] <- "Fines"
        }

        # Void Ratios ----
        if (is.null(DrID)) {
            RANGE <- VoidRatiosUSCS[USCS == UID[k]]
        }
        if (!is.null(DrID)) {
            RANGE <- RelativeDensityRanges[toupper(DensityID) == toupper(DrID) & GroupID == GID[k]]
        }

        #
        if (IgnoreModelIntervals == TRUE) {
            emin[k] <- runif(1, min = RANGE$eminMin, max = RANGE$eminMax)
            emax[k] <- runif(1, min = RANGE$emaxMin, max = RANGE$emaxMax)

            if (UniformDistribution == TRUE) {
                eo[k] <- runif(1, min = emin[k], max = emax[k])
            } else {
                eo[k] <- rtriangle(1, a = emin[k], b = emax[k], c = 1 / 2 * (emin[k] + emax[k]))
            }
        }

        if (IgnoreModelIntervals == FALSE) {
            # Build ranges within the model ranges.
            AvailableModels <- FALSE
            NIT <- 0
            NITmax <- 100
            while (AvailableModels == FALSE & NIT < NITmax) {
                NIT <- NIT + 1
                emin[k] <- runif(1, min = RANGE$eminMin, max = RANGE$eminMax)
                emax[k] <- runif(1, min = RANGE$emaxMin, max = RANGE$emaxMax)

                # No risk of emin==emax since ranges were splitted

                if (UniformDistribution == TRUE) {
                    eo[k] <- runif(1, min = emin[k], max = emax[k])
                } else {
                    eo[k] <- rtriangle(1, a = emin[k], b = emax[k], c = 1 / 2 * (emin[k] + emax[k]))
                }

                # Shear Module ----
                MID <- ShearModelParameters[GroupID == GID[k]][eo[k] >= emin & eo[k] <= emax]$ModelID
                if (length(MID) > 0) {
                    AvailableModels <- TRUE
                } else {
                    AvailableModels <- FALSE
                }
            }
            if (NIT >= NITmax) {
                return(NULL)
            }
        }


        # Relative Density
        Dr[k] <- (emax[k] - eo[k]) / (emax[k] - emin[k])
        # Unit Weight kN/m3
        RANGE <- UnitWeightRanges[USCS == UID[k]]

        # **********************************
        # gsmin[k] <- runif(1,min=RANGE$gsminMin,max=RANGE$gsminMax) #[kN/m3]
        # gsmax[k] <- runif(1,min=RANGE$gsmaxMin,max=RANGE$gsmaxMax) #[kN/m3]
        gsmin[k] <- RANGE$gsatMin # [kN/m3]
        gsmax[k] <- RANGE$gsatMax # [kN/m3]
        # **********************************
        # Saturated wunit weight
        # gsat[k] <- gsmax[k] - (gsmax[k] - gsmin[k]) * (1 - Dr[k]) #[kN/m3]
        gsat[k] <- runif(n = 1, min = gsmin[k], max = gsmax[k])

        # Dry weight
        gd[k] <- max(0, gsat[k] - 10)



        # Submerged weight
        if (zm[k] > zw) {
            gw[k] <- 10
            gs[k] <- max(0, gsat[k] - gw[k]) # some organic clays has gs<10
        } else {
            gs[k] <- gsat[k]
        }
        # Density
        # rho[k] <- gs[k]*1000/9.81 # gs in KN/m3, rho in kg/m3


        # Plasticity Index

        if (UID[k] %in% ValidFines) {
            if (UID[k] %in% c("CH", "MH", "OH")) {
                LL[k] <- runif(n = 1, min = 50, max = 100)
            }
            if (UID[k] %in% c("CL", "ML", "OL")) {
                LL[k] <- runif(n = 1, min = 8, max = 50)
            }
            U_LINE <- 0.9 * (LL[k] - 8) # >8 by def LL
            A_LINE <- ifelse(LL[k] >= 20, 0.73 * (LL[k] - 20), 0)
            if (UID[k] %in% ValidClays) {
                IP[k] <- runif(n = 1, min = A_LINE, max = U_LINE)
            }
            if (UID[k] %in% c(ValidSilts, ValidOrganic)) {
                IP[k] <- runif(n = 1, min = 0, max = A_LINE)
            }
        }

        # Octahedral pressures
        Ko <- 0.5
        if (k > 1) {
            po[k] <- pi[k - 1]
        }


        Dp <- 1 / 3 * (1 + 2 * Ko) * gs[k] * hs[k] # [kN/m3->kPa]
        pi[k] <- po[k] + Dp # kPa
        pm[k] <- (pi[k] + po[k]) / 2 # kPa

        # OverConsolidation Ratio
        OCR[k] <- max((pm[k] + POP) / pm[k], 1)


        # Shear Module ----
        MID <- ShearModelParameters[GroupID == GID[k]]$ModelID
        A <- ShearModelParameters[ModelID %in% MID, A] # kPa
        Ce <- ShearModelParameters[ModelID %in% MID, Ce]
        N <- ShearModelParameters[ModelID %in% MID, n]
        Fe <- sapply(Ce, function(x) {
            (x - eo[k])^2 / (1 + eo[k])
        })


        m1 <- ifelse(
            GID[k] == "Fines",
            approx(x = c(0, 20, 40, 60, 80, 100), y = c(0, 0.18, 0.30, 0.41, 0.48, 0.48), xout = IP[k])$y, 0
        )

        # pref <- 100 #[kPa]
        Gref <- sapply(seq(1, length(Fe)), function(n) {
            A[n] * Fe[n] * (OCR[k]^m1) * 100^(N[n])
        }) # [MPa]
        Go <- sapply(seq(1, length(Fe)), function(n) {
            Gref[n] * (pm[k] / 100)^(N[n])
        }) # [MPa]
        Gm[k] <- mean(Go) # [MPa]
        VSm[k] <- sqrt(9.81 * Gm[k] * 1000 / gsat[k]) # sqrt(kPa / (9.81 m/s2 kN/m3)) [m/s]
    }


    # Get properties ----
    VSa <- (cumsum(hs) / cumsum(hs / VSm)) |> round(digits = 1)
    VS30 <- .setVs30(hs = hs, VSm = VSm, VSo = Vref) |> round(digits = 1)
    SID <- Vs30toSID(VS30)
    # Fit Go,mo model ----

    # if(any(zm/Hs<=0) || any(Gm<=0) ){browser()}
    DATA <- (data.table(zm = zm, Gm = Gm))[Gm > 0][, .(X = log(zm / Hs), Y = log(Gm))]


    LM <- lm(formula = Y ~ X, data = DATA)
    LnGo <- LM$coefficients[1] |> unname()
    mo <- LM$coefficients[2] |>
        round(digits = 3) |>
        unname()
    Go <- exp(LnGo) |> round(digits = 2) # [MPa]


    # Gmax = rho* Vs^2
    # LM <- lm(formula=Z~X+Y,data=DATA)
    # VSo <- predict.lm(LM,newdata = data.table(X=0,Y=log(Go))) |> exp() |> round(digits = 1)
    VSo <- sqrt(9.81 * Go * 1000 / gsat[NL]) |> round(digits = 1)

    Ts <- fitModel.Ts(VSm = VSm, hs = hs, zm = zm) |> round(digits = 3)
    USCS.ID <- UID |>
        unique() |>
        sort() |>
        paste0(collapse = ".")
    # browser()
    Coarse <- round(sum(hs[UID %in% ValidCoarse]) / sum(hs) * 100, digits = 0)
    Gravels <- round(sum(hs[UID %in% ValidGravels]) / sum(hs) * 100, digits = 0)
    Fines <- round(sum(hs[UID %in% ValidFines]) / sum(hs) * 100, digits = 0)
    Sands <- round(sum(hs[UID %in% ValidSands]) / sum(hs) * 100, digits = 0)
    Silts <- round(sum(hs[UID %in% ValidSilts]) / sum(hs) * 100, digits = 0)
    Clays <- round(sum(hs[UID %in% ValidClays]) / sum(hs) * 100, digits = 0)
    Organic <- round(sum(hs[UID %in% ValidOrganic]) / sum(hs) * 100, digits = 0)

    MAP <- .mapLayers(UID = UID, hs = hs)

    Water <- round(100 * Water)

    # if(any(zm/Hs<=0) || any(Gm<=0) ){browser()}
    DATA <- (data.table(zm = zm, VSm = VSm))[VSm > 0][, .(X = log(zm / Hs), Z = log(VSm))]
    # DATA <- data.table(X=log(zm/Hs),Z=log(VSm))

    LM <- lm(formula = X ~ Z, data = DATA)
    Z500 <- predict(LM, newdata = data.table(Z = log(500))) |>
        exp() |>
        round(digits = 1)
    Z1000 <- predict(LM, newdata = data.table(Z = log(1000))) |>
        exp() |>
        round(digits = 1)
    # browser()
    SN <- paste0(
        str_pad(string = (round(10 * Ts, digits = 0)), width = 2, side = "left", pad = 0),
        str_pad(string = Hs, width = 3, side = "left", pad = 0), SID
    )

    # Output -----

    SiteProperties <- data.table(Hs, Hw, NL, Z500, Z1000, SID, Go, mo, Ts, VSo, VS30, UID = USCS.ID, Gravels, Sands, Fines, Clays, Silts, Organic, Water, MAP, Go_Units = "MPa", Vs_Units = "m/s", SN, POP, POP_Units = "kPa")

    if (getSiteLayers) {
        SiteLayers <- data.table(
            USCS = UID[1:NL],
            GroupID = GID[1:NL],
            zm = zm[1:NL],
            hs = hs[1:NL],
            emin = emin[1:NL] |> round(digits = 2),
            eo = eo[1:NL] |> round(digits = 2),
            emax = emax[1:NL] |> round(digits = 2),
            gsmin = gsmin[1:NL] |> round(digits = 1),
            gs = gs[1:NL] |> round(digits = 1),
            gsat = gsat[1:NL] |> round(digits = 1),
            gsmax = gsmax[1:NL] |> round(digits = 1),
            gw = gw[1:NL],
            Dr = Dr[1:NL] |> round(digits = 2),
            IP = IP[1:NL] |> round(digits = 2),
            LL = LL[1:NL] |> round(digits = 2),
            pm = pm[1:NL] |> round(digits = 2),
            OCR = OCR[1:NL] |> round(digits = 1),
            Gm = (Gm[1:NL]) |> round(digits = 2), # MPa
            VSm = VSm[1:NL] |> round(digits = 1),
            VSa = VSa[1:NL] |> round(digits = 1),
            SID = SID
        )

        # browser()
        return(list(SiteLayers = SiteLayers, SiteProperties = SiteProperties))
    } else {
        return(SiteProperties)
    }
}





#
# .tagSites <- function(.SD){
#   sapply(seq(1,nrow(.SD)), function(n){digest::digest(object = .SD[n],algo="crc32")})
# }
#


.setVs30 <- function(hs, VSm, VSo) {
    ZTop <- c(0, head(cumsum(hs), -1)) # depth to top of each layer
    ZBot <- ZTop + hs # depth to bottom of each layer

    h30 <- pmax(pmin(ZBot, 30) - ZTop, 0) # thickness that lies in 0–30 m
    t30 <- h30 / VSm # travel time of logged layers

    if (sum(hs) < 30) { # pad with foundation if needed
        t30 <- c(t30, (30 - sum(hs)) / VSo)
    }

    round(30 / sum(t30), 1)
}
