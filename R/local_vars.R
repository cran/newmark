#' @importFrom stats weighted.mean
NULL

# Work around R CMD check "no visible binding" notes
utils::globalVariables(
  c("..ID", "AF", "w", "pga", "PGA", "N", "depth", "lat", "lon",
    "SaRock", "TnVal", "FValue", "SaF", "LnF", "muLnF", "sdLnF",
    "pnum", "lnSa", "Tn", "p", "Dn", "weight", "muLnD", "sdLnD",
    "LnD", "LnSaF", "muLnSaF", "sdLnSaF", "ID", "Sa", "SaF_SS17",
    "TR", "Vs30", "Vref", ".", "AEP", "IT", "POE", "iml", "imt",
    "ky", "s", "IDg", "IDn", "K", "kmax",
    "k", "mu", "i.mu0", "i.s0", "SiteAmp2020",
    # dsra globals
    "ValidSands", "ValidGravels", "ValidSilts", "ValidClays",
    "ValidOrganic", "ValidFines", "ValidCoarse", "ValidUSCS",
    "ValidGroups", "SiteTable", "VoidRatiosUSCS",
    "RelativeDensityRanges", "UnitWeightRanges", "ParticleSize",
    "USCS", "ShearModelParameters", "CylinderRoots", "SiteClass",
    "DensityID", "GroupID", "ModelID", "lo", "mo", "Hs", "Hw",
    "VSo", "Go", "an", "Fines", "Gravels", "Sands", "Water",
    "Clays", "Silts", "l", "m", "n",
    "Ts", "VS30", "Z500", "Z1000"
  )
)
