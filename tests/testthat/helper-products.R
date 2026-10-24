# Projects for listProjectTables() and queryProjectTable(), whose product files are
# written directly at the destinations of processProject() from the slices
# of the AR-M2V4D products in the kmax fixture: the DnDraws, the PGA table
# read back from its CSV, the shear table and the kmaxTable rows of
# geometry S27 and material U4 at site CFD, under TR 100 and the MCE.
productSlice <- function() readRDS(test_path("fixtures", "kmax-M2V4D-S27-U4.Rds"))

PRODUCTCONTRACT <- list(path = list(tables = "products"), Vref.gmdp = 760, ID.gmdp = "M2V4D")

writeProductContract <- function(root, contract) {
  writeLines(jsonlite::toJSON(contract, digits = I(17), auto_unbox = TRUE), file.path(root, "newmark.json"))
}

# Normalized, because the product queries resolve their root in that
# spelling.
newProductRoot <- function() {
  Root <- tempfile("newmarkProducts")
  dir.create(Root, recursive = TRUE)
  normalizePath(Root, winslash = "/")
}

# The global ShearTable.Rds and kmaxTable.Rds and the DnDraws.Rds,
# PGATable.csv and kmaxTable.Rds of site CFD, written as the steps write
# them; no DnTable.Rds and no DnPlotTable.Rds.
productProject <- function(contract = PRODUCTCONTRACT) {
  Root <- newProductRoot()
  Slice <- productSlice()
  Site <- file.path(Root, "products", "M2V4D", "CFD")
  dir.create(Site, recursive = TRUE)
  saveRDS(Slice$shearTable, file.path(Root, "products", "ShearTable.Rds"))
  saveRDS(Slice$drawTable, file.path(Site, "DnDraws.Rds"))
  data.table::fwrite(Slice$PGATable, file.path(Site, "PGATable.csv"), na = "NA")
  saveRDS(Slice$kmaxTable, file.path(Site, "kmaxTable.Rds"))
  saveRDS(Slice$kmaxTable, file.path(Root, "products", "kmaxTable.Rds"))
  writeProductContract(Root, contract)
  Root
}

# The DnTable.Rds and PGATable.csv of site TST and the global DnTable.Rds
# of the product-pairs case of the dn reference (AM88 and BM19 with their
# trees, geometries G1 and G2, under TR 100 and the MCE), written as the
# DnTable step writes them.
dnProductProject <- function() {
  Root <- newProductRoot()
  Products <- readRDS(test_path("fixtures", "dn-reference.rds"))$pairs$products
  Site <- file.path(Root, "products", "TEST", "TST")
  dir.create(Site, recursive = TRUE)
  saveRDS(Products$DnTable, file.path(Site, "DnTable.Rds"))
  writeLines(Products$PGATable, file.path(Site, "PGATable.csv"))
  saveRDS(Products$GlobalDnTable, file.path(Root, "products", "DnTable.Rds"))
  writeProductContract(Root, list(path = list(tables = "products"), Vref.gmdp = 760, ID.gmdp = "TEST"))
  Root
}

# Every file and directory under a root, with its size and modification
# time, so that a call that writes, removes or touches anything changes it.
productTree <- function(root) {
  FILES <- list.files(root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE, no.. = TRUE)
  data.frame(file = FILES, size = file.size(file.path(root, FILES)), mtime = file.mtime(file.path(root, FILES)))
}
