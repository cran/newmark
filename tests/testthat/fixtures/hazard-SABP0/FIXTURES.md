# Hazard fixtures of project AR-SABP0 (model SABP0, site caldera)

Built on 2026-10-13 by the oracle script fixturesHazard.R of the reconstruction chain, kept in the Git history of this repository, from the OpenQuake exports under the project's `data/hazard` (OpenQuake engine 3.26.2), read only, with R 4.6.1, data.table 1.18.6.1 and the installed newmark 1.7.1. The four archives, cut from output/classical and output/disaggregation of SABP0/caldera, were removed with the newmark OpenQuake readers, whose calculations the hazard product now produces; their cuts and SHA-256 are recorded in this file at Newmark f336176.

| Fixture | Bytes | SHA-256 | Source (under data/hazard) | Cut |
|---|---|---|---|---|
| AEPTable-1.7.1.Rds | 5237 | 8dac9f389e9043aa26c8815b60df6d736ae60bb5f631cee8833fdcdd29526816 | computed from the four archives cut from those exports | newmark 1.7.1 buildGMDP(path = <classical/ archives and the disagg_mag_dist archive staged together>, IDo = "SABP0", vref = 760, TRo = c(100, 200, 500, 1000, 2000, 2500, 5000, 10000))$AEPTable[p == "mean"] with siteID "caldera", 189 rows; the same recipe on the full cell reproduces uhs/SABP0/caldera/760/AEPTable.Rds (all.equal with tolerance 0 and attributes, identical column types) |
