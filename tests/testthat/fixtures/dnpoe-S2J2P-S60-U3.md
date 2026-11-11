# DnPOE fixture

Source: AR-S2J2P, `newmark.json`, `data/newmark/ShearTable.Rds`,
`data/newmark/ZAF/Doornpoort/DnDraws.Rds` and `_local/data/DnPOETable.Rds`.
Extracted 2026-11-01 without modifying the project. Oracle SHA-256:
`d5c8ccc99dda9ad828763324b93a6d2aa75deb4a5ff062be9e4b9ada6e33fcef`.

`dnpoe-S2J2P-S60-U3.Rds` holds the 30000 draws of Doornpoort/MCE/S60/U3,
the shear rows of that pair, NS 200, modelID ZAF, MCE alone (`TR = numeric()`),
the five active unit-weight models, the contract's seven DaH thresholds,
and the corresponding 210 independently produced oracle rows. No draws
were regenerated. The first four ky at DaH 1 give 0.764, 0.740, 0.706,
0.673; all rows have n 1000. Tests compare column values, types and order
exactly; the oracle has no data.table key.
